% Phase 1.3 -- Dubins controller synthesis: constrained + unconstrained.
%
% At the SAME settings (mu, xi0, ds, dv, lambda) reads from the Phase 1.2 outputs,
% synthesizes TWO controllers for later comparison:
%   * CONSTRAINED   -- hard bounded-control SOP over the Phase 1.2 kept samples
%                      (frozen solve_k1_controller_sop.m), input limits enforced.
%                      Solved at the SAVED Phase 1.2 lambda (p12.probe.lambda),
%                      not a freshly re-solved vanilla lambda (deterministic, and
%                      consistent with the lambda the kept samples were derived under).
%   * UNCONSTRAINED -- the vanilla k1 (solve_vanilla_k1_controller_xi); its own
%                      matched lambda, with NO per-sample control-bound constraints.
%                      Baseline reach-avoid law that ignores actuator limits.
%
% Exports (revision/controllers/):
%   dubins_constrained.py     u_opt + certificate (bounded)
%   dubins_unconstrained.py   u_opt + certificate (vanilla, unbounded)
% Meta: revision/data/phase1_3_meta_dubins.json (both controllers).
%
% Self-contained: addpath only revision/matlab_frozen; inputs from
% revision/data/phase1_2_outputs_dubins.json + kept_samples_dubins.csv.

function solve_dubins()
    here = fileparts(mfilename('fullpath')); % revision/controller_synthesis
    revision_dir = fileparts(here); % revision
    addpath(fullfile(revision_dir, 'matlab_frozen'));
    data_dir = fullfile(revision_dir, 'data');
    ctrl_dir = fullfile(revision_dir, 'controllers');
    if ~exist(ctrl_dir, 'dir'); mkdir(ctrl_dir); end

    t_all = tic;

    % ---- Phase 1.2 outputs (the interface) ----------------------------------
    p12 = jsondecode(fileread(fullfile(data_dir, 'phase1_2_outputs_dubins.json')));
    mu_val = p12.probe.mu;
    xi0 = p12.probe.xi0;
    ds = p12.probe.ds;
    dv = p12.probe.dv;
    lambda_saved = p12.probe.lambda; % SAVED Phase 1.2 backstepping scale; the
    % constrained SOP is solved at THIS lambda
    % (deterministic, consistent with the kept
    % samples / u_max_eff) instead of re-solving.
    ub_omega = p12.u_max_eff.omega;
    ub_a = p12.u_max_eff.a;
    fprintf('solve_dubins: mu=%g xi0=%g ds=%d dv=%d lambda_saved=%.6g  u_max_eff=[omega %.4g, a %.4g]\n', ...
        mu_val, xi0, ds, dv, lambda_saved, ub_omega, ub_a);

    % ---- kept samples S_filtered from Phase 1.2 -----------------------------
    Tk = readtable(fullfile(data_dir, 'kept_samples_dubins.csv'));
    x_samples_valid = [Tk.x1, Tk.x2, Tk.th, Tk.v]'; % 4 x N
    n_samples = size(x_samples_valid, 2);
    fprintf('  loaded %d kept samples (S_filtered)\n', n_samples);

    % ---- Dubins system (frozen example def) ---------------------------------
    syms x1 x2 th v y1 y2;
    x_vars = [x1; x2; th; v];
    y_vars = [y1; y2];
    fx_sym = [v * cos(th); v * sin(th); 0; 0];
    gx_sym = [0, 0; 0, 0; 1, 0; 0, 1];
    hx_sym = [x1; x2];
    h_raw =- (y1 ^ 4 + y2 ^ 4 - 16) * (y1 ^ 4 + y2 ^ 4 - 4);
    target_set = (y2 - 0) ^ 2 + ((y1 + 1.7) / 0.5) ^ 2 - 0.4;
    safe_set = 1e-3 * (-target_set + 300) * h_raw;

    t_stage = tic;
    [u, k1, J_k1, mu, lambda, certificate, cert_term_dict, ~, ~, ~, p, r_deg] = ...
        reach_avoid_controller(fx_sym, gx_sym, hx_sym, x_vars, y_vars, safe_set);
    fprintf('[stage] backstepping design done in %.1fs; solving vanilla k1 ...\n', toc(t_stage)); t_stage = tic;
    [k1_y, k1_lambda, ~] = solve_vanilla_k1_controller_xi(y_vars, safe_set, target_set, dv, ds, xi0);
    J_k1_y = jacobian(k1_y, y_vars);
    fprintf('[stage] vanilla k1 done in %.1fs (lambda=%.4g)\n', toc(t_stage), k1_lambda);

    % ---- UNCONSTRAINED controller (vanilla k1, NO per-sample bound constraints)
    u_unc = subs(subs(subs(u, k1, k1_y), J_k1, J_k1_y), y_vars, hx_sym);
    u_unc = sub_lambda(sub_mu(u_unc, mu, mu_val), lambda, k1_lambda);
    cert_unc = sub_lambda(sub_mu(subs(subs(certificate, k1, k1_y), y_vars, hx_sym), mu, mu_val), lambda, k1_lambda);
    params_unc = struct('example', 'dubins', 'controller_type', 'unconstrained_vanilla', ...
        'mu_val', mu_val, 'xi0', xi0, 'ds', ds, 'dv', dv, 'lambda', k1_lambda);
    unc_path = fullfile(ctrl_dir, 'dubins_unconstrained.py');
    export_to_python(u_unc, cert_unc, k1_y, params_unc, unc_path);
    fprintf('  exported %s\n', unc_path);

    % ---- CONSTRAINED controller (HARD bounded-control SOP, no slack) --------
    % Solved at the SAVED Phase 1.2 lambda (lambda_saved), NOT the freshly
    % re-solved vanilla k1_lambda, so the constrained synthesis is deterministic
    % and uses the same backstepping scale the kept samples / u_max_eff were
    % derived under.
    fprintf('[stage] solving HARD bounded SOP over %d samples (lambda=%.6g) ...\n', n_samples, lambda_saved); t_stage = tic;
    lb = [-ub_omega; -ub_a]; ub = [ub_omega; ub_a];
    ux_for_sop = sub_lambda(sub_mu(subs(u, y_vars, hx_sym), mu, mu_val), lambda, lambda_saved);
    [k1_opt, J_k1_opt, k1_delta] = solve_k1_controller_sop( ...
        ux_for_sop, k1, J_k1, cert_term_dict, p, r_deg, ...
        x_vars, y_vars, hx_sym, safe_set, target_set, x_samples_valid, ...
        lb, ub, ds, dv, lambda_saved, mu_val);
    fprintf('[stage] hard SOP done in %.1fs (delta=%.4g)\n', toc(t_stage), k1_delta);

    u_con = subs(subs(subs(u, k1, k1_opt), J_k1, J_k1_opt), y_vars, hx_sym);
    u_con = sub_lambda(sub_mu(u_con, mu, mu_val), lambda, lambda_saved);
    cert_con = sub_lambda(sub_mu(subs(subs(certificate, k1, k1_opt), y_vars, hx_sym), mu, mu_val), lambda, lambda_saved);
    params_con = struct('example', 'dubins', 'controller_type', 'constrained_hard_sop', ...
        'u_max_eff', [ub_omega; ub_a], 'mu_val', mu_val, 'xi0', xi0, 'ds', ds, 'dv', dv, ...
        'n_samples', n_samples, 'k1_delta', k1_delta, 'lambda', lambda_saved);
    con_path = fullfile(ctrl_dir, 'dubins_constrained.py');
    export_to_python(u_con, cert_con, k1_opt, params_con, con_path);
    fprintf('  exported %s\n', con_path);

    % ---- Phase 1.3 synthesis meta (registry for the notebook) ---------------
    meta = struct('example', 'dubins', ...
        'controller_constrained_py', fullfile('revision', 'controllers', 'dubins_constrained.py'), ...
        'controller_unconstrained_py', fullfile('revision', 'controllers', 'dubins_unconstrained.py'), ...
        'u_max_eff', struct('omega', ub_omega, 'a', ub_a), ...
        'mu', mu_val, 'xi0', xi0, 'ds', ds, 'dv', dv, 'n_samples', n_samples, ...
        'k1_delta', k1_delta, 'k1_lambda', k1_lambda, 'lambda_constrained', lambda_saved, ...
        'wallclock_s', toc(t_all));
    fid = fopen(fullfile(data_dir, 'phase1_3_meta_dubins.json'), 'w');
    fwrite(fid, jsonencode(meta, 'PrettyPrint', true));
    fclose(fid);
    fprintf('solve_dubins: done in %.1fs (delta=%.4g, %d samples)\n', toc(t_all), k1_delta, n_samples);
end

% =============================================================================
function expr = sub_mu(expr, mu_cells, mu_val)

    for i = 1:numel(mu_cells)

        if ~isempty(mu_cells{i})
            expr = subs(expr, mu_cells{i}, mu_val);
        end

    end

end

function expr = sub_lambda(expr, lambda_sym, lambda_val)
    expr = subs(expr, lambda_sym, lambda_val);
end
