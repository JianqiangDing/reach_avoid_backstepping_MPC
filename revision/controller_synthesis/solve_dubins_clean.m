% Phase 1.3 -- Dubins final controller synthesis (HARD bounded SOP, no slack).
%
% Reads the Phase 1.2 outputs (the only interface): the per-channel effective
% bound u_max_eff and the kept sample subset S_filtered. Re-solves the
% bounded-control SOP with HARD input limits over S_filtered (frozen
% solve_k1_controller_sop.m), then exports the closed-form controller +
% reach-avoid certificate to revision/controllers/dubins_clean.py.
%
% Self-contained: addpath only revision/matlab_frozen; inputs read from
% revision/data/phase1_2_outputs_dubins.json + kept_samples_dubins.csv; outputs
% to revision/controllers/ and revision/data/. No notebook ENV is needed -- all
% parameters come from the Phase 1.2 outputs (which the notebook produced).

function solve_dubins_clean()
    here = fileparts(mfilename('fullpath'));        % revision/controller_synthesis
    revision_dir = fileparts(here);                 % revision
    addpath(fullfile(revision_dir, 'matlab_frozen'));
    data_dir = fullfile(revision_dir, 'data');
    ctrl_dir = fullfile(revision_dir, 'controllers');
    if ~exist(ctrl_dir, 'dir'); mkdir(ctrl_dir); end

    t_all = tic;

    % ---- Phase 1.2 outputs (the interface) ----------------------------------
    p12 = jsondecode(fileread(fullfile(data_dir, 'phase1_2_outputs_dubins.json')));
    mu_val   = p12.probe.mu;
    xi0      = p12.probe.xi0;
    ds       = p12.probe.ds;
    dv       = p12.probe.dv;
    ub_omega = p12.u_max_eff.omega;
    ub_a     = p12.u_max_eff.a;
    fprintf('solve_dubins_clean: mu=%g xi0=%g ds=%d dv=%d  u_max_eff=[omega %.4g, a %.4g]\n', ...
        mu_val, xi0, ds, dv, ub_omega, ub_a);

    % ---- kept samples S_filtered from Phase 1.2 -----------------------------
    Tk = readtable(fullfile(data_dir, 'kept_samples_dubins.csv'));
    x_samples_valid = [Tk.x1, Tk.x2, Tk.th, Tk.v]';   % 4 x N
    n_samples = size(x_samples_valid, 2);
    fprintf('  loaded %d kept samples (S_filtered)\n', n_samples);

    % ---- Dubins system (frozen example def; identical to solve_dubins_slack) -
    syms x1 x2 th v y1 y2;
    x_vars = [x1; x2; th; v];
    y_vars = [y1; y2];
    fx_sym = [v * cos(th); v * sin(th); 0; 0];
    gx_sym = [0, 0; 0, 0; 1, 0; 0, 1];
    hx_sym = [x1; x2];
    h_raw  = -(y1^4 + y2^4 - 16) * (y1^4 + y2^4 - 4);
    target_set = (y2 - 0)^2 + ((y1 + 1.7) / 0.5)^2 - 0.4;
    safe_set = 1e-3 * (-target_set + 300) * h_raw;

    t_stage = tic;
    [u, k1, J_k1, mu, lambda, certificate, cert_term_dict, ~, ~, ~, p, r_deg] = ...
        reach_avoid_controller(fx_sym, gx_sym, hx_sym, x_vars, y_vars, safe_set);
    fprintf('[stage] backstepping design done in %.1fs; solving vanilla k1 (for lambda) ...\n', toc(t_stage)); t_stage = tic;
    [~, k1_lambda, ~] = solve_vanilla_k1_controller_xi(y_vars, safe_set, target_set, dv, ds, xi0);
    fprintf('[stage] vanilla k1 done in %.1fs (lambda=%.4g); solving HARD bounded SOP over %d samples ...\n', ...
        toc(t_stage), k1_lambda, n_samples); t_stage = tic;

    % ---- HARD bounded-control SOP (no slack) --------------------------------
    lb = [-ub_omega; -ub_a]; ub = [ub_omega; ub_a];
    ux_for_sop = sub_lambda(sub_mu(subs(u, y_vars, hx_sym), mu, mu_val), lambda, k1_lambda);
    [k1_opt, J_k1_opt, k1_delta] = solve_k1_controller_sop( ...
        ux_for_sop, k1, J_k1, cert_term_dict, p, r_deg, ...
        x_vars, y_vars, hx_sym, safe_set, target_set, x_samples_valid, ...
        lb, ub, ds, dv, k1_lambda, mu_val);
    fprintf('[stage] hard SOP done in %.1fs (delta=%.4g)\n', toc(t_stage), k1_delta);

    % ---- assemble final controller + certificate ----------------------------
    u_opt = subs(u, k1, k1_opt);
    u_opt = subs(u_opt, J_k1, J_k1_opt);
    u_opt = subs(u_opt, y_vars, hx_sym);
    u_opt = sub_lambda(sub_mu(u_opt, mu, mu_val), lambda, k1_lambda);

    certificate_opt = subs(certificate, k1, k1_opt);
    certificate_opt = subs(certificate_opt, y_vars, hx_sym);
    certificate_opt = sub_lambda(sub_mu(certificate_opt, mu, mu_val), lambda, k1_lambda);

    % ---- export to revision/controllers/ ------------------------------------
    params = struct();
    params.example = 'dubins';
    params.u_max_eff = [ub_omega; ub_a];
    params.mu_val = mu_val;
    params.ds = ds;
    params.dv = dv;
    params.n_samples = n_samples;
    params.k1_delta = k1_delta;
    py_path = fullfile(ctrl_dir, 'dubins_clean.py');
    export_to_python(u_opt, certificate_opt, k1_opt, params, py_path);
    fprintf('  exported %s\n', py_path);

    % ---- Phase 1.3 synthesis meta (registry for the notebook) ---------------
    meta = struct('example', 'dubins', ...
        'controller_py', fullfile('revision', 'controllers', 'dubins_clean.py'), ...
        'u_max_eff', struct('omega', ub_omega, 'a', ub_a), ...
        'mu', mu_val, 'ds', ds, 'dv', dv, 'n_samples', n_samples, ...
        'k1_delta', k1_delta, 'k1_lambda', k1_lambda, 'wallclock_s', toc(t_all));
    fid = fopen(fullfile(data_dir, 'phase1_3_meta_dubins.json'), 'w');
    fwrite(fid, jsonencode(meta, 'PrettyPrint', true));
    fclose(fid);
    fprintf('solve_dubins_clean: done in %.1fs (delta=%.4g, %d samples)\n', toc(t_all), k1_delta, n_samples);
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
