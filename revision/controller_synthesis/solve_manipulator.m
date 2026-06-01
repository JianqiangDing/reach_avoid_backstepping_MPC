% Phase 1.3 -- manipulator controller synthesis: constrained + unconstrained.
%
% Counterpart of solve_dubins.m. At the SAME settings (from the Phase 1.2 outputs:
% per-channel u_max_eff = [tau1, tau2], kept samples S_filtered, mu/ds/dv, physical
% constants), synthesizes TWO controllers for later comparison:
%   * CONSTRAINED   -- hard bounded-control SOP over S_filtered (per-channel torque
%                      limits enforced).
%   * UNCONSTRAINED -- vanilla k1 at the same lambda, NO per-sample bound
%                      constraints (baseline reach-avoid law ignoring torque limits).
%
% Exports (revision/controllers/):
%   manipulator_constrained.py     u_opt + certificate (bounded)
%   manipulator_unconstrained.py   u_opt + certificate (vanilla, unbounded)
% The 5-D augmented synthesis state x = [q1,q2,dq1,dq2,q1+q2] is collapsed
% (x5 = x1+x2) to a 4-D export for both, as in the main example_manipulator.m.

function solve_manipulator()
    here = fileparts(mfilename('fullpath'));        % revision/controller_synthesis
    revision_dir = fileparts(here);                 % revision
    addpath(fullfile(revision_dir, 'matlab_frozen'));
    data_dir = fullfile(revision_dir, 'data');
    ctrl_dir = fullfile(revision_dir, 'controllers');
    if ~exist(ctrl_dir, 'dir'); mkdir(ctrl_dir); end

    t_all = tic;

    % ---- Phase 1.2 outputs (the interface) ----------------------------------
    p12 = jsondecode(fileread(fullfile(data_dir, 'phase1_2_outputs_manipulator.json')));
    mu_val = p12.probe.mu;
    xi0    = p12.probe.xi0;
    ds     = p12.probe.ds;
    dv     = p12.probe.dv;
    ub1    = p12.u_max_eff.tau1;
    ub2    = p12.u_max_eff.tau2;
    c      = p12.system_constants;   % m1 m2 l1 l2 lc1 lc2 I1 I2 g
    fprintf('solve_manipulator: mu=%g xi0=%g ds=%d dv=%d  u_max_eff=[tau1 %.4g, tau2 %.4g]\n', ...
        mu_val, xi0, ds, dv, ub1, ub2);

    % ---- kept samples S_filtered from Phase 1.2 (5-D on manifold) ------------
    Tk = readtable(fullfile(data_dir, 'kept_samples_manipulator.csv'));
    x_samples_valid = [Tk.x1, Tk.x2, Tk.x3, Tk.x4, Tk.x5]';   % 5 x N
    n_samples = size(x_samples_valid, 2);
    fprintf('  loaded %d kept samples (S_filtered)\n', n_samples);

    % ---- manipulator system (frozen example def; constants from Phase 1.2) ---
    syms x1 x2 x3 x4 x5 y1 y2;
    x_vars = [x1; x2; x3; x4; x5];
    y_vars = [y1; y2];
    m1 = c.m1; m2 = c.m2; l1 = c.l1; l2 = c.l2; lc1 = c.lc1; lc2 = c.lc2; I1 = c.I1; I2 = c.I2; g = c.g;
    M11 = I1 + I2 + m1 * lc1^2 + m2 * (l1^2 + lc2^2 + 2 * l1 * lc2 * cos(x2));
    M12 = m2 * (lc2^2 + l1 * lc2 * cos(x2)) + I2;
    M22 = m2 * lc2^2 + I2;
    M = [M11, M12; M12, M22];
    C_mat = [-m2 * l1 * lc2 * sin(x2) * x4, -m2 * l1 * lc2 * sin(x2) * (x3 + x4); ...
        m2 * l1 * lc2 * sin(x2) * x3, 0];
    G_vec = [(m1 * g * lc1 + m2 * g * l1) * cos(x1) + m2 * g * lc2 * cos(x5); ...
        m2 * g * lc2 * cos(x5)];
    Minv = inv(M);
    qddot_f = simplify(Minv * (-C_mat * [x3; x4] - G_vec));
    fx_sym = [x3; x4; qddot_f; x3 + x4];
    gx_sym = [zeros(2, 2); Minv; zeros(1, 2)];
    hx_sym = [l1 * cos(x1) + l2 * cos(x5); l1 * sin(x1) + l2 * sin(x5)];
    safe_set = -((4 * (y1 - 2) - 2 * y2^3)^2) + 0.8 * y2^3 + 10;
    target_set = ((y1 - 2 - 3.5)^2 / 1.2^2) + ((y2 - 1.8)^2 / 0.4^2) - 2;

    t_stage = tic;
    fprintf('[stage] building symbolic manipulator dynamics done; running backstepping ...\n');
    [u, k1, J_k1, mu, lambda, certificate, cert_term_dict, ~, ~, ~, p, r_deg] = ...
        reach_avoid_controller(fx_sym, gx_sym, hx_sym, x_vars, y_vars, safe_set);
    fprintf('[stage] backstepping design done in %.1fs; solving vanilla k1 ...\n', toc(t_stage)); t_stage = tic;
    [k1_y, k1_lambda, ~] = solve_vanilla_k1_controller_xi(y_vars, safe_set, target_set, dv, ds, xi0);
    J_k1_y = jacobian(k1_y, y_vars);
    fprintf('[stage] vanilla k1 done in %.1fs (lambda=%.4g)\n', toc(t_stage), k1_lambda);

    % ---- UNCONSTRAINED controller (vanilla k1, NO per-sample bound constraints)
    u_unc = subs(subs(subs(u, k1, k1_y), J_k1, J_k1_y), y_vars, hx_sym);
    u_unc = sub_lambda(sub_mu(u_unc, mu, mu_val), lambda, k1_lambda);
    u_unc = subs(u_unc, x5, x1 + x2);
    cert_unc = sub_lambda(sub_mu(subs(subs(certificate, k1, k1_y), y_vars, hx_sym), mu, mu_val), lambda, k1_lambda);
    cert_unc = subs(cert_unc, x5, x1 + x2);
    params_unc = struct('example', 'manipulator', 'controller_type', 'unconstrained_vanilla', ...
        'mu_val', mu_val, 'xi0', xi0, 'ds', ds, 'dv', dv, 'lambda', k1_lambda);
    unc_path = fullfile(ctrl_dir, 'manipulator_unconstrained.py');
    export_to_python(u_unc, cert_unc, k1_y, params_unc, unc_path);
    fprintf('  exported %s\n', unc_path);

    % ---- CONSTRAINED controller (HARD bounded-control SOP, no slack) --------
    fprintf('[stage] solving HARD bounded SOP over %d samples ...\n', n_samples); t_stage = tic;
    lb = [-ub1; -ub2]; ub = [ub1; ub2];
    ux_for_sop = sub_lambda(sub_mu(subs(u, y_vars, hx_sym), mu, mu_val), lambda, k1_lambda);
    [k1_opt, J_k1_opt, k1_delta] = solve_k1_controller_sop( ...
        ux_for_sop, k1, J_k1, cert_term_dict, p, r_deg, ...
        x_vars, y_vars, hx_sym, safe_set, target_set, x_samples_valid, ...
        lb, ub, ds, dv, k1_lambda, mu_val);
    fprintf('[stage] hard SOP done in %.1fs (delta=%.4g)\n', toc(t_stage), k1_delta);

    u_con = subs(subs(subs(u, k1, k1_opt), J_k1, J_k1_opt), y_vars, hx_sym);
    u_con = sub_lambda(sub_mu(u_con, mu, mu_val), lambda, k1_lambda);
    u_con = subs(u_con, x5, x1 + x2);
    cert_con = sub_lambda(sub_mu(subs(subs(certificate, k1, k1_opt), y_vars, hx_sym), mu, mu_val), lambda, k1_lambda);
    cert_con = subs(cert_con, x5, x1 + x2);
    params_con = struct('example', 'manipulator', 'controller_type', 'constrained_hard_sop', ...
        'u_max_eff', [ub1; ub2], 'mu_val', mu_val, 'xi0', xi0, 'ds', ds, 'dv', dv, ...
        'n_samples', n_samples, 'k1_delta', k1_delta, 'lambda', k1_lambda);
    con_path = fullfile(ctrl_dir, 'manipulator_constrained.py');
    export_to_python(u_con, cert_con, k1_opt, params_con, con_path);
    fprintf('  exported %s\n', con_path);

    % ---- Phase 1.3 synthesis meta -------------------------------------------
    meta = struct('example', 'manipulator', ...
        'controller_constrained_py',   fullfile('revision', 'controllers', 'manipulator_constrained.py'), ...
        'controller_unconstrained_py', fullfile('revision', 'controllers', 'manipulator_unconstrained.py'), ...
        'u_max_eff', struct('tau1', ub1, 'tau2', ub2), ...
        'mu', mu_val, 'xi0', xi0, 'ds', ds, 'dv', dv, 'n_samples', n_samples, ...
        'k1_delta', k1_delta, 'k1_lambda', k1_lambda, 'wallclock_s', toc(t_all));
    fid = fopen(fullfile(data_dir, 'phase1_3_meta_manipulator.json'), 'w');
    fwrite(fid, jsonencode(meta, 'PrettyPrint', true));
    fclose(fid);
    fprintf('solve_manipulator: done in %.1fs (delta=%.4g, %d samples)\n', toc(t_all), k1_delta, n_samples);
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
