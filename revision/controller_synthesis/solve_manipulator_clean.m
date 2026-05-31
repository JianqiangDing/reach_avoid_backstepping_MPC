% Phase 1.3 -- manipulator final controller synthesis (HARD bounded SOP, no slack).
%
% Counterpart of solve_dubins_clean.m. Reads the Phase 1.2 outputs (per-channel
% effective bound u_max_eff = [tau1, tau2], kept samples S_filtered, mu/ds/dv and
% the physical constants), re-solves the bounded-control SOP with HARD per-channel
% torque limits over S_filtered, and exports the controller + certificate to
% revision/controllers/manipulator_clean.py. The 5-D augmented synthesis state
% x = [q1,q2,dq1,dq2,q1+q2] is collapsed (x5 = x1+x2) to a 4-D export, as in the
% main example_manipulator.m.

function solve_manipulator_clean()
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
    fprintf('solve_manipulator_clean: mu=%g xi0=%g ds=%d dv=%d  u_max_eff=[tau1 %.4g, tau2 %.4g]\n', ...
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
    fprintf('[stage] backstepping design done in %.1fs; solving vanilla k1 (for lambda) ...\n', toc(t_stage)); t_stage = tic;
    [~, k1_lambda, ~] = solve_vanilla_k1_controller_xi(y_vars, safe_set, target_set, dv, ds, xi0);
    fprintf('[stage] vanilla k1 done in %.1fs (lambda=%.4g); solving HARD bounded SOP over %d samples ...\n', ...
        toc(t_stage), k1_lambda, n_samples); t_stage = tic;

    % ---- HARD bounded-control SOP (no slack) --------------------------------
    lb = [-ub1; -ub2]; ub = [ub1; ub2];
    ux_for_sop = sub_lambda(sub_mu(subs(u, y_vars, hx_sym), mu, mu_val), lambda, k1_lambda);
    [k1_opt, J_k1_opt, k1_delta] = solve_k1_controller_sop( ...
        ux_for_sop, k1, J_k1, cert_term_dict, p, r_deg, ...
        x_vars, y_vars, hx_sym, safe_set, target_set, x_samples_valid, ...
        lb, ub, ds, dv, k1_lambda, mu_val);
    fprintf('[stage] hard SOP done in %.1fs (delta=%.4g)\n', toc(t_stage), k1_delta);

    % ---- assemble final controller + certificate, collapse x5 = x1 + x2 ------
    u_opt = subs(u, k1, k1_opt);
    u_opt = subs(u_opt, J_k1, J_k1_opt);
    u_opt = subs(u_opt, y_vars, hx_sym);
    u_opt = sub_lambda(sub_mu(u_opt, mu, mu_val), lambda, k1_lambda);
    u_opt = subs(u_opt, x5, x1 + x2);

    certificate_opt = subs(certificate, k1, k1_opt);
    certificate_opt = subs(certificate_opt, y_vars, hx_sym);
    certificate_opt = sub_lambda(sub_mu(certificate_opt, mu, mu_val), lambda, k1_lambda);
    certificate_opt = subs(certificate_opt, x5, x1 + x2);

    % ---- export to revision/controllers/ ------------------------------------
    params = struct();
    params.example = 'manipulator';
    params.u_max_eff = [ub1; ub2];
    params.mu_val = mu_val;
    params.ds = ds;
    params.dv = dv;
    params.n_samples = n_samples;
    params.k1_delta = k1_delta;
    py_path = fullfile(ctrl_dir, 'manipulator_clean.py');
    export_to_python(u_opt, certificate_opt, k1_opt, params, py_path);
    fprintf('  exported %s\n', py_path);

    % ---- Phase 1.3 synthesis meta -------------------------------------------
    meta = struct('example', 'manipulator', ...
        'controller_py', fullfile('revision', 'controllers', 'manipulator_clean.py'), ...
        'u_max_eff', struct('tau1', ub1, 'tau2', ub2), ...
        'mu', mu_val, 'ds', ds, 'dv', dv, 'n_samples', n_samples, ...
        'k1_delta', k1_delta, 'k1_lambda', k1_lambda, 'wallclock_s', toc(t_all));
    fid = fopen(fullfile(data_dir, 'phase1_3_meta_manipulator.json'), 'w');
    fwrite(fid, jsonencode(meta, 'PrettyPrint', true));
    fclose(fid);
    fprintf('solve_manipulator_clean: done in %.1fs (delta=%.4g, %d samples)\n', toc(t_all), k1_delta, n_samples);
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
