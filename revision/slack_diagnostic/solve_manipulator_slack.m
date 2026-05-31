% Phase 1.2 -- 2-DoF manipulator per-sample slack diagnostic (self-contained).
%
% Same idea as solve_dubins_slack: request a tight per-channel torque bound and
% solve the per-sample slack SOP (each sample j, channel i gets s_{j,i} >= 0 with
% |u_i(x_j)| <= ub_i + s_{j,i}, reach-avoid kept hard, min delta + w*sum s_{j,i}).
% The per-sample slacks form a distribution we inspect (manipulator_slack_diag.ipynb)
% to pick per-channel effective bounds and the kept subset.
%
% Manipulator specifics:
%   * 5-D augmented state x = [q1; q2; dq1; dq2; q1+q2]; x5 = x1 + x2 is a manifold.
%     Sampling is done on the manifold (sample x1..x4, set x5 = x1+x2) -- consistent
%     with Phase 1.1's 4-D sampling, NOT the independent x5-box of the old test.
%   * channels have very different scales (tau1 ~ thousands, tau2 ~ hundreds), so
%     per-channel ub (UB1, UB2) and per-channel effective bounds.
%
% Self-contained: addpath only revision/matlab_frozen; X_S_eff from
% revision/data/phase1_1_outputs_manipulator.json; outputs to revision/data/.
%
% ALL parameters are supplied by the calling notebook via environment variables
% (this .m sets no parameter values): MU, XI0, UB1, UB2, N_VALID, SLACK_WEIGHT, DS, DV.

function solve_manipulator_slack()
    here = fileparts(mfilename('fullpath'));        % revision/slack_diagnostic
    revision_dir = fileparts(here);                 % revision
    addpath(fullfile(revision_dir, 'matlab_frozen'));
    data_dir = fullfile(revision_dir, 'data');
    if ~exist(data_dir, 'dir'); mkdir(data_dir); end

    mu_val       = read_env_num('MU');
    xi0          = read_env_num('XI0');
    ub1          = read_env_num('UB1');
    ub2          = read_env_num('UB2');
    n_valid      = read_env_num('N_VALID');
    slack_weight = read_env_num('SLACK_WEIGHT');
    ds           = read_env_num('DS');
    dv           = read_env_num('DV');

    % ---- relaxed synthesis region X_S_eff from Phase 1.1 --------------------
    j = jsondecode(fileread(fullfile(data_dir, 'phase1_1_outputs_manipulator.json')));
    lo = j.X_S_eff_def.lower(:);   % [x1; x2; x3; x4; x5]
    hi = j.X_S_eff_def.upper(:);
    box_lo = lo(1:4); box_hi = hi(1:4);   % sample (x1..x4); x5 = x1 + x2 on the manifold
    fprintf('solve_manipulator_slack: mu=%g xi0=%g ub=[%g %g] n_valid=%d slack_weight=%g ds=%d dv=%d\n', ...
        mu_val, xi0, ub1, ub2, n_valid, slack_weight, ds, dv);
    fprintf('  X_S_eff box (x1..x4, from phase1_1_outputs_manipulator.json):\n');
    fprintf('    lower = [%s]\n    upper = [%s]   (x5 = x1 + x2)\n', ...
        strtrim(sprintf('%.4g ', box_lo)), strtrim(sprintf('%.4g ', box_hi)));

    % ---- manipulator system (5-D augmented; x5 independent symbol) ----------
    syms x1 x2 x3 x4 x5 y1 y2;
    x_vars = [x1; x2; x3; x4; x5];
    y_vars = [y1; y2];
    m1 = 1.0; m2 = 1.0; l1 = 4.0; l2 = 4.0; lc1 = 2.0; lc2 = 2.0; I1 = 0.02; I2 = 0.02; g = 9.81;
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

    [u, k1, J_k1, mu, lambda, certificate, cert_term_dict, ~, ~, ~, p, r_deg] = ...
        reach_avoid_controller(fx_sym, gx_sym, hx_sym, x_vars, y_vars, safe_set);

    [k1_y, k1_lambda, ~] = solve_vanilla_k1_controller_xi(y_vars, safe_set, target_set, dv, ds, xi0);

    % ---- sample n_valid reach-avoid states ON THE MANIFOLD x5 = x1 + x2 -----
    cert_vanilla = sub_lambda(sub_mu(subs(subs(certificate, k1, k1_y), y_vars, hx_sym), mu, mu_val), lambda, k1_lambda);
    rng(42);
    x_samples_valid = manifold_sample(n_valid, x_vars, y_vars, hx_sym, safe_set, target_set, ...
        cert_vanilla, box_lo, box_hi);
    valid_count = size(x_samples_valid, 2);
    fprintf('  sampled %d valid reach-avoid states on the manifold\n', valid_count);

    % ---- per-sample slack SOP -----------------------------------------------
    lb = [-ub1; -ub2]; ub = [ub1; ub2];
    ux_for_sop = sub_lambda(sub_mu(subs(u, y_vars, hx_sym), mu, mu_val), lambda, k1_lambda);
    [~, ~, k1_delta, slack_persample, slack_max] = ...
        solve_k1_controller_sop_slack_persample( ...
            ux_for_sop, k1, J_k1, cert_term_dict, p, r_deg, ...
            x_vars, y_vars, hx_sym, safe_set, target_set, x_samples_valid, ...
            lb, ub, ds, dv, k1_lambda, mu_val, slack_weight);

    % ---- write per-sample slack CSV -----------------------------------------
    sample_idx = (1:valid_count)';
    T_ps = table(sample_idx, x_samples_valid(1,:)', x_samples_valid(2,:)', ...
        x_samples_valid(3,:)', x_samples_valid(4,:)', x_samples_valid(5,:)', ...
        slack_persample(:,1), slack_persample(:,2), ...
        'VariableNames', {'sample_idx', 'x1', 'x2', 'x3', 'x4', 'x5', 'slack1', 'slack2'});
    ps_path = fullfile(data_dir, 'slack_persample_manipulator.csv');
    writetable(T_ps, ps_path);
    fprintf('  wrote %s  (%d rows)\n', ps_path, valid_count);

    % ---- meta CSV -----------------------------------------------------------
    T = table(mu_val, xi0, ub1, ub2, valid_count, slack_weight, k1_lambda, k1_delta, ds, dv, ...
        slack_max(1), slack_max(2), ...
        'VariableNames', {'mu', 'xi0', 'ub1_req', 'ub2_req', 'n_valid_samples', 'slack_weight', ...
                          'lambda', 'delta', 'ds', 'dv', 'slack1_max', 'slack2_max'});
    meta_path = fullfile(data_dir, 'slack_meta_manipulator.csv');
    writetable(T, meta_path);
    fprintf('  wrote %s\n', meta_path);

    fprintf('solve_manipulator_slack: done. lambda=%.4g delta=%.4g n_valid=%d\n', ...
        k1_lambda, k1_delta, valid_count);
    ubv = [ub1, ub2];
    for i = 1:2
        s = slack_persample(:, i);
        fprintf('  channel %d slack: max=%.4g p50=%.4g p90=%.4g p99=%.4g nz=%.1f%%  (achievable ub+max=%.4g)\n', ...
            i, max(s), median(s), prctile(s, 90), prctile(s, 99), 100*mean(s > 1e-6), ubv(i) + max(s));
    end
end

% =============================================================================
% manifold sampler: draw (x1..x4) in the box, set x5 = x1 + x2, keep reach-avoid
% (safe >= 0, outside target, certificate >= 0). Returns 5 x n_valid.
function Xv = manifold_sample(n_valid, x_vars, y_vars, hx, safe_set, target_set, cert_x, box_lo, box_hi)
    safe_x = subs(safe_set, y_vars, hx);
    target_x = subs(target_set, y_vars, hx);
    % substitute the manifold x5 = x1 + x2 so everything is a function of x1..x4
    x1 = x_vars(1); x2 = x_vars(2); x5 = x_vars(5);
    safe_x = subs(safe_x, x5, x1 + x2);
    target_x = subs(target_x, x5, x1 + x2);
    cert_x = subs(cert_x, x5, x1 + x2);
    v4 = x_vars(1:4);
    safe_f = matlabFunction(safe_x, 'Vars', v4);
    targ_f = matlabFunction(target_x, 'Vars', v4);
    cert_f = matlabFunction(cert_x, 'Vars', v4);

    Xv = zeros(5, 0);
    batch = max(5000, 20 * n_valid);
    while size(Xv, 2) < n_valid
        S = box_lo(:)' + (box_hi(:)' - box_lo(:)') .* rand(batch, 4);  % batch x 4
        a = num2cell(S, 1);
        sv = safe_f(a{:}); tv = targ_f(a{:}); cv = cert_f(a{:});
        if isscalar(sv); sv = sv * ones(batch, 1); end
        if isscalar(tv); tv = tv * ones(batch, 1); end
        if isscalar(cv); cv = cv * ones(batch, 1); end
        keep = (sv(:) >= 0) & (tv(:) >= 0) & (cv(:) >= 0);
        Sk = S(keep, :);
        x5k = Sk(:, 1) + Sk(:, 2);
        Xv = [Xv, [Sk, x5k]'];   %#ok<AGROW>  5 x m
    end
    Xv = Xv(:, 1:n_valid);
end

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

function v = read_env_num(name)
    % All parameters come from the calling notebook; no defaults are set here.
    raw = getenv(name);
    if isempty(raw)
        error('solve_manipulator_slack:missingEnv', ...
            '%s must be set by the calling notebook (the .m sets no parameter values)', name);
    end
    v = str2double(raw);
    if isnan(v); error('solve_manipulator_slack:badEnv', '%s=%s is not a number', name, raw); end
end
