% Synthesize the two Dubins controllers at the same (mu, xi0):
%   - constrained: per-sample slack SOP
%       (each SOP sample j gets its own slack s_{j,i} per channel; the
%        objective minimizes delta + sum_{j,i} s_{j,i}, giving a *distribution*
%        of per-sample violations that we can inspect to pick a sensible ub)
%   - unconstrained: vanilla k1 at the SAME lambda
%
% Outputs (in scripts/dubins_demo/data/):
%   controller_con_<mu>_<xi0>_<ub>.py     constrained controller + cert
%   controller_unc_<mu>_<xi0>_<ub>.py     unconstrained controller + cert
%   meta_<mu>_<xi0>_<ub>.csv              scalar metadata (lambda, delta, n_valid, slack stats)
%   slack_persample_<mu>_<xi0>_<ub>.csv   per-sample slacks + sample coords
%                                         (columns: sample_idx, x1, x2, th, v, slack1, slack2)

function build_dubins_demo()
    mu_val       = read_env_num('MU',           0.1);
    xi0          = read_env_num('XI0',          10);
    ub_req       = read_env_num('UB',           20);
    n_valid      = read_env_num('N_VALID',      250);
    slack_weight = read_env_num('SLACK_WEIGHT', 1.0);

    matlab_local = fileparts(mfilename('fullpath'));
    demo_dir     = fileparts(matlab_local);
    data_dir     = fullfile(demo_dir, 'data');
    if ~exist(data_dir, 'dir'); mkdir(data_dir); end

    fprintf('build_dubins_demo: mu=%g xi0=%g ub=%g n_valid=%d\n', ...
        mu_val, xi0, ub_req, n_valid);

    % ---- Dubins system ------------------------------------------------------
    syms x1 x2 th v y1 y2;
    x_vars = [x1; x2; th; v];
    y_vars = [y1; y2];
    fx_sym = [v*cos(th); v*sin(th); 0; 0];
    gx_sym = [0, 0; 0, 0; 1, 0; 0, 1];
    hx_sym = [x1; x2];
    h_raw  = -(y1^4 + y2^4 - 16) * (y1^4 + y2^4 - 4);
    target_set = (y2 - 0)^2 + ((y1 + 1.7) / 0.5)^2 - 0.4;
    alpha  = 1e-3 * (-target_set + 300);
    safe_set = alpha * h_raw;

    [u, k1, J_k1, mu, lambda, certificate, cert_term_dict, ~, ~, ~, p, r_deg] = ...
        reach_avoid_controller(fx_sym, gx_sym, hx_sym, x_vars, y_vars, safe_set);

    ds = 4; dv = 4;
    [k1_y, k1_lambda, ~] = solve_vanilla_k1_controller_xi(y_vars, safe_set, target_set, dv, ds, xi0);
    J_k1_y = jacobian(k1_y, y_vars);

    bound_min = [-2; -2; 2*pi/3; 0.1];
    bound_max = [2; 2; 4*pi/3; 1.0];
    cert_vanilla = sub_lambda(sub_mu(subs(subs(certificate, k1, k1_y), y_vars, hx_sym), mu, mu_val), lambda, k1_lambda);
    rng(42);
    x_samples_valid = sample_n_valid(n_valid, x_vars, y_vars, hx_sym, ...
        safe_set, target_set, cert_vanilla, bound_min, bound_max);
    valid_count = size(x_samples_valid, 2);

    % ---- per-sample slack SOP -----------------------------------------------
    lb = [-ub_req; -ub_req]; ub = [ub_req; ub_req];
    ux_for_sop = sub_lambda(sub_mu(subs(u, y_vars, hx_sym), mu, mu_val), lambda, k1_lambda);
    [k1_con, J_k1_con, k1_delta, slack_persample, slack_max] = ...
        solve_k1_controller_sop_slack_persample( ...
            ux_for_sop, k1, J_k1, cert_term_dict, p, r_deg, ...
            x_vars, y_vars, hx_sym, safe_set, target_set, x_samples_valid, ...
            lb, ub, ds, dv, k1_lambda, mu_val, slack_weight);

    % ---- assemble controllers + certificates in x-space ---------------------
    ux_con = sub_lambda(sub_mu(subs(subs(subs(u, k1, k1_con), J_k1, J_k1_con), y_vars, hx_sym), mu, mu_val), lambda, k1_lambda);
    cert_con = sub_lambda(sub_mu(subs(subs(certificate, k1, k1_con), y_vars, hx_sym), mu, mu_val), lambda, k1_lambda);
    ux_uncon = sub_lambda(sub_mu(subs(subs(subs(u, k1, k1_y), J_k1, J_k1_y), y_vars, hx_sym), mu, mu_val), lambda, k1_lambda);
    cert_uncon = sub_lambda(sub_mu(subs(subs(certificate, k1, k1_y), y_vars, hx_sym), mu, mu_val), lambda, k1_lambda);

    suffix = sprintf('%s_%s_%s', fmt_num(mu_val), fmt_num(xi0), fmt_num(ub_req));

    params_con = struct( ...
        'controller_type',  'constrained_slack_persample', ...
        'mu_val', mu_val, 'xi0', xi0, 'ub_req', ub_req, ...
        'lambda', k1_lambda, 'delta', k1_delta, ...
        'slack_weight', slack_weight, ...
        'slack_max1', slack_max(1), 'slack_max2', slack_max(2), ...
        'achievable_max1', ub_req + slack_max(1), ...
        'achievable_max2', ub_req + slack_max(2), ...
        'n_valid_samples', valid_count, 'ds', ds, 'dv', dv);
    export_to_python(ux_con, cert_con, k1_con, params_con, ...
        fullfile(data_dir, sprintf('controller_con_%s.py', suffix)));
    fprintf('  wrote controller_con_%s.py\n', suffix);

    params_unc = struct( ...
        'controller_type', 'unconstrained_vanilla', ...
        'mu_val', mu_val, 'xi0', xi0, 'ub_req', ub_req, ...
        'lambda', k1_lambda, ...
        'n_valid_samples', valid_count, 'ds', ds, 'dv', dv);
    export_to_python(ux_uncon, cert_uncon, k1_y, params_unc, ...
        fullfile(data_dir, sprintf('controller_unc_%s.py', suffix)));
    fprintf('  wrote controller_unc_%s.py\n', suffix);

    % ---- per-sample slack CSV (slack distribution + sample coords) ---------
    sample_idx = (1:valid_count)';
    T_ps = table(sample_idx, x_samples_valid(1,:)', x_samples_valid(2,:)', ...
        x_samples_valid(3,:)', x_samples_valid(4,:)', ...
        slack_persample(:,1), slack_persample(:,2), ...
        'VariableNames', {'sample_idx', 'x1', 'x2', 'th', 'v', 'slack1', 'slack2'});
    ps_path = fullfile(data_dir, sprintf('slack_persample_%s.csv', suffix));
    writetable(T_ps, ps_path);
    fprintf('  wrote slack_persample_%s.csv  (%d rows)\n', suffix, valid_count);

    % ---- meta CSV with summary stats ----------------------------------------
    p50_1 = median(slack_persample(:,1)); p90_1 = prctile(slack_persample(:,1), 90);
    p99_1 = prctile(slack_persample(:,1), 99);
    p50_2 = median(slack_persample(:,2)); p90_2 = prctile(slack_persample(:,2), 90);
    p99_2 = prctile(slack_persample(:,2), 99);
    nz_frac1 = mean(slack_persample(:,1) > 1e-6);
    nz_frac2 = mean(slack_persample(:,2) > 1e-6);
    T = table(mu_val, xi0, ub_req, valid_count, slack_weight, k1_lambda, k1_delta, ...
        slack_max(1), p50_1, p90_1, p99_1, nz_frac1, ...
        slack_max(2), p50_2, p90_2, p99_2, nz_frac2, ...
        'VariableNames', {'mu', 'xi0', 'ub_req', 'n_valid_samples', 'slack_weight', 'lambda', 'delta', ...
                          'slack1_max', 'slack1_p50', 'slack1_p90', 'slack1_p99', 'slack1_nz_frac', ...
                          'slack2_max', 'slack2_p50', 'slack2_p90', 'slack2_p99', 'slack2_nz_frac'});
    meta_path = fullfile(data_dir, sprintf('meta_%s.csv', suffix));
    writetable(T, meta_path);
    fprintf('  wrote meta_%s.csv\n', suffix);

    fprintf('build_dubins_demo: done. lambda=%.4g delta=%.4g n_valid=%d\n', ...
        k1_lambda, k1_delta, valid_count);
    fprintf('  channel 1 slack: max=%.4g  p50=%.4g  p90=%.4g  p99=%.4g  nz=%.1f%%\n', ...
        slack_max(1), p50_1, p90_1, p99_1, 100*nz_frac1);
    fprintf('  channel 2 slack: max=%.4g  p50=%.4g  p90=%.4g  p99=%.4g  nz=%.1f%%\n', ...
        slack_max(2), p50_2, p90_2, p99_2, 100*nz_frac2);
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

function v = read_env_num(name, default_value)
    raw = getenv(name);
    if isempty(raw)
        v = default_value;
    else
        v = str2double(raw);
        if isnan(v)
            error('build_dubins_demo:badEnv', '%s=%s is not a number', name, raw);
        end
    end
end

function s = fmt_num(v)
    if v == fix(v) && abs(v) < 1e6
        s = sprintf('%d', round(v));
    else
        s = strrep(sprintf('%g', v), '.', 'p');
    end
end
