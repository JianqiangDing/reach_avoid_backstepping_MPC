% Subset-experiment driver: load a *user-supplied* subset of SOP samples from a CSV
% and run the per-sample slack SOP with possibly per-channel bounds. Used by the
% notebook to test: "if I keep only the well-behaved samples (slack <= p80) and
% relax the bound to cover them (ub_demo = ub_request + p80), what do con and unc
% look like under the same (mu, xi0, lambda)?"
%
% Env vars (all required unless noted):
%   MU              backstepping coefficient (default 0.1)
%   XI0             SOS floor on lambda  (default 10)
%   UB1, UB2        per-channel bounds   (default 5, 5)
%   SAMPLES_CSV     absolute path to a CSV with columns x1, x2, th, v
%                   (the samples to use as SOP per-sample constraints)
%   TAG             optional filename tag (default 'subset')
%
% Outputs in scripts/dubins_demo/data/ with filenames suffixed by TAG:
%   controller_con_<TAG>_<mu>_<xi0>_<ub1>x<ub2>.py
%   controller_unc_<TAG>_<mu>_<xi0>_<ub1>x<ub2>.py
%   meta_<TAG>_<mu>_<xi0>_<ub1>x<ub2>.csv
%   slack_persample_<TAG>_<mu>_<xi0>_<ub1>x<ub2>.csv
function build_dubins_demo_subset()
    mu_val       = read_env_num('MU',           0.1);
    xi0          = read_env_num('XI0',          10);
    ub1          = read_env_num('UB1',          5);
    ub2          = read_env_num('UB2',          5);
    slack_weight = read_env_num('SLACK_WEIGHT', 1.0);
    tag          = getenv('TAG');
    if isempty(tag); tag = 'subset'; end
    samples_csv = getenv('SAMPLES_CSV');
    if isempty(samples_csv) || ~exist(samples_csv, 'file')
        error('build_dubins_demo_subset:noSamples', ...
            'SAMPLES_CSV env var must point to an existing CSV; got "%s"', samples_csv);
    end

    matlab_local = fileparts(mfilename('fullpath'));
    demo_dir     = fileparts(matlab_local);
    data_dir     = fullfile(demo_dir, 'data');
    if ~exist(data_dir, 'dir'); mkdir(data_dir); end

    fprintf('build_dubins_demo_subset: tag=%s mu=%g xi0=%g ub=[%g %g]\n', ...
        tag, mu_val, xi0, ub1, ub2);
    fprintf('  samples_csv: %s\n', samples_csv);

    % --- load samples from CSV ----------------------------------------------
    T_in = readtable(samples_csv);
    required = {'x1', 'x2', 'th', 'v'};
    for k = 1:numel(required)
        if ~ismember(required{k}, T_in.Properties.VariableNames)
            error('build_dubins_demo_subset:badCSV', 'CSV missing column %s', required{k});
        end
    end
    x_samples_valid = [T_in.x1, T_in.x2, T_in.th, T_in.v]';   % 4 × n_samples
    valid_count = size(x_samples_valid, 2);
    fprintf('  loaded %d samples\n', valid_count);

    % --- Dubins system ------------------------------------------------------
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

    % --- per-sample slack SOP on the loaded samples + per-channel bounds ----
    lb = [-ub1; -ub2]; ub = [ub1; ub2];
    ux_for_sop = sub_lambda(sub_mu(subs(u, y_vars, hx_sym), mu, mu_val), lambda, k1_lambda);
    [k1_con, J_k1_con, k1_delta, slack_persample, slack_max] = ...
        solve_k1_controller_sop_slack_persample( ...
            ux_for_sop, k1, J_k1, cert_term_dict, p, r_deg, ...
            x_vars, y_vars, hx_sym, safe_set, target_set, x_samples_valid, ...
            lb, ub, ds, dv, k1_lambda, mu_val, slack_weight);

    % --- assemble controllers + certificates --------------------------------
    ux_con = sub_lambda(sub_mu(subs(subs(subs(u, k1, k1_con), J_k1, J_k1_con), y_vars, hx_sym), mu, mu_val), lambda, k1_lambda);
    cert_con = sub_lambda(sub_mu(subs(subs(certificate, k1, k1_con), y_vars, hx_sym), mu, mu_val), lambda, k1_lambda);
    ux_uncon = sub_lambda(sub_mu(subs(subs(subs(u, k1, k1_y), J_k1, J_k1_y), y_vars, hx_sym), mu, mu_val), lambda, k1_lambda);
    cert_uncon = sub_lambda(sub_mu(subs(subs(certificate, k1, k1_y), y_vars, hx_sym), mu, mu_val), lambda, k1_lambda);

    suffix = sprintf('%s_%s_%s_%sx%s', tag, fmt_num(mu_val), fmt_num(xi0), fmt_num(ub1), fmt_num(ub2));

    params_con = struct( ...
        'controller_type',  'constrained_slack_persample_subset', ...
        'mu_val', mu_val, 'xi0', xi0, 'ub1', ub1, 'ub2', ub2, 'tag', tag, ...
        'lambda', k1_lambda, 'delta', k1_delta, ...
        'slack_weight', slack_weight, ...
        'slack_max1', slack_max(1), 'slack_max2', slack_max(2), ...
        'n_valid_samples', valid_count, 'ds', ds, 'dv', dv);
    export_to_python(ux_con, cert_con, k1_con, params_con, ...
        fullfile(data_dir, sprintf('controller_con_%s.py', suffix)));
    fprintf('  wrote controller_con_%s.py\n', suffix);

    params_unc = struct( ...
        'controller_type', 'unconstrained_vanilla_same_xi0_lambda', ...
        'mu_val', mu_val, 'xi0', xi0, 'tag', tag, ...
        'lambda', k1_lambda, 'ds', ds, 'dv', dv);
    export_to_python(ux_uncon, cert_uncon, k1_y, params_unc, ...
        fullfile(data_dir, sprintf('controller_unc_%s.py', suffix)));
    fprintf('  wrote controller_unc_%s.py\n', suffix);

    % per-sample slack table for the subset run
    sample_idx = (1:valid_count)';
    T_ps = table(sample_idx, x_samples_valid(1,:)', x_samples_valid(2,:)', ...
        x_samples_valid(3,:)', x_samples_valid(4,:)', ...
        slack_persample(:,1), slack_persample(:,2), ...
        'VariableNames', {'sample_idx', 'x1', 'x2', 'th', 'v', 'slack1', 'slack2'});
    writetable(T_ps, fullfile(data_dir, sprintf('slack_persample_%s.csv', suffix)));

    % scalar summary
    T_meta = table(mu_val, xi0, ub1, ub2, valid_count, slack_weight, k1_lambda, k1_delta, ...
        slack_max(1), slack_max(2), ...
        'VariableNames', {'mu', 'xi0', 'ub1', 'ub2', 'n_valid_samples', 'slack_weight', 'lambda', 'delta', ...
                          'slack1_max', 'slack2_max'});
    writetable(T_meta, fullfile(data_dir, sprintf('meta_%s.csv', suffix)));

    fprintf('build_dubins_demo_subset: done. lambda=%.4g delta=%.4g n_valid=%d slack=[%.4g %.4g]\n', ...
        k1_lambda, k1_delta, valid_count, slack_max(1), slack_max(2));
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
            error('build_dubins_demo_subset:badEnv', '%s=%s is not a number', name, raw);
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
