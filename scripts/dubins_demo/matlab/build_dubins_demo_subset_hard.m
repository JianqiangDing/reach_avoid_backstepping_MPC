% Hard-bound (no per-sample slack) variant of build_dubins_demo_subset.m. Loads a
% user-supplied subset of SOP samples from a CSV, then solves
%
%     min  delta
%     s.t. y-space CBF descent ∈ Σ[y]
%          M_j(k_1) ⪰ 0                  ∀ j (certificate LMI per kept sample)
%          lb_i ≤ u_i(x_j) ≤ ub_i        ∀ j, i   (HARD bound, no slack)
%
% Used together with build_dubins_demo.m's slack distribution to test:
% "after picking ub_demo + the corresponding kept samples, can the SOP enforce
% the bound exactly at those samples?". A `feasible` flag is returned so the
% caller can detect infeasibility instead of crashing.
%
% Env vars:
%   MU             default 0.1
%   XI0            default 10
%   UB1, UB2       per-channel HARD bounds (default 5, 5)
%   SAMPLES_CSV    required: CSV with columns x1, x2, th, v
%   TAG            filename suffix tag (default 'subset_hard')

function build_dubins_demo_subset_hard()
    mu_val  = read_env_num('MU',  0.1);
    xi0     = read_env_num('XI0', 10);
    ub1     = read_env_num('UB1', 5);
    ub2     = read_env_num('UB2', 5);
    tag     = getenv('TAG');
    if isempty(tag); tag = 'subset_hard'; end
    samples_csv = getenv('SAMPLES_CSV');
    if isempty(samples_csv) || ~exist(samples_csv, 'file')
        error('build_dubins_demo_subset_hard:noSamples', ...
            'SAMPLES_CSV env var must point to an existing CSV; got "%s"', samples_csv);
    end

    matlab_local = fileparts(mfilename('fullpath'));
    demo_dir     = fileparts(matlab_local);
    data_dir     = fullfile(demo_dir, 'data');
    if ~exist(data_dir, 'dir'); mkdir(data_dir); end

    fprintf('build_dubins_demo_subset_hard: tag=%s mu=%g xi0=%g ub=[%g %g]\n', ...
        tag, mu_val, xi0, ub1, ub2);
    fprintf('  samples_csv: %s\n', samples_csv);

    % --- load samples from CSV ----------------------------------------------
    T_in = readtable(samples_csv);
    for k = {'x1', 'x2', 'th', 'v'}
        if ~ismember(k{1}, T_in.Properties.VariableNames)
            error('build_dubins_demo_subset_hard:badCSV', 'CSV missing column %s', k{1});
        end
    end
    x_samples_valid = [T_in.x1, T_in.x2, T_in.th, T_in.v]';
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

    % --- hard SOP on the loaded samples -------------------------------------
    lb = [-ub1; -ub2]; ub = [ub1; ub2];
    ux_for_sop = sub_lambda(sub_mu(subs(u, y_vars, hx_sym), mu, mu_val), lambda, k1_lambda);
    [k1_con, J_k1_con, k1_delta, feasible, info] = solve_k1_controller_sop_hard( ...
        ux_for_sop, k1, J_k1, cert_term_dict, p, r_deg, ...
        x_vars, y_vars, hx_sym, safe_set, target_set, x_samples_valid, ...
        lb, ub, ds, dv, k1_lambda, mu_val);

    suffix = sprintf('%s_%s_%s_%sx%s', tag, fmt_num(mu_val), fmt_num(xi0), fmt_num(ub1), fmt_num(ub2));

    if feasible
        ux_con = sub_lambda(sub_mu(subs(subs(subs(u, k1, k1_con), J_k1, J_k1_con), y_vars, hx_sym), mu, mu_val), lambda, k1_lambda);
        cert_con = sub_lambda(sub_mu(subs(subs(certificate, k1, k1_con), y_vars, hx_sym), mu, mu_val), lambda, k1_lambda);
        params_con = struct( ...
            'controller_type', 'constrained_hard_subset', ...
            'mu_val', mu_val, 'xi0', xi0, 'ub1', ub1, 'ub2', ub2, 'tag', tag, ...
            'lambda', k1_lambda, 'delta', k1_delta, ...
            'pinf', info.pinf, 'dinf', info.dinf, 'feasratio', info.feasratio, ...
            'n_valid_samples', valid_count, 'ds', ds, 'dv', dv);
        export_to_python(ux_con, cert_con, k1_con, params_con, ...
            fullfile(data_dir, sprintf('controller_con_%s.py', suffix)));
        fprintf('  wrote controller_con_%s.py\n', suffix);

        % unc controller in matching shape (always vanilla, but tagged with same suffix
        % for convenience when the notebook pairs them)
        ux_uncon = sub_lambda(sub_mu(subs(subs(subs(u, k1, k1_y), J_k1, J_k1_y), y_vars, hx_sym), mu, mu_val), lambda, k1_lambda);
        cert_uncon = sub_lambda(sub_mu(subs(subs(certificate, k1, k1_y), y_vars, hx_sym), mu, mu_val), lambda, k1_lambda);
        params_unc = struct( ...
            'controller_type', 'unconstrained_vanilla_same_xi0_lambda', ...
            'mu_val', mu_val, 'xi0', xi0, 'tag', tag, ...
            'lambda', k1_lambda, 'ds', ds, 'dv', dv);
        export_to_python(ux_uncon, cert_uncon, k1_y, params_unc, ...
            fullfile(data_dir, sprintf('controller_unc_%s.py', suffix)));
        fprintf('  wrote controller_unc_%s.py\n', suffix);
    else
        fprintf('  HARD SOP infeasible — not exporting controllers.\n');
    end

    T_meta = table(mu_val, xi0, ub1, ub2, valid_count, double(feasible), k1_lambda, k1_delta, ...
        info.pinf, info.dinf, info.feasratio, ...
        'VariableNames', {'mu', 'xi0', 'ub1', 'ub2', 'n_valid_samples', 'feasible', 'lambda', 'delta', ...
                          'pinf', 'dinf', 'feasratio'});
    meta_path = fullfile(data_dir, sprintf('meta_%s.csv', suffix));
    writetable(T_meta, meta_path);
    fprintf('  wrote meta_%s.csv (feasible=%d)\n', suffix, double(feasible));
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
            error('build_dubins_demo_subset_hard:badEnv', '%s=%s is not a number', name, raw);
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
