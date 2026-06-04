% §8 driver: HARD-bound SOP on the Dubins demo configuration.
% Uses sample_n_valid (same as build_dubins_demo) to control the exact number of
% per-sample constraints added to the SOS program.

function test_hardbound_dubins()
    mu_val  = read_env_num('MU',      0.1);
    xi0     = read_env_num('XI0',     10);
    ub_req  = read_env_num('UB',      20);
    n_valid = read_env_num('N_VALID', 250);

    matlab_local = fileparts(mfilename('fullpath'));
    demo_dir     = fileparts(matlab_local);
    data_dir     = fullfile(demo_dir, 'data');
    if ~exist(data_dir, 'dir'); mkdir(data_dir); end

    fprintf('test_hardbound_dubins: mu=%g xi0=%g ub=%g n_valid=%d\n', ...
        mu_val, xi0, ub_req, n_valid);

    % --- Dubins system --------------------------------------------------------
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

    bound_min = [-2; -2; 2*pi/3; 0.1];
    bound_max = [2; 2; 4*pi/3; 1.0];
    cert_vanilla = sub_lambda(sub_mu(subs(subs(certificate, k1, k1_y), y_vars, hx_sym), mu, mu_val), lambda, k1_lambda);
    rng(42);
    x_samples_valid = sample_n_valid(n_valid, x_vars, y_vars, hx_sym, ...
        safe_set, target_set, cert_vanilla, bound_min, bound_max);
    valid_count = size(x_samples_valid, 2);

    lb = [-ub_req; -ub_req]; ub = [ub_req; ub_req];
    ux_for_sop = sub_lambda(sub_mu(subs(u, y_vars, hx_sym), mu, mu_val), lambda, k1_lambda);
    [k1_con, J_k1_con, k1_delta, feasible, info] = solve_k1_controller_sop_hard( ...
        ux_for_sop, k1, J_k1, cert_term_dict, p, r_deg, ...
        x_vars, y_vars, hx_sym, safe_set, target_set, x_samples_valid, ...
        lb, ub, ds, dv, k1_lambda, mu_val);

    suffix = sprintf('%s_%s_%s', fmt_num(mu_val), fmt_num(xi0), fmt_num(ub_req));

    if feasible
        ux_con = sub_lambda(sub_mu(subs(subs(subs(u, k1, k1_con), J_k1, J_k1_con), y_vars, hx_sym), mu, mu_val), lambda, k1_lambda);
        cert_con = sub_lambda(sub_mu(subs(subs(certificate, k1, k1_con), y_vars, hx_sym), mu, mu_val), lambda, k1_lambda);
        params = struct( ...
            'controller_type', 'constrained_hardbound', ...
            'mu_val', mu_val, 'xi0', xi0, 'ub_req', ub_req, ...
            'lambda', k1_lambda, 'delta', k1_delta, ...
            'pinf', info.pinf, 'dinf', info.dinf, 'feasratio', info.feasratio, ...
            'n_valid_samples', valid_count, 'ds', ds, 'dv', dv);
        export_to_python(ux_con, cert_con, k1_con, params, ...
            fullfile(data_dir, sprintf('controller_hard_%s.py', suffix)));
        fprintf('  wrote controller_hard_%s.py\n', suffix);
    else
        fprintf('  HARD SOP infeasible — not exporting a controller.\n');
    end

    T = table(mu_val, xi0, ub_req, valid_count, double(feasible), k1_lambda, ...
        info.pinf, info.dinf, info.feasratio, ...
        'VariableNames', {'mu', 'xi0', 'ub_req', 'n_valid_samples', 'feasible', 'lambda', ...
                          'pinf', 'dinf', 'feasratio'});
    meta_path = fullfile(data_dir, sprintf('hardbound_%s.csv', suffix));
    writetable(T, meta_path);
    fprintf('  wrote hardbound_%s.csv (feasible=%d)\n', suffix, double(feasible));
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
            error('test_hardbound_dubins:badEnv', '%s=%s is not a number', name, raw);
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
