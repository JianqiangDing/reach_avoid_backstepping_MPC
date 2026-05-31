% Phase 1.2 -- Dubins per-sample slack diagnostic (self-contained).
%
% Requests a tight input bound ub and solves the per-sample slack SOP: each
% (sample j, channel i) gets its own slack s_{j,i} >= 0 with
%   |u_i(x_j)| <= ub_i + s_{j,i},   minimize  delta + slack_weight * sum s_{j,i},
% while the reach-avoid (descent + certificate LMI) stays HARD. The optimal
% per-sample slacks form a distribution we inspect (in dubins_slack_diag.ipynb)
% to pick a sensible effective bound ub_eff and the kept sample subset.
%
% Self-contained: addpath only revision/matlab_frozen; the synthesis region
% X_S_eff is read from revision/data/phase1_1_outputs_dubins.json (Phase 1.1);
% outputs are written to revision/data/. No controller is exported (diagnostic only).
%
% ALL parameters are supplied by the calling notebook via environment variables
% (this .m sets no parameter values): MU, XI0, UB, N_VALID, SLACK_WEIGHT, DS, DV.

function solve_dubins_slack()
    here = fileparts(mfilename('fullpath'));        % revision/slack_diagnostic
    revision_dir = fileparts(here);                 % revision
    addpath(fullfile(revision_dir, 'matlab_frozen'));
    data_dir = fullfile(revision_dir, 'data');
    if ~exist(data_dir, 'dir'); mkdir(data_dir); end

    mu_val       = read_env_num('MU');
    xi0          = read_env_num('XI0');
    ub_req       = read_env_num('UB');
    n_valid      = read_env_num('N_VALID');
    slack_weight = read_env_num('SLACK_WEIGHT');
    ds           = read_env_num('DS');
    dv           = read_env_num('DV');

    % ---- synthesis region X_S_eff from Phase 1.1 ----------------------------
    j = jsondecode(fileread(fullfile(data_dir, 'phase1_1_outputs_dubins.json')));
    bound_min = j.X_S_eff_def.lower(:);   % [x1; x2; theta; v]
    bound_max = j.X_S_eff_def.upper(:);
    fprintf('solve_dubins_slack: mu=%g xi0=%g ub=%g n_valid=%d slack_weight=%g\n', ...
        mu_val, xi0, ub_req, n_valid, slack_weight);
    fprintf('  X_S_eff box (from phase1_1_outputs_dubins.json):\n');
    fprintf('    lower = [%s]\n    upper = [%s]\n', ...
        strtrim(sprintf('%.4g ', bound_min)), strtrim(sprintf('%.4g ', bound_max)));

    % ---- Dubins system (in-script; self-contained) --------------------------
    syms x1 x2 th v y1 y2;
    x_vars = [x1; x2; th; v];
    y_vars = [y1; y2];
    fx_sym = [v * cos(th); v * sin(th); 0; 0];
    gx_sym = [0, 0; 0, 0; 1, 0; 0, 1];
    hx_sym = [x1; x2];
    h_raw  = -(y1^4 + y2^4 - 16) * (y1^4 + y2^4 - 4);
    target_set = (y2 - 0)^2 + ((y1 + 1.7) / 0.5)^2 - 0.4;
    safe_set = 1e-3 * (-target_set + 300) * h_raw;

    [u, k1, J_k1, mu, lambda, certificate, cert_term_dict, ~, ~, ~, p, r_deg] = ...
        reach_avoid_controller(fx_sym, gx_sym, hx_sym, x_vars, y_vars, safe_set);

    [k1_y, k1_lambda, ~] = solve_vanilla_k1_controller_xi(y_vars, safe_set, target_set, dv, ds, xi0);

    % ---- sample n_valid reach-avoid states over X_S_eff ---------------------
    cert_vanilla = sub_lambda(sub_mu(subs(subs(certificate, k1, k1_y), y_vars, hx_sym), mu, mu_val), lambda, k1_lambda);
    rng(42);
    x_samples_valid = sample_n_valid(n_valid, x_vars, y_vars, hx_sym, ...
        safe_set, target_set, cert_vanilla, bound_min, bound_max);
    valid_count = size(x_samples_valid, 2);

    % ---- per-sample slack SOP -----------------------------------------------
    lb = [-ub_req; -ub_req]; ub = [ub_req; ub_req];
    ux_for_sop = sub_lambda(sub_mu(subs(u, y_vars, hx_sym), mu, mu_val), lambda, k1_lambda);
    [~, ~, k1_delta, slack_persample, slack_max] = ...
        solve_k1_controller_sop_slack_persample( ...
            ux_for_sop, k1, J_k1, cert_term_dict, p, r_deg, ...
            x_vars, y_vars, hx_sym, safe_set, target_set, x_samples_valid, ...
            lb, ub, ds, dv, k1_lambda, mu_val, slack_weight);

    % ---- write per-sample slack CSV (slack distribution + sample coords) ----
    sample_idx = (1:valid_count)';
    T_ps = table(sample_idx, x_samples_valid(1,:)', x_samples_valid(2,:)', ...
        x_samples_valid(3,:)', x_samples_valid(4,:)', ...
        slack_persample(:,1), slack_persample(:,2), ...
        'VariableNames', {'sample_idx', 'x1', 'x2', 'th', 'v', 'slack1', 'slack2'});
    ps_path = fullfile(data_dir, 'slack_persample_dubins.csv');
    writetable(T_ps, ps_path);
    fprintf('  wrote %s  (%d rows)\n', ps_path, valid_count);

    % ---- meta CSV with the probe parameters + summary -----------------------
    T = table(mu_val, xi0, ub_req, valid_count, slack_weight, ds, dv, k1_lambda, k1_delta, ...
        slack_max(1), slack_max(2), ...
        'VariableNames', {'mu', 'xi0', 'ub_req', 'n_valid_samples', 'slack_weight', 'ds', 'dv', ...
                          'lambda', 'delta', 'slack1_max', 'slack2_max'});
    meta_path = fullfile(data_dir, 'slack_meta_dubins.csv');
    writetable(T, meta_path);
    fprintf('  wrote %s\n', meta_path);

    fprintf('solve_dubins_slack: done. lambda=%.4g delta=%.4g n_valid=%d\n', ...
        k1_lambda, k1_delta, valid_count);
    for i = 1:2
        s = slack_persample(:, i);
        fprintf('  channel %d slack: max=%.4g p50=%.4g p90=%.4g p99=%.4g nz=%.1f%%  (achievable bound ub+max=%.4g)\n', ...
            i, max(s), median(s), prctile(s, 90), prctile(s, 99), 100*mean(s > 1e-6), ub_req + max(s));
    end
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

function v = read_env_num(name)
    % All parameters come from the calling notebook; no defaults are set here.
    raw = getenv(name);
    if isempty(raw)
        error('solve_dubins_slack:missingEnv', ...
            '%s must be set by the calling notebook (the .m sets no parameter values)', name);
    end
    v = str2double(raw);
    if isnan(v); error('solve_dubins_slack:badEnv', '%s=%s is not a number', name, raw); end
end
