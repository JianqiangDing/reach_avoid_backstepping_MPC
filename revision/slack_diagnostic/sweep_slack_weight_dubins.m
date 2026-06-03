% Phase 1.2 diagnostic -- slack_weight sweep for the dubins per-sample slack SOP.
%
% Motivation: the per-sample slack SOP returns a descent slack `delta` and a fixed
% backstepping scale `lambda`. The certified reach-avoid funnel of the reduced
% single-integrator system (ydot = k1) is {y : psi(y) := safe(y) >= delta/lambda}
% -- because along the flow d/dt safe = grad(safe).k1 >= lambda*safe - delta, which
% is positive (barrier grows, trajectory driven to the target while staying safe)
% only where safe > delta/lambda. A large delta/lambda shrinks -- and can EMPTY --
% that funnel. This script sweeps slack_weight (default down to 1e-6) and records
% delta/lambda for each, plus the VANILLA (no per-sample constraints) reference.
%
% One system build + ONE set of valid samples (seed-reproducible) is reused for all
% weights, so the only varying knob is slack_weight. Settings (mu, xi0, ub_req,
% n_valid, ds, dv, seed) are read from the saved Phase 1.2 probe -> self-contained.
%
% Output: revision/data/slack_weight_sweep_dubins.csv
%   slack_weight, lambda, delta, delta_over_lambda, slack1_max, slack2_max
% (vanilla row has slack_weight = NaN). Companion: slack_weight_reachset_diag.py
% counts how much of X_S\X_T lies in {safe >= delta/lambda} for each row.
%
% Self-contained: addpath only revision/matlab_frozen; reads only revision/data/.

function sweep_slack_weight_dubins()
    here = fileparts(mfilename('fullpath'));        % revision/slack_diagnostic
    revision_dir = fileparts(here);                 % revision
    addpath(fullfile(revision_dir, 'matlab_frozen'));
    data_dir = fullfile(revision_dir, 'data');

    % ---- Phase 1.2 settings (read from the saved probe; self-contained) -----
    p12 = jsondecode(fileread(fullfile(data_dir, 'phase1_2_outputs_dubins.json')));
    mu_val  = p12.probe.mu;   xi0 = p12.probe.xi0;
    ds      = p12.probe.ds;   dv  = p12.probe.dv;
    ub_req  = p12.probe.ub_req;
    n_valid = p12.probe.n_samples;
    seed    = p12.probe.seed;

    % ---- slack_weight list (env SLACK_WEIGHTS, comma-separated; or default) --
    raw = getenv('SLACK_WEIGHTS');
    if isempty(raw)
        weights = [1, 0.1, 1e-2, 1e-3, 1e-4, 1e-5, 1e-6];
    else
        weights = str2num(raw); %#ok<ST2NM>
    end
    fprintf('sweep_slack_weight_dubins: mu=%g xi0=%g ub=%g n_valid=%d ds=%d dv=%d seed=%d\n', ...
        mu_val, xi0, ub_req, n_valid, ds, dv, seed);
    fprintf('  slack_weights = [%s]\n', strtrim(sprintf('%g ', weights)));

    % ---- X_S_eff box (Phase 1.1) -------------------------------------------
    j = jsondecode(fileread(fullfile(data_dir, 'phase1_1_outputs_dubins.json')));
    bound_min = j.X_S_eff_def.lower(:);
    bound_max = j.X_S_eff_def.upper(:);

    % ---- dubins system (frozen example def) --------------------------------
    syms x1 x2 th v y1 y2;
    x_vars = [x1; x2; th; v];
    y_vars = [y1; y2];
    fx_sym = [v * cos(th); v * sin(th); 0; 0];
    gx_sym = [0, 0; 0, 0; 1, 0; 0, 1];
    hx_sym = [x1; x2];
    h_raw = -(y1^4 + y2^4 - 16) * (y1^4 + y2^4 - 4);
    target_set = (y2 - 0)^2 + ((y1 + 1.7) / 0.5)^2 - 0.4;
    safe_set = 1e-3 * (-target_set + 300) * h_raw;

    [u, k1, J_k1, mu, lambda, certificate, cert_term_dict, ~, ~, ~, p, r_deg] = ...
        reach_avoid_controller(fx_sym, gx_sym, hx_sym, x_vars, y_vars, safe_set);

    % ---- vanilla k1 (NO per-sample constraints): lambda + reference delta ---
    [k1_y, k1_lambda, k1_delta_vanilla] = ...
        solve_vanilla_k1_controller_xi(y_vars, safe_set, target_set, dv, ds, xi0);
    fprintf('[vanilla] lambda=%.6g delta=%.6g  (delta/lambda=%.6g)\n', ...
        k1_lambda, k1_delta_vanilla, k1_delta_vanilla / k1_lambda);

    % ---- ONE set of valid samples reused for every weight (seed-reproducible)
    cert_vanilla = sub_lambda(sub_mu(subs(subs(certificate, k1, k1_y), y_vars, hx_sym), mu, mu_val), lambda, k1_lambda);
    rng(seed);
    x_samples_valid = sample_n_valid(n_valid, x_vars, y_vars, hx_sym, ...
        safe_set, target_set, cert_vanilla, bound_min, bound_max);
    valid_count = size(x_samples_valid, 2);
    fprintf('  %d valid reach-avoid samples (reused for all weights)\n', valid_count);

    ux_for_sop = sub_lambda(sub_mu(subs(u, y_vars, hx_sym), mu, mu_val), lambda, k1_lambda);
    lb = [-ub_req; -ub_req]; ub = [ub_req; ub_req];

    % ---- sweep --------------------------------------------------------------
    wcol = NaN; lamcol = k1_lambda; delcol = k1_delta_vanilla; s1col = 0; s2col = 0; % vanilla row
    for w = weights
        t = tic;
        [~, ~, k1_delta, ~, slack_max] = solve_k1_controller_sop_slack_persample( ...
            ux_for_sop, k1, J_k1, cert_term_dict, p, r_deg, ...
            x_vars, y_vars, hx_sym, safe_set, target_set, x_samples_valid, ...
            lb, ub, ds, dv, k1_lambda, mu_val, w);
        fprintf('[w=%.4g] delta=%.6g  delta/lambda=%.6g  slack_max=[%.4g %.4g]  (%.1fs)\n', ...
            w, k1_delta, k1_delta / k1_lambda, slack_max(1), slack_max(2), toc(t));
        wcol(end+1,1)   = w;            %#ok<AGROW>
        lamcol(end+1,1) = k1_lambda;    %#ok<AGROW>
        delcol(end+1,1) = k1_delta;     %#ok<AGROW>
        s1col(end+1,1)  = slack_max(1); %#ok<AGROW>
        s2col(end+1,1)  = slack_max(2); %#ok<AGROW>
    end

    T = table(wcol, lamcol, delcol, delcol ./ lamcol, s1col, s2col, ...
        'VariableNames', {'slack_weight', 'lambda', 'delta', 'delta_over_lambda', ...
                          'slack1_max', 'slack2_max'});
    out = fullfile(data_dir, 'slack_weight_sweep_dubins.csv');
    writetable(T, out);
    fprintf('wrote %s (%d rows)\n', out, height(T));
end

% =============================================================================
function expr = sub_mu(expr, mu_cells, mu_val)
    for i = 1:numel(mu_cells)
        if ~isempty(mu_cells{i}); expr = subs(expr, mu_cells{i}, mu_val); end
    end
end

function expr = sub_lambda(expr, lambda_sym, lambda_val)
    expr = subs(expr, lambda_sym, lambda_val);
end
