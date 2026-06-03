% Per-sample slack variant of solve_k1_controller_sop_slack.
%
% Each (sample j, channel i) gets its OWN non-negative slack s_{j,i}, so the
% optimal solution returns a *distribution* of slacks across the SOP samples:
%   |u_i(x_j)| <= ub_i + s_{j,i},  s_{j,i} >= 0,
%   minimize  delta + sum_{j,i} s_{j,i}.
% This lets the user inspect where (and how far) the controller exceeds the
% requested bound at each sample — and pick a reasonable ub by reading off the
% per-sample slack distribution (e.g. "95th percentile of slack is 5 → ub_eff ≈ ub+5").
%
% The reach-avoid (certificate non-negativity + descent) constraints remain HARD,
% identical to the single-slack variant.
%
% Returns:
%   k1_opt, J_k1_opt   solved controller polynomial + Jacobian
%   k1_delta           descent slack (scalar)
%   slack_persample    n_samples × m matrix of per-sample slacks
%   slack_max          1 × m max slack per channel (for backward-compat reporting)
function [k1_opt, J_k1_opt, k1_delta, slack_persample, slack_max] = ...
        solve_k1_controller_sop_slack_persample( ...
        ux, k1_sym, J_k1_sym, cert_term_dict, p, r_deg, ...
        x_vars, y_vars, hx, safe_set, target_set, x_samples_valid, ...
        lb, ub, ds, dv, k1_lambda, mu_val, slack_weight)
    % slack_weight (optional, default 1): coefficient on each per-sample slack
    %   in the objective  delta + slack_weight * sum_{j,i} s_{j,i}.
    %   Increasing slack_weight makes the solver prefer keeping s_{j,i} = 0
    %   over reducing delta, which helps when the unweighted problem returns
    %   a numerically-suboptimal solution (Mosek trades off slack for delta).
    if nargin < 19 || isempty(slack_weight)
        slack_weight = 1.0;
    end

    n_samples = size(x_samples_valid, 2);
    m = length(ux);

    % --- pvar variables for x, y --------------------------------------------
    x_vars_pvar = [];

    for i = 1:length(x_vars)
        x_vars_pvar = [x_vars_pvar; pvar(strcat(char(x_vars(i)), '_pvar'))]; %#ok<AGROW>
    end

    y_vars_pvar = [];

    for i = 1:p
        y_vars_pvar = [y_vars_pvar; pvar(strcat(char(y_vars(i)), '_pvar'))]; %#ok<AGROW>
    end

    % --- decision variables: delta + per-sample slacks ----------------------
    dpvar delta;
    monom_const = monomials(y_vars_pvar, 0); % degree-0 monomials → scalar dpvar

    prog = sosprogram(y_vars_pvar, delta);

    % k1 polynomial decision vars
    monom_k1 = monomials(y_vars_pvar, 0:dv);
    monom_sos = monomials(y_vars_pvar, 0:ds);
    k1_poly = [];

    for i = 1:p
        [prog, k1_i] = sospolyvar(prog, monom_k1);
        k1_poly = [k1_poly; k1_i]; %#ok<AGROW>
    end

    % per-sample slack scalars (created as degree-0 polynomial decision vars)
    slack_cells = cell(n_samples, m);

    for j = 1:n_samples

        for i = 1:m
            [prog, sji] = sospolyvar(prog, monom_const);
            slack_cells{j, i} = sji;
            prog = sosineq(prog, sji); % s_{j,i} >= 0
        end

    end

    % --- substitution pvars (k1, J_k1, trig dummies) ------------------------
    k1_pvar = [];

    for i = 1:p
        k1_pvar = [k1_pvar; pvar(strcat(char(k1_sym(i)), '_pvar'))]; %#ok<AGROW>
    end

    J_k1_poly = [];

    for i = 1:p
        J_row = [];

        for j = 1:p
            J_row = [J_row, diff(k1_poly(i), y_vars_pvar(j))]; %#ok<AGROW>
        end

        J_k1_poly = [J_k1_poly; J_row]; %#ok<AGROW>
    end

    J_k1_poly_flat = dpvar2poly(J_k1_poly(:));
    J_k1_sym_flat = J_k1_sym(:);
    J_k1_pvar_flat = [];

    for i = 1:length(J_k1_sym_flat)
        J_k1_pvar_flat = [J_k1_pvar_flat; pvar(strcat(char(J_k1_sym_flat(i)), '_pvar'))]; %#ok<AGROW>
    end

    trig_terms = detect_trigonometric_terms([ux; hx]);
    dummy_trig_vars = sym(zeros(length(trig_terms), 1));
    dummy_trig_vars_pvar = [];

    for i = 1:length(trig_terms)
        trig_term_str = char(trig_terms(i));
        trig_term_str = strrep(trig_term_str, '(', '_');
        trig_term_str = strrep(trig_term_str, ')', '');
        dummy_trig_vars(i) = str2sym(strcat('dummy_trig_var_', trig_term_str));
        dummy_trig_vars_pvar = [dummy_trig_vars_pvar; pvar(strcat('dummy_trig_var_', trig_term_str, '_pvar'))]; %#ok<AGROW>
    end

    hx_sub = hx;

    for i = 1:length(trig_terms)
        hx_sub = subs(hx_sub, trig_terms(i), dummy_trig_vars(i));
    end

    hx_pvar = polynomial(zeros(length(hx), 1));

    for i = 1:length(hx)
        hx_pvar(i) = sym2pvar(hx_sub(i), [x_vars; dummy_trig_vars], [x_vars_pvar; dummy_trig_vars_pvar]);
    end

    safe_set_pvar = sym2pvar(safe_set, y_vars, y_vars_pvar);
    target_set_pvar = sym2pvar(target_set, y_vars, y_vars_pvar);

    [prog, s1] = sospolyvar(prog, monom_sos);
    [prog, s2] = sospolyvar(prog, monom_sos);

    % --- reach-avoid descent (HARD, identical to the original) --------------
    grad_safe_set = [];

    for i = 1:p
        grad_safe_set = [grad_safe_set; diff(safe_set_pvar, y_vars_pvar(i))]; %#ok<AGROW>
    end

    Lie_safe_set = grad_safe_set(1) * k1_poly(1);

    for i = 2:p
        Lie_safe_set = Lie_safe_set + grad_safe_set(i) * k1_poly(i);
    end

    sos_expr = Lie_safe_set - k1_lambda * safe_set_pvar + delta - s1 * safe_set_pvar - s2 * target_set_pvar;
    prog = sosineq(prog, sos_expr);
    prog = sosineq(prog, delta);
    prog = sosineq(prog, s1);
    prog = sosineq(prog, s2);

    % --- per-sample SOFT control bounds ------------------------------------
    prog = add_control_limit_constraints_persample( ...
        prog, ux, x_samples_valid, lb, ub, slack_cells, x_vars, y_vars, ...
        k1_sym, J_k1_sym_flat, x_vars_pvar, y_vars_pvar, k1_pvar, J_k1_pvar_flat, ...
        dummy_trig_vars_pvar, k1_poly, J_k1_poly_flat, hx_pvar, trig_terms, dummy_trig_vars);

    % --- certificate non-negativity (HARD) ---------------------------------
    prog = add_certificate_nonnegativity_constraints( ...
        prog, cert_term_dict, p, r_deg, x_samples_valid, ...
        x_vars, y_vars, x_vars_pvar, y_vars_pvar, k1_sym, J_k1_sym_flat, dummy_trig_vars, ...
        k1_pvar, J_k1_pvar_flat, dummy_trig_vars_pvar, k1_poly, J_k1_poly_flat, hx_pvar, ...
        trig_terms, safe_set_pvar, mu_val);

    % --- objective: delta + slack_weight * sum of per-sample slacks ---------
    total = delta;

    for j = 1:n_samples

        for i = 1:m
            total = total + slack_weight * slack_cells{j, i};
        end

    end

    prog = sossetobj(prog, total);
    fprintf('per-sample slack SOP objective: delta + %.4g * sum(slack)\n', slack_weight);

    solver_opt.solver = 'mosek';
    prog = sossolve(prog, solver_opt);

    % --- extract solution ---------------------------------------------------
    k1_opt = [];

    for i = 1:p
        k1_opt_i = sosgetsol(prog, k1_poly(i));
        k1_opt = [k1_opt; poly2sym(k1_opt_i, y_vars_pvar, y_vars)]; %#ok<AGROW>
    end

    J_k1_opt = jacobian(k1_opt, y_vars);
    k1_delta = double(sosgetsol(prog, delta));

    slack_persample = zeros(n_samples, m);

    for j = 1:n_samples

        for i = 1:m
            slack_persample(j, i) = double(sosgetsol(prog, slack_cells{j, i}));
        end

    end

    slack_persample = max(slack_persample, 0); % numerical safety
    slack_max = max(slack_persample, [], 1);

    if isempty(k1_delta) || isnan(k1_delta) || any(isnan(slack_persample(:)))
        error('solve_k1_controller_sop_slack_persample:solverFailed', ...
        'SOS solve returned no valid solution (delta/slack empty or NaN).');
    end

    fprintf('per-sample slack solve:\n');
    fprintf('  delta = %.4g\n', k1_delta);

    for i = 1:m
        s = slack_persample(:, i);
        fprintf('  channel %d slack:  max=%.4g  mean=%.4g  p50=%.4g  p90=%.4g  p99=%.4g  nz_frac=%.1f%%\n', ...
            i, max(s), mean(s), median(s), prctile(s, 90), prctile(s, 99), 100 * mean(s > 1e-6));
    end

end

% =============================================================================
% certificate non-negativity (identical to the single-slack variant)
% =============================================================================
function prog = add_certificate_nonnegativity_constraints(prog, cert_term_dict, p, r_deg, x_samples_valid, ...
        x_vars, y_vars, x_vars_pvar, y_vars_pvar, k1_sym, J_k1_sym_flat, dummy_trig_vars, ...
        k1_pvar, J_k1_pvar_flat, dummy_trig_vars_pvar, k1_poly, J_k1_poly_flat, hx_pvar, ...
        trig_terms, safe_set_pvar, mu_val)

    term_sum = 0;

    for i = 1:p
        term_sum = term_sum + r_deg(i) - 1;
    end

    for sample_idx = 1:size(x_samples_valid, 2)
        x_sample = x_samples_valid(:, sample_idx);
        n_size = 1 + term_sum;
        M_sample = polynomial(zeros(n_size, n_size));

        hx_sample_pvar = subs(hx_pvar, x_vars_pvar, x_sample);

        for k = 1:length(trig_terms)
            trig_val = double(subs(trig_terms(k), x_vars, x_sample));
            hx_sample_pvar = subs(hx_sample_pvar, dummy_trig_vars_pvar(k), trig_val);
        end

        hx_sample = double(hx_sample_pvar);
        psi_sample = subs(safe_set_pvar, y_vars_pvar, hx_sample);
        M_sample(1, 1) = psi_sample;

        idx = 2;

        for i = 1:p

            for j = 1:(r_deg(i) - 1)
                term_ij = cert_term_dict(sprintf("output_%d_k%d", i, j));
                eta_kj_sym = term_ij{1}(2);

                for k_trig = 1:length(trig_terms)
                    eta_kj_sym = subs(eta_kj_sym, trig_terms(k_trig), dummy_trig_vars(k_trig));
                end

                eta_kj_pvar = sym2pvar(eta_kj_sym, ...
                    [x_vars; y_vars; k1_sym; J_k1_sym_flat; dummy_trig_vars], ...
                    [x_vars_pvar; y_vars_pvar; k1_pvar; J_k1_pvar_flat; dummy_trig_vars_pvar]);
                eta_kj_pvar = subs(eta_kj_pvar, k1_pvar, dpvar2poly(k1_poly));
                eta_kj_pvar = subs(eta_kj_pvar, J_k1_pvar_flat, J_k1_poly_flat);
                eta_kj_pvar = subs(eta_kj_pvar, y_vars_pvar, hx_pvar);
                eta_kj_val = subs(eta_kj_pvar, x_vars_pvar, x_sample);

                for k = 1:length(trig_terms)
                    trig_val = double(subs(trig_terms(k), x_vars, x_sample));
                    eta_kj_val = subs(eta_kj_val, dummy_trig_vars_pvar(k), trig_val);
                end

                M_sample(1, idx) = eta_kj_val;
                M_sample(idx, 1) = eta_kj_val;
                M_sample(idx, idx) = 2 * mu_val;
                idx = idx + 1;
            end

        end

        prog = sosmatrixineq(prog, M_sample);
    end

    disp('Added certificate LMI constraints (per-sample slack variant).');
end

% =============================================================================
% control limits with PER-SAMPLE slack
% =============================================================================
function prog = add_control_limit_constraints_persample(prog, ux, x_samples_valid, lb, ub, slack_cells, ...
        x_vars, y_vars, k1_sym, J_k1_sym_flat, x_vars_pvar, y_vars_pvar, k1_pvar, J_k1_pvar_flat, ...
        dummy_trig_vars_pvar, k1_poly, J_k1_poly_flat, hx_pvar, trig_terms, dummy_trig_vars)

    n_samples = size(x_samples_valid, 2);

    for i = 1:length(ux)
        [num, den] = numden(ux(i));
        den_func = matlabFunction(den, 'Vars', [x_vars]);
        x_samples_valid_cell = num2cell(x_samples_valid, 2);
        den_vals = den_func(x_samples_valid_cell{:});

        if isscalar(den_vals)
            den_vals = den_vals * ones(n_samples, 1);
        end

        for j = 1:length(trig_terms)
            num = subs(num, trig_terms(j), dummy_trig_vars(j));
        end

        num_pvar = sym2pvar(num, ...
            [x_vars; y_vars; k1_sym; J_k1_sym_flat; dummy_trig_vars], ...
            [x_vars_pvar; y_vars_pvar; k1_pvar; J_k1_pvar_flat; dummy_trig_vars_pvar]);
        num_pvar = subs(num_pvar, k1_pvar, dpvar2poly(k1_poly));
        num_pvar = subs(num_pvar, J_k1_pvar_flat, J_k1_poly_flat);
        num_pvar = subs(num_pvar, y_vars_pvar, hx_pvar);

        for j = 1:n_samples
            x_sample = x_samples_valid(:, j);
            den_val = den_vals(j);

            if abs(den_val) < 1e-6
                continue;
            end

            num_pvar_sample = subs(num_pvar, x_vars_pvar, x_sample);

            for k = 1:length(trig_terms)
                trig_term_val = double(subs(trig_terms(k), x_vars, x_sample));
                num_pvar_sample = subs(num_pvar_sample, dummy_trig_vars_pvar(k), trig_term_val);
            end

            sji = slack_cells{j, i};
            % lb(i) - slack <= u_i(x_j) <= ub(i) + slack
            prog = sosineq(prog, num_pvar_sample / den_val - lb(i) + sji);
            prog = sosineq(prog, ub(i) - num_pvar_sample / den_val + sji);
        end

    end

    disp('Added PER-SAMPLE slack control-limit constraints.');
end
