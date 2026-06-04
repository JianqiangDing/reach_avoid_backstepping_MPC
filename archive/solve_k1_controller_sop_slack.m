% SLACK (soft-constraint) variant of solve_k1_controller_sop.
%
% Difference from the original: the per-sample control-input bounds are RELAXED
% with a per-channel slack variable and the slack is minimized, so the program is
% ALWAYS feasible. The optimal slack reports how far the achievable control range
% exceeds the requested bound:  achievable bound_i = ub_i + slack_i.
%
%   |u_i(x_j)| <= ub_i + slack_i,   slack_i >= 0,   min ( delta + sum_i slack_i )
%
% The reach-avoid (certificate descent + Schur-complement non-negativity) remains a
% HARD constraint, so the reported slack is "tightest input bound while keeping the
% reach-avoid property". The original solve_k1_controller_sop.m is left unchanged.
function [k1_opt, J_k1_opt, k1_delta, slack_opt] = solve_k1_controller_sop_slack(ux, k1_sym, J_k1_sym, cert_term_dict, p, r_deg, ...
        x_vars, y_vars, hx, safe_set, target_set, x_samples_valid, ...
        lb, ub, ds, dv, k1_lambda, mu_val)

    m = length(ux);  % number of control inputs (one slack per channel)

    % create pvar variables for the x_vars
    x_vars_pvar = [];
    for i = 1:length(x_vars)
        x_vars_pvar = [x_vars_pvar; pvar(strcat(char(x_vars(i)), '_pvar'))];
    end

    % create pvar variables for the single-integrator system
    y_vars_pvar = [];
    for i = 1:p
        y_vars_pvar = [y_vars_pvar; pvar(strcat(char(y_vars(i)), '_pvar'))];
    end

    % decision variables: delta (descent slack, as in the original) + per-channel
    % control-bound slacks. Only m in {1,2} occur in the examples.
    if m == 1
        dpvar delta slack1;
        slack = slack1; decvars = [delta; slack1];
    else
        dpvar delta slack1 slack2;
        slack = [slack1; slack2]; decvars = [delta; slack1; slack2];
    end

    monom_k1 = monomials(y_vars_pvar, 0:dv);
    monom_sos = monomials(y_vars_pvar, 0:ds);
    prog = sosprogram(y_vars_pvar, decvars);

    % create polynomials for every element of the k1 controller
    k1_poly = [];
    for i = 1:p
        [prog, k1_i] = sospolyvar(prog, monom_k1);
        k1_poly = [k1_poly; k1_i];
    end

    % pvar variables for the k1 controller
    k1_pvar = [];
    for i = 1:p
        k1_pvar = [k1_pvar; pvar(strcat(char(k1_sym(i)), '_pvar'))];
    end

    % polynomial Jacobian of k1 w.r.t y
    J_k1_poly = [];
    for i = 1:p
        J_row = [];
        for j = 1:p
            J_row = [J_row, diff(k1_poly(i), y_vars_pvar(j))];
        end
        J_k1_poly = [J_k1_poly; J_row];
    end
    J_k1_poly_flat = dpvar2poly(J_k1_poly(:));
    J_k1_sym_flat = J_k1_sym(:);
    J_k1_pvar_flat = [];
    for i = 1:length(J_k1_sym_flat)
        J_k1_pvar_flat = [J_k1_pvar_flat; pvar(strcat(char(J_k1_sym_flat(i)), '_pvar'))];
    end

    % trig-term dummy variables (so trig functions become polynomial vars for SOS)
    trig_terms = detect_trigonometric_terms([ux; hx]);
    dummy_trig_vars = sym(zeros(length(trig_terms), 1));
    dummy_trig_vars_pvar = [];
    for i = 1:length(trig_terms)
        trig_term_str = char(trig_terms(i));
        trig_term_str = strrep(trig_term_str, '(', '_');
        trig_term_str = strrep(trig_term_str, ')', '');
        dummy_trig_vars(i) = str2sym(strcat('dummy_trig_var_', trig_term_str));
        dummy_trig_vars_pvar = [dummy_trig_vars_pvar; pvar(strcat('dummy_trig_var_', trig_term_str, '_pvar'))];
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

    % Lie derivative of the safe set along the single-integrator system: grad(S).k1
    grad_safe_set = [];
    for i = 1:p
        grad_safe_set = [grad_safe_set; diff(safe_set_pvar, y_vars_pvar(i))];
    end
    Lie_safe_set = grad_safe_set(1) * k1_poly(1);
    for i = 2:p
        Lie_safe_set = Lie_safe_set + grad_safe_set(i) * k1_poly(i);
    end

    % reach-avoid descent (HARD, unchanged from the original)
    sos_expr = Lie_safe_set - k1_lambda * safe_set_pvar + delta - s1 * safe_set_pvar - s2 * target_set_pvar;
    prog = sosineq(prog, sos_expr);
    prog = sosineq(prog, delta);
    prog = sosineq(prog, s1);
    prog = sosineq(prog, s2);

    % control-input bounds, now SOFT (relaxed by per-channel slack)
    for i = 1:m
        prog = sosineq(prog, slack(i));   % slack_i >= 0
    end
    prog = add_control_limit_constraints_slack(prog, ux, x_samples_valid, lb, ub, slack, x_vars, y_vars, k1_sym, J_k1_sym_flat, ...
        x_vars_pvar, y_vars_pvar, k1_pvar, J_k1_pvar_flat, dummy_trig_vars_pvar, ...
        k1_poly, J_k1_poly_flat, hx_pvar, trig_terms, dummy_trig_vars);

    % certificate non-negativity (HARD, unchanged from the original)
    prog = add_certificate_nonnegativity_constraints(prog, cert_term_dict, p, r_deg, x_samples_valid, ...
        x_vars, y_vars, x_vars_pvar, y_vars_pvar, k1_sym, J_k1_sym_flat, dummy_trig_vars, ...
        k1_pvar, J_k1_pvar_flat, dummy_trig_vars_pvar, k1_poly, J_k1_poly_flat, hx_pvar, ...
        trig_terms, safe_set_pvar, mu_val);

    % minimize the descent slack delta + the total control-bound slack
    prog = sossetobj(prog, delta + sum(slack));

    solver_opt.solver = 'mosek';
    prog = sossolve(prog, solver_opt);

    % extract k1
    k1_opt = [];
    for i = 1:p
        k1_opt_i = sosgetsol(prog, k1_poly(i));
        k1_opt = [k1_opt; poly2sym(k1_opt_i, y_vars_pvar, y_vars)];
    end
    J_k1_opt = jacobian(k1_opt, y_vars);
    k1_delta = double(sosgetsol(prog, delta));
    slack_opt = double(sosgetsol(prog, slack));

    if isempty(k1_delta) || isnan(k1_delta) || any(isnan(slack_opt))
        error('solve_k1_controller_sop_slack:solverFailed', ...
            'SOS solve returned no valid solution (delta/slack empty or NaN); check the Mosek solver status.');
    end
    fprintf('Slack solve: delta=%.4g, slack(achievable-bound margin) = [%s]\n', ...
        k1_delta, strtrim(sprintf('%.4g ', slack_opt)));
end

% ---- certificate non-negativity (identical to the original) ------------------
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
    disp('Added constraints from valid samples for non-negativity of reach-avoid certificate.');
end

% ---- control limits, RELAXED with per-channel slack --------------------------
function prog = add_control_limit_constraints_slack(prog, ux, x_samples_valid, lb, ub, slack, x_vars, y_vars, ...
        k1_sym, J_k1_sym_flat, x_vars_pvar, y_vars_pvar, k1_pvar, J_k1_pvar_flat, dummy_trig_vars_pvar, ...
        k1_poly, J_k1_poly_flat, hx_pvar, trig_terms, dummy_trig_vars)

    for i = 1:length(ux)
        [num, den] = numden(ux(i));
        den_func = matlabFunction(den, 'Vars', [x_vars]);
        x_samples_valid_cell = num2cell(x_samples_valid, 2);
        den_vals = den_func(x_samples_valid_cell{:});
        if isscalar(den_vals)
            den_vals = den_vals * ones(size(x_samples_valid, 2), 1);
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

        for j = 1:size(x_samples_valid, 2)
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
            % SOFT: lb(i) - slack(i) <= u_i <= ub(i) + slack(i)
            prog = sosineq(prog, num_pvar_sample / den_val - lb(i) + slack(i));
            prog = sosineq(prog, ub(i) - num_pvar_sample / den_val + slack(i));
        end
    end
    disp('Added RELAXED (slack) control-limit constraints over the sampled region.');
end
