% Hard-bound variant of solve_k1_controller_sop_slack — no slack variables,
% hard per-sample input bounds |u_i(x_j)| <= ub_i. The problem can be INFEASIBLE
% when ub is too tight at the current (mu, xi0); the caller detects this via the
% returned `feasible` flag and Mosek's pinf/dinf.
%
%   reach-avoid descent + non-negativity:  HARD (unchanged)
%   per-sample input bound:                HARD  (no slack)
%   objective:                             minimize delta
%
% The slack-variant file (matlab/solve_k1_controller_sop_slack.m) is untouched.
function [k1_opt, J_k1_opt, k1_delta, feasible, info] = solve_k1_controller_sop_hard( ...
        ux, k1_sym, J_k1_sym, cert_term_dict, p, r_deg, ...
        x_vars, y_vars, hx, safe_set, target_set, x_samples_valid, ...
        lb, ub, ds, dv, k1_lambda, mu_val)

    % --- create pvar variables --------------------------------------------------
    x_vars_pvar = [];
    for i = 1:length(x_vars)
        x_vars_pvar = [x_vars_pvar; pvar(strcat(char(x_vars(i)), '_pvar'))]; %#ok<AGROW>
    end
    y_vars_pvar = [];
    for i = 1:p
        y_vars_pvar = [y_vars_pvar; pvar(strcat(char(y_vars(i)), '_pvar'))]; %#ok<AGROW>
    end

    % --- decision variable: only delta (no slack) ------------------------------
    dpvar delta;
    decvars = delta;

    monom_k1 = monomials(y_vars_pvar, 0:dv);
    monom_sos = monomials(y_vars_pvar, 0:ds);
    prog = sosprogram(y_vars_pvar, decvars);

    % create polynomials for every element of the k1 controller
    k1_poly = [];
    for i = 1:p
        [prog, k1_i] = sospolyvar(prog, monom_k1);
        k1_poly = [k1_poly; k1_i]; %#ok<AGROW>
    end

    % pvar substitution variables
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

    % trig substitution
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

    % reach-avoid descent (HARD)
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

    % control-input bounds — HARD (no slack)
    prog = add_control_limit_constraints_hard( ...
        prog, ux, x_samples_valid, lb, ub, x_vars, y_vars, k1_sym, J_k1_sym_flat, ...
        x_vars_pvar, y_vars_pvar, k1_pvar, J_k1_pvar_flat, dummy_trig_vars_pvar, ...
        k1_poly, J_k1_poly_flat, hx_pvar, trig_terms, dummy_trig_vars);

    % certificate non-negativity (HARD)
    prog = add_certificate_nonnegativity_constraints( ...
        prog, cert_term_dict, p, r_deg, x_samples_valid, ...
        x_vars, y_vars, x_vars_pvar, y_vars_pvar, k1_sym, J_k1_sym_flat, dummy_trig_vars, ...
        k1_pvar, J_k1_pvar_flat, dummy_trig_vars_pvar, k1_poly, J_k1_poly_flat, hx_pvar, ...
        trig_terms, safe_set_pvar, mu_val);

    prog = sossetobj(prog, delta);

    solver_opt.solver = 'mosek';
    info = struct('pinf', NaN, 'dinf', NaN, 'feasratio', NaN, 'numerr', NaN);
    feasible = false;
    try
        prog = sossolve(prog, solver_opt);
        if isfield(prog, 'solinfo') && isfield(prog.solinfo, 'info')
            si = prog.solinfo.info;
            for fn = {'pinf', 'dinf', 'feasratio', 'numerr'}
                if isfield(si, fn{1}); info.(fn{1}) = si.(fn{1}); end
            end
        end
        feasible = (isnan(info.pinf) || info.pinf == 0) && ...
                   (isnan(info.dinf) || info.dinf == 0);
    catch ME
        fprintf(2, 'hard SOP solver error: %s\n', ME.message);
        feasible = false;
    end

    if ~feasible
        k1_opt = sym(zeros(p, 1));
        J_k1_opt = sym(zeros(p, p));
        k1_delta = NaN;
        fprintf('HARD SOP: INFEASIBLE  (pinf=%g dinf=%g feasratio=%g)\n', ...
            info.pinf, info.dinf, info.feasratio);
        return;
    end

    % extract
    k1_opt = [];
    for i = 1:p
        k1_opt_i = sosgetsol(prog, k1_poly(i));
        k1_opt = [k1_opt; poly2sym(k1_opt_i, y_vars_pvar, y_vars)]; %#ok<AGROW>
    end
    J_k1_opt = jacobian(k1_opt, y_vars);
    k1_delta = double(sosgetsol(prog, delta));
    fprintf('HARD SOP: FEASIBLE  delta=%.4g  (pinf=%g dinf=%g feasratio=%g)\n', ...
        k1_delta, info.pinf, info.dinf, info.feasratio);
end

% =============================================================================
% certificate non-negativity (identical to the slack file)
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
    disp('Added certificate non-negativity LMI constraints (hard).');
end

% =============================================================================
% control limits — HARD (no slack)
% =============================================================================
function prog = add_control_limit_constraints_hard(prog, ux, x_samples_valid, lb, ub, x_vars, y_vars, ...
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
            % HARD: lb(i) <= u_i <= ub(i)
            prog = sosineq(prog, num_pvar_sample / den_val - lb(i));
            prog = sosineq(prog, ub(i) - num_pvar_sample / den_val);
        end
    end
    disp('Added HARD control-limit constraints over the sampled region.');
end
