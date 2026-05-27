% helper function to solve the k1 controller for the single-integrator system considering control input bounds
function [k1_opt, J_k1_opt, k1_delta] = solve_k1_controller_sop(ux, k1_sym, J_k1_sym, cert_term_dict, p, r_deg, ...
        x_vars, y_vars, hx, safe_set, target_set, x_samples_valid, ...
        lb, ub, ds, dv, k1_lambda, mu_val)
    % solve the k1 controller for the single-integrator system considering control input bounds using scenario optimization programming (SOP) with SOS constraints
    % the single-integrator system is defined as y_dot = k1, where y is the output of the original system
    % here we aims to find a k1 controller such that the single-integrator system can satisfy the reach-avoid specification defined by the safe set and target set while ensuring the control input bounds are satisfied

    % INPUTS:
    % ux: original controller expression obtained from backstepping design, symbolic vector of size (m x 1), m is the control input dimension
    % k1_sym: symbolic vector for the k1 controller in the original controller expression, symbolic vector of size (p x 1), p is the output dimension
    % J_k1_sym: symbolic matrix for the Jacobian of k1 w.r.t. y, symbolic matrix of size (p x p)
    % cert_term_dict: dictionary of terms in the reach-avoid certificate
    % x_vars: state variables of the original system, symbolic vector of size (n x 1)
    % y_vars: output variables of the original system, symbolic vector of size (p x 1)
    % hx: output mapping from state variables to output variables, symbolic vector of size (p x 1)
    % safe_set: symbolic expression for the safe set constraint, should be non-negative inside the safe set
    % target_set: symbolic expression for the target set constraint, should be non-positive out of the target set
    % x_samples_valid: random samples for state variables that satisfy the safe set, target set constraint and certificate constraint, matrix of size (n x num_valid_samples)
    % lb: lower bounds for control inputs, vector of size (m x 1)
    % ub: upper bounds for control inputs, vector of size (m x 1)
    % ds: degree of the auxiliary SOS polynomials for the single-integrator system
    % dv: degree of the k1 controller polynomial
    % k1_lambda: the lambda value obtained from the vanilla k1 controller design without considering control input bounds,
    % which will be used as a fixed value for solving the bounded control inputs

    % create pvar variables  for the x_vars
    x_vars_pvar = [];

    for i = 1:length(x_vars)
        x_vars_pvar = [x_vars_pvar; pvar(strcat(char(x_vars(i)), '_pvar'))];
    end

    % create pvar variables for the single-intgegrator system
    y_vars_pvar = [];

    for i = 1:p
        y_vars_pvar = [y_vars_pvar; pvar(strcat(char(y_vars(i)), '_pvar'))];
    end

    % define decision variable for solving k1 controller
    dpvar delta; % only consider delta here, as we will fix lambda to be the same value as the one obtained from the vanilla k1 controller design without considering control input bounds

    % define monomials for k1 controller and auxiliary SOS polynomials
    monom_k1 = monomials(y_vars_pvar, 0:dv);
    monom_sos = monomials(y_vars_pvar, 0:ds);

    % define the SOS program for solving k1 controller
    % prog = sosprogram(y_vars_pvar, [delta; lambda]);
    prog = sosprogram(y_vars_pvar, delta);

    % create polynomials for every single element in k1 controller
    k1_poly = [];

    for i = 1:p
        [prog, k1_i] = sospolyvar(prog, monom_k1);
        k1_poly = [k1_poly; k1_i];
    end

    % create pvar variables for the k1 controller
    k1_pvar = [];

    for i = 1:p
        k1_pvar = [k1_pvar; pvar(strcat(char(k1_sym(i)), '_pvar'))];
    end

    % create polynomial for J_k1_sym, which is the Jacobian of k1 w.r.t. y
    J_k1_poly = [];

    for i = 1:p
        J_row = [];

        for j = 1:p
            J_row = [J_row, diff(k1_poly(i), y_vars_pvar(j))];
        end

        J_k1_poly = [J_k1_poly; J_row];
    end

    J_k1_poly_flat = dpvar2poly(J_k1_poly(:));

    % reshape J_k1_sym into a vector form to make it easier for substitution later
    J_k1_sym_flat = J_k1_sym(:);

    % create pvar variables for the Jacobian of k1 w.r.t. y
    J_k1_pvar_flat = [];

    for i = 1:length(J_k1_sym_flat)
        J_k1_pvar_flat = [J_k1_pvar_flat; pvar(strcat(char(J_k1_sym_flat(i)), '_pvar'))];
    end

    % detect all trigonometric terms in both ux and hx combined so that trig
    % functions appearing in the output map hx are also given dummy variables
    trig_terms = detect_trigonometric_terms([ux; hx]);
    % create dummy symbolic variables for all detected trigonometric terms with corresponding detected names, which will be used for substitution to convert the trigonometric terms into polynomial terms for SOS programming
    dummy_trig_vars = sym(zeros(length(trig_terms), 1));
    dummy_trig_vars_pvar = [];

    for i = 1:length(trig_terms)
        % first replace the "(" and ")" in the detected trig term with "_" to make it a valid variable name, then create a symbolic variable with the processed name
        % for example, replace "sin(x1)" with "sin_x1" to avoid issues when creating symbolic variables with names that contain parentheses
        trig_term_str = char(trig_terms(i));
        trig_term_str = strrep(trig_term_str, '(', '_');
        trig_term_str = strrep(trig_term_str, ')', '');
        dummy_trig_vars(i) = str2sym(strcat('dummy_trig_var_', trig_term_str));
        dummy_trig_vars_pvar = [dummy_trig_vars_pvar; pvar(strcat('dummy_trig_var_', trig_term_str, '_pvar'))];
    end

    % create the pvar polynomials for the system output mapping hx(x)
    % substitute any trig terms with dummy symbolic vars first so that hx is polynomial in x_vars
    hx_sub = hx;

    for i = 1:length(trig_terms)
        hx_sub = subs(hx_sub, trig_terms(i), dummy_trig_vars(i));
    end

    hx_pvar = polynomial(zeros(length(hx), 1));

    for i = 1:length(hx)
        hx_pvar(i) = sym2pvar(hx_sub(i), [x_vars; dummy_trig_vars], [x_vars_pvar; dummy_trig_vars_pvar]);
    end

    % create the pvar version of the safe set and target set by substituting y_vars with y_vars_pvar
    safe_set_pvar = sym2pvar(safe_set, y_vars, y_vars_pvar);
    target_set_pvar = sym2pvar(target_set, y_vars, y_vars_pvar);

    % define auxiliary SOS polynomials
    [prog, s1] = sospolyvar(prog, monom_sos);
    [prog, s2] = sospolyvar(prog, monom_sos);

    % before adding constraints, compute the Lie derivative of the safe_set constraint along the single-integrator system: L_f(safe_set) = ∂safe_set/∂y * k1

    % Compute gradient of safe_set as a pvar polynomial
    grad_safe_set = [];

    for i = 1:p
        grad_safe_set = [grad_safe_set; diff(safe_set_pvar, y_vars_pvar(i))];
    end

    % Compute Lie derivative as dot product: ∇S · k1
    % This keeps it as a proper dpvar polynomial
    Lie_safe_set = grad_safe_set(1) * k1_poly(1);

    for i = 2:p
        Lie_safe_set = Lie_safe_set + grad_safe_set(i) * k1_poly(i);
    end

    % Construct the SOS constraint expression
    sos_expr = Lie_safe_set - k1_lambda * safe_set_pvar + delta - s1 * safe_set_pvar - s2 * target_set_pvar;

    % add the SOS constraints for k1 controller
    prog = sosineq(prog, sos_expr);
    prog = sosineq(prog, delta);
    prog = sosineq(prog, s1);
    prog = sosineq(prog, s2);

    % add additional constraints to ensure the control input bounds are satisfied for the sampled states
    % ===========================================================================
    prog = add_control_limit_constraints(prog, ux, x_samples_valid, lb, ub, x_vars, y_vars, k1_sym, J_k1_sym_flat, ...
        x_vars_pvar, y_vars_pvar, k1_pvar, J_k1_pvar_flat, dummy_trig_vars_pvar, ...
        k1_poly, J_k1_poly_flat, hx_pvar, trig_terms, dummy_trig_vars);
    % ===========================================================================
    % add additional constraints to ensure the certificate non-negativity for the sampled states
    prog = add_certificate_nonnegativity_constraints(prog, cert_term_dict, p, r_deg, x_samples_valid, ...
        x_vars, y_vars, x_vars_pvar, y_vars_pvar, k1_sym, J_k1_sym_flat, dummy_trig_vars, ...
        k1_pvar, J_k1_pvar_flat, dummy_trig_vars_pvar, k1_poly, J_k1_poly_flat, hx_pvar, ...
        trig_terms, safe_set_pvar, mu_val);

    % ===========================================================================

    % set objective to minimize delta
    prog = sossetobj(prog, delta);

    % set solver options and solve the SOS program
    solver_opt.solver = 'mosek';
    prog = sossolve(prog, solver_opt);

    % extract the obtained k1 controller

    k1_opt = [];

    for i = 1:p
        k1_opt_i = sosgetsol(prog, k1_poly(i));
        k1_opt = [k1_opt; poly2sym(k1_opt_i, y_vars_pvar, y_vars)]; % vector of symbolic polynomials for k1 controller
    end

    J_k1_opt = jacobian(k1_opt, y_vars); % Jacobian of the obtained k1 controller w.r.t. y

    % extract the obtained delta
    k1_delta = double(sosgetsol(prog, delta)); % extract delta as a double value
end

% helper function to add constraints from valid samles for certificate non-negativity
function prog = add_certificate_nonnegativity_constraints(prog, cert_term_dict, p, r_deg, x_samples_valid, ...
        x_vars, y_vars, x_vars_pvar, y_vars_pvar, k1_sym, J_k1_sym_flat, dummy_trig_vars, ...
        k1_pvar, J_k1_pvar_flat, dummy_trig_vars_pvar, k1_poly, J_k1_poly_flat, hx_pvar, ...
        trig_terms, safe_set_pvar, mu_val)

    % compute total number of off-diagonal terms across all outputs
    term_sum = 0;

    for i = 1:p
        term_sum = term_sum + r_deg(i) - 1;
    end

    for sample_idx = 1:size(x_samples_valid, 2)
        x_sample = x_samples_valid(:, sample_idx);

        n_size = 1 + term_sum;
        M_sample = polynomial(zeros(n_size, n_size)); % polynomial matrix to hold dpvar/numeric entries

        % --- (1,1) entry: ψ(y(x_sample)) ---
        % evaluate hx at the sample point to get y numerically
        % first substitute state vars, then substitute any trig dummy pvars
        hx_sample_pvar = subs(hx_pvar, x_vars_pvar, x_sample);

        for k = 1:length(trig_terms)
            trig_val = double(subs(trig_terms(k), x_vars, x_sample));
            hx_sample_pvar = subs(hx_sample_pvar, dummy_trig_vars_pvar(k), trig_val);
        end

        hx_sample = double(hx_sample_pvar);
        % then evaluate safe_set_pvar at that y (still a dpvar because k1_poly may appear)
        psi_sample = subs(safe_set_pvar, y_vars_pvar, hx_sample);
        M_sample(1, 1) = psi_sample;

        % --- off-diagonal and diagonal blocks ---
        idx = 2;

        for i = 1:p

            for j = 1:(r_deg(i) - 1)
                term_ij = cert_term_dict(sprintf("output_%d_k%d", i, j));

                eta_kj_sym = term_ij{1}(2); % Lf^j h_i − k_j^i  — symbolic

                % substitute trig terms with dummy symbolic vars before sym2pvar
                % (same pattern as add_control_limit_constraints does for num)
                for k_trig = 1:length(trig_terms)
                    eta_kj_sym = subs(eta_kj_sym, trig_terms(k_trig), dummy_trig_vars(k_trig));
                end

                % Convert to pvar, substitute decision variables, then evaluate at x_sample
                % (same pattern as the control limit constraints)
                eta_kj_pvar = sym2pvar(eta_kj_sym, ...
                    [x_vars; y_vars; k1_sym; J_k1_sym_flat; dummy_trig_vars], ...
                    [x_vars_pvar; y_vars_pvar; k1_pvar; J_k1_pvar_flat; dummy_trig_vars_pvar]);
                eta_kj_pvar = subs(eta_kj_pvar, k1_pvar, dpvar2poly(k1_poly));
                eta_kj_pvar = subs(eta_kj_pvar, J_k1_pvar_flat, J_k1_poly_flat);
                eta_kj_pvar = subs(eta_kj_pvar, y_vars_pvar, hx_pvar);
                eta_kj_val = subs(eta_kj_pvar, x_vars_pvar, x_sample);
                % substitute any trig dummy variables
                for k = 1:length(trig_terms)
                    trig_val = double(subs(trig_terms(k), x_vars, x_sample));
                    eta_kj_val = subs(eta_kj_val, dummy_trig_vars_pvar(k), trig_val);
                end

                % fill symmetric off-diagonal entries (scalar dpvar)
                M_sample(1, idx) = eta_kj_val;
                M_sample(idx, 1) = eta_kj_val;

                % diagonal block: 2*mu_ij * I  (here I is 1×1 since h_i is scalar)
                M_sample(idx, idx) = 2 * mu_val; % 2*mu_ij  — numeric scalar
                % (for simplicity, we use same value for all mu paramters, but in theory, it should not)

                idx = idx + 1;
            end

        end

        % Add the per-sample LMI constraint M(x_sample) ≽ 0
        prog = sosmatrixineq(prog, M_sample);
    end

    disp('Added constraints from valid samples for non-negativity of reach-avoid certificate.');

end

% helper function to add constraints from valid samples for control limits
function prog = add_control_limit_constraints(prog, ux, x_samples_valid, lb, ub, x_vars, y_vars, ...
        k1_sym, J_k1_sym_flat, x_vars_pvar, y_vars_pvar, k1_pvar, J_k1_pvar_flat, dummy_trig_vars_pvar, ...
        k1_poly, J_k1_poly_flat, hx_pvar, trig_terms, dummy_trig_vars)

    for i = 1:length(ux)
        % get the numerator and denominator of the symbolic expression ux(i)
        [num, den] = numden(ux(i));
        % convert the den to a matlab function for later evaluation
        den_func = matlabFunction(den, 'Vars', [x_vars]);
        % evaluate the den at all valied samples at once using the converted matlab function
        x_samples_valid_cell = num2cell(x_samples_valid, 2); % convert the valid samples into a cell array for function evaluation
        den_vals = den_func(x_samples_valid_cell{:});
        % if den_vals is just a scalar, convert it to a vector with the same value for all samples to avoid issues in later processing
        if isscalar(den_vals)
            den_vals = den_vals * ones(size(x_samples_valid, 2), 1);
        end

        % replace all trigonometric terms in num with the corresponding dummy symbolic variables
        for j = 1:length(trig_terms)
            num = subs(num, trig_terms(j), dummy_trig_vars(j));
        end

        % convert the num into a pvar polynomial by substituting relevant symbolic variables with the corresponding pvar variables
        num_pvar = sym2pvar(num, ...
            [x_vars; y_vars; k1_sym; J_k1_sym_flat; dummy_trig_vars], ...
            [x_vars_pvar; y_vars_pvar; k1_pvar; J_k1_pvar_flat; dummy_trig_vars_pvar]);

        % substitute the k1_pvar variables with the k1_poly decision variables
        num_pvar = subs(num_pvar, k1_pvar, dpvar2poly(k1_poly));
        % substitute the J_k1_pvar_flat variables with the corresponding Jacobian polynomials computed from k1_poly
        num_pvar = subs(num_pvar, J_k1_pvar_flat, J_k1_poly_flat);
        % substitute the y_vars with the output mapping hx to make it a proper polynomial in terms of x_vars_pvar
        num_pvar = subs(num_pvar, y_vars_pvar, hx_pvar);

        % for all samples
        for j = 1:size(x_samples_valid, 2)
            x_sample = x_samples_valid(:, j);
            den_val = den_vals(j); % get the pre-evaluated denominator value for the sampled state from the den_vals array

            if abs(den_val) < 1e-6 % if the denominator is close to zero, we skip this sample to avoid numerical issues
                continue;
            end

            % now evaluate the num_pvar for the sampled state by substituting x_vars_pvar with the sampled state
            num_pvar_sample = subs(num_pvar, x_vars_pvar, x_sample);
            % substitute the dummy_trig_vars_pvar variables with the corresponding trigonometric functions evaluated at the sampled state
            for k = 1:length(trig_terms)
                trig_term_val = double(subs(trig_terms(k), x_vars, x_sample));
                num_pvar_sample = subs(num_pvar_sample, dummy_trig_vars_pvar(k), trig_term_val);
            end

            % then add the constraint for the control input bounds: lb(i) <= num_pvar/den_val <= ub(i)
            prog = sosineq(prog, num_pvar_sample / den_val - lb(i)); % ensure num_pvar/den_val >= lb(i) -> num_pvar/den_val - lb(i) >= 0
            prog = sosineq(prog, ub(i) - num_pvar_sample / den_val); % ensure num_pvar/den_val <= ub(i) -> ub(i) - num_pvar/den_val >= 0

        end

    end

    disp('Added constraints from valid samples for control limits.');
end
