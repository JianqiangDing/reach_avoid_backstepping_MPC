% function to solve bounded control for given reach-avoid backstepping controller
function [ux_opt, certificate_opt, valid_count, k1_opt] = solvesop_bounded_control(ux, k1_sym, J_k1_sym, mu, lambda, certificate, cert_term_dict, p, r_deg, ...
        x_vars, y_vars, hx, safe_set, target_set, mu_val, lb, ub, ds, dv, samples_num, bound_min, bound_max)

    % INPUTS:
    % ux: original controller expression obtained from backstepping design, symbolic vector of size (m x 1), m is the control input dimension
    % k1_sym: symbolic vector for the k1 controller in the original controller expression, symbolic vector of size (p x 1), p is the output dimension
    % J_k1: Jacobian of k1 w.r.t. y (symbolic, p x p)
    % certificate: symbolic expression for the certificate obtained from backstepping design, should be non-negative for the valid state space
    % cert_term_dict: dictionary of terms in the reach-avoid certificate
    % p: number of outputs
    % r_deg: vector of relative degrees for each output
    % x_vars: state variables of the original system, symbolic vector of size (n x 1)
    % y_vars: output variables of the original system, symbolic vector of size (p x 1)
    % hx: output mapping from state variables to output variables, symbolic vector of size (p x 1)
    % safe_set: symbolic expression for the safe set constraint, should be non-negative inside the safe set
    % target_set: symbolic expression for the target set constraint, should be non-positive out of the target set
    % Ax: matrix for linear constraints on control inputs, matrix of size (num_constraints x m)
    % lb: lower bounds for control inputs, vector of size (m x 1)
    % ub: upper bounds for control inputs, vector of size (m x 1)
    % ds: degree of the auxiliary SOS polynomials for the single-integrator system
    % dv: degree of the k1 controller polynomial
    % bound_min: lower bounds for sampling the state space for finding valid samples that satisfy the set (safe, target, vanilla reach-avoid certificate) constraints, vector of size (n x 1)
    % bound_max: maximum lower bounds for sampling the state space for finding valid samples that satisfy the set (safe, target, vanilla reach-avoid certificate) constraints, vector of size (n x 1)

    [k1_y, k1_lambda, k1_delta] = solve_vanilla_k1_controller(y_vars, safe_set, target_set, dv, ds);

    % compute the Jacobian of k1 w.r.t. y
    J_k1_y = jacobian(k1_y, y_vars);

    % substitute the obtained k1 controller, Jacobian, and output mapping into the certificate expression
    certificate_subs = subs(certificate, k1_sym, k1_y);
    % certificate_subs = subs(certificate_subs, J_k1_sym, J_k1_y); % not necessary to substitute the Jacobian of k1 for solving the bounded control inputs since the Jacobian only appears in the derivative of the certificate function
    certificate_subs = subs(certificate_subs, y_vars, hx);

    certificate_subs = substitute_mu_lambda(certificate_subs, mu, lambda, mu_val, k1_lambda); % use the computed value from the vanilla k1 controller design as the lambda value for solving the bounded control inputs

    % sample some random samples from the state space that satisfy the safe set, target set constraint and certificate constraint
    x_samples = get_random_samples(samples_num, x_vars, y_vars, hx, safe_set, target_set, certificate_subs, bound_min, bound_max);

    % substitute the obtained k1 controller and output mapping into the original controller expression to get pseudo ux for scenario optimization programming
    ux_pseudo = subs(ux, k1_sym, k1_y);
    ux_pseudo = subs(ux_pseudo, J_k1_sym, J_k1_y);
    ux_pseudo = subs(ux_pseudo, y_vars, hx);

    % also substitute the obtained mu_val and lambda_val into the original controller expression

    ux_pseudo = substitute_mu_lambda(ux_pseudo, mu, lambda, mu_val, k1_lambda); % use the computed value from the vanilla k1 controller design as the lambda value for solving the bounded control inputs

    % >>>>>>>>>>>>>>>>>>>>>> DEBUG <<<<<<<<<<<<<<<<<<<<
    % export the obtained unconstrained controller and the certificate for comparison later
    % If an auxiliary variable x5 = x1+x2 was introduced to work around the
    % compound-trig-argument limitation, substitute it back before exporting
    % so the Python output uses only the natural state variables.
    ux_pseudo_exp = ux_pseudo;
    certificate_exp = certificate_subs;

    if any(arrayfun(@(v) isequal(v, sym('x5')), x_vars))
        ux_pseudo_exp = subs(ux_pseudo_exp, sym('x5'), sym('x1') + sym('x2'));
        certificate_exp = subs(certificate_exp, sym('x5'), sym('x1') + sym('x2'));
    end

    params_for_export = struct();
    params_for_export.lb = lb;
    params_for_export.ub = ub;
    params_for_export.mu_val = mu_val;
    params_for_export.ds = ds;
    params_for_export.dv = dv;
    export_to_python(ux_pseudo_exp, certificate_exp, k1_y, params_for_export, 'sop_bounded_control_unconstrained_controller.py');
    % >>>>>>>>>>>>>>>>>>>>>> DEBUG <<<<<<<<<<<<<<<<<<<<

    % select the samples that satisfy the control input bounds for the pseudo ux

    % convert the samples to a format that can be substituted into the symbolic expression (each column is a sample)
    args = num2cell(x_samples, 2); % convert to cell array of row vectors for substitution
    % convert the ux_pseudo to a MATLAB function for evaluation
    ux_pseudo_func = matlabFunction(ux_pseudo, 'Vars', [x_vars]);
    % evaluate the pseudo ux for all samples in a vectorized manner
    ux_pseudo_values = ux_pseudo_func(args{:}); % this will return a matrix of size (m x num_samples)
    % check which samples satisfy the control input bounds
    % use all(..., 1) to reduce over the control dimension (rows) so that the
    % result is always a (1 x num_samples) logical vector, even when m = 1
    valid_indices = all(ux_pseudo_values >= lb, 1) & all(ux_pseudo_values <= ub, 1);
    % select the valid samples
    x_samples_valid = x_samples(:, valid_indices);
    valid_count = sum(valid_indices);
    fprintf('SOP Feasibility Report: %d Valid / %d Total (%.2f%%)\n', valid_count, samples_num, (valid_count / samples_num) * 100);

    ux_for_sop = subs(ux, y_vars, hx);
    ux_for_sop = substitute_mu_lambda(ux_for_sop, mu, lambda, mu_val, k1_lambda);

    % then we can use the valid samples to solve the scenario optimization programming (SOP) with SOS constraints to find a feasible k1 controller
    [k1_opt, J_k1_opt, k1_delta_opt] = solve_k1_controller_sop(ux_for_sop, k1_sym, J_k1_sym, cert_term_dict, p, r_deg, x_vars, y_vars, ...
        hx, safe_set, target_set, x_samples_valid, lb, ub, ds, dv, k1_lambda, mu_val);

    % substitute the obtainted k1 controller into the original controller expression ux to get the final controller with bounded control inputs
    ux_opt = subs(ux, k1_sym, k1_opt);
    ux_opt = subs(ux_opt, J_k1_sym, J_k1_opt);
    ux_opt = subs(ux_opt, y_vars, hx);
    % substitute the obtained k1_lambda_opt and k1_delta_opt into the final controller expression ux_opt
    ux_opt = substitute_mu_lambda(ux_opt, mu, lambda, mu_val, k1_lambda); % use the computed value from the SOP design as the lambda value for the final controller expression

    % substitute the obtainted k1_opt controller into the original certificate expression to get the final certificate with bounded control inputs
    certificate_opt = subs(certificate, k1_sym, k1_opt);
    certificate_opt = subs(certificate_opt, y_vars, hx);

    certificate_opt = substitute_mu_lambda(certificate_opt, mu, lambda, mu_val, k1_lambda); % use the computed value from the SOP design as the lambda value for the final certificate expression

    % DEBUG
    % check the corresponding certificate value & ux value for the sampled valid states
    % convert the certificate_opt to a MATLAB function for evaluation
    certificate_func = matlabFunction(certificate_opt, 'Vars', [x_vars]);
    ux_opt_func = matlabFunction(ux_opt, 'Vars', [x_vars]);

    % evaluate the certificate and ux_opt for all valid samples in a vectorized manner
    certificate_values = certificate_func(args{:}); % this will return a row vector of size (1 x num_valid_samples)
    ux_opt_values = ux_opt_func(args{:}); % this will return a matrix of size (m x num_valid_samples)

    % count the number of valid samples that satisfy the certificate constraint (certificate value >= 0) and control input bounds for the final controller
    valid_certificate_count = sum(certificate_values >= 0 & all(ux_opt_values >= lb) & all(ux_opt_values <= ub));
    fprintf('Final Controller Report: %d Valid Certificate / %d Total (%.2f%%)\n', ...
        valid_certificate_count, size(x_samples, 2), (valid_certificate_count / size(x_samples, 2)) * 100);
    % DEBUG

end

% helper functikon to substitute mu values and lambda with predefined numerical values in the controller expression
function u_subs = substitute_mu_lambda(expr, mu, lambda, mu_val, lambda_val)
    expr_subs = expr;

    % substitute mu values
    for i = 1:numel(mu)

        if ~isempty(mu{i})
            expr_subs = subs(expr_subs, mu{i}, mu_val);
        end

    end

    % substitute lambda value
    u_subs = subs(expr_subs, lambda, lambda_val);

end

% helper function to get random samples from the state space that satisfy the safe set, target set constraint and certificate constraint
function x_samples = get_random_samples(num_samples, x_vars, y_vars, hx, safe_set, target_set, certificate, bound_min, bound_max)
    % num_samples: number of random samples to generate
    % x_vars: state variables, symbolic vector of size (n x 1)
    % y_vars: output variables, symbolic vector of size (p x 1)
    % safe_set: symbolic expression for the safe set constraint, should be non-negative inside the safe set
    % target_set: symbolic expression for the target set constraint, should be non-positive inside the target set
    % certificate: symbolic expression for the certificate, should be non-negative for the valid state space
    % x_samples: generated random samples for state variables, matrix of size (n x num_samples)
    % y_samples: generated random samples for output variables, matrix of size (p x num_samples)
    % bound_min: lower bounds for sampling the state space, vector of size (n x 1)
    % bound_max: upper bounds for sampling the state space, vector of size (n x 1)

    % convert the given symbolic expressions to MATLAB functions for evaluation
    safe_set_x = subs(safe_set, y_vars, hx);
    target_set_x = subs(target_set, y_vars, hx);
    certificate_x = subs(certificate, y_vars, hx);
    safe_set_func = matlabFunction(safe_set_x, 'Vars', [x_vars]);
    target_set_func = matlabFunction(target_set_x, 'Vars', [x_vars]);
    certificate_func = matlabFunction(certificate_x, 'Vars', [x_vars]);

    n = length(x_vars);
    p = length(y_vars);

    % generate random samples with specific bounds  (result: num_samples × n)
    x_samples = bound_min' + (bound_max - bound_min)' .* rand(num_samples, n);

    % evaluate the safe set, target set, and certificate functions for the generated samples
    % x_samples is (num_samples × n); transpose to (n × num_samples) then split each row
    % into a separate argument to match the matlabFunction signature f(x1, x2, ..., xn)
    x_T = x_samples'; % n × num_samples
    args = num2cell(x_T, 2); % n×1 cell, each cell is a 1×num_samples row vector
    safe_set_values = safe_set_func(args{:});
    target_set_values = target_set_func(args{:});
    certificate_values = certificate_func(args{:});

    % select the samples that satisfy the safe set constraint, target set constraint, and certificate constraint
    valid_indices = find(safe_set_values >= 0 & target_set_values >= 0 & certificate_values >= 0);
    x_samples = x_samples(valid_indices, :)';

end
