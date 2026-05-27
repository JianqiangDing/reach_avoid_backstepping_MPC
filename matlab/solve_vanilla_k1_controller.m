% helper function for vanilla k1 controller solving without considering control input bounds
function [k1_opt, k1_lambda, k1_delta] = solve_vanilla_k1_controller(y_vars, safe_set, target_set, dv, ds)
    % solve the vanilla k1 controller for the single-integrator system without considering control input bounds
    % the single-integrator system is defined as y_dot = k1, where y is the output of the original system
    % here we aims to find a k1 controller such that the single-integrator system can satisfy the reach-avoid specification defined by the safe set and target set

    % INPUTS:
    % y_vars: output variables of the original system, symbolic vector of size (p x 1)
    % safe_set: symbolic expression for the safe set constraint, should be non-negative inside the safe set
    % target_set: symbolic expression for the target set constraint, should be non-positive out of the target set
    % dv: degree of the k1 controller polynomial
    % ds: degree of the auxiliary SOS polynomials for the single-integrator system

    % OUTPUTS:
    % k1_opt: obtained k1 controller for the single-integrator system, symbolic vector of size (p x 1)
    % k1_lambda: obtained lambda scalar for the single-integrator system, scalar decision variable in the SOS program
    % k1_delta: obtained delta scalar for the single-integrator system, scalar decision variable in the SOS program

    echo on;

    % get the number of outputs
    p = length(y_vars);
    % create pvar variables for the single-intgegrator system
    y_vars_pvar = [];

    for i = 1:p
        y_vars_pvar = [y_vars_pvar; pvar(strcat(char(y_vars(i)), '_pvar'))];
    end

    % define decision variable for solving k1 controller
    dpvar delta lambda;

    xi0 = 1e-8; % small positive value to ensure the SOS constraints are strictly feasible

    % define monomials for k1 controller and auxiliary SOS polynomials
    monom_k1 = monomials(y_vars_pvar, 0:dv);
    monom_sos = monomials(y_vars_pvar, 0:ds);

    % define the SOS program for solving k1 controller
    prog = sosprogram(y_vars_pvar, [delta; lambda]);

    % create polynomials for every single element in k1 controller
    k1 = [];

    for i = 1:p
        [prog, k1_i] = sospolyvar(prog, monom_k1);
        k1 = [k1; k1_i];
    end

    % create the pvar version of the safe set and target set by substituting y_vars with y_vars_pvar
    safe_set_pvar = sym2pvar(safe_set, y_vars, y_vars_pvar);
    target_set_pvar = sym2pvar(target_set, y_vars, y_vars_pvar);

    % % disp the safe_set_pvar and target_set_pvar to check the correctness of the substitution
    % disp('Safe set in pvar form:');
    % safe_set_pvar
    % disp('Target set in pvar form:');
    % target_set_pvar

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
    Lie_safe_set = grad_safe_set(1) * k1(1);

    for i = 2:p
        Lie_safe_set = Lie_safe_set + grad_safe_set(i) * k1(i);
    end

    % % disp the Lie derivative to check the correctness of the computation
    % disp('Lie derivative of the safe set constraint along the single-integrator system:');
    % disp(Lie_safe_set);

    % Construct the SOS constraint expression
    sos_expr = Lie_safe_set - lambda * safe_set_pvar + delta - s1 * safe_set_pvar - s2 * target_set_pvar;

    % disp('The expression for the SOS constraint for k1 controller:');
    % disp(sos_expr);

    % add the SOS constraints for k1 controller
    prog = sosineq(prog, sos_expr);
    prog = sosineq(prog, delta);
    prog = sosineq(prog, lambda - xi0); % ensure lambda is positive
    prog = sosineq(prog, s1);
    prog = sosineq(prog, s2);

    % set objective to minimize delta
    prog = sossetobj(prog, delta);

    % set solver options and solve the SOS program
    solver_opt.solver = 'mosek';
    prog = sossolve(prog, solver_opt);

    % extract the obtained k1 controller

    k1_opt = [];

    for i = 1:p
        k1_opt_i = sosgetsol(prog, k1(i));
        k1_opt = [k1_opt; poly2sym(k1_opt_i, y_vars_pvar, y_vars)];
    end

    % extract the obtained lambda and delta
    k1_lambda = double(sosgetsol(prog, lambda));
    k1_delta = double(sosgetsol(prog, delta));

    echo off;

end
