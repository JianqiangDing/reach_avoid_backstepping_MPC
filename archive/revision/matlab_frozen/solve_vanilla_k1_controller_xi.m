% Variant of solve_vanilla_k1_controller with the SOS floor xi0 exposed as an
% argument (original file untouched). Raising xi0 forces lambda >= xi0, so the
% vanilla k1 is RE-SOLVED (co-designed) consistently with a larger backstepping
% scale -- this is the proper test of "would a mu-scale lambda make k1 act?",
% WITHOUT fixing lambda to a constant (it stays a decision variable, floored at xi0).
%
% Returns the solver feasibility info so a sweep can detect where it becomes
% infeasible instead of erroring out.
function [k1_opt, k1_lambda, k1_delta, info] = solve_vanilla_k1_controller_xi(y_vars, safe_set, target_set, dv, ds, xi0)

    p = length(y_vars);
    y_vars_pvar = [];
    for i = 1:p
        y_vars_pvar = [y_vars_pvar; pvar(strcat(char(y_vars(i)), '_pvar'))]; %#ok<AGROW>
    end

    dpvar delta lambda;

    monom_k1 = monomials(y_vars_pvar, 0:dv);
    monom_sos = monomials(y_vars_pvar, 0:ds);

    prog = sosprogram(y_vars_pvar, [delta; lambda]);

    k1 = [];
    for i = 1:p
        [prog, k1_i] = sospolyvar(prog, monom_k1);
        k1 = [k1; k1_i]; %#ok<AGROW>
    end

    safe_set_pvar = sym2pvar(safe_set, y_vars, y_vars_pvar);
    target_set_pvar = sym2pvar(target_set, y_vars, y_vars_pvar);

    [prog, s1] = sospolyvar(prog, monom_sos);
    [prog, s2] = sospolyvar(prog, monom_sos);

    grad_safe_set = [];
    for i = 1:p
        grad_safe_set = [grad_safe_set; diff(safe_set_pvar, y_vars_pvar(i))]; %#ok<AGROW>
    end

    Lie_safe_set = grad_safe_set(1) * k1(1);
    for i = 2:p
        Lie_safe_set = Lie_safe_set + grad_safe_set(i) * k1(i);
    end

    sos_expr = Lie_safe_set - lambda * safe_set_pvar + delta - s1 * safe_set_pvar - s2 * target_set_pvar;

    prog = sosineq(prog, sos_expr);
    prog = sosineq(prog, delta);
    prog = sosineq(prog, lambda - xi0);   % <-- xi0 is now a parameter (mu-scale test)
    prog = sosineq(prog, s1);
    prog = sosineq(prog, s2);

    prog = sossetobj(prog, delta);

    solver_opt.solver = 'mosek';
    prog = sossolve(prog, solver_opt);

    % solver feasibility info (SeDuMi-style fields exposed by SOSTOOLS)
    info = struct('pinf', NaN, 'dinf', NaN, 'feasratio', NaN, 'numerr', NaN);
    if isfield(prog, 'solinfo') && isfield(prog.solinfo, 'info')
        si = prog.solinfo.info;
        for fn = {'pinf', 'dinf', 'feasratio', 'numerr'}
            if isfield(si, fn{1}); info.(fn{1}) = si.(fn{1}); end
        end
    end

    k1_opt = [];
    for i = 1:p
        k1_opt_i = sosgetsol(prog, k1(i));
        k1_opt = [k1_opt; poly2sym(k1_opt_i, y_vars_pvar, y_vars)]; %#ok<AGROW>
    end

    k1_lambda = double(sosgetsol(prog, lambda));
    k1_delta = double(sosgetsol(prog, delta));
end
