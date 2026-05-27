% Compute the lower and upper bounds of a rational polynomial over the zero superlevel set of another polynomial.
%
% The idea here is to determine the bounds separately for denominator positive and negative regions, to avoid unboundedness issues in the SOS program.
%  Specifically, we use the SOS relaxation:
% for denominator positive region:
%   num_poly - lb * den_poly - s1 * superlevel_set_poly  - s2 * den_poly is SOS; sos condition for lower bound certificate on the superlevel set defined by (superlevel_set_poly >= 0 and the denominator den_poly >= 0)
%   ub * den_poly - num_poly - s3 * superlevel_set_poly + s4 * den_poly is SOS; sos condition for upper bound certificate on the superlevel set defined by (superlevel_set_poly >= 0 and the denominator den_poly >= 0)
% for denominator negative region:
%   lb * den_poly - num_poly - s1 * superlevel_set_poly + s2 * den_poly is SOS; sos condition for lower bound certificate on the superlevel set defined by (superlevel_set_poly >= 0 and the denominator den_poly <= 0)
%   num_poly - ub * den_poly - s3 * superlevel_set_poly  - s4 * den_poly is SOS; sos condition for upper bound certificate on the superlevel set defined by (superlevel_set_poly >= 0 and the denominator den_poly <= 0)
% where s1, s2, s3, s4 are SOS multipliers.
% CAUTION: THIS FUNCTION ASSUMES THE ZERO SUPERLEVEL SET OF superlevel_set_poly IS COMPACT (BOUNDED) TO AVOID UNBOUNDEDNESS ISSUES IN THE SOS PROGRAM.

% INPUTS
%   num_poly            – symbolic polynomial (numerator)
%   den_poly            – symbolic polynomial (denominator)
%   superlevel_set_poly – symbolic polynomial defining the feasible region (the set is  {x : superlevel_set_poly(x) >= 0})
%   ds_min              – scalar integer; minimum degree for the SOS multiplier
%   epsilon             – small positive scalar; threshold to separate positive and negative denominator regions (i.e. we consider den_poly >= epsilon as positive region, and den_poly <= -epsilon as negative region)
% OUTPUTS
%   lower_bound – scalar double; the certified lower bound (returns NaN on solver failure)
%   upper_bound – scalar double; the certified upper bound (returns NaN on solver failure)
% NOTE: The feasible region should be compact (bounded) for the SOS problem
%       to be well-posed.  For unbounded sets add bounding-box multipliers following the pattern in demo_find_poly_bounds.m.
function [lower_bound, upper_bound] = compute_poly_bounds(num_poly, den_poly, superlevel_set_poly, ds_min, epsilon)
    % ── Convert sym inputs to pvar if needed ────────────────────────────────
    if isa(num_poly, 'sym') || isa(den_poly, 'sym') || isa(superlevel_set_poly, 'sym')
        % detect trig terms across all three expressions and replace with
        % dummy symbolic vars so that every coefficient is numeric
        all_sym = [num_poly; den_poly; superlevel_set_poly];
        trig_terms = detect_trigonometric_terms(all_sym);
        dummy_trig_syms = sym(zeros(length(trig_terms), 1));

        for i_trig = 1:length(trig_terms)
            trig_str = char(trig_terms(i_trig));
            trig_str = strrep(strrep(trig_str, '(', '_'), ')', '');
            dummy_trig_syms(i_trig) = str2sym(['dummy_trig_sos_' trig_str]);
        end

        for i_trig = 1:length(trig_terms)
            num_poly = subs(num_poly, trig_terms(i_trig), dummy_trig_syms(i_trig));
            den_poly = subs(den_poly, trig_terms(i_trig), dummy_trig_syms(i_trig));
            superlevel_set_poly = subs(superlevel_set_poly, trig_terms(i_trig), dummy_trig_syms(i_trig));
        end

        sym_vars = symvar(num_poly + den_poly + superlevel_set_poly);
        var_names = arrayfun(@char, sym_vars, 'UniformOutput', false);
        eval(['pvar ' strjoin(var_names, ' ')]);
        pvar_vars = eval(['[' strjoin(var_names, ',') ']']).';
        num_poly = sym2pvar(num_poly, sym_vars, pvar_vars);
        den_poly = sym2pvar(den_poly, sym_vars, pvar_vars);
        superlevel_set_poly = sym2pvar(superlevel_set_poly, sym_vars, pvar_vars);
    end

    % ── Extract pvar variable list from the polynomial objects ───────────────
    var_names = unique([num_poly.varname; den_poly.varname; superlevel_set_poly.varname]);
    eval(['pvar ' strjoin(var_names, ' ')]);
    vars = eval(['[' strjoin(var_names, ', ') ']']);

    % Four explicit calls; algebraic sign logic is encapsulated in the helper.
    lb_pos = solve_fraction_bound(vars, num_poly, den_poly, superlevel_set_poly, ds_min, epsilon, 1, 'lower', 1);
    ub_pos = solve_fraction_bound(vars, num_poly, den_poly, superlevel_set_poly, ds_min, epsilon, 1, 'upper', 1);

    lb_neg = solve_fraction_bound(vars, num_poly, den_poly, superlevel_set_poly, ds_min, epsilon, -1, 'lower', 1);
    ub_neg = solve_fraction_bound(vars, num_poly, den_poly, superlevel_set_poly, ds_min, epsilon, -1, 'upper', 1);

    lower_bound = min(lb_pos, lb_neg);
    upper_bound = max(ub_pos, ub_neg);

end

% =========================================================================
% Helper: solve one rational bound (lower or upper) for a fixed sign of Q
% sign_Q =  1 : denominator region  {den_poly >= epsilon}
% sign_Q = -1 : denominator region  {den_poly <= -epsilon}
% =========================================================================
function bound_val = solve_fraction_bound(poly_vars, num_poly, den_poly, S_poly, ds_min, epsilon, sign_Q, bound_type, is_quiet)
    ds = ds_min + mod(ds_min, 2);
    prog = sosprogram(poly_vars);

    dpvar b;
    prog = sosdecvar(prog, b);

    [prog, s_S] = sospolyvar(prog, monomials(poly_vars, 0:ds));
    [prog, s_Q] = sospolyvar(prog, monomials(poly_vars, 0:ds));

    % Archimedean ball: required for compactness so the SOS problem is well-posed
    % (radius 10 is sufficient to cover the unit disk)
    R2 = 10;
    ball_poly = R2 - sum(poly_vars .^ 2);
    [prog, s_ball] = sospolyvar(prog, monomials(poly_vars, 0:ds));

    % Core algebraic case analysis based on the sign of the denominator:
    if sign_Q == 1 % Region: den_poly >= epsilon  (denominator positive)
        Q_expr = den_poly - epsilon;

        if strcmp(bound_type, 'lower')
            % P/Q >= b  =>  P - b*Q >= 0  (Q positive, inequality preserved)
            main_eq = num_poly - b * den_poly;
            prog = sossetobj(prog, -b); % maximize lower bound
        else
            % P/Q <= b  =>  b*Q - P >= 0
            main_eq = b * den_poly - num_poly;
            prog = sossetobj(prog, b); % minimize upper bound
        end

    else % Region: den_poly <= -epsilon  (i.e. -den_poly >= epsilon, denominator negative)
        Q_expr = -den_poly - epsilon;

        if strcmp(bound_type, 'lower')
            % P/Q >= b  =>  b*Q - P >= 0  (Q negative, inequality flips!)
            main_eq = b * den_poly - num_poly;
            prog = sossetobj(prog, -b); % maximize lower bound
        else
            % P/Q <= b  =>  P - b*Q >= 0
            main_eq = num_poly - b * den_poly;
            prog = sossetobj(prog, b); % minimize upper bound
        end

    end

    % Assemble the S-procedure constraint
    prog = sosineq(prog, main_eq - s_S * S_poly - s_Q * Q_expr - s_ball * ball_poly);

    prog = sosineq(prog, s_S);
    prog = sosineq(prog, s_Q);
    prog = sosineq(prog, s_ball);

    solver_opt.solver = 'mosek';

    if exist('is_quiet', 'var') && is_quiet
        solver_opt.verbose = 0;
    else
        solver_opt.verbose = 1;
    end

    prog = sossolve(prog, solver_opt);

    b_sol = sosgetsol(prog, b);

    if isempty(b_sol)
        bound_val = NaN;
    else
        bound_val = double(b_sol);
    end

end
