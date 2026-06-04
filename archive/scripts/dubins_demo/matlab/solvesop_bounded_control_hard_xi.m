% Hard-bound + raised xi0 variant of solvesop_bounded_control_slack_xi.m
% (originals in matlab/ are untouched). Calls solve_k1_controller_sop_hard which
% drops the slack variables and forces |u_i| <= ub_i at samples as HARD constraints.
% Returns a `feasible` flag so callers can detect infeasibility instead of crashing.
function [ux_con, ux_uncon, cert_con, cert_uncon, valid_count, k1_con, feasible, k1_lambda, info] = ...
        solvesop_bounded_control_hard_xi( ...
        ux, k1_sym, J_k1_sym, mu, lambda, certificate, cert_term_dict, p, r_deg, ...
        x_vars, y_vars, hx, safe_set, target_set, mu_val, lb, ub, ds, dv, samples_num, bound_min, bound_max, xi0)

    [k1_y, k1_lambda, ~] = solve_vanilla_k1_controller_xi(y_vars, safe_set, target_set, dv, ds, xi0);
    J_k1_y = jacobian(k1_y, y_vars);

    certificate_subs = subs(certificate, k1_sym, k1_y);
    certificate_subs = subs(certificate_subs, y_vars, hx);
    certificate_subs = substitute_mu_lambda(certificate_subs, mu, lambda, mu_val, k1_lambda);

    x_samples = get_random_samples(samples_num, x_vars, y_vars, hx, safe_set, target_set, ...
        certificate_subs, bound_min, bound_max);
    x_samples_valid = x_samples;
    valid_count = size(x_samples, 2);
    fprintf('  HARD-xi SOP: xi0=%.3g -> lambda=%.4g, %d reach-avoid samples\n', ...
        xi0, k1_lambda, valid_count);

    ux_for_sop = subs(ux, y_vars, hx);
    ux_for_sop = substitute_mu_lambda(ux_for_sop, mu, lambda, mu_val, k1_lambda);

    [k1_con, J_k1_con, ~, feasible, info] = solve_k1_controller_sop_hard( ...
        ux_for_sop, k1_sym, J_k1_sym, cert_term_dict, p, r_deg, ...
        x_vars, y_vars, hx, safe_set, target_set, x_samples_valid, ...
        lb, ub, ds, dv, k1_lambda, mu_val);

    if ~feasible
        ux_con = sym(zeros(numel(ux), 1));
        cert_con = sym(0);
        ux_uncon = substitute_mu_lambda(subs(subs(subs(ux, k1_sym, k1_y), J_k1_sym, J_k1_y), y_vars, hx), ...
            mu, lambda, mu_val, k1_lambda);
        cert_uncon = substitute_mu_lambda(subs(subs(certificate, k1_sym, k1_y), y_vars, hx), ...
            mu, lambda, mu_val, k1_lambda);
        return;
    end

    ux_con = substitute_mu_lambda(subs(subs(subs(ux, k1_sym, k1_con), J_k1_sym, J_k1_con), y_vars, hx), ...
        mu, lambda, mu_val, k1_lambda);
    cert_con = substitute_mu_lambda(subs(subs(certificate, k1_sym, k1_con), y_vars, hx), ...
        mu, lambda, mu_val, k1_lambda);

    ux_uncon = substitute_mu_lambda(subs(subs(subs(ux, k1_sym, k1_y), J_k1_sym, J_k1_y), y_vars, hx), ...
        mu, lambda, mu_val, k1_lambda);
    cert_uncon = substitute_mu_lambda(subs(subs(certificate, k1_sym, k1_y), y_vars, hx), ...
        mu, lambda, mu_val, k1_lambda);
end

% =============================================================================
function u_subs = substitute_mu_lambda(expr, mu, lambda, mu_val, lambda_val)
    expr_subs = expr;
    for i = 1:numel(mu)
        if ~isempty(mu{i})
            expr_subs = subs(expr_subs, mu{i}, mu_val);
        end
    end
    u_subs = subs(expr_subs, lambda, lambda_val);
end

function x_samples = get_random_samples(num_samples, x_vars, y_vars, hx, safe_set, target_set, certificate, bound_min, bound_max)
    safe_set_x = subs(safe_set, y_vars, hx);
    target_set_x = subs(target_set, y_vars, hx);
    certificate_x = subs(certificate, y_vars, hx);
    safe_set_func = matlabFunction(safe_set_x, 'Vars', [x_vars]);
    target_set_func = matlabFunction(target_set_x, 'Vars', [x_vars]);
    certificate_func = matlabFunction(certificate_x, 'Vars', [x_vars]);
    n = length(x_vars);
    x_samples = bound_min' + (bound_max - bound_min)' .* rand(num_samples, n);
    x_T = x_samples';
    args = num2cell(x_T, 2);
    safe_set_values = safe_set_func(args{:});
    target_set_values = target_set_func(args{:});
    certificate_values = certificate_func(args{:});
    valid_indices = find(safe_set_values >= 0 & target_set_values >= 0 & certificate_values >= 0);
    x_samples = x_samples(valid_indices, :)';
end
