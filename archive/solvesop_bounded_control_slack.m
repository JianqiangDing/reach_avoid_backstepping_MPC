% SLACK (soft-constraint) variant of solvesop_bounded_control.
%
% Two changes vs the original (which is left untouched):
%   1. The control-input bounds are enforced over the WHOLE sampled reach-avoid
%      region -- it does NOT pre-filter to samples where the unconstrained
%      controller already satisfies the bound (the original's valid_indices step),
%      so the optimization actually sees the violating region.
%   2. The bounds are relaxed by a per-channel slack (see solve_k1_controller_sop_slack),
%      so the program is always feasible and the optimal slack reports the tightest
%      input bound achievable while keeping the reach-avoid property:
%          achievable bound_i = ub_i + slack_opt_i.
%
% Returns the extra output slack_opt. Does not export controllers (test/research
% variant; the caller evaluates the returned ux_opt directly).
function [ux_opt, certificate_opt, valid_count, k1_opt, slack_opt] = solvesop_bounded_control_slack( ...
        ux, k1_sym, J_k1_sym, mu, lambda, certificate, cert_term_dict, p, r_deg, ...
        x_vars, y_vars, hx, safe_set, target_set, mu_val, lb, ub, ds, dv, samples_num, bound_min, bound_max)

    [k1_y, k1_lambda, ~] = solve_vanilla_k1_controller(y_vars, safe_set, target_set, dv, ds);
    J_k1_y = jacobian(k1_y, y_vars);

    % certificate in x-space with the vanilla k1 (used to define the reach-avoid sampling region)
    certificate_subs = subs(certificate, k1_sym, k1_y);
    certificate_subs = subs(certificate_subs, y_vars, hx);
    certificate_subs = substitute_mu_lambda(certificate_subs, mu, lambda, mu_val, k1_lambda);

    % sample the reach-avoid region: safe >= 0, outside target, certificate >= 0
    x_samples = get_random_samples(samples_num, x_vars, y_vars, hx, safe_set, target_set, ...
        certificate_subs, bound_min, bound_max);

    % KEY: use ALL reach-avoid samples (no already-feasible pre-filter)
    x_samples_valid = x_samples;
    valid_count = size(x_samples, 2);
    fprintf('Slack SOP: %d reach-avoid samples (whole region, no already-feasible pre-filter)\n', valid_count);

    ux_for_sop = subs(ux, y_vars, hx);
    ux_for_sop = substitute_mu_lambda(ux_for_sop, mu, lambda, mu_val, k1_lambda);

    [k1_opt, J_k1_opt, ~, slack_opt] = solve_k1_controller_sop_slack(ux_for_sop, k1_sym, J_k1_sym, ...
        cert_term_dict, p, r_deg, x_vars, y_vars, hx, safe_set, target_set, x_samples_valid, ...
        lb, ub, ds, dv, k1_lambda, mu_val);

    % final controller and certificate with the solved k1
    ux_opt = subs(ux, k1_sym, k1_opt);
    ux_opt = subs(ux_opt, J_k1_sym, J_k1_opt);
    ux_opt = subs(ux_opt, y_vars, hx);
    ux_opt = substitute_mu_lambda(ux_opt, mu, lambda, mu_val, k1_lambda);

    certificate_opt = subs(certificate, k1_sym, k1_opt);
    certificate_opt = subs(certificate_opt, y_vars, hx);
    certificate_opt = substitute_mu_lambda(certificate_opt, mu, lambda, mu_val, k1_lambda);
end

% ---- helpers copied from solvesop_bounded_control.m (file-scoped) ------------
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
