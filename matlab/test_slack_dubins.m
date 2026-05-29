% Test the slack (soft-constraint) SOP on the Dubins car.
% Requests the tight bound +-5; the optimal slack reports how far the achievable
% control range exceeds it while keeping reach-avoid (achievable = ub + slack).
clc; clear; close all;
addpath(fileparts(mfilename('fullpath')));  % matlab/

syms x1 x2 th v y1 y2;
x_vars_sym = [x1; x2; th; v];
fx_sym = [v * cos(th); v * sin(th); 0; 0];
gx_sym = [0, 0; 0, 0; 1, 0; 0, 1];
hx_sym = [x1; x2];
y_vars_sym = [y1; y2];
h_raw = -(y1^4 + y2^4 - 16) * (y1^4 + y2^4 - 4);
target_set_sym = (y2 - 0)^2 + ((y1 + 1.7) / 0.5)^2 - 0.4;
alpha = 1e-3 * (-target_set_sym + 300);
safe_set_sym = alpha * h_raw;

[u, k1, J_k1, mu, lambda, certificate, cert_term_dict, A_matrix, b_vector, ks, p, r_deg] = ...
    reach_avoid_controller(fx_sym, gx_sym, hx_sym, x_vars_sym, y_vars_sym, safe_set_sym);

rng(42);
lb = [-5; -5]; ub = [5; 5];          % desired tight bound
ds = 4; dv = 4; mu_val = 0.1; samples_num = 1000;
bound_min = [-2; -2; 2 * pi / 3; 0.1];  % FL region (v >= 0.1)
bound_max = [2; 2; 4 * pi / 3; 1.0];

[u_opt, certificate_opt, valid_count, k1_opt, slack_opt] = solvesop_bounded_control_slack( ...
    u, k1, J_k1, mu, lambda, certificate, cert_term_dict, p, r_deg, ...
    x_vars_sym, y_vars_sym, hx_sym, safe_set_sym, target_set_sym, ...
    mu_val, lb, ub, ds, dv, samples_num, bound_min, bound_max);

[lb1, ub1] = compute_poly_bounds_sampling(x_vars_sym, u_opt(1), certificate_opt, 10000, bound_min, bound_max);
[lb2, ub2] = compute_poly_bounds_sampling(x_vars_sym, u_opt(2), certificate_opt, 10000, bound_min, bound_max);

fprintf('\n=== DUBINS SLACK RESULT (mu=%.3g) ===\n', mu_val);
fprintf('requested  ub        = [%g, %g]\n', ub(1), ub(2));
fprintf('slack_opt            = [%.4g, %.4g]\n', slack_opt(1), slack_opt(2));
fprintf('achievable bound u+s = [%.4g, %.4g]\n', ub(1) + slack_opt(1), ub(2) + slack_opt(2));
fprintf('slack ctrl range     : u1=[%.4g, %.4g]  u2=[%.4g, %.4g]\n', lb1, ub1, lb2, ub2);
fprintf('(unconstrained reference ~ u1 in [-90, 95])\n');
