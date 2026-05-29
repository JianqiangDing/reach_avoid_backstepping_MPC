% Test the slack (soft-constraint) SOP on the 2-DoF manipulator.
% Requests +-500; the optimal slack reports the achievable torque bound (= ub + slack)
% while keeping reach-avoid. Evaluated on the 5-D augmented state (x5 = x1+x2).
clc; clear; close all;
addpath(fileparts(mfilename('fullpath')));  % matlab/

syms x1 x2 x3 x4 x5 y1 y2;
x_vars_sym = [x1; x2; x3; x4; x5];
y_vars_sym = [y1; y2];
m1 = 1.0; m2 = 1.0; l1 = 4.0; l2 = 4.0; lc1 = 2.0; lc2 = 2.0; I1 = 0.02; I2 = 0.02; g = 9.81;
M11 = I1 + I2 + m1 * lc1^2 + m2 * (l1^2 + lc2^2 + 2 * l1 * lc2 * cos(x2));
M12 = m2 * (lc2^2 + l1 * lc2 * cos(x2)) + I2;
M22 = m2 * lc2^2 + I2;
M = [M11, M12; M12, M22];
C_mat = [-m2 * l1 * lc2 * sin(x2) * x4, -m2 * l1 * lc2 * sin(x2) * (x3 + x4); ...
    m2 * l1 * lc2 * sin(x2) * x3, 0];
G_vec = [(m1 * g * lc1 + m2 * g * l1) * cos(x1) + m2 * g * lc2 * cos(x5); ...
    m2 * g * lc2 * cos(x5)];
Minv = inv(M);
qddot_f = simplify(Minv * (-C_mat * [x3; x4] - G_vec));
fx_sym = [x3; x4; qddot_f; x3 + x4];
gx_sym = [zeros(2, 2); Minv; zeros(1, 2)];
hx_sym = [l1 * cos(x1) + l2 * cos(x5); l1 * sin(x1) + l2 * sin(x5)];
safe_set_sym = -((4 * (y1 - 2) - 2 * y2^3)^2) + 0.8 * y2^3 + 10;
target_set_sym = ((y1 - 2 - 3.5)^2 / 1.2^2) + ((y2 - 1.8)^2 / 0.4^2) - 2;

[u, k1, J_k1, mu, lambda, certificate, cert_term_dict, A_matrix, b_vector, ks, p, r_deg] = ...
    reach_avoid_controller(fx_sym, gx_sym, hx_sym, x_vars_sym, y_vars_sym, safe_set_sym);

rng(42);
lb = [-500; -500]; ub = [500; 500];      % desired bound
ds = 4; dv = 2; mu_val = 15; samples_num = 10000;
bound_min = [-2.0; 0.8; -0.5; -0.5; -1.2];  % FL region (q2 in [0.8, pi-0.8])
bound_max = [0.5; pi - 0.8; 0.5; 0.5; pi - 0.3];

[u_opt, certificate_opt, valid_count, k1_opt, slack_opt] = solvesop_bounded_control_slack( ...
    u, k1, J_k1, mu, lambda, certificate, cert_term_dict, p, r_deg, ...
    x_vars_sym, y_vars_sym, hx_sym, safe_set_sym, target_set_sym, ...
    mu_val, lb, ub, ds, dv, samples_num, bound_min, bound_max);

[lb1, ub1] = compute_poly_bounds_sampling(x_vars_sym, u_opt(1), certificate_opt, 10000, bound_min, bound_max);
[lb2, ub2] = compute_poly_bounds_sampling(x_vars_sym, u_opt(2), certificate_opt, 10000, bound_min, bound_max);

fprintf('\n=== MANIPULATOR SLACK RESULT (mu=%.3g) ===\n', mu_val);
fprintf('requested  ub        = [%g, %g]\n', ub(1), ub(2));
fprintf('slack_opt            = [%.4g, %.4g]\n', slack_opt(1), slack_opt(2));
fprintf('achievable bound u+s = [%.4g, %.4g]\n', ub(1) + slack_opt(1), ub(2) + slack_opt(2));
fprintf('slack ctrl range     : tau1=[%.4g, %.4g]  tau2=[%.4g, %.4g]\n', lb1, ub1, lb2, ub2);
fprintf('(unconstrained reference ~ tau1 in [-5000, 7000])\n');
