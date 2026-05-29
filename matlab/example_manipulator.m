clc;
clear;
close all;

% ── 2-DoF Planar Manipulator (Acrobot-like) ───────────────────────────────────
%   State  x = [q1; q2; dq1; dq2; q1+q2]   (5-D augmented; see note below)
%   Input  u = [tau1; tau2]                  (joint torques)
%   Output y = h(x) = [l1*cos(x1)+l2*cos(x5); l1*sin(x1)+l2*sin(x5)]
%
%   NOTE on augmented state:  x5 = x1 + x2  is introduced as an independent
%   symbol so that every trig argument is a single symbolic variable.
%   The framework's detect_trigonometric_terms maps each trig argument to its
%   variable list; compound arguments like cos(x1+x2) resolve to [x1, x2]
%   (2 columns) while simple arguments like cos(x5) resolve to [x5] (1 column).
%   Mixing these causes a vertcat dimension mismatch.  Using x5 makes all
%   trig arguments single-variable → consistent shapes → vertcat succeeds.
%   Consistency constraint  x5 = x1+x2  is enforced through the ODE:
%     dx5/dt = dx1/dt + dx2/dt = x3 + x4.
% ─────────────────────────────────────────────────────────────────────────────

syms x1 x2 x3 x4 x5 y1 y2;
x_vars_sym = [x1; x2; x3; x4; x5]; % augmented 5-D state
y_vars_sym = [y1; y2];

% ── Physical parameters ───────────────────────────────────────────────────────
m1 = 1.0; m2 = 1.0; % link masses [kg]
l1 = 4.0; l2 = 4.0; % link lengths [m]
lc1 = 2.0; lc2 = 2.0; % centre-of-mass distances [m]
I1 = 0.02; I2 = 0.02; % moments of inertia [kg·m²]
g = 9.81; % gravitational acceleration [m/s²]

% ── Mass (Inertia) Matrix M(q2) — depends only on x2, no compound argument ───
M11 = I1 + I2 + m1 * lc1 ^ 2 + m2 * (l1 ^ 2 + lc2 ^ 2 + 2 * l1 * lc2 * cos(x2));
M12 = m2 * (lc2 ^ 2 + l1 * lc2 * cos(x2)) + I2;
M22 = m2 * lc2 ^ 2 + I2;
M = [M11, M12; M12, M22];

% ── Coriolis and Centrifugal Matrix C(q2, qdot) — depends only on x2 ─────────
C11 = -m2 * l1 * lc2 * sin(x2) * x4;
C12 = -m2 * l1 * lc2 * sin(x2) * (x3 + x4);
C21 = m2 * l1 * lc2 * sin(x2) * x3;
C_mat = [C11, C12; C21, 0];

% ── Gravity Vector G(q1, q1+q2) — uses x1 and x5 (both single-var args) ──────
G1 = (m1 * g * lc1 + m2 * g * l1) * cos(x1) + m2 * g * lc2 * cos(x5);
G2 = m2 * g * lc2 * cos(x5);
G_vec = [G1; G2];

% ── Drift term f(x) = [dq; M^{-1}(-C*dq - G); dq1+dq2] ──────────────────────
qdot = [x3; x4];
Minv = inv(M);
qddot_f = simplify(Minv * (-C_mat * qdot - G_vec));

fx_sym = [x3; x4; qddot_f; x3 + x4]; % 5×1 (last entry: dx5/dt = x3+x4)

% ── Input matrix g(x) = [0_{2×2}; M^{-1}; 0_{1×2}] ──────────────────────────
gx_sym = [zeros(2, 2); Minv; zeros(1, 2)]; % 5×2

% ── Output map: end-effector position — trig args are x1 and x5 (single-var) ─
hx_sym = [l1 * cos(x1) + l2 * cos(x5);
          l1 * sin(x1) + l2 * sin(x5)]; % 2×1

% ── Sets in output (y) space ──────────────────────────────────────────────────
%   safe_set_sym  >= 0  inside the safe region
%   target_set_sym <= 0  inside the target region (ellipse centred at y=(5.8,1.9))
safe_set_sym =- ((4 * (y1 - 2) - 2 * y2 ^ 3) ^ 2) + 0.8 * y2 ^ 3 + 10;
target_set_sym = ((y1 - 2 - 3.5) ^ 2/1.2 ^ 2) + ((y2 - 1.8) ^ 2/0.4 ^ 2) - 2;

% ── Synthesize reach-avoid backstepping controller (symbolic) ─────────────────
t_design = tic;
[u, k1, J_k1, mu, lambda, certificate, cert_term_dict, A_matrix, b_vector, ks, p, r_deg] = ...
    reach_avoid_controller(fx_sym, gx_sym, hx_sym, x_vars_sym, y_vars_sym, safe_set_sym);
fprintf('__TIMING__,%s,design,%.6f\n', mfilename, toc(t_design));

% ── Optimisation parameters ───────────────────────────────────────────────────
rng(42);

% Torque bounds chosen to cover what the mu=15 reach-avoid controller actually
% needs over the FL region (tau1 in ~[-5016,4402], tau2 in ~[-135,541]); the
% manipulator's drift cancellation (gravity/Coriolis) dominates the magnitude,
% so |tau|<=500 is not achievable while keeping the large reach-avoid set.
lb = [-5500; -700]; % lower torque bounds [N·m]
ub = [5500; 700]; % upper torque bounds [N·m]

ds = 4; % degree of auxiliary SOS polynomials for single-integrator
dv = 2; % degree of k1 controller polynomial

mu_val = 15; % mu — matches Python acrobot notebook (mu_1 = 15)

samples_num = 10000; % ~10 % of samples land in safe set with these bounds

% Sampling bounds for [q1; q2; dq1; dq2; q1+q2]
% q1 ∈ [-2, 0.5]; q2 ∈ [0.8, π-0.8] excludes the elbow singularities q2 ∈ {0, π}
% (det A ∝ sin(q2)), so the decoupling matrix stays nonsingular over the region.
% q1+q2 bounds derived from q1+q2 ∈ [-2+0.8, 0.5+(π-0.8)] = [-1.2, π-0.3].
bound_min = [-2.0; 0.8; -0.5; -0.5; -1.2];
bound_max = [0.5; pi - 0.8; 0.5; 0.5; pi - 0.3];

% ── Solve for bounded control inputs via SOP with SOS constraints ─────────────
t_solve = tic;
[u_opt, certificate_opt, valid_count, k1_opt] = solvesop_bounded_control( ...
    u, k1, J_k1, mu, lambda, certificate, cert_term_dict, p, r_deg, ...
    x_vars_sym, y_vars_sym, hx_sym, safe_set_sym, target_set_sym, ...
    mu_val, lb, ub, ds, dv, samples_num, bound_min, bound_max, ...
    'sop_bounded_control_acrobot_unconstrained.py');
fprintf('__TIMING__,%s,sop_solve,%.6f\n', mfilename, toc(t_solve));

% ── Substitute x5 = x1+x2 back → 4-D expressions for Python export ───────────
% The 5-D state was only needed to satisfy the SOS framework's trig detection.
% Python simulation uses the natural 4-D state [q1; q2; dq1; dq2].
u_opt_4d = subs(u_opt, x5, x1 + x2);
certificate_opt_4d = subs(certificate_opt, x5, x1 + x2);
% k1_opt lives in y-space (y1, y2) — no x5 involved, no substitution needed.
x_vars_4d = [x1; x2; x3; x4];
hx_sym_4d = subs(hx_sym, x5, x1 + x2);
fx_sym_4d = subs(fx_sym(1:4), x5, x1 + x2); % drop 5th row (dx5/dt)
gx_sym_4d = subs(gx_sym(1:4, :), x5, x1 + x2); % drop 5th row (zero anyway)

% ── Estimate controller bounds using 4-D expressions ──────────────────────
bound_min_4d = bound_min(1:4);
bound_max_4d = bound_max(1:4);

[estimated_lb_1, estimated_ub_1] = compute_poly_bounds_sampling( ...
    x_vars_4d, u_opt_4d(1), certificate_opt_4d, 10000, bound_min_4d, bound_max_4d);
[estimated_lb_2, estimated_ub_2] = compute_poly_bounds_sampling( ...
    x_vars_4d, u_opt_4d(2), certificate_opt_4d, 10000, bound_min_4d, bound_max_4d);

disp('------------------------------------------------------------------------------------');
disp('Obtained controller (4-D, x5 substituted) after solving with bounded control inputs:');
disp(u_opt_4d);

disp('Obtained certificate (4-D, x5 substituted) after solving with bounded control inputs:');
disp(certificate_opt_4d);

disp(['Estimated bounds for tau1 over cert. zero-superlevel set: [', ...
          num2str(estimated_lb_1), ', ', num2str(estimated_ub_1), ']']);
disp(['Estimated bounds for tau2 over cert. zero-superlevel set: [', ...
          num2str(estimated_lb_2), ', ', num2str(estimated_ub_2), ']']);

% ── Export 4-D results to Python ────────────────────────────────────────────
params_for_export = struct();
params_for_export.fx_sym = fx_sym_4d;
params_for_export.gx_sym = gx_sym_4d;
params_for_export.hx_sym = hx_sym_4d;
params_for_export.x_vars_sym = x_vars_4d;
params_for_export.y_vars_sym = y_vars_sym;
params_for_export.safe_set_sym = safe_set_sym;
params_for_export.target_set_sym = target_set_sym;
params_for_export.lb = lb;
params_for_export.ub = ub;
params_for_export.ds = ds;
params_for_export.dv = dv;
params_for_export.mu_val = mu_val;
params_for_export.samples_num = samples_num;
params_for_export.valid_count = valid_count;
params_for_export.bound_min = bound_min_4d;
params_for_export.bound_max = bound_max_4d;

export_to_python(u_opt_4d, certificate_opt_4d, k1_opt, params_for_export, ...
'sop_bounded_control_acrobot_result.py');
