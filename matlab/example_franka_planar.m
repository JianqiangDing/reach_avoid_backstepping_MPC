clc;
clear;
close all;

% Planar task-space double integrator (Franka EE in a horizontal plane):
%   dpx/dt = vx
%   dpy/dt = vy
%   dvx/dt = ux   (task-space acceleration, |ux| <= a_max)
%   dvy/dt = uy   (task-space acceleration, |uy| <= a_max)
% state x = [px; py; vx; vy], output y = [px; py].
% The decoupling matrix is identity (no singularity), so the sampling box may
% include zero velocity.
syms x1 x2 x3 x4 y1 y2;
x_vars_sym = [x1; x2; x3; x4];

% drift term f(x)
fx_sym = [x3; x4; 0; 0];

% input matrix g(x), input u = [ux; uy]
gx_sym = [0, 0; 0, 0; 1, 0; 0, 1];

% output map to planar position
hx_sym = [x1; x2];

y_vars_sym = [y1; y2];

% scene-specific symbolic sets + parameters (overwritten per-scene by synth.py):
%   safe_set_sym (>=0 inside safe), target_set_sym (<=0 inside target),
%   a_max, ds, dv, mu_val, samples_num, bound_min, bound_max
franka_planar_scene;

% synthesize the reach-avoid backstepping controller (symbolic)
t_total = tic;
t_design = tic;
[u, k1, J_k1, mu, lambda, certificate, cert_term_dict, A_matrix, b_vector, ks, p, r_deg] = reach_avoid_controller(fx_sym, gx_sym, hx_sym, x_vars_sym, y_vars_sym, safe_set_sym);
fprintf('__TIMING__,%s,design,%.6f\n', mfilename, toc(t_design));

% control input bounds |u_i| <= a_max
rng(42);
lb = [-a_max; -a_max];
ub = [a_max; a_max];

% solve the bounded control inputs using SOP with SOS constraints
t_solve = tic;
[u_opt, certificate_opt, valid_count, k1_opt] = solvesop_bounded_control(u, k1, J_k1, mu, lambda, certificate, cert_term_dict, p, r_deg, x_vars_sym, y_vars_sym, ...
    hx_sym, safe_set_sym, target_set_sym, mu_val, lb, ub, ds, dv, samples_num, bound_min, bound_max, 'sop_bounded_control_franka_planar_unconstrained.py');
fprintf('__TIMING__,%s,sop_solve,%.6f\n', mfilename, toc(t_solve));

disp('------------------------------------------------------------------------------------');
disp('Obtained controller after solving with bounded control inputs:');
disp(u_opt);
disp('Obtained certificate after solving with bounded control inputs:');
disp(certificate_opt);

% export results + settings to a python file for verification and testing
params_for_export = struct();
params_for_export.fx_sym = fx_sym;
params_for_export.gx_sym = gx_sym;
params_for_export.hx_sym = hx_sym;
params_for_export.x_vars_sym = x_vars_sym;
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
params_for_export.bound_min = bound_min;
params_for_export.bound_max = bound_max;

t_export = tic;
export_to_python(u_opt, certificate_opt, k1_opt, params_for_export, 'sop_bounded_control_franka_planar_result.py');
fprintf('__TIMING__,%s,result_export,%.6f\n', mfilename, toc(t_export));
fprintf('__TIMING__,%s,total,%.6f\n', mfilename, toc(t_total));
