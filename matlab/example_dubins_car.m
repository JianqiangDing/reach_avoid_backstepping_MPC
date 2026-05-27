clc;
clear;
close all;

% define the system dynamics using symbolic variables
% Dubins car with speed dynamics:
%   dx1/dt    = v*cos(th)
%   dx2/dt    = v*sin(th)
%   dth/dt    = omega
%   dv/dt     = a
syms x1 x2 th v y1 y2;
x_vars_sym = [x1; x2; th; v];

% drift term f(x)
fx_sym = [v * cos(th); v * sin(th); 0; 0];

% input matrix g(x), input u = [omega; a]
gx_sym = [0, 0; 0, 0; 1, 0; 0, 1];

% output map to planar position
hx_sym = [x1; x2];

y_vars_sym = [y1; y2];

% define the safe set h(y) > 0 and target set g(y) < 0 in y-space
h_raw =- (y1 ^ 4 + y2 ^ 4 - 16) * (y1 ^ 4 + y2 ^ 4 - 4);
target_set_sym = (y2 - 0) ^ 2 + ((y1 + 1.7) / 0.5) ^ 2 - 0.4; % g(y): target_set <= 0 inside target
alpha = 1e-3 * (-target_set_sym + 300);
safe_set_sym = alpha * h_raw; % h(y): zsafe_set >= 0 inside safe set

% synthesize the reach-avoid backstepping controller with bounded control inputs (symbolic)
[u, k1, J_k1, mu, lambda, certificate, cert_term_dict, A_matrix, b_vector, ks, p, r_deg] = reach_avoid_controller(fx_sym, gx_sym, hx_sym, x_vars_sym, y_vars_sym, safe_set_sym);
% here u is a symbolic vector [u1; u2], each entry is a symbolic expression of state variables x1,x2,x3,x4, and unknown parameters
% ers (mu values, lambda, k1 controler polynomial)

% define constraint for the control input bounds
% Ax = [1 0; -1 0; 0 1; 0 -1]; % example constraint matrix for control input bounds, here we want to enforce |u1| <= 100 and |u2| <= 100
% set the random seed for reproducibility
rng(42);
lb = [-5; -5]; % lower bounds for control inputs [omega_min; a_min]
ub = [5; 5]; % upper bounds for control inputs [omega_max; a_max]
ds = 4; % degree of the auxiliary SOS polynomials for the single-integrator system
dv = 4; % degree of the k1 controller polynomial

mu_val = 0.1; % example value for mu just for testing, HAVE TO DISCUSS THIS IN THE PAPER

samples_num = 1000; % number of random samples to find the valid samples that satisfy the control input bounds for the pseudo ux

bound_min = [-2; -2; 2 * pi / 3; -1.0]; % lower bounds for sampling state [x1; x2; th; v]
bound_max = [2; 2; 4 * pi / 3; 1.0]; % upper bounds for sampling state [x1; x2; th; v]

% solve the bounded control inputs using scenario optimization programming (SOP) with SOS constraints
[u_opt, certificate_opt, valid_count, k1_opt] = solvesop_bounded_control(u, k1, J_k1, mu, lambda, certificate, cert_term_dict, p, r_deg, x_vars_sym, y_vars_sym, ...
    hx_sym, safe_set_sym, target_set_sym, mu_val, lb, ub, ds, dv, samples_num, bound_min, bound_max);

% compute the bounds of the obtained controller over zero superlevel set of the certificate
% [num_1, den_1] = numden(u_opt(1));
% [lb_1, ub_1] = compute_poly_bounds_sos(num_1, den_1, certificate_opt, ds, 1e-3);
% estimate the bounds of the obtained controller over zero superlevel set of the certificate using sampling (for verification)
[estimated_lb_1, estimated_ub_1] = compute_poly_bounds_sampling(x_vars_sym, u_opt(1), certificate_opt, 10000, bound_min, bound_max);
[estimated_lb_2, estimated_ub_2] = compute_poly_bounds_sampling(x_vars_sym, u_opt(2), certificate_opt, 10000, bound_min, bound_max);

disp('------------------------------------------------------------------------------------');

% disp the obtained controller after solving with bounded control inputs
disp('Obtained controller after solving with bounded control inputs:');
disp(u_opt);

% disp the obtained certificate after solving with bounded control inputs
disp('Obtained certificate after solving with bounded control inputs:');
disp(certificate_opt);

% % disp the number of valid samples that satisfy the control input bounds for the pseudo ux
% disp(['Number of valid samples that satisfy the control input bounds for the pseudo ux: ', num2str(valid_count), ' out of ', num2str(samples_num)]);

% disp(['Bounds for u1 over the zero superlevel set of the certificate: [', ...
%           num2str(lb_1), ', ', num2str(ub_1), ']']);

disp(['Estimated bounds for u1 over the zero superlevel set of the certificate using sampling: [', ...
          num2str(estimated_lb_1), ', ', num2str(estimated_ub_1), ']']);
disp(['Estimated bounds for u2 over the zero superlevel set of the certificate using sampling: [', ...
          num2str(estimated_lb_2), ', ', num2str(estimated_ub_2), ']']);

% export the results, also all setting parameters, to a python file for verification and testing
% first, use a struct to collect all the relevant parameters for exporting
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

export_to_python(u_opt, certificate_opt, k1_opt, params_for_export, 'sop_bounded_control_dubins_car_result.py');
% export the computed controller and certificate to a python file for verification and testing
