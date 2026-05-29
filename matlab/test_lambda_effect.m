% Direct test of the hypothesis: the synthesized k1 contributes little to the
% final control u because the backstepping scale lambda is too small.
%
% lambda enters u ONLY through the term  0.5*lambda*(Lf h_i - k1_i)  in b(x)
% (u = A(x)^{-1} b(x)). This is also the only place k1's *value* (as opposed to
% its Jacobian) is weighted. So:
%
%   - PRIMARY test (what the user asked): with the obtained k1, build u twice,
%       u_full : lambda = computed (k1_lambda from solve_vanilla_k1_controller)
%       u_zero : lambda = 0
%     and compare them over the WHOLE FL box, every state dimension varying
%     (no slicing / no fixed dimensions). |u_full - u_zero| is exactly the
%     contribution of the 0.5*lambda*(Lf h - k1) term to u.
%
%   - SUPPLEMENTARY: u_k0 with k1 = 0 (and J_k1 = 0), lambda = computed, to show
%     k1's TOTAL effect on u (k1 also enters lambda-independently via the
%     Jacobian term sum_j J_k1(i,j) Lf h_j).
%
% No certificate / no reach-avoid sampling is needed here: we evaluate the
% closed-form control laws directly on the full feedback-linearizable box.

clc; clear; close all;
script_dir = fileparts(mfilename('fullpath'));
addpath(script_dir);  % matlab/
data_dir = fullfile(fileparts(script_dir), 'data');
if ~exist(data_dir, 'dir'); mkdir(data_dir); end

% ============================ DUBINS ====================================
syms x1 x2 th v y1 y2;
x_vars_d = [x1; x2; th; v];
fx_d = [v * cos(th); v * sin(th); 0; 0];
gx_d = [0, 0; 0, 0; 1, 0; 0, 1];
hx_d = [x1; x2];
y_vars_d = [y1; y2];
h_raw = -(y1^4 + y2^4 - 16) * (y1^4 + y2^4 - 4);
target_d = (y2 - 0)^2 + ((y1 + 1.7) / 0.5)^2 - 0.4;
alpha = 1e-3 * (-target_d + 300);
safe_d = alpha * h_raw;
[ud, k1d, Jk1d, mud, lamd, ~, ~, ~, ~, ~, ~, ~] = ...
    reach_avoid_controller(fx_d, gx_d, hx_d, x_vars_d, y_vars_d, safe_d);
compare_lambda('DUBINS', ud, k1d, Jk1d, mud, lamd, x_vars_d, y_vars_d, hx_d, ...
    safe_d, target_d, 4, 4, 0.1, [-2; -2; 2*pi/3; 0.1], [2; 2; 4*pi/3; 1.0], ...
    fullfile(data_dir, 'lambda_effect_dubins_car.csv'));

% ========================= MANIPULATOR ==================================
clearvars -except data_dir  % keep workspace tidy (but keep the output dir)
syms x1 x2 x3 x4 x5 y1 y2;
x_vars_m = [x1; x2; x3; x4; x5];
y_vars_m = [y1; y2];
m1 = 1.0; m2 = 1.0; l1 = 4.0; l2 = 4.0; lc1 = 2.0; lc2 = 2.0; I1 = 0.02; I2 = 0.02; g = 9.81;
M11 = I1 + I2 + m1*lc1^2 + m2*(l1^2 + lc2^2 + 2*l1*lc2*cos(x2));
M12 = m2*(lc2^2 + l1*lc2*cos(x2)) + I2;
M22 = m2*lc2^2 + I2;
M = [M11, M12; M12, M22];
C_mat = [-m2*l1*lc2*sin(x2)*x4, -m2*l1*lc2*sin(x2)*(x3 + x4); ...
    m2*l1*lc2*sin(x2)*x3, 0];
G_vec = [(m1*g*lc1 + m2*g*l1)*cos(x1) + m2*g*lc2*cos(x5); m2*g*lc2*cos(x5)];
Minv = inv(M);
qddot_f = simplify(Minv * (-C_mat * [x3; x4] - G_vec));
fx_m = [x3; x4; qddot_f; x3 + x4];
gx_m = [zeros(2, 2); Minv; zeros(1, 2)];
hx_m = [l1*cos(x1) + l2*cos(x5); l1*sin(x1) + l2*sin(x5)];
safe_m = -((4*(y1 - 2) - 2*y2^3)^2) + 0.8*y2^3 + 10;
target_m = ((y1 - 2 - 3.5)^2 / 1.2^2) + ((y2 - 1.8)^2 / 0.4^2) - 2;
[um, k1m, Jk1m, mum, lamm, ~, ~, ~, ~, ~, ~, ~] = ...
    reach_avoid_controller(fx_m, gx_m, hx_m, x_vars_m, y_vars_m, safe_m);
compare_lambda('MANIPULATOR', um, k1m, Jk1m, mum, lamm, x_vars_m, y_vars_m, hx_m, ...
    safe_m, target_m, 4, 2, 15, [-2; 0.8; -0.5; -0.5; -1.2], [0.5; pi - 0.8; 0.5; 0.5; pi - 0.3], ...
    fullfile(data_dir, 'lambda_effect_manipulator.csv'));

% ------------------------------------------------------------------------
function compare_lambda(name, u, k1, J_k1, mu, lambda, x_vars, y_vars, hx, ...
        safe_set, target_set, ds, dv, mu_val, bound_min, bound_max, csv_path)

    rng(42);
    xi0 = 1e-8;  % the SOS floor lambda >= xi0 in solve_vanilla_k1_controller
    % obtained (vanilla) k1 and its computed lambda
    [k1_y, k1_lambda, ~] = solve_vanilla_k1_controller(y_vars, safe_set, target_set, dv, ds);
    J_k1_y = jacobian(k1_y, y_vars);

    % base control law with the obtained k1, in x-space, mu substituted, lambda still free
    u_base = subs(u, k1, k1_y);
    u_base = subs(u_base, J_k1, J_k1_y);
    u_base = subs(u_base, y_vars, hx);
    u_base = sub_mu(u_base, mu, mu_val);

    u_full = subs(u_base, lambda, k1_lambda);   % computed lambda
    u_zero = subs(u_base, lambda, 0);           % lambda = 0

    % supplementary: k1 = 0 (value AND Jacobian), lambda = computed -> total k1 effect
    p = numel(k1);
    u_k0 = subs(u, k1, zeros(p, 1));
    u_k0 = subs(u_k0, J_k1, zeros(p, p));
    u_k0 = subs(u_k0, y_vars, hx);
    u_k0 = sub_mu(u_k0, mu, mu_val);
    u_k0 = subs(u_k0, lambda, k1_lambda);

    % sample the WHOLE feedback-linearizable box (every dimension varies)
    n = numel(x_vars);
    m = numel(u_full);
    N = 20000;
    bmin = bound_min(:)'; bmax = bound_max(:)';
    samples = rand(N, n) .* (bmax - bmin) + bmin;
    cols = num2cell(samples, 1);

    fprintf('\n================= %s : does lambda (hence k1) matter for u? =================\n', name);
    fprintf('computed lambda (k1_lambda) = %.6g    (mu = %.4g)\n', k1_lambda, mu_val);
    fprintf('full-box samples: %d   (all %d state dims varied, no slicing)\n', N, n);

    max_rel_lambda = 0; max_rel_k1 = 0;
    UF = zeros(N, m); UZ = zeros(N, m); UK = zeros(N, m);
    for i = 1:m
        ff = ev(u_full(i), x_vars, cols, N);
        fz = ev(u_zero(i), x_vars, cols, N);
        fk = ev(u_k0(i),  x_vars, cols, N);
        UF(:, i) = ff; UZ(:, i) = fz; UK(:, i) = fk;

        d_lambda = abs(ff - fz);              % effect of the 0.5*lambda*(Lf h - k1) term
        d_k1     = abs(ff - fk);              % total effect of k1 (both pathways)
        scale    = max(abs(ff));              % reference magnitude of u

        rel_lambda = 100 * max(d_lambda) / scale;
        rel_k1     = 100 * max(d_k1) / scale;
        max_rel_lambda = max(max_rel_lambda, rel_lambda);
        max_rel_k1     = max(max_rel_k1, rel_k1);

        fprintf('\n u%d:\n', i);
        fprintf('   range(lambda=computed) = [%.4g, %.4g]\n', min(ff), max(ff));
        fprintf('   range(lambda=0)        = [%.4g, %.4g]\n', min(fz), max(fz));
        fprintf('   range(k1=0)            = [%.4g, %.4g]\n', min(fk), max(fk));
        fprintf('   |u(lam=c) - u(lam=0)| : max=%.4g  mean=%.4g  -> max %.3g%% of |u|   [lambda/k1-value term]\n', ...
            max(d_lambda), mean(d_lambda), rel_lambda);
        fprintf('   |u(lam=c) - u(k1=0)|  : max=%.4g  mean=%.4g  -> max %.3g%% of |u|   [TOTAL k1 effect]\n', ...
            max(d_k1), mean(d_k1), rel_k1);
    end

    fprintf('\n VERDICT (%s):\n', name);
    fprintf('   lambda pathway : u(lambda=computed) vs u(lambda=0) differ by at most %.3g%% of |u|\n', max_rel_lambda);
    fprintf('   total k1       : u(lambda=computed) vs u(k1=0)     differ by at most %.3g%% of |u|\n', max_rel_k1);
    if max_rel_lambda < 1
        fprintf('   => lambda is effectively negligible: the 0.5*lambda*(Lf h - k1) term barely moves u.\n');
    end
    if max_rel_k1 < 1
        fprintf('   => k1 has essentially NO effect on u (both pathways): u is set by drift cancellation + mu*grad.\n');
    elseif max_rel_lambda < 1
        fprintf('   => k1 still affects u, but via the lambda-independent Jacobian term, NOT via lambda.\n');
    end

    % ---- dump per-sample data for the companion notebook -------------------
    state_names = arrayfun(@(s) char(s), x_vars(:), 'UniformOutput', false)';
    u_names = @(suf) arrayfun(@(i) sprintf('u%d_%s', i, suf), 1:m, 'UniformOutput', false);
    col_names = [state_names, u_names('full'), u_names('zero'), u_names('k0'), ...
        {'k1_lambda', 'mu_val', 'xi0'}];
    const_cols = repmat([k1_lambda, mu_val, xi0], N, 1);
    T = array2table([samples, UF, UZ, UK, const_cols], 'VariableNames', col_names);
    writetable(T, csv_path);
    fprintf('   wrote %s\n', csv_path);
end

function expr = sub_mu(expr, mu, mu_val)
    for i = 1:numel(mu)
        if ~isempty(mu{i})
            expr = subs(expr, mu{i}, mu_val);
        end
    end
end

function vals = ev(expr, x_vars, cols, N)
    f = matlabFunction(expr, 'vars', x_vars);
    out = f(cols{:});
    vals = out .* ones(N, 1);   % broadcast in case expr is constant in some/all vars
end
