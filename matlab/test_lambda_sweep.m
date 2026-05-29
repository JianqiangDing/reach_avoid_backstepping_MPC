% Sweep the backstepping scale lambda and measure how much the (constructed) k1
% then shapes the control u.  Mirrors the user's hypothesis: the constructed k1
% barely acts because lambda is pinned at ~1e-8; if lambda were raised toward mu's
% magnitude, would k1 matter?
%
% lambda enters u only via 0.5*lambda*(Lf h_i - k1_i). We hold the obtained
% (vanilla) k1 FIXED and substitute a grid of lambda values, comparing:
%     u_full(lambda) : k1 = obtained
%     u_k0(lambda)   : k1 = 0  (value + Jacobian)
% k1's effect at a given lambda = max|u_full - u_k0| / max|u_full|  (over the box).
%
% NOTE (validity): the obtained k1 was synthesized at the tiny solved lambda, so
% the large-lambda points are a *what-if sensitivity* of the SAME k1, not a re-
% designed controller. That is exactly the user's question ("raise lambda, does
% the constructed k1 act more?"). A proper co-design would re-solve k1 per lambda.
%
% Writes data/lambda_sweep_<sys>.csv for scripts/compare_lambda_effect.ipynb.

clc; clear; close all;
script_dir = fileparts(mfilename('fullpath'));
addpath(script_dir);
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
sweep_lambda('DUBINS', ud, k1d, Jk1d, mud, lamd, x_vars_d, y_vars_d, hx_d, ...
    safe_d, target_d, 4, 4, 0.1, [-2; -2; 2*pi/3; 0.1], [2; 2; 4*pi/3; 1.0], ...
    fullfile(data_dir, 'lambda_sweep_dubins_car.csv'));

% ========================= MANIPULATOR ==================================
clearvars -except data_dir
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
sweep_lambda('MANIPULATOR', um, k1m, Jk1m, mum, lamm, x_vars_m, y_vars_m, hx_m, ...
    safe_m, target_m, 4, 2, 15, [-2; 0.8; -0.5; -0.5; -1.2], [0.5; pi - 0.8; 0.5; 0.5; pi - 0.3], ...
    fullfile(data_dir, 'lambda_sweep_manipulator.csv'));

% ------------------------------------------------------------------------
function sweep_lambda(name, u, k1, J_k1, mu, lambda, x_vars, y_vars, hx, ...
        safe_set, target_set, ds, dv, mu_val, bound_min, bound_max, csv_path)

    rng(42);
    [k1_y, k1_lambda, ~] = solve_vanilla_k1_controller(y_vars, safe_set, target_set, dv, ds);
    J_k1_y = jacobian(k1_y, y_vars);
    p = numel(k1);

    % u with the obtained k1, lambda left FREE (mu substituted, in x-space)
    u_full = subs(u, k1, k1_y);
    u_full = subs(u_full, J_k1, J_k1_y);
    u_full = subs(u_full, y_vars, hx);
    u_full = sub_mu(u_full, mu, mu_val);

    % u with k1 = 0 (value + Jacobian), lambda free
    u_k0 = subs(u, k1, zeros(p, 1));
    u_k0 = subs(u_k0, J_k1, zeros(p, p));
    u_k0 = subs(u_k0, y_vars, hx);
    u_k0 = sub_mu(u_k0, mu, mu_val);

    % lambda grid: 1e-8 .. 1e3, plus the solved lambda and mu_val themselves
    lam_grid = unique([logspace(-8, 3, 34), k1_lambda, mu_val]);

    n = numel(x_vars); m = numel(u_full);
    N = 20000;
    bmin = bound_min(:)'; bmax = bound_max(:)';
    samples = rand(N, n) .* (bmax - bmin) + bmin;

    % vectorized handles in [x_vars; lambda]
    ff = cell(m, 1); fk = cell(m, 1);
    for i = 1:m
        ff{i} = matlabFunction(u_full(i), 'vars', [x_vars(:); lambda]);
        fk{i} = matlabFunction(u_k0(i),  'vars', [x_vars(:); lambda]);
    end

    fprintf('\n================= %s : lambda sweep (mu=%.4g, solved lambda=%.3g) =================\n', ...
        name, mu_val, k1_lambda);
    fprintf('%-12s', 'lambda');
    for i = 1:m; fprintf('  k1_effect(u%d) %%', i); end
    fprintf('\n');

    rows = [];  % lambda, dim, k1_effect_pct, u_scale
    for li = 1:numel(lam_grid)
        lam = lam_grid(li);
        args = num2cell([samples, lam * ones(N, 1)], 1);
        fprintf('%-12.3g', lam);
        for i = 1:m
            uf = ff{i}(args{:}) .* ones(N, 1);
            uk = fk{i}(args{:}) .* ones(N, 1);
            scale = max(abs(uf));
            eff = 100 * max(abs(uf - uk)) / scale;
            rows = [rows; lam, i, eff, scale]; %#ok<AGROW>
            fprintf('   %12.4g', eff);
        end
        fprintf('\n');
    end

    % crossing lambdas (worst control dim): where does k1's effect reach 1/10/50%?
    for thr = [1, 10, 50]
        lam_cross = NaN;
        for li = 1:numel(lam_grid)
            effs = rows(rows(:, 1) == lam_grid(li), 3);
            if max(effs) >= thr; lam_cross = lam_grid(li); break; end
        end
        fprintf('  k1 effect reaches %2d%% of |u| at lambda >= %.3g   (mu = %.3g)\n', ...
            thr, lam_cross, mu_val);
    end

    T = array2table(rows, 'VariableNames', {'lambda', 'dim', 'k1_effect_pct', 'u_scale'});
    T.mu_val = repmat(mu_val, height(T), 1);
    T.k1_lambda = repmat(k1_lambda, height(T), 1);
    writetable(T, csv_path);
    fprintf('  wrote %s\n', csv_path);
end

function expr = sub_mu(expr, mu, mu_val)
    for i = 1:numel(mu)
        if ~isempty(mu{i})
            expr = subs(expr, mu{i}, mu_val);
        end
    end
end
