% Proper co-design test: re-solve the vanilla k1 with the SOS floor xi0 raised
% toward mu's magnitude (lambda stays a decision variable, floored at xi0 -- NOT
% fixed). For each xi0 we report:
%   - feasibility of the SOS solve (Mosek pinf/dinf),
%   - the solved lambda and delta,
%   - k1's effect on the resulting u  ( max|u(k1) - u(k1=0)| / max|u| over the box ),
%   - the control magnitude max|u|.
%
% Question: does a mu-scale lambda (a) stay feasible and (b) make the co-designed
% k1 actually shape u? Writes data/lambda_codesign_<sys>.csv.

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
mu_d = 0.1;
codesign('DUBINS', ud, k1d, Jk1d, mud, lamd, x_vars_d, y_vars_d, hx_d, ...
    safe_d, target_d, 4, 4, mu_d, [-2; -2; 2*pi/3; 0.1], [2; 2; 4*pi/3; 1.0], ...
    [1e-8, 1e-3, 1e-2, mu_d, 1, 10*mu_d, 10, 100], ...
    fullfile(data_dir, 'lambda_codesign_dubins_car.csv'));

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
mu_m = 15;
codesign('MANIPULATOR', um, k1m, Jk1m, mum, lamm, x_vars_m, y_vars_m, hx_m, ...
    safe_m, target_m, 4, 2, mu_m, [-2; 0.8; -0.5; -0.5; -1.2], [0.5; pi - 0.8; 0.5; 0.5; pi - 0.3], ...
    [1e-8, 1e-2, 0.1, 1, mu_m, 10*mu_m, 100*mu_m], ...
    fullfile(data_dir, 'lambda_codesign_manipulator.csv'));

% ------------------------------------------------------------------------
function codesign(name, u, k1, J_k1, mu, lambda, x_vars, y_vars, hx, ...
        safe_set, target_set, ds, dv, mu_val, bound_min, bound_max, xi0_grid, csv_path)

    rng(42);
    p = numel(k1);
    n = numel(x_vars); m = numel(u);
    N = 20000;
    bmin = bound_min(:)'; bmax = bound_max(:)';
    samples = rand(N, n) .* (bmax - bmin) + bmin;
    cols = num2cell(samples, 1);

    fprintf('\n================= %s : co-design (raise xi0 toward mu=%.4g) =================\n', name, mu_val);
    fprintf('%-10s %-6s %-12s %-12s %-12s %-12s\n', ...
        'xi0', 'feas', 'solved_lam', 'k1_eff_%', 'max|u|', 'delta');

    rows = [];
    for xi0 = xi0_grid
        ok = true; lam = NaN; del = NaN; eff_worst = NaN; umax = NaN; feasflag = 0;
        try
            [k1_y, lam, del, info] = solve_vanilla_k1_controller_xi(y_vars, safe_set, target_set, dv, ds, xi0);
            feasflag = double(isfinite(lam) && (isnan(info.pinf) || info.pinf == 0) && ...
                (isnan(info.dinf) || info.dinf == 0));
        catch ME
            ok = false;
            fprintf('%-10.3g %-6s  (solve error: %s)\n', xi0, 'ERR', ME.message);
        end

        if ok && isfinite(lam)
            J_k1_y = jacobian(k1_y, y_vars);
            u_full = sub_mu(subs(subs(subs(u, k1, k1_y), J_k1, J_k1_y), y_vars, hx), mu, mu_val);
            u_full = subs(u_full, lambda, lam);
            u_k0 = sub_mu(subs(subs(subs(u, k1, zeros(p, 1)), J_k1, zeros(p, p)), y_vars, hx), mu, mu_val);
            u_k0 = subs(u_k0, lambda, lam);

            eff_worst = 0; umax = 0;
            for i = 1:m
                uf = ev(u_full(i), x_vars, cols, N);
                uk = ev(u_k0(i),  x_vars, cols, N);
                scale = max(abs(uf));
                eff_worst = max(eff_worst, 100 * max(abs(uf - uk)) / scale);
                umax = max(umax, scale);
            end
            fprintf('%-10.3g %-6d %-12.4g %-12.4g %-12.4g %-12.4g\n', ...
                xi0, feasflag, lam, eff_worst, umax, del);
            rows = [rows; xi0, feasflag, lam, eff_worst, umax, del]; %#ok<AGROW>
        elseif ok
            fprintf('%-10.3g %-6s  (no finite lambda returned)\n', xi0, 'INF');
            rows = [rows; xi0, 0, NaN, NaN, NaN, NaN]; %#ok<AGROW>
        end
    end

    if ~isempty(rows)
        T = array2table(rows, 'VariableNames', ...
            {'xi0', 'feasible', 'solved_lambda', 'k1_effect_pct', 'max_u', 'delta'});
        T.mu_val = repmat(mu_val, height(T), 1);
        writetable(T, csv_path);
        fprintf('  wrote %s\n', csv_path);
    end
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
    vals = f(cols{:}) .* ones(N, 1);
end
