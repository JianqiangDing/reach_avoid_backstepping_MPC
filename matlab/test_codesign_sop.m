% Does raising xi0 (-> larger lambda -> SOP k1 leverage) let the SAMPLE
% constraints pull u toward the bound? For a grid of xi0 we run the whole-region
% slack SOP (solvesop_bounded_control_slack_xi) and report, per control channel:
%   - distinguishability:  max|u_con - u_uncon| / max|u_uncon|  over the reach-avoid region
%   - dense |u| range of the CONSTRAINED controller (the true achievable bound)
%   - dense |u| range of the UNCONSTRAINED controller (same lambda)
%   - slack-achievable bound = ub + slack (what the SOP claims at the samples)
%   - reach-avoid set size (fraction of the box with cert_con >= 0) -- validity check
%
% The requested bound is just a reference; the slack + dense range reveal the
% improvement limit. Writes data/codesign_sop_<sys>.csv.

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
[ud, k1d, Jk1d, mud, lamd, certd, ctd, ~, ~, ~, pd, rd] = ...
    reach_avoid_controller(fx_d, gx_d, hx_d, x_vars_d, y_vars_d, safe_d);
run_codesign('DUBINS', ud, k1d, Jk1d, mud, lamd, certd, ctd, pd, rd, ...
    x_vars_d, y_vars_d, hx_d, safe_d, target_d, 0.1, [-5; -5], [5; 5], 4, 4, 1000, ...
    [-2; -2; 2*pi/3; 0.1], [2; 2; 4*pi/3; 1.0], [1e-8, 0.1, 1, 10], ...
    fullfile(data_dir, 'codesign_sop_dubins_car.csv'));

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
[um, k1m, Jk1m, mum, lamm, certm, ctm, ~, ~, ~, pm, rm] = ...
    reach_avoid_controller(fx_m, gx_m, hx_m, x_vars_m, y_vars_m, safe_m);
run_codesign('MANIPULATOR', um, k1m, Jk1m, mum, lamm, certm, ctm, pm, rm, ...
    x_vars_m, y_vars_m, hx_m, safe_m, target_m, 15, [-500; -500], [500; 500], 4, 2, 10000, ...
    [-2; 0.8; -0.5; -0.5; -1.2], [0.5; pi - 0.8; 0.5; 0.5; pi - 0.3], [1e-8, 0.1, 1, 15], ...
    fullfile(data_dir, 'codesign_sop_manipulator.csv'));

% ------------------------------------------------------------------------
function run_codesign(name, u, k1, J_k1, mu, lambda, certificate, cert_term_dict, p, r_deg, ...
        x_vars, y_vars, hx, safe_set, target_set, mu_val, lb, ub, ds, dv, samples_num, ...
        bound_min, bound_max, xi0_grid, csv_path)

    n = numel(x_vars); m = numel(u);
    N = 20000;
    bmin = bound_min(:)'; bmax = bound_max(:)';

    fprintf('\n================= %s : co-design + slack SOP (mu=%.4g, req ub=%.4g) =================\n', ...
        name, mu_val, ub(1));

    rows = [];
    for xi0 = xi0_grid
        rng(42);
        [ux_con, ux_uncon, cert_con, cert_uncon, ~, ~, slack, lam] = ...
            solvesop_bounded_control_slack_xi(u, k1, J_k1, mu, lambda, certificate, ...
            cert_term_dict, p, r_deg, x_vars, y_vars, hx, safe_set, target_set, ...
            mu_val, lb, ub, ds, dv, samples_num, bound_min, bound_max, xi0);

        % box samples for distinguishability + reach-avoid set size
        samples = rand(N, n) .* (bmax - bmin) + bmin;
        cols = num2cell(samples, 1);
        cu = ev(cert_uncon, x_vars, cols, N) >= 0;     % reach-avoid region (uncon cert)
        cc = ev(cert_con,   x_vars, cols, N) >= 0;
        raset_frac = mean(cc);

        fprintf('  xi0=%-8.3g lam=%-10.4g raset=%.3f   slack=[%.4g %.4g]\n', ...
            xi0, lam, raset_frac, slack(1), slack(2));
        for i = 1:m
            [clo, chi] = compute_poly_bounds_sampling(x_vars, ux_con(i),   cert_con,   10000, bmin, bmax);
            [ulo, uhi] = compute_poly_bounds_sampling(x_vars, ux_uncon(i), cert_uncon, 10000, bmin, bmax);
            if isempty(clo); clo = NaN; chi = NaN; end   % reach-avoid set empty (cert never >= 0)
            if isempty(ulo); ulo = NaN; uhi = NaN; end
            % distinguishability over the (uncon) reach-avoid region
            ufi = ev(ux_con(i),   x_vars, cols, N);
            uui = ev(ux_uncon(i), x_vars, cols, N);
            reg = cu;
            if nnz(reg) < 50; reg = true(N, 1); end
            dist = 100 * max(abs(ufi(reg) - uui(reg))) / max(max(abs(uui(reg))), eps);
            achiev = ub(i) + slack(i);
            fprintf('    u%d: con[% .3g, % .3g] unc[% .3g, % .3g] achiev=%.4g dist=%.3g%%\n', ...
                i, clo, chi, ulo, uhi, achiev, dist);
            rows = [rows; xi0, lam, i, slack(i), achiev, clo, chi, ulo, uhi, dist, raset_frac]; %#ok<AGROW>
        end
    end

    T = array2table(rows, 'VariableNames', {'xi0', 'lambda', 'dim', 'slack', 'achievable', ...
        'con_lo', 'con_hi', 'unc_lo', 'unc_hi', 'dist_pct', 'raset_frac'});
    T.ub_req = repmat(ub(1), height(T), 1);
    T.mu_val = repmat(mu_val, height(T), 1);
    writetable(T, csv_path);
    fprintf('  wrote %s\n', csv_path);
end

function vals = ev(expr, x_vars, cols, N)
    f = matlabFunction(expr, 'vars', x_vars);
    vals = f(cols{:}) .* ones(N, 1);
end
