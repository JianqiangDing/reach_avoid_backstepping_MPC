% Manipulator analog of sweep_mu_dubins.m (§4): sweep mu and report, per mu,
%   - reach-avoid set size as a fraction of the safe set (ra/safe)
%   - DENSE |u| range over the certified reach-avoid set
%   - within the requested bound +-500?
%
% Method: solve vanilla k1 once at the baseline xi0=1e-8 (mu-independent),
% then substitute mu as a scalar parameter into both u(x;mu) and V(x;mu) and
% evaluate over a single box sampling. (mu-independent solve + mu-vectorized eval
% keeps the run cheap — one matlabFunction compilation, fast scan over the grid.)
%
% Writes data/mu_sweep_manipulator.csv (consumed by REVISION_NOTES §8 / notebook).

clc; clear; close all;
script_dir = fileparts(mfilename('fullpath'));
addpath(script_dir);
data_dir = fullfile(fileparts(script_dir), 'data');
if ~exist(data_dir, 'dir'); mkdir(data_dir); end

% ---- manipulator system (matches example_manipulator.m, FL-restricted region) ----
syms x1 x2 x3 x4 x5 y1 y2;
x_vars = [x1; x2; x3; x4; x5];
y_vars = [y1; y2];
m1 = 1.0; m2 = 1.0; l1 = 4.0; l2 = 4.0; lc1 = 2.0; lc2 = 2.0; I1 = 0.02; I2 = 0.02; g = 9.81;
M11 = I1 + I2 + m1*lc1^2 + m2*(l1^2 + lc2^2 + 2*l1*lc2*cos(x2));
M12 = m2*(lc2^2 + l1*lc2*cos(x2)) + I2;
M22 = m2*lc2^2 + I2;
M = [M11, M12; M12, M22];
C_mat = [-m2*l1*lc2*sin(x2)*x4, -m2*l1*lc2*sin(x2)*(x3 + x4); m2*l1*lc2*sin(x2)*x3, 0];
G_vec = [(m1*g*lc1 + m2*g*l1)*cos(x1) + m2*g*lc2*cos(x5); m2*g*lc2*cos(x5)];
Minv = inv(M);
qddot_f = simplify(Minv * (-C_mat * [x3; x4] - G_vec));
fx = [x3; x4; qddot_f; x3 + x4];
gx = [zeros(2, 2); Minv; zeros(1, 2)];
hx = [l1*cos(x1) + l2*cos(x5); l1*sin(x1) + l2*sin(x5)];
safe_set = -((4*(y1 - 2) - 2*y2^3)^2) + 0.8*y2^3 + 10;
target_set = ((y1 - 2 - 3.5)^2 / 1.2^2) + ((y2 - 1.8)^2 / 0.4^2) - 2;

[u, k1, J_k1, mu, lambda, certificate, ~, ~, ~, ~, ~, ~] = ...
    reach_avoid_controller(fx, gx, hx, x_vars, y_vars, safe_set);

% baseline vanilla k1 (xi0=1e-8; mu-independent solve)
ds = 4; dv = 2;
[k1_y, k1_lambda, ~] = solve_vanilla_k1_controller_xi(y_vars, safe_set, target_set, dv, ds, 1e-8);
J_k1_y = jacobian(k1_y, y_vars);

% bind k1, J_k1, y, lambda; keep mu symbolic (replace per-output mu vectors with
% a single shared scalar mu_param so eval can scan a grid without re-compiling)
syms mu_param real;
u_p = subs(u, k1, k1_y);
u_p = subs(u_p, J_k1, J_k1_y);
u_p = subs(u_p, y_vars, hx);
u_p = subs(u_p, lambda, k1_lambda);

cert_p = subs(certificate, k1, k1_y);
cert_p = subs(cert_p, y_vars, hx);
cert_p = subs(cert_p, lambda, k1_lambda);

for i = 1:numel(mu)
    if ~isempty(mu{i})
        u_p = subs(u_p, mu{i}, mu_param);
        cert_p = subs(cert_p, mu{i}, mu_param);
    end
end

% compile once with mu_param as an extra argument
fprintf('compiling matlabFunction(u_p(1), u_p(2), cert_p) — large symbolic expressions...\n');
t = tic;
uf1 = matlabFunction(u_p(1), 'vars', [x_vars; mu_param]);
uf2 = matlabFunction(u_p(2), 'vars', [x_vars; mu_param]);
cf  = matlabFunction(cert_p,  'vars', [x_vars; mu_param]);
safe_func = matlabFunction(subs(safe_set, y_vars, hx), 'vars', x_vars);
fprintf('compile time: %.1fs\n', toc(t));

% one box sampling (FL-restricted region, matches the other manip tests)
N = 30000;
bmin = [-2; 0.8; -0.5; -0.5; -1.2];
bmax = [0.5; pi - 0.8; 0.5; 0.5; pi - 0.3];
rng(42);
samples = bmin' + (bmax - bmin)' .* rand(N, 5);
cols = num2cell(samples, 1);
sv = safe_func(cols{:}) .* ones(N, 1);
safe_mask = sv >= 0;
safe_cnt = nnz(safe_mask);
fprintf('samples in {safe>=0}: %d / %d\n', safe_cnt, N);

mu_list = [5, 10, 15, 30, 50, 100, 200, 500, 1000];
ub_req = 500;

fprintf('\n%-6s  %-10s  %-22s  %-22s  %-12s  %-8s\n', ...
    'mu', 'ra/safe', 'u1 range (RA set)', 'u2 range (RA set)', 'max|u|', 'within?');

rows = [];
for mu_val = mu_list
    args = num2cell([samples, mu_val * ones(N, 1)], 1);
    cv = cf(args{:}) .* ones(N, 1);
    ra_mask = (cv >= 0) & safe_mask;
    ra_cnt = nnz(ra_mask);
    ra_over_safe = ra_cnt / max(safe_cnt, 1);
    uv1 = uf1(args{:}) .* ones(N, 1);
    uv2 = uf2(args{:}) .* ones(N, 1);
    if ra_cnt >= 10
        lb1 = min(uv1(ra_mask)); ub1 = max(uv1(ra_mask));
        lb2 = min(uv2(ra_mask)); ub2 = max(uv2(ra_mask));
        max_u = max(max(abs([lb1, ub1])), max(abs([lb2, ub2])));
    else
        lb1 = NaN; ub1 = NaN; lb2 = NaN; ub2 = NaN; max_u = NaN;
    end
    within = double(isfinite(max_u) && max_u <= ub_req);
    fprintf('%-6g  %-10.4f  %-22s  %-22s  %-12.4g  %-8d\n', ...
        mu_val, ra_over_safe, ...
        sprintf('[%.3g, %.3g]', lb1, ub1), ...
        sprintf('[%.3g, %.3g]', lb2, ub2), ...
        max_u, within);
    rows = [rows; mu_val, ra_cnt, safe_cnt, ra_over_safe, lb1, ub1, lb2, ub2, max_u, within]; %#ok<AGROW>
end

T = array2table(rows, 'VariableNames', {'mu', 'ra_cnt', 'safe_cnt', 'ra_over_safe', ...
    'u1_lo', 'u1_hi', 'u2_lo', 'u2_hi', 'max_abs_u', 'within_ub'});
T.ub_req = repmat(ub_req, height(T), 1);
T.lambda = repmat(k1_lambda, height(T), 1);
writetable(T, fullfile(data_dir, 'mu_sweep_manipulator.csv'));
fprintf('\nwrote %s\n', fullfile(data_dir, 'mu_sweep_manipulator.csv'));
