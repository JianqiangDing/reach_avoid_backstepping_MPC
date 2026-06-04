% Pin down WHY the manipulator reach-avoid set empties as xi0 (-> lambda) rises.
% For rd=2 the certificate is  V = safe(h) - sum_i (1/(2 mu_i)) (Lf h_i - k1_i)^2,
% and lambda enters ONLY through the solved (vanilla) k1. We sweep xi0, re-solve the
% vanilla k1, and over the FL box report:
%   - solved lambda
%   - max |k1| (does k1 grow with lambda?)
%   - median / max of the penalty  sum_i (1/(2 mu_i))(Lf h_i - k1_i)^2  on {safe>=0}
%   - min V  and  raset = fraction of box with V>=0   (where does it cross to 0?)
% Cheap: only the small vanilla SOS solve + certificate eval (no heavy SOP / no A^{-1}b).

clc; clear; close all;
script_dir = fileparts(mfilename('fullpath'));
addpath(script_dir);
data_dir = fullfile(fileparts(script_dir), 'data');
if ~exist(data_dir, 'dir'); mkdir(data_dir); end

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

[~, k1, ~, mu, lambda, certificate, ~, ~, ~, ~, ~, ~] = ...
    reach_avoid_controller(fx, gx, hx, x_vars, y_vars, safe_set);

mu_val = 15; ds = 4; dv = 2;
bound_min = [-2; 0.8; -0.5; -0.5; -1.2]; bound_max = [0.5; pi - 0.8; 0.5; 0.5; pi - 0.3];
N = 20000;
rng(42);
samples = rand(N, 5) .* (bound_max - bound_min)' + bound_min';
cols = num2cell(samples, 1);

safe_x = subs(safe_set, y_vars, hx);
safe_fun = matlabFunction(safe_x, 'vars', x_vars);
sv = safe_fun(cols{:}) .* ones(N, 1);
insafe = sv >= 0;
fprintf('box samples in {safe>=0}: %d / %d\n', nnz(insafe), N);

xi0_grid = [1e-8, 1e-4, 1e-3, 1e-2, 0.03, 0.05, 0.1, 0.5];
fprintf('\n%-10s %-12s %-12s %-12s %-12s %-12s %-10s\n', ...
    'xi0', 'lambda', 'max|k1|', 'med_penalty', 'max_penalty', 'min_V', 'raset');
rows = [];
for xi0 = xi0_grid
    [k1_y, lam, ~, ~] = solve_vanilla_k1_controller_xi(y_vars, safe_set, target_set, dv, ds, xi0);

    % k1 in x-space, magnitude over the box
    k1_x = subs(k1_y, y_vars, hx);
    k1f1 = matlabFunction(k1_x(1), 'vars', x_vars);
    k1f2 = matlabFunction(k1_x(2), 'vars', x_vars);
    k1v = max(abs(k1f1(cols{:}) .* ones(N, 1)), abs(k1f2(cols{:}) .* ones(N, 1)));

    % V(x) with the vanilla k1
    V = subs(certificate, k1, k1_y);
    V = subs(V, y_vars, hx);
    for i = 1:numel(mu)
        if ~isempty(mu{i}); V = subs(V, mu{i}, mu_val); end
    end
    V = subs(V, lambda, lam);
    Vf = matlabFunction(V, 'vars', x_vars);
    Vv = Vf(cols{:}) .* ones(N, 1);

    penalty = sv - Vv;                       % = sum_i (1/2mu)(Lf h_i - k1_i)^2
    raset = mean(Vv >= 0);                   % fraction of box (V>=0 implies safe>=0)
    rows = [rows; xi0, lam, max(k1v(insafe)), median(penalty(insafe)), ...
        max(penalty(insafe)), min(Vv(insafe)), raset]; %#ok<AGROW>
    fprintf('%-10.3g %-12.4g %-12.4g %-12.4g %-12.4g %-12.4g %-10.4f\n', ...
        xi0, lam, max(k1v(insafe)), median(penalty(insafe)), max(penalty(insafe)), ...
        min(Vv(insafe)), raset);
end

T = array2table(rows, 'VariableNames', ...
    {'xi0', 'lambda', 'max_k1', 'med_penalty', 'max_penalty', 'min_V', 'raset'});
writetable(T, fullfile(data_dir, 'manip_collapse.csv'));
fprintf('\nwrote %s\n', fullfile(data_dir, 'manip_collapse.csv'));
