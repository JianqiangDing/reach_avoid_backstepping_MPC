% Driver: check decoupling-matrix A(x) invertibility / feedback-linearizability
% over each example system's region. Writes data/decoupling_<system>.csv for the
% companion notebook scripts/check_decoupling_matrix.ipynb to visualize.
clc; clear;

repo_root = fileparts(fileparts(mfilename('fullpath')));
addpath(fullfile(repo_root, 'matlab'));
data_dir = fullfile(repo_root, 'data');
if ~exist(data_dir, 'dir'); mkdir(data_dir); end

n_samples = 30000;

% ---------------------------------------------------------------- double integrator
syms x1 x2 y1
check_decoupling_invertibility('double_integrator', ...
    [x2; 0], [0; 1], x1, [x1; x2], y1, 1 - y1^2, ...
    [-1.1; -1.1], [1.1; 1.1], n_samples, ...
    fullfile(data_dir, 'decoupling_double_integrator.csv'));
clear x1 x2 y1

% ---------------------------------------------------------------- Dubins car
syms x1 x2 th v y1 y2
h_raw = -(y1^4 + y2^4 - 16) * (y1^4 + y2^4 - 4);
target = (y2)^2 + ((y1 + 1.7) / 0.5)^2 - 0.4;
safe = 1e-3 * (-target + 300) * h_raw;
check_decoupling_invertibility('dubins_car', ...
    [v * cos(th); v * sin(th); 0; 0], [0, 0; 0, 0; 1, 0; 0, 1], [x1; x2], ...
    [x1; x2; th; v], [y1; y2], safe, ...
    [-2; -2; 2 * pi / 3; 0.1], [2; 2; 4 * pi / 3; 1.0], n_samples, ...
    fullfile(data_dir, 'decoupling_dubins_car.csv'));
clear x1 x2 th v y1 y2 h_raw target safe

% ---------------------------------------------------------------- 2-DoF manipulator
syms x1 x2 x3 x4 x5 y1 y2
m1 = 1.0; m2 = 1.0; l1 = 4.0; l2 = 4.0; lc1 = 2.0; lc2 = 2.0;
I1 = 0.02; I2 = 0.02; g = 9.81;
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
fx = [x3; x4; qddot_f; x3 + x4];
gx = [zeros(2, 2); Minv; zeros(1, 2)];
hx = [l1 * cos(x1) + l2 * cos(x5); l1 * sin(x1) + l2 * sin(x5)];
safe = -((4 * (y1 - 2) - 2 * y2^3)^2) + 0.8 * y2^3 + 10;
check_decoupling_invertibility('manipulator', ...
    fx, gx, hx, [x1; x2; x3; x4; x5], [y1; y2], safe, ...
    [-2.0; 0.8; -0.5; -0.5; -1.2], [0.5; pi - 0.8; 0.5; 0.5; pi - 0.3], n_samples, ...
    fullfile(data_dir, 'decoupling_manipulator.csv'), ...
    @(X) [X(:, 1:4), X(:, 1) + X(:, 2)]);  % x5 = x1 + x2 (physical state manifold)

fprintf('\nAll decoupling checks done.\n');
