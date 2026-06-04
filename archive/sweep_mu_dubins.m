% Sweep the backstepping certificate parameter mu_val for the Dubins example and
% report the resulting control-law bounds, to test whether mu affects whether the
% synthesized controller stays within the input limits |u| <= 5.
%
% Rationale: the bound violation is k1-independent (regularizing k1 had no effect);
% the controller magnitude u = A(x)^{-1} b(x) is dominated by drift/certificate terms
% that depend on mu (and lambda). This sweep checks empirically whether mu is a lever.
%
% The symbolic controller is synthesized once (mu stays symbolic); only the SOP
% re-solve + sampling-based bound estimate run per mu value. The unconstrained
% controller export is sent to a throwaway filename so committed controllers are
% left untouched. Run:  matlab -batch "addpath('matlab'); sweep_mu_dubins"
clc; clear;

repo_root = fileparts(fileparts(mfilename('fullpath')));
addpath(fullfile(repo_root, 'matlab'));

% ── Dubins system (matches example_dubins_car.m, FL-restricted region v>=0.1) ──
syms x1 x2 th v y1 y2;
x_vars_sym = [x1; x2; th; v];
fx_sym = [v * cos(th); v * sin(th); 0; 0];
gx_sym = [0, 0; 0, 0; 1, 0; 0, 1];
hx_sym = [x1; x2];
y_vars_sym = [y1; y2];
h_raw = -(y1^4 + y2^4 - 16) * (y1^4 + y2^4 - 4);
target_set_sym = (y2 - 0)^2 + ((y1 + 1.7) / 0.5)^2 - 0.4;
alpha = 1e-3 * (-target_set_sym + 300);
safe_set_sym = alpha * h_raw;

% synthesize the symbolic controller ONCE (mu is symbolic, substituted per sweep)
[u, k1, J_k1, mu, lambda, certificate, cert_term_dict, A_matrix, b_vector, ks, p, r_deg] = ...
    reach_avoid_controller(fx_sym, gx_sym, hx_sym, x_vars_sym, y_vars_sym, safe_set_sym);

% ── fixed settings (match example_dubins_car.m) ──
lb = [-5; -5];
ub = [5; 5];
ds = 4;
dv = 4;
samples_num = 1000;
bound_min = [-2; -2; 2 * pi / 3; 0.1];
bound_max = [2; 2; 4 * pi / 3; 1.0];

% reach-avoid region-size proxy: fraction of the box / of the safe set that the
% solved certificate certifies. Smaller mu shrinks the control magnitude BUT also
% shrinks this certified reach-avoid set, so we report both to see the tradeoff.
safe_set_x = subs(safe_set_sym, y_vars_sym, hx_sym);
safe_fun = matlabFunction(safe_set_x, 'Vars', x_vars_sym);
n_state = numel(x_vars_sym);

% ── mu values to test (current example uses mu_val = 0.1) ──
mu_list = [0.006, 0.008, 0.01, 0.02, 0.05, 0.1];

fprintf('\n==================== mu sweep (Dubins, limit |u|<=5) ====================\n');
for mu_val = mu_list
    rng(42); % identical state sampling across mu values, isolating mu's effect
    try
        [u_opt, cert_opt, valid_count, k1_opt] = solvesop_bounded_control( ...
            u, k1, J_k1, mu, lambda, certificate, cert_term_dict, p, r_deg, ...
            x_vars_sym, y_vars_sym, hx_sym, safe_set_sym, target_set_sym, ...
            mu_val, lb, ub, ds, dv, samples_num, bound_min, bound_max, ...
            'sweep_tmp_dubins_unconstrained.py');

        rng(7); % fixed seed so the bound estimate is comparable across mu values
        [lb1, ub1] = compute_poly_bounds_sampling(x_vars_sym, u_opt(1), cert_opt, 10000, bound_min, bound_max);
        [lb2, ub2] = compute_poly_bounds_sampling(x_vars_sym, u_opt(2), cert_opt, 10000, bound_min, bound_max);

        % reach-avoid set size at this mu: box samples with cert_opt >= 0 (and safe >= 0)
        rng(11); % identical region-sampling grid across mu values
        Nra = 20000;
        Xra = bound_min' + (bound_max - bound_min)' .* rand(Nra, n_state);
        cra = num2cell(Xra, 1);
        cert_opt_fun = matlabFunction(cert_opt, 'Vars', x_vars_sym);
        cv = cert_opt_fun(cra{:}); if isscalar(cv); cv = cv * ones(Nra, 1); end
        sv = safe_fun(cra{:}); if isscalar(sv); sv = sv * ones(Nra, 1); end
        safe_cnt = sum(sv(:) >= 0);
        ra_cnt = sum(cv(:) >= 0 & sv(:) >= 0);
        ra_over_safe = ra_cnt / max(safe_cnt, 1);

        within = (lb1 >= -5) && (ub1 <= 5) && (lb2 >= -5) && (ub2 <= 5);
        fprintf('__SWEEP__ mu=%.4g valid=%d/%d ra_cnt=%d ra/safe=%.3f u1=[%.4g,%.4g] u2=[%.4g,%.4g] within_pm5=%d\n', ...
            mu_val, valid_count, samples_num, ra_cnt, ra_over_safe, lb1, ub1, lb2, ub2, within);
    catch ME
        fprintf('__SWEEP__ mu=%.4g FAILED: %s\n', mu_val, ME.message);
    end
end
fprintf('==========================================================================\n');

% clean up the throwaway export
tmp_export = fullfile(repo_root, 'controllers', 'sweep_tmp_dubins_unconstrained.py');
if exist(tmp_export, 'file'); delete(tmp_export); end
