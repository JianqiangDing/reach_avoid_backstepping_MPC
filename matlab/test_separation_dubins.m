% test_separation_dubins.m
% ------------------------------------------------------------------------
% Demonstration plan: at mu = 0.1, raise the SOS floor xi0 to ~50-100x mu so
% that k1 gains leverage on u, then synthesize TWO controllers at the SAME
% lambda via the vanilla-certificate-region SOP (solvesop_bounded_control_slack_xi):
%     - ux_con   : input-constrained k1 (bound |u| <= ub baked in, slack SOP)
%     - ux_uncon : vanilla (unconstrained) k1
% Build both closed-loop controllers u = A(x)^{-1} b(x) and search over initial
% conditions x0 in the reach-avoid set X_RA for a SEPARATION:
%   (1) two-trajectory: con reaches target with |u|<=ub, uncon's own trajectory
%       violates |u|>ub somewhere (both reach);
%   (2) single-trajectory: along the CON trajectory (which reaches with |u|<=ub),
%       the unconstrained controller WOULD command |u_uncon|>ub somewhere.
%
% Sampling region = vanilla certificate region (the hard certificate LMI cannot
% hold over the whole X_S\X_T, so direct-set sampling is infeasible -- confirmed).
%
% NEW FILE. Does NOT touch controllers/ (no export_to_python). Writes only to
% data/separation_dubins_*.csv and data/sep_best_*.csv.
%
% Run:  matlab -batch "addpath('matlab'); test_separation_dubins"
% ------------------------------------------------------------------------
clc; clear; close all;
script_dir = fileparts(mfilename('fullpath'));
addpath(script_dir);
repo_root = fileparts(script_dir);
data_dir = fullfile(repo_root, 'data');
if ~exist(data_dir, 'dir'); mkdir(data_dir); end

% ---------------- Dubins system (matches example_dubins_car.m) ----------
syms x1 x2 th v y1 y2;
x_vars = [x1; x2; th; v];
fx = [v*cos(th); v*sin(th); 0; 0];
gx = [0, 0; 0, 0; 1, 0; 0, 1];
hx = [x1; x2];
y_vars = [y1; y2];
h_raw = -(y1^4 + y2^4 - 16) * (y1^4 + y2^4 - 4);
target_set_sym = (y2 - 0)^2 + ((y1 + 1.7)/0.5)^2 - 0.4;   % X_T = {target_set < 0}
alpha = 1e-3 * (-target_set_sym + 300);
safe_set_sym = alpha * h_raw;                              % X_S = {safe_set > 0}

% ---------------- fixed settings ----------------------------------------
mu_val = 0.1;
lb = [-5; -5];
ub = [5; 5];
ubnd = 5;                            % demonstration bound (per channel)
ds = 4;
dv = 4;
samples_num = 1000;
bound_min = [-2; -2; 2*pi/3; 0.1];   % v >= 0.1 keeps det A = -v nonsingular
bound_max = [ 2;  2; 4*pi/3; 1.0];
xi0_list = 10;                       % 100x mu

% numeric helpers
f_fn      = matlabFunction(fx,            'Vars', {x_vars});   % f(X) -> 4x1
g_const   = double(gx);
safe_fn   = matlabFunction(subs(safe_set_sym,  y_vars, hx), 'Vars', {x_vars});
target_fn = matlabFunction(subs(target_set_sym, y_vars, hx), 'Vars', {x_vars});

% ---------------- synthesize symbolic reach-avoid controller (once) ------
fprintf('Synthesizing symbolic reach-avoid controller...\n');
[u, k1, J_k1, mu, lambda, certificate, cert_term_dict, A_matrix, b_vector, ks, p, r_deg] = ...
    reach_avoid_controller(fx, gx, hx, x_vars, y_vars, safe_set_sym);

% ======================================================================
for xi0 = xi0_list
    fprintf('\n=================== xi0 = %g  (=%gx mu) ===================\n', xi0, xi0/mu_val);
    rng(42);
    [ux_con, ux_uncon, cert_con, cert_uncon, valid_count, k1_con, slack_opt, lam] = ...
        solvesop_bounded_control_slack_xi(u, k1, J_k1, mu, lambda, certificate, ...
        cert_term_dict, p, r_deg, x_vars, y_vars, hx, safe_set_sym, target_set_sym, ...
        mu_val, lb, ub, ds, dv, samples_num, bound_min, bound_max, xi0);
    fprintf('  lambda=%.4g  slack=[%.4g %.4g]  achievable|u|<=[%.4g %.4g]  valid_samples=%d\n', ...
        lam, slack_opt(1), slack_opt(2), ub(1)+slack_opt(1), ub(2)+slack_opt(2), valid_count);

    ucon_fn   = matlabFunction(ux_con,   'Vars', {x_vars});
    uuncon_fn = matlabFunction(ux_uncon, 'Vars', {x_vars});
    certcon_fn   = matlabFunction(cert_con,   'Vars', {x_vars});
    certuncon_fn = matlabFunction(cert_uncon, 'Vars', {x_vars});

    % -------- sample initial conditions inside X_RA (both certs >= 0) -----
    rng(7);
    Nsamp = 12000;
    Xs = bound_min' + (bound_max - bound_min)' .* rand(Nsamp, 4);
    keep = false(Nsamp, 1);
    for idx = 1:Nsamp
        X = Xs(idx, :)';
        if safe_fn(X) >= 0 && target_fn(X) >= 0 && certcon_fn(X) >= 0 && certuncon_fn(X) >= 0
            keep(idx) = true;
        end
    end
    X0 = Xs(keep, :)';
    nX0 = size(X0, 2);
    if nX0 > 350; X0 = X0(:, 1:350); nX0 = 350; end
    fprintf('  %d candidate initial conditions inside X_RA\n', nX0);

    % -------- simulate + measure both separation criteria ----------------
    Tmax = 30;
    res = nan(nX0, 12);
    % cols: rc sc p1c p2c | ru su p1u p2u | p1u_onCon p2u_onCon | x1 x2
    best1 = struct('score', -inf);   % two-trajectory separation
    best2 = struct('score', -inf);   % single-(con)-trajectory separation
    for j = 1:nX0
        x0 = X0(:, j);
        [rc, sc, p1c, p2c, Tc, Xc, Uc] = sim_cl(f_fn, g_const, ucon_fn,   safe_fn, target_fn, x0, Tmax);
        [ru, su, p1u, p2u, Tu, Xu, Uu] = sim_cl(f_fn, g_const, uuncon_fn, safe_fn, target_fn, x0, Tmax);

        % uncon controller evaluated ALONG the con trajectory (counterfactual)
        p1uc = NaN; p2uc = NaN; Uuc = [];
        if rc && sc && ~isempty(Xc)
            Uuc = zeros(size(Xc, 1), 2);
            for kk = 1:size(Xc, 1)
                Uuc(kk, :) = uuncon_fn(Xc(kk, :)')';
            end
            p1uc = max(abs(Uuc(:, 1))); p2uc = max(abs(Uuc(:, 2)));
        end
        res(j, :) = [rc, sc, p1c, p2c, ru, su, p1u, p2u, p1uc, p2uc, x0(1), x0(2)];

        pc = max(p1c, p2c);
        % criterion (1): both reach safely, con<=ubnd, uncon-own-traj>ubnd
        if rc && sc && ru && su && isfinite(pc)
            pu = max(p1u, p2u);
            if pc <= ubnd && pu > ubnd
                sc1 = pu - pc;
                if sc1 > best1.score
                    best1 = struct('score', sc1, 'x0', x0, 'pc1', p1c, 'pc2', p2c, ...
                        'pu1', p1u, 'pu2', p2u, 'Tc', Tc, 'Xc', Xc, 'Uc', Uc, 'Tu', Tu, 'Xu', Xu, 'Uu', Uu);
                end
            end
        end
        % criterion (2): con reaches safely with con<=ubnd, uncon-on-con-traj>ubnd
        if rc && sc && isfinite(pc) && pc <= ubnd
            puc = max(p1uc, p2uc);
            if isfinite(puc) && puc > ubnd
                sc2 = puc - pc;
                if sc2 > best2.score
                    best2 = struct('score', sc2, 'x0', x0, 'pc1', p1c, 'pc2', p2c, ...
                        'puc1', p1uc, 'puc2', p2uc, 'Tc', Tc, 'Xc', Xc, 'Uc', Uc, 'Uuc', Uuc);
                end
            end
        end
    end

    % -------- summary ----------------------------------------------------
    reach_con = res(:,1)==1 & res(:,2)==1;
    reach_unc = res(:,5)==1 & res(:,6)==1;
    fprintf('  con reached&safe: %d/%d   uncon reached&safe: %d/%d\n', ...
        nnz(reach_con), nX0, nnz(reach_unc), nX0);
    n1 = nnz(reach_con & reach_unc & max(res(:,3),res(:,4))<=ubnd & max(res(:,7),res(:,8))>ubnd);
    n2 = nnz(reach_con & max(res(:,3),res(:,4))<=ubnd & max(res(:,9),res(:,10))>ubnd);
    fprintf('  criterion (1) two-trajectory separations : %d\n', n1);
    fprintf('  criterion (2) con-trajectory separations : %d\n', n2);

    T = array2table(res, 'VariableNames', {'rc','sc','p1c','p2c','ru','su','p1u','p2u', ...
        'p1u_onCon','p2u_onCon','x1_0','x2_0'});
    T.xi0 = repmat(xi0, height(T), 1);
    writetable(T, fullfile(data_dir, sprintf('separation_dubins_xi%g.csv', xi0)));

    if isfinite(best1.score)
        fprintf('  >>> BEST (1) two-traj @ x0=[% .3f % .3f % .3f % .3f]\n', best1.x0);
        fprintf('      peak|u_con|=[%.3g %.3g]  peak|u_uncon(own traj)|=[%.3g %.3g]\n', ...
            best1.pc1, best1.pc2, best1.pu1, best1.pu2);
        writematrix([best1.Tc(:), best1.Xc, best1.Uc], fullfile(data_dir, sprintf('sep1_con_xi%g.csv', xi0)));
        writematrix([best1.Tu(:), best1.Xu, best1.Uu], fullfile(data_dir, sprintf('sep1_unc_xi%g.csv', xi0)));
    else
        fprintf('  >>> NO criterion-(1) separation found.\n');
    end
    if isfinite(best2.score)
        fprintf('  >>> BEST (2) con-traj @ x0=[% .3f % .3f % .3f % .3f]\n', best2.x0);
        fprintf('      peak|u_con|=[%.3g %.3g]  peak|u_uncon ALONG con-traj|=[%.3g %.3g]\n', ...
            best2.pc1, best2.pc2, best2.puc1, best2.puc2);
        writematrix([best2.Tc(:), best2.Xc, best2.Uc, best2.Uuc], fullfile(data_dir, sprintf('sep2_xi%g.csv', xi0)));
    else
        fprintf('  >>> NO criterion-(2) separation found.\n');
    end
end
fprintf('\nDONE.\n');

% ======================================================================
function [reached, safe_ok, peak1, peak2, tout, Xout, Uout] = ...
        sim_cl(f_fn, g, u_fn, safe_fn, target_fn, x0, Tmax)
    reached = 0; safe_ok = 1; peak1 = NaN; peak2 = NaN;
    tout = []; Xout = []; Uout = [];
    odef = @(t, X) f_fn(X) + g * u_fn(X);
    opts = odeset('RelTol', 1e-6, 'AbsTol', 1e-8, 'MaxStep', 0.05, ...
        'Events', @(t, X) ev_fn(t, X, safe_fn, target_fn));
    try
        sol = ode45(odef, [0 Tmax], x0, opts);
    catch
        safe_ok = 0; return;
    end
    tout = sol.x(:); Xout = sol.y';
    Nt = numel(tout);
    U = zeros(Nt, 2);
    for k = 1:Nt
        U(k, :) = u_fn(Xout(k, :)')';
    end
    Uout = U;
    peak1 = max(abs(U(:, 1))); peak2 = max(abs(U(:, 2)));
    if ~isempty(sol.ie)
        if any(sol.ie == 1); reached = 1; end
        if any(sol.ie == 2); safe_ok = 0; end
        if any(sol.ie == 3); safe_ok = 0; end
    end
    sv = arrayfun(@(k) safe_fn(Xout(k, :)'), 1:Nt);
    if any(sv < -1e-6); safe_ok = 0; end
end

function [val, ister, dir] = ev_fn(~, X, safe_fn, target_fn)
    val   = [target_fn(X); X(4) - 0.05; safe_fn(X)];
    ister = [1; 1; 1];
    dir   = [-1; -1; -1];
end
