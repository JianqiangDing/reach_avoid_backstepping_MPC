% test_xi_sweep_dubins.m
% ------------------------------------------------------------------------
% Sweep the SOS floor xi0 over the 50x-100x mu range (mu=0.1 -> xi0 in [5,10])
% to find which lambda gives the best constrained-vs-unconstrained SEPARATION
% for the Dubins demonstration. For each xi0 we synthesize both controllers via
% the vanilla-certificate-region slack SOP and report:
%   - lambda, slack (=> constrained achievable |u| bound = ub + slack)
%   - distinguishability dist% = max|u_con - u_uncon|/max|u_uncon| over X_RA
%   - dense |u| range of con and of uncon over X_RA
%   - a LIGHT trajectory separation count from ~50 initial conditions in X_RA:
%       (1) two-traj : con reaches with |u|<=5, uncon's own traj has |u|>5
%       (2) con-traj : con reaches with |u|<=5, uncon WOULD command |u|>5 on it
%     + the best (most pronounced) margin for each.
%
% NEW FILE. No export_to_python (controllers/ untouched). Writes data/xi_sweep_dubins.csv.
% Run:  matlab -batch "addpath('matlab'); test_xi_sweep_dubins"
% ------------------------------------------------------------------------
clc; clear; close all;
script_dir = fileparts(mfilename('fullpath'));
addpath(script_dir);
data_dir = fullfile(fileparts(script_dir), 'data');
if ~exist(data_dir, 'dir'); mkdir(data_dir); end

syms x1 x2 th v y1 y2;
x_vars = [x1; x2; th; v];
fx = [v*cos(th); v*sin(th); 0; 0];
gx = [0, 0; 0, 0; 1, 0; 0, 1];
hx = [x1; x2];
y_vars = [y1; y2];
h_raw = -(y1^4 + y2^4 - 16) * (y1^4 + y2^4 - 4);
target_set_sym = (y2 - 0)^2 + ((y1 + 1.7)/0.5)^2 - 0.4;
alpha = 1e-3 * (-target_set_sym + 300);
safe_set_sym = alpha * h_raw;

mu_val = 0.1;
lb = [-5; -5]; ub = [5; 5]; ubnd = 5;
ds = 4; dv = 4; samples_num = 1000;
bound_min = [-2; -2; 2*pi/3; 0.1];
bound_max = [ 2;  2; 4*pi/3; 1.0];
xi0_grid = [5, 6, 7, 8, 9, 10];      % 50x .. 100x mu
n_x0 = 50;                           % light search per xi0
Tmax = 25;

f_fn      = matlabFunction(fx, 'Vars', {x_vars});
g_const   = double(gx);
safe_fn   = matlabFunction(subs(safe_set_sym,  y_vars, hx), 'Vars', {x_vars});
target_fn = matlabFunction(subs(target_set_sym, y_vars, hx), 'Vars', {x_vars});

fprintf('Synthesizing symbolic reach-avoid controller...\n');
[u, k1, J_k1, mu, lambda, certificate, cert_term_dict, A_matrix, b_vector, ks, p, r_deg] = ...
    reach_avoid_controller(fx, gx, hx, x_vars, y_vars, safe_set_sym);

rows = [];
fprintf('\n%-5s %-7s %-16s %-7s %-22s %-22s | %-9s %-9s | sep1 sep2  best1   best2\n', ...
    'xi0','lambda','con_achiev|u|','dist%','con_dense u1/u2','unc_dense u1/u2','reachCon','reachUnc');
for xi0 = xi0_grid
    rng(42);
    try
        [ux_con, ux_uncon, cert_con, cert_uncon, valid_count, k1_con, slack_opt, lam] = ...
            solvesop_bounded_control_slack_xi(u, k1, J_k1, mu, lambda, certificate, ...
            cert_term_dict, p, r_deg, x_vars, y_vars, hx, safe_set_sym, target_set_sym, ...
            mu_val, lb, ub, ds, dv, samples_num, bound_min, bound_max, xi0);
    catch ME
        fprintf('xi0=%-4g SOLVE FAILED: %s\n', xi0, ME.message);
        continue;
    end

    ucon_fn   = matlabFunction(ux_con,   'Vars', {x_vars});
    uuncon_fn = matlabFunction(ux_uncon, 'Vars', {x_vars});
    certcon_fn   = matlabFunction(cert_con,   'Vars', {x_vars});
    certuncon_fn = matlabFunction(cert_uncon, 'Vars', {x_vars});

    % dense ranges + distinguishability over X_RA
    [c1lo,c1hi] = compute_poly_bounds_sampling(x_vars, ux_con(1),   cert_con,   10000, bound_min, bound_max);
    [c2lo,c2hi] = compute_poly_bounds_sampling(x_vars, ux_con(2),   cert_con,   10000, bound_min, bound_max);
    [u1lo,u1hi] = compute_poly_bounds_sampling(x_vars, ux_uncon(1), cert_uncon, 10000, bound_min, bound_max);
    [u2lo,u2hi] = compute_poly_bounds_sampling(x_vars, ux_uncon(2), cert_uncon, 10000, bound_min, bound_max);
    con_dense = max(abs([c1lo c1hi]));  con_dense2 = max(abs([c2lo c2hi]));
    unc_dense = max(abs([u1lo u1hi]));  unc_dense2 = max(abs([u2lo u2hi]));

    rng(123); Nd = 8000;
    Xd = bound_min' + (bound_max-bound_min)'.*rand(Nd,4);
    inRA = false(Nd,1); ducon = zeros(Nd,2); duunc = zeros(Nd,2);
    for i=1:Nd
        X = Xd(i,:)';
        if certuncon_fn(X) >= 0; inRA(i)=true; ducon(i,:)=ucon_fn(X)'; duunc(i,:)=uuncon_fn(X)'; end
    end
    reg = inRA; if nnz(reg)<50; reg=true(Nd,1); end
    dist = 100*max(max(abs(ducon(reg,:)-duunc(reg,:))))/max(max(max(abs(duunc(reg,:)))),eps);

    % light trajectory separation search
    rng(7); Ns=12000; Xs = bound_min'+(bound_max-bound_min)'.*rand(Ns,4); keep=false(Ns,1);
    for i=1:Ns
        X=Xs(i,:)';
        if safe_fn(X)>=0 && target_fn(X)>=0 && certcon_fn(X)>=0 && certuncon_fn(X)>=0; keep(i)=true; end
    end
    X0 = Xs(keep,:)'; if size(X0,2)>n_x0; X0=X0(:,1:n_x0); end
    nX0=size(X0,2); rcN=0; ruN=0; sep1=0; sep2=0; best1=0; best2=0;
    for j=1:nX0
        x0=X0(:,j);
        [rc,sc,p1c,p2c,~,Xc,~] = sim_cl(f_fn,g_const,ucon_fn,safe_fn,target_fn,x0,Tmax);
        [ru,su,p1u,p2u,~,~,~]  = sim_cl(f_fn,g_const,uuncon_fn,safe_fn,target_fn,x0,Tmax);
        if rc&&sc; rcN=rcN+1; end
        if ru&&su; ruN=ruN+1; end
        pc=max(p1c,p2c);
        if rc&&sc&&ru&&su&&isfinite(pc)&&pc<=ubnd && max(p1u,p2u)>ubnd
            sep1=sep1+1; best1=max(best1, max(p1u,p2u)-pc);
        end
        if rc&&sc&&isfinite(pc)&&pc<=ubnd && ~isempty(Xc)
            puc=0; for kk=1:size(Xc,1); uu=uuncon_fn(Xc(kk,:)'); puc=max(puc,max(abs(uu))); end
            if puc>ubnd; sep2=sep2+1; best2=max(best2, puc-pc); end
        end
    end

    fprintf('%-5g %-7.3g [%.3g %.3g]%4s %-7.3g [%.3g/%.3g]%6s [%.3g/%.3g]%6s | %3d/%-5d %3d/%-5d | %3d  %3d  %5.2f  %5.2f\n', ...
        xi0, lam, ub(1)+slack_opt(1), ub(2)+slack_opt(2), '', dist, ...
        con_dense, con_dense2, '', unc_dense, unc_dense2, '', rcN, nX0, ruN, nX0, sep1, sep2, best1, best2);
    rows=[rows; xi0, lam, slack_opt(1), slack_opt(2), dist, con_dense, con_dense2, unc_dense, unc_dense2, ...
        rcN, ruN, nX0, sep1, sep2, best1, best2]; %#ok<AGROW>
end

T = array2table(rows, 'VariableNames', {'xi0','lambda','slack1','slack2','dist_pct', ...
    'con_dense1','con_dense2','unc_dense1','unc_dense2','reach_con','reach_unc','nX0', ...
    'sep1_count','sep2_count','sep1_best_margin','sep2_best_margin'});
writetable(T, fullfile(data_dir,'xi_sweep_dubins.csv'));
fprintf('\nwrote %s\nDONE.\n', fullfile(data_dir,'xi_sweep_dubins.csv'));

% ======================================================================
function [reached, safe_ok, peak1, peak2, tout, Xout, Uout] = ...
        sim_cl(f_fn, g, u_fn, safe_fn, target_fn, x0, Tmax)
    reached=0; safe_ok=1; peak1=NaN; peak2=NaN; tout=[]; Xout=[]; Uout=[];
    odef = @(t,X) f_fn(X) + g*u_fn(X);
    opts = odeset('RelTol',1e-6,'AbsTol',1e-8,'MaxStep',0.05,'Events',@(t,X) ev_fn(t,X,safe_fn,target_fn));
    try; sol = ode45(odef,[0 Tmax],x0,opts); catch; safe_ok=0; return; end
    tout=sol.x(:); Xout=sol.y'; Nt=numel(tout); U=zeros(Nt,2);
    for k=1:Nt; U(k,:)=u_fn(Xout(k,:)')'; end
    Uout=U; peak1=max(abs(U(:,1))); peak2=max(abs(U(:,2)));
    if ~isempty(sol.ie)
        if any(sol.ie==1); reached=1; end
        if any(sol.ie==2); safe_ok=0; end
        if any(sol.ie==3); safe_ok=0; end
    end
    sv=arrayfun(@(k) safe_fn(Xout(k,:)'),1:Nt); if any(sv<-1e-6); safe_ok=0; end
end
function [val,ister,dir] = ev_fn(~,X,safe_fn,target_fn)
    val=[target_fn(X); X(4)-0.05; safe_fn(X)]; ister=[1;1;1]; dir=[-1;-1;-1];
end
