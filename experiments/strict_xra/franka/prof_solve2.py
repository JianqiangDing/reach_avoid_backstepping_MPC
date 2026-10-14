import sys, time, numpy as np
import config_franka as cf, mpc_franka as mp
X0 = cf.initial_set(); Y0 = cf.C.y(X0); i = int(np.argmin(Y0[0])); j = int(np.argmin(((Y0[:2] - cf.TARGET[:, None])**2).sum(0)))
for nsub, hess in ((2, "exact"), (2, "limited-memory"), (5, "limited-memory")):
    for idx, lab in ((i, "far"), (j, "near")):
        for term in ("target", "cert"):
            m = mp.FrankaMPC(cf.C, dt=cf.T, N=cf.N, dfrac=cf.DFRAC, u_max=cf.U_MAX, nsub=nsub, terminal=term, hessian=hess); m._retry = True
            opti = m._nlp(cf.N)[0]; t0 = time.time(); ok, U, X = m.solve(X0[:, idx], cf.N); dt = time.time() - t0; st = opti.stats()
            fe = m.feasible(X, U) if ok else False
            print(f"nsub {nsub} {hess:14s} {lab:4s} {term:6s}: ok {ok} (rollout feasible {fe}), {dt:6.2f}s, iterations {st.get('iter_count')}, {dt/max(1, st.get('iter_count', 1)):.3f}s/iter, cost {m.cost(X, U) if ok else np.nan:.4g}", flush=True)
