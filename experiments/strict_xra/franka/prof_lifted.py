import time, numpy as np
import config_franka as cf, mpc_lifted as ml
X0 = cf.initial_set(); Y0 = cf.C.y(X0); i = int(np.argmin(Y0[0])); j = int(np.argmin(((Y0[:2] - cf.TARGET[:, None])**2).sum(0)))
for term in ("target", "cert"):
    m = ml.LiftedMPC(cf.C, dt=cf.T, N=cf.N, dfrac=cf.DFRAC, u_max=cf.U_MAX, nsub=2, terminal=term); m._retry = True; t0 = time.time(); m._solver(); print(f"{term}: build {time.time()-t0:.1f}s", flush=True)
    for idx, lab in ((i, "far"), (j, "near")):
        for Na in ((cf.N, 30) if term == "cert" else (cf.N,)):
            t0 = time.time(); ok, U, X = m.solve(X0[:, idx], Na); dt = time.time() - t0; st = m._S.stats()
            print(f"  {lab} Na={Na}: ok {ok}, {dt:.2f}s, iterations {st['iter_count']}, {dt/max(1, st['iter_count']):.3f}s/iter, status {st['return_status']}, rollout feasible {m.feasible(X, U) if ok else '-'}", flush=True)
