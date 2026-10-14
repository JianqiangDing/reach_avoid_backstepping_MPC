import time, numpy as np
import config_franka as cf, mpc_franka as mp
X0 = cf.initial_set(); Y0 = cf.C.y(X0); i = int(np.argmin(Y0[0])); j = int(np.argmin(((Y0[:2] - cf.TARGET[:, None])**2).sum(0)))
for idx, lab in ((j, "nearest start"), (i, "farthest start")):
    for term in ("target", "cert"):
        m = mp.FrankaMPC(cf.C, dt=cf.T, N=cf.N, dfrac=cf.DFRAC, u_max=cf.U_MAX, nsub=5, terminal=term); m._retry = True      # one solve, no fallbacks
        opti = m._nlp(cf.N)[0]; t0 = time.time(); ok, U, X = m.solve(X0[:, idx], cf.N); dt = time.time() - t0
        st = opti.stats(); print(f"{lab} {term}: ok {ok}, {dt:.1f}s, iterations {st.get('iter_count')}, status {st.get('return_status')}, per iteration {dt/max(1, st.get('iter_count', 1)):.3f}s", flush=True)
