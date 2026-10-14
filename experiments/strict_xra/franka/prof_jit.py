import time, numpy as np, casadi as ca
import config_franka as cf, mpc_franka as mp
from cert_franka import FDYN
psi = cf.C.F["psi"]; x = ca.SX.sym("x", 8); H = ca.hessian(psi(x), x)[0]
print(f"nodes: psi {psi.n_nodes()}, Psi {cf.C.F['Psi'].n_nodes()}, FDYN {FDYN.n_nodes()}, hessian of psi {ca.Function('h', [x], [H]).n_nodes()}", flush=True)
X0 = cf.initial_set(); Y0 = cf.C.y(X0); i = int(np.argmin(Y0[0]))
for jit in (False, True):
    m = mp.FrankaMPC(cf.C, dt=cf.T, N=cf.N, dfrac=cf.DFRAC, u_max=cf.U_MAX, nsub=2, terminal="cert"); m._retry = True
    if jit: m.opts.update({"jit": True, "compiler": "shell", "jit_options": {"flags": ["-O1"], "verbose": False}})
    t0 = time.time(); opti = m._nlp(cf.N)[0]; t1 = time.time(); ok, U, X = m.solve(X0[:, i], cf.N); t2 = time.time(); st = opti.stats()
    ok2, U2, X2 = m.solve(X0[:, i], cf.N, U); t3 = time.time(); st2 = opti.stats()
    print(f"jit {jit}: build {t1-t0:.1f}s, first solve (incl. compile) {t2-t1:.1f}s ok {ok} iters {st.get('iter_count')}, warm solve {t3-t2:.2f}s iters {st2.get('iter_count')} -> {(t3-t2)/max(1, st2.get('iter_count', 1)):.3f}s/iter", flush=True)
