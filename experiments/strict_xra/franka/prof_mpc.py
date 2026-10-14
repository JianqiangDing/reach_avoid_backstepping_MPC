import sys, time, os, numpy as np
os.environ["MPC_HESSIAN"] = sys.argv[1]
import config_franka as cf, mpc_franka as mp
X0 = cf.initial_set(); x = X0[:, 0]
for term in ("target", "cert"):
    t0 = time.time(); m = mp.FrankaMPC(cf.C, dt=cf.T, N=60, terminal=term, u_max=cf.U_MAX); t1 = time.time(); m._nlp(60); t2 = time.time(); ok, U, X = m.solve(x, 60); t3 = time.time(); ok2, U2, X2 = m.solve(x, 60, U); t4 = time.time()
    print(f"{sys.argv[1]} {term}: object {t1-t0:.1f}s, build N=60 {t2-t1:.1f}s, first solve ok {ok} {t3-t2:.1f}s, warm solve ok {ok2} {t4-t3:.2f}s", flush=True)
