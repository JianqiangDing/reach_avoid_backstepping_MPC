import time, sys, numpy as np, casadi as ca
import config_franka as cf, mpc_franka as mp
CL = "/Users/jianqiang/miniforge3/envs/rab_franka/bin/arm64-apple-darwin20.0.0-clang"; SDK = "/Library/Developer/CommandLineTools/SDKs/MacOSX.sdk"
X0 = cf.initial_set(); Y0 = cf.C.y(X0); i = int(np.argmin(Y0[0]))
for N, opt in ((60, "-O1"), (20, "-O1")):
    m = mp.FrankaMPC(cf.C, dt=cf.T, N=N, dfrac=cf.DFRAC, u_max=cf.U_MAX, nsub=2, terminal="cert"); m._retry = True
    m.opts.update({"jit": True, "compiler": "shell", "jit_options": {"compiler": CL, "linker": CL, "compiler_flags": [opt, "-fPIC", "-isysroot", SDK], "linker_flags": ["-isysroot", SDK, "-fPIC"], "verbose": False}})
    t0 = time.time(); opti = m._nlp(N)[0]; t1 = time.time(); ok, U, X = m.solve(X0[:, i], N); t2 = time.time(); it1 = opti.stats()["iter_count"] if ok else -1
    if ok:
        ok2, U2, X2 = m.solve(X0[:, i], N, U); t3 = time.time(); it2 = opti.stats()["iter_count"]
        print(f"N={N} {opt}: build {t1-t0:.1f}s, first solve incl. compile {t2-t1:.1f}s (iters {it1}), warm solve {t3-t2:.3f}s, iters {it2} -> {(t3-t2)/max(1,it2):.4f}s/iter", flush=True)
    else: print(f"N={N}: first solve failed ({t2-t1:.1f}s)", flush=True)
