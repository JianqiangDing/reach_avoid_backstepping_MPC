import time, numpy as np
import config_franka as cf
t0 = time.time(); X0 = cf.initial_set(); cu = cf.C0.scaled(cf.URA_SCALE); print(f"{X0.shape[1]} initial states ({time.time()-t0:.1f}s); in the URA's X_RA {(cu.Psi(X0) >= 0).mean()*100:.0f}%, terminal controller's X_RA {(cf.C.Psi(X0) >= 0).mean()*100:.0f}%", flush=True)
Y0 = cf.C.y(X0); i = int(np.argmin(Y0[0])); print("farthest start y =", np.round(Y0[:, i], 3), flush=True)
t0 = time.time(); u = cf.ura_run(X0[:, i]); print(f"URA: {u['status']} {u['t']:.2f} s, peaks {np.round(u['upeak'], 2)} vs bounds {cf.U_MAX} ({time.time()-t0:.1f}s)", flush=True)
import mpc_franka as mp
for term in ("target", "cert"):
    t0 = time.time(); m = mp.FrankaMPC(cf.C, dt=cf.T, N=cf.N, dfrac=cf.DFRAC, u_max=cf.U_MAX, nsub=5, terminal=term); t1 = time.time(); m._nlp(cf.N); t2 = time.time()
    print(f"{term}: object {t1-t0:.1f}s, build NLP {t2-t1:.1f}s", flush=True)
    L = m.run(X0[:, i], T_max=cf.T_MAX_MPC); print(f"{term}: {L['status']}, steps {len(L['J'])}, t {L['t'][-1]:.2f} s, solve ms mean {np.mean(L['solve_ms']) if len(L['solve_ms']) else np.nan:.0f} (first {L['solve_ms'][0] if len(L['solve_ms']) else np.nan:.0f}), omega {L['omega']:.3g}, eps {L['eps']:.3g} ({time.time()-t2:.0f}s)", flush=True)
