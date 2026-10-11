import numpy as np, sys
import sample2_reconstruction as S
d = np.load("body4_mesh.npz"); V, F, N, W = d["V"], d["F"], d["N"], d["W"]
fr = np.load("s4_frame.npz"); n_up = fr["n_up"]
under = (N @ n_up) < -0.35
print("underside vertices: %d of %d (area share %.2f)" % (under.sum(), len(V), W[under].sum()/W.sum()))
rng = np.random.default_rng(21)
Pu, Nu, Wu = V[under], N[under], W[under]
sel = rng.choice(len(Pu), size=2200, replace=False, p=Wu/Wu.sum())
sols, cen = S.global_axis_search(Pu[sel], Nu[sel], n_dirs=500, keep=14)
print("\nrank  rms(mm)  direction                 angle to plate normal   distance axis-to-sherd centre (mm)")
out = []
for k, (f, p) in enumerate(sols):
    a, c = S.axis_from(p, cen)
    if a @ n_up < 0: a = -a
    ang = np.degrees(np.arccos(np.clip(a @ n_up, -1, 1)))
    dist = np.linalg.norm(np.cross(cen - c, a))
    out.append((f, a, c, ang, dist))
    print("%3d  %6.3f   %s     %5.1f                 %7.1f" % (k+1, f, np.round(a, 3), ang, dist))
np.savez("axis4_search.npz", f=np.array([o[0] for o in out]), a=np.array([o[1] for o in out]), c=np.array([o[2] for o in out]), cen=cen)
# PCA-normal axis for comparison: only the centre free
from scipy.optimize import minimize
a_pc = n_up
th = np.arccos(a_pc[2]); ph = np.arctan2(a_pc[1], a_pc[0])
best = (1e9, None)
for u0 in (-150, -80, 0, 80, 150):
    for v0 in (-150, -80, 0, 80, 150):
        r = minimize(lambda uv: S.objective(np.array([th, ph, uv[0], uv[1]]), Pu[sel], cen), [u0, v0], method="Nelder-Mead", options=dict(xatol=1e-2, fatol=1e-5, maxiter=200))
        if r.fun < best[0]: best = (r.fun, r.x)
print("\naxis fixed to the plate normal, centre optimised: profile scatter %.3f mm" % best[0])
