import numpy as np
from scipy import ndimage as ndi
from scipy.optimize import minimize
import sample2_reconstruction as S
import sample3_base_reconstruction as B
sp = np.array([0.9765625,0.9765625,1.0])
d = np.load("/mnt/user-data/uploads/sample1.npz"); origin = np.array([d["x"].min(), d["y"].min(), d["z"].min()], float)
clean = np.load("brown_clean.npy"); fr = np.load("brown_frame.npz"); n3 = fr["a"]; c0 = fr["c"]
V, F, N, Wt = S.mesh_from_field(clean.astype(float), 0.5, sp, origin, sigma=1.0)
print("cleaned mesh: %d vertices" % len(V))
up = (N @ n3) > 0.35; lo = (N @ n3) < -0.35
print("faces: upper %d, lower %d vertices" % (up.sum(), lo.sum()))
rng = np.random.default_rng(3)
def samp(m, n=5000):
    i = np.where(m)[0]; return V[rng.choice(i, min(n, len(i)), replace=False)]
Pu, Pl = samp(up), samp(lo); cen = np.vstack([Pu, Pl]).mean(0)
# global direction search with a few centre offsets
dirs = S.fib_hemisphere(400); cand = []
for a in dirs:
    th = np.arccos(np.clip(a[2],-1,1)); ph = np.arctan2(a[1], a[0]); e1, e2 = S.basis(S.dir_from(th, ph))
    best = min(((B.two_face_objective(np.array([th, ph, du, dv]), Pu, Pl, cen), (du, dv)) for du in (-15, 0, 15) for dv in (-15, 0, 15)), key=lambda t: t[0])
    cand.append((best[0], np.array([th, ph, best[1][0], best[1][1]])))
cand.sort(key=lambda t: t[0])
picked = []
for f, p in cand:
    if all(np.degrees(np.arccos(min(1, abs(S.dir_from(p[0],p[1]) @ S.dir_from(q[0],q[1]))))) > 8 for _, q in picked): picked.append((f, p))
    if len(picked) >= 6: break
res = []
for f, p in picked:
    r = minimize(B.two_face_objective, p, args=(Pu, Pl, cen), method="Nelder-Mead", options=dict(xatol=5e-4, fatol=1e-5, maxiter=500))
    a, c = S.axis_from(r.x, cen); res.append((r.fun, a, c))
res.sort(key=lambda t: t[0])
print("\nrank  spread(mm)  angle to PCA axis  distance from the piece centroid (mm)  direction")
for k, (f, a, c) in enumerate(res):
    if a @ n3 < 0: a = -a
    print("%3d   %.3f       %5.1f              %6.1f                             %s" % (k+1, f, np.degrees(np.arccos(np.clip(a@n3,-1,1))), np.linalg.norm(np.cross(V.mean(0)-c, a)), np.round(a,3)))
print("\nfor comparison, the PCA axis through the centroid: spread %.3f mm" % B.two_face_objective(np.array([np.arccos(n3[2]), np.arctan2(n3[1], n3[0]), 0.0, 0.0]) , Pu, Pl, V.mean(0)))
np.savez("axis_brown.npz", f=[r[0] for r in res], a=np.array([r[1] for r in res]), c=np.array([r[2] for r in res]))
