import numpy as np
from scipy import ndimage as ndi
from scipy.spatial import cKDTree
from scipy.optimize import minimize
import sample2_reconstruction as S
import sample4_wall_reconstruction as W
HU = np.load("grid1.npy"); sp = np.array([0.9765625, 0.9765625, 1.0])
d = np.load("/mnt/user-data/uploads/sample1.npz"); origin = np.array([d["x"].min(), d["y"].min(), d["z"].min()], float)
seg = np.load("seg1.npy"); body_all = np.load("body1.npy")
main = 1 + int(np.argmax([(seg==j).sum() for j in range(1, seg.max()+1)]))
m = seg == main
present = HU >= 500
crust = present & ~ndi.binary_dilation(body_all, iterations=2)
Vb, Fb, Nb, Wb = S.mesh_from_field(m.astype(float), 0.5, sp, origin, sigma=1.0)
Xb = np.argwhere(m)*sp + origin; cb = Xb.mean(0)
Xc = np.argwhere(crust)*sp + origin
n3 = np.linalg.svd(Xb-cb, full_matrices=False)[2][2]
# crust side of THIS piece: crust voxels within 25 mm of it
near_c = cKDTree(Xb).query(Xc)[0] < 12
n_up = n3*np.sign((Xc[near_c].mean(0)-cb)@n3)
print("main piece: %d vertices; crust within 12 mm of it: %d voxels, centroid offset along the plate normal %.1f mm" % (len(Vb), near_c.sum(), (Xc[near_c].mean(0)-cb)@n_up))
under = (Nb @ n_up) < -0.35
print("crust-free-side vertices: %d (%.0f%% of area)" % (under.sum(), 100*Wb[under].sum()/Wb.sum()))
rng = np.random.default_rng(21)
Pu, Nu, Wu = Vb[under], Nb[under], Wb[under]
sel = rng.choice(len(Pu), size=2200, replace=False, p=Wu/Wu.sum())
sols, cen = W.wall_axis_search(Pu[sel], Nu[sel], n_up, min_angle=0.0, n_dirs=1200, keep=10)
print("\nrank  spread(mm)  angle to plate normal  distance axis-to-piece(mm)  direction")
big = Pu[rng.choice(len(Pu), size=min(6000, len(Pu)), replace=False, p=Wu/Wu.sum())]
out = []
for f, p in sols[:8]:
    r = minimize(W.wall_objective, p, args=(big, cen), method="Nelder-Mead", options=dict(xatol=5e-4, fatol=1e-6, maxiter=600))
    a, c = S.axis_from(r.x, cen)
    if a @ n_up < 0: a = -a
    out.append((r.fun, a, c))
out.sort(key=lambda t: t[0])
for k, (f, a, c) in enumerate(out):
    print("%3d   %.3f       %5.1f                  %6.1f                  %s" % (k+1, f, np.degrees(np.arccos(np.clip(a@n_up,-1,1))), np.linalg.norm(np.cross(Vb.mean(0)-c, a)), np.round(a,3)))
np.savez("axis1_search.npz", f=[o[0] for o in out], a=np.array([o[1] for o in out]), c=np.array([o[2] for o in out]), n_up=n_up)
