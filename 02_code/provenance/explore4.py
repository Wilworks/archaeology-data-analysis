import numpy as np, json
from scipy import ndimage as ndi
from collections import defaultdict, deque
import sample2_reconstruction as S

HU, sp, origin, hu = S.load_grid("/mnt/user-data/uploads/sample4.npz")
present = HU >= 500
F = np.where(present, HU.astype(float), 0.0)
Fs = ndi.gaussian_filter(F, 1.2); w = ndi.gaussian_filter(present.astype(float), 1.2)
Fn = np.where(w > 0.5, Fs / np.maximum(w, 1e-6), 0.0)
T = 1100.0
body = ndi.binary_fill_holes(S.largest_component((Fn >= T) & present))
crust = present & ~ndi.binary_dilation(body, iterations=2)
vox = float(np.prod(sp))
print("body voxels %d (%.1f cm3), crust-class %d (%.1f cm3); body share of present %.1f%%" %
      (body.sum(), body.sum()*vox/1000, crust.sum(), crust.sum()*vox/1000, 100*body.sum()/present.sum()))
Vb, Fb, Nb, Wb = S.mesh_from_field(body.astype(float), 0.5, sp, origin, sigma=1.0)
print("body mesh: %d vertices, surface %.0f mm2" % (len(Vb), Wb.sum()))
np.savez("body4_mesh.npz", V=Vb, F=Fb, N=Nb, W=Wb)
np.save("body4.npy", body); np.save("crust4.npy", crust)

# ---- plate geometry ----
Xb = np.argwhere(body) * sp + origin; cb = Xb.mean(0)
U_, S_, Vt = np.linalg.svd(Xb - cb, full_matrices=False)
sd = S_ / np.sqrt(len(Xb)); Q = (Xb - cb) @ Vt.T
print("PCA std along PC1,PC2,PC3: %s mm ; extents %s mm" % (np.round(sd, 1), np.round(Q.max(0) - Q.min(0), 1)))
Xc = np.argwhere(crust) * sp + origin; cc = Xc.mean(0)
n3 = Vt[2]; sgn = np.sign((cc - cb) @ n3); n_up = n3 * sgn
print("crust centroid offset from the body centroid along the plate normal: %.1f mm (crust is on the +n_up side)" % ((cc - cb) @ n_up))
print("PC directions (scanner xyz): PC1 %s PC2 %s PC3(n_up) %s" % (np.round(Vt[0], 3), np.round(Vt[1], 3), np.round(n_up, 3)))
Qs = Q.copy(); Qs[:, 2] *= sgn
A = np.column_stack([Qs[:, 0] ** 2, Qs[:, 1] ** 2, Qs[:, 0] * Qs[:, 1], Qs[:, 0], Qs[:, 1], np.ones(len(Qs))])
coef, *_ = np.linalg.lstsq(A, Qs[:, 2], rcond=None)
print("mid-surface sag: u3 = %.5f u1^2 + %.5f u2^2 + %.5f u1u2 (+ linear terms) -> curvature radii %.0f mm (along PC1), %.0f mm (along PC2)" %
      (coef[0], coef[1], coef[2], 1 / abs(2 * coef[0]), 1 / abs(2 * coef[1])))
# thickness of the plate
edt = ndi.distance_transform_edt(body, sampling=tuple(sp))
print("plate local thickness (2 x EDT): median of the medial-axis-like voxels %.1f mm, 95th pct %.1f, max %.1f" %
      (2 * np.percentile(edt[body][edt[body] > 1.5], 60), 2 * np.percentile(edt[body], 95), 2 * edt.max()))

# ---- flat patches on the body mesh ----
tri = Vb[Fb]; fn = np.cross(tri[:, 1] - tri[:, 0], tri[:, 2] - tri[:, 0]); farea = 0.5 * np.linalg.norm(fn, axis=1)
fn = fn / (2 * farea[:, None] + 1e-12); fc = tri.mean(1)
em = defaultdict(list)
for i, f in enumerate(Fb):
    for k in range(3): em[tuple(sorted((f[k], f[(k + 1) % 3])))].append(i)
adj = [[] for _ in range(len(Fb))]
for e, fs in em.items():
    if len(fs) == 2: adj[fs[0]].append(fs[1]); adj[fs[1]].append(fs[0])
tol = np.cos(np.radians(7)); label = -np.ones(len(Fb), int); regions = []
for seed in np.argsort(-farea):
    if label[seed] >= 0: continue
    rid = len(regions); q = deque([seed]); label[seed] = rid; nsum = fn[seed] * farea[seed]; mem = [seed]
    while q:
        i = q.popleft()
        for j in adj[i]:
            if label[j] >= 0: continue
            nm = nsum / np.linalg.norm(nsum)
            if fn[j] @ nm > tol: label[j] = rid; q.append(j); mem.append(j); nsum += fn[j] * farea[j]
    regions.append(np.array(mem))
print("\nlargest nearly-flat patches (7 deg tolerance): area mm2, plane RMS mm, normal.n_up, extents mm")
rows = []
for idx in regions:
    A_ = farea[idx].sum()
    if A_ < 150: continue
    Xp = fc[idx]; wp = farea[idx]; mu = (Xp * wp[:, None]).sum(0) / wp.sum()
    _, Sg, Vp = np.linalg.svd((Xp - mu) * np.sqrt(wp)[:, None], full_matrices=False)
    nrm = Vp[2]; rms = np.sqrt((wp * ((Xp - mu) @ nrm) ** 2).sum() / wp.sum())
    rows.append((A_, rms, nrm, np.sqrt(Sg[:2] ** 2 / wp.sum()) * np.sqrt(12), mu))
rows.sort(key=lambda r: -r[0])
for A_, rms, nrm, ext, mu in rows[:10]:
    print("  %7.0f   %.2f   %+.2f   %4.0f x %3.0f   centre %s" % (A_, rms, nrm @ n_up, ext[0], ext[1], np.round(mu, 0)))
np.savez("s4_frame.npz", cb=cb, Vt=Vt, n_up=n_up, sgn=sgn)
