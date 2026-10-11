import numpy as np
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.collections import PolyCollection
from scipy import ndimage as ndi
from skimage import measure
import sample2_reconstruction as S
import sample4_wall_reconstruction as W
HU = np.load("grid1.npy"); sp = np.array([0.9765625, 0.9765625, 1.0])
d = np.load("/mnt/user-data/uploads/sample1.npz"); origin = np.array([d["x"].min(), d["y"].min(), d["z"].min()], float)
present, body, T, p_crust, p_body = W.smoothed_body(HU, sp)
crust = present & ~ndi.binary_dilation(body, iterations=2)
vox = float(np.prod(sp))
print("split at %.0f HU (peaks %.0f / %.0f): body %.1f cm3 (%.0f%% of stored), crust-class %.1f cm3" % (T, p_crust, p_body, body.sum()*vox/1000, 100*body.sum()/present.sum(), (present & ~body).sum()*vox/1000))
np.save("body1.npy", body)
Xb = np.argwhere(body)*sp + origin; cb = Xb.mean(0)
U_, Sg, Vt = np.linalg.svd(Xb-cb, full_matrices=False)
if np.linalg.det(Vt) < 0: Vt[2] *= -1
Q = (Xb-cb)@Vt.T; sd = Sg/np.sqrt(len(Xb))
print("body PCA std %s mm; extents %s mm" % (np.round(sd,1), np.round(Q.max(0)-Q.min(0),1)))
Xc = np.argwhere(crust)*sp + origin
n3 = Vt[2]; sgn = np.sign((Xc.mean(0)-cb)@n3); print("crust centroid offset along the plate normal: %.1f mm" % ((Xc.mean(0)-cb)@(n3*sgn)))
A = np.column_stack([Q[:,0]**2, Q[:,1]**2, Q[:,0]*Q[:,1], Q[:,0], Q[:,1], np.ones(len(Q))]); coef,*_ = np.linalg.lstsq(A, Q[:,2]*sgn, rcond=None)
print("mid-surface sag u3 = %.5f u1^2 %+.5f u2^2 %+.5f u1u2 -> curvature radii %.0f mm (PC1), %.0f mm (PC2)" % (coef[0], coef[1], coef[2], 1/abs(2*coef[0]), 1/abs(2*coef[1])))
edt = ndi.distance_transform_edt(body, sampling=tuple(sp))
print("body thickness (2 x EDT): 60th pct of interior %.1f mm, 95th %.1f, max %.1f" % (2*np.percentile(edt[body][edt[body]>1.5],60), 2*np.percentile(edt[body],95), 2*edt.max()))
# ---- decimated surface for the views ----
F = np.where(present, HU.astype(float), 0.0); Fp = np.pad(ndi.gaussian_filter(F, 0.8), 3)
v, f, n, _ = measure.marching_cubes(Fp, level=500.0, spacing=tuple(sp), step_size=2)
v = v - 3*sp + origin
tri = f; cen = v[tri].mean(1)
bd = ndi.binary_dilation(body, iterations=2)
ijk = np.clip(np.rint((cen-origin)/sp).astype(int), 0, np.array(HU.shape)-1)
is_body = bd[ijk[:,0], ijk[:,1], ijk[:,2]]
P = (v-cb)@Vt.T; P[:,2] *= sgn
T3 = P[tri]; nrm = np.cross(T3[:,1]-T3[:,0], T3[:,2]-T3[:,0]); nrm /= (np.linalg.norm(nrm,axis=1,keepdims=True)+1e-12)
ctr = (P.max(0)+P.min(0))/2
base_col = np.where(is_body[:,None], np.array([0.88,0.82,0.74]), np.array([0.68,0.43,0.20]))
def render(ax, elev, azim, title):
    t_, p_ = np.radians(elev), np.radians(azim); cam = np.array([np.cos(t_)*np.cos(p_), np.cos(t_)*np.sin(p_), np.sin(t_)])
    right = np.cross([0,0,1.0], cam); right = right/np.linalg.norm(right) if np.linalg.norm(right)>1e-6 else np.array([1.0,0,0]); up = np.cross(cam, right)
    Qv = T3-ctr; sx, sy, dz = Qv@right, Qv@up, Qv@cam
    order = np.argsort(dz.mean(1)); lv = 0.6*cam+0.4*np.array([0.3,-0.4,0.85]); lv /= np.linalg.norm(lv)
    shade = 0.35+0.65*np.clip(np.abs(nrm@lv),0,1)
    ax.add_collection(PolyCollection(np.stack([sx,sy],-1)[order], facecolors=np.clip(base_col*shade[:,None],0,1)[order], edgecolors="none"))
    ax.set_xlim(sx.min()-3, sx.max()+3); ax.set_ylim(sy.min()-3, sy.max()+3); ax.set_aspect("equal"); ax.axis("off"); ax.set_title(title, fontsize=10)
fig, axs = plt.subplots(2,3, figsize=(17,10))
for ax,(t,el,az) in zip(axs.ravel(), [("Looking onto the crust side",90,-90),("Looking at the crust-free side",-90,-90),("Side view",0,-90),("Edge-on",0,0),("Oblique, crust side",35,-60),("Oblique, crust-free side",-35,30)]): render(ax, el, az, t)
fig.suptitle("Sample 1 (scan label 2) as stored in sample1.npz (500 HU surface). Pale = ceramic body, brown = crust", fontsize=11)
fig.tight_layout(rect=(0,0,1,0.95)); fig.savefig("s1_surface_views.png", dpi=70)
print("view mesh triangles:", len(tri))
