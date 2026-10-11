import numpy as np
from scipy.spatial import ConvexHull
from scipy.optimize import least_squares
import sample2_reconstruction as S
import sample3_base_reconstruction as B
sp = np.array([0.9765625,0.9765625,1.0])
d = np.load("/mnt/user-data/uploads/sample1.npz"); origin = np.array([d["x"].min(), d["y"].min(), d["z"].min()], float)
clean = np.load("brown_clean.npy"); fr = np.load("brown_frame.npz"); a0 = fr["a"]; c0 = fr["c"]
V, F, N, Wt = S.mesh_from_field(clean.astype(float), 0.5, sp, origin, sigma=1.0)
# side B = the side whose plateau protrudes: try both signs, keep the one with the larger flat cap
best = None
for sgn in (+1, -1):
    ax = sgn*a0; h = (V-c0)@ax
    m = (h > np.percentile(h, 0.5)) & (h < np.percentile(h, 0.5) + 5.0) & (N@ax < -0.6)   # lowest 5 mm on that side, facing outward
    if m.sum() < 200: continue
    P = V[m]; mu = P.mean(0); Vp = np.linalg.svd(P-mu, full_matrices=False)[2]; nrm = Vp[2]
    rms = np.sqrt(np.mean(((P-mu)@nrm)**2)); ext = np.sqrt(np.linalg.svd(P-mu, full_matrices=False)[1][:2]**2/len(P))*np.sqrt(12)
    print("sign %+d: %d lowest-side vertices, plane RMS %.2f mm, extents %.0f x %.0f mm, normal %s" % (sgn, m.sum(), rms, ext[0], ext[1], np.round(nrm,3)))
    if best is None or rms < best[0]: best = (rms, sgn, m, mu, nrm)
rms, sgn, m, mu, nrm = best
if nrm @ (sgn*a0) > 0: nrm = -nrm            # normal pointing out of the foot (down)
axis_up = -nrm                                # axis direction pointing from the foot into the piece
P = V[m]; e1, e2 = S.basis(axis_up)
uv = np.column_stack([(P-mu)@e1, (P-mu)@e2]); hull = ConvexHull(uv); bp = uv[hull.vertices]
def fit(p): return np.hypot(bp[:,0]-p[0], bp[:,1]-p[1]) - p[2]
s = least_squares(fit, [bp[:,0].mean(), bp[:,1].mean(), 40.0]); cu, cv, R = s.x
print("\nfoot (lowest cap): circle through its outline: radius %.1f mm (diameter %.0f mm), RMS %.2f mm; foot plane normal vs the PCA axis: %.1f deg" % (R, 2*R, np.sqrt(np.mean(s.fun**2)), np.degrees(np.arccos(abs(axis_up@a0)))))
centre = mu + cu*e1 + cv*e2
np.savez("foot1.npz", axis_up=axis_up, centre=centre, R=R, mu=mu)
# profile about this axis
rng = np.random.default_rng(1)
Q = V - centre; h = Q@axis_up; ex, ey = S.basis(axis_up); rho = np.hypot(Q@ex, Q@ey)
print("height range %.0f..%.0f mm above the foot plane origin; rho max %.0f mm" % (h.min(), h.max(), rho.max()))
Nf = N; up_f = (Nf@axis_up) > 0.35; lo_f = (Nf@axis_up) < -0.35
Pu, Pl = V[up_f], V[lo_f]; cen = np.vstack([Pu, Pl]).mean(0)
Pu_s = Pu[rng.choice(len(Pu), min(5000,len(Pu)), replace=False)]; Pl_s = Pl[rng.choice(len(Pl), min(5000,len(Pl)), replace=False)]
th = np.arccos(axis_up[2]); ph = np.arctan2(axis_up[1], axis_up[0]); e1c, e2c = S.basis(axis_up)
p0 = np.array([th, ph, (centre-cen)@e1c, (centre-cen)@e2c])
print("two-face spread with the axis fixed by the foot: %.2f mm" % B.two_face_objective(p0, Pu_s, Pl_s, cen))
from scipy.optimize import minimize
r_ = minimize(B.two_face_objective, p0, args=(Pu_s, Pl_s, cen), method="Nelder-Mead", options=dict(xatol=1e-3, fatol=1e-5, maxiter=400))
a_r, c_r = S.axis_from(r_.x, cen)
if a_r @ axis_up < 0: a_r = -a_r
print("refined on both faces: spread %.2f mm, turned by %.1f deg, moved by %.1f mm" % (r_.fun, np.degrees(np.arccos(np.clip(a_r@axis_up,-1,1))), np.linalg.norm(np.cross(c_r-centre, a_r))))
