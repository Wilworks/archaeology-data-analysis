import numpy as np, json, trimesh
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.collections import PolyCollection
from scipy.spatial import cKDTree
import sample2_reconstruction as S
HU, sp, origin, hu = S.load_grid("/mnt/user-data/uploads/sample3.npz")
present = HU >= 500; core = S.largest_component(HU >= 1200)
Vo, Fo, No, Wo = S.mesh_from_field(np.where(present, HU.astype(float), 0.0), 500.0, sp, origin, sigma=0.8)
Xc = np.argwhere(core)*sp + origin
near = cKDTree(Xc).query(Vo)[0] < 1.6
rec = trimesh.load("/mnt/user-data/outputs/sample3_base_scannerLPS.stl")
R = json.load(open("out3d/results.json")); a = np.array(R["axis"]["direction"]); c = np.array(R["axis"]["point_mm"])
e1, e2 = S.centred_basis(a, c, Vo[near])
# ---- distance of the scan surface from the reconstruction (signed: + = the scan sticks out of the reconstruction) ----
idx = np.where(near)[0]; rng = np.random.default_rng(0); sel = idx[rng.choice(len(idx), 14000, replace=False)]
cp, dist, tri = trimesh.proximity.closest_point(rec, Vo[sel])
sgn = np.where(np.einsum("ij,ij->i", Vo[sel]-cp, rec.face_normals[tri]) > 0, 1.0, -1.0)
sd = sgn*dist
print("scan surface vs reconstruction: median |distance| %.2f mm; within 1 mm: %.0f%%; within 2 mm: %.0f%%; within 3 mm: %.0f%%" % (np.median(dist), 100*(dist<1).mean(), 100*(dist<2).mean(), 100*(dist<3).mean()))
print("signed distance: median %+.2f mm, 5th-95th percentile %+.1f to %+.1f mm" % (np.median(sd), *np.percentile(sd,[5,95])))
q = Vo[sel]-c; h = q@a; x = q@e1; y = q@e2; rho = np.hypot(x,y)
print("\nmean |distance| and mean signed distance by distance from the axis:")
for lo_, hi_ in ((0,20),(20,40),(40,50),(50,58),(58,64),(64,72)):
    m = (rho>=lo_)&(rho<hi_)
    if m.sum()>50: print("  rho %2d-%2d mm: %5d points, mean |d| %.2f mm, mean signed %+.2f mm, 90th pct |d| %.2f" % (lo_, hi_, m.sum(), dist[m].mean(), sd[m].mean(), np.percentile(dist[m],90)))
# ---- renderer ----
def tri_arrays(V, F): 
    T3 = V[F]; n = np.cross(T3[:,1]-T3[:,0], T3[:,2]-T3[:,0]); n /= (np.linalg.norm(n,axis=1,keepdims=True)+1e-12); return T3, n
def render(ax, T3, nrm, cols, cam, upref, title):
    cam = cam/np.linalg.norm(cam); right = np.cross(upref, cam); right /= np.linalg.norm(right); up = np.cross(cam, right)
    Q = T3 - c; sx, sy, dz = Q@right, Q@up, Q@cam
    order = np.argsort(dz.mean(1)); lv = 0.6*cam + 0.4*np.array([0.3,-0.4,0.85]); lv /= np.linalg.norm(lv)
    shade = 0.4 + 0.6*np.clip(np.abs(nrm@lv),0,1)
    ax.add_collection(PolyCollection(np.stack([sx,sy],-1)[order], facecolors=np.clip(cols*shade[:,None],0,1)[order], edgecolors="none"))
    ax.set_xlim(sx.min()-3, sx.max()+3); ax.set_ylim(sy.min()-3, sy.max()+3); ax.set_aspect("equal"); ax.axis("off"); ax.set_title(title, fontsize=9.5)
views = [("looking down the axis (from the top)", a, e1), ("side view, axis vertical", e1, a), ("oblique from above", 0.55*e1+0.4*e2+0.75*a, a)]
Ts, ns = tri_arrays(Vo, Fo)
fk = near[Fo].all(1)
col_scan = np.where(fk[:,None], np.array([0.88,0.82,0.74]), np.array([0.68,0.43,0.20]))
Tr, nr = tri_arrays(rec.vertices, rec.faces)
col_rec = np.tile(np.array([0.62,0.71,0.80]), (len(Tr),1))
# scan coloured by distance to the reconstruction: give every scan vertex a distance (nearest sampled vertex)
tree = cKDTree(Vo[sel]); dv = tree.query(Vo)[0]; _, nn_i = tree.query(Vo); sdv = sd[nn_i]
sdf = sdv[Fo].mean(1); cm = plt.get_cmap("coolwarm"); col_d = cm(np.clip((sdf+4)/8, 0, 1))[:, :3]
col_d[~fk] = np.array([0.75,0.75,0.75])
fig, axs = plt.subplots(3, 3, figsize=(16, 15))
for j, (t, cam, upr) in enumerate(views):
    render(axs[0,j], Ts, ns, col_scan, cam, upr, "SCAN, "+t+"\n(brown = crust flecks)")
    render(axs[1,j], Tr, nr, col_rec, cam, upr, "RECONSTRUCTION, "+t)
    render(axs[2,j], Ts, ns, col_d, cam, upr, "SCAN coloured by distance from the reconstruction\n(red = sticks out, blue = inside, grey = crust)")
sm = plt.cm.ScalarMappable(cmap=cm, norm=plt.Normalize(-4,4)); fig.colorbar(sm, ax=axs[2,:], fraction=0.02, label="mm")
fig.suptitle("Scan versus reconstruction, same camera", fontsize=12, fontweight="bold")
fig.subplots_adjust(left=0.01, right=0.93, top=0.94, bottom=0.01, wspace=0.02, hspace=0.12)
fig.savefig("compare3.png", dpi=70)
np.savez("compare3_stats.npz", dist=dist, sd=sd, rho=rho)
