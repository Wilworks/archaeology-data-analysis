import numpy as np
from scipy import ndimage as ndi
from collections import defaultdict, deque
import sample2_reconstruction as S
HU, sp, origin, hu = S.load_grid("/mnt/user-data/uploads/sample3.npz")
present = HU >= 500; core = np.load("core3.npy")
V, F, N, W = S.mesh_from_field(np.where(present, HU.astype(float), 0.0), 500.0, sp, origin, sigma=0.8)
# drop small crust flecks: keep mesh vertices that lie within 1.5 mm of the ceramic core
from scipy.spatial import cKDTree
Xc = np.argwhere(core)*sp+origin
near = cKDTree(Xc).query(V)[0] < 1.6
print("mesh vertices: %d, kept near the ceramic core: %d" % (len(V), near.sum()))
tri = V[F]; fn = np.cross(tri[:,1]-tri[:,0], tri[:,2]-tri[:,0]); fa = 0.5*np.linalg.norm(fn,axis=1); fn = fn/(2*fa[:,None]+1e-12); fc = tri.mean(1)
keepf = near[F].all(1)
em = defaultdict(list)
for i,f in enumerate(F):
    if not keepf[i]: continue
    for k in range(3): em[tuple(sorted((f[k], f[(k+1)%3])))].append(i)
adj = defaultdict(list)
for e,fs in em.items():
    if len(fs)==2: adj[fs[0]].append(fs[1]); adj[fs[1]].append(fs[0])
tol = np.cos(np.radians(4)); lab = -np.ones(len(F), int); regs = []
for seed in np.argsort(-fa):
    if not keepf[seed] or lab[seed]>=0: continue
    rid=len(regs); q=deque([seed]); lab[seed]=rid; ns=fn[seed]*fa[seed]; mem=[seed]
    while q:
        i=q.popleft()
        for j in adj[i]:
            if lab[j]>=0: continue
            nm=ns/np.linalg.norm(ns)
            if fn[j]@nm>tol: lab[j]=rid; q.append(j); mem.append(j); ns+=fn[j]*fa[j]
    regs.append(np.array(mem))
rows=[]
for idx in regs:
    A=fa[idx].sum()
    if A<300: continue
    Xp=fc[idx]; w=fa[idx]; mu=(Xp*w[:,None]).sum(0)/w.sum()
    _,Sg,Vp=np.linalg.svd((Xp-mu)*np.sqrt(w)[:,None], full_matrices=False)
    nrm=Vp[2]; rms=np.sqrt((w*((Xp-mu)@nrm)**2).sum()/w.sum())
    rows.append((A, rms, nrm, mu, np.sqrt(Sg[:2]**2/w.sum())*np.sqrt(12), idx))
rows.sort(key=lambda r:-r[0])
print("largest flat patches (4 deg tolerance): area mm2, plane RMS mm, extents, centre, normal")
for A,rms,nrm,mu,ext,_ in rows[:8]: print("  %6.0f  %.2f  %3.0f x %3.0f  %s  %s" % (A, rms, ext[0], ext[1], np.round(mu,0), np.round(nrm,3)))
# relations between the two biggest flat patches
if len(rows)>=2:
    n1, n2 = rows[0][2], rows[1][2]
    print("angle between the two largest flat patches: %.1f deg; distance between their planes along n1: %.1f mm" % (np.degrees(np.arccos(min(1,abs(n1@n2)))), abs((rows[1][3]-rows[0][3])@n1)))
np.savez("flat3.npz", n1=rows[0][2], mu1=rows[0][3], n2=rows[1][2], mu2=rows[1][3], idx1=rows[0][5], idx2=rows[1][5])
np.savez("mesh3.npz", V=V, F=F, N=N, W=W, near=near)
