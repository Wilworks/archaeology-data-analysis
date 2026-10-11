import numpy as np
from collections import defaultdict, deque
import sample2_reconstruction as S
sp = np.array([0.9765625,0.9765625,1.0])
d = np.load("/mnt/user-data/uploads/sample1.npz"); origin = np.array([d["x"].min(), d["y"].min(), d["z"].min()], float)
clean = np.load("brown_clean.npy")
V, F, N, Wt = S.mesh_from_field(clean.astype(float), 0.5, sp, origin, sigma=1.5)
tri = V[F]; fn = np.cross(tri[:,1]-tri[:,0], tri[:,2]-tri[:,0]); fa = 0.5*np.linalg.norm(fn,axis=1); fn = fn/(2*fa[:,None]+1e-12); fc = tri.mean(1)
em = defaultdict(list)
for i,f in enumerate(F):
    for k in range(3): em[tuple(sorted((f[k], f[(k+1)%3])))].append(i)
adj = defaultdict(list)
for e,fs in em.items():
    if len(fs)==2: adj[fs[0]].append(fs[1]); adj[fs[1]].append(fs[0])
for tol_deg in (6, 10):
    tol = np.cos(np.radians(tol_deg)); lab = -np.ones(len(F), int); regs = []
    for seed in np.argsort(-fa):
        if lab[seed] >= 0: continue
        rid = len(regs); q = deque([seed]); lab[seed] = rid; ns = fn[seed]*fa[seed]; mem = [seed]
        while q:
            i = q.popleft()
            for j in adj[i]:
                if lab[j] >= 0: continue
                nm = ns/np.linalg.norm(ns)
                if fn[j]@nm > tol: lab[j] = rid; q.append(j); mem.append(j); ns += fn[j]*fa[j]
        regs.append(np.array(mem))
    rows = []
    for idx in regs:
        A = fa[idx].sum()
        if A < 1500: continue
        Xp = fc[idx]; w = fa[idx]; mu = (Xp*w[:,None]).sum(0)/w.sum()
        _, Sg, Vp = np.linalg.svd((Xp-mu)*np.sqrt(w)[:,None], full_matrices=False)
        nrm = Vp[2]; rms = np.sqrt((w*((Xp-mu)@nrm)**2).sum()/w.sum())
        rows.append((A, rms, nrm, mu, np.sqrt(Sg[:2]**2/w.sum())*np.sqrt(12), idx))
    rows.sort(key=lambda r: -r[0])
    print("\nflat patches at %d deg tolerance: area mm2, plane RMS mm, extents mm, centre, normal" % tol_deg)
    for A, rms, nrm, mu, ext, _ in rows[:6]: print("  %7.0f  %.2f  %3.0f x %3.0f  %s  %s" % (A, rms, ext[0], ext[1], np.round(mu,0), np.round(nrm,3)))
    if tol_deg == 10:
        np.savez("plateau_patches.npz", n=np.array([r[2] for r in rows[:4]]), mu=np.array([r[3] for r in rows[:4]]), A=[r[0] for r in rows[:4]])
