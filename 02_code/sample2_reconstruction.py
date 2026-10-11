#!/usr/bin/env python3
"""
Sample 2 (Fort Amsterdam ceramics): virtual reconstruction of a rotationally symmetric vessel
from one CT-scanned sherd.

Input  : sample2.npz  (x, y, z in mm; hu as int16; every voxel of the sherd at >= 500 HU)
Output : results.json, profile_sample2.csv, fig_sweep.png, fig_profile.png, fig_meridian.png, fig_lip_check.png, and the CAD/mesh
         exports (*.stl, *.dxf, *.step; DXF needs ezdxf, STEP needs cadquery, both optional)  (in --out)

Stages
  1  voxel grid, HU histogram valley -> body/crust split, body mesh
  2  axis of rotation: (a) global search that makes the (radius, height) profile of the
     smooth underside as thin as possible, (b) independent check/constraint from the flat,
     finished lip (plane fit; its normal is parallel to the axis)
  3  profile (outer + inner) in the vessel frame, wall thickness, half-maximum edge correction
  4  sweep the profile 360 degrees, mark the part that survives, render figures

Requirements: numpy, scipy, scikit-image, matplotlib      (tested: Python 3.12, numpy 2.4.4, scipy 1.17.1,
              scikit-image 0.26.0, matplotlib 3.10.8; run time a few minutes on a laptop)
Run         : python sample2_reconstruction.py --npz sample2.npz --out out
All random choices use fixed seeds, so repeated runs give the same numbers.
"""
import argparse, json, os
import numpy as np
from scipy import ndimage as ndi
from scipy.linalg import eigh
from scipy.optimize import minimize
from scipy.spatial import cKDTree
from scipy.interpolate import splprep, splev
from scipy.signal import medfilt
from scipy.sparse import coo_matrix
from scipy.sparse.csgraph import connected_components
from skimage import measure
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d.art3d import Poly3DCollection

# ----------------------------------------------------------------------------------------
# 0. small helpers
# ----------------------------------------------------------------------------------------
def basis(a):
    """Two unit vectors perpendicular to a (and to each other)."""
    ref = np.array([1.0, 0, 0]) if abs(a[0]) < 0.9 else np.array([0, 1.0, 0])
    e1 = np.cross(a, ref); e1 /= np.linalg.norm(e1)
    return e1, np.cross(a, e1)

def centred_basis(a, c, pts):
    """Basis (e1, e2) perpendicular to a, rotated so the mean azimuth of pts about the axis is 0."""
    e1, e2 = basis(a)
    q = pts - c
    phi0 = np.angle(np.mean(np.exp(1j * np.arctan2(q @ e2, q @ e1))))
    return np.cos(phi0) * e1 + np.sin(phi0) * e2, -np.sin(phi0) * e1 + np.cos(phi0) * e2

def dir_from(th, ph):
    return np.array([np.sin(th) * np.cos(ph), np.sin(th) * np.sin(ph), np.cos(th)])

def largest_component(mask):
    lab, n = ndi.label(mask, structure=np.ones((3, 3, 3)))
    sizes = np.bincount(lab.ravel())[1:]
    return lab == (np.argmax(sizes) + 1)

def largest_vertex_component(mask, faces, nv):
    """Largest connected patch (through mesh faces) of the vertices selected by mask."""
    f = faces[mask[faces].all(1)]
    A = coo_matrix((np.ones(3 * len(f)), (np.r_[f[:, 0], f[:, 1], f[:, 2]], np.r_[f[:, 1], f[:, 2], f[:, 0]])), shape=(nv, nv))
    _, lab = connected_components(A, directed=False)
    keep = np.bincount(lab[mask]).argmax()
    return mask & (lab == keep)

# ----------------------------------------------------------------------------------------
# 1. grid, threshold, meshes
# ----------------------------------------------------------------------------------------
def load_grid(path):
    d = np.load(path)
    x, y, z = (d[k].astype(float) for k in ("x", "y", "z"))
    hu = d["hu"].astype(int)
    def step(v):
        u = np.unique(np.round(v, 3)); return float(np.median(np.diff(u)))
    spacing = np.array([step(x), step(y), step(z)])
    origin = np.array([x.min(), y.min(), z.min()])
    idx = np.rint((np.column_stack([x, y, z]) - origin) / spacing).astype(int)
    shape = tuple(idx.max(0) + 1)
    HU = np.full(shape, -2000, dtype=np.int16)
    HU[idx[:, 0], idx[:, 1], idx[:, 2]] = hu
    return HU, spacing, origin, hu

def histogram_valley(hu, lo=1000, hi=1450):
    """Valley of the 25-HU histogram between the crust peak (~1000) and the body peak (~1330)."""
    edges = np.arange(500, hu.max() + 26, 25)
    h, e = np.histogram(hu, bins=edges); c = 0.5 * (e[:-1] + e[1:])
    hs = np.convolve(h, np.ones(3) / 3, mode="same")
    sel = (c >= lo) & (c <= hi)
    p1 = c[sel & (c < 1150)][np.argmax(hs[sel & (c < 1150)])]
    p2 = c[sel & (c > 1250)][np.argmax(hs[sel & (c > 1250)])]
    between = (c > p1) & (c < p2)
    return float(c[between][np.argmin(hs[between])]), float(p1), float(p2)

def mesh_from_field(F, level, spacing, origin, pad=4, sigma=1.0):
    Fp = np.pad(F.astype(float), pad, mode="constant", constant_values=float(F.min()))
    if sigma > 0:
        Fp = ndi.gaussian_filter(Fp, sigma)
    v, f, n, _ = measure.marching_cubes(Fp, level=level, spacing=tuple(spacing), gradient_direction="descent")
    v = v - pad * spacing + origin
    tri = v[f]
    area = 0.5 * np.linalg.norm(np.cross(tri[:, 1] - tri[:, 0], tri[:, 2] - tri[:, 0]), axis=1)
    wv = np.zeros(len(v))
    for k in range(3):
        np.add.at(wv, f[:, k], area / 3)
    n = n / np.linalg.norm(n, axis=1, keepdims=True)
    return v, f, n, wv

# ----------------------------------------------------------------------------------------
# 2. axis of rotation
# ----------------------------------------------------------------------------------------
DEG = 4
POWERS = [(i, j) for i in range(DEG + 1) for j in range(DEG + 1 - i)]
SC = 50.0

def _design(x, y):
    D = np.column_stack([x ** i * y ** j for i, j in POWERS])
    Dx = np.column_stack([i * x ** max(i - 1, 0) * y ** j if i > 0 else np.zeros_like(x) for i, j in POWERS])
    Dy = np.column_stack([j * x ** i * y ** max(j - 1, 0) if j > 0 else np.zeros_like(x) for i, j in POWERS])
    return D, Dx, Dy

def taubin_rms(rho, h):
    """RMS distance (mm) of the (rho, h) cloud from the best degree-4 implicit curve (Taubin fit).
    Small value  <=>  the cloud is a thin curve  <=>  the trial axis is right."""
    x = (rho - rho.mean()) / SC; y = (h - h.mean()) / SC
    D, Dx, Dy = _design(x, y); n = len(x)
    A = D.T @ D / n
    B = (Dx.T @ Dx + Dy.T @ Dy) / n + 1e-9 * np.eye(len(POWERS))
    w = eigh(A, B, eigvals_only=True, subset_by_index=[0, 0])
    return float(np.sqrt(max(w[0], 0.0)) * SC)

def axis_from(params, cen):
    th, ph, u, v = params
    a = dir_from(th, ph); e1, e2 = basis(a)
    return a, cen + u * e1 + v * e2

def profile_coords(pts, a, c):
    q = pts - c; h = q @ a
    return np.linalg.norm(q - h[:, None] * a, axis=1), h

def objective(params, pts, cen):
    a, c = axis_from(params, cen)
    return taubin_rms(*profile_coords(pts, a, c))

def centre_from_normals(a, pts, nrm):
    """Where the surface normals (projected perpendicular to a) meet: a good start for the axis position."""
    Np = nrm - (nrm @ a)[:, None] * a
    ok = np.linalg.norm(Np, axis=1) > 0.15
    t = np.cross(nrm[ok], a); t /= np.linalg.norm(t, axis=1, keepdims=True)
    M = t[:, :, None] * t[:, None, :]
    return np.linalg.pinv(M.sum(0), rcond=1e-6) @ np.einsum("nij,nj->i", M, pts[ok])

def fib_hemisphere(n):
    k = np.arange(n) + 0.5
    phi = np.arccos(1 - k / n); th = np.pi * (1 + 5 ** 0.5) * k
    return np.column_stack([np.sin(phi) * np.cos(th), np.sin(phi) * np.sin(th), np.cos(phi)])

def global_axis_search(pts, nrm, n_dirs=500, keep=12):
    cen = pts.mean(0); cand = []
    offs = [(0, 0)] + [(du, dv) for du in (-60, 60) for dv in (-60, 60)] + [(0, 90), (0, -90), (90, 0), (-90, 0)]
    for a in fib_hemisphere(n_dirs):
        th = np.arccos(np.clip(a[2], -1, 1)); ph = np.arctan2(a[1], a[0])
        e1, e2 = basis(dir_from(th, ph))
        c0 = centre_from_normals(a, pts - cen, nrm)
        u0, v0 = c0 @ e1, c0 @ e2
        best = min(((objective(np.array([th, ph, u0 + du, v0 + dv]), pts, cen), (du, dv)) for du, dv in offs),
                   key=lambda t: t[0])
        cand.append((best[0], np.array([th, ph, u0 + best[1][0], v0 + best[1][1]])))
    cand.sort(key=lambda t: t[0])
    picked = []
    for f, p in cand:
        if all(np.degrees(np.arccos(min(1, abs(dir_from(p[0], p[1]) @ dir_from(q[0], q[1]))))) > 8 for _, q in picked):
            picked.append((f, p))
        if len(picked) >= keep: break
    out = []
    for f, p in picked:
        r = minimize(objective, p, args=(pts, cen), method="Nelder-Mead",
                     options=dict(xatol=1e-3, fatol=1e-4, maxiter=300))
        out.append((r.fun, r.x))
    out.sort(key=lambda t: t[0])
    return out, cen

def fit_plane(X, w):
    mu = (X * w[:, None]).sum(0) / w.sum()
    _, _, Vt = np.linalg.svd((X - mu) * np.sqrt(w)[:, None], full_matrices=False)
    n = Vt[2]; d = (X - mu) @ n
    return mu, n, float(np.sqrt((w * d ** 2).sum() / w.sum()))

# ----------------------------------------------------------------------------------------
# 3. sweep rendering
# ----------------------------------------------------------------------------------------
LIGHT = np.array([-0.45, 0.35, 0.82]); LIGHT /= np.linalg.norm(LIGHT)
def _rgb(h): return np.array([int(h[i:i + 2], 16) / 255 for i in (1, 3, 5)])

def revolve(rr, zz, phis):
    return np.stack([rr[:, None] * np.cos(phis)[None, :], rr[:, None] * np.sin(phis)[None, :],
                     np.repeat(zz[:, None], len(phis), 1)], axis=-1)

def quads(G, colour_fn=None, base=None, two_sided=True):
    a, b, c, e = G[:-1, :-1], G[1:, :-1], G[1:, 1:], G[:-1, 1:]
    polys = np.stack([a, b, c, e], axis=2).reshape(-1, 4, 3)
    n = np.cross(c - a, e - b).reshape(-1, 3); n /= (np.linalg.norm(n, axis=1, keepdims=True) + 1e-12)
    shade = 0.55 + 0.45 * np.abs(n @ LIGHT)
    cols = colour_fn(polys.mean(1)) if colour_fn is not None else np.tile(base, (len(polys), 1))
    return polys, np.clip(cols * shade[:, None], 0, 1)


# ----------------------------------------------------------------------------------------
# 5. export: meshes (STL), a 2-D profile (DXF) and a solid (STEP) that open in 3D Slicer / SolidWorks
# ----------------------------------------------------------------------------------------

def cad_pieces(poly, tol=0.10, corner_deg=35.0, n_fit=28):
    """Split a closed (rho, z) polygon at its sharp corners and describe each piece as either a straight line or a
    smooth cubic spline (approximating, tolerance tol mm, so measurement noise does not become wiggles).
    Returns a list of ('line', p0, p1) or ('spline', fit_points) covering the closed polygon in order."""
    P = np.asarray(poly, float); n = len(P)
    def ang(v): return np.arctan2(v[1], v[0])
    turn = np.zeros(n)
    for i in range(n):
        d1 = P[i] - P[(i - 2) % n]; d2 = P[(i + 2) % n] - P[i]
        turn[i] = abs((ang(d2) - ang(d1) + np.pi) % (2 * np.pi) - np.pi)
    cand = [i for i in range(n) if np.degrees(turn[i]) > corner_deg]
    corners = []
    for i in sorted(cand, key=lambda k: -turn[k]):                       # strongest first, at least 3 mm apart in space
        if all(np.linalg.norm(P[i] - P[j]) >= 3.0 for j in corners): corners.append(i)
    corners = sorted(corners)
    if not corners: corners = [0]
    pieces = []
    for a_, b_ in zip(corners, corners[1:] + [corners[0] + n]):
        idx = [k % n for k in range(a_, b_ + 1)]
        Q = P[idx]
        chord = Q[-1] - Q[0]; L = np.linalg.norm(chord)
        dev = np.abs(np.cross(chord, Q - Q[0])) / max(L, 1e-9) if L > 1e-9 else np.zeros(len(Q))
        if len(Q) < 6 or dev.max() < 0.05:
            pieces.append(("line", Q[0], Q[-1])); continue
        w = np.ones(len(Q)); w[0] = w[-1] = 1e3
        tck, _ = splprep([Q[:, 0], Q[:, 1]], w=w, s=len(Q) * tol ** 2, k=3)
        u = np.linspace(0, 1, 600); x_, y_ = splev(u, tck)
        seg = np.concatenate([[0], np.cumsum(np.hypot(np.diff(x_), np.diff(y_)))])
        q = np.linspace(0, seg[-1], n_fit)
        F = np.column_stack([np.interp(q, seg, x_), np.interp(q, seg, y_)])
        F[0], F[-1] = Q[0], Q[-1]
        pieces.append(("spline", F))
    return pieces

def revolve_mesh(poly, nseg=144):
    """Triangle mesh of the solid made by revolving a closed (rho, z) polygon 360 degrees about the z axis.
    Points with rho = 0 become a single vertex. Returns vertices (n,3) and faces (m,3), outward oriented."""
    poly = np.asarray(poly, float); n = len(poly)
    th = np.linspace(0, 2 * np.pi, nseg, endpoint=False)
    verts, vid = [], -np.ones((n, nseg), int)
    for j, (r, z) in enumerate(poly):
        if r < 1e-9:
            verts.append((0.0, 0.0, z)); vid[j, :] = len(verts) - 1
        else:
            for k, t in enumerate(th):
                verts.append((r * np.cos(t), r * np.sin(t), z)); vid[j, k] = len(verts) - 1
    faces = []
    for j in range(n):
        j2 = (j + 1) % n
        for k in range(nseg):
            k2 = (k + 1) % nseg
            p, q, r_, s_ = vid[j, k], vid[j2, k], vid[j2, k2], vid[j, k2]
            for tri in ((p, q, r_), (p, r_, s_)):
                if len(set(tri)) == 3: faces.append(tri)
    V = np.array(verts); F = np.array(faces)
    vol = signed_volume(V, F)
    if vol < 0: F = F[:, ::-1]
    return V, F

def signed_volume(V, F):
    t = V[F]
    return float(np.einsum("ij,ij->i", t[:, 0], np.cross(t[:, 1], t[:, 2])).sum() / 6.0)

def is_watertight(F):
    e = np.sort(np.concatenate([F[:, [0, 1]], F[:, [1, 2]], F[:, [2, 0]]]), axis=1)
    _, cnt = np.unique(e, axis=0, return_counts=True)
    return bool((cnt == 2).all())

def write_stl(path, V, F, name="sample2"):
    t = V[F].astype(np.float32)
    n = np.cross(t[:, 1] - t[:, 0], t[:, 2] - t[:, 0]); n /= (np.linalg.norm(n, axis=1, keepdims=True) + 1e-20)
    rec = np.zeros(len(F), dtype=[("n", "<f4", 3), ("v", "<f4", (3, 3)), ("a", "<u2")])
    rec["n"], rec["v"] = n.astype(np.float32), t
    with open(path, "wb") as fh:
        fh.write(name.encode()[:79].ljust(80, b" ")); fh.write(np.uint32(len(F)).tobytes()); fh.write(rec.tobytes())

# ----------------------------------------------------------------------------------------
def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--npz", default="sample2.npz"); ap.add_argument("--out", default="out")
    ap.add_argument("--dirs", type=int, default=500, help="directions in the global axis search")
    args = ap.parse_args(); os.makedirs(args.out, exist_ok=True)
    R = {}                                                     # results collected for results.json

    # ---------------- stage 1: grid, split, meshes ----------------
    HU, sp, origin, hu = load_grid(args.npz)
    R["voxels"] = int(len(hu)); R["voxel_size_mm"] = [round(float(s), 4) for s in sp]
    T, p1, p2 = histogram_valley(hu)
    R["hu_valley"] = T; R["hu_peaks"] = [p1, p2]
    T = 1212.0 if abs(T - 1212) <= 25 else T                  # 1212 = value used in the analysis (valley 1200-1225)
    present = HU >= 500
    body = ndi.binary_fill_holes(largest_component(HU >= T))
    crust = present & ~ndi.binary_dilation(body, iterations=2)
    Vb, Fb, Nb, Wb = mesh_from_field(body.astype(float), 0.5, sp, origin, sigma=1.0)      # body, 1212 level
    Fo = np.where(present, HU.astype(float), 0.0)                                         # unknown (<500) voxels set to 0
    Vo, Fo_, No, Wo = mesh_from_field(Fo, 500.0, sp, origin, sigma=0.6)                   # outer surface, 500 HU level
    R["body_voxels"] = int(body.sum()); R["body_volume_cm3"] = float(body.sum() * sp.prod() / 1000)
    R["crust_class_volume_cm3"] = float((present & ~body).sum() * sp.prod() / 1000)

    # side of the plate that carries the crust = vessel interior ("up")
    Xb = np.argwhere(body) * sp + origin; cb = Xb.mean(0)
    Xc = np.argwhere(crust) * sp + origin
    n3 = np.linalg.svd(Xb - cb, full_matrices=False)[2][2]
    n_up = n3 * np.sign((Xc.mean(0) - cb) @ n3)
    under_b = (Nb @ n_up) < -0.35                                                         # underside of the body mesh
    print("[1] voxels %d, valley %.0f HU, body %.1f cm3" % (len(hu), R["hu_valley"], R["body_volume_cm3"]))

    # ---------------- stage 2: axis ----------------
    rng = np.random.default_rng(21)
    Pu, Nu, Wu = Vb[under_b], Nb[under_b], Wb[under_b]
    sel = rng.choice(len(Pu), size=2200, replace=False, p=Wu / Wu.sum())
    sols, cen = global_axis_search(Pu[sel], Nu[sel], n_dirs=args.dirs, keep=12)
    p_best = sols[0][1]
    # re-optimise on a larger sample
    big = Pu[rng.choice(len(Pu), size=min(6000, len(Pu)), replace=False, p=Wu / Wu.sum())]
    best = None
    for f, p in sols[:5]:
        r = minimize(objective, p, args=(big, cen), method="Nelder-Mead", options=dict(xatol=5e-4, fatol=1e-5, maxiter=500))
        if best is None or r.fun < best[0]: best = (r.fun, r.x)
    f0, p0 = best
    a0, c0 = axis_from(p0, cen)
    if a0 @ n_up < 0: a0 = -a0
    R["axis_profile_only"] = dict(direction=np.round(a0, 4).tolist(), profile_rms_mm=round(float(f0), 3))
    R["alternative_axes"] = [dict(rms_mm=round(float(f), 3), direction=np.round((axis_from(p, cen)[0]), 3).tolist())
                             for f, p in sols[:6]]

    # independent check: the finished lip is flat and perpendicular to the axis
    e1, e2 = centred_basis(a0, c0, Vb[under_b])
    q = Vb - c0; hh = q @ a0; xx = q @ e1; yy = q @ e2
    hh = hh - hh[under_b].min(); rho = np.hypot(xx, yy); phi = np.degrees(np.arctan2(yy, xx))
    top_zone = (rho > 95) & (hh > 50)
    phi_mid = np.median(phi[top_zone]) if top_zone.any() else 0.0
    hcrest = np.percentile(hh[(rho > 95) & (np.abs(phi - phi_mid) < 30)], 99)
    Nl_list = []
    for (dlo, dhi, plo, phi_hi, nc) in [(2.5, -1.0, -28, 30, 0.8), (3.5, -1.0, -28, 30, 0.7), (1.5, -1.0, -25, 28, 0.85),
                                        (2.5, -1.0, -15, 25, 0.8), (3.0, -1.0, -28, 30, 0.6)]:
        m = ((hh > hcrest - dlo) & (hh < hcrest - dhi) & (phi > phi_mid + plo) & (phi < phi_mid + phi_hi)
             & (Nb @ a0 > nc) & (rho > 95))
        mu, nl, rms = fit_plane(Vb[m], Wb[m])
        if nl @ a0 < 0: nl = -nl
        Nl_list.append((nl, rms, int(m.sum())))
    Nl = np.mean([n for n, _, _ in Nl_list], axis=0); Nl /= np.linalg.norm(Nl)
    lip_spread = max(np.degrees(np.arccos(np.clip(n @ Nl, -1, 1))) for n, _, _ in Nl_list)
    lip_vs_profile = float(np.degrees(np.arccos(np.clip(Nl @ a0, -1, 1))))
    R["lip"] = dict(plane_rms_mm=[round(r, 2) for _, r, _ in Nl_list], normal=np.round(Nl, 4).tolist(),
                    angle_to_profile_axis_deg=round(lip_vs_profile, 2), selection_spread_deg=round(float(lip_spread), 2))
    # final axis: direction from the lip, position from the profile fit
    th, ph = np.arccos(Nl[2]), np.arctan2(Nl[1], Nl[0])
    bestc = (1e9, None)
    for u0 in (-120, -60, 0, 60, 120):
        for v0 in (-120, -60, 0, 60, 120):
            r = minimize(lambda uv: objective(np.array([th, ph, uv[0], uv[1]]), big, cen), [u0, v0],
                         method="Nelder-Mead", options=dict(xatol=1e-2, fatol=1e-5, maxiter=250))
            if r.fun < bestc[0]: bestc = (r.fun, r.x)
    a, c = axis_from(np.array([th, ph, *bestc[1]]), cen)
    if a @ n_up < 0: a = -a
    e1, e2 = centred_basis(a, c, Vb[under_b])
    ref = Vb.mean(0)
    axis_pt = c + ((ref - c) @ a) * a
    R["axis"] = dict(direction=np.round(a, 4).tolist(), point_mm=np.round(axis_pt, 2).tolist(),
                     tilt_from_scanner_z_deg=round(float(np.degrees(np.arccos(abs(a[2])))), 1),
                     profile_rms_mm=round(float(bestc[0]), 3))
    print("[2] axis", np.round(a, 4), "profile scatter %.3f mm; lip vs profile axis %.2f deg" % (bestc[0], lip_vs_profile))

    # ---------------- stage 3: profiles ----------------
    def frame(P):
        q = P - c; h = q @ a; x = q @ e1; y = q @ e2
        return np.hypot(x, y), h, np.arctan2(y, x)
    rb, hb, pb = frame(Vb); ro, ho, po = frame(Vo)
    # underside of the 500 HU mesh: faces away from the crust and lies close to the body
    dist_to_body = cKDTree(Vb).query(Vo)[0]
    near = dist_to_body < 2.5                                   # excludes crust lumps that do not touch the plate
    under_o = largest_vertex_component(((No @ n_up) < -0.35) & near, Fo_, len(Vo))
    # half-maximum edge correction. Slope across the edge from the distance between the 1212 and 500 HU surfaces.
    d12_500 = float(np.median(cKDTree(Vb[under_b]).query(Vo[under_o])[0]))
    grad = (T - 500.0) / d12_500                                # HU per mm
    HALF = (-1000.0 + 1350.0) / 2.0                             # half-way between air (-1000) and the body plateau (~1350)
    delta500 = (500.0 - HALF) / grad                            # 500 HU surface -> half-max surface (outwards)
    delta_T = (T - HALF) / grad                                 # 1212 HU surface -> half-max surface
    R["edge"] = dict(distance_1212_to_500_mm=round(d12_500, 2), gradient_HU_per_mm=int(round(grad)),
                     half_max_HU=HALF, outward_shift_from_500_mm=round(delta500, 2), shift_from_1212_mm=round(delta_T, 2))
    hmin = ho[under_o].min()                                    # datum of this stage: lowest underside point (500 HU level)
    wrap = lambda p: (p - np.angle(np.mean(np.exp(1j * po[under_o]))) + np.pi) % (2 * np.pi) - np.pi

    # ---- outer profile: medians in bins of polar angle in the (rho, h) plane, then a smoothing spline ----
    ru_, hu_ = ro[under_o], ho[under_o] - hmin
    cpr, cph = 0.5 * (np.percentile(ru_, 1) + np.percentile(ru_, 99)), 0.8 * hu_.max()     # a point inside the concavity
    ang = np.arctan2(hu_ - cph, ru_ - cpr)
    uu = np.linspace(0, 1, 900)
    def fit_spline(param):
        """medians of the underside points in 70 equal-count bins of 'param', then a smoothing spline (0.25 mm)"""
        edges = np.quantile(param, np.linspace(0, 1, 71)); pts2 = []
        for k in range(70):
            m = (param >= edges[k]) & (param <= edges[k + 1])
            if m.sum() >= 20: pts2.append((np.median(ru_[m]), np.median(hu_[m])))
        pts2 = np.array(pts2)
        tck, _ = splprep([pts2[:, 0], pts2[:, 1]], s=len(pts2) * 0.25 ** 2, k=3)
        pr, ph = splev(uu, tck)
        return tck, pr, ph
    tck, prho, ph_ = fit_spline(np.arctan2(hu_ - cph, ru_ - cpr))          # first pass: order by polar angle
    for _ in range(2):                                                    # then order by arc length along the fitted curve
        s_ = np.concatenate([[0], np.cumsum(np.hypot(np.diff(prho), np.diff(ph_)))])
        keep = cKDTree(np.column_stack([prho, ph_])).query(np.column_stack([ru_, hu_]))[1]
        tck, prho, ph_ = fit_spline(s_[keep])
    # the spline ends are medians over the first/last bin: extend the curve to the true tip of the sherd and past the lip plane
    P2 = np.column_stack([prho, ph_]); n2 = len(P2)
    tau0 = P2[0] - P2[3]; tau0 /= np.linalg.norm(tau0)
    tau1 = P2[-1] - P2[-4]; tau1 /= np.linalg.norm(tau1)
    pts_u = np.column_stack([ru_, hu_]); nn = cKDTree(P2).query(pts_u)[1]
    ext0 = float(np.percentile((pts_u[nn == 0] - P2[0]) @ tau0, 99.5)) if (nn == 0).sum() > 5 else 0.0
    ext1 = (float(np.percentile((pts_u[nn == n2 - 1] - P2[-1]) @ tau1, 99.5)) if (nn == n2 - 1).sum() > 5 else 0.0)
    ext0, ext1 = max(ext0, 0.0), max(ext1, 0.0) + 6.0
    st = P2[0] + np.arange(ext0, 0.0, -0.2)[:, None] * tau0
    en = P2[-1] + np.arange(0.2, ext1 + 1e-9, 0.2)[:, None] * tau1
    P2 = np.vstack([st, P2, en]); prho, ph_ = P2[:, 0], P2[:, 1]
    L = len(prho)
    dpr, dph = np.gradient(prho), np.gradient(ph_)
    s_arc = np.concatenate([[0], np.cumsum(np.hypot(np.diff(prho), np.diff(ph_)))])
    tx, ty = dpr / np.hypot(dpr, dph), dph / np.hypot(dpr, dph)
    nur, nuh = -ty, tx                                          # left normal
    mid = L // 2
    if nur[mid] * (cpr - prho[mid]) + nuh[mid] * (cph - ph_[mid]) < 0:   # make it point into the vessel interior
        nur, nuh = -nur, -nuh
    R["profile_extension_mm"] = dict(tip=round(ext0, 1), top=round(ext1 - 6.0, 1))
    # ---- correct the outer profile against the CT field itself ----
    # The mesh classification can be biased (e.g. by crust touching the wall).  Along the outward normal of the current
    # profile, at many azimuths, find where the field falls through 500 HU; the median over azimuth moves the profile.
    Ffield = np.where(HU >= 500, HU.astype(float), 0.0)
    us = np.arange(-3.0, 6.0001, 0.1)
    az_grid = np.radians(np.arange(-36, 36.1, 3.0))
    samp = np.arange(0, L, 8)
    E = np.full((len(samp), len(az_grid)), np.nan)
    for j, phj in enumerate(az_grid):
        rr_ = prho[samp][:, None] - us[None, :] * nur[samp][:, None]
        hh_ = ph_[samp][:, None] - us[None, :] * nuh[samp][:, None] + hmin
        Pw = (c + (rr_ * np.cos(phj))[..., None] * e1 + (rr_ * np.sin(phj))[..., None] * e2 + hh_[..., None] * a)
        Fv = ndi.map_coordinates(Ffield, np.moveaxis((Pw - origin) / sp, -1, 0), order=1, mode="constant", cval=0.0)
        inside0 = Fv[:, 0] >= 900.0
        below = Fv < 500.0
        k1 = np.where(below.any(1), below.argmax(1), -1)
        ok = inside0 & (k1 > 0)
        rows = np.where(ok)[0]
        f0 = Fv[rows, k1[rows] - 1]; f1 = Fv[rows, k1[rows]]
        E[rows, j] = us[k1[rows] - 1] + (f0 - 500.0) / (f0 - f1) * 0.1
    e_med = np.nanmedian(np.where(np.isfinite(E), E, np.nan), axis=1)
    okr = np.isfinite(e_med)
    e_full = np.interp(np.arange(L), samp[okr], e_med[okr])
    e_full = ndi.gaussian_filter1d(medfilt(e_full, 9), 12, mode="nearest")           # ~1 mm smoothing along the profile
    R["outer_profile_correction_mm"] = dict(median=round(float(np.median(e_full)), 2), min=round(float(e_full.min()), 2),
                                            max=round(float(e_full.max()), 2),
                                            top_third_median=round(float(np.median(e_full[int(0.7 * L):])), 2))
    prho = prho - nur * e_full; ph_ = ph_ - nuh * e_full
    s_arc = np.concatenate([[0], np.cumsum(np.hypot(np.diff(prho), np.diff(ph_)))])
    tree2 = cKDTree(np.column_stack([prho, ph_]))
    res_d = tree2.query(np.column_stack([ru_, hu_]))[0]
    R["outer_profile_scatter_mm"] = dict(median=round(float(np.median(res_d)), 2), p90=round(float(np.percentile(res_d, 90)), 2),
                                        p95=round(float(np.percentile(res_d, 95)), 2))

    # ---- is the vessel round? the edge radius at each azimuth, relative to the azimuth-median profile ----
    s_samp = s_arc[samp]
    m_base = (s_samp > 5) & (s_samp < 0.45 * s_arc[-1]); m_wall = s_samp > 0.7 * s_arc[-1]
    dE = E - np.nanmedian(E, axis=1, keepdims=True)
    def az_range(mask):
        med = np.array([np.nanmedian(dE[mask, j]) if np.isfinite(dE[mask, j]).sum() > 8 else np.nan for j in range(len(az_grid))])
        okj = np.isfinite(med)
        return (round(float(np.ptp(med[okj])), 1) if okj.sum() > 2 else None), np.degrees(az_grid[okj]), med[okj]
    rb_, _, _ = az_range(m_base); rw_, az_w, med_w = az_range(m_wall)
    R["roundness_mm"] = dict(base_departure_range=rb_, wall_departure_range=rw_,
                             note="peak-to-peak over azimuth of the median edge position (500 HU) relative to the azimuth-median profile")

    # ---- inner face (1212 HU body surface facing the interior) and wall thickness ----
    # inner face: body-mesh vertices whose normal points into the vessel interior, seen in the (rho, h) meridian plane
    rho_hat = np.column_stack([np.cos(pb), np.sin(pb)])                  # radial unit vector (x', y') of each body vertex
    n_rho = (Nb @ e1) * rho_hat[:, 0] + (Nb @ e2) * rho_hat[:, 1]; n_h = Nb @ a
    near_i = cKDTree(np.column_stack([prho, ph_])).query(np.column_stack([rb, hb - hmin]))[1]
    dot_in = n_rho * nur[near_i] + n_h * nuh[near_i]
    off_in = (rb - prho[near_i]) * nur[near_i] + (hb - hmin - ph_[near_i]) * nuh[near_i]
    inner_b = (dot_in > 0.6) & (off_in > 1.5) & (off_in < 18.0)
    inner2 = np.column_stack([rb[inner_b], hb[inner_b] - hmin])
    # crust lumps only: opening removes the thin partial-volume shell (about 1-2 voxels) that surrounds every free face
    crust_lumps = ndi.binary_opening(present & ~body, structure=ndi.generate_binary_structure(3, 1), iterations=2)
    covered = cKDTree(np.argwhere(crust_lumps) * sp + origin).query(Vb[inner_b])[0] < 2.0     # crust lies against this face
    R["inner_face_crust_covered_fraction"] = round(float(covered.mean()), 2)
    ti = cKDTree(inner2)
    t_meas = np.full(L, np.nan); cov_frac = np.zeros(L)
    for i in range(L):
        idx = np.array(ti.query_ball_point([prho[i], ph_[i]], r=12.0))
        if len(idx) < 8: continue
        dv = inner2[idx] - np.array([prho[i], ph_[i]])
        along = dv @ np.array([tx[i], ty[i]]); nrm_d = dv @ np.array([nur[i], nuh[i]])
        m = (np.abs(along) < 1.5) & (nrm_d > 1.0)
        if m.sum() < 5: continue
        t_meas[i] = np.median(nrm_d[m]); cov_frac[i] = covered[idx[m]].mean()
    good = np.isfinite(t_meas)
    t_est = t_meas + delta500 + np.where(cov_frac > 0.5, 0.0, delta_T)     # half-max thickness
    t_f = np.interp(s_arc, s_arc[good], t_est[good])
    k = int(10.0 / (s_arc[-1] / len(s_arc))) | 1
    t_s = medfilt(t_f, kernel_size=max(k, 5))

    # ---- lip crest on the 500 HU mesh (flat top, perpendicular to the axis) ----
    hrel_o = ho - hmin
    in_az = np.abs(wrap(po)) < np.radians(30)
    wall_o = (ro > 95) & in_az & near
    crest500 = float(np.percentile(hrel_o[wall_o], 99))
    lipo = wall_o & ((No @ a) > 0.85) & (hrel_o > crest500 - 2.5)
    for _ in range(3):                                          # tighten around the flat top
        med = np.median(hrel_o[lipo]); lipo = wall_o & ((No @ a) > 0.85) & (np.abs(hrel_o - med) < 1.2)
    lip_h500 = float(np.median(hrel_o[lipo])); lip_r = ro[lipo]
    R["lip_crest_500HU_mesh"] = dict(n_vertices=int(lipo.sum()), height_sd_mm=round(float(np.std(hrel_o[lipo])), 2),
                                     radial_extent_mm=[round(float(np.percentile(lip_r, 5)), 1), round(float(np.percentile(lip_r, 95)), 1)])

    # ---- half-max outer / inner curves, cut at the lip plane ----
    o_rho_f = prho - nur * delta500; o_h_f = ph_ - nuh * delta500              # outer, half-max (uncut)
    hshift = o_h_f.min()
    H_top = lip_h500 + delta500                                                # lip plane, half-max, same datum as o_h_f
    i_rho_f = o_rho_f + nur * t_s; i_h_f = o_h_f + nuh * t_s
    # a parallel offset folds over itself where the wall is thicker than the local radius of curvature (the shoulder):
    # drop the points of that loop, which leaves the sharp inner corner a real wall would have
    d_to_outer = cKDTree(np.column_stack([o_rho_f, o_h_f])).query(np.column_stack([i_rho_f, i_h_f]))[0]
    keep_i = d_to_outer > (t_s - 0.6)
    i_rho_f, i_h_f = i_rho_f[keep_i], i_h_f[keep_i]
    seg_i = np.concatenate([[0], np.cumsum(np.hypot(np.diff(i_rho_f), np.diff(i_h_f)))])
    q_i = np.arange(0, seg_i[-1], 0.25)                                        # resample every 0.25 mm and smooth lightly (~1 mm)
    i_rho_f = ndi.gaussian_filter1d(np.interp(q_i, seg_i, i_rho_f), 5, mode="nearest")
    i_h_f = ndi.gaussian_filter1d(np.interp(q_i, seg_i, i_h_f), 5, mode="nearest")
    def cut(rr, hh_, H):
        j = np.where(hh_ >= H)[0]
        if len(j) == 0: return rr, hh_
        j0 = j[0]; f = (H - hh_[j0 - 1]) / (hh_[j0] - hh_[j0 - 1])
        return np.append(rr[:j0], rr[j0 - 1] + f * (rr[j0] - rr[j0 - 1])), np.append(hh_[:j0], H)
    o_rho, o_h = cut(o_rho_f, o_h_f, H_top); i_rho, i_h = cut(i_rho_f, i_h_f, H_top)
    o_h_d, i_h_d = o_h - hshift, i_h - hshift                                  # heights above the lowest surviving point
    rim_out, rim_in = float(o_rho[-1]), float(i_rho[-1])
    i_max = int(np.argmax(o_rho))
    R["vessel"] = dict(
        lowest_surviving_point_radius_mm=round(float(o_rho[0]), 1),
        rim_height_above_lowest_point_mm=round(float(o_h_d[-1]), 1),
        rim_outer_radius_mm=round(rim_out, 1), rim_outer_diameter_mm=round(2 * rim_out),
        rim_inner_radius_mm=round(rim_in, 1), lip_radial_width_mm=round(rim_out - rim_in, 1),
        max_outer_radius_mm=round(float(o_rho.max()), 1), max_outer_diameter_mm=round(2 * float(o_rho.max())),
        height_of_max_diameter_mm=round(float(o_h_d[i_max]), 1),
        thickness_mm=dict(median=round(float(np.median(t_s)), 1), min=round(float(t_s.min()), 1), max=round(float(t_s.max()), 1)),
        surviving_azimuth_deg=round(float(np.degrees(np.percentile(wrap(po[under_o]), 99.5) - np.percentile(wrap(po[under_o]), 0.5)))),
        lip_level_azimuth_deg=round(float(np.degrees(np.ptp(wrap(po[lipo]))))))
    print("[3] rim diameter %.0f mm, max diameter %.0f mm, rim height %.1f mm, thickness median %.1f mm, lip width %.1f mm" %
          (2 * rim_out, 2 * o_rho.max(), o_h_d[-1], np.median(t_s), rim_out - rim_in))

    # ---- footprint of the surviving sherd along the profile (azimuth range per arc-length bin) ----
    tree_d = cKDTree(np.column_stack([o_rho_f, o_h_f - hshift]))
    pts_under_d = np.column_stack([ro[under_o], ho[under_o] - hmin - hshift])
    s_pt = s_arc[tree_d.query(pts_under_d)[1]]
    nb = 60; sb = np.linspace(0, s_arc[-1], nb + 1)
    ph_lo = np.full(nb, np.nan); ph_hi = np.full(nb, np.nan); pw = wrap(po[under_o])
    for k in range(nb):
        m = (s_pt >= sb[k]) & (s_pt <= sb[k + 1])
        if m.sum() > 10: ph_lo[k], ph_hi[k] = np.percentile(pw[m], 1), np.percentile(pw[m], 99)
    okb = np.isfinite(ph_lo)
    def surviving(pt):
        """True where a point of the swept surface lies in the part that was actually scanned."""
        r_ = np.hypot(pt[:, 0], pt[:, 1]); z_ = pt[:, 2]; p_ = np.arctan2(pt[:, 1], pt[:, 0])
        s_here = s_arc[tree_d.query(np.column_stack([r_, z_]))[1]]
        kk = np.clip(np.searchsorted(sb, s_here) - 1, 0, nb - 1)
        return okb[kk] & (p_ >= ph_lo[kk]) & (p_ <= ph_hi[kk])

    # table for reuse
    tt = np.linspace(0, 1, 200)
    def rs(v): return np.interp(tt, np.linspace(0, 1, len(v)), v)
    np.savetxt(os.path.join(args.out, "profile_sample2.csv"),
               np.column_stack([tt, rs(o_rho), rs(o_h_d), rs(i_rho), rs(i_h_d)]), delimiter=",",
               header="t_along_profile,outer_radius_mm,outer_height_mm,inner_radius_mm,inner_height_mm", comments="")

    # ---- validation figure: real CT sections through the axis with the extracted profiles on top ----
    Ffield = np.where(HU >= 500, HU.astype(float), 0.0)
    rr = np.linspace(65, 130, 260); yy = np.linspace(-8, 68, 300); RRg, YYg = np.meshgrid(rr, yy)
    def section(phi_c):
        acc = []
        for dphi in np.radians([-1.5, -0.75, 0.0, 0.75, 1.5]):
            ph = phi_c + dphi
            Pg = (c + (RRg * np.cos(ph))[..., None] * e1 + (RRg * np.sin(ph))[..., None] * e2
                  + (YYg + hmin + hshift)[..., None] * a)
            acc.append(ndi.map_coordinates(Ffield, np.moveaxis((Pg - origin) / sp, -1, 0), order=1, mode="constant", cval=0.0))
        return np.mean(acc, axis=0)
    fig, axs = plt.subplots(1, 3, figsize=(17, 6.4))
    phi_mid_ = np.angle(np.mean(np.exp(1j * po[under_o])))
    for axm, dphi_deg in zip(axs, (-20, 0, 20)):
        img = section(phi_mid_ + np.radians(dphi_deg))
        axm.imshow(img, origin="lower", extent=[rr[0], rr[-1], yy[0], yy[-1]], cmap="gray_r", vmin=0, vmax=1600, aspect="equal")
        axm.contour(RRg, YYg, img, levels=[500], colors="cyan", linewidths=0.7)
        axm.contour(RRg, YYg, img, levels=[T], colors="red", linewidths=0.6)
        axm.plot(o_rho, o_h_d, color="yellow", lw=1.4, label="outer profile (half-max)")
        axm.plot(i_rho, i_h_d, color="orange", lw=1.4, label="inner profile (half-max)")
        axm.set_title("CT section through the axis, azimuth %+d deg\n(cyan 500 HU, red %d HU)" % (dphi_deg, int(T)), fontsize=9)
        axm.set_xlabel("distance from the axis (mm)"); axm.set_ylabel("height above the lowest surviving point (mm)")
    axs[0].legend(fontsize=7, loc="upper left")
    fig.tight_layout(); fig.savefig(os.path.join(args.out, "fig_meridian.png"), dpi=100); plt.close(fig)

    # ---- lip check figure: the finished lip is flat and level around the axis ----
    fig, axl = plt.subplots(1, 3, figsize=(17, 6))
    hrel_all = ho - hmin
    axl[0].scatter(ro[near], hrel_all[near], s=0.5, color="#9aa8b8", alpha=0.3, label="all near-body surface")
    axl[0].scatter(ro[under_o], hrel_all[under_o], s=0.8, color="#3a5f8a", alpha=0.5, label="underside")
    axl[0].scatter(ro[lipo], hrel_all[lipo], s=4, color="#d03030", label="finished lip (flat top)")
    axl[0].axhline(lip_h500, color="#d03030", lw=0.8, ls="--"); axl[0].set_aspect("equal"); axl[0].set_xlabel("distance from the axis (mm)")
    axl[0].set_ylabel("height above the lowest underside point (mm)"); axl[0].set_title("Profile with the lip highlighted", fontsize=10); axl[0].legend(fontsize=8, markerscale=6, loc="upper left")
    azd = np.degrees(wrap(po)); wall_z = near & (ro > 95) & (hrel_all > 40)
    bins = np.linspace(azd[wall_z].min(), azd[wall_z].max(), 45); ix = np.digitize(azd[wall_z], bins)
    crest_az = [hrel_all[wall_z][ix == kk].max() if (ix == kk).sum() > 5 else np.nan for kk in range(1, len(bins))]
    axl[1].plot(0.5 * (bins[1:] + bins[:-1]), crest_az, "o-"); axl[1].axhspan(lip_h500 - 1, lip_h500 + 1, color="#d03030", alpha=0.15, label="+/- 1 mm")
    axl[1].set_ylim(30, 70); axl[1].set_xlabel("azimuth (deg)"); axl[1].set_ylabel("highest point of the wall (mm)"); axl[1].legend(fontsize=8)
    axl[1].set_title("Crest height around the axis (level = flat lip)", fontsize=10)
    axl[2].scatter(ro[under_o] * np.cos(wrap(po)[under_o]), ro[under_o] * np.sin(wrap(po)[under_o]), s=0.6, color="#3a5f8a", alpha=0.5)
    axl[2].scatter(ro[lipo] * np.cos(wrap(po)[lipo]), ro[lipo] * np.sin(wrap(po)[lipo]), s=3, color="#d03030"); axl[2].plot(0, 0, "k+", ms=14, mew=2)
    axl[2].set_aspect("equal"); axl[2].set_title("Plan view down the axis", fontsize=10); axl[2].set_xlabel("mm")
    fig.suptitle("Sample 2: the finished lip confirms the axis", fontsize=12, fontweight="bold"); fig.tight_layout(rect=(0, 0, 1, 0.95))
    fig.savefig(os.path.join(args.out, "fig_lip_check.png"), dpi=100); plt.close(fig)

    # ---------------- stage 4: sweep ----------------
    # ILLUSTRATIVE closure of the unmeasured centre: sphere cap tangent to the outer profile at its lowest point
    n0 = 8
    slope = float(np.arctan2(o_h_d[n0] - o_h_d[0], o_rho[n0] - o_rho[0]))
    Rs = o_rho[0] / np.sin(slope)
    thc = np.linspace(0, slope, 40)
    clos_o_rho = Rs * np.sin(thc); clos_o_h = Rs * (np.cos(slope) - np.cos(thc))
    t0 = float(t_s[0])
    clos_i_rho = (Rs - t0) * np.sin(thc)
    clos_i_h = (Rs - t0) * (np.cos(slope) - np.cos(thc)) + t0 * np.cos(slope)
    R["closure_illustrative"] = dict(sphere_radius_mm=round(float(Rs)), drop_below_lowest_surviving_point_mm=round(float(-clos_o_h[0]), 1),
                                     tangent_slope_deg=round(float(np.degrees(slope)), 1),
                                     note="NOT measured. Smooth cap tangent to the surviving base, for display only.")
    def thin(rr, zz, n=140):
        s_ = np.concatenate([[0], np.cumsum(np.hypot(np.diff(rr), np.diff(zz)))])
        q = np.linspace(0, s_[-1], n)
        return np.interp(q, s_, rr), np.interp(q, s_, zz)
    so_rho, so_h = thin(o_rho, o_h_d); si_rho, si_h = thin(i_rho, i_h_d)
    C_OUT, C_IN = _rgb("#c9a07a"), _rgb("#b58c66"); C_INF_O, C_INF_I = _rgb("#9db6cc"), _rgb("#b9cde0")
    C_GHOST_O, C_GHOST_I = _rgb("#d0d0d0"), _rgb("#e2e2e2"); C_EDGE = _rgb("#8f8f8f")
    def build(phis, cut_edges=False):
        P, C = [], []
        for rr, zz, c_s, c_i in ((so_rho, so_h, C_OUT, C_INF_O), (si_rho, si_h, C_IN, C_INF_I)):
            pl, cl = quads(revolve(rr, zz, phis), colour_fn=lambda cn: np.where(surviving(cn)[:, None], c_s[None, :], c_i[None, :]))
            P.append(pl); C.append(cl)
        G = np.stack([revolve(so_rho[-1:], so_h[-1:], phis)[0], revolve(si_rho[-1:], si_h[-1:], phis)[0]], axis=0)
        pl, cl = quads(G, base=C_EDGE); P.append(pl); C.append(cl)
        for rr, zz, col in ((clos_o_rho, clos_o_h, C_GHOST_O), (clos_i_rho, clos_i_h, C_GHOST_I)):
            pl, cl = quads(revolve(rr, zz, phis), base=col); P.append(pl); C.append(cl)
        if cut_edges:
            for col_i in (0, -1):
                Go = revolve(so_rho, so_h, phis)[:, col_i]; Gi = revolve(si_rho, si_h, phis)[:, col_i]
                pl, cl = quads(np.stack([Go, Gi], 0), base=C_EDGE); P.append(pl); C.append(cl)
        return np.concatenate(P), np.concatenate(C)

    N_PHI = 170
    phis_full = np.linspace(-np.pi, np.pi, N_PHI + 1)
    lim = float(o_rho.max()) * 1.02; Ht = float(o_h_d.max()); zlo = float(clos_o_h.min())
    def add3(ax, P, C, elev, azim, title):
        ax.add_collection3d(Poly3DCollection(P, facecolors=C, edgecolors="none", linewidths=0))
        ax.set_xlim(-lim, lim); ax.set_ylim(-lim, lim); ax.set_zlim(zlo - 2, Ht + 2)
        ax.set_box_aspect((2 * lim, 2 * lim, Ht - zlo + 4), zoom=1.2)
        ax.view_init(elev=elev, azim=azim); ax.set_title(title, fontsize=10)
        ax.tick_params(labelsize=6); ax.set_xlabel("x", fontsize=7); ax.set_ylabel("y", fontsize=7)
    fig = plt.figure(figsize=(16, 10.6))
    P_, C_ = build(phis_full)
    add3(fig.add_subplot(2, 2, 1, projection="3d"), P_, C_, 24, -30,
         "Whole vessel from outside\nclay = scanned sherd, blue = repeated by symmetry, grey = centre not measured")
    ph_cut = np.linspace(np.radians(-120), np.radians(150), int(N_PHI * 0.75) + 1)
    Pc_, Cc_ = build(ph_cut, cut_edges=True)
    add3(fig.add_subplot(2, 2, 2, projection="3d"), Pc_, Cc_, 42, 195, "Cutaway (quarter removed) to show the inside")
    add3(fig.add_subplot(2, 2, 3, projection="3d"), P_, C_, 88, -90, "Looking down into the vessel")
    ax4 = fig.add_subplot(2, 2, 4)
    for sg in (1, -1):
        ax4.fill(sg * np.concatenate([o_rho, i_rho[::-1]]), np.concatenate([o_h_d, i_h_d[::-1]]), color="#d9bfa2", alpha=0.9 if sg > 0 else 0.55)
        ax4.plot(sg * o_rho, o_h_d, color="#5a3d22", lw=1.8); ax4.plot(sg * i_rho, i_h_d, color="#8a6a44", lw=1.4)
        ax4.plot(sg * clos_o_rho, clos_o_h, color="#999999", lw=1.4, ls="--"); ax4.plot(sg * clos_i_rho, clos_i_h, color="#bbbbbb", lw=1.1, ls="--")
    ax4.axvline(0, color="gray", ls="-.", lw=0.8)
    ax4.annotate("", (-rim_out, o_h_d[-1] + 7), (rim_out, o_h_d[-1] + 7), arrowprops=dict(arrowstyle="<->", lw=0.9))
    ax4.text(0, o_h_d[-1] + 9, "rim outer diameter ~%.0f mm" % (2 * rim_out), ha="center", fontsize=8)
    hm_ = float(o_h_d[i_max])
    ax4.annotate("", (-o_rho.max(), hm_), (o_rho.max(), hm_), arrowprops=dict(arrowstyle="<->", lw=0.9, color="#3a5f8a"))
    ax4.text(0, hm_ + 2, "max diameter ~%.0f mm" % (2 * o_rho.max()), ha="center", fontsize=8, color="#3a5f8a")
    ax4.text(0, zlo - 9, "centre of the base not measured (dashed = illustrative closure)", ha="center", fontsize=8, color="#666666")
    ax4.set_aspect("equal"); ax4.set_xlabel("mm"); ax4.set_ylabel("height above the lowest surviving point (mm)")
    ax4.set_xlim(-lim - 12, lim + 12); ax4.set_ylim(zlo - 14, o_h_d[-1] + 18)
    ax4.set_title("Cross-section through the axis", fontsize=10); ax4.grid(alpha=0.25)
    fig.suptitle("Sample 2: rotationally symmetric vessel swept from the scanned sherd", fontsize=12, fontweight="bold")
    fig.subplots_adjust(left=0.03, right=0.98, top=0.92, bottom=0.05, wspace=0.05, hspace=0.2)
    fig.savefig(os.path.join(args.out, "fig_sweep.png"), dpi=110); plt.close(fig)

    # ---- profile / thickness figure ----
    fig, axs = plt.subplots(1, 3, figsize=(17, 6))
    axs[0].scatter(inner2[:, 0], inner2[:, 1] - hshift, s=1, color="#999999", alpha=0.25, label="inner face (1212 HU surface)")
    axs[0].scatter(ro[under_o], ho[under_o] - hmin - hshift, s=1, color="#3a5f8a", alpha=0.5, label="underside (500 HU surface)")
    axs[0].plot(o_rho, o_h_d, color="k", lw=1.3, label="outer profile (half-max)"); axs[0].plot(i_rho, i_h_d, color="#c04020", lw=1.3, label="inner profile (half-max)")
    axs[0].set_aspect("equal"); axs[0].set_xlabel("distance from the axis (mm)"); axs[0].set_ylabel("height above the lowest surviving point (mm)")
    axs[0].legend(fontsize=7, loc="upper left"); axs[0].set_title("Extracted profiles", fontsize=10)
    axs[1].plot(s_arc, t_meas, ".", ms=2, color="#bbbbbb", label="measured (500 HU outer to 1212 HU inner)")
    axs[1].plot(s_arc, t_s, color="#c04020", lw=2, label="half-max estimate (used)")
    axs[1].set_xlabel("distance along the outer profile from the lowest point (mm)"); axs[1].set_ylabel("wall thickness (mm)")
    axs[1].legend(fontsize=8); axs[1].set_title("Wall thickness in the meridian plane", fontsize=10); axs[1].grid(alpha=0.3)
    axs[2].hist(res_d, bins=60, color="#3a5f8a"); axs[2].set_xlabel("distance of underside points from the fitted profile (mm)")
    axs[2].set_title("How thin is the underside profile?", fontsize=10)
    fig.tight_layout(); fig.savefig(os.path.join(args.out, "fig_profile.png"), dpi=105); plt.close(fig)


    # ---------------- stage 5: export ----------------
    def thin_(rr, zz, n):
        s_ = np.concatenate([[0], np.cumsum(np.hypot(np.diff(rr), np.diff(zz)))]); q = np.linspace(0, s_[-1], n)
        return np.column_stack([np.interp(q, s_, rr), np.interp(q, s_, zz)])
    Po = thin_(o_rho, o_h_d, 90); Pi = thin_(i_rho, i_h_d, 90)                     # outer / inner curves, tip -> lip
    Co = np.column_stack([clos_o_rho, clos_o_h]); Ci = np.column_stack([clos_i_rho, clos_i_h])   # caps, pole -> tip
    poly_meas = np.vstack([Po, Pi[::-1]])                                          # closed: outer up, lip, inner down, tip
    poly_full = np.vstack([Co[:-1], Po, Pi[::-1], Ci[::-1][1:]])                   # with the ILLUSTRATIVE base cap
    exp = {}
    def to_scanner(V, zshift=0.0):
        z_raw = V[:, 2] - zshift + hshift + hmin
        return c + V[:, [0]] * e1 + V[:, [1]] * e2 + z_raw[:, None] * a
    for tag, poly in (("measured", poly_meas), ("with_illustrative_base", poly_full)):
        V, Fc = revolve_mesh(poly)
        z0 = -V[:, 2].min() if tag != "measured" else 0.0                          # local files rest on z = 0
        write_stl(os.path.join(args.out, "sample2_vessel_%s_scannerLPS.stl" % tag), to_scanner(V), Fc)
        Vl = V.copy(); Vl[:, 2] += z0
        write_stl(os.path.join(args.out, "sample2_vessel_%s_local.stl" % tag), Vl, Fc)
        # Pappus check: volume of a solid of revolution = 2 pi * integral of rho dA (shoelace on the polygon)
        x_, z_ = poly[:, 0], poly[:, 1]; x2, z2 = np.roll(x_, -1), np.roll(z_, -1)
        cross = x_ * z2 - x2 * z_; A_ = 0.5 * cross.sum()
        Vp = abs(np.pi / 3.0 * ((x_ ** 2 + x_ * x2 + x2 ** 2) * (z2 - z_)).sum())   # exact for straight-edged polygon
        exp[tag] = dict(triangles=int(len(Fc)), watertight=is_watertight(Fc), mesh_volume_cm3=round(abs(signed_volume(V, Fc)) / 1000, 1),
                        polygon_volume_cm3=round(float(Vp) / 1000, 1))
    R["export"] = exp
    # smooth CAD description of the same profile (splines between the sharp corners, lines elsewhere)
    Po_d = thin_(o_rho, o_h_d, 240); Pi_d = thin_(i_rho, i_h_d, 240)
    poly_meas_d = np.vstack([Po_d, Pi_d[::-1]])
    Co_d = thin_(clos_o_rho, clos_o_h, 60); Ci_d = thin_(clos_i_rho, clos_i_h, 60)
    poly_full_d = np.vstack([Co_d[:-1], Po_d, Pi_d[::-1], Ci_d[::-1][1:]])
    pieces_meas = cad_pieces(poly_meas_d); pieces_full = cad_pieces(poly_full_d)
    exp["cad_pieces"] = dict(measured=len(pieces_meas), with_base=len(pieces_full))
    # 2-D profile for a native CAD revolve (X = radius, Y = height above the lowest surviving point; units mm)
    try:
        import ezdxf
        doc = ezdxf.new("R2010", setup=True); doc.units = ezdxf.units.MM; msp = doc.modelspace()
        for nm, col in (("PROFILE_MEASURED", 7), ("AXIS", 1), ("ILLUSTRATIVE_BASE_CAP", 8)):
            doc.layers.add(nm, color=col)
        for pc in pieces_meas:
            if pc[0] == "line": msp.add_line(tuple(pc[1]), tuple(pc[2]), dxfattribs={"layer": "PROFILE_MEASURED"})
            else: msp.add_spline([tuple(p) for p in pc[1]], degree=3, dxfattribs={"layer": "PROFILE_MEASURED"})
        ymin, ymax = float(min(poly_full[:, 1].min(), 0)) - 10, float(o_h_d.max()) + 10
        msp.add_line((0, ymin), (0, ymax), dxfattribs={"layer": "AXIS", "linetype": "CENTER"})
        for cap in (Co_d, Ci_d):
            msp.add_spline([tuple(p) for p in cap[::6]] + [tuple(cap[-1])], degree=3, dxfattribs={"layer": "ILLUSTRATIVE_BASE_CAP", "linetype": "DASHED"})
        doc.saveas(os.path.join(args.out, "sample2_vessel_profile.dxf")); exp["dxf"] = True
    except Exception as ex:                                                            # optional dependency
        exp["dxf"] = "skipped: %s" % ex
    # exact solids (splines revolved 360 degrees) for SolidWorks / any CAD system
    try:
        import cadquery as cq
        def solid(pieces, dz=0.0):
            V3 = lambda p: cq.Vector(float(p[0]), 0.0, float(p[1]) + dz)
            edges = []
            for pc in pieces:
                if pc[0] == "line":
                    if np.linalg.norm(pc[1] - pc[2]) > 1e-6: edges.append(cq.Edge.makeLine(V3(pc[1]), V3(pc[2])))
                else:
                    edges.append(cq.Edge.makeSpline([V3(p) for p in pc[1]]))
            return cq.Solid.revolve(cq.Wire.assembleEdges(edges), [], 360.0, cq.Vector(0, 0, 0), cq.Vector(0, 0, 1))
        sol = solid(pieces_meas)
        cq.exporters.export(cq.Workplane("XY").add(sol), os.path.join(args.out, "sample2_vessel_measured.step"))
        exp["step_measured"] = dict(volume_cm3=round(sol.Volume() / 1000, 1), valid=bool(sol.isValid()))
        sol2 = solid(pieces_full, dz=-float(poly_full[:, 1].min()))
        cq.exporters.export(cq.Workplane("XY").add(sol2), os.path.join(args.out, "sample2_vessel_with_illustrative_base.step"))
        exp["step_with_base"] = dict(volume_cm3=round(sol2.Volume() / 1000, 1), valid=bool(sol2.isValid()))
    except Exception as ex:
        exp["step"] = "skipped: %s" % ex
    print("[5] export:", json.dumps(exp))

    with open(os.path.join(args.out, "results.json"), "w") as fh:
        json.dump(R, fh, indent=2)
    print(json.dumps(R, indent=2))

if __name__ == "__main__":
    main()
