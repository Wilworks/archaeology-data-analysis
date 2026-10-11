#!/usr/bin/env python3
"""
make_figures.py - regenerates the figures that are NOT written by the reconstruction scripts themselves.

Figures produced (file names in --out)
  hu_histograms_samples1to4.png            HU histograms with the split thresholds, and the high-density tails
  surface_views_sample{1,2,3,4}.png        six shaded views of each specimen (500 HU surface)
  sample1_parts.png                        Sample 1 cluster split into its ceramic pieces (and the brown piece)
  sample1_piece_isolated_vs_cleaned.png    the brown piece alone, before/after removing debris
  sample1_piece_on_its_base.png            the cleaned brown piece standing on its foot (from the reconstruction STL)
  sample3_height_maps.png                  lower face, upper face and extent maps of the Sample 3 plate
  sample3_scan_vs_reconstruction.png       scan against reconstruction from identical viewpoints, with distance map
  sample3_orientation_A_vs_B.png           the two possible orientations of the Sample 3 base
  cross_sample_profiles.png                half cross-sections of all reconstructions on one scale

Usage
  python make_figures.py --npz-dir ../01_inputs/point_files --work WORK --out OUT [--only NAME ...]
  WORK is the folder where run_all.py placed the reconstruction outputs (subfolders s1, s2, s3, s3_flip, s4).
Requirements: numpy, scipy, scikit-image, matplotlib, trimesh (for the Sample 3 comparison only).
"""
import argparse, json, os
import numpy as np
from scipy import ndimage as ndi
from scipy.spatial import cKDTree
from skimage import measure
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.collections import PolyCollection
import sample2_reconstruction as S
import sample4_wall_reconstruction as W

PALE, BROWN = np.array([0.88, 0.82, 0.74]), np.array([0.68, 0.43, 0.20])


def paint(ax, T3, nrm, cols, cam, upref, title, ctr):
    cam = cam / np.linalg.norm(cam); right = np.cross(upref, cam)
    right = right / np.linalg.norm(right) if np.linalg.norm(right) > 1e-6 else np.array([1.0, 0, 0]); up = np.cross(cam, right)
    Q = T3 - ctr; sx, sy, dz = Q @ right, Q @ up, Q @ cam
    order = np.argsort(dz.mean(1)); lv = 0.6 * cam + 0.4 * np.array([0.3, -0.4, 0.85]); lv /= np.linalg.norm(lv)
    shade = 0.35 + 0.65 * np.clip(np.abs(nrm @ lv), 0, 1)
    ax.add_collection(PolyCollection(np.stack([sx, sy], -1)[order], facecolors=np.clip(cols * shade[:, None], 0, 1)[order], edgecolors="none"))
    ax.set_xlim(sx.min() - 3, sx.max() + 3); ax.set_ylim(sy.min() - 3, sy.max() + 3); ax.set_aspect("equal"); ax.axis("off"); ax.set_title(title, fontsize=9.5)


def tri_normals(T3):
    n = np.cross(T3[:, 1] - T3[:, 0], T3[:, 2] - T3[:, 0]); return n / (np.linalg.norm(n, axis=1, keepdims=True) + 1e-12)


def load(npz_dir, k):
    HU, sp, origin, hu = S.load_grid(os.path.join(npz_dir, "sample%d.npz" % k)); return HU, sp, origin, hu


def ceramic_mask(k, HU, sp):
    """The 'pale' (ceramic) part used for colouring each specimen."""
    if k == 2: return ndi.binary_fill_holes(S.largest_component(HU >= 1212))
    if k == 3: return S.largest_component(HU >= 1200)
    return W.smoothed_body(HU, sp)[1]


# ---------------------------------------------------------------------------------------------------------------------
def fig_hu_histograms(a):
    thr = {2: 1212.0, 3: 1200.0, 4: 1112.0, 1: 1088.0}
    try: thr[4] = json.load(open(os.path.join(a.work, "s4", "results.json")))["split"]["threshold_HU"]
    except Exception: pass
    try: thr[1] = json.load(open(os.path.join(a.work, "s1", "results.json")))["isolation"]["split_threshold_HU"]
    except Exception: pass
    fig, axs = plt.subplots(2, 4, figsize=(20, 8))
    for j, k in enumerate((1, 2, 3, 4)):
        hu = np.load(os.path.join(a.npz_dir, "sample%d.npz" % k))["hu"].astype(int)
        e = np.arange(500, 3125, 25); h, _ = np.histogram(hu, bins=e); c = 0.5 * (e[:-1] + e[1:])
        axs[0, j].bar(c, h, width=25, color="#3a5f8a"); axs[0, j].axvline(thr[k], color="#d03030", ls="--", lw=1)
        axs[0, j].set_xlim(500, 2000); axs[0, j].set_title("Sample %d: %d voxels, mean %.0f HU" % (k, len(hu), hu.mean()), fontsize=10)
        axs[0, j].text(thr[k] + 15, h.max() * 0.9, "split %.0f HU" % thr[k], color="#d03030", fontsize=8); axs[0, j].set_xlabel("HU"); axs[0, j].set_ylabel("voxels per 25 HU")
        hc = np.maximum(h, 0.5); axs[1, j].bar(c, hc, width=25, color="#6b6b6b"); axs[1, j].set_yscale("log"); axs[1, j].set_xlim(1400, 3125); axs[1, j].set_ylim(0.8, hc.max() * 1.5)
        axs[1, j].axvline(3071, color="#d03030", ls="--", lw=1); axs[1, j].text(3060, hc.max() * 0.5, "12-bit cap\n%d voxels" % int((hu == 3071).sum()), ha="right", fontsize=8, color="#d03030")
        axs[1, j].set_xlabel("HU"); axs[1, j].set_ylabel("voxels (log)"); axs[1, j].set_title("high-density tail", fontsize=10)
    fig.suptitle("HU distributions of the four scanned specimens (voxels at or above 500 HU)", fontsize=12, fontweight="bold")
    fig.tight_layout(rect=(0, 0, 1, 0.95)); fig.savefig(os.path.join(a.out, "hu_histograms_samples1to4.png"), dpi=100); plt.close(fig)


def fig_surface_views(a, k):
    HU, sp, origin, hu = load(a.npz_dir, k); present = HU >= 500; body = ceramic_mask(k, HU, sp)
    step = 2 if k == 1 else 1
    F = np.where(present, HU.astype(float), 0.0); Fp = np.pad(ndi.gaussian_filter(F, 0.8), 3)
    v, f, n, _ = measure.marching_cubes(Fp, level=500.0, spacing=tuple(sp), step_size=step); v = v - 3 * sp + origin
    Xb = np.argwhere(body) * sp + origin; cb = Xb.mean(0); Vt = np.linalg.svd(Xb - cb, full_matrices=False)[2]
    if np.linalg.det(Vt) < 0: Vt[2] *= -1
    sgn = 1.0
    if k != 3:
        Xc = np.argwhere(present & ~ndi.binary_dilation(body, iterations=2)) * sp + origin
        if len(Xc): sgn = float(np.sign((Xc.mean(0) - cb) @ Vt[2])) or 1.0
    P = (v - cb) @ Vt.T; P[:, 2] *= sgn
    cen = v[f].mean(1); ijk = np.clip(np.rint((cen - origin) / sp).astype(int), 0, np.array(HU.shape) - 1)
    is_body = ndi.binary_dilation(body, iterations=2)[ijk[:, 0], ijk[:, 1], ijk[:, 2]]
    cols = np.where(is_body[:, None], PALE, BROWN); T3 = P[f]; nrm = tri_normals(T3); ctr = 0.5 * (P.max(0) + P.min(0))
    fig, axs = plt.subplots(2, 3, figsize=(17, 10))
    for ax, (t, el, az) in zip(axs.ravel(), [("plan, one face", 90, -90), ("plan, other face", -90, -90), ("side view", 0, -90), ("edge-on", 0, 0), ("oblique, one side", 35, -60), ("oblique, other side", -35, 30)]):
        t_, p_ = np.radians(el), np.radians(az); cam = np.array([np.cos(t_) * np.cos(p_), np.cos(t_) * np.sin(p_), np.sin(t_)])
        paint(ax, T3, nrm, cols, cam, np.array([0, 0, 1.0]), t, ctr)
    lab = {1: "Sample 1 (scan label 2; pale = ceramic group, brown = the brown piece and crust)", 2: "Sample 2 (scan label 424; pale = ceramic body, brown = crust)",
           3: "Sample 3 (scan label 62; almost no crust)", 4: "Sample 4 (scan label 440; pale = ceramic body, brown = crust)"}[k]
    fig.suptitle(lab + ": 500 HU surface", fontsize=11); fig.tight_layout(rect=(0, 0, 1, 0.95))
    fig.savefig(os.path.join(a.out, "surface_views_sample%d.png" % k), dpi=70); plt.close(fig)


def fig_s3_height_maps(a):
    HU, sp, origin, hu = load(a.npz_dir, 3); core = S.largest_component(HU >= 1200)
    X = np.argwhere(core) * sp + origin; cb = X.mean(0); Vt = np.linalg.svd(X - cb, full_matrices=False)[2]
    if np.linalg.det(Vt) < 0: Vt[2] *= -1
    Q = (X - cb) @ Vt.T; g = 4.0
    i1 = np.floor((Q[:, 0] - Q[:, 0].min()) / g).astype(int); i2 = np.floor((Q[:, 1] - Q[:, 1].min()) / g).astype(int)
    lo = np.full((i1.max() + 1, i2.max() + 1), np.nan); hi = lo.copy()
    for aa in range(lo.shape[0]):
        m = i1 == aa
        for b in range(lo.shape[1]):
            mm = m & (i2 == b)
            if mm.sum() >= 8: lo[aa, b] = np.percentile(Q[mm, 2], 2); hi[aa, b] = np.percentile(Q[mm, 2], 98)
    ext = [Q[:, 0].min(), Q[:, 0].max(), Q[:, 1].min(), Q[:, 1].max()]
    fig, ax = plt.subplots(1, 3, figsize=(18, 5.5))
    for a_, arr, t in zip(ax, (lo, hi, hi - lo), ("lower face height u3 (mm)", "upper face height u3 (mm)", "extent across the plate (mm)")):
        im = a_.imshow(arr.T, origin="lower", extent=ext, cmap="viridis"); plt.colorbar(im, ax=a_, fraction=0.046); a_.set_title(t, fontsize=10); a_.set_xlabel("u1 (mm)"); a_.set_ylabel("u2 (mm)")
    fig.suptitle("Sample 3: height maps of the two faces in the plate's own frame (u3 = across the plate)", fontsize=11)
    fig.tight_layout(rect=(0, 0, 1, 0.94)); fig.savefig(os.path.join(a.out, "sample3_height_maps.png"), dpi=90); plt.close(fig)


def fig_s3_scan_vs_recon(a):
    import trimesh
    HU, sp, origin, hu = load(a.npz_dir, 3); present = HU >= 500; core = S.largest_component(HU >= 1200)
    Vo, Fo, No, Wo = S.mesh_from_field(np.where(present, HU.astype(float), 0.0), 500.0, sp, origin, sigma=0.8)
    near = cKDTree(np.argwhere(core) * sp + origin).query(Vo)[0] < 1.6
    rec = trimesh.load(os.path.join(a.work, "s3", "sample3_base_scannerLPS.stl")); R0 = json.load(open(os.path.join(a.work, "s3", "results.json")))
    ax_ = np.array(R0["axis"]["direction"]); c = np.array(R0["axis"]["point_mm"]); e1, e2 = S.centred_basis(ax_, c, Vo[near])
    idx = np.where(near)[0]; rng = np.random.default_rng(0); sel = idx[rng.choice(len(idx), 14000, replace=False)]
    cp, dist, tri = trimesh.proximity.closest_point(rec, Vo[sel]); sgn = np.where(np.einsum("ij,ij->i", Vo[sel] - cp, rec.face_normals[tri]) > 0, 1.0, -1.0); sd = sgn * dist
    tree = cKDTree(Vo[sel]); nn = tree.query(Vo)[1]; sdf = sd[nn][Fo].mean(1); fk = near[Fo].all(1)
    cols_d = plt.get_cmap("coolwarm")(np.clip((sdf + 4) / 8, 0, 1))[:, :3]; cols_d[~fk] = 0.75
    cols_s = np.where(fk[:, None], PALE, BROWN); Ts = Vo[Fo]; ns = tri_normals(Ts); Tr = rec.vertices[rec.faces]; nr = tri_normals(Tr)
    cols_r = np.tile(np.array([0.62, 0.71, 0.80]), (len(Tr), 1))
    fig, axs = plt.subplots(3, 3, figsize=(16, 15))
    for j, (t, cam, upr) in enumerate((("looking down the axis", ax_, e1), ("side view, axis vertical", e1, ax_), ("oblique from above", 0.55 * e1 + 0.4 * e2 + 0.75 * ax_, ax_))):
        paint(axs[0, j], Ts, ns, cols_s, cam, upr, "SCAN, " + t + "\n(brown = crust flecks)", c); paint(axs[1, j], Tr, nr, cols_r, cam, upr, "RECONSTRUCTION, " + t, c)
        paint(axs[2, j], Ts, ns, cols_d, cam, upr, "SCAN coloured by distance from the reconstruction\n(red = sticks out, blue = inside, grey = crust)", c)
    sm = plt.cm.ScalarMappable(cmap="coolwarm", norm=plt.Normalize(-4, 4)); fig.colorbar(sm, ax=axs[2, :], fraction=0.02, label="mm")
    fig.suptitle("Sample 3: scan versus reconstruction, same camera (median distance %.2f mm; within 1/2/3 mm: %.0f/%.0f/%.0f%%)" %
                 (np.median(dist), 100 * (dist < 1).mean(), 100 * (dist < 2).mean(), 100 * (dist < 3).mean()), fontsize=12, fontweight="bold")
    fig.subplots_adjust(left=0.01, right=0.93, top=0.94, bottom=0.01, wspace=0.02, hspace=0.12); fig.savefig(os.path.join(a.out, "sample3_scan_vs_reconstruction.png"), dpi=70); plt.close(fig)
    json.dump(dict(median_distance_mm=round(float(np.median(dist)), 2), within_1mm_percent=round(100 * float((dist < 1).mean())), within_2mm_percent=round(100 * float((dist < 2).mean())), within_3mm_percent=round(100 * float((dist < 3).mean())),
              signed_median_mm=round(float(np.median(sd)), 2), n_points=int(len(dist))), open(os.path.join(a.out, 'sample3_scan_vs_reconstruction_stats.json'), 'w'), indent=1)


def fig_s3_orientation(a):
    pa = np.loadtxt(os.path.join(a.work, "s3_flip", "profile_base.csv"), delimiter=",", skiprows=1)      # the rejected orientation (--flip)
    pb = np.loadtxt(os.path.join(a.work, "s3", "profile_base.csv"), delimiter=",", skiprows=1)            # the accepted orientation
    fig, axs = plt.subplots(1, 2, figsize=(15, 4.6))
    for ax, p, ttl in zip(axs, (pa, pb), ("A (rejected): stands on the wide flat face; the peaked crest is on top", "B (accepted): stands on the rounded ridge; the wall root bends upward at the outer edge")):
        r, up, lo = p[:, 0], p[:, 1], p[:, 2]
        for sg in (1, -1):
            ax.fill(sg * np.concatenate([r, r[::-1]]), np.concatenate([up, lo[::-1]]), color="#d9bfa2"); ax.plot(sg * r, up, "k", lw=1.2); ax.plot(sg * r, lo, "k", lw=1.2)
        ax.axhline(0, color="#3a5f8a", lw=2); ax.axvline(0, color="gray", ls="-.", lw=0.8); ax.set_aspect("equal"); ax.set_ylim(-6, 32); ax.set_xlabel("mm"); ax.set_title(ttl, fontsize=9.5)
    fig.tight_layout(); fig.savefig(os.path.join(a.out, "sample3_orientation_A_vs_B.png"), dpi=95); plt.close(fig)


def fig_s1_parts(a):
    from skimage.segmentation import watershed
    HU, sp, origin, hu = load(a.npz_dir, 1); present, body, T, _, _ = W.smoothed_body(HU, sp)
    st = ndi.generate_binary_structure(3, 1); er = ndi.binary_erosion(body, structure=st, iterations=3)
    lab, n = ndi.label(er, structure=np.ones((3, 3, 3))); sizes = np.bincount(lab.ravel())[1:]; markers = np.zeros_like(lab)
    for j, kk in enumerate(np.where(sizes >= 3000)[0] + 1, 1): markers[lab == kk] = j
    seg = watershed(-ndi.distance_transform_edt(body, sampling=tuple(sp)), markers, mask=body)
    F = np.where(present, HU.astype(float), 0.0); v, f, nn, _ = measure.marching_cubes(np.pad(ndi.gaussian_filter(F, 0.8), 3), level=500.0, spacing=tuple(sp), step_size=2); v = v - 3 * sp + origin
    cen = v[f].mean(1); ijk = np.clip(np.rint((cen - origin) / sp).astype(int), 0, np.array(HU.shape) - 1); labt = ndi.grey_dilation(seg, size=7)[ijk[:, 0], ijk[:, 1], ijk[:, 2]]
    counts = [(seg == j).sum() for j in range(1, seg.max() + 1)]; order = np.argsort(counts)[::-1]; rank = {int(j + 1): r for r, j in enumerate(order)}
    pal = {0: BROWN, 1: np.array([0.30, 0.52, 0.78]), 2: np.array([0.35, 0.72, 0.42]), 3: np.array([0.90, 0.65, 0.20])}
    cols = np.array([pal[0] if l == 0 else pal[min(rank[int(l)] + 1, 3)] for l in labt])
    Xb = np.argwhere(body) * sp + origin; cb = Xb.mean(0); Vt = np.linalg.svd(Xb - cb, full_matrices=False)[2]
    if np.linalg.det(Vt) < 0: Vt[2] *= -1
    P = (v - cb) @ Vt.T; T3 = P[f]; nrm = tri_normals(T3); ctr = 0.5 * (P.max(0) + P.min(0))
    fig, axs = plt.subplots(1, 4, figsize=(20, 5.6))
    for ax, (t, el, az) in zip(axs, [("plan, one side", 90, -90), ("plan, other side", -90, -90), ("side view", 0, -90), ("oblique", 35, -60)]):
        t_, p_ = np.radians(el), np.radians(az); cam = np.array([np.cos(t_) * np.cos(p_), np.cos(t_) * np.sin(p_), np.sin(t_)]); paint(ax, T3, nrm, cols, cam, np.array([0, 0, 1.0]), t, ctr)
    vol = {rank[j]: (seg == j).sum() * np.prod(sp) / 1000 for j in range(1, seg.max() + 1)}
    fig.suptitle("Sample 1 (scan label 2) split into its parts: blue main slab %.0f cm3, green %.0f cm3, orange %.0f cm3; brown = the brown piece and debris" % (vol[0], vol[1], vol[2]), fontsize=11)
    fig.tight_layout(rect=(0, 0, 1, 0.93)); fig.savefig(os.path.join(a.out, "sample1_parts.png"), dpi=70); plt.close(fig)


def fig_s1_isolation(a):
    sys_path_ok = True
    import sample1_base_reconstruction as B1
    HU, sp, origin, hu = load(a.npz_dir, 1); present, body, rec, clean, T = B1.isolate_piece(HU, sp)
    def mesh(mask):
        v, f, n, _ = measure.marching_cubes(np.pad(ndi.gaussian_filter(mask.astype(float), 1.0), 3), level=0.5, spacing=tuple(sp), step_size=2); return v - 3 * sp + origin, f
    Vr, Fr = mesh(rec); Vc, Fc = mesh(clean); X = np.argwhere(clean) * sp + origin; c = X.mean(0); Vt = np.linalg.svd(X - c, full_matrices=False)[2]
    e1, e2, ax_ = Vt[0], Vt[1], Vt[2]
    views = [("along the disc axis, side A", ax_, e1), ("along the disc axis, side B", -ax_, e1), ("side view", e1, ax_), ("oblique", 0.6 * e1 + 0.4 * e2 + 0.7 * ax_, ax_)]
    fig, axs = plt.subplots(2, 4, figsize=(20, 9.6))
    for j, (t, cam, upr) in enumerate(views):
        for i, (V_, F_, nm) in enumerate(((Vr, Fr, "ISOLATED (debris still on)"), (Vc, Fc, "CLEANED (lumps under about 10 mm removed)"))):
            T3 = V_[F_]; paint(axs[i, j], T3, tri_normals(T3), np.tile(np.array([0.72, 0.47, 0.24]), (len(T3), 1)), cam, upr, nm + ", " + t, c)
    fig.suptitle("Sample 1: the brown piece alone (the ceramic pieces removed)", fontsize=12, fontweight="bold")
    fig.tight_layout(rect=(0, 0, 1, 0.95)); fig.savefig(os.path.join(a.out, "sample1_piece_isolated_vs_cleaned.png"), dpi=65); plt.close(fig)


def fig_s1_on_base(a):
    import trimesh
    m = trimesh.load(os.path.join(a.work, "s1", "sample1_piece_cleaned_on_its_base_local.stl")); X = np.asarray(m.vertices); F = np.asarray(m.faces); T3 = X[F]; nrm = tri_normals(T3)
    cols = np.tile(np.array([0.84, 0.72, 0.58]), (len(T3), 1)); ctr = np.array([0, 0, 25.0]); z = np.array([0, 0, 1.0])
    fig, axs = plt.subplots(1, 4, figsize=(20, 5.4))
    for ax, (t, cam, upr) in zip(axs, (("side view, standing on its round base", np.array([1.0, 0, 0]), z), ("side view from the other direction", np.array([0, 1.0, 0]), z),
                                        ("oblique from above", np.array([0.5, 0.4, 0.75]), z), ("looking up at the base", np.array([0, 0, -1.0]), np.array([0, 1.0, 0])))):
        paint(ax, T3, nrm, cols, cam, upr, t, ctr)
    axs[0].axhline(-ctr[2], color="#3a5f8a", lw=2)
    fig.suptitle("Sample 1: the brown piece, debris removed, standing on its round foot", fontsize=12, fontweight="bold")
    fig.tight_layout(rect=(0, 0, 1, 0.93)); fig.savefig(os.path.join(a.out, "sample1_piece_on_its_base.png"), dpi=72); plt.close(fig)


def fig_cross_profiles(a):
    fig, axs = plt.subplots(1, 4, figsize=(21, 5.4), sharey=False)
    p2 = np.loadtxt(os.path.join(a.work, "s2", "profile_sample2.csv"), delimiter=",", skiprows=1)
    axs[0].fill(np.concatenate([p2[:, 1], p2[::-1, 3]]), np.concatenate([p2[:, 2], p2[::-1, 4]]), color="#d9bfa2"); axs[0].plot(p2[:, 1], p2[:, 2], "k", lw=1.2); axs[0].plot(p2[:, 3], p2[:, 4], "k", lw=1.2)
    axs[0].set_title("Sample 2: rim to shoulder (bowl)", fontsize=10)
    p3 = np.loadtxt(os.path.join(a.work, "s3", "profile_base.csv"), delimiter=",", skiprows=1)
    axs[1].fill(np.concatenate([p3[:, 0], p3[::-1, 0]]), np.concatenate([p3[:, 1], p3[::-1, 2]]), color="#d9bfa2"); axs[1].plot(p3[:, 0], p3[:, 1], "k", lw=1.2); axs[1].plot(p3[:, 0], p3[:, 2], "k", lw=1.2)
    axs[1].set_title("Sample 3: ring foot (base)", fontsize=10)
    p1 = np.loadtxt(os.path.join(a.work, "s1", "profile_base.csv"), delimiter=",", skiprows=1)
    axs[2].fill(np.concatenate([p1[:, 0], p1[::-1, 0]]), np.concatenate([p1[:, 1], p1[::-1, 2]]), color="#d9bfa2"); axs[2].plot(p1[:, 0], p1[:, 1], "k", lw=1.2); axs[2].plot(p1[:, 0], p1[:, 2], "k", lw=1.2)
    axs[2].set_title("Sample 1: pedestal dish (base)", fontsize=10)
    p4 = np.loadtxt(os.path.join(a.work, "s4", "profile_wall.csv"), delimiter=",", skiprows=1)
    axs[3].fill(np.concatenate([p4[:, 1], p4[::-1, 2]]), np.concatenate([p4[:, 0], p4[::-1, 0]]), color="#d9bfa2"); axs[3].plot(p4[:, 1], p4[:, 0], "k", lw=1.2); axs[3].plot(p4[:, 2], p4[:, 0], "k", lw=1.2)
    axs[3].set_title("Sample 4: wall band (up direction assumed)", fontsize=10)
    for ax in axs: ax.set_aspect("equal"); ax.set_xlabel("distance from the axis (mm)"); ax.set_ylabel("height (mm)"); ax.grid(alpha=0.3)
    fig.suptitle("Measured half cross-sections (each panel 1:1, distance from its own axis)", fontsize=12, fontweight="bold")
    fig.tight_layout(rect=(0, 0, 1, 0.94)); fig.savefig(os.path.join(a.out, "cross_sample_profiles.png"), dpi=95); plt.close(fig)


JOBS = {"hu": fig_hu_histograms, "views1": lambda a: fig_surface_views(a, 1), "views2": lambda a: fig_surface_views(a, 2), "views3": lambda a: fig_surface_views(a, 3),
        "views4": lambda a: fig_surface_views(a, 4), "heightmaps3": fig_s3_height_maps, "scanrecon3": fig_s3_scan_vs_recon, "orient3": fig_s3_orientation,
        "parts1": fig_s1_parts, "isolation1": fig_s1_isolation, "onbase1": fig_s1_on_base, "profiles": fig_cross_profiles}

if __name__ == "__main__":
    ap = argparse.ArgumentParser(); ap.add_argument("--npz-dir", required=True); ap.add_argument("--work", required=True); ap.add_argument("--out", required=True)
    ap.add_argument("--only", nargs="*", help="subset of: " + ", ".join(JOBS)); a = ap.parse_args(); os.makedirs(a.out, exist_ok=True)
    for name in (a.only or JOBS):
        try:
            JOBS[name](a); print("done:", name)
        except Exception as ex:
            print("FAILED:", name, "->", repr(ex))
