#!/usr/bin/env python3
"""
Show the REAL scanned sherd (crust removed, broken edges intact) inside its symmetric completion.

For each sample it writes
  * <sample>_real_sherd_cleaned_scannerLPS.stl : the ceramic body only (surface of the body mask grown by one voxel, which puts
    it near the true edge); where crust was attached, the surface is the ceramic-crust interface. An open/closed shell in
    scanner coordinates, to be loaded together with the vessel STL.
  * fig_<sample>_real_sherd_in_vessel.png      : clay = real sherd, blue = symmetric completion, grey = not measured

Usage: python show_real_sherd_in_vessel.py --sample 2 --npz sample2.npz --vessel sample2_vessel_with_illustrative_base_scannerLPS.stl --results sample2_results.json
       python show_real_sherd_in_vessel.py --sample 4 --npz sample4.npz --vessel sample4_wall_band_scannerLPS.stl --results sample4_wall_results.json
Needs sample2_reconstruction.py and sample4_wall_reconstruction.py in the same folder. Requirements: numpy, scipy, scikit-image, matplotlib, trimesh.
"""
import argparse, json, os
import numpy as np
from scipy import ndimage as ndi
from scipy.spatial import cKDTree
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.collections import PolyCollection
import trimesh
import sample2_reconstruction as S
import sample4_wall_reconstruction as W


def body_mask(sample, HU, sp):
    if sample == 2:
        return ndi.binary_fill_holes(S.largest_component(HU >= 1212))
    if sample == 4:
        return W.smoothed_body(HU, sp)[1]
    if sample == 3:
        return S.largest_component(HU >= 1200)
    raise ValueError("sample must be 2, 3 or 4")


def paint(ax, T3, nrm, cols, cam, upref, ttl, ctr):
    cam = cam / np.linalg.norm(cam); right = np.cross(upref, cam); right /= np.linalg.norm(right); upv = np.cross(cam, right)
    Q = T3 - ctr; sx, sy, dz = Q @ right, Q @ upv, Q @ cam
    order = np.argsort(dz.mean(1)); lv = 0.6 * cam + 0.4 * np.array([0.3, -0.4, 0.85]); lv /= np.linalg.norm(lv)
    shade = 0.4 + 0.6 * np.clip(np.abs(nrm @ lv), 0, 1)
    ax.add_collection(PolyCollection(np.stack([sx, sy], -1)[order], facecolors=np.clip(cols * shade[:, None], 0, 1)[order], edgecolors="none"))
    ax.set_xlim(sx.min() - 4, sx.max() + 4); ax.set_ylim(sy.min() - 4, sy.max() + 4); ax.set_aspect("equal"); ax.axis("off"); ax.set_title(ttl, fontsize=9.5)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--sample", type=int, required=True); ap.add_argument("--npz", required=True)
    ap.add_argument("--vessel", required=True); ap.add_argument("--results", required=True); ap.add_argument("--out", default="out_real")
    ap.add_argument("--unmeasured-radius", type=float, default=None, help="radius (mm) inside which the vessel is NOT measured (drawn grey)")
    args = ap.parse_args(); os.makedirs(args.out, exist_ok=True)
    HU, sp, origin, hu = S.load_grid(args.npz)
    body = body_mask(args.sample, HU, sp)
    grown = ndi.binary_dilation(body, iterations=1)
    V, F, N, Wt = S.mesh_from_field(grown.astype(float), 0.5, sp, origin, sigma=1.0)
    name = "sample%d" % args.sample
    S.write_stl(os.path.join(args.out, "%s_real_sherd_cleaned_scannerLPS.stl" % name), V, F)
    R0 = json.load(open(args.results)); a = np.array(R0["axis"]["direction"]); c = np.array(R0["axis"]["point_mm"])
    vs = trimesh.load(args.vessel); Vv, Fv = np.asarray(vs.vertices), np.asarray(vs.faces)
    e1, e2 = S.centred_basis(a, c, V)
    q = V - c; phi_s = np.angle(np.mean(np.exp(1j * np.arctan2(q @ e2, q @ e1))))
    Tv = Vv[Fv]; cen = Tv.mean(1)
    tree = cKDTree(V); hide = tree.query(cen)[0] < 4.0
    qv = cen - c; hv = qv @ a; rv = np.linalg.norm(qv - hv[:, None] * a, axis=1); phi_v = np.arctan2(qv @ e2, qv @ e1)
    cut = np.abs((phi_v - phi_s + np.pi) % (2 * np.pi) - np.pi) > np.radians(130)
    nv = np.cross(Tv[:, 1] - Tv[:, 0], Tv[:, 2] - Tv[:, 0]); nv /= np.linalg.norm(nv, axis=1, keepdims=True) + 1e-12
    Ts = V[F]; ns = np.cross(Ts[:, 1] - Ts[:, 0], Ts[:, 2] - Ts[:, 0]); ns /= np.linalg.norm(ns, axis=1, keepdims=True) + 1e-12
    col_v = np.tile(np.array([0.62, 0.71, 0.80]), (len(Tv), 1))
    if args.unmeasured_radius:
        col_v[rv < args.unmeasured_radius] = np.array([0.72, 0.72, 0.72])
    col_s = np.tile(np.array([0.83, 0.68, 0.52]), (len(Ts), 1)); ctr = 0.5 * (Tv.reshape(-1, 3).max(0) + Tv.reshape(-1, 3).min(0))
    away = -(np.cos(phi_s) * e1 + np.sin(phi_s) * e2)          # direction pointing away from the sherd
    toward = -away
    views = [("oblique from above", 0.55 * toward + 0.45 * e2 * 0 + 0.75 * a, a, False),
             ("side view, axis vertical", toward, a, False),
             ("cutaway from above (far side removed)", 0.35 * toward + 0.9 * a, a, True)]
    fig, axs = plt.subplots(1, 3, figsize=(18, 6.2))
    for ax_, (ttl, cam, upr, use_cut) in zip(axs, views):
        keep = ~hide & (~cut if use_cut else True)
        T3 = np.concatenate([Tv[keep], Ts]); nn = np.concatenate([nv[keep], ns]); cc = np.concatenate([col_v[keep], col_s])
        paint(ax_, T3, nn, cc, cam, upr, ttl + "\n(clay = REAL scanned sherd, blue = symmetric completion%s)" % (", grey = not measured" if args.unmeasured_radius else ""), ctr)
    fig.suptitle("Sample %d: the real broken sherd shown inside the completed vessel" % args.sample, fontsize=12, fontweight="bold")
    fig.tight_layout(rect=(0, 0, 1, 0.94)); fig.savefig(os.path.join(args.out, "fig_%s_real_sherd_in_vessel.png" % name), dpi=85); plt.close(fig)
    print("sample %d: cleaned sherd %d vertices; vessel triangles hidden behind the real sherd: %.1f%%" % (args.sample, len(V), 100 * hide.mean()))


if __name__ == "__main__":
    main()
