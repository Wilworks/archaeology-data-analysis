#!/usr/bin/env python3
"""
Reconstruction of a rotationally symmetric vessel BASE (foot ring on a thin recessed floor) from one CT-scanned base sherd
(used for the piece with scan label 62, called Sample 3 in the manuscript's Table 2).

Input  : sample3.npz  (x, y, z in mm; hu int16; every voxel of the piece at >= 500 HU)
Output : results.json, profile_base.csv, fig_base_*.png and CAD/mesh exports (*.stl, *.step, *.dxf) in --out
Needs sample2_reconstruction.py in the same folder (shared helper functions).
Requirements: numpy, scipy, scikit-image, matplotlib (+ optional ezdxf, cadquery). Fixed seeds: repeatable.

Method in one paragraph: the ceramic core (>= 1200 HU) is meshed at the 500 HU edge. The two faces of the plate are
separated by their normals. For a candidate axis every face point becomes (rho, h); when the axis is right, each face is a
single thin curve h(rho), so the axis is found by minimising the spread of h inside 1 mm bins of rho (both faces). The
profile is then the bin medians, moved to the half-maximum edge, and swept 360 degrees.
Orientation: the default stands the base on its foot ring with the wall root bending upward at the outer edge (confirmed by
the owner of the sherd); --flip reverses it (the earlier, rejected orientation).
"""
import argparse, json, os
import numpy as np
from scipy import ndimage as ndi
from scipy.optimize import minimize
from scipy.spatial import cKDTree
from scipy.signal import medfilt
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.collections import PolyCollection
import sample2_reconstruction as S


def bin_rms(rho, h, width=1.0):
    """RMS spread of h inside bins of rho (a curve h(rho) has zero spread apart from noise)."""
    k = np.floor(rho / width).astype(int); k -= k.min()
    n = np.bincount(k); s1 = np.bincount(k, h); s2 = np.bincount(k, h * h)
    ok = n >= 5
    var = s2[ok] / n[ok] - (s1[ok] / n[ok]) ** 2
    return float(np.sqrt(np.sum(n[ok] * np.maximum(var, 0)) / np.sum(n[ok])))


def two_face_objective(params, P_up, P_lo, cen):
    a, c = S.axis_from(params, cen)
    r1, h1 = S.profile_coords(P_up, a, c); r2, h2 = S.profile_coords(P_lo, a, c)
    return float(np.sqrt(0.5 * (bin_rms(r1, h1) ** 2 + bin_rms(r2, h2) ** 2)))


def to_params(a, c, cen):
    a = a / np.linalg.norm(a); th = np.arccos(np.clip(a[2], -1, 1)); ph = np.arctan2(a[1], a[0]); e1, e2 = S.basis(a)
    return np.array([th, ph, (c - cen) @ e1, (c - cen) @ e2])


def find_axis(P_up, P_lo, n3, c_guess, cen, rng):
    """Nelder-Mead from a grid of starts: directions within a few degrees of the plate normal, centres around the guess."""
    e1, e2 = S.basis(n3)
    starts = []
    for tilt in (0.0, 4.0, 8.0):
        for az in ((0.0,) if tilt == 0 else (0.0, 90.0, 180.0, 270.0)):
            d = n3 * np.cos(np.radians(tilt)) + np.sin(np.radians(tilt)) * (np.cos(np.radians(az)) * e1 + np.sin(np.radians(az)) * e2)
            for du in (-12.0, 0.0, 12.0):
                for dv in (-12.0, 0.0, 12.0):
                    starts.append(to_params(d, c_guess + du * e1 + dv * e2, cen))
    scored = sorted(((two_face_objective(p, P_up, P_lo, cen), i) for i, p in enumerate(starts)), key=lambda t: t[0])[:8]
    out = []
    for f, i in scored:
        r = minimize(two_face_objective, starts[i], args=(P_up, P_lo, cen), method="Nelder-Mead",
                     options=dict(xatol=5e-4, fatol=1e-5, maxiter=500))
        out.append((r.fun, r.x))
    out.sort(key=lambda t: t[0])
    return out


def synthetic_base(seed, noise=0.25, n=6000):
    """Synthetic base with a ring, known axis, D-shaped footprint: to test find_axis."""
    rng = np.random.default_rng(seed)
    rho_max = 60.0
    def h_up(r): return 2.0 + 10.0 * np.exp(-((r - 50.0) / 5.0) ** 2)
    def h_lo(r): return -4.0 - 9.0 * np.exp(-((r - 50.0) / 5.0) ** 2)
    pts_u, pts_l = [], []
    while len(pts_u) < n:
        r = rho_max * np.sqrt(rng.uniform(0, 1, 4000)); ph = rng.uniform(-np.pi, np.pi, 4000)
        x, y = r * np.cos(ph), r * np.sin(ph); ok = (y < 35)          # straight break edge
        for xi, yi, ri in zip(x[ok], y[ok], r[ok]):
            pts_u.append([xi, yi, h_up(ri)]); pts_l.append([xi, yi, h_lo(ri)])
    pu = np.array(pts_u[:n]) + rng.normal(0, noise, (n, 3)) * np.array([0, 0, 1]); pl = np.array(pts_l[:n]) + rng.normal(0, noise, (n, 3)) * np.array([0, 0, 1])
    q = rng.normal(size=4); q /= np.linalg.norm(q); w, x, y, z = q
    R = np.array([[1 - 2 * (y * y + z * z), 2 * (x * y - z * w), 2 * (x * z + y * w)],
                  [2 * (x * y + z * w), 1 - 2 * (x * x + z * z), 2 * (y * z - x * w)],
                  [2 * (x * z - y * w), 2 * (y * z + x * w), 1 - 2 * (x * x + y * y)]])
    t = np.array([80.0, -30.0, 1140.0])
    return pu @ R.T + t, pl @ R.T + t, R @ np.array([0, 0, 1.0]), t


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--npz", default="sample3.npz"); ap.add_argument("--out", default="out3")
    ap.add_argument("--flip", action="store_true", help="reverse the default orientation"); ap.add_argument("--selftest", action="store_true")
    args = ap.parse_args(); os.makedirs(args.out, exist_ok=True)
    rng = np.random.default_rng(21)
    if args.selftest:
        for seed in range(3):
            Pu_, Pl_, a_t, c_t = synthetic_base(seed)
            cen_ = np.vstack([Pu_, Pl_]).mean(0)
            n3_ = a_t + np.radians(6) * np.array([0.6, -0.7, 0.3]); n3_ /= np.linalg.norm(n3_)     # deliberately 6 degrees off
            sol = find_axis(Pu_, Pl_, n3_, cen_ + np.array([7.0, -6.0, 4.0]), cen_, rng)
            a, c = S.axis_from(sol[0][1], cen_)
            ang = np.degrees(np.arccos(min(1, abs(a @ a_t))))
            at = lambda pt, ad, cc: pt + ((cc - pt) @ ad) * ad
            print("synthetic base seed %d: spread %.3f mm, direction error %.2f deg, axis position error %.2f mm" %
                  (seed, sol[0][0], ang, np.linalg.norm(at(c, a, cen_) - at(c_t, a_t, cen_))))
        return
    R = {}
    # ---------------- 1. data, ceramic core, edge mesh ----------------
    HU, sp, origin, hu = S.load_grid(args.npz)
    present = HU >= 500
    core = S.largest_component(HU >= 1200)
    vox = float(np.prod(sp))
    R["input"] = dict(voxels=int(len(hu)), voxel_size_mm=[round(float(s), 4) for s in sp],
                      core_volume_cm3=round(core.sum() * vox / 1000, 1), stored_volume_cm3=round(len(hu) * vox / 1000, 1))
    Vo, Fo, No, Wo = S.mesh_from_field(np.where(present, HU.astype(float), 0.0), 500.0, sp, origin, sigma=0.8)
    Xc = np.argwhere(core) * sp + origin; cb = Xc.mean(0)
    near = cKDTree(Xc).query(Vo)[0] < 1.6                               # drops crust flecks that do not touch the ceramic
    n3 = np.linalg.svd(Xc - cb, full_matrices=False)[2][2]
    if not args.flip: n3 = -n3          # default = the corrected orientation (stands on the foot ring); --flip gives the opposite
    faces_up = near & ((No @ n3) > 0.35); faces_lo = near & ((No @ n3) < -0.35)
    R["faces"] = dict(up_vertices=int(faces_up.sum()), lower_vertices=int(faces_lo.sum()))
    # ---------------- 2. axis ----------------
    def sample(mask, n=5000):
        idx = np.where(mask)[0]; return Vo[rng.choice(idx, min(n, len(idx)), replace=False)]
    P_up, P_lo = sample(faces_up), sample(faces_lo)
    cen = np.vstack([P_up, P_lo]).mean(0)
    # centre guess: middle of the D-shaped footprint (mean of the plate, moved towards the arc side by the plan-view centroid)
    Qc = (Xc - cb); U, Sg, Vt = np.linalg.svd(Qc, full_matrices=False)
    sols = find_axis(P_up, P_lo, n3, cb, cen, rng)
    f0, p0 = sols[0]; a, c = S.axis_from(p0, cen)
    if a @ n3 < 0: a = -a
    R["axis_candidates"] = []
    for f, p in sols[:6]:
        aa, cc = S.axis_from(p, cen); aa = aa if aa @ a > 0 else -aa
        at = lambda pt, ad, c0: pt + ((c0 - pt) @ ad) * ad
        R["axis_candidates"].append(dict(spread_mm=round(float(f), 3), angle_to_best_deg=round(float(np.degrees(np.arccos(np.clip(aa @ a, -1, 1)))), 2),
                                         centre_shift_mm=round(float(np.linalg.norm(at(cc, aa, cb) - at(c, a, cb))), 2)))
    # bootstrap by spatial blocks of the faces
    Pall = np.vstack([P_up, P_lo]); tree_u = cKDTree(P_up); tree_l = cKDTree(P_lo); bd, bs = [], []
    ctr_u = rng.choice(len(P_up), len(P_up) // 40, replace=False); ctr_l = rng.choice(len(P_lo), len(P_lo) // 40, replace=False)
    for b in range(20):
        iu = np.unique(np.concatenate(tree_u.query_ball_point(P_up[rng.choice(ctr_u, len(ctr_u))], r=5.0)))
        il = np.unique(np.concatenate(tree_l.query_ball_point(P_lo[rng.choice(ctr_l, len(ctr_l))], r=5.0)))
        r = minimize(two_face_objective, p0, args=(P_up[iu], P_lo[il], cen), method="Nelder-Mead", options=dict(xatol=1e-3, fatol=1e-5, maxiter=300))
        ab, cbx = S.axis_from(r.x, cen); ab = ab if ab @ a > 0 else -ab
        bd.append(np.degrees(np.arccos(np.clip(ab @ a, -1, 1))))
        bs.append(np.linalg.norm(np.cross(cbx - c, ab)))                 # sideways shift of the bootstrap axis from the main one
    e1, e2 = S.centred_basis(a, c, Vo[near])
    R["axis"] = dict(direction=np.round(a, 4).tolist(), point_mm=np.round(c + ((cb - c) @ a) * a, 2).tolist(),
                     spread_mm=round(float(f0), 3), angle_to_plate_normal_deg=round(float(np.degrees(np.arccos(np.clip(a @ n3, -1, 1)))), 2),
                     bootstrap_direction_deg_median=round(float(np.median(bd)), 2), bootstrap_direction_deg_p90=round(float(np.percentile(bd, 90)), 2),
                     bootstrap_centre_mm_median=round(float(np.median(bs)), 2))
    print("[2] axis %s spread %.3f mm; %.1f deg from the plate normal" % (np.round(a, 3), f0, R["axis"]["angle_to_plate_normal_deg"]))

    def frame(P):
        q = P - c; h = q @ a; x = q @ e1; y = q @ e2
        return np.hypot(x, y), h, np.arctan2(y, x)
    ro, ho, po = frame(Vo)
    # ---------------- 3. edge gradient and half-max shift (raw HU across the outer faces) ----------------
    Ff = np.where(present, HU.astype(float), 0.0)
    iu = np.where(near & (faces_up | faces_lo))[0]; iu = iu[rng.choice(len(iu), min(3000, len(iu)), replace=False)]
    us = np.arange(-6.0, 1.5, 0.1)
    Pr = Vo[iu][:, None, :] + us[None, :, None] * No[iu][:, None, :]
    pr = ndi.map_coordinates(Ff, np.moveaxis((Pr - origin) / sp, -1, 0), order=1, mode="constant", cval=0.0)
    def cross(pv, lvl):
        k = np.where((pv[:-1] >= lvl) & (pv[1:] < lvl))[0]
        return np.nan if len(k) == 0 else us[k[-1]] + (pv[k[-1]] - lvl) / (pv[k[-1]] - pv[k[-1] + 1]) * 0.1
    u1000 = np.array([cross(p_, 1000.0) for p_ in pr]); u600 = np.array([cross(p_, 600.0) for p_ in pr])
    okg = np.isfinite(u1000) & np.isfinite(u600) & (u600 > u1000)
    d_lvl = float(np.median(u600[okg] - u1000[okg])); grad = 400.0 / d_lvl
    plateau = float(np.median(HU[core])); HALF = (-1000.0 + plateau) / 2.0
    d500 = (500.0 - HALF) / grad
    R["edge"] = dict(distance_1000_to_600HU_mm=round(d_lvl, 2), gradient_HU_per_mm=int(round(grad)), plateau_HU=plateau,
                     half_max_HU=round(HALF, 1), outward_shift_mm=round(d500, 2))
    # ---------------- 4. profile h(rho) of each face ----------------
    def curve(mask, sign):
        r_, h_ = ro[mask], ho[mask]; k = np.floor(r_).astype(int)
        rr, hh, nn = [], [], []
        for kk in range(0, k.max() + 1):
            m = k == kk
            if m.sum() >= 8: rr.append(np.median(r_[m])); hh.append(np.median(h_[m])); nn.append(m.sum())
        rr, hh = np.array(rr), np.array(hh)
        hh = ndi.gaussian_filter1d(medfilt(hh, 3), 0.8, mode="nearest")
        return rr, hh + sign * d500                                # faces move outward (up for the upper face, down for the lower)
    r_up, h_up = curve(near & faces_up, +1); r_lo, h_lo = curve(near & faces_lo, -1)
    r_edge = float(np.percentile(ro[near], 99.7))
    r_max_face = float(min(r_up.max(), r_lo.max()))
    h0 = float(h_lo.min()); zdatum = h0
    up_z = h_up - zdatum; lo_z = h_lo - zdatum
    thick = np.interp(r_up, r_lo, lo_z)
    ring_az = np.sort(po[near & (ro > 45.0) & (ro < 62.0)]); gaps = np.diff(np.concatenate([ring_az, [ring_az[0] + 2 * np.pi]]))
    kmax = int(np.argmax(gaps)); az_start = float(ring_az[(kmax + 1) % len(ring_az)]); az_span = float(2 * np.pi - gaps[kmax])
    cover = float(np.degrees(az_span))
    # how well do the face points sit on the two profile curves? (2-D distance to the curves)
    cu = np.column_stack([r_up, up_z + zdatum]); cl = np.column_stack([r_lo, lo_z + zdatum])
    d_up = cKDTree(cu).query(np.column_stack([ro[near & faces_up], ho[near & faces_up]]))[0]
    d_lo = cKDTree(cl).query(np.column_stack([ro[near & faces_lo], ho[near & faces_lo]]))[0]
    dd_all = np.concatenate([d_up, d_lo])
    flat = r_lo[lo_z <= 1.5]
    R["profile"] = dict(floor_radius_covered_mm=[0.0, round(r_max_face, 1)], outer_edge_radius_mm=round(r_edge, 1),
                        outer_edge_diameter_mm=round(2 * r_edge),
                        upper_face_height_range_mm=[round(float(up_z.min()), 1), round(float(up_z.max()), 1)],
                        lower_face_height_range_mm=[round(float(lo_z.min()), 1), round(float(lo_z.max()), 1)],
                        total_height_mm=round(float(up_z.max()), 1),
                        centre_thickness_mm=round(float(np.median((up_z - np.interp(r_up, r_lo, lo_z))[r_up < 25])), 1),
                        ring_crest_upper=dict(radius_mm=round(float(r_up[np.argmax(up_z)]), 1), height_mm=round(float(up_z.max()), 1)),
                        ring_low_point=dict(radius_mm=round(float(r_lo[np.argmin(lo_z)]), 1), height_mm=round(float(lo_z.min()), 1)),
                        azimuth_covered_deg=round(cover, 0),
                        distance_of_face_points_from_profile_mm=dict(median=round(float(np.median(dd_all)), 2), p90=round(float(np.percentile(dd_all, 90)), 2)),
                        flat_foot_radius_range_mm=[round(float(flat.min()), 1), round(float(flat.max()), 1)] if len(flat) else None,
                        crest_width_mm=round(float(np.ptp(r_up[up_z >= up_z.max() - 2.0])), 1))
    np.savetxt(os.path.join(args.out, "profile_base.csv"),
               np.column_stack([r_up, up_z, np.interp(r_up, r_lo, lo_z)]), delimiter=",", header="radius_mm,upper_face_height_mm,lower_face_height_mm", comments="")
    print("[4] outer edge radius %.1f mm (diameter %.0f); upper crest %.1f mm at r=%.1f; lower low point %.1f mm at r=%.1f; azimuth covered %.0f deg" %
          (r_edge, 2 * r_edge, up_z.max(), r_up[np.argmax(up_z)], lo_z.min(), r_lo[np.argmin(lo_z)], cover))
    # ---------------- 5. sweep polygon + exports ----------------
    rmax = r_max_face
    ru = np.linspace(0, rmax, 160); hu_c = np.interp(ru, r_up, up_z); hl_c = np.interp(ru, r_lo, lo_z)
    poly = np.vstack([np.column_stack([ru, hu_c]), np.column_stack([ru, hl_c])[::-1]])
    Vl, Fl = S.revolve_mesh(poly)
    zshift = -float(lo_z.min())                                       # local files rest on z = 0
    to_scanner = lambda V_: c + V_[:, [0]] * e1 + V_[:, [1]] * e2 + (V_[:, 2] + zdatum)[:, None] * a
    S.write_stl(os.path.join(args.out, "sample3_base_scannerLPS.stl"), to_scanner(Vl), Fl)
    Vloc = Vl.copy(); Vloc[:, 2] += zshift
    S.write_stl(os.path.join(args.out, "sample3_base_local.stl"), Vloc, Fl)
    exp = dict(triangles=int(len(Fl)), watertight=S.is_watertight(Fl), mesh_volume_cm3=round(abs(S.signed_volume(Vl, Fl)) / 1000, 1))
    try:
        pieces = S.cad_pieces(np.vstack([np.column_stack([np.linspace(0, rmax, 240), np.interp(np.linspace(0, rmax, 240), r_up, up_z)]),
                                         np.column_stack([np.linspace(0, rmax, 240), np.interp(np.linspace(0, rmax, 240), r_lo, lo_z)])[::-1]]))
        import ezdxf
        doc = ezdxf.new("R2010", setup=True); doc.units = ezdxf.units.MM; msp = doc.modelspace()
        for nm, col in (("PROFILE_MEASURED", 7), ("AXIS", 1)): doc.layers.add(nm, color=col)
        for pc in pieces:
            if pc[0] == "line": msp.add_line(tuple(pc[1]), tuple(pc[2]), dxfattribs={"layer": "PROFILE_MEASURED"})
            else: msp.add_spline([tuple(p) for p in pc[1]], degree=3, dxfattribs={"layer": "PROFILE_MEASURED"})
        msp.add_line((0, float(lo_z.min()) - 8), (0, float(up_z.max()) + 8), dxfattribs={"layer": "AXIS", "linetype": "CENTER"})
        doc.saveas(os.path.join(args.out, "sample3_base_profile.dxf")); exp["dxf"] = True
    except Exception as ex:
        exp["dxf"] = "skipped: %s" % ex
    try:
        import cadquery as cq
        V3 = lambda p: cq.Vector(float(p[0]), 0.0, float(p[1]) + zshift)
        edges = [cq.Edge.makeLine(V3(pc[1]), V3(pc[2])) if pc[0] == "line" else cq.Edge.makeSpline([V3(p) for p in pc[1]]) for pc in pieces]
        sol = cq.Solid.revolve(cq.Wire.assembleEdges(edges), [], 360.0, cq.Vector(0, 0, 0), cq.Vector(0, 0, 1))
        cq.exporters.export(cq.Workplane("XY").add(sol), os.path.join(args.out, "sample3_base.step"))
        exp["step"] = dict(volume_cm3=round(sol.Volume() / 1000, 1), valid=bool(sol.isValid()))
    except Exception as ex:
        exp["step"] = "skipped: %s" % ex
    R["frame"] = dict(e1=np.round(e1, 6).tolist(), e2=np.round(e2, 6).tolist(), zdatum_mm=round(float(zdatum), 4), axis_point_mm=np.round(c, 4).tolist())   # local -> scanner: c + x e1 + y e2 + (z + zdatum) a
    R["export"] = exp; print("[5] export:", json.dumps(exp))
    # ---------------- figures ----------------
    fig, axs = plt.subplots(1, 3, figsize=(18, 6.2))
    axs[0].scatter(ro[near & faces_up], ho[near & faces_up] - zdatum, s=0.6, color="#b58c66", alpha=0.4, label="upper face")
    axs[0].scatter(ro[near & faces_lo], ho[near & faces_lo] - zdatum, s=0.6, color="#3a5f8a", alpha=0.4, label="lower face")
    axs[0].plot(r_up, up_z, "k", lw=1.4); axs[0].plot(r_lo, lo_z, "k", lw=1.4)
    axs[0].set_aspect("equal"); axs[0].set_xlabel("distance from the axis (mm)"); axs[0].set_ylabel("height (mm)"); axs[0].legend(fontsize=8, markerscale=8)
    axs[0].set_title("Both faces plotted as (radius, height) with the fitted axis", fontsize=10)
    phi_mid = float(np.angle(np.mean(np.exp(1j * po[near]))))
    rr = np.linspace(0, r_edge + 8, 260); yy = np.linspace(lo_z.min() - 6, up_z.max() + 6, 200); RRg, YYg = np.meshgrid(rr, yy)
    for ax_, dphi in zip(axs[1:], (-25, 25)):
        acc = []
        for dd in np.radians([-1.5, 0, 1.5]):
            ph = phi_mid + np.radians(dphi) + dd
            Pg = c + (RRg * np.cos(ph))[..., None] * e1 + (RRg * np.sin(ph))[..., None] * e2 + (YYg + zdatum)[..., None] * a
            acc.append(ndi.map_coordinates(Ff, np.moveaxis((Pg - origin) / sp, -1, 0), order=1, mode="constant", cval=0.0))
        img = np.mean(acc, axis=0)
        ax_.imshow(img, origin="lower", extent=[rr[0], rr[-1], yy[0], yy[-1]], cmap="gray_r", vmin=0, vmax=1500, aspect="equal")
        ax_.contour(RRg, YYg, img, levels=[500], colors="cyan", linewidths=0.7)
        ax_.plot(r_up, up_z, color="yellow", lw=1.2); ax_.plot(r_lo, lo_z, color="orange", lw=1.2)
        ax_.set_title("CT section through the axis at azimuth %+d deg (cyan = 500 HU)" % dphi, fontsize=9); ax_.set_xlabel("distance from the axis (mm)")
    fig.tight_layout(); fig.savefig(os.path.join(args.out, "fig_base_profile_sections.png"), dpi=95); plt.close(fig)
    # swept base, painter's algorithm
    def surf(rr_, zz_, phis):
        Xg = rr_[:, None] * np.cos(phis)[None, :]; Yg = rr_[:, None] * np.sin(phis)[None, :]; Zg = np.repeat(zz_[:, None], len(phis), 1)
        G = np.stack([Xg, Yg, Zg], -1); a_, b_, c_, d_ = G[:-1, :-1], G[1:, :-1], G[1:, 1:], G[:-1, 1:]
        return np.stack([a_, b_, c_, d_], 2).reshape(-1, 4, 3)
    phis = np.linspace(-np.pi, np.pi, 181); sub = slice(None, None, 2)
    polys3 = np.concatenate([surf(ru[sub], hu_c[sub], phis), surf(ru[sub], hl_c[sub], phis)])
    rim = np.stack([np.column_stack([rmax * np.cos(phis[:-1]), rmax * np.sin(phis[:-1]), np.full(len(phis) - 1, hu_c[-1])]),
                    np.column_stack([rmax * np.cos(phis[1:]), rmax * np.sin(phis[1:]), np.full(len(phis) - 1, hu_c[-1])]),
                    np.column_stack([rmax * np.cos(phis[1:]), rmax * np.sin(phis[1:]), np.full(len(phis) - 1, hl_c[-1])]),
                    np.column_stack([rmax * np.cos(phis[:-1]), rmax * np.sin(phis[:-1]), np.full(len(phis) - 1, hl_c[-1])])], axis=1)
    polys3 = np.concatenate([polys3, rim])
    cph = np.arctan2(polys3[:, :, 1].mean(1), polys3[:, :, 0].mean(1)); on = ((cph - az_start) % (2 * np.pi)) <= az_span
    colr = np.where(on[:, None], np.array([0.79, 0.63, 0.48]), np.array([0.62, 0.71, 0.80]))
    fig, axs = plt.subplots(1, 3, figsize=(18, 6.2))
    for ax_, (el, azd, ttl) in zip(axs[:2], ((30, -35, "Swept base from the top"), (-30, -35, "Swept base from below"))):
        t_, p_ = np.radians(el), np.radians(azd); cam = np.array([np.cos(t_) * np.cos(p_), np.cos(t_) * np.sin(p_), np.sin(t_)])
        right = np.cross([0, 0, 1.0], cam); right /= np.linalg.norm(right); upv = np.cross(cam, right)
        Q = polys3.reshape(-1, 3); sx, sy, dz = (Q @ right).reshape(-1, 4), (Q @ upv).reshape(-1, 4), (Q @ cam).reshape(-1, 4)
        order = np.argsort(dz.mean(1)); nn = np.cross(polys3[:, 1] - polys3[:, 0], polys3[:, 2] - polys3[:, 0]); nn /= np.linalg.norm(nn, axis=1, keepdims=True) + 1e-12
        shade = 0.5 + 0.5 * np.abs(nn @ (0.6 * cam + 0.4 * np.array([0.3, -0.4, 0.85])))
        ax_.add_collection(PolyCollection(np.stack([sx, sy], -1)[order], facecolors=np.clip(colr * shade[:, None], 0, 1)[order], edgecolors="none"))
        ax_.set_xlim(sx.min() - 3, sx.max() + 3); ax_.set_ylim(sy.min() - 3, sy.max() + 3); ax_.set_aspect("equal"); ax_.axis("off")
        ax_.set_title(ttl + "\n(clay = scanned arc, blue = repeated by symmetry)", fontsize=10)
    axs[2].fill(np.concatenate([ru, ru[::-1]]), np.concatenate([hu_c, hl_c[::-1]]), color="#d9bfa2"); axs[2].fill(-np.concatenate([ru, ru[::-1]]), np.concatenate([hu_c, hl_c[::-1]]), color="#d9bfa2", alpha=0.6)
    for sg in (1, -1): axs[2].plot(sg * ru, hu_c, "k", lw=1.3); axs[2].plot(sg * ru, hl_c, "k", lw=1.3)
    axs[2].axvline(0, color="gray", ls="-.", lw=0.8); axs[2].set_aspect("equal"); axs[2].set_xlabel("mm"); axs[2].set_title("Cross-section through the axis", fontsize=10)
    fig.suptitle("Base reconstructed from the scanned base sherd (label 62, Table 2 'Sample 3')", fontsize=12, fontweight="bold")
    fig.tight_layout(rect=(0, 0, 1, 0.94)); fig.savefig(os.path.join(args.out, "fig_base_sweep.png"), dpi=100); plt.close(fig)
    with open(os.path.join(args.out, "results.json"), "w") as fh: json.dump(R, fh, indent=2)
    print(json.dumps(R, indent=2))


if __name__ == "__main__":
    main()
