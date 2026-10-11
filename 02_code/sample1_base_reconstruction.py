#!/usr/bin/env python3
"""
Reconstruction of the wide dish on a pedestal foot that is the brown piece inside scan label 2 ("Sample 1" in Table 2).

Input  : sample1.npz (x, y, z in mm; hu int16; every voxel of label 2 at >= 500 HU)
Stages : 1 isolate the piece (remove the touching ceramic pieces and the debris), 2 find its round foot -> base plane,
         axis and centre, 3 profile (radius, height) of both faces and the sphere of the inner face, 4 wall pulled up along
         that sphere to chosen heights z, 5 exports and figures (real piece shown inside the completed vessel).
Needs sample2_reconstruction.py, sample4_wall_reconstruction.py and sample3_vessel_spherical_wall.py in the same folder.
Requirements: numpy, scipy, scikit-image, matplotlib (+ ezdxf, cadquery, trimesh optional). Fixed seeds: repeatable.
Usage  : python sample1_base_reconstruction.py --npz sample1.npz --out out1 --z 70 --z 90 --z 110
"""
import argparse, json, os
import numpy as np
from scipy import ndimage as ndi
from scipy.spatial import ConvexHull, cKDTree
from scipy.optimize import least_squares
from scipy.signal import medfilt
from skimage.morphology import ball
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.collections import PolyCollection
import sample2_reconstruction as S
import sample4_wall_reconstruction as W
from sample3_vessel_spherical_wall import fit_circle, arc_points


def isolate_piece(HU, sp):
    """The brown piece = everything stored that is not part of the (larger) ceramic group, opened to drop thin bridges and debris."""
    present, body, T, p_crust, p_body = W.smoothed_body(HU, sp)
    rest = present & ~ndi.binary_dilation(body, structure=ball(2))
    opened = ndi.binary_opening(rest, structure=ball(3))
    lab, _ = ndi.label(opened, structure=np.ones((3, 3, 3))); core = lab == (np.argmax(np.bincount(lab.ravel())[1:]) + 1)
    rec = core.copy()
    for _ in range(3): rec = ndi.binary_dilation(rec, structure=ball(1)) & rest        # give the edge back
    clean = ndi.binary_opening(rec, structure=ball(5))                                  # lumps under about 10 mm
    lab2, _ = ndi.label(clean, structure=np.ones((3, 3, 3))); clean = ndi.binary_fill_holes(lab2 == (np.argmax(np.bincount(lab2.ravel())[1:]) + 1))
    return present, body, rec, clean, T


def base_plane(V, N, mask_vox, sp, origin):
    """Round foot: the flat cap of the piece. Both faces are tried; the one with the larger flat cap is the base."""
    X = np.argwhere(mask_vox) * sp + origin; c0 = X.mean(0)
    a0 = np.linalg.svd(X - c0, full_matrices=False)[2][2]
    best = None
    for sgn in (1.0, -1.0):
        ax = sgn * a0; h = (V - c0) @ ax
        cap = (h < np.percentile(h, 0.5) + 6.0) & (N @ ax < -0.5)
        if best is None or cap.sum() > best[0]: best = (int(cap.sum()), ax, cap)
    _, ax, cap = best; P = V[cap]
    for _ in range(4):
        mu = P.mean(0); nrm = np.linalg.svd(P - mu, full_matrices=False)[2][2]
        if nrm @ ax > 0: nrm = -nrm
        P = V[cap][np.abs((V[cap] - mu) @ nrm) < 1.5]
    mu = P.mean(0); nrm = np.linalg.svd(P - mu, full_matrices=False)[2][2]
    if nrm @ ax > 0: nrm = -nrm
    up = -nrm; e1, e2 = S.basis(up)
    uv = np.column_stack([(P - mu) @ e1, (P - mu) @ e2]); bp = uv[ConvexHull(uv).vertices]
    s = least_squares(lambda p: np.hypot(bp[:, 0] - p[0], bp[:, 1] - p[1]) - p[2], [bp[:, 0].mean(), bp[:, 1].mean(), 40.0])
    cen = mu + s.x[0] * e1 + s.x[1] * e2
    return up, mu, cen, float(s.x[2]), float(np.sqrt(np.mean(s.fun ** 2))), len(P), (float(np.ptp(uv[:, 0])), float(np.ptp(uv[:, 1])))


def profile_faces(V, N, cen, mu, up, e1, e2, rmax=100.0):
    q = V - cen; rho = np.hypot(q @ e1, q @ e2); h = (V - mu) @ up
    res = {}
    for name, m in (("up", N @ up > 0.35), ("lo", N @ up < -0.35)):
        r_, h_ = rho[m], h[m]; k = np.floor(r_).astype(int); rr, hh = [], []
        for kk in range(0, int(min(k.max(), rmax)) + 1):
            s = k == kk
            if s.sum() >= 8: rr.append(np.median(r_[s])); hh.append(np.median(h_[s]))
        rr, hh = np.array(rr), np.array(hh)
        res[name] = (rr, ndi.gaussian_filter1d(medfilt(hh, 3), 0.8, mode="nearest"))
    return res, rho, h


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--npz", default="sample1.npz"); ap.add_argument("--out", default="out1"); ap.add_argument("--z", type=float, action="append", help="wall height above the base plane (mm), repeatable"); ap.add_argument("--join-radius", type=float, default=50.0); ap.add_argument("--prev-height", type=float, default=110.0, help="wall height already shown in blue; anything above it is coloured green"); ap.add_argument("--thickness", type=float, default=None, help="wall thickness in mm (default: from the measured outer skirt)")
    args = ap.parse_args(); os.makedirs(args.out, exist_ok=True)
    R = {}
    HU, sp, origin, hu = S.load_grid(args.npz); vox = float(np.prod(sp))
    present, body, rec, clean, T = isolate_piece(HU, sp)
    R["isolation"] = dict(ceramic_group_volume_cm3=round(body.sum() * vox / 1000, 1), piece_isolated_cm3=round(rec.sum() * vox / 1000, 1),
                          piece_cleaned_cm3=round(clean.sum() * vox / 1000, 1), split_threshold_HU=T,
                          piece_mean_HU=round(float(HU[rec].astype(float).mean())), ceramic_group_mean_HU=round(float(HU[body].astype(float).mean())))
    print("[1] piece isolated %.0f cm3, cleaned %.0f cm3 (mean HU %.0f vs ceramic group %.0f)" % (R["isolation"]["piece_isolated_cm3"], R["isolation"]["piece_cleaned_cm3"], R["isolation"]["piece_mean_HU"], R["isolation"]["ceramic_group_mean_HU"]))
    V, F, N, Wt = S.mesh_from_field(clean.astype(float), 0.5, sp, origin, sigma=1.0)
    up, mu, cen, Rfoot, rms_foot, nfoot, ext_foot = base_plane(V, N, clean, sp, origin)
    e1, e2 = S.basis(up)
    R["foot"] = dict(circle_radius_mm=round(Rfoot, 1), circle_diameter_mm=round(2 * Rfoot, 1), circle_rms_mm=round(rms_foot, 2), cap_vertices=nfoot,
                     cap_extents_mm=[round(ext_foot[0]), round(ext_foot[1])], axis=np.round(up, 4).tolist(), centre_mm=np.round(cen, 2).tolist())
    print("[2] foot: circle radius %.1f mm (RMS %.2f), axis %s" % (Rfoot, rms_foot, np.round(up, 3)))
    # cut debris below the base plane
    idx = np.argwhere(clean); Xv = idx * sp + origin
    below = (Xv - mu) @ up < -1.0
    c2 = clean.copy(); c2[idx[below, 0], idx[below, 1], idx[below, 2]] = False
    lab, _ = ndi.label(c2, structure=np.ones((3, 3, 3))); c2 = lab == (np.argmax(np.bincount(lab.ravel())[1:]) + 1)
    V, F, N, Wt = S.mesh_from_field(c2.astype(float), 0.5, sp, origin, sigma=1.0)
    prof, rho, h = profile_faces(V, N, cen, mu, up, e1, e2)
    r_up, h_up = prof["up"]; r_lo, h_lo = prof["lo"]
    # debris cut: material sitting more than 4 mm above the (robust, azimuth-median) inner-face profile is lumps of crust, not the piece
    idx = np.argwhere(c2); Xv = idx * sp + origin; q_ = Xv - cen; rho_v = np.hypot(q_ @ e1, q_ @ e2); h_v = (Xv - mu) @ up
    debris = (h_v > np.interp(rho_v, r_up, h_up, right=h_up[-1]) + 4.0) & (rho_v <= 85.0)
    c3 = c2.copy(); c3[idx[debris, 0], idx[debris, 1], idx[debris, 2]] = False
    lab, _ = ndi.label(c3, structure=np.ones((3, 3, 3))); c3 = lab == (np.argmax(np.bincount(lab.ravel())[1:]) + 1)
    R["isolation"]["lumps_above_inner_face_removed_cm3"] = round((c2.sum() - c3.sum()) * vox / 1000, 1); R["isolation"]["piece_final_cm3"] = round(c3.sum() * vox / 1000, 1)
    print("    lumps above the inner face removed: %.1f cm3; final piece %.0f cm3" % (R["isolation"]["lumps_above_inner_face_removed_cm3"], R["isolation"]["piece_final_cm3"]))
    V, F, N, Wt = S.mesh_from_field(c3.astype(float), 0.5, sp, origin, sigma=1.0)
    rho = np.hypot((V - cen) @ e1, (V - cen) @ e2); h = (V - mu) @ up
    # sphere of the inner face: circle with its centre ON the axis (a surface of revolution about the foot axis); the free-centre fit is kept for comparison
    m = (r_up >= 5) & (r_up <= 55)
    (xc_f, yc_f, Ri_f), rms_f = fit_circle(r_up[m], h_up[m])
    on = least_squares(lambda p: np.sqrt(r_up[m] ** 2 + (h_up[m] - p[0]) ** 2) - p[1], [150.0, 150.0]); yc, Ri = float(on.x[0]), float(on.x[1]); xc = 0.0
    rms_i = float(np.sqrt(np.mean(on.fun ** 2)))
    R["interior_sphere"] = dict(radius_mm=round(Ri, 1), centre_height_mm=round(yc, 1), fit_rms_mm=round(rms_i, 3),
                                free_centre_fit=dict(radius_mm=round(float(Ri_f), 1), centre_radius_mm=round(float(xc_f), 1), fit_rms_mm=round(float(rms_f), 3)))
    print("[3] interior sphere (centre on the axis): radius %.1f mm, centre %.1f mm above the base plane, fit RMS %.3f mm | free-centre circle: R %.0f mm, RMS %.3f" % (Ri, yc, rms_i, Ri_f, rms_f))
    cu = cKDTree(np.column_stack([r_up, h_up])); cl = cKDTree(np.column_stack([r_lo, h_lo]))
    mu_, ml_ = N @ up > 0.35, N @ up < -0.35
    d_all = np.concatenate([cu.query(np.column_stack([rho[mu_], h[mu_]]))[0], cl.query(np.column_stack([rho[ml_], h[ml_]]))[0]])
    R["profile_fit"] = dict(median_distance_mm=round(float(np.median(d_all)), 2), p90_mm=round(float(np.percentile(d_all, 90)), 2),
                            note="face points (including remaining lumps and the non-round outline) against the two profile curves")
    # exterior sphere through the measured outer skirt, concentric with the interior one
    r_x = 70.0; h_x = float(np.interp(r_x, r_lo, h_lo)); Re = float(np.hypot(r_x - xc, h_x - yc)); yc_e = yc
    if args.thickness is not None: Re = Ri + args.thickness; yc_e = h_x + float(np.sqrt(Re ** 2 - (r_x - xc) ** 2))       # exterior sphere through the measured skirt point
    t = Re - Ri
    R["wall"] = dict(exterior_sphere_radius_mm=round(Re, 1), wall_thickness_mm=round(t, 1), exterior_join=dict(radius_mm=r_x, height_mm=round(h_x, 1)),
                     interior_gap_below_sphere_at_r70_mm=round(float((yc - np.sqrt(Ri ** 2 - 70.0 ** 2)) - np.interp(70.0, r_up, h_up)), 1))
    print("    exterior sphere %.1f mm -> wall %.1f mm; measured interior at r=70 is %.1f mm below the sphere" % (Re, t, R["wall"]["interior_gap_below_sphere_at_r70_mm"]))
    np.savetxt(os.path.join(args.out, "profile_base.csv"), np.column_stack([r_up, h_up, np.interp(r_up, r_lo, h_lo)]), delimiter=",", header="radius_mm,upper_face_height_mm,lower_face_height_mm", comments="")
    to_scanner = lambda X_: cen + X_[:, [0]] * e1 + X_[:, [1]] * e2 + X_[:, [2]] * up
    S.write_stl(os.path.join(args.out, "sample1_piece_cleaned_scannerLPS.stl"), V, F)
    Xloc = np.column_stack([(V - cen) @ e1, (V - cen) @ e2, (V - mu) @ up]); S.write_stl(os.path.join(args.out, "sample1_piece_cleaned_on_its_base_local.stl"), Xloc, F)
    # ---------------- 4/5. wall pulled up to each height z, exports, figures ----------------
    zs = args.z or [70.0, 90.0, 110.0]
    j0 = int(np.argmin(np.abs(r_up - args.join_radius))); r_j = float(r_up[j0])
    lo_keep = r_lo <= r_x
    ext_meas = np.column_stack([r_lo[lo_keep], h_lo[lo_keep]]); ext_meas[-1] = [r_x, h_x]
    tree_piece = cKDTree(V); Tp = V[F]; ns_ = np.cross(Tp[:, 1] - Tp[:, 0], Tp[:, 2] - Tp[:, 0]); ns_ /= np.linalg.norm(ns_, axis=1, keepdims=True) + 1e-12
    fig, axs = plt.subplots(len(zs), 3, figsize=(16, 5.2 * len(zs))); axs = np.atleast_2d(axs); info = {}
    def paint(ax_, T3, nrm, cols, cam, upref, ttl, ctr):
        cam = cam / np.linalg.norm(cam); right = np.cross(upref, cam); right /= np.linalg.norm(right); upv = np.cross(cam, right)
        Q = T3 - ctr; sx, sy, dz = Q @ right, Q @ upv, Q @ cam; order = np.argsort(dz.mean(1)); lv = 0.6 * cam + 0.4 * np.array([0.3, -0.4, 0.85]); lv /= np.linalg.norm(lv)
        shade = 0.4 + 0.6 * np.clip(np.abs(nrm @ lv), 0, 1)
        ax_.add_collection(PolyCollection(np.stack([sx, sy], -1)[order], facecolors=np.clip(cols * shade[:, None], 0, 1)[order], edgecolors="none"))
        ax_.set_xlim(sx.min() - 4, sx.max() + 4); ax_.set_ylim(sy.min() - 4, sy.max() + 4); ax_.set_aspect("equal"); ax_.axis("off"); ax_.set_title(ttl, fontsize=9.5)
    for k, z in enumerate(zs):
        int_meas = np.column_stack([r_up[:j0 + 1], h_up[:j0 + 1]]); int_arc = arc_points(xc, yc, Ri, r_j, z); ext_arc = arc_points(xc, yc_e, Re, r_x, z)
        poly = np.vstack([int_meas, int_arc[1:], ext_arc[::-1], ext_meas[::-1][1:]])
        Vl, Fl = S.revolve_mesh(poly, nseg=180); tag = "z%03d" % int(z)
        S.write_stl(os.path.join(args.out, "sample1_vessel_%s_local.stl" % tag), Vl, Fl); S.write_stl(os.path.join(args.out, "sample1_vessel_%s_scannerLPS.stl" % tag), to_scanner(Vl), Fl)
        x_, z_ = poly[:, 0], poly[:, 1]; x2, z2 = np.roll(x_, -1), np.roll(z_, -1); pv = abs(np.pi / 3.0 * ((x_ ** 2 + x_ * x2 + x2 ** 2) * (z2 - z_)).sum()) / 1000
        ent = dict(rim_inner_diameter_mm=round(2 * float(int_arc[-1, 0])), rim_outer_diameter_mm=round(2 * float(ext_arc[-1, 0])), total_height_mm=z,
                   watertight=S.is_watertight(Fl), mesh_volume_cm3=round(abs(S.signed_volume(Vl, Fl)) / 1000), polygon_volume_cm3=round(float(pv)))
        try:
            import ezdxf, cadquery as cq
            def dense(P_, n_):
                s_ = np.concatenate([[0], np.cumsum(np.hypot(*np.diff(P_, axis=0).T))]); qq = np.linspace(0, s_[-1], n_)
                return np.column_stack([np.interp(qq, s_, P_[:, 0]), np.interp(qq, s_, P_[:, 1])])
            up_c = np.vstack([int_meas, int_arc[1:]]); lo_c = np.vstack([ext_meas, ext_arc[1:]])
            pieces = S.cad_pieces(np.vstack([dense(up_c, 700), dense(lo_c, 700)[::-1]]), tol=0.04, n_fit=110)
            doc = ezdxf.new("R2010", setup=True); doc.units = ezdxf.units.MM; msp = doc.modelspace()
            for nm, col in (("PROFILE", 7), ("AXIS", 1)): doc.layers.add(nm, color=col)
            for pc in pieces:
                if pc[0] == "line": msp.add_line(tuple(pc[1]), tuple(pc[2]), dxfattribs={"layer": "PROFILE"})
                else: msp.add_spline([tuple(p) for p in pc[1]], degree=3, dxfattribs={"layer": "PROFILE"})
            msp.add_line((0, -8), (0, z + 8), dxfattribs={"layer": "AXIS", "linetype": "CENTER"}); doc.saveas(os.path.join(args.out, "sample1_vessel_%s_profile.dxf" % tag))
            V3 = lambda p: cq.Vector(float(p[0]), 0.0, float(p[1]))
            edges = [cq.Edge.makeLine(V3(pc[1]), V3(pc[2])) if pc[0] == "line" else cq.Edge.makeSpline([V3(p) for p in pc[1]]) for pc in pieces]
            sol = cq.Solid.revolve(cq.Wire.assembleEdges(edges), [], 360.0, cq.Vector(0, 0, 0), cq.Vector(0, 0, 1))
            cq.exporters.export(cq.Workplane("XY").add(sol), os.path.join(args.out, "sample1_vessel_%s.step" % tag))
            ent["step"] = dict(volume_cm3=round(sol.Volume() / 1000), valid=bool(sol.isValid()), volume_vs_polygon_percent=round(100 * (sol.Volume() / 1000 - pv) / pv, 1))
        except Exception as ex:
            ent["cad"] = "skipped: %s" % ex
        # scan-versus-vessel distance for the base region (real piece vertices against the vessel surface)
        try:
            import trimesh
            mv = trimesh.Trimesh(vertices=to_scanner(Vl), faces=Fl, process=False)
            sel = np.where(rho < 70.0)[0]; sel = sel[np.random.default_rng(0).choice(len(sel), min(6000, len(sel)), replace=False)]
            dist = trimesh.proximity.closest_point(mv, V[sel])[1]
            ent["real_piece_vs_vessel_mm"] = dict(median=round(float(np.median(dist)), 2), within_1mm_percent=round(100 * float((dist < 1).mean())), within_2mm_percent=round(100 * float((dist < 2).mean())),
                                                  within_3mm_percent=round(100 * float((dist < 3).mean())))
        except Exception as ex:
            ent["compare"] = "skipped: %s" % ex
        info[tag] = ent
        print("z = %3.0f mm: rim %d mm outside / %d inside; %s; STEP %s; real piece vs vessel %s" % (z, ent["rim_outer_diameter_mm"], ent["rim_inner_diameter_mm"], "watertight" if ent["watertight"] else "NOT watertight", json.dumps(ent.get("step")), json.dumps(ent.get("real_piece_vs_vessel_mm"))))
        Vs = to_scanner(Vl); Tv = Vs[Fl]; cenv = Tv.mean(1); hide = tree_piece.query(cenv)[0] < 4.0
        qv = Vl[Fl].mean(1); phi_v = np.arctan2(qv[:, 1], qv[:, 0]); cut = np.abs(phi_v) < np.radians(65)
        nv = np.cross(Tv[:, 1] - Tv[:, 0], Tv[:, 2] - Tv[:, 0]); nv /= np.linalg.norm(nv, axis=1, keepdims=True) + 1e-12
        col_v = np.tile(np.array([0.62, 0.71, 0.80]), (len(Tv), 1)); col_s = np.tile(np.array([0.83, 0.68, 0.52]), (len(Tp), 1)); ctr = Vs.mean(0)
        col_v[qv[:, 2] > args.prev_height] = np.array([0.45, 0.74, 0.55])          # new projection above the height we already had (green)
        for jv, (ttl, cam, use_cut) in enumerate((("oblique from above", 0.55 * e1 + 0.45 * e2 + 0.75 * up, False), ("side view, axis vertical", e1, False), ("cutaway from above (a wedge removed)", 0.3 * e1 + 0.9 * up, True))):
            keep = ~hide & (~cut if use_cut else True)
            paint(axs[k, jv], np.concatenate([Tv[keep], Tp]), np.concatenate([nv[keep], ns_]), np.concatenate([col_v[keep], col_s]), cam, up,
                  "z = %.0f mm, %s\n(clay = REAL piece, blue = wall to %.0f mm, green = new projection)" % (z, ttl, args.prev_height), ctr)
    fig.suptitle("Sample 1 vessel: real piece on its foot, wall pulled up along the interior sphere (R = %.0f mm)" % Ri, fontsize=12, fontweight="bold")
    fig.tight_layout(rect=(0, 0, 1, 0.97)); fig.savefig(os.path.join(args.out, "fig_vessel_with_real_piece.png"), dpi=75); plt.close(fig)
    R["vessels"] = info
    # profile figure
    fig, ax = plt.subplots(1, 2, figsize=(15, 5.6))
    ax[0].scatter(rho[mu_], h[mu_], s=0.5, color="#b58c66", alpha=0.35, label="upper (interior) face"); ax[0].scatter(rho[ml_], h[ml_], s=0.5, color="#3a5f8a", alpha=0.35, label="lower (base side) face")
    ax[0].plot(r_up, h_up, "k", lw=1.4); ax[0].plot(r_lo, h_lo, "k", lw=1.4)
    th_ = np.linspace(0, 0.85, 200); ax[0].plot(Ri * np.sin(th_), yc - Ri * np.cos(th_), "r--", lw=1, label="interior sphere (centre on axis)"); ax[0].plot(Re * np.sin(th_), yc_e - Re * np.cos(th_), "m--", lw=1, label="exterior sphere")
    ax[0].set_aspect("equal"); ax[0].set_xlim(0, 110); ax[0].set_ylim(-5, 80); ax[0].legend(fontsize=8, markerscale=8); ax[0].set_xlabel("distance from the foot axis (mm)"); ax[0].set_ylabel("height above the base plane (mm)")
    ax[0].set_title("Both faces as (radius, height); black = the profile used", fontsize=10)
    phi = np.arctan2(((V - cen) @ e2), ((V - cen) @ e1)); sc = ax[1].scatter(phi[mu_] * 180 / np.pi, rho[mu_], s=1, c=h[mu_], cmap="viridis"); ax[1].set_xlabel("azimuth (deg)"); ax[1].set_ylabel("radius (mm)")
    ax[1].set_title("Upper face by azimuth (colour = height): how round is it?", fontsize=10); fig.tight_layout(); fig.savefig(os.path.join(args.out, "fig_profile.png"), dpi=80); plt.close(fig)
    with open(os.path.join(args.out, "results.json"), "w") as fh: json.dump(R, fh, indent=2)
    print(json.dumps(R, indent=1))


if __name__ == "__main__":
    main()
