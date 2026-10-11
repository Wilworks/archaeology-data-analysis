#!/usr/bin/env python3
"""
Reconstruction of a rotationally symmetric vessel WALL BAND from one CT-scanned wall sherd
(used for the piece with scan label 440, called Sample 4 in the manuscript's Table 2).

Input  : sample4.npz  (x, y, z in mm; hu int16; every voxel of the piece at >= 500 HU)
Output : results.json, profile_wall.csv, fig_wall_sweep.png, fig_wall_ct_sections.png and CAD/mesh exports
         (*.stl, *.step, *.dxf) in --out

Needs sample2_reconstruction.py in the same folder (it provides the shared helper functions).
Requirements: numpy, scipy, scikit-image, matplotlib (+ optional ezdxf, cadquery).  Fixed seeds: repeatable.

Differences from the Sample 2 script: (1) the crust/body split works on a smoothed field because the two HU peaks
overlap; (2) no finished lip exists, so the axis comes from the profile-thinness criterion alone, restricted to
axes roughly parallel to the wall (--min-angle); (3) the result is a wall BAND (only the scanned height range),
not a whole vessel; (4) which end is 'up' is not known from the data (see --flip).
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


def smoothed_body(HU, sp):
    """Body/crust split on a smoothed, edge-normalised field (the raw HU peaks overlap too much)."""
    present = HU >= 500
    F = np.where(present, HU.astype(float), 0.0)
    w = ndi.gaussian_filter(present.astype(float), 1.2)
    Fn = np.where(w > 0.5, ndi.gaussian_filter(F, 1.2) / np.maximum(w, 1e-6), 0.0)
    core = ndi.binary_erosion(present, iterations=2)
    h, e = np.histogram(Fn[core], bins=np.arange(500, 1700, 25)); c = 0.5 * (e[:-1] + e[1:])
    hs = np.convolve(h, np.ones(3) / 3, mode="same")
    p_crust = c[(c < 1100)][np.argmax(hs[c < 1100])]
    p_body = c[(c > 1150)][np.argmax(hs[c > 1150])]
    mid = (c > p_crust) & (c < p_body)
    T = float(c[mid][np.argmin(hs[mid])])
    body = ndi.binary_fill_holes(S.largest_component((Fn >= T) & present))
    return present, body, T, float(p_crust), float(p_body)


def wall_objective(params, pts, cen, deg=3):
    """RMS (mm) of the radius about a cubic function of height: rho = f(h). For a wall piece the profile is a
    single-valued curve, so this is a plain least-squares residual (the Taubin curve fit used for Sample 2 is
    unreliable for nearly straight profiles because it can pick a degenerate multiple-root curve)."""
    a, c = S.axis_from(params, cen)
    rho, h = S.profile_coords(pts, a, c)
    hh = (h - h.mean()) / 50.0
    A = np.vander(hh, deg + 1)
    coef, *_ = np.linalg.lstsq(A, rho, rcond=None)
    return float(np.sqrt(np.mean((rho - A @ coef) ** 2)))


def wall_axis_search(pts, nrm, n_up, min_angle, n_dirs=1500, keep=10):
    """Global search for the axis, restricted to directions at least min_angle degrees away from the plate normal."""
    cen = pts.mean(0); cand = []
    offs = [(0, 0)] + [(du, dv) for du in (-80, 80) for dv in (-80, 80)] + [(0, 120), (0, -120), (120, 0), (-120, 0)]
    dirs = S.fib_hemisphere(n_dirs)
    dirs = dirs[np.abs(dirs @ n_up) < np.cos(np.radians(min_angle))]
    for a in dirs:
        th = np.arccos(np.clip(a[2], -1, 1)); ph = np.arctan2(a[1], a[0])
        e1, e2 = S.basis(S.dir_from(th, ph))
        c0 = S.centre_from_normals(a, pts - cen, nrm)
        u0, v0 = c0 @ e1, c0 @ e2
        best = min(((wall_objective(np.array([th, ph, u0 + du, v0 + dv]), pts, cen), (du, dv)) for du, dv in offs), key=lambda t: t[0])
        cand.append((best[0], np.array([th, ph, u0 + best[1][0], v0 + best[1][1]])))
    cand.sort(key=lambda t: t[0]); picked = []
    for f, p in cand:
        if all(np.degrees(np.arccos(min(1, abs(S.dir_from(p[0], p[1]) @ S.dir_from(q[0], q[1]))))) > 6 for _, q in picked):
            picked.append((f, p))
        if len(picked) >= keep: break
    out = []
    for f, p in picked:
        r = minimize(wall_objective, p, args=(pts, cen), method="Nelder-Mead", options=dict(xatol=1e-3, fatol=1e-5, maxiter=400))
        out.append((r.fun, r.x))
    out.sort(key=lambda t: t[0])
    return out, cen


def synthetic_wall(seed, radius=100.0, lean_deg=8.0, height=110.0, arc_deg=46.0, noise=0.25, nnoise_deg=2.0, n=2200):
    """Synthetic wall sherd with a known axis, to test how well the axis search recovers it."""
    rng = np.random.default_rng(seed); lean = np.radians(lean_deg)
    h = rng.uniform(0, height, 6000); ph = np.radians(rng.uniform(-arc_deg / 2, arc_deg / 2, 6000))
    rho = radius + np.tan(lean) * (h - height / 2)
    P = np.column_stack([rho * np.cos(ph), rho * np.sin(ph), h])
    nrm = np.column_stack([np.cos(ph) * np.cos(lean), np.sin(ph) * np.cos(lean), -np.sin(lean) * np.ones_like(ph)])
    sel = rng.choice(len(P), n, replace=False); P, nrm = P[sel], nrm[sel]
    P = P + nrm * rng.normal(0, noise, (n, 1))
    nrm = nrm + np.radians(nnoise_deg) * rng.normal(size=(n, 3)); nrm /= np.linalg.norm(nrm, axis=1, keepdims=True)
    q = rng.normal(size=4); q /= np.linalg.norm(q); w, x, y, z = q
    R = np.array([[1 - 2 * (y * y + z * z), 2 * (x * y - z * w), 2 * (x * z + y * w)],
                  [2 * (x * y + z * w), 1 - 2 * (x * x + z * z), 2 * (y * z - x * w)],
                  [2 * (x * z - y * w), 2 * (y * z + x * w), 1 - 2 * (x * x + y * y)]])
    t = np.array([50.0, -40.0, 1380.0])
    return P @ R.T + t, nrm @ R.T, R @ np.array([0, 0, 1.0]), t, R @ np.array([1.0, 0, 0])


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--npz", default="sample4.npz"); ap.add_argument("--out", default="out4")
    ap.add_argument("--min-angle", type=float, default=55.0, help="axis must be at least this many degrees from the plate normal (wall reading)")
    ap.add_argument("--flip", action="store_true", help="flip the up direction")
    ap.add_argument("--selftest", action="store_true", help="only run the synthetic wall test of the axis search")
    args = ap.parse_args(); os.makedirs(args.out, exist_ok=True)
    R = {}
    if args.selftest:
        for seed in range(4):
            P, Nn, a_t, c_t, r_t = synthetic_wall(seed)
            n_up_t = np.cross(a_t, np.cross(a_t, r_t)); n_up_t = r_t / np.linalg.norm(r_t)      # plate normal points away from the axis
            sols, cen = wall_axis_search(P, Nn, -n_up_t, args.min_angle, n_dirs=800, keep=6)
            a, c = S.axis_from(sols[0][1], cen)
            ang = np.degrees(np.arccos(min(1, abs(a @ a_t))))
            def at(pt, ad, cc): return pt + ((cc - pt) @ ad) * ad
            dr = np.linalg.norm(at(c, a, P.mean(0)) - at(c_t, a_t, P.mean(0)))
            print("synthetic seed %d: rms %.3f mm, direction error %.1f deg, axis position error %.1f mm" % (seed, sols[0][0], ang, dr))
        return

    # ---------------- 1. grid, body/crust split, meshes ----------------
    HU, sp, origin, hu = S.load_grid(args.npz)
    present, body, T, p_crust, p_body = smoothed_body(HU, sp)
    crust = present & ~ndi.binary_dilation(body, iterations=2)
    vox = float(np.prod(sp))
    R["input"] = dict(voxels=int(len(hu)), voxel_size_mm=[round(float(s), 4) for s in sp], thresholded_volume_cm3=round(len(hu) * vox / 1000, 1))
    R["split"] = dict(threshold_HU=T, crust_peak_HU=p_crust, body_peak_HU=p_body, body_volume_cm3=round(body.sum() * vox / 1000, 1),
                      crust_class_volume_cm3=round((present & ~body).sum() * vox / 1000, 1))
    Vb, Fb, Nb, Wb = S.mesh_from_field(body.astype(float), 0.5, sp, origin, sigma=1.0)
    Vo, Fo, No, Wo = S.mesh_from_field(np.where(present, HU.astype(float), 0.0), 500.0, sp, origin, sigma=0.6)
    Xb = np.argwhere(body) * sp + origin; cb = Xb.mean(0); Xc = np.argwhere(crust) * sp + origin
    n3 = np.linalg.svd(Xb - cb, full_matrices=False)[2][2]
    n_up = n3 * np.sign((Xc.mean(0) - cb) @ n3)                       # towards the crust = the concave (inner) face
    under_b = (Nb @ n_up) < -0.35
    print("[1] body %.1f cm3, crust %.1f cm3, split at %.0f HU (peaks %.0f / %.0f)" % (R["split"]["body_volume_cm3"], R["split"]["crust_class_volume_cm3"], T, p_crust, p_body))

    # ---------------- 2. axis (wall reading) ----------------
    rng = np.random.default_rng(21)
    Pu, Nu, Wu = Vb[under_b], Nb[under_b], Wb[under_b]
    sel = rng.choice(len(Pu), size=2200, replace=False, p=Wu / Wu.sum())
    sols, cen = wall_axis_search(Pu[sel], Nu[sel], n_up, args.min_angle)
    big = Pu[rng.choice(len(Pu), size=min(6000, len(Pu)), replace=False, p=Wu / Wu.sum())]
    ref = []
    for f, p in sols[:6]:
        r = minimize(wall_objective, p, args=(big, cen), method="Nelder-Mead", options=dict(xatol=5e-4, fatol=1e-6, maxiter=600)); ref.append((r.fun, r.x))
    ref.sort(key=lambda t: t[0]); f0, p0 = ref[0]
    a, c = S.axis_from(p0, cen)
    R["axis_candidates"] = []
    for f, p in ref:
        aa, cc = S.axis_from(p, cen); aa = aa if aa @ a > 0 else -aa
        R["axis_candidates"].append(dict(rms_mm=round(float(f), 3), angle_to_best_deg=round(float(np.degrees(np.arccos(np.clip(aa @ a, -1, 1)))), 1),
                                         distance_from_sherd_mm=round(float(np.linalg.norm(np.cross(Vb.mean(0) - cc, aa))), 1)))
    # bootstrap (spatial blocks) for the spread of the axis direction and of the radius
    tree = cKDTree(Pu); centres = rng.choice(len(Pu), size=len(Pu) // 60, replace=False); bd, bdist = [], []
    for b in range(24):
        pick = rng.choice(centres, size=len(centres), replace=True)
        idx = np.unique(np.concatenate(tree.query_ball_point(Pu[pick], r=6.0)))
        if len(idx) > 4000: idx = rng.choice(idx, 4000, replace=False)
        r = minimize(wall_objective, p0, args=(Pu[idx], cen), method="Nelder-Mead", options=dict(xatol=1e-3, fatol=1e-5, maxiter=300))
        ab, cbx = S.axis_from(r.x, cen); ab = ab if ab @ a > 0 else -ab
        bd.append(np.degrees(np.arccos(np.clip(ab @ a, -1, 1)))); bdist.append(np.linalg.norm(np.cross(Vb.mean(0) - cbx, ab)))
    dist0 = float(np.linalg.norm(np.cross(Vb.mean(0) - c, a)))
    R["axis_bootstrap"] = dict(direction_spread_deg_median=round(float(np.median(bd)), 1), direction_spread_deg_p90=round(float(np.percentile(bd, 90)), 1),
                               radius_sd_mm=round(float(np.std(bdist)), 1))
    # is the crust on the axis side? (it must be if the crust is on the inside of the vessel)
    e1, e2 = S.centred_basis(a, c, Vb[under_b])
    def frame(P):
        q = P - c; h = q @ a; x = q @ e1; y = q @ e2
        return np.hypot(x, y), h, np.arctan2(y, x)
    rc_, hc_, _ = frame(Xc); rb_, hb_, pb_ = frame(Xb)
    crust_inside = bool(np.median(rc_) < np.median(rb_))
    # orientation convention: +h is the direction in which the wall leans outward (open bowl); --flip reverses it
    ro_, ho_, _ = frame(Vb[under_b])
    slope = np.polyfit(ho_, ro_, 1)[0]
    if (slope < 0) != args.flip: a = -a
    e1, e2 = S.centred_basis(a, c, Vb[under_b])
    R["axis"] = dict(direction=np.round(a, 4).tolist(), point_mm=np.round(c + ((Vb.mean(0) - c) @ a) * a, 2).tolist(),
                     profile_rms_mm=round(float(f0), 3), angle_to_plate_normal_deg=round(float(np.degrees(np.arccos(abs(a @ n_up)))), 1),
                     distance_from_sherd_centre_mm=round(dist0, 1), crust_on_axis_side=crust_inside,
                     up_convention="wall leans outward going up (assumption)" if not args.flip else "flipped by --flip")
    print("[2] axis %s rms %.3f mm; %.0f deg from plate normal; %.0f mm from sherd; crust on the axis side: %s" %
          (np.round(a, 3), f0, R["axis"]["angle_to_plate_normal_deg"], dist0, crust_inside))

    # ---------------- 3. profile: outer face, inner face, thickness ----------------
    rb, hb, pb = frame(Vb); ro, ho, po = frame(Vo)
    near = cKDTree(Vb).query(Vo)[0] < 2.5
    under_o = S.largest_vertex_component(((No @ n_up) < -0.35) & near, Fo, len(Vo))
    Ff_ = np.where(present, HU.astype(float), 0.0)
    iu = np.where(under_o)[0]; iu = iu[rng.choice(len(iu), min(3000, len(iu)), replace=False)]
    us_ = np.arange(-6.0, 1.5, 0.1)                                   # along the OUTWARD normal (negative = inside the wall)
    Pray = Vo[iu][:, None, :] + us_[None, :, None] * No[iu][:, None, :]
    prof = ndi.map_coordinates(Ff_, np.moveaxis((Pray - origin) / sp, -1, 0), order=1, mode="constant", cval=0.0)
    def cross_level(pr, lvl):                                          # outermost crossing of lvl going outward
        k = np.where((pr[:-1] >= lvl) & (pr[1:] < lvl))[0]
        if len(k) == 0: return np.nan
        k = k[-1]; return us_[k] + (pr[k] - lvl) / (pr[k] - pr[k + 1]) * 0.1
    u1000 = np.array([cross_level(p_, 1000.0) for p_ in prof]); u600 = np.array([cross_level(p_, 600.0) for p_ in prof])
    okg = np.isfinite(u1000) & np.isfinite(u600) & (u600 > u1000)
    d_lvl = float(np.median(u600[okg] - u1000[okg]))
    grad = 400.0 / d_lvl; HALF = (-1000.0 + p_body) / 2.0
    d500 = (500.0 - HALF) / grad; dT = (T - HALF) / grad
    R["edge"] = dict(distance_1000_to_600HU_mm=round(d_lvl, 2), gradient_HU_per_mm=int(round(grad)), half_max_HU=HALF,
                     outward_shift_500_mm=round(d500, 2), shift_split_level_mm=round(dT, 2),
                     note="gradient measured on the raw HU profile across the outer face")
    hmin = float(ho[under_o].min()); hmax = float(ho[under_o].max())
    step = 0.5; hg = np.arange(np.ceil(hmin), np.floor(hmax), step)
    rho_hat = np.column_stack([np.cos(pb), np.sin(pb)])
    n_rho_b = (Nb @ e1) * rho_hat[:, 0] + (Nb @ e2) * rho_hat[:, 1]
    inner_b = (n_rho_b < -0.5) & (rb > 0)
    crust_lumps = ndi.binary_opening(present & ~body, structure=ndi.generate_binary_structure(3, 1), iterations=2)
    covered = cKDTree(np.argwhere(crust_lumps) * sp + origin).query(Vb[inner_b])[0] < 2.0     # crust lies against this face
    ri_all, hi_all = rb[inner_b], hb[inner_b]
    out_r = np.full(len(hg), np.nan); in_r = np.full(len(hg), np.nan); cov = np.zeros(len(hg)); nout = np.zeros(len(hg), int)
    for k, h0 in enumerate(hg):
        m = under_o & (np.abs(ho - h0) < 0.75)
        if m.sum() >= 30: out_r[k] = np.median(ro[m]); nout[k] = m.sum()
        m2 = (np.abs(hi_all - h0) < 1.5)
        if m2.sum() >= 30: in_r[k] = np.median(ri_all[m2]); cov[k] = covered[m2].mean()
    ok = np.isfinite(out_r) & np.isfinite(in_r)
    hg, out_r, in_r, cov = hg[ok], out_r[ok], in_r[ok], cov[ok]
    sm = lambda v: ndi.gaussian_filter1d(medfilt(v, 9), 3, mode="nearest")
    out_r, in_r = sm(out_r), sm(in_r)
    out_hm = out_r + d500                                            # outer face: half-max is further out
    in_hm = in_r - np.where(cov > 0.5, 0.0, dT)                      # bare inner face: half-max is further in
    lean = np.arctan(np.gradient(out_hm, hg))
    t_norm = (out_hm - in_hm) * np.cos(lean)
    h_rel = hg - hg.min()
    # roundness: does the outer radius change with azimuth at a given height?
    az_o = np.degrees(po - np.angle(np.mean(np.exp(1j * po[under_o])))); az_o = (az_o + 180) % 360 - 180      # azimuth of every outer-mesh vertex
    az = az_o[under_o]
    dev = {}
    for b0 in range(-30, 30, 10):
        vals = []
        for h0, ref_r in zip(hg[::4], out_r[::4]):
            m = under_o & (np.abs(ho - h0) < 1.0) & (az_o >= b0) & (az_o < b0 + 10)
            if m.sum() >= 15: vals.append(np.median(ro[m]) - ref_r)
        if len(vals) >= 5: dev[b0] = float(np.median(vals))
    rr = max(dev.values()) - min(dev.values()) if len(dev) > 1 else None
    R["wall"] = dict(radius_outer_mm=dict(bottom=round(float(out_hm[0]), 1), middle=round(float(out_hm[len(out_hm) // 2]), 1), top=round(float(out_hm[-1]), 1)),
                     outer_diameter_mm=dict(bottom=round(float(2 * out_hm[0])), top=round(float(2 * out_hm[-1]))),
                     lean_from_axis_deg=round(float(np.degrees(np.median(lean))), 1),
                     band_height_mm=round(float(h_rel.max()), 1),
                     thickness_mm=dict(median=round(float(np.median(t_norm)), 1), min=round(float(t_norm.min()), 1), max=round(float(t_norm.max()), 1)),
                     surviving_azimuth_deg=round(float(np.ptp(np.percentile(az, [0.5, 99.5]))), 0),
                     roundness_range_mm=None if rr is None else round(float(rr), 1),
                     inner_face_crust_covered_fraction=round(float(np.mean(cov > 0.5)), 2))
    # ---- symmetry check straight on the CT field: does the 500 HU outer edge sit at the profile radius at other azimuths? ----
    Ff2 = np.where(present, HU.astype(float), 0.0)
    phi_mid0 = float(np.angle(np.mean(np.exp(1j * po[under_o]))))
    chk = {}
    for dg in (-20, -14, -7, 0, 7, 14, 20):
        ph_ = phi_mid0 + np.radians(dg); dif = []
        for h0, r_prof in zip(hg[::3], out_r[::3]):
            rr_ = np.arange(r_prof - 8.0, r_prof + 8.0, 0.05)
            Pq = c + (rr_ * np.cos(ph_))[:, None] * e1 + (rr_ * np.sin(ph_))[:, None] * e2 + h0 * a
            f_ = ndi.map_coordinates(Ff2, ((Pq - origin) / sp).T, order=1, mode="constant", cval=0.0)
            i_in = int(np.argmin(np.abs(rr_ - (r_prof - 3.0))))
            if f_[i_in] < 900.0: continue
            below = np.where(f_[i_in:] < 500.0)[0]
            if len(below) == 0 or below[0] == 0: continue
            k = i_in + below[0]
            dif.append(rr_[k - 1] + (f_[k - 1] - 500.0) / (f_[k - 1] - f_[k]) * 0.05 - r_prof)
        if len(dif) >= 8: chk[dg] = (float(np.median(dif)), float(np.percentile(np.abs(dif), 90)), len(dif))
    R["symmetry_check_mm"] = {str(k): dict(median_offset=round(v[0], 2), p90_abs=round(v[1], 2), n=v[2]) for k, v in chk.items()}
    R["datum_raw_height_mm"] = round(float(hg.min()), 3)
    np.savetxt(os.path.join(args.out, "profile_wall.csv"), np.column_stack([h_rel, out_hm, in_hm, t_norm]), delimiter=",",
               header="height_above_lowest_mm,outer_radius_mm,inner_radius_mm,thickness_normal_mm", comments="")
    print("[3] outer radius %.1f -> %.1f mm over %.0f mm of height; lean %.1f deg; thickness median %.1f mm; arc %.0f deg" %
          (out_hm[0], out_hm[-1], h_rel.max(), np.degrees(np.median(lean)), np.median(t_norm), R["wall"]["surviving_azimuth_deg"]))

    # ---------------- 4. the swept band + exports ----------------
    poly = np.vstack([np.column_stack([out_hm, h_rel]), np.column_stack([in_hm, h_rel])[::-1]])
    def thin(P, n):
        s_ = np.concatenate([[0], np.cumsum(np.hypot(*np.diff(P, axis=0).T))]); q = np.linspace(0, s_[-1], n)
        return np.column_stack([np.interp(q, s_, P[:, 0]), np.interp(q, s_, P[:, 1])])
    Po = thin(np.column_stack([out_hm, h_rel]), 90); Pi = thin(np.column_stack([in_hm, h_rel]), 90)
    mesh_poly = np.vstack([Po, Pi[::-1]])
    Vl, Fl = S.revolve_mesh(mesh_poly)
    h_raw0 = hg.min()                                                # raw height (along the axis from c) of the bottom of the band
    to_scanner = lambda V_: c + V_[:, [0]] * e1 + V_[:, [1]] * e2 + (V_[:, 2] + h_raw0)[:, None] * a
    S.write_stl(os.path.join(args.out, "sample4_wall_band_scannerLPS.stl"), to_scanner(Vl), Fl)
    S.write_stl(os.path.join(args.out, "sample4_wall_band_local.stl"), Vl, Fl)
    exp = dict(triangles=int(len(Fl)), watertight=S.is_watertight(Fl), mesh_volume_cm3=round(abs(S.signed_volume(Vl, Fl)) / 1000, 1))
    try:
        pieces = S.cad_pieces(np.vstack([thin(np.column_stack([out_hm, h_rel]), 240), thin(np.column_stack([in_hm, h_rel]), 240)[::-1]]))
        import ezdxf
        doc = ezdxf.new("R2010", setup=True); doc.units = ezdxf.units.MM; msp = doc.modelspace()
        for nm, col in (("PROFILE_MEASURED", 7), ("AXIS", 1)): doc.layers.add(nm, color=col)
        for pc in pieces:
            if pc[0] == "line": msp.add_line(tuple(pc[1]), tuple(pc[2]), dxfattribs={"layer": "PROFILE_MEASURED"})
            else: msp.add_spline([tuple(p) for p in pc[1]], degree=3, dxfattribs={"layer": "PROFILE_MEASURED"})
        msp.add_line((0, -10), (0, float(h_rel.max()) + 10), dxfattribs={"layer": "AXIS", "linetype": "CENTER"})
        doc.saveas(os.path.join(args.out, "sample4_wall_band_profile.dxf")); exp["dxf"] = True
    except Exception as ex:
        exp["dxf"] = "skipped: %s" % ex
    try:
        import cadquery as cq
        V3 = lambda p: cq.Vector(float(p[0]), 0.0, float(p[1]))
        edges = [cq.Edge.makeLine(V3(pc[1]), V3(pc[2])) if pc[0] == "line" else cq.Edge.makeSpline([V3(p) for p in pc[1]]) for pc in pieces]
        sol = cq.Solid.revolve(cq.Wire.assembleEdges(edges), [], 360.0, cq.Vector(0, 0, 0), cq.Vector(0, 0, 1))
        cq.exporters.export(cq.Workplane("XY").add(sol), os.path.join(args.out, "sample4_wall_band.step"))
        exp["step"] = dict(volume_cm3=round(sol.Volume() / 1000, 1), valid=bool(sol.isValid()))
    except Exception as ex:
        exp["step"] = "skipped: %s" % ex
    R["export"] = exp
    print("[4] export:", json.dumps(exp))

    # ---------------- figures ----------------
    # (a) swept band, painter's algorithm; clay = scanned arc, blue = repeated by symmetry
    fp_lo, fp_hi = np.radians(np.percentile(az, [0.5, 99.5]))
    nphi = 200; phis = np.linspace(-np.pi, np.pi, nphi + 1)
    def surf(rr, zz):
        X = rr[:, None] * np.cos(phis)[None, :]; Y = rr[:, None] * np.sin(phis)[None, :]; Z = np.repeat(zz[:, None], len(phis), 1)
        G = np.stack([X, Y, Z], axis=-1)
        a_, b_, c_, d_ = G[:-1, :-1], G[1:, :-1], G[1:, 1:], G[:-1, 1:]
        return np.stack([a_, b_, c_, d_], axis=2).reshape(-1, 4, 3)
    sub = slice(None, None, 6)
    polys3 = np.concatenate([surf(out_hm[sub], h_rel[sub]), surf(in_hm[sub], h_rel[sub])])
    cen_ph = np.arctan2(polys3[:, :, 1].mean(1), polys3[:, :, 0].mean(1))
    on_sherd = (cen_ph >= fp_lo) & (cen_ph <= fp_hi)
    colr = np.where(on_sherd[:, None], np.array([0.79, 0.63, 0.48]), np.array([0.62, 0.71, 0.80]))
    fig, axs = plt.subplots(1, 3, figsize=(18, 6.4))
    for ax_, (el, azd, ttl) in zip(axs[:2], ((25, -35, "Wall band from outside"), (60, -60, "Wall band from above"))):
        th_, ph_ = np.radians(el), np.radians(azd)
        cam = np.array([np.cos(th_) * np.cos(ph_), np.cos(th_) * np.sin(ph_), np.sin(th_)])
        right = np.cross([0, 0, 1.0], cam); right /= np.linalg.norm(right); up = np.cross(cam, right)
        Q = polys3.reshape(-1, 3); sx, sy, dz = Q @ right, Q @ up, Q @ cam
        sx, sy, dz = sx.reshape(-1, 4), sy.reshape(-1, 4), dz.reshape(-1, 4)
        order = np.argsort(dz.mean(1))
        nn = np.cross(polys3[:, 1] - polys3[:, 0], polys3[:, 2] - polys3[:, 0]); nn /= np.linalg.norm(nn, axis=1, keepdims=True) + 1e-12
        shade = 0.5 + 0.5 * np.abs(nn @ (0.6 * cam + 0.4 * np.array([0.3, -0.4, 0.85])))
        ax_.add_collection(PolyCollection(np.stack([sx, sy], -1)[order], facecolors=np.clip(colr * shade[:, None], 0, 1)[order], edgecolors="none"))
        ax_.set_xlim(sx.min() - 3, sx.max() + 3); ax_.set_ylim(sy.min() - 3, sy.max() + 3); ax_.set_aspect("equal"); ax_.axis("off"); ax_.set_title(ttl + "\n(clay = scanned arc, blue = repeated by symmetry)", fontsize=10)
    axs[2].fill(np.concatenate([out_hm, in_hm[::-1]]), np.concatenate([h_rel, h_rel[::-1]]), color="#d9bfa2")
    axs[2].plot(out_hm, h_rel, "k", lw=1.6, label="outer face"); axs[2].plot(in_hm, h_rel, color="#c04020", lw=1.4, label="inner face")
    axs[2].axvline(0, color="gray", ls="-.", lw=0.8); axs[2].set_xlim(0, out_hm.max() + 15); axs[2].set_aspect("equal")
    axs[2].set_xlabel("distance from the axis (mm)"); axs[2].set_ylabel("height above the lowest scanned point (mm)"); axs[2].legend(fontsize=8)
    axs[2].set_title("Cross-section through the axis (only the scanned height is known)", fontsize=10)
    fig.suptitle("Wall band reconstructed from the scanned wall sherd (label 440, Table 2 'Sample 4')", fontsize=12, fontweight="bold")
    fig.tight_layout(rect=(0, 0, 1, 0.94)); fig.savefig(os.path.join(args.out, "fig_wall_sweep.png"), dpi=100); plt.close(fig)
    # (b) real CT sections through the axis with the extracted profile on top
    Ff = np.where(present, HU.astype(float), 0.0)
    rr = np.linspace(out_hm.min() - 25, out_hm.max() + 15, 260); yy = np.linspace(-6, h_rel.max() + 6, 320); RRg, YYg = np.meshgrid(rr, yy)
    phi_mid = float(np.angle(np.mean(np.exp(1j * po[under_o]))))
    fig, axs = plt.subplots(1, 3, figsize=(17, 6.6))
    for ax_, dphi in zip(axs, (-14, 0, 14)):
        acc = []
        for dd in np.radians([-1.5, 0, 1.5]):
            ph = phi_mid + np.radians(dphi) + dd
            Pg = c + (RRg * np.cos(ph))[..., None] * e1 + (RRg * np.sin(ph))[..., None] * e2 + (YYg + h_raw0)[..., None] * a
            acc.append(ndi.map_coordinates(Ff, np.moveaxis((Pg - origin) / sp, -1, 0), order=1, mode="constant", cval=0.0))
        img = np.mean(acc, axis=0)
        ax_.imshow(img, origin="lower", extent=[rr[0], rr[-1], yy[0], yy[-1]], cmap="gray_r", vmin=0, vmax=1500, aspect="equal")
        ax_.contour(RRg, YYg, img, levels=[500], colors="cyan", linewidths=0.7); ax_.contour(RRg, YYg, img, levels=[T], colors="red", linewidths=0.6)
        ax_.plot(out_hm, h_rel, color="yellow", lw=1.3, label="outer face (half-max)"); ax_.plot(in_hm, h_rel, color="orange", lw=1.3, label="inner face (half-max)")
        ax_.set_title("CT section through the axis, azimuth %+d deg\n(cyan 500 HU, red %d HU)" % (dphi, int(T)), fontsize=9)
        ax_.set_xlabel("distance from the axis (mm)"); ax_.set_ylabel("height above the lowest scanned point (mm)")
    axs[0].legend(fontsize=7, loc="upper left")
    fig.tight_layout(); fig.savefig(os.path.join(args.out, "fig_wall_ct_sections.png"), dpi=95); plt.close(fig)
    with open(os.path.join(args.out, "results.json"), "w") as fh: json.dump(R, fh, indent=2)
    print(json.dumps(R, indent=2))


if __name__ == "__main__":
    main()
