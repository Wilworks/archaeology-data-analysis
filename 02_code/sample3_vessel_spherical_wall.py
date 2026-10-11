#!/usr/bin/env python3
"""
Sample 3 vessel: measured base + a wall "pulled up" to a chosen height z along the radius of curvature of the inner surface.

Finding used: the interior surface of the base is a circular arc in the (radius, height) plane (fit RMS about 0.05 mm over
r = 5-52 mm) with its centre almost on the axis, i.e. a spherical bowl of radius R_i (about 120 mm). The wall simply
continues on that sphere. The outer surface continues on a concentric sphere of radius R_e = R_i + t, where t follows from
the measured outer flank of the foot ring (or from --thickness). Nothing above the broken edge was scanned, so z is a choice.

Inputs : profile_base.csv, results.json (both from sample3_base_reconstruction.py), sample3.npz
Outputs: per height z: local STL, scanner-frame STL, STEP, DXF; the cleaned real sherd (scanner-frame STL);
         fig_vessel_with_real_sherd.png; walls.json
Needs sample2_reconstruction.py in the same folder. Requirements: numpy, scipy, scikit-image, matplotlib (+ ezdxf, cadquery, trimesh optional).
Usage  : python sample3_vessel_spherical_wall.py --profile profile_base.csv --results results.json --npz sample3.npz --z 60 --z 80 --z 100
"""
import argparse, json, os
import numpy as np
from scipy import ndimage as ndi
from scipy.optimize import least_squares
from scipy.spatial import cKDTree
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.collections import PolyCollection
import sample2_reconstruction as S


def fit_circle(x, y):
    A = np.column_stack([2 * x, 2 * y, np.ones(len(x))]); b = x ** 2 + y ** 2
    c, *_ = np.linalg.lstsq(A, b, rcond=None)
    p0 = [c[0], c[1], np.sqrt(c[2] + c[0] ** 2 + c[1] ** 2)]
    s = least_squares(lambda p: np.hypot(x - p[0], y - p[1]) - p[2], p0)
    return s.x, float(np.sqrt(np.mean(s.fun ** 2)))


def arc_points(xc, yc, R, r_from, z_to, n=160):
    """Points on the circle (centre xc, yc; radius R) from radius r_from up to height z_to, lower half of the circle."""
    th0 = np.arcsin((r_from - xc) / R); th1 = np.arccos(np.clip((yc - z_to) / R, -1, 1))
    th = np.linspace(th0, th1, n)
    return np.column_stack([xc + R * np.sin(th), yc - R * np.cos(th)])


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--profile", default="profile_base.csv"); ap.add_argument("--results", default="results.json")
    ap.add_argument("--npz", default="sample3.npz"); ap.add_argument("--out", default="out_w")
    ap.add_argument("--z", type=float, action="append", help="wall height above the ground in mm (repeatable)")
    ap.add_argument("--join-radius", type=float, default=50.0, help="radius where the measured interior hands over to the sphere")
    ap.add_argument("--thickness", type=float, default=None, help="wall thickness in mm (default: from the measured outer flank)")
    args = ap.parse_args(); os.makedirs(args.out, exist_ok=True)
    zs = args.z or [60.0, 80.0, 100.0]
    P = np.loadtxt(args.profile, delimiter=",", skiprows=1); r, up, lo = P[:, 0], P[:, 1], P[:, 2]
    R0 = json.load(open(args.results)); a = np.array(R0["axis"]["direction"]); fr = R0["frame"]
    e1, e2, zdat, c = np.array(fr["e1"]), np.array(fr["e2"]), fr["zdatum_mm"], np.array(fr["axis_point_mm"])
    m = (r >= 5) & (r <= args.join_radius + 2)
    (xc, yc, Ri), rms_i = fit_circle(r[m], up[m])
    i_f = int(np.where((r > 55) & (lo >= 20.0))[0][0]); r_f, h_f = float(r[i_f]), float(lo[i_f])       # top of the outer flank of the foot ring
    Re = float(np.hypot(r_f - xc, h_f - yc)) if args.thickness is None else Ri + args.thickness
    t = Re - Ri
    info = dict(interior_circle=dict(radius_mm=round(Ri, 1), centre_radius_mm=round(xc, 1), centre_height_mm=round(yc, 1), fit_rms_mm=round(rms_i, 3)),
                exterior_sphere_radius_mm=round(Re, 1), wall_thickness_mm=round(t, 1), flank_top=dict(radius_mm=round(r_f, 1), height_mm=round(h_f, 1)),
                measured_interior_slope_at_r69_deg=round(float(np.degrees(np.arctan(np.gradient(up, r)[-1]))), 1),
                sphere_slope_at_r69_deg=round(float(np.degrees(np.arcsin((69.0 - xc) / Ri))), 1), heights={})
    print("interior sphere: radius %.1f mm, centre at r=%.1f, h=%.1f mm (fit RMS %.3f mm); exterior sphere %.1f mm; wall %.1f mm" % (Ri, xc, yc, rms_i, Re, t))
    # ---- the real sherd (ceramic core surface with the crust flecks left out) ----
    HU, sp, origin, hu = S.load_grid(args.npz)
    present = HU >= 500; core = S.largest_component(HU >= 1200)
    Vo, Fo, No, Wo = S.mesh_from_field(np.where(present, HU.astype(float), 0.0), 500.0, sp, origin, sigma=0.8)
    Xc = np.argwhere(core) * sp + origin; near = cKDTree(Xc).query(Vo)[0] < 1.6
    face_ok = near[Fo].all(1)
    S.write_stl(os.path.join(args.out, "sample3_scanned_sherd_cleaned_scannerLPS.stl"), Vo, Fo[face_ok])
    to_scanner = lambda V_: c + V_[:, [0]] * e1 + V_[:, [1]] * e2 + (V_[:, 2] + zdat)[:, None] * a
    q = Vo[near] - c; phi_s = np.angle(np.mean(np.exp(1j * np.arctan2(q @ e2, q @ e1))))
    tree_scan = cKDTree(Vo[near])
    # ---- one vessel per height ----
    fig, axs = plt.subplots(len(zs), 3, figsize=(16, 5.2 * len(zs))); axs = np.atleast_2d(axs)
    Tscan = Vo[Fo[face_ok]]; nscan = np.cross(Tscan[:, 1] - Tscan[:, 0], Tscan[:, 2] - Tscan[:, 0]); nscan /= np.linalg.norm(nscan, axis=1, keepdims=True) + 1e-12
    def paint(ax, T3, nrm, cols, cam, upref, ttl, ctr):
        cam = cam / np.linalg.norm(cam); right = np.cross(upref, cam); right /= np.linalg.norm(right); upv = np.cross(cam, right)
        Q = T3 - ctr; sx, sy, dz = Q @ right, Q @ upv, Q @ cam
        order = np.argsort(dz.mean(1)); lv = 0.6 * cam + 0.4 * np.array([0.3, -0.4, 0.85]); lv /= np.linalg.norm(lv)
        shade = 0.4 + 0.6 * np.clip(np.abs(nrm @ lv), 0, 1)
        ax.add_collection(PolyCollection(np.stack([sx, sy], -1)[order], facecolors=np.clip(cols * shade[:, None], 0, 1)[order], edgecolors="none"))
        ax.set_xlim(sx.min() - 4, sx.max() + 4); ax.set_ylim(sy.min() - 4, sy.max() + 4); ax.set_aspect("equal"); ax.axis("off"); ax.set_title(ttl, fontsize=9.5)
    for k, z in enumerate(zs):
        if z < 32.0: print("z = %.0f mm is below the measured crest; skipped" % z); continue
        j0 = int(np.argmin(np.abs(r - args.join_radius)))
        # interior: measured up to the join, then the sphere up to z
        int_meas = np.column_stack([r[:j0 + 1], up[:j0 + 1]])
        int_arc = arc_points(xc, yc, Ri, float(r[j0]), z)
        # exterior: measured up to the top of the flank, then the concentric sphere up to z
        ext_meas = np.column_stack([r[:i_f + 1], lo[:i_f + 1]]); ext_meas[-1] = [xc + np.sqrt(Re ** 2 - (yc - h_f) ** 2), h_f]
        ext_arc = arc_points(xc, yc, Re, float(ext_meas[-1, 0]), z)
        poly = np.vstack([int_meas, int_arc[1:], ext_arc[::-1], ext_meas[::-1][1:]])
        Vl, Fl = S.revolve_mesh(poly, nseg=180)
        tag = "z%03d" % int(z)
        S.write_stl(os.path.join(args.out, "sample3_vessel_%s_local.stl" % tag), Vl, Fl)
        S.write_stl(os.path.join(args.out, "sample3_vessel_%s_scannerLPS.stl" % tag), to_scanner(Vl), Fl)
        rim_i, rim_e = float(int_arc[-1, 0]), float(ext_arc[-1, 0])
        x_, z_ = poly[:, 0], poly[:, 1]; x2, z2 = np.roll(x_, -1), np.roll(z_, -1)
        pv = abs(np.pi / 3.0 * ((x_ ** 2 + x_ * x2 + x2 ** 2) * (z2 - z_)).sum()) / 1000
        ent = dict(rim_inner_diameter_mm=round(2 * rim_i), rim_outer_diameter_mm=round(2 * rim_e), total_height_mm=round(z, 1),
                   triangles=int(len(Fl)), watertight=S.is_watertight(Fl), mesh_volume_cm3=round(abs(S.signed_volume(Vl, Fl)) / 1000, 0),
                   polygon_volume_cm3=round(float(pv), 0))
        try:
            import ezdxf, cadquery as cq
            def dense(P_, n): 
                s_ = np.concatenate([[0], np.cumsum(np.hypot(*np.diff(P_, axis=0).T))]); qq = np.linspace(0, s_[-1], n)
                return np.column_stack([np.interp(qq, s_, P_[:, 0]), np.interp(qq, s_, P_[:, 1])])
            up_c = np.vstack([int_meas, int_arc[1:]]); lo_c = np.vstack([ext_meas, ext_arc[1:]])
            pieces = S.cad_pieces(np.vstack([dense(up_c, 700), dense(lo_c, 700)[::-1]]), tol=0.04, n_fit=110)
            doc = ezdxf.new("R2010", setup=True); doc.units = ezdxf.units.MM; msp = doc.modelspace()
            for nm, col in (("PROFILE", 7), ("AXIS", 1)): doc.layers.add(nm, color=col)
            for pc in pieces:
                if pc[0] == "line": msp.add_line(tuple(pc[1]), tuple(pc[2]), dxfattribs={"layer": "PROFILE"})
                else: msp.add_spline([tuple(p) for p in pc[1]], degree=3, dxfattribs={"layer": "PROFILE"})
            msp.add_line((0, -8), (0, z + 8), dxfattribs={"layer": "AXIS", "linetype": "CENTER"})
            doc.saveas(os.path.join(args.out, "sample3_vessel_%s_profile.dxf" % tag)); ent["dxf"] = True
            V3 = lambda p: cq.Vector(float(p[0]), 0.0, float(p[1]))
            edges = [cq.Edge.makeLine(V3(pc[1]), V3(pc[2])) if pc[0] == "line" else cq.Edge.makeSpline([V3(p) for p in pc[1]]) for pc in pieces]
            sol = cq.Solid.revolve(cq.Wire.assembleEdges(edges), [], 360.0, cq.Vector(0, 0, 0), cq.Vector(0, 0, 1))
            cq.exporters.export(cq.Workplane("XY").add(sol), os.path.join(args.out, "sample3_vessel_%s.step" % tag))
            ent["step"] = dict(volume_cm3=round(sol.Volume() / 1000, 0), valid=bool(sol.isValid()), volume_vs_polygon_percent=round(100 * (sol.Volume() / 1000 - pv) / pv, 1))
        except Exception as ex:
            ent["cad"] = "skipped: %s" % ex
        info["heights"][tag] = ent
        print("z = %3.0f mm: rim diameter %d mm outside / %d mm inside; %s; STEP %s" % (z, ent["rim_outer_diameter_mm"], ent["rim_inner_diameter_mm"],
              "watertight" if ent["watertight"] else "NOT watertight", json.dumps(ent.get("step"))))
        # ---- figure row: completed vessel with the real sherd in place ----
        Vs = to_scanner(Vl); Tv = Vs[Fl]; cen = Tv.mean(1)
        hide = tree_scan.query(cen)[0] < 3.0                                # where the real sherd exists, show the real sherd instead
        qv = Vl[Fl].mean(1); phi_v = np.arctan2(qv[:, 1], qv[:, 0])
        cut = np.abs((phi_v - phi_s + np.pi) % (2 * np.pi) - np.pi) > np.radians(130)     # wedge on the far side, for the cutaway
        nv = np.cross(Tv[:, 1] - Tv[:, 0], Tv[:, 2] - Tv[:, 0]); nv /= np.linalg.norm(nv, axis=1, keepdims=True) + 1e-12
        col_v = np.tile(np.array([0.62, 0.71, 0.80]), (len(Tv), 1)); col_s = np.tile(np.array([0.83, 0.68, 0.52]), (len(Tscan), 1))
        ctr = Vs.mean(0)
        for jv, (ttl, cam, upr, use_cut) in enumerate((("oblique from above", 0.55 * e1 + 0.4 * e2 + 0.75 * a, a, False),
                                                       ("side view, axis vertical", -np.cos(phi_s) * e1 - np.sin(phi_s) * e2, a, False),
                                                       ("cutaway from above (far side removed)", 0.35 * (np.cos(phi_s) * e1 + np.sin(phi_s) * e2) * -1 + 0.9 * a, a, True))):
            keep = ~hide & (~cut if use_cut else True)
            T3 = np.concatenate([Tv[keep], Tscan]); nn = np.concatenate([nv[keep], nscan]); cc = np.concatenate([col_v[keep], col_s])
            paint(axs[k, jv], T3, nn, cc, cam, upr, "z = %.0f mm, %s\n(clay = the REAL scanned sherd, blue = assumed completion)" % (z, ttl), ctr)
    fig.suptitle("Sample 3 vessel: real sherd in place, wall pulled up along the interior sphere (R = %.0f mm)" % Ri, fontsize=12, fontweight="bold")
    fig.tight_layout(rect=(0, 0, 1, 0.97)); fig.savefig(os.path.join(args.out, "fig_vessel_with_real_sherd.png"), dpi=75); plt.close(fig)
    with open(os.path.join(args.out, "walls.json"), "w") as fh: json.dump(info, fh, indent=2)
    print(json.dumps(info, indent=2))


if __name__ == "__main__":
    main()
