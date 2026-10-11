#!/usr/bin/env python3
"""make_final_numbers.py - builds FINAL_NUMBERS.csv and FINAL_NUMBERS.md from the result files of the reconstruction scripts.
Every row: sample, quantity, value, unit, uncertainty, how obtained, status (MEASURED / ESTIMATE / ASSUMED).
Usage: python make_final_numbers.py --numbers ../03_numbers --out ../03_numbers"""
import argparse, csv, json, os

def main():
    ap = argparse.ArgumentParser(); ap.add_argument("--numbers", required=True); ap.add_argument("--out", required=True); a = ap.parse_args()
    J = lambda n: json.load(open(os.path.join(a.numbers, n)))
    s1, s1p, s2, s3, s3w, s4 = J("sample1_results.json"), J("sample1_pot_results.json"), J("sample2_results.json"), J("sample3_base_results.json"), J("sample3_vessel_walls.json"), J("sample4_wall_results.json")
    R = []
    def add(*row): R.append(row)
    v = s2["vessel"]
    add("2", "axis tilt from the scanner z axis", s2["axis"]["tilt_from_scanner_z_deg"], "deg", "+/-1-2", "profile-thinness search, confirmed by the lip plane", "MEASURED")
    add("2", "lip-plane normal versus profile axis", s2["lip"]["angle_to_profile_axis_deg"], "deg", "selection spread %.2f deg" % s2["lip"]["selection_spread_deg"], "independent lip plane fit", "MEASURED")
    add("2", "profile scatter about the axis (Taubin)", s2["axis"]["profile_rms_mm"], "mm", "noise level in synthetic tests 0.27-0.33", "implicit degree-4 curve", "MEASURED")
    add("2", "rim outer diameter", v["rim_outer_diameter_mm"], "mm", "+/-3 (wall out of round by %.1f mm)" % s2["roundness_mm"]["wall_departure_range"], "profile intersected with the lip plane", "MEASURED")
    add("2", "widest outer diameter", v["max_outer_diameter_mm"], "mm", "+/-3", "profile maximum", "MEASURED")
    add("2", "height of the widest point above the lowest surviving point", v["height_of_max_diameter_mm"], "mm", "datum is a broken tip", "profile", "MEASURED")
    add("2", "rim height above the lowest surviving point", v["rim_height_above_lowest_point_mm"], "mm", "minimum vertical extent, not a vessel height", "lip plane", "MEASURED")
    add("2", "lip radial width", v["lip_radial_width_mm"], "mm", "+/-1", "outer and inner curve at the lip plane", "ESTIMATE")
    add("2", "wall thickness (median; range)", "%.1f (%.1f-%.1f)" % (v["thickness_mm"]["median"], v["thickness_mm"]["min"], v["thickness_mm"]["max"]), "mm", "+/-0.7", "500 HU outer to 1212 HU inner, half-maximum corrected", "ESTIMATE")
    add("2", "scanned share of the circumference", v["surviving_azimuth_deg"], "deg", "", "azimuth range of the underside", "MEASURED")
    add("2", "base inside this radius", v["lowest_surviving_point_radius_mm"], "mm", "", "not measured; grey cap in figures is illustrative", "ASSUMED (illustrative cap)")
    add("2", "roundness: base / wall departure range", "%.1f / %.1f" % (s2["roundness_mm"]["base_departure_range"], s2["roundness_mm"]["wall_departure_range"]), "mm", "", "edge position at 500 HU against azimuth", "MEASURED")
    p = s3["profile"]; ax3 = s3["axis"]
    add("3", "axis: angle to the plate normal", ax3["angle_to_plate_normal_deg"], "deg", "bootstrap %.2f deg (90th pct %.2f)" % (ax3["bootstrap_direction_deg_median"], ax3["bootstrap_direction_deg_p90"]), "two-face bin criterion; synthetic recovery 0.01 deg", "MEASURED")
    add("3", "foot ring: radius of the lowest point", p["ring_low_point"]["radius_mm"], "mm", "flat bottom spans %.1f-%.1f mm" % tuple(p["flat_foot_radius_range_mm"]), "profile", "MEASURED")
    add("3", "standing circle diameter (centreline of the flat bottom)", round(2 * p["ring_low_point"]["radius_mm"]), "mm", "%d-%d across the flat bottom" % (round(2 * p["flat_foot_radius_range_mm"][0]), round(2 * p["flat_foot_radius_range_mm"][1])), "profile", "MEASURED")
    add("3", "outer edge diameter of the surviving flange", p["outer_edge_diameter_mm"], "mm", "flange thin, least certain", "profile", "MEASURED")
    add("3", "floor thickness at the centre", p["centre_thickness_mm"], "mm", "+/-0.7", "vertical, upper minus lower face", "ESTIMATE")
    add("3", "highest surviving material above the standing plane", p["total_height_mm"], "mm", "broken edge", "profile", "MEASURED")
    add("3", "scanned share of the circumference", p["azimuth_covered_deg"], "deg", "", "largest empty azimuth gap", "MEASURED")
    add("3", "interior sphere radius (circle in r,h; centre near the axis)", s3w["interior_circle"]["radius_mm"], "mm", "fit RMS %.3f mm" % s3w["interior_circle"]["fit_rms_mm"], "circle fit to the inner-face profile, r = 5-52 mm", "MEASURED")
    add("3", "exterior sphere radius / wall thickness of the completion", "%.1f / %.1f" % (s3w["exterior_sphere_radius_mm"], s3w["wall_thickness_mm"]), "mm", "thickness +/-2", "from the top of the outer flank of the ring", "ASSUMED")
    for k, e in s3w["heights"].items(): add("3", "completed vessel %s: rim outer / inner diameter" % k, "%d / %d" % (e["rim_outer_diameter_mm"], e["rim_inner_diameter_mm"]), "mm", "wall height is a choice", "sphere continuation", "ASSUMED")
    i1 = s1["isolation"]; f1 = s1["foot"]; sp1 = s1["interior_sphere"]
    add("1", "brown piece: isolated / cleaned / final volume", "%d / %d / %d" % (i1["piece_isolated_cm3"], i1["piece_cleaned_cm3"], i1["piece_final_cm3"]), "cm3", "", "isolation, opening, debris cut", "MEASURED")
    add("1", "brown piece mean HU versus the ceramic group", "%d vs %d" % (i1["piece_mean_HU"], i1["ceramic_group_mean_HU"]), "HU", "", "voxel statistics", "MEASURED")
    add("1", "foot circle diameter", f1["circle_diameter_mm"], "mm", "circle RMS %.2f mm" % f1["circle_rms_mm"], "flat cap plane and outline circle", "MEASURED")
    add("1", "interior sphere radius (centre on the axis)", sp1["radius_mm"], "mm", "free-centre circle gives %.0f mm" % sp1["free_centre_fit"]["radius_mm"], "circle fit to the inner-face profile", "ESTIMATE (125-175)")
    add("1", "wall thickness of the completion", s1["wall"]["wall_thickness_mm"], "mm", "measured wall at r = 70 mm is about 7", "concentric exterior sphere", "ASSUMED")
    add("1", "median distance of the real piece from the vessel", s1["vessels"]["z090"]["real_piece_vs_vessel_mm"]["median"], "mm", "within 1/2/3 mm: %d/%d/%d %%" % (s1["vessels"]["z090"]["real_piece_vs_vessel_mm"]["within_1mm_percent"], s1["vessels"]["z090"]["real_piece_vs_vessel_mm"]["within_2mm_percent"], s1["vessels"]["z090"]["real_piece_vs_vessel_mm"]["within_3mm_percent"]), "closest-point distance", "MEASURED")
    for tag, src in (("bowls", s1), ("pot projection", s1p)):
        for k, e in src["vessels"].items(): add("1", "%s %s: rim outer / inner diameter" % (tag, k), "%d / %d" % (e["rim_outer_diameter_mm"], e["rim_inner_diameter_mm"]), "mm", "wall height is a choice", "sphere continuation", "ASSUMED")
    w4 = s4["wall"]; ax4 = s4["axis"]
    add("4", "wall band: outer diameter at the bottom / belly / top", "%d / %d / %d" % (w4["outer_diameter_mm"]["bottom"], 2 * w4["radius_outer_mm"]["middle"], w4["outer_diameter_mm"]["top"]), "mm", "+/-2", "profile about the wall-type axis", "ESTIMATE")
    add("4", "band height", w4["band_height_mm"], "mm", "", "scanned extent along the axis", "MEASURED")
    add("4", "wall thickness (median; range)", "%.1f (%.1f-%.1f)" % (w4["thickness_mm"]["median"], w4["thickness_mm"]["min"], w4["thickness_mm"]["max"]), "mm", "+/-1", "meridian thickness, half-maximum corrected", "ESTIMATE")
    add("4", "axis profile spread (least squares, rho against h)", ax4["profile_rms_mm"], "mm", "alternatives at 0.539 (56 deg away) and 0.588 (35 deg away)", "wall-restricted global search", "ESTIMATE (weakly determined)")
    add("4", "axis direction bootstrap spread", s4["axis_bootstrap"]["direction_spread_deg_median"], "deg", "90th pct %.1f deg" % s4["axis_bootstrap"]["direction_spread_deg_p90"], "spatial-block bootstrap", "ESTIMATE")
    add("4", "surviving arc", w4["surviving_azimuth_deg"], "deg", "", "azimuth range", "MEASURED")
    add("4", "up direction", "not determined", "", "band is nearly symmetric top to bottom", "convention: wall leans outward going up", "ASSUMED")
    with open(os.path.join(a.out, "FINAL_NUMBERS.csv"), "w", newline="") as fh:
        w = csv.writer(fh); w.writerow(["sample", "quantity", "value", "unit", "uncertainty_or_note", "method", "status"]); w.writerows(R)
    with open(os.path.join(a.out, "FINAL_NUMBERS.md"), "w") as fh:
        fh.write("# Final numbers (generated by make_final_numbers.py from the result files)\n\nStatus: MEASURED = directly from the scan; ESTIMATE = derived with stated assumptions; ASSUMED = a choice, not a measurement.\n\n")
        fh.write("| Sample | Quantity | Value | Unit | Uncertainty / note | Method | Status |\n|---|---|---|---|---|---|---|\n")
        for r in R: fh.write("| " + " | ".join(str(x) for x in r) + " |\n")
    print("wrote", len(R), "rows")

if __name__ == "__main__":
    main()
