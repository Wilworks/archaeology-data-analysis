#!/usr/bin/env python3
"""
run_all.py - rebuilds every result, figure and 3D file of the package from the four point files.

    python run_all.py --data ../01_inputs/point_files --work ../work            # everything (about 30-40 minutes)
    python run_all.py --data ../01_inputs/point_files --work ../work --only s2 s4
    python run_all.py --list

Steps (in this order; later steps need the earlier ones)
  s2      Sample 2 rim-to-shoulder bowl (axis, lip check, profile, sweep, STL/STEP/DXF)
  s4      Sample 4 wall band
  s3      Sample 3 ring-footed base (accepted orientation)
  s3flip  Sample 3 base in the rejected orientation (only needed for the orientation figure)
  s3v     Sample 3 vessels with the wall pulled up along the interior sphere (z = 60, 80, 100 mm)
  s1      Sample 1 brown piece: isolation, foot, profile, bowls z = 70, 90, 110 mm
  s1pot   Sample 1 pot projection z = 170, 210, 250 mm (wall above 110 mm drawn green)
  real    real sherd shown inside the completed vessel (Samples 2 and 4)
  figs    all remaining figures (make_figures.py)
  numbers FINAL_NUMBERS.csv / .md
  verify  compares the fresh results with the reference results shipped in ../03_numbers
Every step writes to WORK/<step>/ and a log to WORK/logs/<step>.log. Random choices use fixed seeds, so repeated runs agree.
"""
import argparse, os, subprocess, sys, time

def main():
    ap = argparse.ArgumentParser(); ap.add_argument("--data", default="../01_inputs/point_files"); ap.add_argument("--work", default="../work")
    ap.add_argument("--only", nargs="*"); ap.add_argument("--list", action="store_true"); a = ap.parse_args()
    D, W = os.path.abspath(a.data), os.path.abspath(a.work); here = os.path.dirname(os.path.abspath(__file__)); py = sys.executable
    P = lambda *s: os.path.join(W, *s)
    steps = {
        "s2": [py, "sample2_reconstruction.py", "--npz", D + "/sample2.npz", "--out", P("s2")],
        "s4": [py, "sample4_wall_reconstruction.py", "--npz", D + "/sample4.npz", "--out", P("s4")],
        "s3": [py, "sample3_base_reconstruction.py", "--npz", D + "/sample3.npz", "--out", P("s3")],
        "s3flip": [py, "sample3_base_reconstruction.py", "--npz", D + "/sample3.npz", "--out", P("s3_flip"), "--flip"],
        "s3v": [py, "sample3_vessel_spherical_wall.py", "--profile", P("s3", "profile_base.csv"), "--results", P("s3", "results.json"), "--npz", D + "/sample3.npz", "--out", P("s3v"), "--z", "60", "--z", "80", "--z", "100"],
        "s1": [py, "sample1_base_reconstruction.py", "--npz", D + "/sample1.npz", "--out", P("s1"), "--z", "70", "--z", "90", "--z", "110"],
        "s1pot": [py, "sample1_base_reconstruction.py", "--npz", D + "/sample1.npz", "--out", P("s1_pot"), "--z", "170", "--z", "210", "--z", "250"],
        "real2": [py, "show_real_sherd_in_vessel.py", "--sample", "2", "--npz", D + "/sample2.npz", "--vessel", P("s2", "sample2_vessel_with_illustrative_base_scannerLPS.stl"), "--results", P("s2", "results.json"), "--unmeasured-radius", "74.4", "--out", P("real")],
        "real4": [py, "show_real_sherd_in_vessel.py", "--sample", "4", "--npz", D + "/sample4.npz", "--vessel", P("s4", "sample4_wall_band_scannerLPS.stl"), "--results", P("s4", "results.json"), "--out", P("real")],
        "figs": [py, "make_figures.py", "--npz-dir", D, "--work", W, "--out", P("figures")],
        "numbers": [py, "make_final_numbers.py", "--numbers", os.path.join(here, "..", "03_numbers"), "--out", P("numbers")],
        "verify": [py, "verify_results.py", "--work", W, "--reference", os.path.join(here, "..", "03_numbers")],
    }
    if a.list: print(" ".join(steps)); return
    os.makedirs(P("logs"), exist_ok=True)
    for name in (a.only or steps):
        t0 = time.time(); os.makedirs(P(name), exist_ok=True) if name in ("numbers",) else None
        with open(P("logs", name + ".log"), "w") as lg: rc = subprocess.call(steps[name], cwd=here, stdout=lg, stderr=subprocess.STDOUT)
        print("%-8s %s  (%.0f s)" % (name, "ok" if rc == 0 else "FAILED (see logs/%s.log)" % name, time.time() - t0)); sys.stdout.flush()

if __name__ == "__main__":
    main()
