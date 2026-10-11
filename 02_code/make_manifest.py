#!/usr/bin/env python3
"""make_manifest.py - writes MANIFEST.csv (path, size in bytes, SHA-256, description) for every file of the package.
Usage: python make_manifest.py --root ..   (run from 02_code; the manifest is written to 00_START_HERE/MANIFEST.csv)"""
import argparse, csv, hashlib, os

RULES = [
 ("00_START_HERE/README", "start here: contents, specimens, what to read before quoting a number"), ("00_START_HERE/RUN_ALL", "how to rebuild everything"),
 ("00_START_HERE/requirements", "Python packages"), ("00_START_HERE/ENVIRONMENT", "tested software versions, run times, seeds"),
 ("01_inputs/point_files/", "point file: x, y, z (mm, float32) and hu (int16) of every voxel >= 500 HU of one scan label"),
 ("01_inputs/extraction/extract_pieces", "vectorised extraction script with self-test (DICOM series -> point files)"), ("01_inputs/extraction/DICOM", "DICOM tag report and world-coordinate formula"),
 ("01_inputs/extraction/SAMPLE_LABELS", "how specimens were matched to scan labels and to Table 2"), ("01_inputs/extraction/ORIGINAL", "note on the three original extraction scripts"), ("01_inputs/SAMPLE5", "note on Sample 5 (a stone)"),
 ("02_code/provenance/", "exploratory script that produced an early figure or check (paths are hard-coded; not part of run_all)"),
 ("02_code/run_all", "rebuilds every result, figure and 3D file"), ("02_code/verify_results", "compares a fresh run with the reference results"), ("02_code/make_figures", "data views, histograms, comparisons and summary figures"),
 ("02_code/make_final_numbers", "builds FINAL_NUMBERS from the result files"), ("02_code/make_manifest", "writes this manifest"), ("02_code/show_real_sherd_in_vessel", "real sherd drawn inside the completed vessel (Samples 2 and 4)"),
 ("02_code/sample1_", "Sample 1 (brown piece) pipeline: isolation, foot, profile, sphere, bowls and pots"), ("02_code/sample2_", "Sample 2 pipeline and shared library (axis search, surfaces, sweep, exports)"),
 ("02_code/sample3_base", "Sample 3 ring-footed base pipeline"), ("02_code/sample3_vessel", "Sample 3 wall pulled up along the interior sphere"), ("02_code/sample4_", "Sample 4 wall-band pipeline"),
 ("03_numbers/FINAL_NUMBERS", "every reported number with uncertainty and status"), ("03_numbers/ASSUMPTIONS", "register of assumptions"), ("03_numbers/VALIDATION", "validation summary"),
 ("03_numbers/CORRECTIONS_LOG", "errors found and fixed"), ("03_numbers/FIGURE_INDEX", "which script makes which figure"),
 ("03_numbers/", "result file or profile table of one pipeline (see the file name)"),
 ("04_figures/", "figure (see 03_numbers/FIGURE_INDEX.md)"),
 ("05_3d_models/OPEN_IN", "how to open the 3D files"), ("05_3d_models/real_sherds_cleaned", "crust-free scanned sherd surface in scanner coordinates (open shell)"),
 ("05_3d_models/", "3D model or CAD file: _scannerLPS.stl overlays the CT, _local.stl rests on z = 0, .step exact solid, .dxf half cross-section"),
 ("06_methods_documents/", "methods note (Word) and hand-out text for the Sample 2 pipeline"),
 ("07_manuscript/main.tex", "manuscript master file"), ("07_manuscript/main.pdf", "compiled manuscript"), ("07_manuscript/sections/", "manuscript section (LaTeX)"), ("07_manuscript/references.bib", "bibliography (entries to be checked)"),
 ("07_manuscript/figures/", "figure as included in the manuscript (copy of a file in 04_figures)"), ("07_manuscript/TABLE2", "corrections needed in the first draft's Table 2 and abstract"), ("07_manuscript/TO_VERIFY", "items the preparer could not verify"),
]

def main():
    ap = argparse.ArgumentParser(); ap.add_argument("--root", default=".."); a = ap.parse_args(); root = os.path.abspath(a.root); rows = []
    for d, _, fs in os.walk(root):
        for f in sorted(fs):
            p = os.path.join(d, f); rel = os.path.relpath(p, root).replace(os.sep, "/")
            if rel.endswith("MANIFEST.csv"): continue
            h = hashlib.sha256()
            with open(p, "rb") as fh:
                for chunk in iter(lambda: fh.read(1 << 20), b""): h.update(chunk)
            desc = next((t for pre, t in RULES if rel.startswith(pre)), "")
            rows.append((rel, os.path.getsize(p), h.hexdigest(), desc))
    rows.sort()
    with open(os.path.join(root, "00_START_HERE", "MANIFEST.csv"), "w", newline="") as fh:
        w = csv.writer(fh); w.writerow(["path", "bytes", "sha256", "description"]); w.writerows(rows)
    print("manifest:", len(rows), "files,", sum(r[1] for r in rows) // 1000000, "MB")

if __name__ == "__main__":
    main()
