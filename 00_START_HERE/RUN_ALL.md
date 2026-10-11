# Rebuilding everything

```
cd 02_code
pip install -r ../00_START_HERE/requirements.txt
python run_all.py --data ../01_inputs/point_files --work ../work
```
* `run_all.py --list` shows the steps; `--only s3 s3v` runs a subset.
* Each step writes into `work/<step>/` and a log into `work/logs/<step>.log`.
* The last step, `verify`, compares the six result files with the references in `03_numbers/` and prints "identical to the reference" for each.
* Figures are written to `work/figures/` (data views, comparisons, cross-sample summaries) and inside each step folder (profile, sweep, CT-section figures).
* `05_3d_models/` holds the files exactly as shipped; a fresh run reproduces them in `work/<step>/` (STL/STEP/DXF).

Where each figure of the paper comes from: `03_numbers/FIGURE_INDEX.md`.
Where each number of the paper comes from: `03_numbers/FINAL_NUMBERS.md` (built from the result files by `make_final_numbers.py`).

The chain before the point files (DICOM series -> point files) is in `01_inputs/extraction/`; the DICOM series itself is not part of this package.
