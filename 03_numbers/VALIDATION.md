# Validation summary

## Synthetic recovery tests (known axis, known noise)
| Pipeline | Test | Result |
|---|---|---|
| Sample 2 (profile-thinness, Taubin) | 4 synthetic sherds, floor + wall with a rim, noise 0.4 mm, 2 deg normal noise | direction error 0.2-6 deg, position error 0.5-9 mm; noise floor 0.27-0.33 mm |
| Sample 4 (wall, least squares) | 4 synthetic wall sherds (radius 100 mm, lean 8 deg, 46 deg arc), noise 0.25 mm | direction error 0.4-1.0 deg, position 0.7-1.7 mm (`python sample4_wall_reconstruction.py --selftest`) |
| Sample 3 (two-face bin criterion) | 3 synthetic bases with a ring, started 6 deg and about 10 mm off | direction error 0.01-0.02 deg, position 0.01-0.02 mm (`python sample3_base_reconstruction.py --selftest`) |
| Extraction (world coordinates) | vectorised versus per-voxel formula, rotated geometry | max difference 0.0 mm (`python extract_pieces.py --selftest`) |

## Checks against the scan itself
* Sample 2: the lip plane normal agrees with the profile axis to 0.6 deg (selection spread 1.65 deg); crest level within about 1 mm over about 50 deg; the plate normal of the body is 49.7 deg from the axis found, so slices levelled by PCA are tilted sections.
* Sample 2: 96 % of the scanned ceramic voxels lie inside the exported measured solid.
* Sample 4: outer edge within 0.3 mm (median) of the profile at every azimuth from -20 to +20 deg (90th percentile 0.4-0.8 mm); 99 % of the scanned ceramic voxels lie inside the exported band.
* Sample 3: same ring shape in the CT sections at -25 and +25 deg; face points a median 0.57 mm from the profile; 97.5 % of the ceramic core inside the exported base; scan versus reconstruction median distance about 0.5 mm (69/85/91 % within 1/2/3 mm).
* Sample 1: interior circle fit 0.045-0.28 mm; scan versus vessel median 0.99 mm (50/74/83 % within 1/2/3 mm) - poorer than Sample 3 because of the non-round outline and remaining lumps.

## Solids
* STL meshes are closed (watertight) and their volumes agree with the polygon geometry to 0.05 %.
* STEP solids are valid (OpenCascade) and agree with the polygon volumes to within 0.3 % (Sample 3 vessels within 0.8 %).

## Determinism
`run_all.py` twice gives identical result files; `verify_results.py` prints "identical to the reference" for the six result files.

## Not validated
* No comparison with an intact vessel or another reconstruction method.
* Sample 4's axis is weakly determined (near-tied alternatives), Sample 1's sphere radius depends on the fit (126.5 versus 176 mm).
* The DICOM reading in extract_pieces.py was not run on the real series by the preparer.
