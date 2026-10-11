# Corrections log - errors found during the work and what replaced them

Read this before trusting any number quoted in an earlier message or an earlier file. The files in this package are the corrected ones.

| # | What was wrong | How it was found | What replaced it |
|---|---|---|---|
| 1 | Sample 2 lip width quoted as 4.6 mm (flat top only) | comparison with the full lip cross-section | 6.7 mm from the outer and inner curves at the lip plane |
| 2 | Sample 2 test for "crust lies against the inner face" almost never fired (crust set excluded a 2-voxel shell) | coverage fraction printed 7.6 % | opening-based crust lumps; coverage 17 %; lip width 7.2 -> 6.7 mm, minimum thickness 6.8 -> 5.9 mm, volumes < 1.5 % |
| 3 | The implicit-curve (Taubin) axis criterion gives a meaningless near-zero scatter on nearly straight profiles | synthetic wall test returned 0.000 mm at wrong axes | least-squares rho(h) criterion for wall pieces (recovers a known axis to <= 1 deg); the earlier "base-like 0.164 mm against wall-like 0.231 mm" comparison for Sample 4 is discarded |
| 4 | Sample 4 wall profile said to be 2-3 mm outside the CT edge | judged by eye from an overlay | measured directly: within 0.3 mm (median) at all azimuths -20..+20 deg |
| 5 | Sample 4: half-maximum shifts derived from a split surface built on a smoothed field | gradient of 3661 HU/mm was not plausible | gradient measured on the raw HU profile across the outer face: 974 HU/mm |
| 6 | Sample 3 first delivered in the wrong orientation (peaked crest up, wide flat face down) | owner: "it is inverted"; CT sections show the interior rising into the wall | corrected orientation (rounded ridge down); the rejected orientation is kept only as figure A of the orientation comparison |
| 7 | Sample 3 shape called "nearly symmetric between the faces" | CT sections | the two faces are asymmetric (flat wide face versus peaked crest) |
| 8 | Sample 1 brown piece called crust (546 cm3) | owner: it is a base; density check: mean 1214 HU against 1254 HU for the ceramic group | the piece is ceramic (414 cm3 isolated); only the lumps on it are crust |
| 9 | Sample 1 first stood on the wrong side (rough side down, debris under the base) | owner: "this is not the base" | stood on the round flat foot; debris below the base plane cut |
| 10 | Sample 3 angle-based wall scenarios (STEP volumes off by 12-19 %) | STEP versus mesh volume | withdrawn; walls now follow the interior sphere; STEP within 0.8 % |
| 11 | Sample 2 STEP built from splines through dense noisy points (3-8 % volume error) | STEP versus polygon volume | corner-aware approximating splines; within 0.3 % |
| 12 | Early sketches of Sample 2 with a flaring wall | CT sections | the wall leans inward toward the lip; the sketches are excluded from this package |
| 13 | Quick outline circle fit of Sample 3 (54.5 mm) | only 53 of 350 outline points fitted | not used; the standing ring radius is measured from the profile (53.5 mm) |
| 14 | First-draft abstract: Pearson r > 0.97 for the Main Group | heatmap values | 0.915-0.996 (Figure 3 of the draft) |
| 15 | First-draft Table 2 "Volume" read as vessel volume | volume = voxel count x 0.9537 mm3 | it is the volume of everything at or above the threshold (ceramic plus crust) |
| 16 | First-draft Methods: "1.00 mm isotropic" | DICOM tags | voxels are 0.977 x 0.977 x 1.0 mm |
| 17 | Draft of Results 4.2 said r > 0.97 holds for "six of the ten pairs" | reading the heatmap again while compiling | the Main Group has six pairs; four exceed 0.97, two (0.942, 0.915) do not |
| 18 | A garbled sentence in Results 4.4.7 (units) and a "50 deg" lip-tilt limit stated without basis | proof-reading | reworded; limit stated as 30-50 deg (the normal-alignment cut-offs 0.6-0.85 of the lip selections) |
| 19 | VALIDATION.md quoted a PCA-levelled scatter of 0.73 mm for Sample 2 that was not re-verified | proof-reading | removed; replaced by the verified angle between the plate normal and the axis (49.7 deg) |
