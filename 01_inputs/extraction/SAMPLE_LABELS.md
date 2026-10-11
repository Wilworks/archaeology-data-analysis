# How the specimens were matched to scan labels

Pieces were matched to the rows of Table 2 of the first draft **by voxel count only** (largest remaining piece first); no physical
record of the match exists in the extraction scripts. The table below records what the point files contain.

| Sample | Scan label | Voxels (>= 500 HU) | Table 2 voxels | Difference | Bounding box here (mm) | Bounding box in Table 2 (mm) | Remark |
|---|---|---|---|---|---|---|---|
| 1 | 2 | 977,093 | 1,017,293 | -4.0 % | 188.5 x 125.0 x 179.0 | 188.5 x 125.0 x 179.0 | identical box; the piece is a **cluster** (see below) |
| 2 | 424 | 217,925 | 232,328 | -6.2 % | 142.6 x 80.1 x 92.0 | 142.6 x 80.1 x 92.0 | identical box |
| 3 | 62 | 158,229 | 165,240 | -4.2 % | 124.0 x 51.8 x 116.0 | 124.0 x 53.7 x 113.0 | **boxes differ in y and z** |
| 4 | 440 | 137,436 | 150,759 | -8.8 % | 130.9 x 61.5 x 94.0 | 130.9 x 61.5 x 94.0 | identical box; Table 2 calls it "Base", the owner of the sherd identifies it as an outer wall |
| 5 | 491 | 67,438 | 75,587 | -10.8 % | not extracted | not extracted | identified by the owner as a **stone**; not reconstructed |

Every voxel count is 4-11 % below Table 2, consistently: Table 2 was evidently made with a lower HU threshold than 500.

Sample 1 is a cluster: after a watershed split the ceramic group contains a main slab (about 294 cm3), a second piece (about 80 cm3) and
a small piece (about 10 cm3), all in contact with a brown piece of about 414 cm3 that is a ceramic piece in its own right (mean 1214 HU
against 1254 HU for the ceramic group). The blue slab looks like Sample 3 in form but is not scan label 62.
