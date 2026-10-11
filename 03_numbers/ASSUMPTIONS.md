# Register of assumptions

Everything the reconstructions rely on that is a choice or an inference, not a measurement. Effects are given where they were tested.

| # | Assumption | Where used | Value / rule | Effect if wrong |
|---|---|---|---|---|
| A1 | Each vessel is a surface of revolution (perfectly round about one axis) | all samples | - | Wall departures from round: 4.1 mm (Sample 2 wall), 0.6 mm (Sample 2 base), 0.4 mm (Sample 4); Sample 1 piece is visibly not round (scan-versus-vessel median 0.99 mm, 50 % within 1 mm) |
| A2 | The scanned sherd is representative of the whole circumference | all | 79 deg (S2), 232 deg (S3), 52 deg (S4), 360 deg (S1 piece) scanned | Only the scanned arc is measured; the rest is repeated |
| A3 | Unstored voxels (below 500 HU) are set to 0 HU when building the 500 HU surface | outer surfaces | field = HU where >= 500, else 0 | Surface position uncertain by a fraction of a voxel; covered by the half-maximum correction |
| A4 | Half-maximum edge: true surface lies halfway between air (-1000 HU) and the ceramic plateau | thickness, edges | plateau 1350 (S2), 1238 (S4), 1407 (S3); slope measured on the raw HU profile | Edge shifts 0.29-0.39 mm outward from the 500 HU surface, 1.0-1.2 mm from the split level; thickness uncertain by a few tenths of a mm |
| A5 | Body/crust split threshold | S2 1212 HU, S4 1112 HU, S1 1088 HU | valley of the (smoothed) histogram | Split is shallow (valley 15-56 % below the crust peak): connectivity is used as well |
| A6 | Sample 3 orientation: stands on the rounded ridge; wall root bends upward at the outer edge | S3 | confirmed by the owner of the sherd | Opposite orientation (kept as figure A) makes a peaked crest with steep drops, not a plausible vessel wall |
| A7 | The missing centre of the Sample 2 base | S2 | grey cap: sphere tangent to the surviving base (radius 144 mm, 20.7 mm drop) | illustrative only; no base diameter, height or capacity may be quoted |
| A8 | Wall above the broken edge continues the interior sphere fitted to the base (radius 120 mm S3, 126.5 mm S1) | S3, S1 | exterior sphere concentric, thickness from the outer skirt (10.7 mm S3, 13.6 mm S1) | Any wall height or thickness is a choice; free-centre circle for S1 gives 176 mm instead of 126.5 mm |
| A9 | Wall heights of the completed vessels | S3 z = 60/80/100 mm; S1 bowls 70/90/110, pot 170/210/250 | choices | none measured; shown as scenarios |
| A10 | Up direction of the Sample 4 band | S4 | wall leans outward going up | nearly symmetric top-bottom; orientation is arbitrary |
| A11 | Sample 4 axis | S4 | best of near-tied solutions (0.514 mm; 0.539 mm at 56 deg; 0.588 mm at 35 deg) | radius 87-101 mm across the candidates |
| A12 | Sample 1 brown piece: which material belongs to it | S1 | remaining stored material after removing the ceramic group, opening (ball 3), largest component; debris opening (ball 5); cut below the base plane; cut above the inner-face profile + 4 mm | 106 cm3 + 24 cm3 removed as debris; relief above the inner face of more than 4 mm would be lost |
| A13 | Sample 1 axis comes from the flat foot (plane + outline circle), not from a profile search | S1 | foot circle radius 47.1 mm, RMS 1.18 mm | The two-face criterion gives a poor spread (about 8 mm) for this piece: it is not round enough for a search |
| A14 | Sample identification by voxel count against Table 2 | all | see 01_inputs/extraction/SAMPLE_LABELS.md | Sample 3's bounding box does not match Table 2 in y and z |
