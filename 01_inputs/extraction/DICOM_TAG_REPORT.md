# DICOM tag report (series_0000E3E7)

Source: reported by the author's `02_extract_sample2.py`, Step 1 (transcribed from the author's reproducibility README).

| Tag | Value |
|---|---|
| Rows x Columns | 512 x 512 |
| Number of slices / series | 482 / 1 |
| PixelSpacing | [0.9765625, 0.9765625] mm |
| SliceThickness / SpacingBetweenSlices | 1 mm / 1 mm |
| GantryDetectorTilt | 0 deg |
| ImageOrientationPatient | [1, 0, 0, 0, 1, 0] (axis-aligned) |
| ImagePositionPatient, first slice (sorted along the slice normal) | [-273.74670028686, -307.85266113281, 985.27001953125] mm |
| ImagePositionPatient, last slice | [-273.74670028686, -307.85266113281, 1466.27001953125] mm |
| Slice gaps (min / max / mean) | 1.0000 / 1.0000 / 1.0000 mm |
| RescaleSlope / RescaleIntercept | 1 / -1024 |
| BitsStored / PixelRepresentation | 12 / 0 (unsigned) |
| KVP | 120 |
| ConvolutionKernel | UB |

Consequences used in the paper
* The largest stored value is 4095 - 1024 = 3071 HU: every specimen shows a ceiling at 3071 HU that is a storage limit, not a density.
* Voxels are 0.9766 x 0.9766 x 1.0 mm (voxel volume 0.9537 mm3): slightly anisotropic, not "1.00 mm isotropic".

## World coordinates
For every voxel, using the position tag of its own slice:

    P = IPP_slice + col * PixelSpacing[1] * rowDir + row * PixelSpacing[0] * colDir

with rowDir = IOP[0:3] and colDir = IOP[3:6]. Slices are sorted by the projection of IPP on cross(rowDir, colDir), not by InstanceNumber.

## Segmentation
Threshold 500 HU (inclusive); 26-connected components (`scipy.ndimage.generate_binary_structure(3, 3)`); no smoothing, resampling, PCA or rotation.

Tags that the paper's Methods also quotes but that are NOT in this report and must be checked against the DICOM headers before submission:
scanner manufacturer and model, tube current (the first draft says 424 mA), acquisition site, reconstruction field of view.
