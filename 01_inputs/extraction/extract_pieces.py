#!/usr/bin/env python3
"""
extract_pieces.py - clean, vectorised version of the three extraction scripts used for this study
(02_extract_sample2.py, 03_extract_sample4.py, 04_extract_samples1_3.py), with the DICOM folder and the output folder given
on the command line instead of being hard-coded.

What it does (same rules as the originals: report only, no fitting, no PCA, no smoothing, no resampling, no rotation)
  1. reads a single-series CT folder, sorts the slices by the projection of ImagePositionPatient on the slice normal
     (normal = rowDir x colDir from ImageOrientationPatient), converts to Hounsfield units (slope * stored + intercept);
  2. thresholds at >= 500 HU and labels 26-connected pieces;
  3. for every requested piece, prints the statistics, the HU histogram (25 HU bins) and the world-coordinate bounding box,
     and saves <name>.npz with x, y, z (float32, mm) and hu (int16) for every voxel of that piece;
  4. reloads each file and prints its voxel count and coordinate ranges (they must match step 3).

World coordinates (per slice, using that slice's own ImagePositionPatient):
     P = IPP_slice + col * PixelSpacing[1] * rowDir + row * PixelSpacing[0] * colDir

Usage
  python extract_pieces.py --dicom-dir DIR --out-dir OUT --piece sample2=232328 --piece sample4=137436
  python extract_pieces.py --dicom-dir DIR --out-dir OUT --piece sample1=label:2      # choose by scan label instead of size
  python extract_pieces.py --selftest                                                  # no DICOM needed (see below)
Pieces given as name=N are matched to the label whose voxel count is closest to N; the script stops if the match is off by
more than 10 percent. Requirements: numpy, scipy, pydicom.

Status of this file: the world-coordinate function is tested against the per-voxel formula of the original scripts on
synthetic geometry (--selftest). The DICOM reading part follows the original scripts line by line but was NOT run on the real
CT series by the author of this file.
"""
import argparse, os, sys
import numpy as np


def world_coords(ipps, iop, ps_row, ps_col, zi, yi, xi):
    """Vectorised P = IPP_slice + col*ps_col*rowDir + row*ps_row*colDir for voxel indices (zi, yi, xi)."""
    row_dir = np.asarray(iop[:3], float); col_dir = np.asarray(iop[3:], float)
    return ipps[zi] + xi[:, None] * ps_col * row_dir[None, :] + yi[:, None] * ps_row * col_dir[None, :]


def selftest():
    rng = np.random.default_rng(0)
    n_sl = 60
    iop = [0.99, 0.10, 0.0, -0.10, 0.99, 0.0]                                    # a slightly rotated (non-axis-aligned) geometry
    ipps = np.column_stack([-270 + 0.3 * np.arange(n_sl), -300 + 0.1 * np.arange(n_sl), 985.27 + 1.0 * np.arange(n_sl)])
    zi, yi, xi = rng.integers(0, n_sl, 5000), rng.integers(0, 512, 5000), rng.integers(0, 512, 5000)
    fast = world_coords(ipps, iop, 0.9765625, 0.9765625, zi, yi, xi)
    row_dir_arr = np.array(iop[:3]); col_dir_arr = np.array(iop[3:])
    slow = np.empty_like(fast)
    for vi in range(len(zi)):                                                     # the per-voxel loop of the original scripts
        slow[vi] = ipps[zi[vi]] + xi[vi] * 0.9765625 * row_dir_arr + yi[vi] * 0.9765625 * col_dir_arr
    err = np.abs(fast - slow).max()
    print("selftest: max difference between the vectorised and the per-voxel formula = %.3e mm (%s)" % (err, "PASS" if err < 1e-9 else "FAIL"))
    return err < 1e-9


def load_series(dicom_dir):
    import pydicom
    datasets = []
    for f in sorted(os.listdir(dicom_dir)):
        p = os.path.join(dicom_dir, f)
        if os.path.isfile(p):
            try: datasets.append((p, pydicom.dcmread(p, stop_before_pixels=True)))
            except Exception: pass
    iop = np.array(datasets[0][1].ImageOrientationPatient, float); normal = np.cross(iop[:3], iop[3:])
    datasets.sort(key=lambda t: float(np.dot(np.array(t[1].ImagePositionPatient, float), normal)))
    ds0 = datasets[0][1]; vol = np.empty((len(datasets), int(ds0.Rows), int(ds0.Columns)), dtype=np.int16)
    for i, (p, _) in enumerate(datasets):
        full = pydicom.dcmread(p)
        vol[i] = (full.pixel_array.astype(np.float32) * float(full.RescaleSlope) + float(full.RescaleIntercept)).astype(np.int16)
    ipps = np.array([[float(v) for v in d.ImagePositionPatient] for _, d in datasets], float)
    return vol, ipps, iop, float(ds0.PixelSpacing[0]), float(ds0.PixelSpacing[1])


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--dicom-dir"); ap.add_argument("--out-dir", default="."); ap.add_argument("--threshold", type=float, default=500.0)
    ap.add_argument("--piece", action="append", default=[], help="name=voxel_count  or  name=label:N (repeatable)")
    ap.add_argument("--selftest", action="store_true")
    a = ap.parse_args()
    if a.selftest: sys.exit(0 if selftest() else 1)
    from scipy import ndimage
    os.makedirs(a.out_dir, exist_ok=True)
    vol, ipps, iop, ps_row, ps_col = load_series(a.dicom_dir)
    mask = vol >= a.threshold
    labeled, n = ndimage.label(mask, structure=ndimage.generate_binary_structure(3, 3))
    sizes = np.array(ndimage.sum(mask, labeled, range(1, n + 1)), dtype=np.int64)
    print("volume %s; %d connected pieces at >= %.0f HU; five largest: %s" % (vol.shape, n, a.threshold, [(int(i) + 1, int(sizes[i])) for i in np.argsort(sizes)[::-1][:5]]))
    for spec in a.piece:
        name, val = spec.split("=")
        if val.startswith("label:"): lbl = int(val.split(":")[1])
        else:
            lbl = int(np.argmin(np.abs(sizes - int(val)))) + 1
            diff = 100.0 * abs(int(sizes[lbl - 1]) - int(val)) / int(val); print("%s: closest label %d (%d voxels, %.2f%% from the target)" % (name, lbl, sizes[lbl - 1], diff))
            if diff > 10.0: sys.exit("difference above 10 percent - stopping; check the target")
        m = labeled == lbl; zi, yi, xi = np.where(m); hu = vol[m]
        P = world_coords(ipps, iop, ps_row, ps_col, zi, yi, xi).astype(np.float32)
        print("%s: %d voxels, HU min %d max %d (at max %d), mean %.2f sd %.2f" % (name, len(zi), hu.min(), hu.max(), int((hu == hu.max()).sum()), hu.mean(), hu.std()))
        h, e = np.histogram(hu, bins=np.arange(a.threshold, hu.max() + 26, 25))
        print("  histogram (bin start, count):", [(int(lo), int(c)) for lo, c in zip(e[:-1], h)][:80], "...")
        print("  bounding box (mm): x %.3f..%.3f, y %.3f..%.3f, z %.3f..%.3f" % (P[:, 0].min(), P[:, 0].max(), P[:, 1].min(), P[:, 1].max(), P[:, 2].min(), P[:, 2].max()))
        out = os.path.join(a.out_dir, name + ".npz"); np.savez(out, x=P[:, 0], y=P[:, 1], z=P[:, 2], hu=hu.astype(np.int16))
        v = np.load(out); print("  reloaded %s: %d voxels; x %.3f..%.3f, y %.3f..%.3f, z %.3f..%.3f" % (out, len(v["x"]), v["x"].min(), v["x"].max(), v["y"].min(), v["y"].max(), v["z"].min(), v["z"].max()))


if __name__ == "__main__":
    main()
