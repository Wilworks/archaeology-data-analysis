# 17th-Century Broken Ceramics from Fort Amsterdam: A 3D Reconstruction

[![Python 3.10+](https://img.shields.io/badge/python-3.10+-blue.svg)](https://www.python.org/downloads/)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

This repository contains the complete open-source computational pipeline, CT point-cloud datasets, numerical outputs, and reconstructed 3D digital CAD models accompanying the research paper:

> **17th-Century Broken Ceramics from Fort Amsterdam: A 3D Reconstruction**  
> *Department of Biomedical Engineering, School of Engineering Sciences, University of Ghana, Legon, Ghana.*

---

## 🏛️ Archaeological Context & Assemblage

Excavations at Fort Amsterdam (established in 1638 at Abandze/Kormantin, Central Region, Ghana) recovered five heavily encrusted, high-radiodensity ($1037\text{--}1505\,\text{HU}$) material culture specimens from subterranean vaulted chambers. This repository provides the deterministic, non-destructive geometric reconstruction pipeline that resolves their physical identities:

* **Sample 1 (`sample1.npz`):** Agglomerated earthenware cluster containing a pedestal-footed dish ($\varnothing 94.2\,\text{mm}$ standing foot; interior sphere radius $R \approx 126.5\,\text{mm}$).
* **Sample 2 (`sample2.npz`):** Akan carinated serving bowl fragment ($\varnothing 228\,\text{mm}$ rim diameter, $\varnothing 243\,\text{mm}$ widest shoulder, $9.7\,\text{mm}$ median wall thickness).
* **Sample 3 (`sample3.npz`):** Ring-footed vessel base standing on an annular ring ($\varnothing \approx 107\,\text{mm}$; interior sphere radius $R \approx 120.2\,\text{mm}$).
* **Sample 4 (`sample4.npz`):** Convex outer wall band ($112.5\,\text{mm}$ vertical height; $\varnothing 189\text{--}218\,\text{mm}$ outer diameter) from a globular storage jar.
* **Sample 5 (`sample5.npz`):** Water-worn smooth quartz pebble (mean $1505\,\text{HU}$) curated as an Akan potter's burnishing implement (*mogyemogye*, non-ceramic; non-rotational).

---

## 📁 Repository Structure

```
archaeology-data-analysis/
├── 00_START_HERE/              # Quickstart guide, environment notes, requirements, and SHA-256 manifest
│   ├── README.txt              # Detailed overview of the package
│   ├── RUN_ALL.md              # Instructions for full automated pipeline rebuild
│   ├── MANIFEST.csv            # File manifest with cryptographic SHA-256 hashes
│   └── requirements.txt        # Python dependency specifications
├── 01_inputs/                  # Input datasets and extraction provenance
│   ├── point_files/            # NumPy point clouds (sample1.npz - sample5.npz)
│   └── extraction/             # DICOM extraction scripts and tag reports
├── 02_code/                    # Full Python reconstruction codebase
│   ├── sample1_base_reconstruction.py
│   ├── sample2_reconstruction.py
│   ├── sample3_base_reconstruction.py
│   ├── sample3_vessel_spherical_wall.py
│   ├── sample4_wall_reconstruction.py
│   ├── run_all.py              # Master execution pipeline
│   ├── verify_results.py       # Deterministic verification against references
│   ├── make_figures.py         # Figure generation script
│   └── make_final_numbers.py   # Numerical table compilation script
├── 03_numbers/                 # Numerical results and methodological registers
│   ├── FINAL_NUMBERS.md        # Reference dimensions with measurement status (M, E, A)
│   ├── ASSUMPTIONS.md          # Methodological register of geometric assumptions (A1-A14)
│   ├── VALIDATION.md           # Multi-stage validation tolerances and consistency checks
│   └── *_results.json          # Machine-readable output metrics for each specimen
├── 04_figures/                 # Generated diagnostic plots and multi-view renders
├── 05_3d_models/               # Reconstructed 3D CAD assets (open formats)
│   ├── real_sherds_cleaned/    # Cleaned surface meshes of authentic scanned sherds (STL)
│   ├── sample1_bowls/          # Reconstructed dish models (STL, STEP, DXF)
│   ├── sample1_pot_projection/ # Extended globular pot models (STL, STEP, DXF)
│   ├── sample2/                # Carinated bowl solid reconstructions (STL, STEP, DXF)
│   ├── sample3_base/           # Ring base solids (STL, STEP, DXF)
│   ├── sample3_vessel_spherical_wall/ # Reconstructed base-and-wall vessels (STL, STEP, DXF)
│   └── sample4_wall_band/      # Swept cylindrical wall band models (STL, STEP, DXF)
└── 06_methods_documents/       # Mathematical documentation on profile thinness & axis fitting
```

---

## 🚀 Quickstart & Pipeline Reproduction

### 1. Environment Setup

Clone this repository and create a Python virtual environment:

```bash
git clone https://github.com/Wilworks/archaeology-data-analysis.git
cd archaeology-data-analysis

python3 -m venv .venv
source .venv/bin/activate
pip install -r 00_START_HERE/requirements.txt
```

### 2. Run the Full Reconstruction Pipeline

To run the complete pipeline and regenerate all models, figures, and numerical tables:

```bash
cd 02_code
python run_all.py --data ../01_inputs/point_files --work ../work
```

* To inspect available steps: `python run_all.py --list`
* To execute specific specimens: `python run_all.py --data ../01_inputs/point_files --work ../work --only s1 s2`

### 3. Verify Numerical Reproducibility

To verify that the freshly computed results match the peer-reviewed reference numbers:

```bash
python verify_results.py --work ../work --reference ../03_numbers
```

All algorithms utilize fixed random seeds, ensuring strict numerical determinism.

---

## 🖥️ Viewing the 3D Models

All reconstructed digital assets are provided in open, non-proprietary CAD and mesh formats:
* **STL (`.stl`):** Watertight 2-manifold surface meshes compatible with 3D Slicer, MeshLab, Blender, and 3D printing slicers (PrusaSlicer, Cura).
* **STEP (`.step`):** Exact boundary-representation (B-rep) solids compatible with FreeCAD, SolidWorks, Rhino, and Autodesk Fusion 360.
* **DXF (`.dxf`):** Standard 2D CAD meridian directrix profiles with explicit metric dimensions.

See [`05_3d_models/OPEN_IN_SLICER_AND_SOLIDWORKS.md`](05_3d_models/OPEN_IN_SLICER_AND_SOLIDWORKS.md) for step-by-step viewing instructions.

---

## 🔬 Data Archiving Statement

* **CT Point Clouds & 3D Models:** Available directly in this repository under `01_inputs/point_files/` and `05_3d_models/`.
* **DICOM CT Scans & INAA Spectral Data:** The original medical DICOM volume series ($482$ slices) and raw INAA multichannel gamma spectra are archived at the University of Ghana Medical Centre (UGMC) and Ghana Atomic Energy Commission (GAEC), accessible upon academic request.

---

## 📄 Citation & Correspondence

For inquiries regarding the computational pipeline, please contact:
* **Department of Biomedical Engineering**, University of Ghana, Legon: `cfadams@ug.edu.gh` / `waasumboya@st.ug.edu.gh`
