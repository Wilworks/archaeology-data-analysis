# Opening the 3D files

Folders
* `sample2/`, `sample3_base/`, `sample3_vessel_spherical_wall/`, `sample4_wall_band/`, `sample1_bowls/`, `sample1_pot_projection/`, `real_sherds_cleaned/`
* File name endings: `_scannerLPS.stl` = in the CT scan's own coordinates (millimetres) - overlays the CT in 3D Slicer; `_local.stl` = axis along z, resting on z = 0 - for CAD, printing and viewing;
  `.step` = exact solid (revolved splines) for SolidWorks and other CAD; `_profile.dxf` = the half cross-section for a native Revolved Boss/Base.

3D Slicer: load the CT first, then File > Add Data and pick a `_scannerLPS.stl` (Slicer reads STL as LPS, which is what these files use; if the model looks mirrored or far from the scan,
reload through Show Options and choose LPS). Models module > Opacity to see the sherd through the model. The `real_sherds_cleaned` files are the crust-free scanned surfaces (open shells) for the same frame.
SolidWorks: File > Open the `.step` (solid body, mm). DXF: import as a 2D sketch, turn the AXIS line into a centreline, Revolved Boss/Base 360 deg. STL: choose "solid body" in the import options if wanted (slow for large meshes).
Other software: MeshLab, Blender, FreeCAD, Fusion 360 open STL/STEP.

Colour convention used in all figures: clay = the real scanned sherd; blue = symmetric completion / assumed wall; green = further projection above the height already shown; grey = not measured (Sample 2 base centre).

What is measured and what is not: `03_numbers/FINAL_NUMBERS.md` and `03_numbers/ASSUMPTIONS.md`. The files were checked by reading them back with trimesh, OpenCascade (CadQuery) and ezdxf; they were not opened in SolidWorks or 3D Slicer by the preparer.
Sample 3 vessel files: z = 60, 80, 100 mm. Sample 1 bowls: 70, 90, 110 mm; Sample 1 pot projection: 170, 210, 250 mm (all wall heights are choices).
