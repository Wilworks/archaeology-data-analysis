SAMPLE 2 - RECONSTRUCTED VESSEL: 3D FILES AND HOW TO USE THEM
================================================================

Prepared: 19 September 2026
Prepared by / contact: [add name and email here]
Related manuscript (draft): "17th-Century Ceramics from Fort Amsterdam:
A Digital 3D and Geochemical Perspective", University of Ghana.


1. WHAT THIS IS
---------------
Sample 2 is one broken ceramic sherd from Fort Amsterdam (Abandze/Kormantin,
Ghana). It was scanned with a clinical CT scanner (512 x 512 pixels, 482
slices, voxels 0.977 x 0.977 x 1.0 mm, 120 kVp). From the CT data we worked
out where the pot's axis of rotation is, extracted the shape (profile) of
the sherd's wall, and spun that profile 360 degrees to build a model of the
whole vessel.

These files are that model. Only about a fifth of the vessel (roughly 79
degrees of the circumference) was actually scanned. The rest is the same
shape repeated by symmetry. The centre of the base was not on the sherd and
was NOT measured.

All lengths are in millimetres.


2. WHICH FILE SHOULD I OPEN?
----------------------------
To see the complete-looking pot:
   3D Slicer   ->  sample2_vessel_with_illustrative_base_scannerLPS.stl
   SolidWorks  ->  sample2_vessel_with_illustrative_base.step

To see only what the scan supports (open hole in the middle):
   3D Slicer   ->  sample2_vessel_measured_scannerLPS.stl
   SolidWorks  ->  sample2_vessel_measured.step


3. THE FILES
------------
sample2_vessel_measured_scannerLPS.stl
   The measured vessel: a ring-shaped bowl with a hollow centre. Triangle
   mesh, 51,840 triangles, closed (watertight). In the same coordinates as
   the CT scan, so it can be laid over the scan in 3D Slicer.

sample2_vessel_with_illustrative_base_scannerLPS.stl
   The same vessel plus a smooth rounded bottom that closes the centre.
   73,728 triangles, closed. Same coordinates as the CT scan.
   THE BOTTOM IS AN ILLUSTRATION, NOT A MEASUREMENT (see section 6).

sample2_vessel_measured_local.stl
sample2_vessel_with_illustrative_base_local.stl
   The same two shapes, but moved so the vessel's axis is the vertical z axis
   through the origin. Use these for CAD, measuring, viewing in MeshLab or
   Blender, or 3D printing.
      measured:  z = 0 at the lowest surviving point of the sherd
                 (a broken edge), rim at z = 59.7 mm.
      with base: z = 0 at the bottom of the illustrative base, so the model
                 stands on the z = 0 plane; total height 80.4 mm.

sample2_vessel_measured.step
sample2_vessel_with_illustrative_base.step
   Exact CAD solids (curved surfaces, not triangles): 5 faces and 4 faces.
   Local coordinates as for the "local" STL files. Best choice for
   SolidWorks, Fusion 360, FreeCAD and other CAD programs.

sample2_vessel_profile.dxf
   The half cross-section (the shape that is spun to make the vessel) as a
   2D drawing. X = distance from the axis, Y = height above the lowest
   surviving point, units mm. Layers:
      PROFILE_MEASURED        the measured wall (curved lines and two
                              straight ends)
      AXIS                    the vertical axis line at X = 0
      ILLUSTRATIVE_BASE_CAP   the made-up bottom (dashed); ignore or delete it
                              if you want only the measured part
   Use it to make a Revolved Boss/Base in SolidWorks.

Supporting files (if they were sent with this package):
   sample2_reconstruction.py      the program that produces everything
   sample2_results.json           all numbers in one file
   sample2_profile.csv            the profile curves as numbers
   sample2_vessel_sweep.png       pictures of the whole vessel
   sample2_ct_sections_with_profile.png   real CT slices with the profile drawn on
   sample2_profile_thickness.png  the extracted profile and wall thickness
   sample2_methods_math.docx      the mathematics, step by step


4. HOW TO OPEN THEM
-------------------
3D Slicer
   1. Load the CT scan first (the DICOM folder), as usual.
   2. File > Add Data (or drag the STL into the window) and choose one of the
      "scannerLPS" STL files.
   3. The vessel should sit exactly on the scanned sherd. Slicer normally
      reads STL files in LPS coordinates, which is what these files use. If
      the model looks mirrored or far from the scan, load it again through
      Add Data > Show Options and set the coordinate system to LPS.
   4. To see the sherd through the model: Modules > Models, select the model,
      and lower its Opacity.

SolidWorks
   STEP (recommended): File > Open, choose the .step file. It opens as a
      solid part in millimetres.
   DXF: File > Open, choose the .dxf, and import it into a new part as a 2D
      sketch (units mm). Turn the AXIS line into a centreline, then use
      Revolved Boss/Base with 360 degrees.
   STL: File > Open, file type STL, and choose "solid body" in the import
      options if you want a solid. Large meshes are slow, so prefer STEP.

Other programs
   MeshLab, Blender, FreeCAD and Fusion 360 open the STL and STEP files.
   For 3D printing, use a "local" STL file. The models are closed, so slicing
   software should accept them.

Note: these files were checked by reading them back with other software
(mesh, CAD and DXF libraries). They have not yet been opened in SolidWorks or
3D Slicer by the person who prepared them. If one fails to open, try the other
format in the same row of section 2.


5. WHAT THE FILES CAN HELP WITH
-------------------------------
   - Seeing the reconstructed vessel next to the original CT scan.
   - Measuring the wall thickness, the rim and the shape of the profile.
   - Making a CAD model that can be edited, drawn or compared with other
     vessels.
   - Producing a 3D-printed replica for teaching or display.
   - Drawing figures for the manuscript.


6. WHAT IS MEASURED, WHAT IS NOT
--------------------------------
Measured from the scan (Sample 2, with rough uncertainty):
   Rim outer diameter                 about 228 mm     (+/- 3 mm)
   Widest outer diameter (shoulder)   about 243 mm     (+/- 3 mm), at about
                                      37 mm above the lowest surviving point
   Rim height above lowest point      about 59.7 mm    (a minimum; that lowest
                                      point is a broken edge, so it is not the
                                      true base)
   Lip (rim top) radial width         about 6.7 mm
   Wall thickness                     median 9.7 mm, range 5.9 to 12.7 mm
                                      (about +/- 0.7 mm)
   Distance of the lowest surviving edge from the axis: about 74 mm
   Slope of the base at that edge: about 31 degrees from horizontal
   Axis of rotation (scanner coordinates, mm): passes through
      (145.0, -92.6, 1285.9), direction (-0.315, -0.773, 0.550),
      tilted 56.6 degrees from the scanner z axis.

NOT measured:
   - The centre of the base (inside about 74 mm from the axis). The rounded
     bottom in the "with_illustrative_base" files is a smooth cap that just
     touches the measured base (sphere radius about 144 mm, dropping about
     21 mm below the lowest surviving point). It is for display only. Do not
     use it to quote a base diameter, a total height or a capacity.
   - The true vessel height.
   - Anything about the parts of the vessel that were not on the sherd. They
     are copies of the scanned part, assuming a perfectly round pot.

Other cautions:
   - The wall is not perfectly round about the fitted axis: it varies by up
     to about 4 mm around the scanned arc (the base varies by 0.6 mm). This
     is either a small error in the axis or a slightly out-of-round pot, and
     the data cannot tell which. That is why the diameters carry +/- 3 mm.
   - Wall thickness and edge positions use a standard CT edge assumption
     (the true surface lies half way between air and the ceramic). It depends
     on assumed values, so treat thickness as accurate to a few tenths of a
     millimetre at best.
   - The outer surface is placed at the CT edge; the inner surface is
     the wall thickness measured inwards from it.


7. WHERE THE NUMBERS CAME FROM
------------------------------
Input: the CT voxels of the sherd at or above 500 HU (sample2.npz). The
adhering crust was separated from the ceramic body at the valley between
the two density peaks (about 1212 HU). The axis was found by choosing the
axis that makes the wall's (distance, height) profile as thin as possible,
and it was checked independently with a plane fitted to the flat, finished
rim (the two agree to 0.6 degrees). The profile was then measured, corrected
to the half-way edge, given a wall thickness, and spun 360 degrees. The full
method and all equations are in sample2_methods_math.docx. The program
sample2_reconstruction.py repeats the whole process and gives identical
numbers each time it is run.


8. CHECKS ALREADY DONE
----------------------
   - Meshes are closed (no holes) and their volumes match the geometry to
     0.05 percent.
   - The CAD solids are valid and their volumes match the meshes to about
     0.3 percent.
   - About 96 percent of the scanned ceramic voxels lie inside the measured
     solid, which confirms the "scannerLPS" files line up with the CT data.
   - Modelled from one sherd only. Nothing here has been compared with
     another vessel or with an intact one.

End of file.
