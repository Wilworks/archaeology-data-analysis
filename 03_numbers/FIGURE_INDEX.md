# Figure index (files in 04_figures/, and the copy used in 07_manuscript/figures/)

| Folder / file | Shows | Produced by |
|---|---|---|
| 00_data_views/hu_histograms_samples1to4.png | HU histograms with the split thresholds, and the high-density tails | make_figures.py (hu) |
| 00_data_views/surface_views_sample{1,2,3,4}.png | six shaded views of each specimen (500 HU surface; pale = ceramic, brown = crust / the brown piece) | make_figures.py (views1-4) |
| 00_data_views/sample3_height_maps.png | lower face, upper face and extent maps of the Sample 3 plate | make_figures.py (heightmaps3) |
| 00_data_views/sample1_parts.png | Sample 1 cluster split into three ceramic pieces and the brown piece | make_figures.py (parts1) |
| 00_data_views/sample1_piece_isolated_vs_cleaned.png | the brown piece alone, before and after removing the debris | make_figures.py (isolation1) |
| 00_data_views/sample5_stone_screenshot.png | Sample 5 (a stone) in 3D Slicer, supplied by the owner | screenshot; no code |
| 01_sample2/sample2_lip_check.png | the finished lip is flat and level around the axis | sample2_reconstruction.py |
| 01_sample2/sample2_ct_sections_with_profile.png | CT sections through the axis with the extracted profiles | sample2_reconstruction.py |
| 01_sample2/sample2_profile_thickness.png | outer and inner profiles, wall thickness, scatter | sample2_reconstruction.py |
| 01_sample2/sample2_vessel_sweep.png | swept vessel, cutaway, top view, cross-section | sample2_reconstruction.py |
| 01_sample2/sample2_real_sherd_in_vessel.png | real sherd (clay) inside the completed bowl (blue), unmeasured centre grey | show_real_sherd_in_vessel.py (real2) |
| 02_sample3/sample3_orientation_A_vs_B.png | the two orientations of the base (A rejected, B accepted) | make_figures.py (orient3) |
| 02_sample3/sample3_base_ct_sections_with_profile.png | CT sections at -25 and +25 deg with the profile, and both faces plotted as (radius, height) | sample3_base_reconstruction.py |
| 02_sample3/sample3_base_sweep.png | swept base, cross-section | sample3_base_reconstruction.py |
| 02_sample3/sample3_vessel_spherical_wall_with_real_sherd.png | real sherd inside vessels of three wall heights (sphere continuation) | sample3_vessel_spherical_wall.py |
| 02_sample3/sample3_scan_vs_reconstruction.png | scan and reconstruction from the same camera, and a distance map | make_figures.py (scanrecon3) |
| 03_sample1/sample1_piece_on_its_base.png | the cleaned brown piece standing on its foot | make_figures.py (onbase1) |
| 03_sample1/sample1_profile_both_faces.png | both faces as (radius, height), the interior and exterior spheres, upper face by azimuth | sample1_base_reconstruction.py |
| 03_sample1/sample1_vessel_bowls_with_real_piece.png | real piece inside bowls of three heights | sample1_base_reconstruction.py |
| 03_sample1/sample1_pot_projection.png | the same projected to pot heights; green = above 110 mm | sample1_base_reconstruction.py |
| 04_sample4/sample4_wall_ct_sections_with_profile.png | CT sections at -14, 0, +14 deg with the profile | sample4_wall_reconstruction.py |
| 04_sample4/sample4_wall_band_sweep.png | swept wall band | sample4_wall_reconstruction.py |
| 04_sample4/sample4_real_sherd_in_wall_band.png | real sherd inside the band | show_real_sherd_in_vessel.py (real4) |
| 05_cross_sample/cross_sample_profiles.png | half cross-sections of all reconstructions | make_figures.py (profiles) |
| 06_draft_figures/draft_fig{1..4}_*.png | the four figures of the first draft (CT overview, INAA profiles, correlation heatmap, Sample 2 HU distribution), extracted from the draft PDF | not reproducible from this package (INAA data are not included) |

The scripts named "provenance" in 02_code/provenance produced early exploratory pictures (parts, plateau, foot fits); the final versions of those figures are made by the scripts above.
