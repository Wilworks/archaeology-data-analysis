FORT AMSTERDAM CERAMICS - FULL PACKAGE
=======================================

Contents (open in this order)
  00_START_HERE     this file, MANIFEST.csv (every file with its SHA-256), requirements, environment, how to rerun
  01_inputs         point files (sample1-4.npz), extraction record (DICOM tag report, label matching, cleaned extraction script with a self-test), Sample 5 note
  02_code           every script (run_all.py rebuilds everything), provenance scripts of early exploratory figures
  03_numbers        result files and profile CSVs of every pipeline, FINAL_NUMBERS (csv/md), ASSUMPTIONS, VALIDATION, CORRECTIONS_LOG, FIGURE_INDEX
  04_figures        all figures grouped by specimen (see 03_numbers/FIGURE_INDEX.md for what each shows and which script makes it)
  05_3d_models      STL / STEP / DXF of every reconstruction, cleaned real sherds, how to open them
  06_methods_documents  the Word note for the Sample 2 pipeline (the paper's methods are in the LaTeX manuscript)
  07_manuscript     main.tex (abstract, methods and results written; other sections keep the first draft's writing guides), sections/, references.bib, figures/, main.pdf,
                    TABLE2_AND_ABSTRACT_FIXES.md, TO_VERIFY.md

The five specimens (scan label -> what it turned out to be)
  Sample 1  label 2    a cluster: three ceramic pieces plus a brown ceramic piece that is a pedestal-footed dish (foot about 94 mm); the dish is reconstructed
  Sample 2  label 424  a rim-to-shoulder fragment of a bowl with a finished lip (rim diameter about 228 mm)
  Sample 3  label 62   a ring-footed base (standing circle about 107 mm), its interior a sphere of radius about 120 mm
  Sample 4  label 440  an outer wall band, 112 mm tall (barrel-shaped, diameter 190-218 mm)
  Sample 5  label 491  a stone (not reconstructed)

Read these three files before quoting any number
  03_numbers/FINAL_NUMBERS.md     every number with its uncertainty and status (MEASURED / ESTIMATE / ASSUMED)
  03_numbers/ASSUMPTIONS.md       everything that is a choice, not a measurement (all walls above the broken edges, wall thicknesses, the base centre of Sample 2)
  03_numbers/CORRECTIONS_LOG.md   errors found and fixed along the way: older numbers or files from earlier messages must not be used

Colour convention in the figures: clay = real scanned sherd; blue = symmetric completion or assumed wall; green = further projection; grey = not measured.
Prepared: 19 September 2026. Contact: [add name and email].

Manuscript note: table and figure numbers of the LaTeX manuscript follow order of appearance, so the INAA table is Table 2 (the specimen inventory is Table 1);
the first draft's morphometric table is compared in Table 5. Items the preparer could not verify are printed in red in main.pdf and listed in 07_manuscript/TO_VERIFY.md.
Sections 1, 2, 5 and 6 of the manuscript still hold the first draft's writing guides; the abstract, methods, results and the data/code statement were written from this work.
