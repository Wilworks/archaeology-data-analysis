from pptx import Presentation
from pptx.util import Inches, Pt
from pptx.enum.text import PP_ALIGN
from pptx.dml.color import RGBColor
import pandas as pd
import re
import os

def create_presentation():
    base_dir = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
    prs = Presentation()
    
    # slide layout 0 is title, 1 is content, 5 is blank, 6 is title only
    
    # --- Slide 1: Title ---
    slide_layout = prs.slide_layouts[0]
    slide = prs.slides.add_slide(slide_layout)
    title = slide.shapes.title
    subtitle = slide.placeholders[1]
    
    title.text = "Archeological Technical Report: Fort Amsterdam NAA Analysis"
    subtitle.text = "Elemental Analysis of Ceramic Specimens\nJanuary 2026"

    # --- Slide 2: Executive Summary ---
    slide_layout = prs.slide_layouts[1]
    slide = prs.slides.add_slide(slide_layout)
    title = slide.shapes.title
    title.text = "Executive Summary"
    
    content = slide.placeholders[1]
    text_frame = content.text_frame
    text_frame.text = "Context: ~500 Years Post-Deposition (17th Century Colonial)"
    
    p = text_frame.add_paragraph()
    p.text = "Findings:"
    p.level = 0
    
    bullets = [
        "Specimens are confirmed Ceramic Pottery, not biological remains.",
        "High Al, Fe, K signatures indicate terrestrial clay matrix.",
        "High Ca levels consistent with shell tempering and coastal burial.",
        "Sodium (Na) infusion confirms long-term proximity to sea salt."
    ]
    
    for bullet in bullets:
        p = text_frame.add_paragraph()
        p.text = bullet
        p.level = 1

    # --- NEW Slide 3: Methodology ---
    slide_layout = prs.slide_layouts[1]
    slide = prs.slides.add_slide(slide_layout)
    title = slide.shapes.title
    title.text = "Methodology"
    
    content = slide.placeholders[1]
    text_frame = content.text_frame
    text_frame.text = "Neutron Activation Analysis (NAA):"
    
    p = text_frame.add_paragraph()
    p.text = "Non-destructive technique measuring elemental composition via neutron bombardment."
    p.level = 1
    
    p = text_frame.add_paragraph()
    p.text = "Provides precise quantification of trace elements (ppm-level accuracy)."
    p.level = 1
    
    p = text_frame.add_paragraph()
    p.text = ""
    
    p = text_frame.add_paragraph()
    p.text = "Computational Analysis Workflow:"
    p.level = 0
    
    p = text_frame.add_paragraph()
    p.text = "Data extraction: Parsed NAA intensity values from Excel, removing uncertainty ranges (±)."
    p.level = 1
    
    p = text_frame.add_paragraph()
    p.text = "Visualization: Generated comparative histograms using Python (matplotlib/seaborn)."
    p.level = 1
    
    p = text_frame.add_paragraph()
    p.text = "Statistical correlation: Calculated Pearson correlation coefficients between sample elemental profiles to identify geochemical consistency."
    p.level = 1

    # --- Slide 4: Historical Context ---
    slide_layout = prs.slide_layouts[1]
    slide = prs.slides.add_slide(slide_layout)
    title = slide.shapes.title
    title.text = "Historical Context & Speculation"
    
    content = slide.placeholders[1]
    text_frame = content.text_frame
    text_frame.text = "Original Speculation (Slice Analysis):"
    
    p = text_frame.add_paragraph()
    p.text = "High Al, Ca, K, and Na were initially attributed to human remains (bones, coffins, grave goods)."
    p.level = 1
    
    p = text_frame.add_paragraph()
    p.text = "Heavy metals (Ti, Fe) suggested agricultural/industrial contamination."
    p.level = 1
    
    p = text_frame.add_paragraph()
    p.text = ""
    
    p = text_frame.add_paragraph()
    p.text = "NAA Revelation:"
    p.level = 0
    
    p = text_frame.add_paragraph()
    p.text = "3D rendering confirmed pottery shapes, not biological material."
    p.level = 1
    
    p = text_frame.add_paragraph()
    p.text = "Chemical profile matches coastal ceramic production, not burial remains."
    p.level = 1

    # --- Slide 4: Element Intensities ---
    slide_layout = prs.slide_layouts[5] # Blank
    slide = prs.slides.add_slide(slide_layout)
    
    # Add title manually since blank layout
    title_box = slide.shapes.add_textbox(Inches(0.5), Inches(0.2), Inches(9), Inches(1))
    tf = title_box.text_frame
    tf.text = "Chemical Profile: Element Intensities"
    tf.paragraphs[0].font.size = Pt(40)
    tf.paragraphs[0].font.bold = True
    
    img_path = os.path.join(base_dir, 'figures', 'element_intensities_histogram.png')
    try:
        slide.shapes.add_picture(img_path, Inches(1), Inches(1.5), width=Inches(8))
    except FileNotFoundError:
        print(f"Warning: {img_path} not found. Skipping image.")

    textbox = slide.shapes.add_textbox(Inches(0.5), Inches(6.5), Inches(9), Inches(1))
    tf = textbox.text_frame
    tf.text = "Consistent Al/Fe/K/Ca profile across samples indicates shared provenance."

    # --- Slide 6: Data Table ---
    slide_layout = prs.slide_layouts[5] # Blank
    slide = prs.slides.add_slide(slide_layout)
    
    title_box = slide.shapes.add_textbox(Inches(0.5), Inches(0.2), Inches(9), Inches(0.8))
    tf = title_box.text_frame
    tf.text = "Intensity Data (NAA Results)"
    tf.paragraphs[0].font.size = Pt(36)
    tf.paragraphs[0].font.bold = True
    
    # Read the data
    try:
        df = pd.read_excel(os.path.join(base_dir, 'data', 'Project@Archeology.xlsx'))
        target_elements = ['Al', 'Fe', 'Na', 'Ti', 'Ca', 'K']
        samples = ['Sample 1', 'Sample 2', 'Sample 3', 'Sample 4', 'Sample 5']
        
        extracted_rows = []
        for el in target_elements:
            row = df[df['Element'].str.startswith(el, na=False)]
            if not row.empty:
                extracted_rows.append(row)
        
        if extracted_rows:
            filtered_df = pd.concat(extracted_rows)
            
            # Create table
            rows = len(target_elements) + 1  # +1 for header
            cols = 6  # Element + 5 samples
            
            left = Inches(0.5)
            top = Inches(1.2)
            width = Inches(9)
            height = Inches(4.5)
            
            table = slide.shapes.add_table(rows, cols, left, top, width, height).table
            
            # Set column widths
            table.columns[0].width = Inches(1.8)
            for i in range(1, 6):
                table.columns[i].width = Inches(1.4)
            
            # Header row
            table.cell(0, 0).text = "Element"
            for i, sample in enumerate(samples):
                table.cell(0, i+1).text = sample
            
            # Data rows
            for i, el in enumerate(target_elements):
                row_data = filtered_df[filtered_df['Element'].str.startswith(el, na=False)]
                if not row_data.empty:
                    element_name = row_data.iloc[0]['Element'].split('(')[0].strip()
                    table.cell(i+1, 0).text = element_name
                    
                    for j, sample in enumerate(samples):
                        val = row_data.iloc[0][sample]
                        # Clean the value
                        if pd.notna(val):
                            match = re.search(r'([\d\.]+)', str(val))
                            if match:
                                table.cell(i+1, j+1).text = match.group(1)
            
            # Format table
            for row in table.rows:
                for cell in row.cells:
                    cell.text_frame.paragraphs[0].font.size = Pt(14)
                    cell.text_frame.paragraphs[0].alignment = PP_ALIGN.CENTER
                    
    except Exception as e:
        print(f"Warning: Could not create data table: {e}")
        textbox = slide.shapes.add_textbox(Inches(1), Inches(2), Inches(8), Inches(2))
        tf = textbox.text_frame
        tf.text = "Data table could not be generated. Please refer to the Excel file."

    # --- Slide 7: Correlation Analysis ---
    slide_layout = prs.slide_layouts[5] # Blank
    slide = prs.slides.add_slide(slide_layout)
    
    title_box = slide.shapes.add_textbox(Inches(0.5), Inches(0.2), Inches(9), Inches(1))
    tf = title_box.text_frame
    tf.text = "Geochemical Correlation Matrix"
    tf.paragraphs[0].font.size = Pt(40)
    tf.paragraphs[0].font.bold = True
    
    img_path = os.path.join(base_dir, 'figures', 'inter_sample_correlation.png')
    try:
        # Scale to fit
        slide.shapes.add_picture(img_path, Inches(2), Inches(1.5), height=Inches(4.5))
    except FileNotFoundError:
        print(f"Warning: {img_path} not found. Skipping image.")

    # Add text explaining Sample 2
    textbox = slide.shapes.add_textbox(Inches(0.5), Inches(6.2), Inches(9), Inches(1))
    tf = textbox.text_frame
    tf.word_wrap = True
    p = tf.add_paragraph()
    p.text = "Sample 2 is a significant outlier (Correlation 0.22 - 0.45)."
    p.font.bold = True
    p.font.color.rgb = RGBColor(255, 0, 0) # Red
    
    p2 = tf.add_paragraph()
    p2.text = "This suggests it originates from a different vessel, workshop, or has undergone isolated contamination compared to the main group (1, 3, 4, 5)."

    # --- Slide 8: Conclusion ---
    slide_layout = prs.slide_layouts[1]
    slide = prs.slides.add_slide(slide_layout)
    title = slide.shapes.title
    title.text = "Conclusion & Archaeological Implications"
    
    content = slide.placeholders[1]
    text_frame = content.text_frame
    
    conclusions = [
        "Interpretation: The artifacts are a single set of colonial-era utility ware (Main Group), with one anomalous piece (Sample 2).",
        "Coastal Context: Validated by ubiquitous Sodium presence.",
        "Material: Local or imported volcanic-clay pottery (High Fe/Al) with shell temper (Ca).",
        "Historical Significance: Reflects Fort Amsterdam's role as a 17th-century hub with diverse pottery sources (local production + imports).",
        "Result: Speculations of human remains are chemically refuted."
    ]
    
    for line in conclusions:
        p = text_frame.add_paragraph()
        p.text = line
        p.space_after = Pt(12)

    output_path = os.path.join(base_dir, 'reports', 'Fort_Amsterdam_Analysis.pptx')
    prs.save(output_path)
    print(f"Presentation saved as {output_path}")

if __name__ == "__main__":
    create_presentation()
