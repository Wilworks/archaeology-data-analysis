# Archaeological Data Analysis: Fort Amsterdam Pottery

**Author:** Wilfred

This repository contains the data, analysis scripts, and final reports for the Neutron Activation Analysis (NAA) performed on pottery samples found at Fort Amsterdam. The study aimed to determine the elemental composition of these samples, validating hypotheses about their coastal origin and investigating signs of environmental interaction over a 500-year burial period.

## Project Structure

```text
├── data/                  # Raw datasets (NAA intensity results)
│   ├── Project@Archeology.xlsx
│   └── sample_corr.csv
├── src/                   # Python source code for analysis and reporting
│   ├── analyze_archeology.py
│   └── generate_report_pptx.py
├── figures/               # Generated statistical and geochemical visualizations
│   └── *.png
├── reports/               # Formal documentation and presentations
│   ├── Fort_Amsterdam_Analysis.pptx
│   └── Archaeological_Report_Fort_Amsterdam.md
├── requirements.txt       # Python dependencies
└── README.md              # This file
```

## Setup and Installation

1. Clone the repository to your local machine.
2. Ensure you have Python installed.
3. Install the required dependencies:
   ```bash
   pip install -r requirements.txt
   ```

## Running the Analysis

To generate the analytical figures (saved to the `figures/` directory):
```bash
python src/analyze_archeology.py
```

To generate the PowerPoint presentation summarizing the results (saved to the `reports/` directory):
```bash
python src/generate_report_pptx.py
```

## Key Findings

- Specimens are confirmed Ceramic Pottery (not biological remains).
- Chemical signatures strongly indicate a terrestrial clay matrix with shell tempering, matching expected local coastal ceramic profiles.
- Elemental consistencies across samples indicate a shared provenance, with the exception of one notable outlier (Sample 2).

For detailed methodology and historical context, refer to the documents in the `reports/` directory.
