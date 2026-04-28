import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import os
import re
from scipy import stats

def clean_value(val):
    if pd.isna(val):
        return 0.0
    if isinstance(val, (int, float)):
        return float(val)
    # Extract only the number before the ± sign
    match = re.search(r"([\d\.]+)", str(val))
    if match:
        return float(match.group(1))
    return 0.0

def analyze_data():
    base_dir = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
    file_path = os.path.join(base_dir, 'data', 'Project@Archeology.xlsx')
    
    # 1. Load the data
    try:
        df = pd.read_excel(file_path)
    except Exception as e:
        print(f"Error loading Excel file: {e}")
        return

    print("Columns found:", df.columns.tolist())
    
    # Precise element matching
    target_elements = ['Al', 'Fe', 'Na', 'Ti', 'Ca', 'K']
    samples = ['Sample 1', 'Sample 2', 'Sample 3', 'Sample 4', 'Sample 5']
    
    # 2. Extract and Clean Data
    # Find rows where the Element column starts with any of our targets
    extracted_rows = []
    for el in target_elements:
        row = df[df['Element'].str.startswith(el, na=False)]
        if not row.empty:
            extracted_rows.append(row)
    
    if not extracted_rows:
        print("No elements found! Check the 'Element' column content.")
        return
        
    filtered_df = pd.concat(extracted_rows).copy()
    
    for sample in samples:
        filtered_df[sample] = filtered_df[sample].apply(clean_value)
    
    # Clean the element names for legend/axis
    filtered_df['Element_Clean'] = filtered_df['Element'].apply(lambda x: x.split('(')[0].strip())
    filtered_df.set_index('Element_Clean', inplace=True)
    plot_data = filtered_df[samples]
    
    if plot_data.empty:
        print("Plot data is empty after filtering!")
        return

    # 3. Bar Chart (Comparison)
    plt.figure(figsize=(14, 8))
    plot_data.plot(kind='bar', ax=plt.gca())
    plt.title('Elementary Comparison: NAA Results for 5 Samples')
    plt.xlabel('Elements')
    plt.ylabel('Intensity')
    plt.legend(title='Samples')
    plt.xticks(rotation=0)
    plt.grid(axis='y', linestyle='--', alpha=0.7)
    plt.tight_layout()
    plt.savefig(os.path.join(base_dir, 'figures', 'element_intensities_histogram.png'))
    plt.close()

    # 4. Correlation Plots (Aluminum Focus)
    corr_df = plot_data.transpose()
    
    # NEW: Inter-sample correlation
    # We correlate the entire elemental profile of one sample against another
    sample_correlation = plot_data.corr()
    print("\nInter-sample Correlation Matrix (How similar are the samples?):")
    print(sample_correlation)
    
    plt.figure(figsize=(10, 8))
    sns.heatmap(sample_correlation, annot=True, cmap='YlGnBu', vmin=0.9, vmax=1.0)
    plt.title('Inter-sample Geochemical Correlation (Elemental Fingerprint)')
    plt.savefig(os.path.join(base_dir, 'figures', 'inter_sample_correlation.png'))
    plt.close()

    # Create Al vs K Correlation
    if 'Al' in corr_df.columns and 'K' in corr_df.columns:
        plt.figure(figsize=(10, 6))
        sns.regplot(x='Al', y='K', data=corr_df, ci=None, line_kws={"color": "red"})
        plt.title('Geochemical Correlation: Aluminum (Al) vs Potassium (K)')
        plt.xlabel('Al (Aluminum) Intensity')
        plt.ylabel('K (Potassium) Intensity')
        plt.grid(True, alpha=0.3)
        plt.savefig(os.path.join(base_dir, 'figures', 'correlation_al_k.png'))
        plt.close()

    # Create Al vs Fe Correlation
    if 'Al' in corr_df.columns and 'Fe' in corr_df.columns:
        plt.figure(figsize=(10, 6))
        sns.regplot(x='Al', y='Fe', data=corr_df, ci=None, line_kws={"color": "red"})
        plt.title('Geochemical Correlation: Aluminum (Al) vs Iron (Fe)')
        plt.xlabel('Al (Aluminum) Intensity')
        plt.ylabel('Fe (Iron) Intensity')
        plt.grid(True, alpha=0.3)
        plt.savefig(os.path.join(base_dir, 'figures', 'correlation_al_fe.png'))
        plt.close()

    print("Analysis Complete. Generated graphs including Inter-sample Fingerprint.")

if __name__ == "__main__":
    analyze_data()
