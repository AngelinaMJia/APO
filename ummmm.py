import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import interp1d

# 1. Define Common Time Grid
t_common = np.arange(0, 120.5, 0.5)

# 2. Define Labels
file_labels = {
    "Book3.xlsx - LAS 0.1%.csv": "LAS 0.1%",
    "Book3.xlsx - SAF 0.1%.csv": "SAF 0.1%",
    "Book3.xlsx - SDS 0.1%.csv": "SDS 0.1%",
    "Book3.xlsx - LAS 3x CMC.csv": "LAS 3x CMC",
    "Book3.xlsx - SAF 3x CMC.csv": "SAF 3x CMC",
    "Book3.xlsx - SDS 3x CMC.csv": "SDS 3x CMC",
    "Book3.xlsx - DI water.csv": "DI Water"
}

# 3. Define Robust Processing Function
def process_surfactant_data_robust(file_name, common_time_grid):
    try:
        df = pd.read_csv(file_name)
    except:
        return None, None

    time_col = 'Elapsed time [s]'
    ca_col = 'CA(m) [°]'
    
    if time_col not in df.columns or ca_col not in df.columns:
        return None, None

    df = df.dropna(subset=[time_col, ca_col]).reset_index(drop=True)
    
    # Identify separate runs based on time reset
    elapsed = df[time_col].values
    diffs = np.diff(elapsed)
    split_indices = np.where(diffs < -5.0)[0] + 1
    split_indices = np.concatenate(([0], split_indices, [len(df)]))
    
    sets_interpolated = []
    
    for i in range(len(split_indices) - 1):
        start = split_indices[i]
        end = split_indices[i+1]
        subset = df.iloc[start:end]
        if len(subset) < 5: continue
            
        t_vals = subset[time_col].values
        ca_vals = subset[ca_col].values
        
        # Clean duplicates and sort
        _, unique_indices = np.unique(t_vals, return_index=True)
        t_vals = t_vals[unique_indices]
        ca_vals = ca_vals[unique_indices]
        if len(t_vals) < 2: continue

        f = interp1d(t_vals, ca_vals, kind='linear', bounds_error=False, fill_value=np.nan)
        sets_interpolated.append(f(common_time_grid))
        
    if not sets_interpolated:
        return None, None
        
    data_matrix = np.vstack(sets_interpolated)
    
    # --- ROBUST FILTERING ---
    # 1. Median of runs at each time point
    median_curve = np.nanmedian(data_matrix, axis=0)
    
    # 2. Deviation
    abs_diffs = np.abs(data_matrix - median_curve)
    
    # 3. MAD
    mad_curve = np.nanmedian(abs_diffs, axis=0)
    
    # 4. Threshold: max(3*MAD, 8.0 deg)
    threshold_curve = np.maximum(3 * mad_curve, 8.0)
    
    # 5. Mask outliers
    mask = abs_diffs <= threshold_curve
    cleaned_matrix = np.where(mask, data_matrix, np.nan)
    
    # Compute stats
    with np.errstate(invalid='ignore'):
        mean_curve = np.nanmean(cleaned_matrix, axis=0)
        std_curve = np.nanstd(cleaned_matrix, axis=0)
    
    return mean_curve, std_curve

# 4. Loop, Process, and Generate INDIVIDUAL Plots
generated_files = []

for fname, label in file_labels.items():
    mean_curve, std_curve = process_surfactant_data_robust(fname, t_common)
    
    if mean_curve is not None:
        plt.figure(figsize=(8, 5))
        
        mask = ~np.isnan(mean_curve)
        if np.any(mask):
            # Plot
            plt.plot(t_common[mask], mean_curve[mask], label=f'{label} Mean', color='blue', linewidth=2)
            plt.fill_between(t_common[mask], 
                             mean_curve[mask] - std_curve[mask], 
                             mean_curve[mask] + std_curve[mask], 
                             color='blue', alpha=0.2, label='Std Dev')
            
            plt.title(f"Average Contact Angle: {label}")
            plt.xlabel("Time [s]")
            plt.ylabel("Contact Angle [°]")
            plt.legend()
            plt.grid(True, alpha=0.3)
            
            # Set consistent limits if desired, or auto
            plt.xlim(0, 120)
            
            # Save
            safe_label = label.replace(" ", "_").replace("%", "pct").replace(".", "pt")
            out_file = f"graph_{safe_label}.png"
            plt.savefig(out_file)
            generated_files.append(out_file)
            plt.close() # Important to close so figures don't pile up in memory/display

print("Generated individual graphs:", generated_files)