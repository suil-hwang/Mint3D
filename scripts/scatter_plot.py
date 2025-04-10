#!/usr/bin/env python3

import csv
import matplotlib.pyplot as plt
import sys
import os
import numpy as np # Needed for checking finite values

# --- Configuration ---
# Input merged CSV file
MERGED_CSV_FILE = 'merged_sorted_data.csv'
# MERGED_CSV_FILE = 'merged_solids_sorted_data.csv'



# Output plot file
OUTPUT_PLOT_FILE = 'merged_log_scatter_styled_plot.png'

# Column indices (0-based) in the merged CSV
# Assuming: MeshName, EnergyF1, CalcValF1, EnergyF2, CalcValF2
MESH_NAME_IDX = 0
ENERGY_F1_IDX = 1
CALCVAL_F1_IDX = 2
ENERGY_F2_IDX = 3
CALCVAL_F2_IDX = 4

# --- Styling Configuration ---
# Colors similar to example (adjust as needed)
COLOR_F1 = '#2ca02c' # Green
COLOR_F2 = '#5e3c99' # Dark Purple/Blue

# Plot area background color
# AXES_FACE_COLOR = '#f0f0f0' # Light gray
AXES_FACE_COLOR = '#ececec'


# Marker settings
MARKER_SIZE = 75
MARKER_ALPHA = 0.7
# --- End Configuration ---

def plot_styled_scatter(csv_filepath, plot_filepath):
    """
    Reads the merged CSV and creates a styled log-log scatter plot comparing
    (CalcVal vs Energy) for File 1 and File 2 data.

    Args:
        csv_filepath (str): Path to the merged input CSV file.
        plot_filepath (str): Path to save the output plot image.
    """
    energy_f1 = []
    calcval_f1 = []
    energy_f2 = []
    calcval_f2 = []
    # mesh_names = [] # Not strictly needed for scatter plot

    required_indices = [ENERGY_F1_IDX, CALCVAL_F1_IDX, ENERGY_F2_IDX, CALCVAL_F2_IDX]
    max_idx = max(required_indices)

    print(f"Reading merged data from: {csv_filepath}")
    try:
        with open(csv_filepath, 'r', newline='') as infile:
            reader = csv.reader(infile)
            try:
                header = next(reader) # Read header
                if len(header) <= max_idx:
                     print(f"ERROR: CSV file '{csv_filepath}' has fewer columns than expected.", file=sys.stderr)
                     return False
                print(f"CSV Headers: {header}")
            except StopIteration:
                print(f"ERROR: CSV file '{csv_filepath}' is empty or has no header.", file=sys.stderr)
                return False

            rows_read = 0
            rows_plotted = 0
            for i, row in enumerate(reader):
                rows_read += 1
                if len(row) <= max_idx:
                    print(f"Warning: Skipping row {i+2} due to insufficient columns.", file=sys.stderr)
                    continue

                try:
                    # Extract data using configured indices
                    e1 = float(row[ENERGY_F1_IDX])
                    c1 = float(row[CALCVAL_F1_IDX])
                    e2 = float(row[ENERGY_F2_IDX])
                    c2 = float(row[CALCVAL_F2_IDX])

                    # Log plots require positive, finite values
                    if np.isfinite(e1) and e1 > 0 and \
                       np.isfinite(c1) and c1 > 0 and \
                       np.isfinite(e2) and e2 > 0 and \
                       np.isfinite(c2) and c2 > 0:
                        energy_f1.append(e1)
                        calcval_f1.append(c1)
                        energy_f2.append(e2)
                        calcval_f2.append(c2)
                        # mesh_names.append(row[MESH_NAME_IDX])
                        rows_plotted += 1
                    else:
                        mesh = row[MESH_NAME_IDX]
                        print(f"Warning: Skipping row {i+2} for mesh '{mesh}' due to non-positive or non-finite value(s) required for log plot.", file=sys.stderr)

                except (ValueError, TypeError):
                    mesh = row[MESH_NAME_IDX]
                    print(f"Warning: Skipping row {i+2} for mesh '{mesh}'. Could not convert one or more values to float.", file=sys.stderr)
                    continue
                except IndexError:
                     print(f"Warning: Skipping row {i+2} due to unexpected index error (check column indices).", file=sys.stderr)
                     continue

    except FileNotFoundError:
        print(f"ERROR: Input CSV file not found: {csv_filepath}", file=sys.stderr)
        return False
    except Exception as e:
        print(f"ERROR: Failed to read {csv_filepath}: {e}", file=sys.stderr)
        return False

    print(f"Read {rows_read} data rows. Prepared {rows_plotted} rows for plotting.")

    if not energy_f1: # Check if any valid data was collected
        print("No valid data points found for plotting.", file=sys.stderr)
        return False

    # --- Create and Style Plot ---
    print("Generating styled plot...")
    # Set overall font style (optional)
    # plt.rcParams['font.family'] = 'sans-serif'
    # plt.rcParams['font.sans-serif'] = ['DejaVu Sans'] # Or 'Arial', 'Helvetica' if available

    fig, ax = plt.subplots(figsize=(8, 7)) # Adjust figure size as needed

    # Set plot area background color
    ax.set_facecolor(AXES_FACE_COLOR)

    # --- Plot Data using Scatter ---
    # Plot File 1 data
    ax.scatter(energy_f1, calcval_f1,
               s=MARKER_SIZE,           # Marker size
               c=COLOR_F1,              # Face color
               edgecolor='grey',       # Optional: darker edge for definition
               linewidth=0.5,           # Optional: edge line width
               alpha=MARKER_ALPHA,      # Transparency
               label=f'{header[CALCVAL_F1_IDX]} vs {header[ENERGY_F1_IDX]}' # Legend label
              )

    # Plot File 2 data
    ax.scatter(energy_f2, calcval_f2,
               s=MARKER_SIZE,
               c=COLOR_F2,
               marker='^',              # Use a different marker (e.g., triangle)
               edgecolor='grey',
               linewidth=0.5,
               alpha=MARKER_ALPHA,
               label=f'{header[CALCVAL_F2_IDX]} vs {header[ENERGY_F2_IDX]}'
              )

    # --- Apply Styling ---
    # Set log scales
    ax.set_xscale('log')
    ax.set_yscale('log')

    # Grid configuration
    ax.grid(True, which='major', color='white', linestyle='-', linewidth=1.0)
    ax.grid(True, which='minor', color='white', linestyle='--', linewidth=0.5)
    ax.set_axisbelow(True) # Ensure grid is behind data

    # Tick configuration
    ax.tick_params(axis='both', which='both', direction='out', top=False, right=False, labelsize=14) # Ticks inward, also on top/right
    # plt.xticks(fontsize=14)
    # plt.yticks(fontsize=14)
    ax.tick_params(axis='both', which='major', length=20, width=3, color='grey')
    ax.tick_params(axis='both', which='minor', length=10, width=1.5)    

    # ax.patch.set_linewidth(10)
    ax.spines['top'].set_color("grey")
    ax.spines['bottom'].set_color("grey")
    ax.spines['left'].set_color("grey")
    ax.spines['right'].set_color("grey")
    # plt.axis('off') 
    # ax.spines['top'].set_visible(False)
    # ax.spines['right'].set_visible(False)

    # Add labels and title
    ax.set_xlabel("Q_int (Log Scale)")
    ax.set_ylabel("E_mint (Log Scale)")
    ax.set_title("Improvement of ")
    ax.set_aspect(.75)        

    # Add legend
    # ax.legend()

    # plt.tight_layout()
    # plt.figure(figsize=(5, 10))

    # --- Save Plot ---
    try:
        plt.savefig(plot_filepath, dpi=150, facecolor=fig.get_facecolor()) # Ensure figure background color is saved
        print(f"Plot saved successfully to: {plot_filepath}")
    except Exception as e:
        print(f"Error saving plot to {plot_filepath}: {e}", file=sys.stderr)
        return False

    # Optionally display the plot
    # plt.show()

    return True

# --- Main Execution ---
if __name__ == "__main__":
    # Verify input file exists
    if not os.path.exists(MERGED_CSV_FILE):
        print(f"ERROR: Merged CSV file not found: {MERGED_CSV_FILE}", file=sys.stderr)
        sys.exit(1)

    plot_styled_scatter(MERGED_CSV_FILE, OUTPUT_PLOT_FILE)