#!/usr/bin/env python3

import csv
import matplotlib.pyplot as plt
import sys
import os
import datetime # Import datetime to add timestamp to filename

# --- Configuration ---
# Input merged CSV file (should be sorted by the column at X_AXIS_IDX)
MERGED_CSV_FILE = 'merged_solids_sorted_data.csv'

# --- Define Column Indices based on the merged file header ---
# Example Header: MeshName, EnergyPrimalIntegrability_F1, CalculatedValue_F1, EnergyPrimalIntegrability_F2, CalculatedValue_F2
MESH_NAME_IDX = 0
ENERGY_F1_IDX = 1
CALCVAL_F1_IDX = 2 # This will be our X-axis
ENERGY_F2_IDX = 3
CALCVAL_F2_IDX = 4

# --- Define which columns to use for Axes ---
X_AXIS_IDX = CALCVAL_F1_IDX # Set X-axis to the column used for sorting in merge script

# Define columns to plot on Y-axis against the new X-axis
# We'll plot EnergyF1, CalcValF1 (same as X), EnergyF2, CalcValF2
# ADD CALCVAL_F1_IDX BACK TO THIS LIST
Y_AXIS_INDICES = [ENERGY_F1_IDX, CALCVAL_F1_IDX, ENERGY_F2_IDX, CALCVAL_F2_IDX]

# Output plot file (Added timestamp)
timestamp = datetime.datetime.now().strftime("%Y%m%d_%H%M%S")
OUTPUT_PLOT_FILE = f'merged_semi_log_plot_X{X_AXIS_IDX+1}_Y{[i+1 for i in Y_AXIS_INDICES]}_{timestamp}.png'


# --- End Configuration ---

def plot_semi_log_lines(csv_filepath, plot_filepath):
    """
    Reads the merged CSV, ensures data is sorted by the chosen X_AXIS_IDX,
    and creates a semi-log line plot of specified Y-axis streams against
    the chosen X-axis column.

    Args:
        csv_filepath (str): Path to the merged input CSV file.
        plot_filepath (str): Path to save the output plot image.
    """
    # Store data points temporarily before plotting
    # Structure: (x_val, [y1_val, y2_val, y3_val, y4_val]) <- Now 4 Y values
    data_points = []

    # Find highest index needed from all configured columns
    all_indices = Y_AXIS_INDICES + [X_AXIS_IDX, MESH_NAME_IDX]
    max_idx = max(all_indices)

    print(f"Reading merged data from: {csv_filepath}")
    try:
        with open(csv_filepath, 'r', newline='') as infile:
            reader = csv.reader(infile)
            try:
                header = next(reader) # Read header
                if len(header) <= max_idx:
                     print(f"ERROR: CSV file '{csv_filepath}' has fewer columns ({len(header)}) than expected based on column indices (max needed: {max_idx}).", file=sys.stderr)
                     print(f"Header found: {header}", file=sys.stderr)
                     return False
                print(f"CSV Headers: {header}")
                # Extract labels using the indices
                y_labels = [header[idx] for idx in Y_AXIS_INDICES]
                x_label = header[X_AXIS_IDX]
                print(f"Using X-Axis: '{x_label}' (Column {X_AXIS_IDX+1})")
                print(f"Using Y-Axis Streams: {[f'{lbl} (Col {idx+1})' for lbl, idx in zip(y_labels, Y_AXIS_INDICES)]}")

            except StopIteration:
                print(f"ERROR: CSV file '{csv_filepath}' is empty or has no header.", file=sys.stderr)
                return False

            rows_read = 0
            valid_rows_collected = 0
            for i, row in enumerate(reader):
                rows_read += 1
                if len(row) <= max_idx:
                    print(f"Warning: Skipping row {i+2} due to insufficient columns ({len(row)}). Expected at least {max_idx + 1}.", file=sys.stderr)
                    continue

                try:
                    # Extract X value (based on X_AXIS_IDX)
                    x_val = float(row[X_AXIS_IDX])

                    # Extract Y values (based on Y_AXIS_INDICES) - now extracts 4 values
                    y_vals = [float(row[idx]) for idx in Y_AXIS_INDICES]

                    # Crucial for log Y-axis: Ensure ALL Y values AND the X value are positive
                    # (since X value is also one of the Y values now)
                    if all(v > 0 for v in y_vals): # This implicitly checks x_val too now
                        # Store the x value and the list of corresponding y values
                        data_points.append((x_val, y_vals))
                        valid_rows_collected += 1
                    else:
                        # Check X value positivity explicitly for warning clarity if needed
                        x_pos = x_val > 0
                        y_pos_check = [y > 0 for y in y_vals]
                        mesh = row[MESH_NAME_IDX] # Get mesh name for warning
                        print(f"Warning: Skipping data point for mesh '{mesh}' (row {i+2}) due to non-positive value(s) required for log scale. [X positive: {x_pos}, Ys positive: {y_pos_check}, Y values: {y_vals}]", file=sys.stderr)


                except ValueError:
                    mesh = row[MESH_NAME_IDX]
                    print(f"Warning: Skipping data point for mesh '{mesh}' (row {i+2}). Could not convert one or more values to float.", file=sys.stderr)
                    continue
                except IndexError:
                     print(f"Warning: Skipping row {i+2} due to unexpected index error (check column indices configuration).", file=sys.stderr)
                     continue

    except FileNotFoundError:
        print(f"ERROR: Input CSV file not found: {csv_filepath}", file=sys.stderr)
        return False
    except Exception as e:
        print(f"ERROR: Failed to read {csv_filepath}: {e}", file=sys.stderr)
        return False

    print(f"Read {rows_read} data rows. Collected {valid_rows_collected} valid data points for plotting (with positive values).")

    if not data_points: # Check if any valid data was collected
        print("No valid data points found for plotting.", file=sys.stderr)
        return False

    # --- Sort the collected data based on the X value (robustness check) ---
    print(f"Sorting {len(data_points)} data points by X-axis value ('{x_label}')...")
    data_points.sort(key=lambda point: point[0]) # Sort by the first element (x_val)

    # --- Unpack the sorted data into separate lists for plotting ---
    x_axis_sorted = [p[0] for p in data_points]
    # Create lists for each Y stream, extracting from the nested list p[1]
    y_streams_sorted = []
    # The number of Y streams is now len(Y_AXIS_INDICES)
    for i in range(len(Y_AXIS_INDICES)):
        y_streams_sorted.append([p[1][i] for p in data_points])

    # --- Create Plot ---
    print("Generating plot...")
    try:
        # Try to use a modern style first
        plt.style.use('seaborn-v0_8-darkgrid')
    except:
        # Fallback for older matplotlib versions
        plt.style.use('seaborn-darkgrid')


    fig, ax = plt.subplots(figsize=(12, 7)) # Wider figure for line plots

    # Plot the Y data streams using the SORTED data
    # The loop now iterates over all indices in Y_AXIS_INDICES
    for i in range(len(y_labels)):
        ax.plot(x_axis_sorted, y_streams_sorted[i], label=y_labels[i], marker='.', markersize=4)

    # Set Y axis to log scale, X axis remains linear (default)
    ax.set_yscale('log')

    # Add labels and title
    ax.set_xlabel(f"{x_label} (Linear Scale)")
    ax.set_ylabel("Values (Log Scale)")
    ax.set_title(f"Data Streams vs {x_label} (Semi-Log Plot)")

    # Add legend and grid
    ax.legend()
    ax.grid(True, which='major', axis='both', linestyle='--', linewidth=0.5, alpha=0.7)
    ax.grid(True, which='minor', axis='y', linestyle=':', linewidth=0.3, alpha=0.5)

    plt.tight_layout()

    # --- Save Plot ---
    try:
        plt.savefig(plot_filepath, dpi=150)
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

    plot_semi_log_lines(MERGED_CSV_FILE, OUTPUT_PLOT_FILE)