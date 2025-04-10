#!/usr/bin/env python3

import os
import glob
import re
import json
import matplotlib.pyplot as plt
import sys
import csv

# --- Configuration ---
# IMPORTANT: Replace this with the actual top-level path on your system
# 
# TOP_LEVEL_DIRECTORY = '/Users/fool/Documents/mint/condor_files/hexme_remesh_i_nocurl'
TOP_LEVEL_DIRECTORY = '/Users/fool/Documents/mint/condor_files/solids_overflow'


# This pattern assumes the structure is tld/simulation_dir/inner_iters/state_step
STATE_STEP_GLOB_PATTERN = '*/inner_iters/state_step'
FILE_PATTERN = 'conf_1000*.json'
OUTPUT_LOG_FILE = 'multi_dir_output_log.txt'
OUTPUT_PLOT_FILE = 'multi_dir_energy_plot.png'
# OUTPUT_CSV_FILE = 'hexme_remesh_i_nocurl_output_data.csv'
OUTPUT_CSV_FILE = 'solids_overflow_output_data.csv'

# --- End Configuration ---

def find_largest_id_file(directory, pattern):
    """
    Finds the file matching the pattern with the largest numeric ID in its name.
    (Same as before)
    """
    full_pattern = os.path.join(directory, pattern)
    json_files = glob.glob(full_pattern)

    largest_num = -1
    largest_file_path = None
    largest_filename = None

    for file_path in json_files:
        filename = os.path.basename(file_path)
        match = re.search(r'conf_1000(\d{2})\.json', filename)

        if match:
            try:
                num = int(match.group(1))
                if num > largest_num:
                    largest_num = num
                    largest_file_path = file_path
                    largest_filename = filename
            except ValueError:
                print(f"Warning: Invalid number format in {filename} in dir {directory}", file=sys.stderr)
                continue
        # else: # Optional warning
            # print(f"Warning: Filename {filename} in dir {directory} does not match expected pattern.", file=sys.stderr)

    return largest_file_path, largest_filename, largest_num

def get_mesh_name_from_parent_dir(step_dir, log_file_handle):
    """
    Finds the mesh name by looking for a .mesh file in ../.. relative to step_dir.
    """
    mesh_name = "UNKNOWN_MESH" # Default value
    try:
        simulation_dir_path = os.path.dirname(os.path.dirname(step_dir))
        if not os.path.isdir(simulation_dir_path):
             warning_msg = f"Simulation directory '{simulation_dir_path}' not found relative to '{step_dir}'."
             print(f"  Warning: {warning_msg}", file=sys.stderr)
             log_file_handle.write(f"  [WARN] {warning_msg} Cannot find mesh file.\n")
             return mesh_name # Return default early

        mesh_files = glob.glob(os.path.join(simulation_dir_path, '*.mesh'))

        if len(mesh_files) == 1:
            mesh_filename = os.path.basename(mesh_files[0])
            mesh_name = os.path.splitext(mesh_filename)[0]
            # print(f"    Debug: Found mesh file: {mesh_files[0]}, extracted name: {mesh_name}") # Optional debug
        elif len(mesh_files) == 0:
            warning_msg = f"No '.mesh' file found in directory '{simulation_dir_path}'."
            print(f"  Warning: {warning_msg}", file=sys.stderr)
            log_file_handle.write(f"  [WARN] {warning_msg} Using default '{mesh_name}'.\n")
        else: # More than one .mesh file
            warning_msg = f"Multiple '.mesh' files found in '{simulation_dir_path}': {mesh_files}."
            print(f"  Warning: {warning_msg}", file=sys.stderr)
            log_file_handle.write(f"  [WARN] {warning_msg} Ambiguous mesh name. Using default '{mesh_name}'.\n")

    except Exception as e:
        # Catch potential errors during path manipulation/glob
        warning_msg = f"Error finding mesh file relative to '{step_dir}': {e}"
        print(f"  Warning: {warning_msg}", file=sys.stderr)
        log_file_handle.write(f"  [WARN] {warning_msg} Using default '{mesh_name}'.\n")

    return mesh_name


def process_simulation_dirs(tld, state_step_pattern, file_pattern_in_step, log_file, csv_file):
    """
    Processes the largest ID file in each state_step directory, logs info,
    finds mesh name from ../../ *.mesh file, and saves data to CSV.
    """
    state_step_dirs = glob.glob(os.path.join(tld, state_step_pattern))

    energy_primal_list = []
    calculated_value_list = []

    print(f"Scanning TLD: {tld}")
    print(f"Found {len(state_step_dirs)} potential 'state_step' directories.")

    if not os.path.exists(tld):
        print(f"ERROR: Top Level Directory not found: {tld}", file=sys.stderr)
        return [], []
    if not state_step_dirs:
         print(f"Warning: No directories matching '{state_step_pattern}' found under {tld}.", file=sys.stderr)

    # --- Prepare Log and CSV files ---
    try:
        # Ensure we can open both files before starting the main loop
        with open(log_file, 'w') as logfile, \
             open(csv_file, 'w', newline='') as csvfile:

            logfile.write(f"Log for processing largest ID file in state_step dirs under: {tld}\n")
            logfile.write("Using mesh name from '../../ *.mesh' file.\n")
            logfile.write("="*40 + "\n")

            csv_writer = csv.writer(csvfile)
            header = ['MeshName', 'EnergyPrimalIntegrability', 'CalculatedValue']
            csv_writer.writerow(header)
            print(f"Opened Log file: {log_file}")
            print(f"Opened CSV file for writing: {csv_file}")

            # Process each state_step directory found
            for step_dir in sorted(state_step_dirs):
                # --- Extract Mesh Name using the new function ---
                mesh_name = get_mesh_name_from_parent_dir(step_dir, logfile)
                # --- End Extract Mesh Name ---

                print(f"\nProcessing Directory: {step_dir} (Mesh: {mesh_name})")
                logfile.write(f"\n--- Directory: {step_dir} (Mesh: {mesh_name}) ---\n")

                largest_file_path, largest_filename, largest_num = find_largest_id_file(step_dir, file_pattern_in_step)

                if largest_file_path is None:
                    print("  No valid 'conf_1000XX.json' files found.")
                    logfile.write("  No valid 'conf_1000XX.json' files found.\n")
                    continue

                print(f"  Largest ID file found: {largest_filename} (ID={largest_num})")

                # --- Apply conditional logic ---
                if largest_num < 50:
                    logfile.write(f"  [ID < 50] Logged filename: {largest_filename}\n")
                    print(f"  Logged filename (ID={largest_num} < 50)")
                else:
                    # Parse JSON for ID >= 50
                    try:
                        with open(largest_file_path, 'r') as f:
                            data = json.load(f)

                        energy_primal = data.get('energy_primal_integrability')
                        rescale_param = data.get('w_global_rescale_param')
                        energy_mint = data.get('energy_mint')

                        # Check validity (same as before)
                        if energy_primal is not None and rescale_param is not None and energy_mint is not None and \
                           isinstance(energy_primal, (int, float)) and \
                           isinstance(rescale_param, (int, float)) and \
                           isinstance(energy_mint, (int, float)):

                            calculated_value = energy_mint * rescale_param

                            # Append data for plotting
                            energy_primal_list.append(energy_primal)
                            calculated_value_list.append(calculated_value)

                            # Write data row including mesh name to CSV file
                            csv_writer.writerow([mesh_name, energy_primal, calculated_value])

                            # Log extracted and calculated values
                            logfile.write(f"  [ID >= 50] File: {largest_filename}, "
                                          f"Mesh: {mesh_name}, "
                                          f"EnergyPrimal: {energy_primal:.6e}, "
                                          f"CalculatedValue: {calculated_value:.6e}, "
                                          f"Data written to CSV.\n")
                            print(f"  Processed data (ID={largest_num} >= 50), wrote to CSV.")
                        else:
                            # Handle missing/invalid data logging (same as before)
                            missing_keys = [k for k in ['energy_primal_integrability', 'w_global_rescale_param', 'energy_mint'] if data.get(k) is None]
                            non_numeric = not (isinstance(energy_primal, (int, float)) and isinstance(rescale_param, (int, float)) and isinstance(energy_mint, (int, float)))
                            warn_msg = f"  Warning: Skipped data point from {largest_filename}."
                            log_warn_msg = f"  [WARN] Skipped data point from {largest_filename} -"
                            if missing_keys: warn_msg += f" Missing keys: {missing_keys}."; log_warn_msg += f" Missing keys ({', '.join(missing_keys)})."
                            if non_numeric: warn_msg += f" Non-numeric value found."; log_warn_msg += f" Non-numeric value found."
                            warn_msg += " Not written to CSV."; log_warn_msg += " Not written to CSV."
                            print(warn_msg, file=sys.stderr); logfile.write(log_warn_msg + "\n")

                    # Exception handling for JSON parsing etc. (same as before)
                    except FileNotFoundError:
                        print(f"  Error: File not found during processing: {largest_file_path}", file=sys.stderr); logfile.write(f"  [ERROR] File not found: {largest_filename}\n")
                    except json.JSONDecodeError:
                        print(f"  Error: Could not decode JSON from file: {largest_file_path}", file=sys.stderr); logfile.write(f"  [ERROR] Invalid JSON: {largest_filename}\n")
                    except Exception as e:
                        print(f"  An unexpected error occurred processing {largest_file_path}: {e}", file=sys.stderr); logfile.write(f"  [ERROR] Unexpected error ({type(e).__name__}) processing: {largest_filename}\n")

    except IOError as e:
         # Handle errors opening log or CSV file initially
         print(f"ERROR: Could not open/write to output files ({log_file}, {csv_file}): {e}", file=sys.stderr)
         return [], [] # Return empty lists if setup fails

    print(f"\nProcessing complete. Log saved to: {log_file}, Data saved to: {csv_file}")
    return energy_primal_list, calculated_value_list

def plot_data(x_data, y_data, plot_file):
    """
    Generates and saves a plot of y_data vs x_data, sorted by x_data.
    (Same as before)
    """
    if not x_data or not y_data or len(x_data) != len(y_data):
        print("\nNo valid data pairs collected across directories for plotting.")
        return

    print(f"\nGenerating plot with {len(x_data)} data points...")
    # ... (rest of plotting code is identical) ...
    paired_data = list(zip(x_data, y_data))
    paired_data.sort(key=lambda pair: pair[0])
    sorted_x, sorted_y = zip(*paired_data)
    plt.figure(figsize=(12, 7))
    plt.plot(sorted_x, sorted_y, marker='o', linestyle='-', markersize=4, label='Calculated Value')
    plt.xlabel("Energy Primal Integrability (Sorted)")
    plt.ylabel("Energy Mint * W Global Rescale Param")
    plt.title("Calculated Value vs Sorted Energy Primal (Largest ID file per simulation)")
    plt.grid(True); plt.legend(); plt.tight_layout()
    try:
        plt.savefig(plot_file); print(f"Plot saved successfully to: {plot_file}")
    except Exception as e: print(f"Error saving plot to {plot_file}: {e}", file=sys.stderr)
    # plt.show()

# --- Main Execution ---
if __name__ == "__main__":
    # Process directories, using the new mesh name logic, log, save CSV
    primal_energies, calc_values = process_simulation_dirs(
        TOP_LEVEL_DIRECTORY,
        STATE_STEP_GLOB_PATTERN,
        FILE_PATTERN,
        OUTPUT_LOG_FILE,
        OUTPUT_CSV_FILE
    )

    # Plot the numeric data collected
    plot_data(primal_energies, calc_values, OUTPUT_PLOT_FILE)