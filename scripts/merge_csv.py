#!/usr/bin/env python3

import csv
import sys
import os

# --- Configuration ---
# Specify the paths to your two input CSV files
# INPUT_CSV_FILE_1 = 'hexme_remesh_i_output_data.csv' # Example: Output from the previous script
# INPUT_CSV_FILE_2 = 'hexme_remesh_i_nocurl_output_data.csv'  # Example: A second file with the same format
INPUT_CSV_FILE_1 = 'solids_output_data.csv' # Example: Output from the previous script
INPUT_CSV_FILE_2 = 'solids_nocurl_output_data.csv'  # Example: A second file with the same format



# Specify the path for the merged and sorted output CSV file
# OUTPUT_MERGED_CSV_FILE = 'merged_sorted_data.csv'
OUTPUT_MERGED_CSV_FILE = 'merged_solids_sorted_data.csv'



# Column index (0-based) of the mesh name in both files
MESH_NAME_COL_IDX = 0
# Column index (0-based) of the value to sort by (from file 1)
SORT_COL_IDX = 1
# --- End Configuration ---

def merge_and_sort_csv(file1_path, file2_path, output_path):
    """
    Merges two CSV files based on MeshName, sorts by the first numeric column
    from file 1, and writes the output.

    Args:
        file1_path (str): Path to the first input CSV file.
        file2_path (str): Path to the second input CSV file.
        output_path (str): Path to the output merged CSV file.
    """
    data_file1 = {}
    headers1 = []

    # --- Read File 1 ---
    print(f"Reading data from: {file1_path}")
    try:
        with open(file1_path, 'r', newline='') as infile1:
            reader1 = csv.reader(infile1)
            try:
                headers1 = next(reader1) # Read header
                if len(headers1) <= max(MESH_NAME_COL_IDX, SORT_COL_IDX):
                     print(f"ERROR: File 1 '{file1_path}' has fewer columns than expected.", file=sys.stderr)
                     return False
            except StopIteration:
                print(f"ERROR: File 1 '{file1_path}' is empty or has no header.", file=sys.stderr)
                return False

            for i, row in enumerate(reader1):
                if len(row) <= max(MESH_NAME_COL_IDX, SORT_COL_IDX):
                    print(f"Warning: Skipping row {i+2} in {file1_path} due to insufficient columns.", file=sys.stderr)
                    continue
                mesh_name = row[MESH_NAME_COL_IDX]
                # Store the *entire row* (excluding mesh name if desired, but easier to keep it)
                # to handle potential variations in number of data columns later
                data_file1[mesh_name] = row
    except FileNotFoundError:
        print(f"ERROR: Input file not found: {file1_path}", file=sys.stderr)
        return False
    except Exception as e:
        print(f"ERROR: Failed to read {file1_path}: {e}", file=sys.stderr)
        return False

    if not data_file1:
        print(f"Warning: No data read from {file1_path}.", file=sys.stderr)
        # Continue if file2 might have data, but merging won't happen.

    print(f"Read {len(data_file1)} data rows from {file1_path}.")

    # --- Read File 2 and Merge ---
    print(f"Reading data from: {file2_path}")
    merged_data = []
    headers2 = []
    try:
        with open(file2_path, 'r', newline='') as infile2:
            reader2 = csv.reader(infile2)
            try:
                headers2 = next(reader2) # Read header
                if len(headers2) <= MESH_NAME_COL_IDX:
                     print(f"ERROR: File 2 '{file2_path}' has fewer columns than expected for MeshName.", file=sys.stderr)
                     return False
            except StopIteration:
                print(f"ERROR: File 2 '{file2_path}' is empty or has no header.", file=sys.stderr)
                return False

            rows_processed_f2 = 0
            rows_merged = 0
            for i, row2 in enumerate(reader2):
                rows_processed_f2 += 1
                if len(row2) <= MESH_NAME_COL_IDX:
                    print(f"Warning: Skipping row {i+2} in {file2_path} due to insufficient columns.", file=sys.stderr)
                    continue

                mesh_name = row2[MESH_NAME_COL_IDX]

                # Check if mesh_name exists in data from file 1
                if mesh_name in data_file1:
                    row1 = data_file1[mesh_name]

                    # Combine data: MeshName, DataColsFile1, DataColsFile2
                    # Exclude mesh name from row1 and row2 when combining data values
                    data_cols1 = row1[MESH_NAME_COL_IDX+1:]
                    data_cols2 = row2[MESH_NAME_COL_IDX+1:]

                    # Attempt conversion of the sort column value from file 1 to float
                    try:
                        sort_value = float(row1[SORT_COL_IDX])
                    except ValueError:
                         print(f"Warning: Skipping merge for '{mesh_name}'. Cannot convert sort value '{row1[SORT_COL_IDX]}' from {file1_path} to float.", file=sys.stderr)
                         continue # Skip if the sort value isn't numeric

                    # Create the merged row structure including the sort value for sorting
                    # Store values as strings initially, handle conversion during write/sort if needed
                    # Store as [sort_value_float, mesh_name, data_cols1_list, data_cols2_list]
                    merged_row_data = [sort_value, mesh_name] + data_cols1 + data_cols2
                    merged_data.append(merged_row_data)
                    rows_merged += 1

    except FileNotFoundError:
        print(f"ERROR: Input file not found: {file2_path}", file=sys.stderr)
        return False
    except Exception as e:
        print(f"ERROR: Failed to read or process {file2_path}: {e}", file=sys.stderr)
        return False

    print(f"Processed {rows_processed_f2} data rows from {file2_path}. Found {rows_merged} matching rows to merge.")

    if not merged_data:
        print("No common rows found between the two files. Output file will not be created.", file=sys.stderr)
        return False

    # --- Sort Merged Data ---
    print(f"Sorting {len(merged_data)} merged rows based on column {SORT_COL_IDX+1} (from {file1_path})...")
    # Sort based on the first element, which we stored as the float sort_value
    merged_data.sort(key=lambda row: row[0])

    # --- Write Output File ---
    print(f"Writing merged and sorted data to: {output_path}")
    try:
        with open(output_path, 'w', newline='') as outfile:
            csv_writer = csv.writer(outfile)

            # Create new headers
            # Headers from File 1 (excluding mesh name)
            header_f1_data = headers1[MESH_NAME_COL_IDX+1:]
            # Headers from File 2 (excluding mesh name)
            header_f2_data = headers2[MESH_NAME_COL_IDX+1:]
            # Add suffixes to distinguish origins
            output_header = [headers1[MESH_NAME_COL_IDX]] + \
                            [f"{h}_F1" for h in header_f1_data] + \
                            [f"{h}_F2" for h in header_f2_data]
            csv_writer.writerow(output_header)

            # Write the sorted data rows
            # Skip the first element (the sort key float) when writing
            for row_data in merged_data:
                csv_writer.writerow(row_data[1:]) # Write from mesh_name onwards

    except IOError as e:
        print(f"ERROR: Could not write to output file {output_path}: {e}", file=sys.stderr)
        return False
    except Exception as e:
        print(f"ERROR: An unexpected error occurred during writing: {e}", file=sys.stderr)
        return False

    print("Merge and sort process completed successfully.")
    return True

# --- Main Execution ---
if __name__ == "__main__":
    # Verify input files exist before calling merge function
    if not os.path.exists(INPUT_CSV_FILE_1):
        print(f"ERROR: Input file 1 not found: {INPUT_CSV_FILE_1}", file=sys.stderr)
        sys.exit(1)
    if not os.path.exists(INPUT_CSV_FILE_2):
        print(f"ERROR: Input file 2 not found: {INPUT_CSV_FILE_2}", file=sys.stderr)
        sys.exit(1)

    merge_and_sort_csv(INPUT_CSV_FILE_1, INPUT_CSV_FILE_2, OUTPUT_MERGED_CSV_FILE)