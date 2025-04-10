import io
import csv
import sys # Import sys to write directly to standard output

# --- File Paths ---
file_path1 = '/Users/fool/Documents/mint/gpgpt/scripts/solids_overflow_output_data.csv'
file_path2 = '/Users/fool/Documents/mint/gpgpt/scripts/solids_output_data.csv'

# --- Data Structures ---
existing_names = set()
data1_rows = [] # To store rows from the first file (as dictionaries)
rows_to_add = [] # To store unique rows from the second file (as dictionaries)
header = [] # To store the header row

# --- Process the first file ---
try:
    with open(file_path1, 'r', newline='', encoding='utf-8') as infile1:
        reader1 = csv.DictReader(infile1)
        header = reader1.fieldnames # Get header from the first file
        if not header:
             print(f"Error: Could not read header from {file_path1}", file=sys.stderr)
             sys.exit(1)
        # Ensure 'MeshName' is in the header
        if 'MeshName' not in header:
             print(f"Error: 'MeshName' column not found in header of {file_path1}", file=sys.stderr)
             print(f"Header found: {header}", file=sys.stderr)
             sys.exit(1)

        for row in reader1:
             # Check if row is not empty or malformed
             if row and 'MeshName' in row and row['MeshName'] is not None:
                 mesh_name = row['MeshName']
                 existing_names.add(mesh_name)
                 data1_rows.append(row)
             else:
                 print(f"Warning: Skipping potentially malformed row in {file_path1}: {row}", file=sys.stderr)

except FileNotFoundError:
    print(f"Error: File not found at {file_path1}", file=sys.stderr)
    sys.exit(1)
except Exception as e:
    print(f"An error occurred while processing {file_path1}: {e}", file=sys.stderr)
    sys.exit(1)


# --- Process the second file ---
try:
    with open(file_path2, 'r', newline='', encoding='utf-8') as infile2:
        reader2 = csv.DictReader(infile2)
        # Check if header matches (optional but good practice)
        if reader2.fieldnames != header:
             print(f"Warning: Headers in {file_path1} and {file_path2} do not match.", file=sys.stderr)
             print(f"Header 1: {header}", file=sys.stderr)
             print(f"Header 2: {reader2.fieldnames}", file=sys.stderr)
             # Decide if you want to exit or continue based on mismatch
             # sys.exit(1) # Uncomment to exit if headers must match

        # Ensure 'MeshName' is available in the second file's processing logic
        if 'MeshName' not in reader2.fieldnames:
             print(f"Error: 'MeshName' column not found in header of {file_path2}", file=sys.stderr)
             sys.exit(1)

        for row in reader2:
            # Check if row is not empty or malformed
            if row and 'MeshName' in row and row['MeshName'] is not None:
                mesh_name = row['MeshName']
                if mesh_name not in existing_names:
                    rows_to_add.append(row)
                    existing_names.add(mesh_name) # Add to set to avoid duplicates if present multiple times in file 2
            else:
                 print(f"Warning: Skipping potentially malformed row in {file_path2}: {row}", file=sys.stderr)


except FileNotFoundError:
    print(f"Error: File not found at {file_path2}", file=sys.stderr)
    sys.exit(1)
except Exception as e:
    print(f"An error occurred while processing {file_path2}: {e}", file=sys.stderr)
    sys.exit(1)


# --- Output the combined data to standard output (console) ---
# Use sys.stdout as the output file
writer = csv.DictWriter(sys.stdout, fieldnames=header, lineterminator='\n')

writer.writeheader()       # Write the header row
writer.writerows(data1_rows) # Write the original rows from file 1
writer.writerows(rows_to_add)  # Write the new unique rows from file 2

# --- Optional: Write to a new file instead of printing ---
# output_filename = 'merged_output.csv'
# try:
#     with open(output_filename, 'w', newline='', encoding='utf-8') as outfile:
#         writer = csv.DictWriter(outfile, fieldnames=header)
#         writer.writeheader()
#         writer.writerows(data1_rows)
#         writer.writerows(rows_to_add)
#     print(f"\nMerged data successfully written to {output_filename}")
# except Exception as e:
#     print(f"\nAn error occurred while writing to {output_filename}: {e}", file=sys.stderr)