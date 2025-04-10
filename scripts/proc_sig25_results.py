import os
import shutil

# Define the source and target base directories
results_base = "/Users/fool/Downloads/proc/other_baselines_12_models_rendered"
proc_base = "proc"

# Define the mapping of file paths
file_mappings = {
    "metricguided_seamless/isolines_with_singular_curves.png": "solids_final/isolines_with_singular_curves.png",
    "metricguided_seamless/singular_lines.png": "solids_final/singular_lines.png",
    "metricguided_intgrid/dual_stream_iso_lines.png": "solids_final/dual_stream_iso_lines.png",
    "metricguided_intgrid/interior_err_coolwarm_sjac_sliced.png": "solids_final/interior_err_coolwarm_sjac_sliced.png",
    "metricguided_intgrid/interior_err_coolwarm_int_err_sliced.png": "solids_final/interior_err_coolwarm_int_err_sliced.png",
}

# Traverse the directory structure
for animal in os.listdir(results_base):
    animal_path = os.path.join(results_base, animal)
    if not os.path.isdir(animal_path):
        continue  # Skip non-directory entries

    for source_rel_path, target_rel_path in file_mappings.items():
        source_path = os.path.join(animal_path, source_rel_path)
        target_path = os.path.join(proc_base, animal, target_rel_path)

        if os.path.exists(source_path):
            # Create target directories if they don't exist
            os.makedirs(os.path.dirname(target_path), exist_ok=True)

            # Move or copy the file
            shutil.copy2(source_path, target_path)  # Use shutil.move() to move instead of copying

            print(f"Copied: {source_path} -> {target_path}")
