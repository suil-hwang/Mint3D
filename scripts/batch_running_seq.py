import os
import subprocess
import argparse
import time
import csv
import platform


# Function to convert Windows path to WSL path
def windows_to_wsl_path(windows_path):
    """Convert a Windows path to a WSL-compatible path."""
    return windows_path.replace("C:/", "/mnt/c/").replace("\\", "/")


# Function to detect if the script is running on WSL
def is_wsl():
    """Check if the script is running on WSL (Windows Subsystem for Linux)."""
    return 'microsoft' in platform.uname().release.lower()


# Function to get mesh files from the input folder
def get_mesh_files(input_folder):
    """Get a list of mesh files from the input folder."""
    return [f for f in os.listdir(input_folder) if f.endswith('.mesh')]


# Function to get file sizes for sorting
def get_file_sizes(input_folder, mesh_files):
    """Get the size of each mesh file and return a sorted list."""
    mesh_file_sizes = [
        (f, os.path.getsize(os.path.join(input_folder, f))) for f in mesh_files
    ]
    return sorted(mesh_file_sizes, key=lambda x: x[1])


def run_program_with_env(command, time_limit):
    """Run the C++ program with a timeout and environment variables."""
    try:
        subprocess.run(command, check=True, timeout=time_limit)
        return True
    except subprocess.TimeoutExpired:
        print(f"Timeout expired for {command[1]}. Skipping.")
        return False
    except subprocess.CalledProcessError as e:
        print(f"Error processing {command[1]}: {e}. Skipping.")
        return False
    except Exception as e:
        print(f"Unexpected error with {command[1]}: {e}. Skipping.")
        return False


def process_mesh_file(args, results):
    """Process each mesh file with runtime logging."""
    exe_path, mesh_file, input_folder, output_folder, time_limit = args
    mesh_path = os.path.join(input_folder, mesh_file)
    output_path = os.path.join(output_folder, f"output_{mesh_file}")

    command = [
        exe_path,
        "-m", mesh_path,
        "-o", output_path,
        "--headless"
    ]

    print(f"Processing {mesh_file}...")
    start_time = time.time()
    success = run_program_with_env(command, time_limit)
    end_time = time.time()

    runtime = end_time - start_time
    results.append({"mesh_file": mesh_file, "runtime_seconds": runtime, "status": "success" if success else "failed"})


def write_results_to_csv(results, output_file):
    """Write the results to a CSV file."""
    with open(output_file, mode='w', newline='') as file:
        writer = csv.DictWriter(file, fieldnames=["mesh_file", "runtime_seconds", "status"])
        writer.writeheader()
        writer.writerows(results)


# Function to handle command-line arguments
def parse_arguments():
    parser = argparse.ArgumentParser(description="Batch process mesh files.")
    parser.add_argument('-i', '--input_folder', type=str, required=True, help='Input folder containing mesh files')
    parser.add_argument('-o', '--output_folder', type=str, required=True, help='Output folder for processed files')
    parser.add_argument('-t', '--time_limit', type=int, default=3600 * 24, help='Time limit for each process (in seconds)')
    parser.add_argument('-e', '--exe_path', type=str, default="/home/chopsticks/Projects/gpgpt/build/bin/09_mint3d_v2", help='Path to the executable')
    return parser.parse_args()


def main():
    args = parse_arguments()

    # Check if we're on WSL and convert paths accordingly
    if is_wsl():
        input_folder = windows_to_wsl_path(args.input_folder)
        output_folder = windows_to_wsl_path(args.output_folder)
    else:
        input_folder = args.input_folder
        output_folder = args.output_folder

    exe_path = args.exe_path

    # Create output folder if it does not exist
    if not os.path.exists(output_folder):
        os.makedirs(output_folder)

    # Get and sort mesh files
    mesh_files = get_mesh_files(input_folder)
    sorted_mesh_files = get_file_sizes(input_folder, mesh_files)

    results = []

    print("Starting sequential processing...")
    for mesh_file, _ in sorted_mesh_files:
        process_mesh_file((exe_path, mesh_file, input_folder, output_folder, args.time_limit), results)

    # Save results to CSV
    write_results_to_csv(results, os.path.join(output_folder, "time_log.csv"))

    # Print summary
    print("\nProcessing Summary:")
    successes = [r for r in results if r["status"] == "success"]
    failures = [r for r in results if r["status"] == "failed"]
    print(f"Successfully processed: {len(successes)} files")
    print(f"Failed to process: {len(failures)} files")
    if failures:
        print("Failed files:")
        for failure in failures:
            print(f"  - {failure['mesh_file']}")


if __name__ == "__main__":
    main()
