import os
import subprocess
import threading

def process_directories(root_dir, frame_path=None):
    """
    Processes directories under root_dir, finding the frame with the largest ID 
    or using a specified frame, and computes statistics using fields_viewer.

    Args:
        root_dir: The root directory containing subdirectories to process.
        frame_path: Optional. Path to a specific frame file to use instead of 
                    finding the largest ID frame.
    """

    for dir_name in os.listdir(root_dir):
        dir_path = os.path.join(root_dir, dir_name)
        
        # Only process if it's a directory
        if os.path.isdir(dir_path):
            print(f"Processing directory: {dir_path}")
            
            # Determine mesh path (assuming it's named like the directory)
            mesh_name = dir_name.split("__")[1] + ".mesh" # extract mesh name from dir name
            mesh_path = os.path.join(dir_path, mesh_name)

            # Determine output path for statistics
            out_path = os.path.join(dir_path, "statistics.txt")

            # Determine frames path
            frames_dir = os.path.join(dir_path, "inner_iters", "frames")

            if frame_path:
                # Use the specified frame path
                curr_frames = frame_path
            else:
                # Find the frame with the largest ID
                if os.path.isdir(frames_dir):
                    largest_frame = find_largest_frame(frames_dir)
                    if largest_frame:
                        curr_frames = os.path.join(frames_dir, largest_frame)
                    else:
                        print(f"No frames found in {frames_dir}. Skipping directory.")
                        continue
                else:
                    print(f"Frames directory not found: {frames_dir}. Skipping directory.")
                    continue
            
            # Construct the command for fields_viewer
            command = f"./bin/10_fields_viewer -m {mesh_path} -f {curr_frames} -s {out_path} --statistics"
            print(command)
            
            # Execute the command in a separate thread (with optional stack size adjustment for Linux)
            thread = threading.Thread(target=execute_command, args=(command,))
            thread.start()
            thread.join() # Wait for the thread to finish before processing the next directory


def find_largest_frame(frames_dir):
    """
    Finds the frame file with the largest numerical ID in the given directory.

    Args:
        frames_dir: The directory containing frame files.

    Returns:
        The filename of the frame with the largest ID, or None if no frames are found.
    """
    largest_frame = None
    largest_id = -1

    for filename in os.listdir(frames_dir):
        if filename.startswith("f_") and filename.endswith(".bfra"):
            try:
                frame_id = int(filename.split("_")[1].split(".")[0])
                if frame_id > largest_id:
                    largest_id = frame_id
                    largest_frame = filename
            except ValueError:
                pass  # Ignore files that don't match the expected format

    return largest_frame


def execute_command(command):
    """
    Executes a shell command, optionally adjusting stack size on Linux.
    """
#ifdef __linux__
    # Set stack size to 8GB using ulimit (example value, adjust as needed)
    subprocess.run("ulimit -s 8096000", shell=True)  
#endif
    subprocess.run(command, shell=True)


# Example usage:
root_directory = "/Users/fool/mint/results/test"  # Replace with your root directory

# Process all directories, using the largest ID frame by default:
process_directories(root_directory)

# Process all directories, but use a specific frame for all:
# specific_frame = "/Users/fool/mint/results/test/12_30_6_30__oloid_4kT_ss_1.00e+00_is_1.00e-05_o_5.00e+02_u_2.50e+00_vsc_2.50e-07_bsdf/inner_iters/frames/f_100050.bfra" 
# process_directories(root_directory, frame_path=specific_frame)