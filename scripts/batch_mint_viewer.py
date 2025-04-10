import os
import subprocess
import threading
import shutil

def process_directories(root_dir, frame_path=None):
    """
    Processes directories under root_dir, finding the frame with the largest ID 
    or using a specified frame, and computes statistics using fields_viewer.

    Args:
        root_dir: The root directory containing subdirectories to process.
        frame_path: Optional. Path to a specific frame file to use instead of 
                    finding the largest ID frame.
    """

    # Get the absolute path of the root directory
    abs_root_dir = os.path.abspath(root_dir)
    # Get the base name of the root directory
    # root_dir_name = os.path.basename(abs_root_dir)
    path_components = abs_root_dir.split(os.sep)

    # out_path_seamless = os.path.expanduser("/mnt/4TBSSD/mint_results_siggraph_2025/" + path_components[-2] + "/" + path_components[-1] + "_seamless/") 
    # out_path_intgrid = os.path.expanduser("/mnt/4TBSSD/mint_results_siggraph_2025/" + path_components[-2] + "/" + path_components[-1] +"_intgrid/") 
    out_path_seamless = os.path.expanduser("/Users/fool/Documents/mint/results/render_out/" + path_components[-2] + "/" + path_components[-1] + "_seamless/") 
    out_path_intgrid = os.path.expanduser("/Users/fool/Documents/mint/results/render_out/" + path_components[-2] + "/" + path_components[-1] +"_intgrid/") 



    os.makedirs(out_path_seamless, exist_ok=True)
    os.makedirs(out_path_intgrid, exist_ok=True)

    print("out directory")
    print(out_path_seamless)

    os.makedirs(out_path_seamless + "/screenshots", exist_ok=True)
    os.makedirs(out_path_intgrid + "/screenshots", exist_ok=True)

    for dir_name in os.listdir(root_dir):
        dir_path = os.path.join(root_dir, dir_name)
        
        # Only process if it's a directory
        if os.path.isdir(dir_path):
            # print(f"Processing directory: {dir_path}")
            
            # Determine mesh path (now assuming it's directly in the directory and has a general name)
            mesh_name = ""
            for file_name in os.listdir(dir_path):
              if ".mesh" in file_name:
                mesh_name = file_name
                print("aaaaaaaaaaaaaaaa")
                print(mesh_name)
            if mesh_name == "":
            #   print("Mesh file not found in this directory. Skipping.")
              continue
            mesh_path = os.path.join(dir_path, mesh_name)

            # Determine output path for statistics
            out_path = os.path.join(dir_path, "param_data")

            # Determine frames path
            frames_dir = os.path.join(dir_path, "inner_iters", "frames")

            largest_frame_id = -1

            if frame_path:
                # Use the specified frame path
                curr_frames = frame_path
            else:
                # Find the frame with the largest ID
                if os.path.isdir(frames_dir):
                    largest_frame, largest_frame_id = find_largest_frame(frames_dir)
                    if largest_frame:
                        if largest_frame_id > 100010:
                            curr_frames = os.path.join(frames_dir, largest_frame)
                        else:
                            print(f"Frame ID too low: {largest_frame}. Skipping directory.")
                            continue
                    else:
                        print(f"No frames found in {frames_dir}. Skipping directory.")
                        continue
                else:
                    print(f"Frames directory not found: {frames_dir}. Skipping directory.")
                    continue

            
            # print("*********************************************")
            # print(f"Processing mesh: {mesh_path}")
            # print(f"Using frame: {curr_frames} (ID: {largest_frame_id})")
            # if (largest_frame_id < 10):
            #     print(f"Frame ID too low: {largest_frame_id}. Skipping directory.")
            #     continue
            # print("*********************************************")

            save_path_seamless = out_path_seamless + mesh_name[:-5]
            save_path_intgrid = out_path_intgrid + "/" + mesh_name[:-5]
            command = f"~/Documents/mint/gpgpt/build/bin/10_fields_viewer -m {mesh_path} -f {curr_frames} -s {save_path_intgrid} --regenerate-type=1 --subdivision --without-gui" #  --statistics" 
            command2 = f"~/Documents/mint/gpgpt/build/bin/10_fields_viewer -m {mesh_path} -f {curr_frames} -s {save_path_seamless} --regenerate-type=2 --subdivision --without-gui" #  --statistics" --without-gui 
            print(command)
                    
            # Execute the command in a separate thread (with optional stack size adjustment for Linux)
            thread = threading.Thread(target=execute_command, args=(command,))
            thread.start()
            thread.join() # Wait for the thread to finish before processing the next directory

            thread = threading.Thread(target=execute_command, args=(command2,))
            thread.start()
            thread.join() # Wait for the thread to finish before processing the next directory

            # copy save_path_intgrid + "/screenshot.png" to out_path_intgrid + "/screenshots/" + mesh_name[:-5] + ".png"


            try:
                # Copy the screenshot from save_path_intgrid to out_path_intgrid
                shutil.copy(
                    os.path.join(save_path_intgrid, "screenshot.png"),
                    os.path.join(out_path_intgrid, "screenshots", mesh_name[:-5] + ".png")
                )
            except FileNotFoundError as e:
                print(f"Error: {e}")

            # copy save_path_seamless + "/screenshot.png" to out_path_seamless + "/screenshots/" + mesh_name[:-5] + ".png"

            try:
                # Copy the screenshot from save_path_seamless to out_path_seamless
                shutil.copy(
                    os.path.join(save_path_seamless, "screenshot.png"),
                    os.path.join(out_path_seamless, "screenshots", mesh_name[:-5] + ".png")
                )
            except FileNotFoundError as e:
                print(f"Error: {e}")



            #### save the crease data 
            # command = f"~/Documents/mint/gpgpt/build/bin/09_mint3d_v2  -m {mesh_path} -c {save_path_intgrid} --headless" #  --statistics" 
            command3 = f"~/Documents/mint/gpgpt/build/bin/09_mint3d_v2 -m {mesh_path} -c {save_path_seamless} --headless" #  --statistics" --without-gui 
            print(command)
                    
            thread = threading.Thread(target=execute_command, args=(command3,))
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

    return (largest_frame, largest_id)


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
# root_directory = os.path.expanduser("~/Documents/mint/results/test/")   # Replace with your root directory
# root_directory = os.path.expanduser("/mnt/4TBSSD/mint_results_siggraph_2025/diag_2.5e-2_500_2.5_start_1e-5/solids/orthog")   # Replace with your root directory
# root_directory = os.path.expanduser("/mnt/4TBSSD/mint_results_siggraph_2025/diag_1e-3_750_2.5_start_1e-5/gallery")   # Replace with your root directory
# root_directory = os.path.expanduser("/mnt/4TBSSD/mint_results_siggraph_2025/diag_1e-3_750_2.5_start_1e-5/gallery")   # Replace with your root directory
# root_directory = os.path.expanduser("/mnt/4TBSSD/mint_results_siggraph_2025/fix_proj_newton_diag_1e-3_100_2.5_start_1e-5/solids")   # Replace with your root directory

# root_directory = os.path.expanduser("/mnt/4TBSSD/mint_results_siggraph_2025/mesh_is1e-5_vsc_5e-3/solids")
# root_directory = os.path.expanduser("/mnt/4TBSSD/mint_results_siggraph_2025/mesh_is1e-8_vsc_1e-1/solids")
# root_directory = os.path.expanduser("/mnt/4TBSSD/mint_results_siggraph_2025/mesh_is1e-8_vsc_1e-2_halfstep/solids")
# root_directory = os.path.expanduser("/mnt/4TBSSD/mint_results_siggraph_2025/smoother_others/hexme")
# root_directory = os.path.expanduser("/mnt/4TBSSD/mint_results_siggraph_2025/smoother_solids/u_1e-4_o_1e-1_feat_5e+1")
# root_directory = os.path.expanduser("/mnt/4TBSSD/mint_results_siggraph_2025/smoother_solids/u_1e-4_o_1e-1_feat_1e+3")
# root_directory = os.path.expanduser("/mnt/4TBSSD/mint_results_siggraph_2025/smoother_solids/mint_o_1e-1_u1e-4")
# root_directory = os.path.expanduser("/mnt/4TBSSD/mint_results_siggraph_2025/smoother_solids/mint_o_1e2_u_2.5_crs_50+hard")
# root_directory = os.path.expanduser("/mnt/4TBSSD/mint_results_siggraph_2025/smoother_solids/mint_o_5e-1_u_1e-2_crs_1e4")


# root_directory = os.path.expanduser("/mnt/4TBSSD/mint_results_siggraph_2025/smoother_solids/mint_0_5e-1_u_1e-2_feat_1e4_maybe")
# root_directory = os.path.expanduser("/mnt/4TBSSD/mint_results_siggraph_2025/smoother_solids/mint_o_5e-1_u_1e-2_feat_5e1")

# root_directory = os.path.expanduser("/mnt/4TBSSD/mint_results_siggraph_2025/mesh_is1e-8_vsc_1e-2_halfstep/solids")
# root_directory = os.path.expanduser("/mnt/4TBSSD/mint_results_siggraph_2025/smoother_solids/mint_o_1e2_u_2.5_crs_1e3")
# root_directory = os.path.expanduser("/mnt/4TBSSD/mint_results_siggraph_2025/mint_solids/to_proc/")

# root_directory = os.path.expanduser("/mnt/4TBSSD/mint_results_siggraph_2025/mint_solids/to_proc/mint_mesh_with_crease_regularity_term")
# root_directory = os.path.expanduser("/mnt/4TBSSD/mint_results_siggraph_2025/mint_solids/to_proc")
root_directory = os.path.expanduser("/Users/fool/Documents/mint/results/to_render/sphere_inset")

# root_directory = os.path.expanduser("/mnt/4TBSSD/mint_results_siggraph_2025/to_proc/baseline_solids")

dirs = os.listdir(root_directory)

for d in dirs:
    if not os.path.isdir(root_directory + "/" + d):
        continue
    # print(root_directory + "/" + d)
    print("blah")
    process_directories(root_directory + "/" + d)
    # print(root_directory + "/" + d)

# Process all directories, using the largest ID frame by default:
# process_directories(root_directory)

# Process all directories, but use a specific frame for all:
# specific_frame = "/Users/fool/mint/results/test/12_30_6_30__oloid_4kT_ss_1.00e+00_is_1.00e-05_o_5.00e+02_u_2.50e+00_vsc_2.50e-07_bsdf/inner_iters/frames/f_100050.bfra" 
# process_directories(root_directory, frame_path=specific_frame)