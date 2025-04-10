#include "SaveForRendering.h"

#include "CubeCover/readMeshFixed.h"

#include "CommonFunctions.h"

#include <filesystem>

bool SaveMeshForRendering(const std::string& save_folder, const Eigen::MatrixXd& V, const Eigen::MatrixXi& T,
                          const Eigen::RowVector3d& center, double scaling_ratio) {
    // save the mesh
    std::string mesh_path = save_folder + "/model.mesh";
    Eigen::MatrixXd V_normalized = NormalizePtsGivenCenterAndScalingRatio(V, center, scaling_ratio);
    CubeCover::exportMESH(mesh_path, V_normalized, T);

    return true;
}

bool SaveScalarFieldForRendering(const std::string& save_folder, const std::string& save_filename,
                                 const Eigen::VectorXd& tet_vals) {   // save the scalar function values
    std::string tet_vals_path = save_folder + "/" + save_filename + ".txt";
    std::ofstream tet_vals_file(tet_vals_path);
    if (!tet_vals_file.is_open()) {
        std::cerr << "Failed to open file " << tet_vals_path << std::endl;
        return false;
    }
    for (int i = 0; i < tet_vals.size(); i++) {
        tet_vals_file << tet_vals(i) << std::endl;
    }
    tet_vals_file.close();
    return true;
}

bool copyFileToOutputFolder(const std::string& input_file, const std::string& output_folder_path) {
    // Define the output file path

    std::cout << input_file << std::endl;
    std::cout << output_folder_path << std::endl;

    std::filesystem::path output_file =
        std::filesystem::path(output_folder_path) / std::filesystem::path(input_file).filename();

    std::cout << output_file << std::endl;

    try {
        // Copy the file
        std::filesystem::copy(input_file, output_file, std::filesystem::copy_options::overwrite_existing);
        std::cout << "File copied to " << output_file << std::endl;
        return true;
    } catch (std::filesystem::filesystem_error& e) {
        std::cerr << "Error copying file: " << e.what() << std::endl;
        return false;
    }
}

bool SaveSegments(const std::string& save_name, const Eigen::MatrixXd& p_start, const Eigen::MatrixXd& p_end,
                  const Eigen::MatrixXd& seg_colors, const Eigen::RowVector3d& center, double scaling_ratio) {
    std::ofstream out(save_name);
    if (!out.is_open()) {
        std::cerr << "Error opening file for writing: " << save_name << std::endl;
        return false;
    }

    if (p_start.rows() != p_end.rows() || p_start.rows() != seg_colors.rows() || p_end.rows() != seg_colors.rows()) {
        std::cerr << "Error: the input matrices should have the same number of rows" << std::endl;
        return false;
    }

    // Iterate over the matrix elements and write them to the file
    for (int i = 0; i < p_start.rows(); ++i) {
        out << (p_start.row(i) - center) / scaling_ratio << " " << (p_end.row(i) - center) / scaling_ratio << " "
            << seg_colors.row(i) << std::endl;
    }

    out.close();   // Close the file
    return true;
}

bool SaveScalars(const std::string& save_name, const Eigen::VectorXd& scalars) {
    std::ofstream out(save_name);
    if (!out.is_open()) {
        std::cerr << "Error opening file for writing: " << save_name << std::endl;
        return false;
    }

    // Iterate over the matrix elements and write them to the file
    for (int i = 0; i < scalars.size(); ++i) {
        out << scalars(i) << std::endl;
    }

    out.close();   // Close the file
    return true;
}