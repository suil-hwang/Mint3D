#pragma once

#include <Eigen/Dense>

/*
 * Save a mesh for rendering
 * Input:
 * - save_folder: the folder to save the mesh
 * - V: the vertices of the mesh
 * - T: the faces of the mesh
 * - tet_vals: the scalar function values defined on each tet
 * - center: the center of the mesh
 * - scaling_ratio: the scaling ratio of the mesh
 * Remark:
 * - the mesh file is saved as save_folder/model.mesh
 * - the scalar function values are saved as save_folder/tet_err.txt
 * - Save the mesh position as p_new = (p - center) / scaling_ratio.
 *
 */
bool SaveMeshForRendering(const std::string& save_folder, const Eigen::MatrixXd& V, const Eigen::MatrixXi& T,
                          const Eigen::RowVector3d& center, double scaling_ratio);

/*
 * Save a scalar field
 *
 */
bool SaveScalarFieldForRendering(const std::string& save_folder, const std::string& save_filename,
                                 const Eigen::VectorXd& tet_vals);

/*
 * Copy a file to an output folder
 * Input:
 * - input_file: the path to the input file
 * - output_folder: the path to the output folder
 * Output:
 * - Returns true if the file was copied successfully, false otherwise
 * Remark:
 * - The file is copied to output_folder with the same filename as input_file
 * - If the file already exists in the output folder, it will be overwritten
 *
 */
bool copyFileToOutputFolder(const std::string& input_file, const std::string& output_folder);

/*
 * Save a set of segments
 * Input:
 * - save_name: the file name to save the segments
 * - p_start: the start points of the segments
 * - p_end: the end points of the segments
 * - seg_colors: the colors of the segments, should be rgba, or rgb or scalar
 * - center: the center of the segments
 * - scaling_ratio: the scaling ratio of the segments
 * Remark:
 * - Save the segment position as p_new = (p - center) / scaling_ratio.
 */
bool SaveSegments(const std::string& save_name, const Eigen::MatrixXd& p_start, const Eigen::MatrixXd& p_end,
                  const Eigen::MatrixXd& seg_colors, const Eigen::RowVector3d& center, double scaling_ratio);

/*
 * Save a set of scalar values
 * Input:
 * - save_name: the file name to save the scalar values
 * - scalars: the scalar values to save
 * Remark:
 * - Save the scalar values as a column vector in the file.
 */
bool SaveScalars(const std::string& save_name, const Eigen::VectorXd& scalars);