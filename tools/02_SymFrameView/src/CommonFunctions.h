#pragma once

#include <vector>

#include <Eigen/Dense>

#include "CubeCover/StreamlinesExtraction.h"
#include "CubeCover/TetMeshConnectivity.h"

/*
 * Convert a tetrahedral mesh to a soup of tetrahedra.
 * The input is a list of vertices and a list of tetrahedra, where each tetrahedron is represented by a list of 4 vertex
 * indices.
 */
std::pair<Eigen::MatrixXd, Eigen::MatrixXi> GetTetSoup(const Eigen::MatrixXd& V, const Eigen::MatrixXi& T);

/*
 * Convert a tetrahedral mesh to a soup of tetrahedra.
 * The input is a list of vertices and a tetrahedral mesh connectivity.
 */
std::pair<Eigen::MatrixXd, Eigen::MatrixXi> GetTetSoup(const Eigen::MatrixXd& V,
                                                       const CubeCover::TetMeshConnectivity& mesh);

/*
 * Get the surface mesh from a tetrahedral mesh with interior faces, that is, each teterahedron is converted to 4
 * triangles.
 */
std::pair<Eigen::MatrixXd, Eigen::MatrixXi> GetSurfaceMeshFromTetMeshWithInterior(
    const Eigen::MatrixXd& V, const CubeCover::TetMeshConnectivity& mesh);

/*
 * Get the boundary surface mesh from a tetrahedral mesh, that is, only the surface triangles on the boundary are
 * returned.
 */
std::pair<Eigen::MatrixXd, Eigen::MatrixXi> GetBoundarySurfaceMeshFromTetMesh(
    const Eigen::MatrixXd& V, const CubeCover::TetMeshConnectivity& mesh);

void findSharpFeatures(const Eigen::MatrixXd& V, const Eigen::MatrixXi& F, double sharp_feature_threshold,
                       std::vector<Eigen::Vector3d>& nodes, std::vector<Eigen::Vector2i>& edges);

/*
 * Normalize a set of points to [-1/2, 1/2] x [-1/2, 1/2] x [-1/2, 1/2] and record the center and scaling ratio.
 */
Eigen::MatrixXd NormalizePts(const Eigen::MatrixXd& pos, Eigen::RowVector3d& center, double& scaling_ratio);

/*
 * Normalize a set of points given the center and scaling ratio: p_new = (p - center) / scaling_ratio.
 */
Eigen::MatrixXd NormalizePtsGivenCenterAndScalingRatio(const Eigen::MatrixXd& pos, const Eigen::RowVector3d& center,
                                                       double scaling_ratio);

/*
 * Convert a vector of Eigen::Vector3d to an Eigen::MatrixXd.
 */
Eigen::MatrixXd VectorPts2Matrix(const std::vector<Eigen::Vector3d>& pts);

/*
 * Compute the gradient fileds of a scalar function defined on a tetrahedral mesh.
 * Input:
 * - V: the list of vertices.
 * - mesh: the tetrahedral mesh connectivity.
 * - values: the scalar function values defined on the tetrahedral SOUP (Not on the mesh).
 *
 * Return:
 * A ntets x (k * 3) matrix G, where the G.row(i).segment<3>(3 * j) is the gradient fields on i-th tet and j-th value
 * fucntion.
 */
Eigen::MatrixXd ComputeGradient(const Eigen::MatrixXd& V, const CubeCover::TetMeshConnectivity& mesh,
                                const Eigen::MatrixXd& values);

/*
 * Compute the difference between two frames.
 * Input:
 * - frame0: the first frame.
 * - frame1: the second frame.
 * - nvecs: the number of vectors in each frame.
 *
 * Return:
 * A vector of error vectors (err), where err[i][j] is the difference between the i-th vector on j-th element in frame0
 * and frame1.
 *
 * Assumption:
 * frame0 and frame1 has 3 * nvecs columns.
 */
std::vector<Eigen::VectorXd> GetFrameDifference(const Eigen::MatrixXd& frame0, const Eigen::MatrixXd& frame1,
                                                int nvecs);

/*
 * Compute the best match frames in a list of frames for a given frame.
 * Input:
 * - frame_list: a list of frames.
 * - frame: the given frames.
 *
 * Return:
 * A list of frames that are the best match for the given frame.
 *
 * Assumption:
 * - frame_list and frame has 3 * nvecs columns, where nvecs = frame_list.size()
 * - frame_list[i].cols() = 3
 */
std::vector<Eigen::MatrixXd> GetBestMatchFrames(const std::vector<Eigen::MatrixXd>& frame_list,
                                                const Eigen::MatrixXd& frame);

/*
 * Extract the frame vectors from a list of frames, where frames is ntet x (3 * n_frames_per_tet)
 * Return: a list of frame vectors (of size n_frames_per_tet), where each element is a matrix of ntet x 3.
 */
std::vector<Eigen::MatrixXd> ExtractFrameVectors(const Eigen::MatrixXd& frames);

/*
 * Convert hsv color to rgb color.
 */
void HSV2RGB(double h, double s, double v, double& r, double& g, double& b);

/*
 * Convert a vector of S1 to rgb colors.
 */
Eigen::MatrixXd PaintPhi(const Eigen::VectorXd& phi, Eigen::VectorXd* brightness = nullptr);

/*
 * Convert streamlines to start points, end points and colors.
 */
void GetStreamlines(const std::vector<CubeCover::Streamline>& traces, Eigen::MatrixXd& p_start, Eigen::MatrixXd& p_end,
                    Eigen::MatrixXd& colors);

/*
 * Convert streamlines and error values to start points, end points and scalar quantity.
 */
void GetStreamlines(const std::vector<CubeCover::Streamline>& traces, const std::vector<Eigen::VectorXd>& errs,
                    Eigen::MatrixXd& p_start, Eigen::MatrixXd& p_end, Eigen::VectorXd& scalar_quatity,
                    Eigen::MatrixXd* colors = nullptr);

/*
 * Convert the (P, E) streamlines to the (start, end, color) format for rendering.
 */
std::tuple<Eigen::MatrixXd, Eigen::MatrixXd, Eigen::MatrixXd> ConvertCurveNetWorkForRender(
    const Eigen::MatrixXd& P, const Eigen::MatrixXi& E, const Eigen::RowVector4d& const_color);

/*
 * Compute the dihedral angles of two adjacent faces.
 */
double DihedralAngle(const Eigen::MatrixXd& V, const Eigen::MatrixXi& F, int f0, int f1);

/*
 * Check if a tet mesh is valid, that is, it is oriented correctly and is manifold
 */
bool IsValidMesh(Eigen::MatrixXd& tet_pos, CubeCover::TetMeshConnectivity& tet_mesh, bool verbose = false);

/*
 * Check if a frame is valid, that is each tet only contains at most one singular edge
 */
bool IsValidFrame(const CubeCover::TetMeshConnectivity& tet_mesh, const Eigen::MatrixXd& tet_frames,
                  bool verbose = false);

/*
 * Get the tetrahedra matrix from a tetrahedral mesh connectivity
 */
Eigen::MatrixXi GetTetMatrix(const CubeCover::TetMeshConnectivity& tet_mesh);

Eigen::MatrixXd ExpandFramePermutation(const Eigen::MatrixXd& perm, const Eigen::MatrixXd& sgn);
double ComputeCombedL2wrtBasis(const Eigen::VectorXd& f, const Eigen::VectorXd& g, const Eigen::MatrixXd& basis);

void ComputePerTetFacetBasis(const Eigen::MatrixXd& V, const CubeCover::TetMeshConnectivity& mesh,
                             std::vector<std::vector<Eigen::MatrixXd>>& tet_facet_basis,
                             std::vector<std::vector<Eigen::MatrixXd>>& tet_facet_dual_basis,
                             std::vector<std::vector<Eigen::MatrixXd>>& tet_facet_full_basis);

// void GetBestMatchFrames(const std::vector<Eigen::MatrixXd>& frame_list, const Eigen::MatrixXd& frame);