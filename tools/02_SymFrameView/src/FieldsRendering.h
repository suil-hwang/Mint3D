#pragma once

#include <Eigen/Dense>
#include <vector>

#include "polyscope/volume_mesh.h"

#include "CubeCover/TetMeshConnectivity.h"
#include "CubeCover/StreamlinesExtraction.h"

/*
 * Render streamlines on a mesh.
 * Input:
 * - traces: a vector of streamlines
 * - V: vertices of the mesh
 * - T: faces of the mesh
 * - err: error values for each streamline. If provided, the streamline will be colored based on the error values,
 * else the streamline will be colored based on the streamline direction.
 * - name: name of the streamline
 * - set_enabled: whether to set the streamline enabled for visualization
 */
void RenderStreamlines(const std::vector<CubeCover::Streamline>& traces, const Eigen::MatrixXd& V,
                       const Eigen::MatrixXi& T, std::vector<Eigen::VectorXd>* err = nullptr, std::string name = "",
                       bool set_enabled = true);

/*
* Render the isolines of a scalar field on a mesh.
* Input:
* - V: vertices of the mesh
* - mesh: mesh connectivity
* - values: scalar field values
* - set_enable: whether to set the isolines enabled for visualization
*/
void RenderIsolines(const Eigen::MatrixXd& V, const CubeCover::TetMeshConnectivity& mesh, const Eigen::MatrixXd& values,
                    bool set_enable = false);


/*
* Render the isolines of a scalar field on a mesh.
* Input:
* - P: vertices of the good isolines
* - E: edges of the good isolines
* - P2: vertices of the bad isolines
* - E2: edges of the bad isolines
* - name: name of the isolines
* - set_enable: whether to set the isolines enabled for visualization
*/
void RenderIsolines(const Eigen::MatrixXd& P, const Eigen::MatrixXi& E, const Eigen::MatrixXd& P2,
                    const Eigen::MatrixXi& E2, std::string name = "Isolines", bool set_enable = false);

/*
* Render the singular lines of a frame field on a mesh.
* Input:
* - V: vertices of the mesh
* - mesh: mesh connectivity
* - frames: frame field
* - assignments: assignments of the frame field
* - name_prefix: prefix of the name of the singular lines
* - set_enable: whether to set the singular lines enabled for visualization
*/
void RenderSingularLines(const Eigen::MatrixXd& V, const CubeCover::TetMeshConnectivity& mesh,
                         const Eigen::MatrixXd& frames, const Eigen::MatrixXi& assignments,
                         std::string name_prefix = "", bool set_enable = false);

/*
* Render the error between the frame field (scaled by globale_rescaling) and the gradient of the values on a mesh.
* Input:
* - V: vertices of the mesh
* - mesh: mesh connectivity
* - values: scalar field values
* - frames: frame field
* - global_rescaling: global rescaling factor for the frame field
*/
void RenderError(const Eigen::MatrixXd& V, const CubeCover::TetMeshConnectivity& mesh, const Eigen::MatrixXd& values,
                 const Eigen::MatrixXd& frames, double global_rescaling, bool set_enable = false);

/*
* Render scalr fields on a mesh.
*/
void RenderScalarFields(polyscope::VolumeMesh* tet_mesh, const Eigen::MatrixXd& values);