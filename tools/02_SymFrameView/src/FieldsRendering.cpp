#include "FieldsRendering.h"

#include <memory>

#include "polyscope/curve_network.h"
#include "polyscope/point_cloud.h"
#include "polyscope/surface_mesh.h"
#include "polyscope/volume_mesh.h"

#include "CubeCover/FrameField.h"
#include "CubeCover/IsolinesExtraction.h"
#include "CubeCover/SingularCurveNetwork.h"

#include "BoysColoring.h"
#include "CommonFunctions.h"

void RenderStreamlines(const std::vector<CubeCover::Streamline>& traces, const Eigen::MatrixXd& V,
                       const Eigen::MatrixXi& T, std::vector<Eigen::VectorXd>* err, std::string name,
                       bool set_enabled) {
    using namespace CubeCover;
    if (traces.empty()) {
        return;
    }
    int ntraces = traces.size();

    int iter = 0;
    std::vector<Eigen::Vector2d> cur_line_edges;
    std::vector<Eigen::Vector3d> cur_points;
    std::vector<Eigen::Vector3d> cur_colors;
    std::vector<double> cur_err;

    for (int tid = 0; tid < ntraces; tid++) {
        int nsteps = traces.at(tid).stream_pts_.size();

        for (int i = 0; i < nsteps - 1; i++) {
            Eigen::Vector3d edge =
                traces.at(tid).stream_pts_[i].start_pt_ - traces.at(tid).stream_pts_[i + 1].start_pt_;
            int tet_id = traces.at(tid).tet_ids_[i];
            cur_points.push_back(traces.at(tid).stream_pts_[i].start_pt_);

            if (err) {
                double max_err = 0;
                for (int j = 0; j < err->size(); j++) {
                    max_err = std::max(max_err, (*err)[j][tet_id]);
                }
                cur_err.push_back(max_err);
            } else {
                Eigen::Vector3d rgb_color = Boys2RGB(edge);
                cur_colors.push_back(rgb_color);
            }
        }
        cur_points.push_back(traces.at(tid).stream_pts_[nsteps - 1].start_pt_);

        for (int i = 0; i < nsteps - 1; i++) {
            cur_line_edges.push_back(Eigen::Vector2d(iter + i, iter + i + 1));
        }
        iter = iter + nsteps;
    }

    polyscope::CurveNetwork* streamlines;

    streamlines = polyscope::registerCurveNetwork(name + "streamline", cur_points, cur_line_edges);
    streamlines->setTransparency(1);
    streamlines->setRadius(0.002);
    streamlines->setMaterial("flat");
    streamlines->setEnabled(set_enabled);

    if (!cur_colors.empty()) {
        auto color_viewer = streamlines->addEdgeColorQuantity("direction colors", cur_colors);
        // color_viewer->setEnabled(false);
        color_viewer->setEnabled(true);
    }

    if (err) {
        auto color_viewer = streamlines->addEdgeScalarQuantity("Error colors", cur_err);
        color_viewer->setEnabled(true);
    }
}

void RenderIsolines(const Eigen::MatrixXd& P, const Eigen::MatrixXi& E, const Eigen::MatrixXd& P2,
                    const Eigen::MatrixXi& E2, std::string name, bool set_enable) {
    auto* psCurves = polyscope::registerCurveNetwork(name, P, E);
    psCurves->setEnabled(set_enable);
    std::vector<Eigen::Vector3d> isoline_colors;
    for (int i = 0; i < E.rows(); i++) {
        Eigen::Vector3d dir = P.row(E(i, 0)) - P.row(E(i, 1));
        isoline_colors.push_back(Boys2RGB(dir));
    }
    // psCurves->setColor({0, 0, 0});
    // psCurves->setColor(isoline_colors);
    auto color_viewer = psCurves->addEdgeColorQuantity("isoline colors", isoline_colors);
    color_viewer->setEnabled(true);
    psCurves->setRadius(0.005);

    isoline_colors.clear();
    auto* psCurves2 = polyscope::registerCurveNetwork("Bad " + name, P2, E2);
    psCurves2->setEnabled(set_enable);
    for (int i = 0; i < E2.rows(); i++) {
        Eigen::Vector3d dir = P2.row(E2(i, 0)) - P2.row(E2(i, 1));
        isoline_colors.push_back(Boys2RGB(dir));
    }
    psCurves2->setColor({1, 0, 0});
    // psCurves2->setColor(isoline_colors);
    // auto color_viewer = psCurves->addEdgeColorQuantity("direction colors", isoline_colors);

    psCurves2->setRadius(0.002);
}

// void RenderIsolines(const Eigen::MatrixXd& V, const CubeCover::TetMeshConnectivity& mesh, const Eigen::MatrixXd&
// values,
//                     bool set_enable) {
//     Eigen::MatrixXd P;
//     Eigen::MatrixXi E;

//     Eigen::MatrixXd P2;
//     Eigen::MatrixXi E2;

//     CubeCover::ExtractIsolines(V, mesh, values, P, E, P2, E2);

//     auto* psCurves = polyscope::registerCurveNetwork("Isolines", P, E);
//     psCurves->setEnabled(set_enable);

//     psCurves->setColor({0, 0, 0});
//     psCurves->setRadius(0.002);

//     auto* psCurves2 = polyscope::registerCurveNetwork("Bad Isolines", P2, E2);
//     psCurves2->setEnabled(set_enable);

//     psCurves2->setColor({0, 0, 0});
//     psCurves2->setRadius(0.002);
// }

// void RenderIsolines(const Eigen::MatrixXd& P, const Eigen::MatrixXi& E, const Eigen::MatrixXd& P2,
//                     const Eigen::MatrixXi& E2, std::string name, bool set_enable) {
//     auto* psCurves = polyscope::registerCurveNetwork(name, P, E);
//     psCurves->setEnabled(set_enable);
//     std::vector<Eigen::Vector3d> isoline_colors;
//     for (int i = 0; i < E.rows(); i++) {
//         Eigen::Vector3d dir = P.row(E(i, 0)) - P.row(E(i, 1));
//         isoline_colors.push_back(Boys2RGB(dir));
//     }
//     psCurves->setColor({0, 0, 0});
//     psCurves->setRadius(0.002);

//     isoline_colors.clear();
//     auto* psCurves2 = polyscope::registerCurveNetwork("Bad " + name, P2, E2);
//     psCurves2->setEnabled(set_enable);
//     for (int i = 0; i < E2.rows(); i++) {
//         Eigen::Vector3d dir = P2.row(E2(i, 0)) - P2.row(E2(i, 1));
//         isoline_colors.push_back(Boys2RGB(dir));
//     }
//     psCurves2->setColor({0, 0, 0});
//     psCurves2->setRadius(0.002);
// }

void RenderSingularLines(const Eigen::MatrixXd& V, const CubeCover::TetMeshConnectivity& mesh,
                         const Eigen::MatrixXd& frames, const Eigen::MatrixXi& assignments, std::string name_prefix,
                         bool set_enable) {
    Eigen::MatrixXd Pblack;
    Eigen::MatrixXi Eblack;
    Eigen::MatrixXd Pblue;
    Eigen::MatrixXi Eblue;
    Eigen::MatrixXd Pgreen;
    Eigen::MatrixXi Egreen;

    std::unique_ptr<CubeCover::FrameField> field(CubeCover::fromFramesAndAssignments(mesh, frames, assignments, true));
    if (field) {
        std::cout << "No face assignments provided, recomputing: ";
        std::cout.flush();
        field->computeLocalAssignments();
        std::cout << "found " << field->nSingularEdges() << " singular edges" << std::endl;
        field->combAssignments();
        extractSingularCurveNetwork(V, mesh, *field, Pgreen, Egreen, Pblue, Eblue, Pblack, Eblack);

        auto* green = polyscope::registerCurveNetwork(name_prefix + "Singular Curves (+1/4)", Pgreen, Egreen);
        green->setColor({0.0, 1.0, 0.0});
        green->setTransparency(0.8);
        green->setRadius(0.01);
        green->setEnabled(set_enable);

        auto* blue = polyscope::registerCurveNetwork(name_prefix + "Singular Curves (-1/4)", Pblue, Eblue);
        blue->setColor({0.0, 0.0, 1.0});
        blue->setTransparency(0.8);
        blue->setRadius(0.01);

        blue->setEnabled(set_enable);

        auto* black = polyscope::registerCurveNetwork(name_prefix + "Singular Curves (irregular)", Pblack, Eblack);
        black->setColor({1.0, 0.0, 0.0});
        black->setTransparency(0.8);
        black->setRadius(0.005);

        black->setEnabled(set_enable);
    }
}

void RenderError(const Eigen::MatrixXd& V, const CubeCover::TetMeshConnectivity& mesh, const Eigen::MatrixXd& values,
                 const Eigen::MatrixXd& frames, double global_rescaling, bool set_enable) {
    Eigen::MatrixXd grad = ComputeGradient(V, mesh, values);
    grad /= global_rescaling;
    std::vector<Eigen::MatrixXd> grad_vec = ExtractFrameVectors(grad);

    std::vector<Eigen::VectorXd> errs = GetFrameDifference(frames, grad, grad_vec.size());

    if (!polyscope::hasVolumeMesh("tet soup mesh")) {
        Eigen::MatrixXd soup_V;
        Eigen::MatrixXi soup_T;
        std::tie(soup_V, soup_T) = GetTetSoup(V, mesh);
        polyscope::registerVolumeMesh("tet soup mesh", soup_V, soup_T);
    }
    if (!polyscope::hasPointCloud("centroid pc")) {
        std::vector<Eigen::RowVector3d> pts;
        for (int i = 0; i < mesh.nTets(); i++) {
            Eigen::RowVector3d centroid = (V.row(mesh.tetVertex(i, 0)) + V.row(mesh.tetVertex(i, 1)) +
                                           V.row(mesh.tetVertex(i, 2)) + V.row(mesh.tetVertex(i, 3))) /
                                          4;
            pts.push_back(centroid);
        }
        polyscope::registerPointCloud("centroid pc", pts);
    }

    auto tet_mesh = polyscope::getVolumeMesh("tet soup mesh");
    auto pc_mesh = polyscope::getPointCloud("centroid pc");
    for (int i = 0; i < 3; i++) {
        auto grad_vec_viewer = pc_mesh->addVectorQuantity("gradient " + std::to_string(i), grad_vec[i]);
        grad_vec_viewer->setEnabled(false);
        //		tet_mesh->addCellScalarQuantity("error " + std::to_string(i), errs[i]);
    }
    auto err_viewer = tet_mesh->addCellScalarQuantity("error", (errs[0] + errs[1] + errs[2]) / 3);

    // compute combed vector diff operator metrics (combed smoothness, integrability, div)

    std::vector<std::vector<Eigen::MatrixXd>> tet_facet_basis;
    std::vector<std::vector<Eigen::MatrixXd>> tet_facet_dual_basis;
    std::vector<std::vector<Eigen::MatrixXd>> tet_facet_full_basis;

    ComputePerTetFacetBasis(V, mesh, tet_facet_basis, tet_facet_dual_basis, tet_facet_full_basis);

    std::vector<double> combed_smoothness;
    std::vector<double> combed_smoothness_grad;
    std::vector<double> combed_integrability;
    std::vector<double> combed_div;

    double ntets = mesh.nTets();
    for (int i = 0; i < ntets; i++) {
        double n_opp_tet = 0;
        double smoothness = 0;
        double smoothness_grad = 0;
        double integrability = 0;
        double div = 0;
        for (int j = 0; j < 4; j++) {
            Eigen::MatrixXd basis = tet_facet_basis[i][j];
            Eigen::MatrixXd dual_basis = tet_facet_dual_basis[i][j];
            Eigen::MatrixXd full_basis = tet_facet_full_basis[i][j];

            Eigen::VectorXd tf = frames.row(i);
            int opp_tet_idx = mesh.tetOppositeVertex(i, j);
            if (opp_tet_idx == -1) {
                n_opp_tet = n_opp_tet + 1.;
                continue;
            }
            Eigen::VectorXd njf = frames.row(mesh.tetOppositeVertex(i, j));

            integrability += ComputeCombedL2wrtBasis(tf, njf, basis);
            div += ComputeCombedL2wrtBasis(tf, njf, dual_basis);
            smoothness += ComputeCombedL2wrtBasis(tf, njf, full_basis);
            smoothness_grad += ComputeCombedL2wrtBasis(grad.row(i), grad.row(mesh.tetOppositeVertex(i, j)), full_basis);
        }
        double n_neighbor = 4. - n_opp_tet;
        combed_smoothness.push_back(smoothness / n_neighbor);
        combed_smoothness_grad.push_back(smoothness_grad / n_neighbor);
        combed_integrability.push_back(integrability / n_neighbor);
        combed_div.push_back(div / n_neighbor);
    }

    tet_mesh->addCellScalarQuantity("combed_smooth", combed_smoothness);
    tet_mesh->addCellScalarQuantity("combed_smooth_grad", combed_smoothness_grad);
    tet_mesh->addCellScalarQuantity("combed_int", combed_integrability);
    tet_mesh->addCellScalarQuantity("combed_div", combed_div);

    // compute octahedral quality metrics

    std::vector<double> det_vals;
    std::vector<double> sjac_vals;

    for (int i = 0; i < mesh.nTets(); i++) {
        Eigen::Vector3d f0 = grad_vec.at(0).row(i);
        Eigen::Vector3d f1 = grad_vec.at(1).row(i);
        Eigen::Vector3d f2 = grad_vec.at(2).row(i);

        double det = f0.dot(f1.cross(f2));
        double orthog_vol = f0.norm() * f1.norm() * f2.norm();

        det_vals.push_back(det);
        sjac_vals.push_back(std::abs(det) / orthog_vol);
    }

    tet_mesh->addCellScalarQuantity("det", det_vals);
    tet_mesh->addCellScalarQuantity("sjac", sjac_vals);

    err_viewer->setEnabled(set_enable);

    Eigen::MatrixXd surf_V, surf_V_sliced;
    Eigen::MatrixXi surf_F, surf_F_sliced;
    Eigen::VectorXd face_err_sliced;
    std::tie(surf_V, surf_F) = GetSurfaceMeshFromTetMeshWithInterior(V, mesh);

    Eigen::VectorXd face_err(surf_F.rows());

    for (int i = 0; i < surf_F.rows(); i++) {
        int tet_id = i / 4;
        face_err[i] = (errs[0][tet_id] + errs[1][tet_id] + errs[2][tet_id]) / 3;
    }

    /// Boundary alignment error computation

    Eigen::VectorXd face_align_err(mesh.nBoundaryElements());

    for (int i = 0; i < mesh.nBoundaryElements(); i++) {
        int tid = mesh.boundaryElementTet(i);
        int fid = mesh.boundaryFace(i);

        Eigen::Vector3d f0 = grad_vec.at(0).row(tid).normalized();
        Eigen::Vector3d f1 = grad_vec.at(1).row(tid).normalized();
        Eigen::Vector3d f2 = grad_vec.at(2).row(tid).normalized();

        // compute face normals

        Eigen::Vector3d v0 = V.row(mesh.faceVertex(fid, 0));
        Eigen::Vector3d v1 = V.row(mesh.faceVertex(fid, 1));
        Eigen::Vector3d v2 = V.row(mesh.faceVertex(fid, 2));

        Eigen::Vector3d e1 = v1 - v0;
        Eigen::Vector3d e2 = v2 - v0;

        Eigen::Vector3d n = e1.cross(e2);
        n.normalize();

        // max of n.dot(f0), n.dot(f1), n.dot(f2)
        double max_align = std::max(std::max(std::abs(n.dot(f0)), std::abs(n.dot(f1))), std::abs(n.dot(f2)));
        face_align_err[i] = 1 - max_align;
    }

    polyscope::getSurfaceMesh("Boundary Mesh")->addFaceScalarQuantity("alignment_err", face_align_err);

    // boundaryElementTet(int boundaryElement)

    auto surf_mesh = polyscope::registerSurfaceMesh("Surface mesh", surf_V, surf_F);
    surf_mesh->addFaceScalarQuantity("err", face_err);
    surf_mesh->setEnabled(set_enable);
}

void RenderScalarFields(polyscope::VolumeMesh* tet_mesh, const Eigen::MatrixXd& values) {
    for (int i = 0; i < 3; i++) {
        Eigen::MatrixXd vertex_color = PaintPhi(values.col(i));
        tet_mesh->addVertexColorQuantity("color " + std::to_string(i), vertex_color);
    }
}