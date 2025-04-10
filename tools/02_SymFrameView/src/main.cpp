#include <filesystem>
#include <iostream>
#include <regex>

#include <igl/file_dialog_open.h>
#include <igl/file_dialog_save.h>
#include <igl/hsv_to_rgb.h>
#include <igl/rgb_to_hsv.h>
#include <igl/writeOBJ.h>

#include <nlohmann/json.hpp>

#include "args.hxx"

#include "polyscope/curve_network.h"
#include "polyscope/point_cloud.h"
#include "polyscope/surface_mesh.h"
#include "polyscope/volume_mesh.h"

#include "CubeCover/CubeCover.h"
#include "CubeCover/FrameField.h"
#include "CubeCover/IsolinesExtraction.h"
#include "CubeCover/MeshSubdivide.h"
#include "CubeCover/ReadFrameField.h"
#include "CubeCover/ReadHexEx.h"
#include "CubeCover/SingularCurveNetwork.h"
#include "CubeCover/SurfaceExtraction.h"
#include "CubeCover/WriteHexEx.h"
#include "CubeCover/readMeshFixed.h"

#include "BoysColoring.h"
#include "CommonFunctions.h"
#include "FieldsRendering.h"
#include "PointGraph.h"
#include "SaveForRendering.h"

enum ParametrizationType {
    kSeamless = 0,
    kIntegerGrid = 1,
};

enum StreamLineType {
    kInit = 0,
    kGradient = 1,
    kInitPerp = 2,
    kInitBestMatchGrad = 3,
    // kInputField = 4,
    kInputFieldPerp = 4,
};

enum StreamLineTracingType {
    kRandomCentroid = 0,   // randomly sample the tet ids and use its centroid as the starting point
    kRandom = 1,           // randomly sample the tet ids and use a random point inside tets for starting point
    kGridPt = 2,           // use integer grid points
};

static std::pair<Eigen::Vector3d, Eigen::Vector3d> GetCutPlaneFromPolyscopeSlice(polyscope::SlicePlane* pl) {
    glm::mat4x4 mat = pl->getTransform();
    Eigen::Vector3d pl_normal(mat[0][0], mat[0][1], mat[0][2]);
    Eigen::Vector3d pl_center(mat[3][0], mat[3][1], mat[3][2]);

    return {pl_normal, pl_center};
}

static bool SerializeMatrix(const Eigen::MatrixXd& mat, const std::string& filepath, int vector_per_element) {
    try {
        std::filesystem::path path(filepath);
        std::ofstream outFile;
        if (path.extension() == ".bfra") {
            outFile = std::ofstream(filepath, std::ios::binary);
        } else {
            outFile = std::ofstream(filepath);
        }

        if (!outFile.is_open()) {
            std::cerr << "Error: Unable to open file for writing: " << filepath << std::endl;
            return false;
        }

        // Write matrix rows and cols
        int rows = static_cast<int>(mat.rows());
        int cols = static_cast<int>(mat.cols());
        int vpe = static_cast<int>(vector_per_element);

        outFile.write("FRA 2", sizeof("FRA 2"));
        outFile.write(reinterpret_cast<char*>(&rows), sizeof(int));
        outFile.write(reinterpret_cast<char*>(&cols), sizeof(int));
        outFile.write(reinterpret_cast<char*>(&vpe), sizeof(int));

        // Write matrix data
        outFile.write(reinterpret_cast<const char*>(mat.data()), rows * cols * sizeof(double));
        outFile.close();
        return true;
    } catch (const std::exception& e) {
        std::cerr << "Error: Unable to serialize matrix: " << e.what() << std::endl;
        return false;
    }
}

static bool IntegrateFrames(const Eigen::MatrixXd& frames, const Eigen::MatrixXd& V, const Eigen::MatrixXi& T,
                            const ParametrizationType& param_type, Eigen::MatrixXi& assignments,
                            Eigen::MatrixXd& values, double global_rescaling = 1.0) {
    CubeCover::CubeCoverOptions opt;
    if (param_type == kSeamless) {
        opt.parameterizationType = CubeCover::CubeCoverOptions::ParameterizationType::PT_SEAMLESS;
        opt.assignmentHandling = CubeCover::CubeCoverOptions::AssignmentHandling::AH_RECOMPUTE;
        opt.boundaryConditions = CubeCover::CubeCoverOptions::BoundaryConditions::BC_FREE;

        opt.solver = CubeCover::CubeCoverOptions::MIPSolver::MS_COMISO;
        // opt.solver = CubeCover::CubeCoverOptions::MIPSolver::MS_GUROBI;

        opt.verbose = true;
    } else {
        opt.parameterizationType = CubeCover::CubeCoverOptions::ParameterizationType::PT_INTEGERGRID;
        opt.assignmentHandling = CubeCover::CubeCoverOptions::AssignmentHandling::AH_RECOMPUTE;
        opt.boundaryConditions = CubeCover::CubeCoverOptions::BoundaryConditions::BC_FORCEINTEGER;
        opt.solver = CubeCover::CubeCoverOptions::MIPSolver::MS_COMISO;
        // opt.solver = CubeCover::CubeCoverOptions::MIPSolver::MS_GUROBI;

        // set to something non-zero if you want curl-correction. 1.0 == 100% change in the input frames allowed.
        opt.curlCorrection = 0.0;

        opt.verbose = true;
    }
    Eigen::MatrixXd to_round = global_rescaling * frames;
    return CubeCover::cubeCover(V, T, to_round, assignments, values, opt);
}

Eigen::MatrixXd V, P1_iso, P2_iso;
Eigen::MatrixXi T, E1_iso, E2_iso;
CubeCover::TetMeshConnectivity mesh;

Eigen::MatrixXd frames, dual_frames, ref_frames, ref_dual_frames;
Eigen::MatrixXi assignments;
Eigen::MatrixXd values;   // the integrated value

double global_rescaling = 15;
int sample_density = 100;
double stream_pt_eps = 1e-2;

double angle_eps = 2.;
double perturb_eps = 0;

// ParametrizationType param_type = kIntegerGrid;
ParametrizationType param_type = kSeamless;
StreamLineType streamline_type = kInitPerp;
StreamLineTracingType streamline_tracing_type = kGridPt;

std::string input_path = "";
std::string output_path = "";
std::string frame_file_path = "";

std::vector<CubeCover::Streamline> input_traces, grad_traces, dual_traces, best_match_traces;

double dihedral_threshold = 0;

polyscope::SlicePlane* pl = nullptr;
Eigen::Vector3d pl_normal, pl_center;

void Save(const std::string& folder) {
    if (!std::filesystem::exists(folder)) {
        std::filesystem::create_directory(folder);
    }

    // save the global rescaling
    std::ofstream rescale_file(folder + "/frames_curr_config.json");
    nlohmann::json j;
    j["global_rescaling"] = global_rescaling;
    rescale_file << j.dump(4);

    Eigen::MatrixXd enlarged_V = 1.005 * V;
    double scaling_ratio;
    Eigen::RowVector3d center;
    NormalizePts(enlarged_V, center, scaling_ratio);

    // log global integrability error

    Eigen::MatrixXd grad = ComputeGradient(V, mesh, values);
    grad /= global_rescaling;
    std::vector<Eigen::MatrixXd> grad_vec = ExtractFrameVectors(grad);
    std::vector<Eigen::VectorXd> errs = GetFrameDifference(frames, grad, grad_vec.size());

    Eigen::VectorXd tet_scalars = errs[0];
    for (int i = 1; i < errs.size(); i++) {
        tet_scalars += errs[i];
    }
    tet_scalars /= errs.size();

    // log frame field diff operator metrics

    std::vector<std::vector<Eigen::MatrixXd>> tet_facet_basis;
    std::vector<std::vector<Eigen::MatrixXd>> tet_facet_dual_basis;
    std::vector<std::vector<Eigen::MatrixXd>> tet_facet_full_basis;

    ComputePerTetFacetBasis(V, mesh, tet_facet_basis, tet_facet_dual_basis, tet_facet_full_basis);

    Eigen::VectorXd combed_smoothness = Eigen::VectorXd::Zero(mesh.nTets());
    Eigen::VectorXd combed_smoothness_grad = Eigen::VectorXd::Zero(mesh.nTets());
    Eigen::VectorXd combed_integrability = Eigen::VectorXd::Zero(mesh.nTets());
    Eigen::VectorXd combed_div = Eigen::VectorXd::Zero(mesh.nTets());

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
        combed_smoothness(i) = (smoothness / n_neighbor);
        combed_smoothness_grad(i) = (smoothness_grad / n_neighbor);
        combed_integrability(i) = (integrability / n_neighbor);
        combed_div(i) = (div / n_neighbor);
    }

    // Log octahedral quality metrics

    Eigen::VectorXd det_vals = Eigen::VectorXd::Zero(mesh.nTets());
    Eigen::VectorXd sjac_vals = Eigen::VectorXd::Zero(mesh.nTets());
    Eigen::VectorXd aniso_vals = Eigen::VectorXd::Zero(mesh.nTets());

    for (int i = 0; i < mesh.nTets(); i++) {
        Eigen::Vector3d f0 = grad_vec.at(0).row(i);
        Eigen::Vector3d f1 = grad_vec.at(1).row(i);
        Eigen::Vector3d f2 = grad_vec.at(2).row(i);

        double det = f0.dot(f1.cross(f2));
        double orthog_vol = f0.norm() * f1.norm() * f2.norm();

        det_vals(i) = det;
        if (orthog_vol < 1e-10) {
            sjac_vals(i) = 0;
        } else {
            sjac_vals(i) = std::abs(det) / orthog_vol;
        }

        double min_val = std::min(std::min(f0.norm(), f1.norm()), f2.norm());
        double max_val = std::max(std::max(f0.norm(), f1.norm()), f2.norm());

        aniso_vals(i) = min_val / max_val;
    }

    SaveMeshForRendering(folder, enlarged_V, T, center, scaling_ratio);
    SaveScalarFieldForRendering(folder, "int_err", tet_scalars);

    SaveScalarFieldForRendering(folder, "combed_smooth", combed_smoothness);
    SaveScalarFieldForRendering(folder, "combed_smooth_grad", combed_smoothness_grad);
    SaveScalarFieldForRendering(folder, "combed_int", combed_integrability);
    SaveScalarFieldForRendering(folder, "combed_div", combed_div);

    SaveScalarFieldForRendering(folder, "det", det_vals);
    SaveScalarFieldForRendering(folder, "sjac", sjac_vals);
    SaveScalarFieldForRendering(folder, "aniso", aniso_vals);

    // Boundary field to log

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

    SaveScalarFieldForRendering(folder, "balign", face_align_err);

    copyFileToOutputFolder(frame_file_path, folder);

    // slicing plane
    if (pl) {
        std::tie(pl_normal, pl_center) = GetCutPlaneFromPolyscopeSlice(pl);
    }
    pl_center = (pl_center - center.transpose()) / scaling_ratio;

    nlohmann::json jval;

    jval["center"] = {pl_center[0], pl_center[1], pl_center[2]};
    jval["normal"] = {-pl_normal[0], -pl_normal[1], -pl_normal[2]};
    std::ofstream ofs(folder + "/slice.json");
    ofs << jval.dump(4) << std::endl;

    // dual stream lines
    dual_traces.clear();
    CubeCover::TraceStreamlines(V, mesh, dual_frames, values, 700, dual_traces, stream_pt_eps, perturb_eps);

    for (int n = 0; n < dual_traces.size(); n++) {
        CubeCover::Streamline merged_s =
            CubeCover::MergeStreamline(dual_traces[n], std::min(1e-2 * stream_pt_eps, 1e-5), angle_eps * M_PI / 180);
        dual_traces[n] = std::move(merged_s);
    }
    Eigen::MatrixXd p_start, p_end, seg_colors;
    Eigen::VectorXd err_scalar;
    GetStreamlines(dual_traces, errs, p_start, p_end, err_scalar, &seg_colors);
    SaveSegments(folder + "/dual_stream_lines.txt", p_start, p_end, seg_colors, center, scaling_ratio);
    SaveScalars(folder + "/dual_stream_lines_scalar.txt", err_scalar);

    // isolines
    Eigen::MatrixXd P;
    Eigen::MatrixXi E;

    Eigen::MatrixXd P2;
    Eigen::MatrixXi E2;
    CubeCover::ExtractIsolines(V, mesh, values, P, E, P2, E2);

    // merge isolines
    PointGraph pt_graph(P, E), pt_graph2(P2, E2);
    pt_graph.BuildConnectivity();
    // merge close points
    pt_graph.MergeClosePoints(1e-6);
    // merge parallel segments
    pt_graph.MergeConsecutiveParallelSegments(angle_eps * M_PI / 180);
    std::tie(P, E) = pt_graph.GetPosEdges();

    pt_graph2.BuildConnectivity();
    // merge close points
    pt_graph2.MergeClosePoints(1e-6);
    // merge parallel segments
    pt_graph2.MergeConsecutiveParallelSegments(angle_eps * M_PI / 180);
    std::tie(P2, E2) = pt_graph2.GetPosEdges();

    Eigen::MatrixXd iso_start_pt((E.rows() + E2.rows()), 3), iso_end_pt((E.rows() + E2.rows()), 3),
        iso_colors((E.rows() + E2.rows()), 4);
    Eigen::MatrixXd iso_start_pt_sliced, iso_end_pt_sliced, iso_colors_sliced;

    for (int i = 0; i < E.rows(); i++) {
        int vid0 = E(i, 0);
        iso_start_pt.row(i) = P.row(vid0);
        int vid1 = E(i, 1);
        iso_end_pt.row(i) = P.row(vid1);
        Eigen::Vector3d dir = P.row(E(i, 0)) - P.row(E(i, 1));
        Eigen::Vector3d c = Boys2RGB(dir);

        // c.transposeInPlace();
        // igl::rgb_to_hsv(c, c);
        // c(0) += 90;   // rotate hue by 90 degrees
        // igl::hsv_to_rgb(c(0), c(1), c(2), c(0), c(1), c(2));

        iso_colors.row(i) << c(0), c(1), c(2), .5;
    }

    for (int i = 0; i < E2.rows(); i++) {
        int vid0 = E2(i, 0);
        iso_start_pt.row(E.rows() + i) = P2.row(vid0);
        int vid1 = E2(i, 1);
        iso_end_pt.row(E.rows() + i) = P2.row(vid1);
        iso_colors.row(E.rows() + i) << 0, 0, 0, 1;
    }
    SaveSegments(folder + "/isolines.txt", iso_start_pt, iso_end_pt, iso_colors, center, scaling_ratio);

    // singular lines
    Eigen::MatrixXd Pblack;
    Eigen::MatrixXi Eblack;
    Eigen::MatrixXd Pblue;
    Eigen::MatrixXi Eblue;
    Eigen::MatrixXd Pgreen;
    Eigen::MatrixXi Egreen;

    CubeCover::FrameField* field = CubeCover::fromFramesAndAssignments(mesh, frames, assignments, true);
    field->computeLocalAssignments();
    std::cout << "found " << field->nSingularEdges() << " singular edges" << std::endl;
    field->combAssignments();
    extractSingularCurveNetwork(V, mesh, *field, Pgreen, Egreen, Pblue, Eblue, Pblack, Eblack);

    Eigen::MatrixXd black_start_pts, black_end_pts, black_color, green_start_pts, green_end_pts, green_color,
        blue_start_pts, blue_end_pts, blue_color;

    std::tie(black_start_pts, black_end_pts, black_color) =
        ConvertCurveNetWorkForRender(Pblack, Eblack, Eigen::RowVector4d(0, 0, 0, 1));
    SaveSegments(folder + "/irregular_singular_curves.txt", black_start_pts, black_end_pts, black_color, center,
                 scaling_ratio);

    std::tie(green_start_pts, green_end_pts, green_color) =
        ConvertCurveNetWorkForRender(Pgreen, Egreen, Eigen::RowVector4d(0, 1, 0, 1));
    SaveSegments(folder + "/one_quarter_singular_curves.txt", green_start_pts, green_end_pts, green_color, center,
                 scaling_ratio);

    std::tie(blue_start_pts, blue_end_pts, blue_color) =
        ConvertCurveNetWorkForRender(Pblue, Eblue, Eigen::RowVector4d(0, 1, 0, 1));
    SaveSegments(folder + "/neg_one_quarter_singular_curves.txt", blue_start_pts, blue_end_pts, blue_color, center,
                 scaling_ratio);

    Eigen::MatrixXd bdryV;
    Eigen::MatrixXi bdryF;
    std::tie(bdryV, bdryF) = GetBoundarySurfaceMeshFromTetMesh(V, mesh);
    igl::writeOBJ(folder + "/boundary_mesh.obj", bdryV, bdryF);

    // Eigen::MatrixXd sharp_start_pts, sharp_end_pts, sharp_color;
    // double sharp_feature_threshold = 25;
    // std::vector<Eigen::Vector3d> sharp_nodes;
    // std::vector<Eigen::Vector2i> sharp_edges;

    // findSharpFeatures(bdryV, bdryF, sharp_feature_threshold, sharp_nodes, sharp_edges);

    // Eigen::MatrixXi sharp_edges_mat(sharp_edges.size(), 2);
    // for (int i = 0; i < sharp_edges.size(); i++) {
    //     sharp_edges_mat.row(i) << sharp_edges[i].x(), sharp_edges[i].y();
    // }
    // Eigen::MatrixXd sharp_nodes_mat(sharp_nodes.size(), 3);
    // for (int i = 0; i < sharp_nodes.size(); i++) {
    //     sharp_nodes_mat.row(i) = sharp_nodes[i];
    // }

    // std::tie(sharp_start_pts, sharp_end_pts, sharp_color) =
    //     ConvertCurveNetWorkForRender(sharp_nodes_mat, sharp_edges_mat, Eigen::RowVector4d(.7, .7, .7, 1));

    // SaveSegments(folder + "/sharp_crease.txt", sharp_start_pts, sharp_end_pts, sharp_color, center, scaling_ratio);
}

void callback() {
    ImGui::PushItemWidth(100);
    ImGui::Text(input_path.c_str());
    if (ImGui::CollapsingHeader("Save Options")) {
        if (ImGui::Button("Save Hexex")) {
            std::string file = igl::file_dialog_save();
            CubeCover::writeHexEx(file, V, T, values);
        }
        if (ImGui::Button("Save Boundary Mesh")) {
            std::string file = igl::file_dialog_save();
            Eigen::MatrixXd bdryV;
            Eigen::MatrixXi bdryF;
            std::tie(bdryV, bdryF) = GetBoundarySurfaceMeshFromTetMesh(V, mesh);
            igl::writeOBJ(file, bdryV, bdryF);
        }
        if (ImGui::Button("Save for Rendering")) {
            if (output_path.empty()) {
                std::string file = igl::file_dialog_save();
                std::string folder = std::filesystem::path(file).parent_path().string();
                Save(folder);
            } else {
                Save(output_path);
            }
            // std::string file Eigen::MatrixXd bdryV;
            // Eigen::MatrixXi bdryF;
            // std::tie(bdryV, bdryF) = GetBoundarySurfaceMeshFromTetMesh(V, mesh);
            // igl::writeOBJ(file, bdryV, bdryF);
        }
    }

    if (ImGui::CollapsingHeader("Integration Options", ImGuiTreeNodeFlags_DefaultOpen)) {
        ImGui::Combo("Paramaterization Type", (int*)&param_type, "Seamless\0Integer grid\0");
        ImGui::InputDouble("Global Rescale", &global_rescaling);
        if (ImGui::Button("Integrate Frames")) {
            if (!IntegrateFrames(frames, V, T, param_type, assignments, values, global_rescaling)) {
                std::cout << "cube cover failed!" << std::endl;
            } else {
                auto tet_mesh = polyscope::getVolumeMesh("tet soup mesh");
                RenderScalarFields(tet_mesh, values);
            }
        }
    }

    ImGui::InputDouble("Parallel Angle Eps", &angle_eps);

    if (ImGui::CollapsingHeader("Streamlines tracing", ImGuiTreeNodeFlags_DefaultOpen)) {
        ImGui::Combo("Streamline Type", (int*)&streamline_type, "Init\0Gradient\0Init Perp\0Best Match Init\0");
        ImGui::Combo("Tracing Pt Type", (int*)&streamline_tracing_type,
                     "Random Sample Centroid\0Random\0Grid Points\0");
        ImGui::InputInt("Sample Density", &sample_density);
        ImGui::InputDouble("Stream Point Eps", &stream_pt_eps);
        if (ImGui::InputDouble("Perturb Eps", &perturb_eps)) {
            if (perturb_eps < 0 && perturb_eps > 1) {
                perturb_eps = 0.1;
            }
        }

        if (ImGui::Button("Grid Points")) {
            std::vector<Eigen::Vector3d> pts;
            for (int i = 0; i < T.rows(); i++) {
                std::vector<std::pair<Eigen::Vector3i, Eigen::Vector3d>> tmp_res =
                    CubeCover::ComputeGridPts(V, mesh, values, i, perturb_eps);
                for (auto& pt : tmp_res) {
                    pts.push_back(pt.second);
                }
            }
            polyscope::registerPointCloud("grid points", pts);
        }

        if (ImGui::Button("Trace Stream lines")) {
            std::vector<CubeCover::Streamline> traces;
            int max_iter_per_trace = 700;

            Eigen::MatrixXd frame_vecs;
            std::string frame_name = "input ";
            std::vector<Eigen::MatrixXd> frame_list;

            auto pc_mesh = polyscope::getPointCloud("centroid pc");

            switch (streamline_type) {
                case kInit: {
                    frame_vecs = frames;
                    frame_list = ExtractFrameVectors(frame_vecs);
                    if (streamline_tracing_type == kGridPt) {
                        CubeCover::TraceStreamlines(V, mesh, frame_vecs, values, max_iter_per_trace, traces,
                                                    stream_pt_eps, perturb_eps);
                    } else {
                        CubeCover::TraceStreamlines(V, mesh, frame_vecs, max_iter_per_trace, traces, sample_density,
                                                    streamline_tracing_type == kRandom, stream_pt_eps);
                    }

                    // Get the error
                    Eigen::MatrixXd grad = ComputeGradient(V, mesh, values);
                    grad /= global_rescaling;
                    std::vector<Eigen::MatrixXd> grad_vec = ExtractFrameVectors(grad);
                    std::vector<Eigen::VectorXd> errs = GetFrameDifference(frames, grad, grad_vec.size());

                    RenderStreamlines(traces, V, T, &errs, frame_name);
                    input_traces = std::move(traces);
                    break;
                }
                case kGradient: {
                    frame_vecs = ComputeGradient(V, mesh, values);
                    frame_vecs /= global_rescaling;
                    frame_name = "gradient ";
                    frame_list = ExtractFrameVectors(frame_vecs);
                    if (streamline_tracing_type == kGridPt) {
                        CubeCover::TraceStreamlines(V, mesh, frame_vecs, values, max_iter_per_trace, traces,
                                                    stream_pt_eps, perturb_eps);
                    } else {
                        CubeCover::TraceStreamlines(V, mesh, frame_vecs, max_iter_per_trace, traces, sample_density,
                                                    streamline_tracing_type == kRandom, stream_pt_eps);
                    }
                    RenderStreamlines(traces, V, T, nullptr, frame_name);
                    grad_traces = std::move(traces);
                    break;
                }
                case kInitPerp: {
                    frame_vecs = dual_frames;
                    frame_name = "dual input ";
                    frame_list = ExtractFrameVectors(frame_vecs);
                    if (streamline_tracing_type == kGridPt) {
                        CubeCover::TraceStreamlines(V, mesh, frame_vecs, values, max_iter_per_trace, traces,
                                                    stream_pt_eps, perturb_eps);
                    } else {
                        CubeCover::TraceStreamlines(V, mesh, frame_vecs, max_iter_per_trace, traces, sample_density,
                                                    streamline_tracing_type == kRandom, stream_pt_eps);
                    }
                    // Get the error
                    Eigen::MatrixXd grad = ComputeGradient(V, mesh, values);
                    grad /= global_rescaling;
                    std::vector<Eigen::MatrixXd> grad_vec = ExtractFrameVectors(grad);
                    std::vector<Eigen::VectorXd> errs = GetFrameDifference(frames, grad, grad_vec.size());

                    RenderStreamlines(traces, V, T, &errs, frame_name);
                    dual_traces = std::move(traces);
                    break;
                }
                case kInitBestMatchGrad: {
                    Eigen::MatrixXd grad_vec = ComputeGradient(V, mesh, values);
                    grad_vec /= global_rescaling;
                    std::vector<Eigen::MatrixXd> grad_list = ExtractFrameVectors(grad_vec);
                    frame_vecs = dual_frames;
                    frame_list = GetBestMatchFrames(grad_list, frames);
                    frame_name = "best match input ";
                    if (streamline_tracing_type == kGridPt) {
                        CubeCover::TraceStreamlines(V, mesh, frame_vecs, values, max_iter_per_trace, traces,
                                                    stream_pt_eps, perturb_eps);
                    } else {
                        CubeCover::TraceStreamlines(V, mesh, frame_vecs, max_iter_per_trace, traces, sample_density,
                                                    streamline_tracing_type == kRandom, stream_pt_eps);
                    }
                    RenderStreamlines(traces, V, T, nullptr, frame_name);
                    best_match_traces = std::move(traces);
                }
                case kInputFieldPerp: {
                    frame_vecs = dual_frames;
                    frame_name = "dual input ";
                    frame_list = ExtractFrameVectors(frame_vecs);
                    if (streamline_tracing_type == kGridPt) {
                        CubeCover::TraceStreamlines(V, mesh, frame_vecs, values, max_iter_per_trace, traces,
                                                    stream_pt_eps, perturb_eps);
                    } else {
                        CubeCover::TraceStreamlines(V, mesh, frame_vecs, max_iter_per_trace, traces, sample_density,
                                                    streamline_tracing_type == kRandom, stream_pt_eps);
                    }
                    RenderStreamlines(traces, V, T, nullptr, frame_name);
                    dual_traces = std::move(traces);
                    break;
                }
            }

            for (int i = 0; i < frame_list.size(); i++) {
                pc_mesh->addVectorQuantity(frame_name + std::to_string(i), frame_list[i]);
            }
        }

        if (ImGui::Button("Merge Stream lines")) {
            std::vector<CubeCover::Streamline> traces;
            std::string frame_name = "input merged ";

            switch (streamline_type) {
                case kInit: {
                    traces = input_traces;
                    break;
                }
                case kGradient: {
                    frame_name = "gradient merged ";
                    traces = grad_traces;
                    break;
                }
                case kInitPerp: {
                    frame_name = "dual input merged ";
                    traces = dual_traces;
                    break;
                }
                case kInitBestMatchGrad: {
                    frame_name = "best match input merged ";
                    traces = best_match_traces;
                }
            }

            int total_count = 0;
            std::cout << "traces size: " << traces.size() << std::endl;
            for (int n = 0; n < traces.size(); n++) {
                std::cout << "merge stream line: " << n << " ";
                CubeCover::Streamline s = traces[n];
                CubeCover::Streamline merged_s =
                    CubeCover::MergeStreamline(s, std::min(1e-2 * stream_pt_eps, 1e-5), angle_eps * M_PI / 180);
                std::cout << ", #merged segments: " << s.stream_pts_.size() - merged_s.stream_pts_.size() << std::endl;
                total_count += s.stream_pts_.size() - merged_s.stream_pts_.size();
                traces[n] = std::move(merged_s);
            }
            std::cout << "total merged size: " << total_count << std::endl;

            RenderStreamlines(traces, V, T, nullptr, frame_name);
        }
    }

    if (ImGui::CollapsingHeader("Visualization Options", ImGuiTreeNodeFlags_DefaultOpen)) {
        if (ImGui::Button("Iso Lines extraction")) {
            Eigen::MatrixXd P;
            Eigen::MatrixXi E;

            Eigen::MatrixXd P2;
            Eigen::MatrixXi E2;

            CubeCover::ExtractIsolines(V, mesh, values, P, E, P2, E2);
            RenderIsolines(P, E, P2, E2, "Isolines", true);

            P1_iso = std::move(P);
            E1_iso = std::move(E);

            P2_iso = std::move(P2);
            E2_iso = std::move(E2);
        }

        if (ImGui::Button("Merge Iso Lines")) {
            PointGraph pt_graph(P1_iso, E1_iso), pt_graph2(P2_iso, E2_iso);
            pt_graph.BuildConnectivity();
            // merge close points
            pt_graph.MergeClosePoints(1e-6);
            // merge parallel segments
            pt_graph.MergeConsecutiveParallelSegments(angle_eps * M_PI / 180);
            std::tie(P1_iso, E1_iso) = pt_graph.GetPosEdges();

            pt_graph2.BuildConnectivity();
            // merge close points
            pt_graph2.MergeClosePoints(1e-6);
            // merge parallel segments
            pt_graph2.MergeConsecutiveParallelSegments(angle_eps * M_PI / 180);
            std::tie(P2_iso, E2_iso) = pt_graph2.GetPosEdges();

            RenderIsolines(P1_iso, E1_iso, P2_iso, E2_iso, "Merged Isolines", true);
        }

        if (ImGui::Button("Iso Surfaces Extraction")) {
            Eigen::MatrixXd isoV;
            Eigen::MatrixXi isoF;

            CubeCover::isosurfaceSoup(V, mesh, values, isoV, isoF);

            auto* ps = polyscope::registerSurfaceMesh("Isosurfaces", isoV, isoF);
        }

        if (ImGui::Button("Compute Error")) {
            RenderError(V, mesh, values, frames, global_rescaling, true);
        }

        if (ImGui::Button("Get Singular lines")) {
            RenderSingularLines(V, mesh, frames, assignments, "Singular lines", true);
        }

        if (ImGui::Button("Draw tet frames around singular lines")) {
            std::unique_ptr<CubeCover::FrameField> field(
                CubeCover::fromFramesAndAssignments(mesh, frames, assignments, true));
            if (field) {
                std::cout << "No face assignments provided, recomputing: ";
                std::cout.flush();
                field->computeLocalAssignments();
                std::cout << "found " << field->nSingularEdges() << " singular edges" << std::endl;
                field->combAssignments();
                std::vector<int> labels = getSingularEdgeLabels(V, mesh, *field);

                for (int i = 0; i < labels.size(); i++) {
                    int eid = field->singularEdge(i);
                    // green curves
                    std::vector<Eigen::Vector3d> pts;
                    std::vector<std::vector<Eigen::Vector3d>> vec_lists(3);
                    CubeCover::AssignmentGroup o(3);
                    for (int j = 0; j < mesh.nEdgeTets(eid); j++) {
                        int tet_id = mesh.edgeTet(eid, j);
                        Eigen::Vector3d pt = V.row(mesh.tetVertex(tet_id, 0)) + V.row(mesh.tetVertex(tet_id, 1)) +
                                             V.row(mesh.tetVertex(tet_id, 2)) + V.row(mesh.tetVertex(tet_id, 3));
                        pt /= 4;
                        pts.emplace_back(pt);

                        int tet = mesh.edgeTet(eid, j);
                        int faceidx = mesh.edgeTetFaceIndex(eid, j, 1);
                        int face = mesh.tetFace(tet, faceidx);
                        int orient = mesh.tetFaceOrientation(tet, faceidx);

                        for (int f1 = 0; f1 < 3; f1++) {
                            Eigen::Vector3d f1_vec = field->tetFrame(tet).row(o.targetVector(f1)).transpose();
                            f1_vec *= o.targetSign(f1);
                            vec_lists[f1].emplace_back(f1_vec);
                        }

                        CubeCover::AssignmentGroup faceo = field->faceAssignment(face);

                        if (orient == 1) {
                            faceo = faceo.inverse();
                        }
                        o = faceo * o;
                    }

                    auto local_frames_viz = polyscope::registerPointCloud(
                        "local frames " + std::to_string(i) + ", edge id " + std::to_string(eid), pts);
                    local_frames_viz->setPointRadius(0);
                    for (int j = 0; j < 3; j++) {
                        auto frame_viz =
                            local_frames_viz->addVectorQuantity("frame " + std::to_string(j), vec_lists[j]);
                        frame_viz->setEnabled(true);
                    }
                }
            }
        }

        ImGui::InputDouble("Angle Threshold", &dihedral_threshold);
    }

    if (ImGui::CollapsingHeader("Subdivision", ImGuiTreeNodeFlags_DefaultOpen)) {
        if (ImGui::Button("Sudivide frames")) {
            bool is_valid = IsValidFrame(mesh, frames, false);
            std::cout << "Is valid frame before subdivision: " << is_valid << std::endl;
            Eigen::MatrixXd V_sub, frames_sub;
            CubeCover::TetMeshConnectivity mesh_sub;
            CubeCover::MeshSubdivision(V, mesh, frames, V_sub, mesh_sub, frames_sub);
            CubeCover::TestMeshSubdivision(V, mesh, frames, V_sub, mesh_sub, frames_sub);
            std::cout << "Is valid mesh after subdivision: " << IsValidMesh(V_sub, mesh_sub) << std::endl;
            std::cout << "Is valid frame after subdivision:  " << IsValidFrame(mesh_sub, frames_sub, false)
                      << std::endl;
            Eigen::MatrixXi assignments_sub;
            RenderSingularLines(V_sub, mesh_sub, frames_sub, assignments_sub, "Singular lines after subdivision", true);
        }
    }
    ImGui::PopItemWidth();
}

// Default arguments
struct {
    std::string mesh_file = "";     // the input mesh file, required input
    std::string frame_file = "";    // the input frame file, required input
    std::string hexex_file = "";    // the input hexex file, optional
    std::string config_json = "";   // the input config json file, optional
    std::string slice_json = "";    // the slice plane json file, optional
    std::string save_folder = "";   // the save folder, optional
    bool no_gui = false;            // without gui option, default is false
    int regenerate = 0;   // regeneration option, 0 for false, 1 for integer-grid, 2 for seamless, default is 0
    bool is_smoothest_combing = true;   // is smoothest combing option, default is true
    bool is_subdivision = false;        // is subdivision option, default is false. If true, will subdivide the mesh and
                                        // frames to make sure that every tet has at most one singular line
    bool renormalize_scale = false;     // renormalize scale option, default is false
} exe_args;

int main(int argc, char** argv) {
    args::ArgumentParser parser("Frame Fields Viewer", "This goes after the options.");
    args::HelpFlag help(parser, "help", "Display this help menu", {'h', "help"});

    // the input mesh and frame file
    args::ValueFlag<std::string> input_mesh_filename(parser, "input mesh", "Input volume mesh file",
                                                     {'m', "input-mesh"}, "");
    args::ValueFlag<std::string> input_frame_filename(parser, "input frames", "Input frame file", {'f', "input-frame"},
                                                      "");
    // the hexex options, we only need at most one of them
    args::ValueFlag<std::string> input_hexex_filename(parser, "input hexex", "Input hexex file", {"hexex"}, "");

    //  args::ValueFlag<std::string> input_intgrid_hexex_filename(parser, "input intgrid hexex",
    //                                                            "Input intgrid hexex file", {"intgrid-hexex"}, "");
    //  args::ValueFlag<std::string> input_seamless_hexex_filename(parser, "input seamless hexex",
    //                                                             "Input seamless hexex file", {"seamless-hexex"}, "");

    // the config json file
    args::ValueFlag<std::string> input_config_json(parser, "input config json", "Input config json file",
                                                   {"config-json"}, "");
    // the slice json file
    args::ValueFlag<std::string> input_slice_json(parser, "input slice json", "Input slice json file", {"slice-json"},
                                                  "");

    // the save folder
    args::ValueFlag<std::string> input_save_folder(parser, "input save folder", "Input save folder",
                                                   {'s', "save-folder"}, "");

    // the with gui flag
    args::Flag without_gui_flag(parser, "without gui", "Without GUI", {"without-gui"});

    // the regenerate flag
    args::ValueFlag<int> regenerate_type(parser, "regenerate", "Is regenerate", {'r', "regenerate-type"}, 0);

    // the smoothest combing flag
    args::Flag is_smoothest_combing_flag(parser, "is smoothest combing", "Is smoothest combing", {"smoothest-combing"});

    // the subdivision flag
    args::Flag is_subdivision_flag(parser, "is subdivision", "Is subdivision", {"subdivision"});

    // the subdivision flag
    args::Flag renormalize_scale(parser, "renormalize_scale", "renormalize_scale", {"renormalize-scale"});

    // Parse args
    try {
        parser.ParseCLI(argc, argv);
    } catch (args::Help&) {
        std::cout << parser;
        return EXIT_SUCCESS;
    } catch (args::ParseError& e) {
        std::cerr << e.what() << std::endl;
        std::cerr << parser;
        return EXIT_FAILURE;
    }

    // parse to the exe_args
    exe_args.mesh_file = args::get(input_mesh_filename);
    exe_args.frame_file = args::get(input_frame_filename);
    exe_args.hexex_file = args::get(input_hexex_filename);
    exe_args.config_json = args::get(input_config_json);
    exe_args.save_folder = args::get(input_save_folder);
    exe_args.no_gui = args::get(without_gui_flag);
    exe_args.regenerate = args::get(regenerate_type);
    exe_args.is_smoothest_combing = args::get(is_smoothest_combing_flag);
    exe_args.is_subdivision = args::get(is_subdivision_flag);
    exe_args.slice_json = args::get(input_slice_json);
    exe_args.renormalize_scale = args::get(renormalize_scale);

    output_path = exe_args.save_folder;
    frame_file_path = exe_args.frame_file;

    if (exe_args.mesh_file == "") {
        std::cerr << "Please specify an input mesh file" << std::endl;
        return EXIT_FAILURE;
    }

    if (exe_args.mesh_file.find("no_curl") != std::string::npos) {
        input_path = "No Curl";
    } else {
        input_path = "With Curl";
    }

    Eigen::MatrixXi tet_F;
    if (!CubeCover::readMESH(exe_args.mesh_file, V, T, tet_F)) {
        std::cerr << "could not read .mesh file " << exe_args.mesh_file << std::endl;
        return EXIT_FAILURE;
    }

    mesh = CubeCover::TetMeshConnectivity(T);

    if (!CubeCover::readFrameField(exe_args.frame_file, "", T, frames, assignments, true)) {
        std::cerr << "could not read frames/permutations" << std::endl;
        return EXIT_FAILURE;
    }

    int vpe = frames.cols() / 3;
    if (vpe != 3) {
        std::cout << "should always be 3 frames" << std::endl;
        return EXIT_FAILURE;
    }

    if (exe_args.no_gui) {
        if (exe_args.save_folder == "") {
            std::cerr << "Without gui, please provide the save folder" << std::endl;
            return EXIT_FAILURE;
        }

        if (!std::filesystem::exists(exe_args.save_folder)) {
            std::filesystem::create_directory(exe_args.save_folder);
        }
    }

    std::string save_file_suffix = "";
    if (exe_args.is_subdivision) {
        bool is_valid = IsValidFrame(mesh, frames, false);
        std::cout << "Is valid frame before subdivision: " << is_valid << std::endl;
        if (!is_valid) {
            std::cout << "frame is not valid, applying subdivision" << std::endl;
            Eigen::MatrixXd V_sub, frames_sub;
            CubeCover::TetMeshConnectivity mesh_sub;
            CubeCover::MeshSubdivision(V, mesh, frames, V_sub, mesh_sub, frames_sub);
            // CubeCover::TestMeshSubdivision(V, mesh, frames, V_sub, mesh_sub, frames_sub);
            save_file_suffix = "_subdiv";
            V = std::move(V_sub);
            mesh = std::move(mesh_sub);
            frames = std::move(frames_sub);
            T = GetTetMatrix(mesh);
            std::cout << "Is valid mesh after subdivision: " << IsValidMesh(V, mesh) << std::endl;
            std::cout << "Is valid frame after subdivision:  " << IsValidFrame(mesh, frames, false) << std::endl;

            std::filesystem::path meshPath(exe_args.mesh_file);
            std::filesystem::path filePath(exe_args.frame_file);

            std::string cur_file_name = exe_args.save_folder + "/" + meshPath.filename().stem().string() + "_subd.mesh";
            std::cout << cur_file_name << std::endl;
            // cur_file_name = std::regex_replace(cur_file_name, std::regex("\\.mesh$"), save_file_suffix + ".mesh");
            CubeCover::exportMESH(cur_file_name, V, T);

            // cur_file_name = exe_args.frame_file;
            cur_file_name = exe_args.save_folder + "/" + filePath.filename().stem().string() + "_subd.bfra";
            exe_args.frame_file = cur_file_name;
            std::cout << cur_file_name << std::endl;
            // cur_file_name =
            //     std::regex_replace(cur_file_name, std::regex("\\.(bfra|_ascii\\.ff3)$"), save_file_suffix + ".bfra");
            std::cout << cur_file_name << std::endl;
            SerializeMatrix(frames, cur_file_name, 3);
        }
    }

    dual_frames = frames;
    double ave_det = 0;
    for (int i = 0; i < T.rows(); i++) {
        ave_det /= T.rows();
        double fdet = frames.row(i).segment<3>(0).dot(frames.row(i).segment<3>(3).cross(frames.row(i).segment<3>(6)));
        ave_det += fdet;
    }
    std::cout << "average inp det: " << ave_det << std::endl;
    double frames_rescale = pow(ave_det, 1. / 3.);
    frames *= 1. / frames_rescale;
    for (int i = 0; i < T.rows(); i++) {
        for (int j = 0; j < 3; j++) {
            Eigen::Vector3d v1 = CubeCover::GetSingleFrame(frames, vpe, i, (j + 1) % 3);
            Eigen::Vector3d v2 = CubeCover::GetSingleFrame(frames, vpe, i, (j + 2) % 3);
            Eigen::Vector3d v0 = v1.cross(v2);
            dual_frames.row(i).segment<3>(3 * j) = v0.transpose();
        }
    }

    Eigen::RowVector3d center = Eigen::RowVector3d::Zero();
    for (int i = 0; i < V.rows(); i++) {
        center += V.row(i);
    }
    center /= V.rows();

    bool global_scaling_set_from_file = false;

    std::ifstream config_file(exe_args.config_json);
    using json = nlohmann::json;
    json jval;

    if (config_file) {
        config_file >> jval;
    }

    if (jval.contains(std::string_view{"global_rescaling"})) {
        global_scaling_set_from_file = true;
        global_rescaling = jval["global_rescaling"];
    }

    std::cout << "global scaling: " << global_rescaling << std::endl;

    // round
    if ((exe_args.hexex_file == "") || exe_args.regenerate != 0) {
        param_type = kIntegerGrid;
        std::string file_suffix = "_intgrid.hexex";
        if (exe_args.regenerate == 2) {
            param_type = kSeamless;
            file_suffix = "_seamless.hexex";
        }

        if (!global_scaling_set_from_file && exe_args.renormalize_scale) {
            global_rescaling = 100.0;
        }

        // if (param_type == kIntegerGrid) {
        //     global_rescaling = 17.5;
        // }

        IntegrateFrames(frames, V, T, param_type, assignments, values, global_rescaling);

        double jumps = values.col(0).maxCoeff() - values.col(0).minCoeff();
        jumps = std::min(jumps, values.col(1).maxCoeff() - values.col(1).minCoeff());
        jumps = std::min(jumps, values.col(2).maxCoeff() - values.col(2).minCoeff());
        std::string file_name = std::regex_replace(exe_args.frame_file, std::regex("\\.fra$"), ".bfra");

        std::cout << "#jumps: " << jumps << std::endl;

        if (jumps > 1e-6 && !global_scaling_set_from_file && exe_args.renormalize_scale) {
            global_rescaling = global_rescaling * 5. / (jumps);
            // if (param_type == kIntegerGrid) {
            //     global_rescaling /= 2;
            // }
            IntegrateFrames(frames, V, T, param_type, assignments, values, global_rescaling);
        }
        file_name = std::regex_replace(file_name, std::regex("\\.bfra$"), save_file_suffix + file_suffix);
        CubeCover::writeHexEx(file_name, V, T, values);
    } else {
        Eigen::MatrixXd V1;
        Eigen::MatrixXi T1;

        std::string hexex_file_path = exe_args.hexex_file;
        param_type = hexex_file_path.find("intgrid") != std::string::npos ? kIntegerGrid : kSeamless;

        if (!CubeCover::readHexEx(hexex_file_path, V1, T1, values)) {
            std::cerr << "error reading the .hexex file in: " << hexex_file_path << std::endl;
            return EXIT_FAILURE;
        } else {
            if (V1.rows() != V.rows() || (T - T1).norm() != 0) {
                std::cout << "mesh mismatched!" << std::endl;
                if (exe_args.is_subdivision) {
                    std::cout
                        << "force to do subdivision. Possible mismatching is that not providing the subdivided hexex"
                        << std::endl;
                }
                return EXIT_FAILURE;
            }
        }
    }

    // dual stream lines
    Eigen::Vector3d min_corner = V.colwise().minCoeff();
    Eigen::Vector3d max_corner = V.colwise().maxCoeff();
    Eigen::Vector3d diag = max_corner - min_corner;
    stream_pt_eps = diag.maxCoeff() / 20;

    // slice json file
    if (exe_args.slice_json != "" && std::filesystem::path(exe_args.slice_json).extension() == ".json") {
        std::ifstream slice_file(exe_args.slice_json);
        using json = nlohmann::json;
        json js;

        if (slice_file) {
            slice_file >> js;
        }
        if (js.contains("normal")) {
            std::vector<double> vec_n = js["normal"].get<std::vector<double>>();
            pl_normal << vec_n[0], vec_n[1], vec_n[2];
        } else {
            Eigen::Vector3d min_corner = V.colwise().minCoeff();
            Eigen::Vector3d max_corner = V.colwise().maxCoeff();

            double min =
                std::min({max_corner[0] - min_corner[0], max_corner[1] - min_corner[1], max_corner[2] - min_corner[2]});

            if (min == max_corner[0] - min_corner[0]) {
                pl_normal << 1, 0, 0;
            } else if (min == max_corner[1] - min_corner[1]) {
                pl_normal << 0, 1, 0;
            } else {
                pl_normal << 0, 0, 1;
            }
        }
        if (js.contains("center")) {
            std::vector<double> vec_p0 = js["center"].get<std::vector<double>>();
            pl_center << vec_p0[0], vec_p0[1], vec_p0[2];
        } else {
            pl_center = V.colwise().mean();
        }
    } else {
        Eigen::Vector3d min_corner = V.colwise().minCoeff();
        Eigen::Vector3d max_corner = V.colwise().maxCoeff();

        double min =
            std::min({max_corner[0] - min_corner[0], max_corner[1] - min_corner[1], max_corner[2] - min_corner[2]});

        if (min == max_corner[0] - min_corner[0]) {
            pl_normal << 1, 0, 0;
        } else if (min == max_corner[1] - min_corner[1]) {
            pl_normal << 0, 1, 0;
        } else {
            pl_normal << 0, 0, 1;
        }

        pl_center = V.colwise().mean();
    }

    // if (exe_args.no_gui) {
    //     if (exe_args.save_folder == "") {
    //         std::cerr << "Without gui, please provide the save folder" << std::endl;
    //         return EXIT_FAILURE;
    //     }
    //     Save(exe_args.save_folder);
    // } else {
    std::vector<Eigen::MatrixXd> frame_vec = ExtractFrameVectors(frames);
    std::vector<Eigen::Vector3d> pts;
    for (int i = 0; i < T.rows(); i++) {
        std::vector<std::pair<Eigen::Vector3i, Eigen::Vector3d>> tmp_res =
            CubeCover::ComputeGridPts(V, mesh, values, i);
        for (auto& pt : tmp_res) {
            pts.push_back(pt.second);
        }
    }

    Eigen::MatrixXd bdryV;
    Eigen::MatrixXi bdryF;
    std::tie(bdryV, bdryF) = GetBoundarySurfaceMeshFromTetMesh(V, mesh);

    double sharp_feature_threshold = 17;
    std::vector<Eigen::Vector3d> sharp_nodes;
    std::vector<Eigen::Vector2i> sharp_edges;

    findSharpFeatures(bdryV, bdryF, sharp_feature_threshold, sharp_nodes, sharp_edges);

    // for (int i = 0; i < sharp_edges.size(); i += 2) {
    //     std::cout << sharp_nodes[i].transpose() << std::endl;
    //     std::cout << sharp_nodes[i + 1].transpose() << std::endl;
    //     std::cout << sharp_edges[i / 2].transpose() << std::endl;
    // }

    polyscope::init();

    auto* psSharp = polyscope::registerCurveNetwork("sharp_features", sharp_nodes, sharp_edges);
    psSharp->setTransparency(1.);
    psSharp->setColor({.5, .5, .5});
    psSharp->setRadius(0.0005);
    psSharp->setEnabled(true);

    auto* psMesh = polyscope::registerSurfaceMesh("Boundary Mesh", bdryV, bdryF);
    psMesh->setTransparency(0.2);
    psMesh->setSurfaceColor({.9, .85, .7});
    psMesh->setEnabled(false);

    Eigen::MatrixXd soup_V;
    Eigen::MatrixXi soup_T;

    std::tie(soup_V, soup_T) = GetTetSoup(V, T);

    auto tet_mesh = polyscope::registerTetMesh("tet soup mesh", soup_V, soup_T);
    tet_mesh->setEnabled(false);

    std::vector<Eigen::Vector3d> centroids;
    for (int i = 0; i < T.rows(); i++) {
        Eigen::Vector3d c;
        c.setZero();
        for (int j = 0; j < 4; j++) {
            c += V.row(T(i, j));
        }
        c /= 4;
        centroids.emplace_back(c);
    }

    auto pc_mesh = polyscope::registerPointCloud("centroid pc", centroids);
    pc_mesh->setPointRadius(0);
    pc_mesh->setEnabled(false);

    std::vector<int> face_ids = {};
    for (int i = 0; i < 3; i++) {
        pc_mesh->addVectorQuantity("frame " + std::to_string(i), frame_vec[i]);
    }

    RenderScalarFields(tet_mesh, values);

    // iso-lines
    Eigen::MatrixXd P;
    Eigen::MatrixXi E;

    Eigen::MatrixXd P2;
    Eigen::MatrixXi E2;

    CubeCover::ExtractIsolines(V, mesh, values, P, E, P2, E2);

    RenderIsolines(P, E, P2, E2, "Isolines", true);
    P1_iso = std::move(P);
    E1_iso = std::move(E);

    P2_iso = std::move(P2);
    E2_iso = std::move(E2);

    RenderSingularLines(V, mesh, frames, assignments, "Singular lines", true);

    CubeCover::TraceStreamlines(V, mesh, dual_frames, values, 700, dual_traces, stream_pt_eps, perturb_eps);

    // Get the error
    Eigen::MatrixXd grad = ComputeGradient(V, mesh, values);
    grad /= global_rescaling;
    std::vector<Eigen::MatrixXd> grad_vec = ExtractFrameVectors(grad);

    std::vector<Eigen::VectorXd> errs = GetFrameDifference(frames, grad, grad_vec.size());

    RenderStreamlines(dual_traces, V, T, &errs, "dual input err", false);
    RenderStreamlines(dual_traces, V, T, nullptr, "dual input ", true);

    CubeCover::TraceStreamlines(V, mesh, frames, 700, input_traces, sample_density, streamline_tracing_type == kRandom,
                                stream_pt_eps);

    RenderStreamlines(input_traces, V, T, &errs, "input ", false);
    // Render the error
    RenderError(V, mesh, values, frames, global_rescaling, false);

    // Add the callback
    polyscope::state::userCallback = callback;
    polyscope::options::autocenterStructures = false;

    pl = polyscope::addSceneSlicePlane();
    pl->setPose({pl_center[0], pl_center[1], pl_center[2]}, {pl_normal[0], pl_normal[1], pl_normal[2]});
    pl->setActive(false);

    // visualize!

    if (exe_args.no_gui) {
        if (exe_args.save_folder == "") {
            std::cerr << "Without gui, please provide the save folder" << std::endl;
            return EXIT_FAILURE;
        }
        Save(exe_args.save_folder);

        Eigen::MatrixXd isoV;
        Eigen::MatrixXi isoF;
        CubeCover::isosurfaceSoup(V, mesh, values, isoV, isoF);
        igl::writeOBJ(exe_args.save_folder + "/isosurfaces.obj", isoV, isoF);
        auto* ps = polyscope::registerSurfaceMesh("Isosurfaces", isoV, isoF);
        ps->setSurfaceColor({50. / 255., 190. / 255., 150. / 255.});
        ps->setEnabled(true);
        ps->setTransparency(0.3);
        polyscope::options::transparencyRenderPasses = 32;
        // polyscope::options::groundPlaneMode = polyscope::GroundPlaneMode::ShadowOnly;
        polyscope::options::groundPlaneMode = polyscope::GroundPlaneMode::None;

        // if (exe_args.renormalize_scale) {
        auto* dis = polyscope::getStructure(polyscope::CurveNetwork::structureTypeName, "dual input streamline");
        dis->setEnabled(false);
        auto* isolines = polyscope::getStructure(polyscope::CurveNetwork::structureTypeName, "Isolines");
        // isolines->setRadius(0.005);
        polyscope::screenshot(exe_args.save_folder + "/screenshot.png", false);
        // }
        std::cout << "done saving to " << exe_args.save_folder << std::endl;

    } else {
        polyscope::show();
    }
    // }

    return EXIT_SUCCESS;
}