/*
 * This file is the entry point into the reference implementation of Mint3D
 * Author: Josh Vekhter, cleaned up by Zhen Chen
 */
#include <thread>

#include "polyscope/point_cloud.h"
#include "polyscope/polyscope.h"
#include "polyscope/volume_mesh.h"

#include "MiNT_mesh.h"
#include "Mint3DHook.h"

#include "args.hxx"

#include "SymmetricTensor/SymmetricKruskalTensor.h"
#include "SymmetricTensor/SymmetricMomentTensor.h"
#include "SymmetricTensor/SymmetricTensorCommon.h"

#include "MiNT3D/MiNTCommonFunc.h"
#include "MiNT3D/MiNTEnergy.h"
#include "MiNT3D/MiNTModel.h"
#include "MiNT3D/MiNTModel_UnitTests.h"

#include "CubeCover/TetMeshConnectivity.h"
#include "CubeCover/readMesh.h"

#include <filesystem>
#include <iostream>
#include <sstream>
#include "../mint3D_hook/ImGuiWidgets.h"
#include "../mint3D_hook/date.h"

#include <memory>

// #include <Eigen/SPQRSupport>

#include <igl/per_face_normals.h>

// #include <igl/min_quad_with_fixed.h>

#if __linux
#    include <cblas.h>
#endif

// static Mint3DHook *hook = NULL;
static std::unique_ptr<Mint3DHook> hook = nullptr;

void toggleSimulation() {
    hook->appState->shouldLogData = true;
    if (!hook) {
        std::cout << "hook is null" << std::endl;
        return;
    }

    if (hook->isPaused()) {
        hook->run();
    } else {
        hook->pause();
        hook->reset();
    }
}

void resetSimulation() {
    if (!hook) return;

    std::cout << "try to reset" << std::endl;
    hook->reset();
    hook->appState->solveStatus = "reset state";
}

// // Solver function
// Eigen::VectorXd solvePoissonEquationWithNeumannBC(const Eigen::MatrixXd& V, const Eigen::MatrixXi& T, const
// Eigen::MatrixXi& F) {
//     // Define the Laplacian matrix (L)
//     Eigen::SparseMatrix<double> L;
//     igl::cotmatrix(V, T, L);

//     // Define the right-hand side vector (b)
//     Eigen::VectorXd b = Eigen::VectorXd::Ones(V.rows());

//     // Compute boundary facets
//     Eigen::MatrixXi boundaryFacets;
//     igl::boundary_facets(T, boundaryFacets);

//     // Compute face normals
//     Eigen::MatrixXd faceNormals;
//     igl::per_face_normals(V, boundaryFacets, faceNormals);

//     // Modify the right-hand side vector for Neumann boundary conditions
//     for (int i = 0; i < boundaryFacets.rows(); ++i) {
//         for (int j = 0; j < boundaryFacets.cols(); ++j) {
//             int vid = boundaryFacets(i, j);
//             b(vid) += faceNormals.row(i).norm(); // Add contribution from Neumann condition
//         }
//     }

//     // Set up the solver
//     Eigen::SparseLU<Eigen::SparseMatrix<double>> solver;
//     solver.analyzePattern(L);
//     solver.factorize(L);

//     // Solve for v
//     Eigen::VectorXd v = solver.solve(b);

//     if(solver.info() != Eigen::Success) {
//         std::cerr << "Decomposition failed" << std::endl;
//         return Eigen::VectorXd();
//     }

//     return v;
// }
std::tuple<Eigen::MatrixXd, Eigen::MatrixXd, Eigen::MatrixXd> ConvertCurveNetWorkForRender(
    const Eigen::MatrixXd &P, const Eigen::MatrixXi &E, const Eigen::RowVector4d &const_color) {
    Eigen::MatrixXd start_pt(E.rows(), 3), end_pt(E.rows(), 3), colors(E.rows(), 4);

    for (int i = 0; i < E.rows(); i++) {
        int vid0 = E(i, 0);
        start_pt.row(i) = P.row(vid0);
        int vid1 = E(i, 1);
        end_pt.row(i) = P.row(vid1);
        colors.row(i) = const_color;
    }
    return std::tuple<Eigen::MatrixXd, Eigen::MatrixXd, Eigen::MatrixXd>(start_pt, end_pt, colors);
}

bool SaveSegments(const std::string &save_name, const Eigen::MatrixXd &p_start, const Eigen::MatrixXd &p_end,
                  const Eigen::MatrixXd &seg_colors, const Eigen::RowVector3d &center, double scaling_ratio) {
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

void drawGUICallback() {
    ImGui::PushItemWidth(100);   // Make ui elements 100 pixels wide,
                                 // instead of full width. Must have
                                 // matching PopItemWidth() below.

    if (hook->showSimButtons() || true) {
        if (ImGui::CollapsingHeader("Start Simulation.", ImGuiTreeNodeFlags_DefaultOpen)) {
            if (ImGui::Button("Run/Pause Sim")) {
                toggleSimulation();
                // hook->appState->turn logging on
            }

            constexpr int MAX_LENGTH = 128;
            char buffer[MAX_LENGTH];
            strncpy(buffer, hook->appState->experiment_name.c_str(), MAX_LENGTH);
            buffer[MAX_LENGTH - 1] = '\0';   // Ensure null-termination

            if (ImGui::InputText("log folder prefix", buffer, MAX_LENGTH)) {
                // Update the std::string only if InputText returns true, indicating a change
                hook->appState->experiment_name = std::string(buffer);
            }

            if (ImGui::Button("Reset Sim")) {
                resetSimulation();
            }

            // if (ImGui::Button("Take one step")) {
            //     if (hook->isPaused()) {
            //         hook->simulateOneStep();
            //         hook->updateRenderGeometry();
            //     }
            // }

            // ImGui::InputDouble( "identity_weight", &hook->appState->identity_weight_log);

            // file explorer

            int *fileIndex = &hook->appState->fcurr;
            int prevFileIndex = *fileIndex;
            ImGui::PushItemWidth(250);
            if (ImGui::SliderInt("Load File Index", fileIndex, hook->appState->fmin, hook->appState->fmax)) {
                // std::cout << "file index changed" << std::endl;
                std::cout << "prev file index " << prevFileIndex << " | current file index " << *fileIndex << std::endl;
                hook->appState->currentFileID = *fileIndex;

                hook->loadPrimaryData();

                hook->updateRenderGeometry();

                hook->appState->restartFromCurrentFrame = true;

                // hacky, should save model on app state
                // and refactor so that solve doesn't block gui
                Eigen::MatrixXd V = hook->appState->V;
                Eigen::MatrixXi T = hook->appState->T;

                CubeCover::TetMeshConnectivity mesh(T);
                MiNT3D::MiNTModel model;
                model.SetMesh(V, mesh);

                Eigen::MatrixXd frames = hook->appState->frames;
                Eigen::MatrixXd bnd_frames = hook->appState->boundary_frames;

                Eigen::MatrixXd frames_out, bnd_frames_out;

                // this is for reproduce the crashed case
                //                hook->appState->solve_params.w_outer_step = 100;

                model.ShowDensity(hook->appState->solve_params, frames, frames_out, &bnd_frames, &bnd_frames_out);

                Eigen::VectorXd cur_scalar_quantity = model.cur_smooth_density;
                cur_scalar_quantity = cur_scalar_quantity.segment(0, T.rows());
                auto cur_scalar_field =
                    polyscope::getVolumeMesh("c")->addCellScalarQuantity("smoothness", cur_scalar_quantity);

                cur_scalar_quantity = model.cur_asap_smoothness_density;
                cur_scalar_quantity = cur_scalar_quantity.segment(0, T.rows());
                cur_scalar_field =
                    polyscope::getVolumeMesh("c")->addCellScalarQuantity("asap_smoothness", cur_scalar_quantity);

                cur_scalar_quantity = model.cur_aiap_smoothness_density;
                cur_scalar_quantity = cur_scalar_quantity.segment(0, T.rows());
                cur_scalar_field =
                    polyscope::getVolumeMesh("c")->addCellScalarQuantity("aiap_smoothness", cur_scalar_quantity);

                cur_scalar_quantity = model.cur_mint_density;
                cur_scalar_quantity = cur_scalar_quantity.segment(0, T.rows());
                cur_scalar_field =
                    polyscope::getVolumeMesh("c")->addCellScalarQuantity("integrability", cur_scalar_quantity);

                cur_scalar_quantity = model.cur_neohookean_density;
                cur_scalar_quantity = cur_scalar_quantity.segment(0, T.rows());
                cur_scalar_field =
                    polyscope::getVolumeMesh("c")->addCellScalarQuantity("neohookean", cur_scalar_quantity);

                cur_scalar_quantity = model.flagged_frames;
                cur_scalar_quantity = cur_scalar_quantity.segment(0, T.rows());
                cur_scalar_field = polyscope::getVolumeMesh("c")->addCellScalarQuantity(
                    "smoothness and integrability combings disagree", cur_scalar_quantity);

                // // cur_scalar_field->setEnabled(true);
                // double maxbound = cur_scalar_quantity.maxCoeff();
                // double minbound = cur_scalar_quantity.minCoeff();
                // double diff = maxbound - minbound;

                // cur_scalar_field->setMapRange(std::make_pair(1e-8, 1e-4));
                // cur_scalar_field->setEnabled(true);

                // FieldBounds cur_bounds = appState->fieldBounds[appState->current_element];
                // double abs_max = cur_bounds.upper * diff + minbound;
                // double abs_min = cur_bounds.lower * diff + minbound;

                // if (appState->override_bounds_active)
                // {
                //     abs_max = appState->override_bounds.upper;
                //     abs_min = appState->override_bounds.lower;
                // }

                // cur_scalar_field->setMapRange(std::make_pair(abs_min, abs_max));

                // std::cout << density << std::endl;

                // hook->outputData->frames = hook->appState->frames;
                // outputData->frames
                // appState.shouldReload = true;
            }

            bool *load_inner = &hook->appState->loadInner;
            ImGui::Checkbox("load inner iterations", load_inner);

            bool *start_from_current = &hook->appState->restartFromCurrentFrame;
            ImGui::Checkbox("start_from_current", start_from_current);

            if (ImGui::Button("View Integrated Fields")) {
                hook->launchFieldViewerFromCurrentState();
            }

            ImGui::PopItemWidth();

            ///////////////////////////////////
            ////////////// TESTS //////////////
            ///////////////////////////////////

            ImGui::Begin("Solver Tests");
            ImGui::PushItemWidth(150);

            // if (ImGui::Button("Two Tet check")) {
            //     Eigen::MatrixXd frames = Eigen::MatrixXd::Zero(2, 9);
            //     frames << -.690558672, -.458358347, .559496522, .707926512, .10121631, .698996723, .2019715,
            //     -.716825,
            //         .66735899, -.48258, .048645, -.881716, 1.0211516, -.1707977, -.0005055, .484414, .3568398,
            //         -.43527472;
            //     Eigen::MatrixXd V = hook->appState->V;
            //     Eigen::MatrixXi T = hook->appState->T;

            //     CubeCover::TetMeshConnectivity mesh(T);

            //     MiNT3D::MiNTModel_UnitTests model;
            //     model.SetMesh(V, mesh);

            //     MiNT3D::SolveParams solve_params;
            //     Eigen::MatrixXd out_frames;

            //     bool succeed = model.SolveTwoTetTest(solve_params, frames, out_frames);

            //     for (int i = 0; i < out_frames.rows(); i++) {
            //         std::cout << out_frames.row(i) << std::endl;
            //     }

            //     // hook->opt->set_current_x(x);
            // }

            if (ImGui::Button("cylinder fit test case")) {
                hook->setFramesToInvertingOnCylinder(hook->appState->invert_frames_in_cylinder_test,
                                                     hook->appState->cylinder_multiplier,
                                                     hook->appState->noise_multiplier);
                hook->appState->restartFromCurrentFrame = true;
                hook->appState->exp_name = "fit";
                // hook->appState->solve_params.setFitSolve();
                hook->appState->solve_params.setSoftAlignSolve();

                hook->appState->show_frames = true;

                hook->appState->solve_params.boundary_hard_constraints =
                    MiNT3D::BoundaryHardConstraints::NoHardConstraints;
                // hook->appState->solve_params.w_unit = 1e-1;
                // hook->appState->solve_params.w_fit = 1e3;

                Eigen::MatrixXd perturb = Eigen::MatrixXd::Random(hook->appState->frames.rows(), 9);
                hook->appState->frames += perturb * 0.00001;
            }

            //

            // checkbox invert_frames_in_cylinder_test

            // input double cylinder_multiplier
            ImGui::InputDouble("cylinder multiplier", &hook->appState->cylinder_multiplier);
            ImGui::InputDouble("noise multiplier", &hook->appState->noise_multiplier);
            // ImGui::SameLine();
            if (ImGui::Checkbox("invert frames in cylinder test", &hook->appState->invert_frames_in_cylinder_test)) {
                hook->appState->cylinder_multiplier = .25;
            }

            // if (ImGui::Button("Check Correctness")) {
            //     Eigen::MatrixXd frames = hook->appState->frames;
            //     Eigen::MatrixXd V = hook->appState->V;
            //     Eigen::MatrixXi T = hook->appState->T;

            //     CubeCover::TetMeshConnectivity mesh(T);

            //     double tinyAD_energy = hook->opt->eval_func_local(hook->opt->get_current_x());

            //     MiNT3D::MiNTModel model;
            //     model.SetMesh(V, mesh);

            //     Eigen::MatrixXd bnd_frames = hook->appState->boundary_frames;

            //     double mint_energy = model.ComputeEnergy(frames, &bnd_frames);

            //     std::cout << "tiny AD energy: " << tinyAD_energy << ", mint energy: " << mint_energy
            //               << ", difference: " << mint_energy - tinyAD_energy << std::endl;
            // }

            if (ImGui::Button("Set Current Field to Dual")) {
                hook->takeFrameFieldDual();
            }

            ImGuiWidgets::ShowCombedFramesPerFacet(*hook->appState);

            if (ImGui::Button("Test Scaling Invariance")) {
                Eigen::MatrixXd frames = hook->appState->frames;
                Eigen::MatrixXd bnd_frames = hook->appState->boundary_frames;
                Eigen::MatrixXd V = hook->appState->V;
                Eigen::MatrixXi T = hook->appState->T;

                CubeCover::TetMeshConnectivity mesh(T);

                MiNT3D::MiNTModel model;
                model.SetMesh(V, mesh);
                double energy = model.ComputeEnergy(frames, &bnd_frames);

                double rescaling_ratio = 2;

                MiNT3D::MiNTModel model_rescaled;
                Eigen::MatrixXd V_rescaled = V * rescaling_ratio;

                model_rescaled.SetMesh(V_rescaled, mesh);
                double energy_rescaled = model_rescaled.ComputeEnergy(frames, &bnd_frames);

                std::cout << "rescaling factor: " << rescaling_ratio << ", before rescaling: " << energy
                          << ", after rescaling: " << energy_rescaled
                          << ", difference: " << (energy_rescaled / energy - rescaling_ratio) << std::endl;
            }

            if (ImGui::Button("Test Energy and Derivatives")) {
                Eigen::MatrixXd frames = hook->appState->frames;
                Eigen::MatrixXd tmpV = hook->appState->V;
                Eigen::MatrixXi tmpT = hook->appState->T;

                CubeCover::TetMeshConnectivity mesh(tmpT);

                Eigen::MatrixXd tmp_fields(tmpT.rows(), 9), tmp_fields1(tmpT.rows(), 9);
                tmp_fields.setRandom();
                tmp_fields1.setRandom();

                MiNT3D::MiNTEnergy energy(tmpV, mesh, mesh);
                energy.Precompute();
                energy.ComputeGuidedTensors(tmp_fields1);

                //  int tet_id = std::rand() % mesh.nTets();
                //  energy.TestComputeKruskalTensorsPerTet(tmp_fields, tet_id);
                //  energy.TestComputeFaceProjectedTensors(tmp_fields, tet_id);

                int face_id = std::rand() % mesh.nFaces();

                while (mesh.faceTet(face_id, 0) == -1 || mesh.faceTet(face_id, 1) == -1) {
                    face_id = std::rand() % mesh.nFaces();
                }

                //  energy.TestComputeSmoothnessEnergyPerFace(tmp_fields, face_id);
                energy.TestComputeIntegrabilityEnergyPerFace(tmp_fields, face_id);

                //
                //          energy.TestPerFaceSmoothnessHessianProjection(tmp_fields, face_id);
                //
                //          energy.TestComputeSmoothnessEnergy(tmp_fields);

                // energy.TestComputeSmoothnessEnergy(tmp_fields);
                // energy.TestComputeIntegrabilityEnergy(tmp_fields);

                //          energy.TestComputeDeviationEnergy(tmp_fields);
                //  energy.TestComputeUnitNormPenalty(tmp_fields);
                //  energy.TestComputeUnitNormalBoundaryPenalty(tmp_fields);
                // energy.TestComputeUnitNormPenaltyPerTet(tmp_fields, 0);

                // std::cout << energy.ComputeDeviationEnergy(tmp_fields1)
                //           << std::endl; // should be 0
            }

            ImGui::PopItemWidth();

            ImGui::End();
        }
    }
    hook->drawGUI();
    ImGui::PopItemWidth();
    hook->render();
}

int main(int argc, char **argv) {
    // #ifdef __APPLE__

    // This get's rid of a weird issue on linux where you run multple jobs at once, but disabling for condor compile.
    // #elif __linux
    //     openblas_set_num_threads(1);
    // #endif

    // Configure the argument parser
    args::ArgumentParser parser(
        ("This is an example optimization that studies integrability in 2d.\n "
         "Built on top of an LLM assisted codebase called gpgpt. \n\n\n"
         ""));

    args::Flag headless(parser, "headless_mode", "This will avoid all of the polyscope calls and will immediately run",
                        {"headless"});

    args::Flag nomint(parser, "disable_integrability_constraint",
                      "this will disable the integrability term from being enforced as a constraint", {"nomint"});

    args::ValueFlag<std::string> filepath(parser, "input_file_path", "This is the input mesh file path",
                                          {'f', 'm', "meshpath"});
    args::ValueFlag<std::string> inDir(parser, "reload_dir", "load a directory to restart optimization", {'d', "dir"});
    args::ValueFlag<std::string> outDir(parser, "out_dir", "choose directory to write output to", {'o', "out_dir"});

    args::ValueFlag<std::string> parseMetricDrivenFrame(parser, "metric_driven_frame",
                                                        "load frames from metric driven frame 3d", {"ff3"});
    args::ValueFlag<std::string> parseFRA(
        parser, "load_fra", "load frames stored in FRA or BFRA format, requires also passing in the mesh", {"fra"});

    args::Flag subdiv_boundary(parser, "subdivide_boundary",
                               "this will subdivide every tet on the boundary of the input mesh",
                               {'s', "subdivide_boundary"});

    args::ValueFlag<std::string> creaseSavePath(parser, "crease_save_path",
                                                "this is a deadline rendering hack, call with -m and --headless and it "
                                                "will write the creases here, fix the surface in the fields viewer...",
                                                {'c', "crease_save_path"});

    // args::Positional<std::string> inDir(parser, "out_dir", "load a directory to
    // restart optimization");
    args::Positional<std::string> inObj(parser, "in_obj", "load a surface mesh to solve if out_dir is empty");

    // Parse args
    try {
        parser.ParseCLI(argc, argv);
    } catch (args::Help &) {
        std::cout << parser;
        return 0;
    } catch (args::ParseError &e) {
        std::cerr << e.what() << std::endl;
        std::cerr << parser;
        return 1;
    }

    if (!headless) {
        // Options
        polyscope::options::autocenterStructures = false;
        polyscope::view::windowWidth = 1024;
        polyscope::view::windowHeight = 1024;
        polyscope::init();

    } else {
        std::cout << "running MiNT 3D in headless mode" << std::endl;
    }

    // hook = static_cast<Mint3DHook *>(new MiNT_mesh());
    hook = std::make_unique<MiNT_mesh>();

    hook->appState->headless_mode = headless;
    // hook->appState->subdivide_boundary_tets = subdiv_boundary;

    // Should move this out of main later
    ///////////////////
    /// Load Mesh ////
    ///////////////////

    // default path
    hook->appState->meshInputFilePath = "../tools/shared/";

    // hook->appState->meshInputFilePath += "octahedron_8.mesh";
    hook->appState->meshInputFilePath += "cube_lil.mesh";

    if (filepath) {
        std::string path = args::get(filepath);
        std::cout << "loading mesh " << path << std::endl;
        hook->appState->meshInputFilePath = path;
        hook->appState->solve_params.setSoftOrthogSolve();
        hook->appState->exp_name = "orthog";
    } else if (inDir) {
        std::cout << "loading from directory " << hook->appState->directoryPath << std::endl;
        FileParser fp = FileParser(args::get(inDir));
        std::cout << "minid: " << fp.minID << " maxid: " << fp.maxID << std::endl;
        hook->appState->fmin = fp.minID;
        hook->appState->fmax = fp.maxID;
        hook->appState->fcurr = hook->appState->fmax;

        hook->appState->meshInputFilePath = fp.meshFilePath;
        // hook->reset();
        hook->appState->directoryPath = args::get(inDir);
        hook->fileParser = std::make_unique<FileParser>(hook->appState->directoryPath);

        hook->appState->loadLastFile = false;   // need to clean this up....
        hook->appState->shouldReload = true;
        hook->appState->loadedPreviousRun = false;
        hook->appState->subdivide_boundary_tets = false;
        hook->appState->restartFromCurrentFrame = true;

        hook->loadPrimaryData();
        hook->loadSecondaryData();

    } else {
        std::cout << "loading default mesh" << std::endl;
    }

    if (outDir) {
        hook->appState->out_dir = args::get(outDir);
    }

    // if (!inDir) {
    // create directory

    if (parseMetricDrivenFrame || parseFRA) {
        hook->appState->subdivide_boundary_tets = false;
        hook->appState->restartFromCurrentFrame = true;
        hook->appState->exp_name = "align";
    }

    hook->initSimulation();
    hook->reset();

    if (parseMetricDrivenFrame) {
        std::string path = args::get(parseMetricDrivenFrame);
        std::cout << "loading metric driven frame " << path << std::endl;
        hook->loadMetricDrivenFrameField(path);
    } else if (parseFRA) {
        std::string path = args::get(parseFRA);
        std::cout << "loading frames from fra " << path << std::endl;
        hook->loadSpecificFrameField(path);

    } else if (!headless) {
        hook->initializeFieldState();
        hook->updateRenderGeometry();
    }

    if (nomint) {
        hook->appState->solve_params.b_int_combed = false;
        hook->appState->solve_params.b_int_sym = false;
    }

    // try {
    //     std::filesystem::copy_file(source_mesh, targ_mesh);
    // } catch (...) {
    //     std::cout << "failed to copy " << source_mesh << " to directory " << targ_mesh << std::endl;
    // }

    // hook->reset();
    // }

    // std::cout << "nvars in opt: after reset" << hook->opt->get_num_vars() << std::endl;

    if (!headless) {
        polyscope::state::userCallback = drawGUICallback;
        polyscope::options::programName = "gpgpt - MINT3D";
        polyscope::options::verbosity = 1;
        polyscope::options::transparencyRenderPasses = 8;
        polyscope::view::resetCameraToHomeView();
        polyscope::show();

        if (parseMetricDrivenFrame) {
            polyscope::getPointCloud("surface_vecs")->setEnabled(false);
        }

    } else {
        hook->appState->headless_mode = true;
        // hook->run();

        if (creaseSavePath) {
            Eigen::MatrixXd sharp_start_pts, sharp_end_pts, sharp_color;
            // double sharp_feature_threshold = 25;
            std::vector<Eigen::Vector3d> sharp_nodes = hook->appState->sharp_nodes;
            std::vector<Eigen::Vector2i> sharp_edges = hook->appState->sharp_edges;

            Eigen::MatrixXi sharp_edges_mat(sharp_edges.size(), 2);
            for (int i = 0; i < sharp_edges.size(); i++) {
                sharp_edges_mat.row(i) << sharp_edges[i].x(), sharp_edges[i].y();
            }
            Eigen::MatrixXd sharp_nodes_mat(sharp_nodes.size(), 3);
            for (int i = 0; i < sharp_nodes.size(); i++) {
                sharp_nodes_mat.row(i) = sharp_nodes[i];
            }

            std::tie(sharp_start_pts, sharp_end_pts, sharp_color) =
                ConvertCurveNetWorkForRender(sharp_nodes_mat, sharp_edges_mat, Eigen::RowVector4d(.7, .7, .7, 1));

            SaveSegments(args::get(creaseSavePath) + "/sharp_crease.txt", sharp_start_pts, sharp_end_pts, sharp_color,
                         Eigen::RowVector3d(0, 0, 0), 1);
        } else {
            while (hook->appState->keepSolving) {
                // hook->warmup();
                hook->simulateOneStep();
                // hook->updateRenderGeometry();
            }
        }
    }

    return 0;
}
