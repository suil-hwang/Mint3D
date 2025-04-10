#include "ImGuiWidgets.h"
#include "imgui.h"
// #include "Polyscope.h"
#include "FieldView.h"
#include "polyscope/point_cloud.h"
#include "polyscope/polyscope.h"
#include "polyscope/surface_mesh.h"

#include "implot/implot.h"

#include "../MiNT3D/MiNTEnergy.h"

namespace ImGuiWidgets {

// Add these in to each function

// static void HelpMarker(const char* desc)
// {
//     ImGui::TextDisabled("(?)");
//     if (ImGui::BeginItemTooltip())
//     {
//         ImGui::PushTextWrapPos(ImGui::GetFontSize() * 35.0f);
//         ImGui::TextUnformatted(desc);
//         ImGui::PopTextWrapPos();
//         ImGui::EndTooltip();
//     }
// }

void ShowSolverParams(AppState& appState) {
    MiNT3D::SolveParams& w = appState.solve_params;

    // MiNT3D::SolveParams& w = *(appState.solve_params.get());

    // ImGui::InputDouble("S Perp Weight", &c.w_s_perp, 0.0f, 0.0f, "%e");
    ImGui::PushItemWidth(150);
    ImGui::InputInt("solve steps", &w.inner_iter_max_steps);
    ImGui::InputDouble("outer step size (mutiplicative)", &w.w_outer_step, 0.0f, 0.0f, "%e");

    // ImGui::InputDouble("Curl Weight", &w.w_mint, 0.0f, 0.0f, "%e");
    ImGui::InputDouble("max lambda", &w.w_max_lambda, 0.0f, 0.0f, "%e");
    ImGui::InputDouble("cur lambda", &w.lambda_penalty, 0.0f, 0.0f, "%e");

    ImGui::Text("presets:");
    if (ImGui::Button("set symmetric dirichlet")) {
        appState.solve_params = MiNT3D::SolveParams();
        appState.exp_name = "sym_dirichlet";
    }

    // ImGui::SameLine();

    if (ImGui::Button("odeco (precon)")) {
        appState.solve_params.setOdecoSolve();
        appState.exp_name = "odeco";
    }

    if (ImGui::Button("mesh (default)")) {
        appState.solve_params.setSoftOrthogSolve();
        appState.exp_name = "mesh";
    }

    if (ImGui::Button("fit current")) {
        appState.solve_params.setFitSolve();
        appState.exp_name = "fit";
    }

    if (ImGui::Button("align current")) {
        appState.solve_params.setSoftAlignSolve();
        appState.exp_name = "align";
    }

    if (ImGui::Button("align w SDF normal bound cond")) {
        appState.solve_params.setAlignWithNormalBoundarySolve();
        appState.exp_name = "align_w_dirch_boundary";
    }

    // ImGui::InputInt("total step number", &w.total_step);

    // use combed smoothness
    // use sym combing vs integrable combing
    // use sym smoothness

    // ImGui::NewLine();
    // ImGui::InputDouble("Boundary Penalty", &w.w_bound, 0.0f, 0.0f, "%e");
    // ImGui::InputDouble("Odeco (scaled jacobian)", &w.w_scaled_jacobian, 0.0f, 0.0f, "%e");
    // ImGui::InputDouble("Unit Norm", &w.w_unit_norm, 0.0f, 0.0f, "%e");
    // ImGui::InputDouble("Unit Barrier", &w.w_unit_barrier, 0.0f, 0.0f, "%e");

    // ImGui::InputDouble("Smoothness Weight", &w.w_smooth, 0.0f, 0.0f, "%e");
    ImGui::PopItemWidth();

    // bool* show_frames = &appState.show_frames;
    bool* b_smooth_combed = &w.b_smooth_combed;
    bool* b_smooth_sym = &w.b_smooth_sym;
    bool* b_smooth_asap_combing = &w.b_smooth_asap_combing;
    bool* b_int_combed = &w.b_int_combed;
    bool* b_int_sym = &w.b_int_sym;
    bool* b_orthog = &w.b_orthog;
    bool* b_unit = &w.b_unit;
    bool* b_use_kruskal_tensors_for_sym_smoothness = &w.b_use_kruskal_tensors_for_sym_smoothness;

    // std::string show_frames_checkbox = ("combed##cb");

    if (!ImGui::CollapsingHeader("Init")) {
        ImGui::Text("boundary conditions and :");

        // int* boundary_condition = &w.boundary_condition;   // 0 for "exact", 1 for "random"
        int boundary_condition = static_cast<int>(w.boundary_condition);   // Convert to int

        ImGui::RadioButton("free", &boundary_condition, 0);
        ImGui::SameLine();
        ImGui::RadioButton("SDF", &boundary_condition, 1);
        ImGui::SameLine();
        ImGui::RadioButton("poisson", &boundary_condition, 2);

        w.boundary_condition = static_cast<MiNT3D::BoundaryCondition>(boundary_condition);

        int boundary_hard_constraints = static_cast<int>(w.boundary_hard_constraints);   // Convert to int

        ImGui::RadioButton("none", &boundary_hard_constraints, 0);
        ImGui::SameLine();
        ImGui::RadioButton("normal", &boundary_hard_constraints, 1);
        ImGui::SameLine();
        ImGui::RadioButton("fixed_boundary", &boundary_hard_constraints, 2);

        w.boundary_hard_constraints = static_cast<MiNT3D::BoundaryHardConstraints>(boundary_hard_constraints);

        int init_state = static_cast<int>(w.init_state);

        // // Create radio buttons
        ImGui::PushItemWidth(150);
        ImGui::RadioButton("exact", &init_state, 0);
        ImGui::SameLine();
        ImGui::RadioButton("random", &init_state, 1);
        ImGui::SameLine();
        ImGui::InputDouble("scale", &w.w_init_scale, 0.0f, 0.0f, "%e");
        ImGui::PopItemWidth();

        w.init_state = static_cast<MiNT3D::InitState>(init_state);
    }

    if (!ImGui::CollapsingHeader("Smoothness")) {
        ImGui::Text("Run Info:");

        ImGui::Checkbox("combed", b_smooth_combed);
        ImGui::SameLine();
        ImGui::Checkbox("as_smooth_as_possible", b_smooth_asap_combing);
        ImGui::PushItemWidth(150);
        ImGui::InputDouble("w_smooth_combed", &w.w_smooth_combed, 0.0f, 0.0f, "%e");
        ImGui::SameLine();
        ImGui::InputDouble("w_smooth_combed_precon", &w.w_smooth_combed_precon, 0.0f, 0.0f, "%e");

        ImGui::PopItemWidth();

        ImGui::Checkbox("symmetric", b_smooth_sym);
        ImGui::SameLine();
        ImGui::Checkbox("use_kruskal_lift", b_use_kruskal_tensors_for_sym_smoothness);
        ImGui::PushItemWidth(150);
        ImGui::InputDouble("w_smooth_sym", &w.w_smooth_sym, 0.0f, 0.0f, "%e");
        ImGui::SameLine();
        ImGui::InputDouble("w_smooth_sym_precon", &w.w_smooth_sym_precon, 0.0f, 0.0f, "%e");

        ImGui::PopItemWidth();
    }

    if (!ImGui::CollapsingHeader("Integrability")) {
        ImGui::Text("Run Info:");
        ImGui::Checkbox("combed##int", b_int_combed);
        ImGui::PushItemWidth(150);
        ImGui::InputDouble("w_int_combed", &w.w_int_combed, 0.0f, 0.0f, "%e");
        ImGui::SameLine();
        ImGui::InputDouble("w_int_combed_fixed_weight", &w.w_int_combed_fixed_weight, 0.0f, 0.0f, "%e");
        ImGui::PopItemWidth();

        ImGui::Checkbox("symmetric##int", b_int_sym);
        ImGui::PushItemWidth(150);
        ImGui::InputDouble("w_int_sym", &w.w_int_sym, 0.0f, 0.0f, "%e");
        ImGui::SameLine();
        ImGui::InputDouble("w_int_sym_fixed_weight", &w.w_int_sym_fixed_weight, 0.0f, 0.0f, "%e");
        ImGui::PopItemWidth();
    }

    if (!ImGui::CollapsingHeader("Orthonormality")) {
        ImGui::Checkbox("init_as_odeco", &w.b_init_as_odeco);

        ImGui::Checkbox("orthogonality", b_orthog);

        ImGui::PushItemWidth(150);
        ImGui::InputDouble("w_orthog", &w.w_orthog, 0.0f, 0.0f, "%e");
        ImGui::SameLine();
        ImGui::InputDouble("w_orthog_precon", &w.w_orthog_precon, 0.0f, 0.0f, "%e");
        ImGui::SameLine();
        ImGui::InputDouble("w_orthog_constraint_weight", &w.w_orthog_constraint_weight, 0.0f, 0.0f, "%e");
        ImGui::PopItemWidth();

        ImGui::Checkbox("unit", b_unit);

        ImGui::PushItemWidth(150);
        ImGui::InputDouble("w_unit", &w.w_unit, 0.0f, 0.0f, "%e");
        ImGui::SameLine();
        ImGui::InputDouble("w_unit_precon", &w.w_unit_precon, 0.0f, 0.0f, "%e");
        ImGui::PopItemWidth();
    }

    if (!ImGui::CollapsingHeader("Fitting to initial frame field")) {
        ImGui::Checkbox("fit", &w.b_fit);

        const char* fit_type_items[] = {"Moments", "Direction", "Kruskal"};
        int current_fit_type = static_cast<int>(w.fit_type);   // Variable to hold the current selected FitType

        // Create the dropdown
        if (ImGui::Combo("Fit Type", &current_fit_type, fit_type_items, IM_ARRAYSIZE(fit_type_items))) {
            // Update the FitType based on the selected item
            w.fit_type = static_cast<MiNT3D::FitType>(current_fit_type);
        }

        ImGui::PushItemWidth(150);
        ImGui::InputDouble("w_fit", &w.w_fit, 0.0f, 0.0f, "%e");
        ImGui::SameLine();
        ImGui::InputDouble("w_fit_precon", &w.w_fit_precon, 0.0f, 0.0f, "%e");
        ImGui::PopItemWidth();

        ImGui::Checkbox("self_align", &w.b_self_align);
        ImGui::InputDouble("w_self_align", &w.w_self_align, 0.0f, 0.0f, "%e");

        ImGui::Checkbox("viscocity", &w.b_viscosity);
        ImGui::InputDouble("w_viscocity", &w.w_viscocity, 0.0f, 0.0f, "%e");

        ImGui::Checkbox("feat_align", &w.b_feat_align);
        ImGui::InputDouble("w_feat_align", &w.w_feat_align, 0.0f, 0.0f, "%e");
        // expose control on the order of the
    }
}

void ShowCombedFramesPerFacet(AppState& appState) {
    bool showOnlyFramesWithDiff = true;

    if (ImGui::Button("Show Fields with discrepancy")) {
        Eigen::MatrixXd frames = appState.frames;
        Eigen::MatrixXd bnd_frames = appState.boundary_frames;
        Eigen::MatrixXd V = appState.V;
        Eigen::MatrixXi T = appState.T;

        CubeCover::TetMeshConnectivity mesh(T);

        MiNT3D::MiNTEnergy energy(V, mesh, mesh);
        energy.Precompute();
        energy.ComputePermutations(frames);

        double rescaling_ratio = 2;

        Eigen::MatrixXd face_centroids = Eigen::MatrixXd::Zero(mesh.nFaces(), 3);
        Eigen::VectorXi incompatible_facets = Eigen::VectorXi::Zero(mesh.nFaces());

        for (int i = 0; i < mesh.nFaces(); i++) {
            Eigen::Vector3d v0 = V.row(mesh.faceVertex(i, 0));
            Eigen::Vector3d v1 = V.row(mesh.faceVertex(i, 1));
            Eigen::Vector3d v2 = V.row(mesh.faceVertex(i, 2));
            face_centroids.row(i) = (v0 + v1 + v2) / 3;
        }

        Eigen::MatrixXd fb0 = Eigen::MatrixXd::Zero(mesh.nFaces(), 3);
        Eigen::MatrixXd fb1 = Eigen::MatrixXd::Zero(mesh.nFaces(), 3);

        Eigen::MatrixXd t0v0 = Eigen::MatrixXd::Zero(mesh.nFaces(), 3);
        Eigen::MatrixXd t0v1 = Eigen::MatrixXd::Zero(mesh.nFaces(), 3);
        Eigen::MatrixXd t0v2 = Eigen::MatrixXd::Zero(mesh.nFaces(), 3);

        Eigen::MatrixXd t1v0_asap = Eigen::MatrixXd::Zero(mesh.nFaces(), 3);
        Eigen::MatrixXd t1v1_asap = Eigen::MatrixXd::Zero(mesh.nFaces(), 3);
        Eigen::MatrixXd t1v2_asap = Eigen::MatrixXd::Zero(mesh.nFaces(), 3);

        Eigen::MatrixXd t1v0_aiap = Eigen::MatrixXd::Zero(mesh.nFaces(), 3);
        Eigen::MatrixXd t1v1_aiap = Eigen::MatrixXd::Zero(mesh.nFaces(), 3);
        Eigen::MatrixXd t1v2_aiap = Eigen::MatrixXd::Zero(mesh.nFaces(), 3);

        for (int face_id = 0; face_id < mesh.nFaces(); face_id++) {
            int tet_id0 = mesh.faceTet(face_id, 0);
            int tet_id1 = mesh.faceTet(face_id, 1);
            if (tet_id0 == -1 || tet_id1 == -1) {
                continue;
            }

            // Extract the frames corresponding to the two tetrahedra
            Eigen::VectorXd v = frames.row(tet_id0);
            Eigen::VectorXd u = frames.row(tet_id1);

            Eigen::MatrixXd P_asap = energy.permutation_as_smooth_as_possible_[face_id];
            Eigen::MatrixXd P_aiap = energy.permutation_as_integrable_as_possible_[face_id];

            Eigen::VectorXd u_asap = P_asap * u;
            Eigen::VectorXd u_aiap = P_aiap * u;

            // Retrieve and expand the face basis
            Eigen::MatrixXd face_basis0 = energy.tet_facet_basis_.at(tet_id0).at(mesh.faceTetIndex(face_id, 0));

            double diff = (P_asap - P_aiap).norm();
            if (diff < 1e-6 && showOnlyFramesWithDiff) {
                continue;
            }

            incompatible_facets(face_id) = 1;

            fb0.row(face_id) = face_basis0.row(0);
            fb1.row(face_id) = face_basis0.row(1);

            t0v0.row(face_id) = v.segment(0, 3);
            t0v1.row(face_id) = v.segment(3, 3);
            t0v2.row(face_id) = v.segment(6, 3);

            t1v0_asap.row(face_id) = u_asap.segment(0, 3);
            t1v1_asap.row(face_id) = u_asap.segment(3, 3);
            t1v2_asap.row(face_id) = u_asap.segment(6, 3);

            t1v0_aiap.row(face_id) = u_aiap.segment(0, 3);
            t1v1_aiap.row(face_id) = u_aiap.segment(3, 3);
            t1v2_aiap.row(face_id) = u_aiap.segment(6, 3);
        }

        int nincompfacets = incompatible_facets.sum();
        Eigen::MatrixXd incompV(3 * nincompfacets, 3);
        Eigen::MatrixXi incompF(nincompfacets, 3);

        int iter = 0;

        for (int face_id = 0; face_id < mesh.nFaces(); face_id++) {
            if (incompatible_facets(face_id) == 0) {
                continue;
            }

            Eigen::Vector3d v0 = V.row(mesh.faceVertex(face_id, 0));
            Eigen::Vector3d v1 = V.row(mesh.faceVertex(face_id, 1));
            Eigen::Vector3d v2 = V.row(mesh.faceVertex(face_id, 2));

            incompV.row(3 * iter + 0) = v0;
            incompV.row(3 * iter + 1) = v1;
            incompV.row(3 * iter + 2) = v2;

            incompF.row(iter) = Eigen::Vector3i(3 * iter + 0, 3 * iter + 1, 3 * iter + 2);

            iter++;
        }

        auto* incompmesh = polyscope::registerSurfaceMesh("Incompatible Facets", incompV, incompF);
        incompmesh->setSurfaceColor({0.0, 0.0, 0.0});
        incompmesh->setTransparency(0.5);
        incompmesh->setEnabled(true);

        glm::vec3 c_r(0.9, 0.1, 0.1);
        glm::vec3 c_g(0.1, 0.9, 0.1);
        glm::vec3 c_b(0.1, 0.1, 0.9);

        polyscope::registerPointCloud("dual_edge_basis", face_centroids)->setPointRadius(0.0)->setEnabled(false);
        polyscope::getPointCloud("dual_edge_basis")
            ->addVectorQuantity("fb0", fb0)   //  polyscope::VectorType::AMBIENT)
            ->setEnabled(true);
        polyscope::getPointCloud("dual_edge_basis")
            ->addVectorQuantity("fb1", fb1)   // polyscope::VectorType::AMBIENT)
            ->setEnabled(true);

        polyscope::registerPointCloud("dual_edge_t0", face_centroids)->setPointRadius(0.0)->setEnabled(true);

        polyscope::getPointCloud("dual_edge_t0")
            ->addVectorQuantity("t0v0", t0v0)   //  polyscope::VectorType::AMBIENT)
            ->setVectorColor(c_r)
            ->setEnabled(true);
        polyscope::getPointCloud("dual_edge_t0")
            ->addVectorQuantity("t0v1", t0v1)   // polyscope::VectorType::AMBIENT)
            ->setVectorColor(c_r)
            ->setEnabled(true);
        polyscope::getPointCloud("dual_edge_t0")
            ->addVectorQuantity("t0v2", t0v2)   // polyscope::VectorType::AMBIENT)
            ->setVectorColor(c_r)
            ->setEnabled(true);

        polyscope::registerPointCloud("dual_edge_t1_asap", face_centroids)->setPointRadius(0.0)->setEnabled(true);

        polyscope::getPointCloud("dual_edge_t1_asap")
            ->addVectorQuantity("t1v0", t1v0_asap)
            ->setVectorColor(c_g)
            ->setEnabled(true);
        polyscope::getPointCloud("dual_edge_t1_asap")
            ->addVectorQuantity("t1v1", t1v1_asap)
            ->setVectorColor(c_g)
            ->setEnabled(true);
        polyscope::getPointCloud("dual_edge_t1_asap")
            ->addVectorQuantity("t1v2", t1v2_asap)
            ->setVectorColor(c_g)
            ->setEnabled(true);

        polyscope::registerPointCloud("dual_edge_t1_aiap", face_centroids)->setPointRadius(0.0)->setEnabled(true);

        polyscope::getPointCloud("dual_edge_t1_aiap")
            ->addVectorQuantity("t1v0", t1v0_aiap)
            ->setVectorColor(c_b)
            ->setEnabled(true);
        polyscope::getPointCloud("dual_edge_t1_aiap")
            ->addVectorQuantity("t1v1", t1v1_aiap)
            ->setVectorColor(c_b)
            ->setEnabled(true);
        polyscope::getPointCloud("dual_edge_t1_aiap")
            ->addVectorQuantity("t1v2", t1v2_aiap)
            ->setVectorColor(c_b)
            ->setEnabled(true);
    }
}

void ShowFrames(AppState& appState) {
    bool* show_frames_as_lines = &appState.show_frames_as_lines;
    bool* show_frames = &appState.show_frames;

    std::string show_frames_checkbox = ("draw vectors as lines##cb");
    ImGui::Checkbox("show frames", show_frames);
    ImGui::SameLine();
    ImGui::Checkbox(show_frames_checkbox.c_str(), show_frames_as_lines);

    // ImGui::InputDouble("Vec Size", &appState.gui_vec_size);

    ImGui::SliderFloat("Vec Size", &appState.gui_vec_size, 1e-4f, 10., "%.16f", ImGuiSliderFlags_Logarithmic);

    bool* recompute_singularities = &appState.recompute_singularities;
    ImGui::Checkbox("recompute singularities", recompute_singularities);

    ImGui::SameLine();

    bool* use_asap_combing_for_singularities = &appState.use_asap_combing_for_singularities;
    ImGui::Checkbox("use_asap_combing", use_asap_combing_for_singularities);
    // ImGui::InputDouble("L4 alpha", &appState.L4_alpha);
}

// Function to display a file scrubber in ImGui
// void ShowFileScrubber(int& fileIndex, int minIndex, int maxIndex) {
void ShowFileScrubber(AppState& appState) {
    ImGui::Text("File Scrubber:");

    // Display a slider to select the current file index

    int* fileIndex = &appState.currentFileID;
    int minIndex = 0;                                        // appState.minFileIndex;
    int maxIndex = std::max(appState.max_saved_index, 10);   // appState.fileList.size() - 1;
    int prevFileIndex = *fileIndex;

    // appState.currentFileIndex, 0, appState.fileList.size() - 1
    ImGui::PushItemWidth(250);
    if (ImGui::SliderInt("File Index", fileIndex, minIndex, maxIndex)) {
        std::cout << "file index changed" << std::endl;
        std::cout << "prev file index " << prevFileIndex << " | current file index " << *fileIndex << std::endl;
        appState.shouldReload = true;
    }
    ImGui::PopItemWidth();

    ImGui::SameLine();

    ImGui::PushButtonRepeat(true);
    if (ImGui::ArrowButton("##file_scrub_der", ImGuiDir_Left) && *fileIndex > minIndex) {
        appState.currentFileID--;
        appState.shouldReload = true;
    }
    ImGui::SameLine(0.0f, 5.0f);
    if (ImGui::ArrowButton("##file_scrub_incr", ImGuiDir_Right) && *fileIndex < maxIndex) {
        appState.currentFileID++;
        appState.shouldReload = true;
    }
    ImGui::PopButtonRepeat();

    if (appState.shouldReload && !appState.loadedPreviousRun) {
        appState.shouldLogData = false;
    }

    // if (ImGui::InputInt("bump cur file id",fileIndex))
    // {
    //     appState.shouldReload = true;
    // }
}

// Function to display checkboxes for field views in ImGui
void ShowFieldViewCheckboxes(AppState& appState) {
    ImGui::Text("File Logging Controls:");
    bool* log_state_to_file = &appState.shouldLogData;
    ImGui::Checkbox("log stuff to file", log_state_to_file);
    //       ImGui::SameLine();

    if (ImGui::CollapsingHeader("Select Which Views To Log:")) {
        for (int i = 0; i < (int)Views::Field_View::Element_COUNT - 1; ++i) {
            Field_View field = static_cast<Field_View>(i);
            // const char* fieldName = fieldViewToString(field).c_str();
            bool* active = &appState.fieldViewActive[(int)field];

            std::string field_str = fieldViewToFileStub(field);
            std::string checkbox = (field_str + "##cb");

            // std::cout << fieldName << std::endl;
            ImGui::Checkbox(checkbox.c_str(), active);
            // ImGui::SameLine();
        }
    }
}

// Function to display run information in ImGui
void ShowRunInfo(AppState& appState) {
    ImGui::Text("Run Info:");
    ImGui::Text("Experiment Type: %s", appState.exp_name.c_str());
    ImGui::Text("Logging Dir", appState.logFolderPath.c_str());

    // ImGui::Text(appState.solveDescription.c_str());
    ImGui::Text("Status: ");
    ImGui::SameLine();
    ImGui::Text(appState.solveStatus.c_str());

    ImGui::Text("");

    // // ImGui::Text("Current File: %s", appState.fileList[appState.currentFileIndex].c_str());
    // // ImGui::Text("Current Element: %s", fieldViewToString(appState.current_element).c_str());
    // ImGui::Text("Current Iteration: %d", appState.currentIteration);
    // ImGui::Text("Current File ID: %d", appState.currentFileID);

    // ImGui::Text("TODO");
    // ImGui::Text("Fix curl operator");
    // // ImGui::Text("refactor optzoo to use stencils");
    // // ImGui::Text("reset/load from folder");
    // //  ImGui::Text("file scrubber");
    // // ImGui::Text("offby1finished issue");
    // ImGui::Text("log energy in vector, graph here");

    // // ImGui::Text("Current Frame: %d", appState.currentFrameIndex);
    // // ImGui::Text("Current Moment: %d", appState.currentMomentIndex);

    // // Add other run information as needed
}

// Function to draw a file loader widget in ImGui
void DrawFileLoader(AppState& appState) {
    ImGui::Text("File Loader:");

    // ImGui::InputText("Directory", appState.directoryPath.c_str());
    if (ImGui::Button("Load Files")) {
        appState.refreshFileLists();
        // appState.setDefaultFileIndices();
    }
}

// Function to add field view scalars to Polyscope
void AddFieldViewScalarsToPolyscope(AppState& appState) {
    for (int i = 0; i < (int)Views::Field_View::Element_COUNT; ++i) {
        Field_View field = static_cast<Field_View>(i);
        std::string fieldName = fieldViewToFileStub(field);
        // std::vector<double> scalarValues(appState.fieldData[field].begin(), appState.fieldData[field].end());

        // Calculate bounds (10th and 90th percentiles)
        float minBound = appState.fieldBounds[field].lower;
        float maxBound = appState.fieldBounds[field].upper;

        // Add the scalar quantity with bounds to Polyscope
        // auto scalarQ = polyscope::curls_primal("my mesh")->addVertexScalarQuantity(fieldName, scalarValues);
        // scalarQ->setMapRange({minBound, maxBound});
    }
}

// Function to display field view checkboxes with sliders in ImGui
void ShowFieldViewCheckboxesWithSliders(AppState& appState) {
    ImGui::Text("Field Views with Sliders:");

    for (int i = 0; i < (int)Views::Field_View::Element_COUNT - 1; ++i) {
        Field_View field = static_cast<Field_View>(i);
        const char* fieldName = fieldViewToFileStub(field).c_str();
        bool* active = &appState.fieldViewActive[(int)field];
        float* minVal = &appState.fieldBounds[field].lower;
        float* maxVal = &appState.fieldBounds[field].upper;

        std::string field_str = fieldViewToFileStub(field);
        // const char* lower = ("l" + field_str).c_str();
        // const char* upper = ("u" + field_str).c_str();
        // const char* checkbox = ("cb" + field_str).c_str();

        // std::string lower = ("l" + field_str).c_str();
        // std::string upper = ("u" + field_str).c_str();
        // std::string checkbox = ("cb" + field_str).c_str();

        std::string lower = ("lower ##" + field_str);
        std::string upper = ("upper##" + field_str);
        std::string checkbox = (field_str + "##cb");

        // std::cout << "lower*****" << lower.c_str() << upper.c_str() << checkbox.c_str() << std::endl;

        ImGui::PushItemWidth(150);
        ImGui::InputFloat(lower.c_str(), minVal, 0.01f, .10f, "%.3f");
        ImGui::SameLine();
        ImGui::InputFloat(upper.c_str(), maxVal, 0.01f, .10f, "%.3f");
        ImGui::PopItemWidth();
        ImGui::SameLine();
        ImGui::Checkbox(checkbox.c_str(), active);

        //                    std::cout << fieldName << std::endl;

        // ImGui::SliderFloat((std::string(fieldName) + " (Lower Bound)").c_str(), minVal, 0.0f, 1.0f);
        // ImGui::SameLine();
        // ImGui::SliderFloat((std::string(fieldName) + " (Upper Bound)").c_str(), maxVal, 0.0f, 1.0f);
    }
}

void FieldViewSelector(AppState& appState, Field_View& currentField) {
    // Display a combo box to select the current field view
    if (ImGui::BeginCombo("Select Field View", fieldViewToString(currentField).c_str())) {
        for (int i = 0; i < (int)Views::Field_View::Element_COUNT; ++i) {
            Field_View field = static_cast<Field_View>(i);
            bool is_selected = (currentField == field);
            if (ImGui::Selectable(fieldViewToString(field).c_str(), is_selected)) {
                currentField = field;   // Update the current field view
            }
            if (is_selected) {
                ImGui::SetItemDefaultFocus();   // Set the default selection
            }
        }
        ImGui::EndCombo();
    }
    // ImGui::Text("Field String: %s", fieldViewToFileStub(currentField).c_str());

    switch (currentField) {
        case Field_View::vec_norms:
            // cur_scalar_quantity = outputData->norms_vec;
            break;
        case Field_View::delta_norms:
            // cur_scalar_quantity = outputData->norms_delta;
            break;
        case Field_View::vec_dirch:
            // cur_scalar_quantity = outputData->smoothness_primal;
            break;
        case Field_View::moment_dirch:
            // List box
            if (ImGui::TreeNode("Choose symmetric dirichlet component")) {
                Views::Sym_Moment_View selected = appState.cur_moment_view;

                if (ImGui::Selectable("L2_dirch", selected == Views::Sym_Moment_View::L2)) {
                    selected = Views::Sym_Moment_View::L2;
                }
                if (ImGui::Selectable("L4_dirch", selected == Views::Sym_Moment_View::L4)) {
                    selected = Views::Sym_Moment_View::L4;
                }
                if (ImGui::Selectable("L2 + L4_dirch", selected == Views::Sym_Moment_View::Total)) {
                    selected = Views::Sym_Moment_View::Total;
                }
                if (ImGui::Selectable("L6_dirch", selected == Views::Sym_Moment_View::L6)) {
                    selected = Views::Sym_Moment_View::L6;
                }
                if (appState.cur_moment_view != selected) {
                    appState.updateRenderGeometryNextFrameIfPaused = true;
                }
                appState.cur_moment_view = selected;

                ImGui::TreePop();
            }

            // cur_scalar_quantity = outputData->smoothness_sym;
            // cur_scalar_quantity = outputData->smoothness_L2;
            // cur_scalar_quantity = outputData->smoothness_L4;
            break;
        case Field_View::primal_curl_residual:
            // cur_scalar_quantity = outputData->curls_primal;
            break;
        case Field_View::sym_curl_residual:
            // cur_scalar_quantity = outputData->curls_sym;
            if (ImGui::TreeNode("Choose symmetric dirichlet component")) {
                Views::Sym_Curl_View selected = appState.cur_curl_view;
                if (ImGui::Selectable("Total Symmetric Curl Residual", selected == Views::Sym_Curl_View::Total)) {
                    selected = Views::Sym_Curl_View::Total;
                }
                if (ImGui::Selectable("L2_curl", selected == Views::Sym_Curl_View::L2)) {
                    selected = Views::Sym_Curl_View::L2;
                }
                if (ImGui::Selectable("L4_curl", selected == Views::Sym_Curl_View::L4)) {
                    selected = Views::Sym_Curl_View::L4;
                }
                if (ImGui::Selectable("L6_curl", selected == Views::Sym_Curl_View::L6)) {
                    selected = Views::Sym_Curl_View::L6;
                }
                if (appState.cur_curl_view != selected) {
                    appState.updateRenderGeometryNextFrameIfPaused = true;
                }
                appState.cur_curl_view = selected;

                ImGui::TreePop();
            }
        case Field_View::odeco:
            // cur_scalar_quantity = outputData->norms_vec;
            break;
        case Field_View::fradeco:
            // cur_scalar_quantity = outputData->norms_vec;
            break;

            break;
        case Field_View::gui_free:
            // Implement logic for gui_free if required
            break;
        default:
            std::cerr << "Unknown Field_View option selected in AppState." << std::endl;
            break;
    }
}

// Function to display a field view scrubber in ImGui
void ShowFieldViewScrubber(AppState& appState, Field_View& currentField) {
    // ImGui::Text("Field View Scrubber:");

    ImGui::Text("Scalar View Bounds:");

    bool* active = &appState.fieldViewActive[(int)currentField];
    float* minVal = &appState.fieldBounds[currentField].lower;
    float* maxVal = &appState.fieldBounds[currentField].upper;

    std::string field_str = fieldViewToFileStub(currentField);

    std::string percentile_lower = ("lower_bound##percentile" + field_str);
    std::string percentile_upper = ("upper_bound##percentile" + field_str);

    ImGui::PushItemWidth(150);
    ImGui::InputFloat(percentile_lower.c_str(), minVal, 0.01f, .10f, "%.3f");
    ImGui::SameLine();
    ImGui::InputFloat(percentile_upper.c_str(), maxVal, 0.01f, .10f, "%.3f");
    ImGui::PopItemWidth();

    std::string override_lower = ("lower_bound##ovr" + field_str);
    std::string overide_upper = ("upper_bound##ovr" + field_str);

    bool* ovr_active = &appState.override_bounds_active;
    float* ovr_minVal = &appState.override_bounds.lower;
    float* ovr_maxVal = &appState.override_bounds.upper;

    std::string checkbox = ("abs override##ovr");

    ImGui::PushItemWidth(150);

    if (*ovr_active) {
        ImGui::Text("vvv Absolute Bounds vvv");
    } else {
        ImGui::Text("^^^ Percentile Bounds ^^^");
    }

    ImGui::SameLine();
    // std::cout << fieldName << std::endl;
    ImGui::Checkbox(checkbox.c_str(), ovr_active);
    ImGui::InputFloat(override_lower.c_str(), ovr_minVal, 0.001f, .10f, "%.8f");
    ImGui::SameLine();
    ImGui::InputFloat(overide_upper.c_str(), ovr_maxVal, 0.001f, .10f, "%.8f");
    ImGui::PopItemWidth();

    ImGui::PushItemWidth(300);
    ImGui::SliderFloat("override max (log slider)", ovr_maxVal, 1e-16f, 1e-3f, "%.16f", ImGuiSliderFlags_Logarithmic);
    ImGui::PopItemWidth();

    // //  const char* element_names[Field_View::Element_COUNT] = { "Vector Norms", "Delta Norms", "Vector
    // Dirichlet", "Symmetric Dirichlet", "Vector Curl", "Symmetric Curl", "free" };
    //         // const char* current_element_name = (current_element >= 0 && current_element <
    //         Field_View::Element_COUNT) ? element_names[current_element] : "Unknown"; const char*
    //         current_element_name = fieldViewToString(currentField).c_str(); ImGui::PushItemWidth(300);
    //         ImGui::SliderInt("Shading Mode", (int *) currentField, 0, Views::Field_View::Element_COUNT - 1,
    //         current_element_name); ImGui::PopItemWidth();
}

void SetSeabornStyle() {
    ImPlot::CreateContext();
    ImPlotStyle& style = ImPlot::GetStyle();

    // Set the background color to white
    style.Colors[ImPlotCol_PlotBg] = ImVec4(1.0f, 1.0f, 1.0f, 1.0f);

    // Set the plot border color to a light gray
    style.Colors[ImPlotCol_PlotBorder] = ImVec4(0.8f, 0.8f, 0.8f, 1.0f);

    // Set the colormap to ImPlotColormap_Deep
    // ImPlot::SetColormap(ImPlotColormap_Deep);

    // // Set the line color to a seaborn-like blue
    // style.Colors[ImPlotCol_Line] = ImVec4(0.42f, 0.67f, 0.84f, 1.0f);

    // style.Colors[ImPlotCol_TitleBg] = ImVec4(1.0f, 1.0f, 1.0f, 1.0f);  // Opaque white

    // Set the line weight
    style.LineWeight = 2.0f;

    // Set the marker size
    style.MarkerSize = 4.0f;

    // Set the color of the x and y axis to black
    style.Colors[ImPlotCol_AxisBg] = ImVec4(0.0f, 0.0f, 0.0f, 1.0f);

    // Set the color of the grid to a light gray
    style.Colors[ImPlotCol_AxisGrid] = ImVec4(0.8f, 0.8f, 0.8f, 1.0f);

    // Set the color of the legend background to white
    style.Colors[ImPlotCol_LegendBg] = ImVec4(1.0f, 1.0f, 1.0f, 1.0f);

    // Set the color of the legend border to a light gray
    style.Colors[ImPlotCol_LegendBorder] = ImVec4(0.8f, 0.8f, 0.8f, 1.0f);

    // Set the color of the legend text to black
    style.Colors[ImPlotCol_LegendText] = ImVec4(0.0f, 0.0f, 0.0f, 1.0f);
}

// void SetLightStyle() {
//     ImPlot::CreateContext();
//     ImPlotStyle& style = ImPlot::GetStyle();

//     // Set the background color to white
//     style.Colors[ImPlotCol_PlotBg] = ImVec4(1.0f, 1.0f, 1.0f, 1.0f);

//     // Set the plot border color to a light gray
//     style.Colors[ImPlotCol_PlotBorder] = ImVec4(0.8f, 0.8f, 0.8f, 1.0f);

//     // Set the line weight
//     style.LineWeight = 2.0f;

//     // Set the marker size
//     style.MarkerSize = 4.0f;

//     // Set the color of the x and y axis to black
//     style.Colors[ImPlotCol_XAxis] = ImVec4(0.0f, 0.0f, 0.0f, 1.0f);
//     style.Colors[ImPlotCol_YAxis] = ImVec4(0.0f, 0.0f, 0.0f, 1.0f);

//     // Set the color of the grid to a light gray
//     style.Colors[ImPlotCol_XGrid] = ImVec4(0.8f, 0.8f, 0.8f, 1.0f);
//     style.Colors[ImPlotCol_YGrid] = ImVec4(0.8f, 0.8f, 0.8f, 1.0f);

//     // Set the color of the legend background to white
//     style.Colors[ImPlotCol_LegendBg] = ImVec4(1.0f, 1.0f, 1.0f, 1.0f);

//     // Set the color of the legend border to a light gray
//     style.Colors[ImPlotCol_LegendBorder] = ImVec4(0.8f, 0.8f, 0.8f, 1.0f);

//     // Set the color of the legend text to black
//     style.Colors[ImPlotCol_LegendText] = ImVec4(0.0f, 0.0f, 0.0f, 1.0f);

//         ImPlot::SetColormap(ImPlotColormap_Deep);
// }

void ShowPlots(AppState& appState) {
    ImPlot::CreateContext();

    ImPlotStyle& style = ImPlot::GetStyle();
    // style.Colors[ImPlotCol_PlotBg] = ImVec4(1.0f, 1.0f, 1.0f, 1.0f);  // Opaque white
    // style.Colors[ImPlotCol_Line] = ImVec4(0.42f, 0.67f, 0.84f, 1.0f);  // Change line color
    style.LineWeight = 1.75f;   // Change line weight

    // For generating paper screenshots:
    // SetSeabornStyle();

    // static double xs[1001], ys1[1001], ys2[1001], ys3[1001];
    //     for (int i = 0; i < 1001; ++i) {
    //         xs[i]  = i*0.1f;
    //         ys1[i] = sin(xs[i]) + 1;
    //         ys2[i] = log(xs[i]);
    //         ys3[i] = pow(10.0, xs[i]);
    //     }

    if (ImPlot::BeginPlot("Top Level", ImVec2(-1, 0))) {
        // conf.values.ys = ;
        // conf.values.count = appState.energy_trace.size();

        int iter_count = appState.total_energy_log.size();
        // ImPlot::SetupAxisScale(ImAxis_X1, ImPlotScale_Log10);
        ImPlot::SetupLegend(ImPlotLocation_NorthEast, ImPlotLegendFlags_Outside | ImPlotLegendFlags_Horizontal);
        ImPlot::SetupAxisScale(ImAxis_Y1, ImPlotScale_Log10);

        ImPlot::SetupAxesLimits(0, iter_count, 1e-23, 1e17);

        // ImPlot::PlotLine("total_time", &(appState.total_time[0]), iter_count);
        // ImPlot::PlotLine("assembly_time", &(appState.assembly_time[0]), iter_count);
        // ImPlot::PlotLine("solve_time", &(appState.solve_time[0]), iter_count);
        // ImPlot::PlotLine("line_search_time", &(appState.line_search_time[0]), iter_count);

        // ImPlot::PlotLine("energy_diff", &(appState.energy_diff_log[0]), iter_count);
        // ImPlot::PlotLine("solve_residual", &(appState.solve_residual_log[0]), iter_count);

        ImPlot::PlotLine("total_energy", &(appState.total_energy_log[0]), iter_count);

        ImPlot::PlotLine("primal_integrability_log", &(appState.primal_integrability_log[0]), iter_count);

        // ImPlot::PlotLine("curl_weight_log", &(appState.curl_weight_log[0]), iter_count);

        ImPlot::PlotLine("global_scale_outiter_log", &(appState.global_scale_outiter_log[0]), iter_count);

        if (appState.solve_params.b_fit) {
            ImPlot::PlotLine("fit_alignment_penalty", &(appState.identity_weight_log[0]), iter_count);
        }

        // ImPlot::PlotLine("objective val", &(appState.energy_trace[0]), iter_count);
        // ImPlot::PlotLine("smoothness part", &(appState.energy_smoothness_part_trace[0]), iter_count);
        // ImPlot::PlotLine("curl part (w/o attenuation)", &(appState.energy_curl_part_trace[0]), iter_count);
        // ImPlot::PlotLine("energy_odeco_part", &(appState.energy_odeco_part_trace[0]), iter_count);
        // ImPlot::PlotLine("energy_fradeco_part", &(appState.energy_fradeco_part_trace[0]), iter_count);
        // ImPlot::PlotLine("energy_bound_vol_part", &(appState.energy_bound_vol_part_trace[0]), iter_count);
        // ImPlot::PlotLine("relative solver residual", &(appState.solve_rel_residual_trace[0]), iter_count);
        // ImPlot::PlotLine("max gradient norm", &(appState.cur_max_gradient_norm_trace[0]), iter_count);
        // ImPlot::PlotLine("Identity Weight", &(appState.identity_weight_trace[0]), iter_count);
        // // ImPlot::PlotLine("Other Parts Energy", &(appState.energy_other_parts_trace[0]), iter_count);
        // ImPlot::PlotLine("Timing", &(appState.timing_trace[0]), iter_count);

        // ImPlot::PlotLine("f(x) = x", &(appState.energy_trace[0]), iter_count);
        // ImPlot::PlotLine("f(x) = sin(x)+1", xs, ys1, 1001);
        // ImPlot::PlotLine("f(x) = log(x)",   xs, ys2, 1001);
        // ImPlot::PlotLine("f(x) = 10^x",     xs, ys3, 21);
        ImPlot::EndPlot();
    }

    if (ImPlot::BeginPlot("Log Plot", ImVec2(-1, 0))) {
        // conf.values.ys = ;
        // conf.values.count = appState.energy_trace.size();

        int inner_iter_count = appState.total_time.size();
        // ImPlot::SetupAxisScale(ImAxis_X1, ImPlotScale_Log10);
        ImPlot::SetupLegend(ImPlotLocation_NorthEast, ImPlotLegendFlags_Outside | ImPlotLegendFlags_Horizontal);
        ImPlot::SetupAxisScale(ImAxis_Y1, ImPlotScale_Log10);

        ImPlot::SetupAxesLimits(0, inner_iter_count, 1e-23, 1e8);

        ImPlot::PlotLine("gradient_norm_log", &(appState.gradient_norm_log[0]), inner_iter_count);
        ImPlot::PlotLine("gradient_norm_step_start_log", &(appState.gradient_norm_step_start_log[0]), inner_iter_count);

        ImPlot::PlotLine("energy_diff", &(appState.energy_diff_log[0]), inner_iter_count);

        ImPlot::PlotLine("solve_residual", &(appState.solve_residual_log[0]), inner_iter_count);

        ImPlot::PlotLine("global_scale_log", &(appState.global_scale_log[0]), inner_iter_count);
        ImPlot::PlotLine("identity_weight", &(appState.identity_weight_log[0]), inner_iter_count);

        ImPlot::EndPlot();
    }

    if (ImPlot::BeginPlot("Combed Smoothness", ImVec2(-1, 0))) {
        // conf.values.ys = ;
        // conf.values.count = appState.energy_trace.size();

        int iter_count = appState.total_energy_log.size();
        // ImPlot::SetupAxisScale(ImAxis_X1, ImPlotScale_Log10);
        ImPlot::SetupLegend(ImPlotLocation_NorthEast, ImPlotLegendFlags_Outside | ImPlotLegendFlags_Horizontal);
        ImPlot::SetupAxisScale(ImAxis_Y1, ImPlotScale_Log10);

        ImPlot::SetupAxesLimits(0, iter_count, 1e-8, 1e3);

        ImPlot::PlotLine("aiap_smoothness_log", &(appState.aiap_combed_smoothness_log[0]), iter_count);
        ImPlot::PlotLine("asap_smoothness_log", &(appState.asap_combed_smoothness_log[0]), iter_count);

        ImPlot::PlotLine("primal_integrability_log", &(appState.primal_integrability_log[0]), iter_count);

        ImPlot::EndPlot();
    }

    if (ImPlot::BeginPlot("Other stuff Plot", ImVec2(-1, 0))) {
        // conf.values.ys = ;
        // conf.values.count = appState.energy_trace.size();

        int iter_count = appState.total_energy_log.size();
        // ImPlot::SetupAxisScale(ImAxis_X1, ImPlotScale_Log10);
        ImPlot::SetupLegend(ImPlotLocation_NorthEast, ImPlotLegendFlags_Outside | ImPlotLegendFlags_Horizontal);
        ImPlot::SetupAxisScale(ImAxis_Y1, ImPlotScale_Log10);

        ImPlot::SetupAxesLimits(0, iter_count, 1e-10, 1e8);

        ImPlot::PlotLine("total_energy", &(appState.total_energy_log[0]), iter_count);
        ImPlot::PlotLine("sym_smooth", &(appState.smoothness_log[0]), iter_count);
        ImPlot::PlotLine("sym_int", &(appState.symmetric_integrability_log[0]), iter_count);

        if (appState.solve_params.b_fit) {
            ImPlot::PlotLine("fit_alignment_penalty", &(appState.identity_weight_log[0]), iter_count);
        }

        ImPlot::PlotLine("unit_penalty_log", &(appState.unit_penalty_log[0]), iter_count);
        ImPlot::PlotLine("orthogonality_log", &(appState.scaled_jacobian_log[0]), iter_count);

        ImPlot::EndPlot();
    }

    // if (ImPlot::BeginPlot("Curl Plot", ImVec2(-1, 0))) {
    //     // conf.values.ys = ;
    //     // conf.values.count = appState.energy_trace.size();

    //     int iter_count = appState.energy_trace.size();
    //     // ImPlot::SetupAxisScale(ImAxis_X1, ImPlotScale_Log10);
    //     ImPlot::SetupLegend(ImPlotLocation_SouthWest, ImPlotLegendFlags_Outside);
    //     ImPlot::SetupAxisScale(ImAxis_Y1, ImPlotScale_Log10);

    //     ImPlot::SetupAxesLimits(0, iter_count, 1e-20, 1e5);
    //     // ImPlot::PlotLine("curl part (w/o attenuation)", &(appState.energy_curl_part_trace[0]), iter_count);

    //     ImPlot::EndPlot();
    // }

    ImPlot::DestroyContext();

    // ImGui::PlotConfig conf;
    // // conf.values.xs = x_data; // this line is optional
    // // conf.values.ys = y_data;
    // conf.values.ys = &(appState.energy_trace[0]);
    // conf.values.count = appState.energy_trace.size();

    // conf.scale.min = 0;
    // conf.scale.max = 1e1;
    // // conf.scale.type = conf.scale.Log10;

    // conf.tooltip.show = true;
    // conf.tooltip.format = "x=%.2f, y=%.2f";
    // conf.grid_x.show = true;
    // conf.grid_y.show = true;
    // conf.frame_size = ImVec2(400, 400);
    // conf.line_thickness = 2.f;

    // ImGui::Plot("plot", conf);
}

// // Function to display a plot in ImGui
// void ShowPlot(const char* label, const std::vector<float>& values, float minY, float maxY) {
//     ImGui::Text("Plot: %s", label);
//     ImGui::PlotLines(label, &values[0], static_cast<int>(values.size()), 0, NULL, minY, maxY, ImVec2(0, 80));
// }
// ImGui::PlotLines("Frame Times", arr, IM_ARRAYSIZE(arr));

// void if(appState.current_element == )

void ShowMainGUI(AppState& appState) {
    // Begin the main ImGui window
    // ImGui::Begin("Main Interface");

    // // Display checkboxes with min and max sliders for field views
    // // ShowFieldViewCheckboxesWithSliders(appState);

    // // Display file scrubber for selecting files
    // ShowFileScrubber(appState);

    // ImGui::Begin("Select Current Scalar Field");
    // FieldViewSelector(appState, appState.current_element);
    // ImGui::End();

    // // Display a field view scrubber
    // ShowFieldViewScrubber(appState, appState.current_element);

    // ShowFieldViewCheckboxes(appState);

    // DrawFileLoader(appState);

    ImGui::Begin("Solver Params");
    ShowSolverParams(appState);
    ImGui::End();

    ImGui::Begin("Frame Display Options");
    ShowFrames(appState);
    ImGui::End();

    // Display run information
    ImGui::Begin("Optimization State");
    ShowRunInfo(appState);
    ImGui::End();

    ImGui::Begin("Convergence Plots");
    ShowPlots(appState);
    ImGui::End();

    // ShowPlot

    ImGui::ShowDemoWindow();   // TODO remove this

    // Display a plot for run_step_times
    // ShowPlot(appState.run_step_times, "Step Times", 0.0f, 100.0f);

    // Display a plot for run_energy
    // ShowPlot(appState.run_energy, "Energy", 0.0f, 100.0f);

    // Additional GUI elements or widgets can be added here

    // End the main ImGui window
    // ImGui::End();
}

}   // namespace ImGuiWidgets