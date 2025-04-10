
///////////////////////
/////   init simulation
///////////////////////

#include "Mint3DHook.h"
#include <chrono>
#include <filesystem>
#include <iostream>
#include <sstream>
#include "AppState.h"
#include "FileParser.h"
#include "Serialization.h"

#include "ImGuiWidgets.h"

#include "polyscope/curve_network.h"
#include "polyscope/point_cloud.h"
#include "polyscope/polyscope.h"
#include "polyscope/surface_mesh.h"
#include "polyscope/volume_mesh.h"

#include "../MiNT3D/MiNTModel.h"
#include "../MiNT3D/MiNTSolveParams.h"

// #include <igl/readOBJ.h>
// #include "CubeCover/readMesh.h"
#include "CubeCover/ReadFrameField.h"
#include "CubeCover/ReadVTK.h"
#include "CubeCover/readMeshFixed.h"

#include <igl/on_boundary.h>
#include <igl/writeOBJ.h>
// #include <Eigen/SparseLU>
#include <igl/boundary_facets.h>
#include <igl/cotmatrix.h>
#include <igl/massmatrix.h>
#include <igl/remove_unreferenced.h>
#include <igl/signed_distance.h>

#include <Eigen/CholmodSupport>

#include "CubeCover/FrameField.h"
#include "CubeCover/SingularCurveNetwork.h"
#include "CubeCover/WriteFrameField.h"
#include "Surface.h"
#include "UtilsMisc.h"
#include "date.h"

#include <thread>

#include <cstdlib>

#include <functional>   // Add this line

#if defined(_WIN32) || defined(_WIN64)
// Windows-specific includes
#    include <winsock.h>
#elif defined(__linux__)
// Linux-specific includes
#    include <unistd.h>
#endif

namespace fs = std::filesystem;

void Mint3DHook::drawGUI() {
    if (!appState->headless_mode) {
        ImGuiWidgets::ShowMainGUI(*appState);
    }
}

void Mint3DHook::updateRenderGeometry() {
    std::cout << "update render geometry" << std::endl;

    if (appState->shouldReload) {
        std::cout << "do from file reload" << std::endl;
        pause();
        if (!this->isPaused()) {
            std::cout << "failed to pause" << std::endl;
        }
        initSimulation();
    }

    if (appState->headless_mode) {
        return;
    }

    int frows = appState->frames.rows();
    int bfrows = appState->boundary_frames.rows();

    int vec_dofs = 3;

    if (vec_dofs != 3) std::cout << "******** wrong number of dofs for mint3d ********" << std::endl;

    // outputData->frames = appState->frames;

    int frame_rank = 3;
    for (int vid = 0; vid < frame_rank; vid++) {
        Eigen::MatrixXd vec_cur = Eigen::MatrixXd::Zero(appState->frames.rows(), 3);

        // vec_cur << appState->frames.block(0, vid * vec_dofs, frows, vec_dofs); // ,  Eigen::MatrixXd::Zero(frows, 1);
        // outputData->frames.push_back(vec_cur);
    }

    for (int vid = 0; vid < frame_rank; vid++) {
        Eigen::MatrixXd vec_cur = Eigen::MatrixXd::Zero(appState->boundary_frames.rows(), 3);
        // vec_cur << appState->boundary_frames.block(0, vid * vec_dofs, bfrows, vec_dofs); // ,
        // Eigen::MatrixXd::Zero(frows, 1); outputData->boundary_frames.push_back(vec_cur);
    }

    if (appState->shouldReload) {
        appState->currentFileID--;
        appState->shouldReload = false;
    }

    // appState->LogToFile("curr");

    // appState->zeroPassiveVars();

    // Request a redraw in Polyscope to update the visualization

    if (!appState->headless_mode && appState->frames.rows() > 0) {
        if (appState->recompute_singularities) {
            this->recomputeSingularStructures();
        }

        polyscope::requestRedraw();
    }
}

void Mint3DHook::recomputeSingularStructures() {
    Eigen::MatrixXi assignments;

    // Call the function to get the raw pointer
    CubeCover::FrameField *rawPointer =
        CubeCover::fromFramesAndAssignments(*(appState->cur_tet_mesh), appState->frames, assignments, true);

    appState->field = std::unique_ptr<CubeCover::FrameField>(rawPointer);

    // if (!field)
    //     return -1;

    // if (recomputeperms)
    // {
    //     std::cout << "No face assignments provided, recomputing: ";
    //     std::cout.flush();

    if (appState->use_asap_combing_for_singularities) {
        appState->field->computeLocalAssignments();
    } else {
        appState->field->computeLocalAssignments(&appState->V);
    }

    std::cout << "found " << appState->field->nSingularEdges() << " singular edges" << std::endl;
    // }
    appState->field->combAssignments();

    extractSingularCurveNetwork(appState->V, *(appState->cur_tet_mesh), *(appState->field.get()), appState->Pgreen,
                                appState->Egreen, appState->Pblue, appState->Eblue, appState->Pblack, appState->Eblack);
}

void Mint3DHook::updateFrameVisualization() {
    int num_vecs = 3;
    for (int v = 0; v < num_vecs; v++) {
        Eigen::MatrixXd cur_vec = appState->frames.block(0, v * 3, appState->frames.rows(), 3);
        cur_vec *= appState->gui_vec_size;

        double color_shift = (v + 1.) * 1.0 / num_vecs;

        auto vectorField = polyscope::getPointCloud("c_vecs")->addVectorQuantity(
            "Vector Field " + std::to_string(v), cur_vec, polyscope::VectorType::AMBIENT);
        auto vectorFieldNeg = polyscope::getPointCloud("c_vecs")->addVectorQuantity(
            "Vector Field (negative) " + std::to_string(v), (-1.) * cur_vec, polyscope::VectorType::AMBIENT);

        if (appState->show_frames && appState->show_frames_as_lines) {
            vectorField->setEnabled(true);
            vectorFieldNeg->setEnabled(true);
            vectorField->setVectorLengthScale(appState->gui_vec_size);
            vectorFieldNeg->setVectorLengthScale(appState->gui_vec_size);

        } else if (appState->show_frames) {
            vectorField->setEnabled(true);
            vectorFieldNeg->setEnabled(false);
        } else {
            vectorField->setEnabled(false);
            vectorFieldNeg->setEnabled(false);
        }
    }

    for (int v = 0; v < num_vecs; v++) {
        Eigen::MatrixXd cur_vec = appState->boundary_frames.block(0, v * 3, appState->boundary_frames.rows(), 3);
        cur_vec *= appState->gui_vec_size;

        double color_shift = (v + 1.) * 1.0 / num_vecs;

        auto vectorField =
            polyscope::getPointCloud("surface_vecs")
                ->addVectorQuantity("Vector Field " + std::to_string(v), cur_vec, polyscope::VectorType::AMBIENT);
        auto vectorFieldNeg = polyscope::getPointCloud("surface_vecs")
                                  ->addVectorQuantity("Vector Field (negative) " + std::to_string(v), (-1.) * cur_vec,
                                                      polyscope::VectorType::AMBIENT);

        if (appState->show_frames && appState->show_frames_as_lines) {
            vectorField->setEnabled(true);
            vectorFieldNeg->setEnabled(true);
            vectorField->setVectorLengthScale(appState->gui_vec_size);
            vectorFieldNeg->setVectorLengthScale(appState->gui_vec_size);

        } else if (appState->show_frames) {
            vectorField->setEnabled(true);
            vectorFieldNeg->setEnabled(false);
        } else {
            vectorField->setEnabled(false);
            vectorFieldNeg->setEnabled(false);
        }
    }
}

void Mint3DHook::updateSingularityVisualiation() {
    auto *green = polyscope::registerCurveNetwork("Singular Curves (+1/4)", appState->Pgreen, appState->Egreen);
    green->setColor({0.0, 1.0, 0.0});

    auto *blue = polyscope::registerCurveNetwork("Singular Curves (-1/4)", appState->Pblue, appState->Eblue);
    blue->setColor({0.0, 0.0, 1.0});

    auto *black = polyscope::registerCurveNetwork("Singular Curves (irregular)", appState->Pblack, appState->Eblack);
    black->setColor({0.0, 0.0, 0.0});

    if (appState->recompute_singularities) {
        green->setEnabled(true);
        blue->setEnabled(true);
        black->setEnabled(true);
    } else {
        green->setEnabled(false);
        blue->setEnabled(false);
        black->setEnabled(false);
    }

    auto *features = polyscope::registerCurveNetwork("Sharp Features", appState->sharp_nodes, appState->sharp_edges);
    features->setColor({0.7, 0.1, 0.1});
    features->setEnabled(true);

    // auto *ps =
    //     polyscope::registerSurfaceMesh("Isosurfaces", appState->cur_surf->data().V, appState->cur_surf->data().F);
}

// Need to fill out viewer for each of: Field_View { vec_dirch, moment_dirch, sym_curl_residual, primal_curl_residual,
void Mint3DHook::renderRenderGeometry() {
    // still do this I guess but return as a no-op in headless mode
    // if (appState->headless_mode)
    // {
    //     return;
    // }

    if ((appState->shouldReload || appState->updateRenderGeometryNextFrameIfPaused) && this->isPaused()) {
        updateRenderGeometry();
        appState->updateRenderGeometryNextFrameIfPaused = false;
    }

    // Depending on the current element view, render different quantities
    const std::string cur_field = fieldViewToFileStub(appState->current_element);

    if (appState->current_element == Field_View::gui_free) {
        // noop
    } else if (appState->showVectorField && appState->frames.rows() > 0) {
        updateFrameVisualization();
        updateSingularityVisualiation();
    }

    polyscope::requestRedraw();
}

void Mint3DHook::pause() {
    PhysicsHook::pause();
    appState->keepSolving = true;
    appState->solveStatus = "paused";
    std::cout << "paused" << std::endl;
}
// Pause the simulation

void Mint3DHook::initSimulation() {
    appState->keepSolving = true;
    appState->solveStatus = "init simulation";

    // bool create_new_dir = false;
    // if (appState->directoryPath.empty()) {
    //     std::cerr << "No directory path provided. Creating new directory." << std::endl;
    //     create_new_dir = true;
    //     appState->directoryPath = "../../results/BLAH"; // switch to ../shared/mint2d_testsequence
    //     // return;
    // }

    Eigen::MatrixXd V;   // Temporary storage for vertices
    Eigen::MatrixXi T;   // Temporary storage for tetrahedra
    Eigen::MatrixXi F;   // Temporary storage for faces

    if (!CubeCover::readMESH(appState->meshInputFilePath, V, T, F)) {
        std::cerr << "Failed to load mesh from " << appState->meshInputFilePath << std::endl;
        return;
    }
    // else
    // if (!CubeCover::readVTK(appState->meshInputFilePath, V, T)) {
    //     std::cerr << "Failed to load mesh from " << appState->meshInputFilePath << std::endl;
    //     return;
    // }

    appState->cur_tet_mesh = std::make_unique<CubeCover::TetMeshConnectivity>(T);

    // Resize to unit bounding box

    Eigen::Vector3d bbox_min = V.row(0);
    Eigen::Vector3d bbox_max = V.row(0);
    int nverts = V.rows();

    for (int i = 0; i < nverts; i++) {
        Eigen::Vector3d c = V.row(i);
        for (int j = 0; j < 3; j++) {
            if (c(j) < bbox_min(j)) bbox_min(j) = c(j);

            if (c(j) > bbox_max(j)) bbox_max(j) = c(j);
        }
    }

    Eigen::Vector3d diff = bbox_max - bbox_min;
    double max_diff = diff.maxCoeff();

    // Calculate the center of the bounding box
    Eigen::Vector3d bbox_center = (bbox_min + bbox_max) / 2.0;

    // Calculate the shift needed to move the center to the origin
    Eigen::Vector3d shift = -bbox_center;

    // Scale factor to resize to unit bounding box
    double scale = 1.0 / max_diff;

    // shift << 1., 2., 3.;
    for (int i = 0; i < nverts; i++) {
        V.row(i) += shift;
        V.row(i) *= scale;
    }

    std::cout << "V.rows() " << V.rows() << " T.rows() " << T.rows() << " F.rows() " << F.rows() << std::endl;

    // Set mesh data to AppState
    appState->V = V;
    appState->T = T;

    // hack
    // appState->subdivide_boundary_tets = false;
    if (appState->subdivide_boundary_tets) {
        subdivideMeshBoundary();
        subdivideLockedFacets();
        T = appState->T;
        V = appState->V;
    }

    // get surface faces
    // MUST COME AFTER SUBDIVISION!!!
    int nfaces = appState->cur_tet_mesh->nBoundaryElements();
    F = Eigen::MatrixXi(nfaces, 3);
    for (int i = 0; i < nfaces; i++) {
        for (int j = 0; j < 3; j++) {
            F(i, j) = appState->cur_tet_mesh->faceVertex(appState->cur_tet_mesh->boundaryFace(i), j);
        }
    }
    appState->F = F;

    auto surface_data = GetBoundarySurfaceMeshFromTetMesh(V, *appState->cur_tet_mesh);
    appState->cur_surf = std::make_unique<Surface>(surface_data.first, surface_data.second);

    appState->nelem = T.rows();
    if (appState->useBoundaryFrames) {
        appState->nelem += appState->cur_tet_mesh->nBoundaryElements();
    }
    appState->ntets = T.rows();
    appState->nbelem = appState->cur_tet_mesh->nBoundaryElements();

    findSharpFeatures();

    if (appState->currentFileID == -1 && appState->logFolderPath != "") {
        //         appState->readAllLogFiles();
        //         appState->currentFileID = appState->max_saved_index;
        //         std::cout << "appState->currentFileID " << appState->currentFileID << std::endl;
        //         std::cout << "appState->max_saved_index " << appState->max_saved_index << std::endl;)
    }

    this->resetAppState();

    if (appState->currentFileID == -1) {
        // call init stuff here?
        // maybe add some flags for this.
        appState->frames = Eigen::MatrixXd::Random(appState->T.rows(), 9) * 1e-1;   //* 1e-1;
        appState->boundary_frames = Eigen::MatrixXd::Random(appState->cur_tet_mesh->nBoundaryElements(), 9);
    } else {
        loadPrimaryData();
    }
}

void Mint3DHook::subdivideMeshBoundary() {
    Eigen::MatrixXd V = appState->V;
    Eigen::MatrixXi T = appState->T;
    Eigen::MatrixXi F = appState->F;

    int nverts = V.rows();
    int ntets = T.rows();
    int nboundfaces = appState->cur_tet_mesh->nBoundaryElements();

    std::set<int> bound_verts;

    // Collect all boundary vertices
    for (int i = 0; i < nboundfaces; i++) {
        int fid = appState->cur_tet_mesh->boundaryFace(i);
        for (int j = 0; j < 3; j++) {
            bound_verts.insert(appState->cur_tet_mesh->faceVertex(fid, j));
        }
    }

    std::vector<int> to_subdiv;   // list of tets to subdivide

    // Identify tets with exactly 4 boundary vertices
    for (int i = 0; i < ntets; i++) {
        int nboundverts = 0;
        for (int j = 0; j < 4; j++) {
            if (bound_verts.find(T(i, j)) != bound_verts.end()) {
                nboundverts++;
            }
        }

        if (nboundverts == 4) {
            to_subdiv.push_back(i);
        }
    }

    // Calculate the size of the new V and T matrices
    int newVerts = nverts + to_subdiv.size();
    int newTets = ntets + to_subdiv.size() * 3;

    Eigen::MatrixXd V_subd(newVerts, 3);
    Eigen::MatrixXi T_subd(newTets, 4);

    V_subd.topRows(nverts) = V;

    int curtetidx = 0;
    int curbvertidx = nverts;

    for (int i = 0; i < ntets; i++) {
        if (!to_subdiv.empty() && to_subdiv.front() == i) {
            to_subdiv.erase(to_subdiv.begin());

            // Add a new vertex in the center of each boundary tet
            Eigen::Vector3d newVert = Eigen::Vector3d::Zero();
            for (int j = 0; j < 4; j++) {
                newVert += V.row(T(i, j));
            }
            newVert /= 4.0;
            V_subd.row(curbvertidx) = newVert;

            // Subdivide each boundary tet into 4 tets
            for (int j = 0; j < 4; j++) {
                Eigen::Vector4i newTet;
                newTet << T(i, j), T(i, (j + 1) % 4), T(i, (j + 2) % 4), curbvertidx;

                // Ensure the determinant is positive
                Eigen::Vector3d v0 = V_subd.row(newTet(0));
                Eigen::Vector3d v1 = V_subd.row(newTet(1));
                Eigen::Vector3d v2 = V_subd.row(newTet(2));
                Eigen::Vector3d v3 = V_subd.row(newTet(3));

                Eigen::Vector3d e1 = v1 - v0;
                Eigen::Vector3d e2 = v2 - v0;
                Eigen::Vector3d e3 = v3 - v0;

                Eigen::Vector3d n = e1.cross(e2);
                double det = n.dot(e3);

                if (det < 0) {
                    std::swap(newTet(2), newTet(3));
                }

                T_subd.row(curtetidx) = newTet;
                curtetidx++;
            }

            curbvertidx++;
        } else {
            T_subd.row(curtetidx) = T.row(i);
            curtetidx++;
        }
    }

    V_subd.conservativeResize(curbvertidx, 3);
    T_subd.conservativeResize(curtetidx, 4);

    std::cout << "Subdivided the boundary of the input mesh, added " << curtetidx - ntets << " tets and "
              << curbvertidx - nverts << " verts to the mesh" << std::endl;

    appState->V = V_subd;
    appState->T = T_subd;

    appState->cur_tet_mesh = std::make_unique<CubeCover::TetMeshConnectivity>(T_subd);
    if (!appState->cur_tet_mesh->isFaceConnected()) {
        std::cout << "Subdivided mesh is not face connected" << std::endl;
        appState->cur_tet_mesh->isManifold(true);
    }
}

// NOTE: This function *MUST* be called after subdivideMeshBoundary in order to guarentee that
//       every tet has at most one such facet.  If this condition is not met this code *will* be bugged!!!
void Mint3DHook::subdivideLockedFacets() {
    Eigen::MatrixXd V = appState->V;
    Eigen::MatrixXi T = appState->T;
    Eigen::MatrixXi F = appState->F;

    int nverts = V.rows();
    int ntets = T.rows();
    int nfaces = F.rows();

    int nboundfaces = appState->cur_tet_mesh->nBoundaryElements();

    std::set<int> bound_verts;

    // Collect all boundary vertices
    for (int i = 0; i < nboundfaces; i++) {
        int fid = appState->cur_tet_mesh->boundaryFace(i);
        for (int j = 0; j < 3; j++) {
            bound_verts.insert(appState->cur_tet_mesh->faceVertex(fid, j));
        }
    }

    std::vector<int> to_subdiv;   // list of facets to subdivide

    // Identify facets with exactly 3 boundary vertices
    for (int i = 0; i < nfaces; i++) {
        if (appState->cur_tet_mesh->isBoundaryFace(i)) continue;

        int nboundverts = 0;
        for (int j = 0; j < 3; j++) {
            int cur_vtx = appState->cur_tet_mesh->faceVertex(i, j);
            if (bound_verts.find(cur_vtx) != bound_verts.end()) {
                nboundverts++;
            }
        }

        if (nboundverts == 3) {
            to_subdiv.push_back(i);
        }
    }

    // Calculate the size of the new V and T matrices
    int newVerts = nverts + to_subdiv.size();
    int newTets = ntets + to_subdiv.size() * 4;

    Eigen::MatrixXd V_subd(newVerts, 3);
    Eigen::MatrixXi T_subd(newTets, 4);

    V_subd.topRows(nverts) = V;

    int curtetidx = 0;
    int curbvertidx = nverts;

    Eigen::VectorXi tetsToCopy = Eigen::VectorXi::Zero(ntets);

    for (int i = 0; i < nfaces; i++) {
        if (!to_subdiv.empty() && to_subdiv.front() == i) {
            to_subdiv.erase(to_subdiv.begin());

            // Add a new vertex in the center of each boundary tet
            Eigen::Vector3d newVert = Eigen::Vector3d::Zero();
            for (int j = 0; j < 3; j++) {
                int vtx = appState->cur_tet_mesh->faceVertex(i, j);
                newVert += V.row(vtx);
            }
            newVert /= 3.0;
            V_subd.row(curbvertidx) = newVert;

            // Subdivide both tets adjacent to this face into 6 tets
            for (int ft = 0; ft < 2; ft++) {
                int cft = appState->cur_tet_mesh->faceTet(i, ft);
                tetsToCopy(cft) = 1;

                int vid0 = appState->cur_tet_mesh->faceTetVertexIndex(i, ft, 0);
                int vid1 = appState->cur_tet_mesh->faceTetVertexIndex(i, ft, 1);
                int vid2 = appState->cur_tet_mesh->faceTetVertexIndex(i, ft, 2);
                int foppidx = -1;
                for (int iter = 0; iter < 4; iter++) {
                    if (iter != vid0 && iter != vid1 && iter != vid2) {
                        foppidx = iter;
                        break;
                    }
                }

                for (int j = 0; j < 3; j++) {
                    int vid0 = appState->cur_tet_mesh->faceTetVertexIndex(i, ft, j);
                    int vid1 = appState->cur_tet_mesh->faceTetVertexIndex(i, ft, (j + 1) % 3);

                    Eigen::Vector4i newTet;
                    newTet << T(cft, vid0), T(cft, vid1), T(cft, foppidx), curbvertidx;

                    // Ensure the determinant is positive
                    Eigen::Vector3d v0 = V_subd.row(newTet(0));
                    Eigen::Vector3d v1 = V_subd.row(newTet(1));
                    Eigen::Vector3d v2 = V_subd.row(newTet(2));
                    Eigen::Vector3d v3 = V_subd.row(newTet(3));

                    Eigen::Vector3d e1 = v1 - v0;
                    Eigen::Vector3d e2 = v2 - v0;
                    Eigen::Vector3d e3 = v3 - v0;

                    Eigen::Vector3d n = e1.cross(e2);
                    double det = n.dot(e3);

                    if (det < 0) {
                        std::swap(newTet(2), newTet(3));
                    }

                    T_subd.row(curtetidx) = newTet;
                    curtetidx++;
                }
            }
            curbvertidx++;
        }
    }

    for (int i = 0; i < ntets; i++) {
        if (tetsToCopy(i) == 0) {
            T_subd.row(curtetidx) = T.row(i);
            curtetidx++;
        }
    }

    V_subd.conservativeResize(curbvertidx, 3);
    T_subd.conservativeResize(curtetidx, 4);

    std::cout << "Subdivided the locked facets in the mesh, added " << curtetidx - ntets << " tets and "
              << curbvertidx - nverts << " verts to the mesh" << std::endl;
    // std::cout << "Subdivided the boundary of the input mesh" << std::endl;

    appState->V = V_subd;
    appState->T = T_subd;

    appState->cur_tet_mesh = std::make_unique<CubeCover::TetMeshConnectivity>(T_subd);
    if (!appState->cur_tet_mesh->isFaceConnected()) {
        std::cout << "Subdivided mesh is not face connected" << std::endl;
        appState->cur_tet_mesh->isManifold(true);
    }
}

void Mint3DHook::findSharpFeatures() {
    double nBoundElem = appState->cur_surf->nFaces();

    Eigen::MatrixXd V = appState->cur_surf->data().V;
    Eigen::MatrixXi F = appState->cur_surf->data().F;

    // std::cout << F << std::endl;

    appState->sharp_feature_edges.clear();

    double thresh = appState->sharp_feature_threshold * M_PI / 180.0;

    std::vector<Eigen::Vector3d> nodes;
    std::vector<Eigen::Vector2i> edges;

    int cedgeid = 0;

    for (int i = 0; i < nBoundElem; i++) {
        Eigen::Vector3d v0 = V.row(F(i, 0));
        Eigen::Vector3d v1 = V.row(F(i, 1));
        Eigen::Vector3d v2 = V.row(F(i, 2));

        Eigen::Vector3d e1 = v1 - v0;
        Eigen::Vector3d e2 = v2 - v0;

        Eigen::Vector3d nf = e1.cross(e2).normalized();

        std::vector<Eigen::Vector3d> face_sharp_edges;

        // face_sharp_edges.push_back(nf);
        // face_sharp_edges.push_back(nf);

        for (int j = 0; j < 3; j++) {
            double cedge = appState->cur_surf->data().faceEdges(i, j);
            int g = appState->cur_surf->data().faceNeighbors(i, j);
            if (g < 0) continue;

            Eigen::Vector3d u0 = V.row(F(g, 0));
            Eigen::Vector3d u1 = V.row(F(g, 1));
            Eigen::Vector3d u2 = V.row(F(g, 2));

            e1 = u1 - u0;
            e2 = u2 - u0;

            Eigen::Vector3d ng = e1.cross(e2).normalized();

            double dot_prod = std::clamp(nf.dot(ng), -.999999, .999999);
            double angle = std::acos(dot_prod);
            // double angle = std::acos(nf.dot(ng) * .999999);

            // std::cout << "cedge" << cedge << " g " << g << " angle " << angle << " nf.dot(ng) " << nf.dot(ng)
            //   << std::endl;
            if (angle > thresh) {
                int faceEdge = appState->cur_surf->data().faceEdges(i, j);

                // std::cout << "faceEdge " << faceEdge << std::endl;

                Eigen::Vector2i e = appState->cur_surf->data().edgeVerts.row(faceEdge);
                Eigen::Vector3d p0 = V.row(e(0));
                Eigen::Vector3d p1 = V.row(e(1));
                Eigen::Vector3d sharp_edge = (p0 - p1).normalized();
                // removed this as a hack to do tangent plane alignment
                face_sharp_edges.push_back(sharp_edge);
                Eigen::Vector2i ve = Eigen::Vector2i(cedgeid, cedgeid + 1);
                edges.push_back(ve);
                cedgeid += 2;

                // edges.push_back(appState->cur_surf->data().edgeVerts.row(faceEdge));
                nodes.push_back(p0);
                nodes.push_back(p1);
                // std::cout << "Face " << i << " is sharp" << std::endl;
            }
        }

        appState->sharp_feature_edges.push_back(face_sharp_edges);
    }

    appState->sharp_nodes = nodes;
    appState->sharp_edges = edges;
}

void Mint3DHook::launchFieldViewerFromCurrentState() {
    // std::string curr_frames = fileParser->getFileWithID("frames", "f_", ".bfra", appState->currentFileID);

    std::string curr_frames = "f_" + appState->meshName + "_" + std::to_string(appState->currentFileID) + ".bfra";
    Serialization::serializeMatrix(appState->frames, curr_frames);
    curr_frames = "frames_curr.bfra";
    Serialization::serializeMatrix(appState->frames, curr_frames);
    std::string mesh_path = appState->meshInputFilePath;

    std::string out_path;

    if (appState->directoryPath.empty() && appState->logFolderPath.empty()) {
        initializeLogFolder();
    }

    if (appState->logFolderPath.empty()) {
        out_path = appState->directoryPath + "/" + appState->meshName + "_f_" + std::to_string(appState->currentFileID);
    } else {
        out_path = appState->logFolderPath + "/" + appState->meshName + "_f_" + std::to_string(appState->currentFileID);
    }

    std::string command = "./bin/02_SymFrameView  -m " + mesh_path + " -f " + curr_frames + " -s " + out_path +
                          " --subdivision --regenerate-type=2 --renormalize-scale";
    std::cout << command << std::endl;
    std::thread t([command]() {
#ifdef __linux__
        // Set stack size to 2048MB using ulimit
        std::system("ulimit -s 8096000");   // 2048000");
#endif
        std::system(command.c_str());
    });
    t.detach();
    // std::system(command.c_str());
}

bool Mint3DHook::loadPrimaryData() {
    bool success = true;

    std::vector<FileTypeInfo> fileInfo = {
        {"frames", "f_", ".bfra", appState->frames}, {"boundary_frames", "bf_", ".bfra", appState->boundary_frames},
        // {"deltas_", ".bmom", appState->deltas},
        // {"moments_", ".bmom", appState->moments}
    };

    for (const auto &info : fileInfo) {
        // std::cout << "info.prefix " << info.prefix << std::endl;
        std::string file = fileParser->getFileWithID(info.folder, info.prefix, info.extension, appState->currentFileID);
        std::cout << file << std::endl;
        if (!file.empty()) {
            if (!Serialization::deserializeMatrix(info.targetMatrix, file)) {
                std::cerr << "Failed to load data for " << info.prefix << " from file: " << file << std::endl;
                success = false;
            }
            std::cout << "loaded " << info.prefix << " from file: " << file << std::endl;
            // std::cout << info.targetMatrix << std::endl;
        } else {
            std::cerr << "File not found for " << info.prefix << " with ID: " << appState->currentFileID << std::endl;
            success = false;
        }
    }

    return success;
}

bool Mint3DHook::loadSecondaryData() {
    bool success = true;

    std::string dataLogPath = appState->directoryPath + "/run_data";

    // Instead of Eigen matrices, directly use std::vector<double> for loading data
    std::vector<double> total_time, assembly_time, solve_time, line_search_time, identity_weight;
    std::vector<double> total_energy, energy_diff, solve_residual, gradient_norm, gradient_norm_step_start;
    std::vector<double> global_scale_log;
    std::vector<double> smoothness, unit_penalty, symmetric_integrability, primal_integrability, scaled_jacobian;
    std::vector<double> unit_barrier, curl_weight, asap_combed_smoothness, aiap_combed_smoothness;
    std::vector<double> combed_smoothness, combed_integrability, total_energy_unscale, global_scale_outiter;
    std::vector<double> fit_penalty;

    // Updated FileTypeInfo to work with std::vector<double>
    struct FileTypeInfo {
        std::string key;
        std::string subfolder;
        std::string filename;
        std::vector<double> &targetVector;   // Reference to the std::vector where data will be loaded

        FileTypeInfo(std::string key, std::string subfolder, std::string filename, std::vector<double> &targetVector)
            : key(key), subfolder(subfolder), filename(filename), targetVector(targetVector) {}
    };

    std::vector<FileTypeInfo> csvFileInfo = {
        {"total_time", "run_data", "total_time.csv", total_time},
        {"assembly_time", "run_data", "assembly_time.csv", assembly_time},
        {"solve_time", "run_data", "solve_time.csv", solve_time},
        {"line_search_time", "run_data", "line_search_time.csv", line_search_time},
        {"identity_weight", "run_data", "identity_weight.csv", identity_weight},
        {"total_energy", "run_data", "total_energy.csv", total_energy},
        {"energy_diff", "run_data", "energy_diff.csv", energy_diff},
        {"solve_residual", "run_data", "solve_residual.csv", solve_residual},
        {"gradient_norm", "run_data", "gradient_norm.csv", gradient_norm},
        {"gradient_norm_step_start", "run_data", "gradient_norm_step_start.csv", gradient_norm_step_start},
        {"global_scale_log", "run_data", "global_scale_log.csv", global_scale_log},
        {"smoothness", "run_data", "smoothness.csv", smoothness},
        {"unit_penalty", "run_data", "unit_penalty.csv", unit_penalty},
        {"symmetric_integrability", "run_data", "symmetric_integrability.csv", symmetric_integrability},
        {"primal_integrability", "run_data", "primal_integrability.csv", primal_integrability},
        {"scaled_jacobian", "run_data", "scaled_jacobian.csv", scaled_jacobian},
        {"unit_barrier", "run_data", "unit_barrier.csv", unit_barrier},
        {"curl_weight", "run_data", "curl_weight.csv", curl_weight},
        {"asap_combed_smoothness", "run_data", "asap_combed_smoothness.csv", asap_combed_smoothness},
        {"aiap_combed_smoothness", "run_data", "aiap_combed_smoothness.csv", aiap_combed_smoothness},
        {"combed_smoothness", "run_data", "combed_smoothness.csv", combed_smoothness},
        {"combed_integrability", "run_data", "combed_integrability.csv", combed_integrability},
        {"total_energy_unscale", "run_data", "total_energy_unscale.csv", total_energy_unscale},
        {"global_scale_outer_iter", "run_data", "global_scale_outer_iter.csv", global_scale_outiter},
        {"fit_penalty", "run_data", "fit_penalty.csv", fit_penalty}};

    for (const auto &info : csvFileInfo) {
        std::string fullFilePath = appState->directoryPath + "/" + info.subfolder + "/" + info.filename;
        if (std::filesystem::exists(fullFilePath)) {
            // Use Serialization::readCSV to directly load data into the std::vector
            if (!Serialization::readCSV(info.targetVector, fullFilePath)) {
                std::cerr << "Failed to load data for " << info.key << " from file: " << fullFilePath << std::endl;
                success = false;
            } else {
                std::cout << "Loaded " << info.key << " from file: " << fullFilePath << std::endl;

                // Assign the loaded vector to the appropriate member in appState
                if (info.key == "total_time") {
                    appState->total_time = info.targetVector;
                } else if (info.key == "assembly_time") {
                    appState->assembly_time = info.targetVector;
                } else if (info.key == "solve_time") {
                    appState->solve_time = info.targetVector;
                } else if (info.key == "line_search_time") {
                    appState->line_search_time = info.targetVector;
                } else if (info.key == "identity_weight") {
                    appState->identity_weight_log = info.targetVector;
                } else if (info.key == "total_energy") {
                    appState->total_energy_log = info.targetVector;
                } else if (info.key == "energy_diff") {
                    appState->energy_diff_log = info.targetVector;
                } else if (info.key == "solve_residual") {
                    appState->solve_residual_log = info.targetVector;
                } else if (info.key == "gradient_norm") {
                    appState->gradient_norm_log = info.targetVector;
                } else if (info.key == "gradient_norm_step_start") {
                    appState->gradient_norm_step_start_log = info.targetVector;
                } else if (info.key == "global_scale_log") {
                    appState->global_scale_log = info.targetVector;
                } else if (info.key == "smoothness") {
                    appState->smoothness_log = info.targetVector;
                } else if (info.key == "unit_penalty") {
                    appState->unit_penalty_log = info.targetVector;
                } else if (info.key == "symmetric_integrability") {
                    appState->symmetric_integrability_log = info.targetVector;
                } else if (info.key == "primal_integrability") {
                    appState->primal_integrability_log = info.targetVector;
                } else if (info.key == "scaled_jacobian") {
                    appState->scaled_jacobian_log = info.targetVector;
                } else if (info.key == "unit_barrier") {
                    // appState->unit_barrier_log = info.targetVector; // Uncomment and add to AppState if needed
                } else if (info.key == "curl_weight") {
                    appState->curl_weight_log = info.targetVector;
                } else if (info.key == "asap_combed_smoothness") {
                    appState->asap_combed_smoothness_log = info.targetVector;
                } else if (info.key == "aiap_combed_smoothness") {
                    appState->aiap_combed_smoothness_log = info.targetVector;
                } else if (info.key == "combed_smoothness") {
                    appState->combed_smoothness_log = info.targetVector;
                } else if (info.key == "combed_integrability") {
                    appState->combed_integrability_log = info.targetVector;
                } else if (info.key == "total_energy_unscale") {
                    appState->total_energy_log = info.targetVector;
                } else if (info.key == "global_scale_outer_iter") {
                    appState->global_scale_outiter_log = info.targetVector;
                } else if (info.key == "fit_penalty") {
                    appState->fit_penalty_log = info.targetVector;
                }
            }
        } else {
            std::cerr << "File not found: " << fullFilePath << std::endl;

            if (info.key == "global_scale_outer_iter") {
                appState->global_scale_outiter_log =
                    std::vector<double>(appState->primal_integrability_log.size(), 1.0);
            }
            success = false;
        }
    }

    return success;
}

/*

bool Mint3DHook::loadGuiState() {
    bool success = true;
    for (int i = 0; i < (int) Views::Field_View::gui_free; ++i) {
        Field_View view = static_cast<Field_View>(i);
        std::string fieldStub = Views::fieldViewToFileStub(view) + "_";
        std::string fieldFile = fileParser->getFileWithID(fieldStub, ".bdat", appState->currentFileID);

        if (!fieldFile.empty()) {
            Eigen::VectorXd data;
            if (!Serialization::deserializeVector(data, fieldFile)) {
                std::cerr << "Failed to load data for " << fieldStub << " from file: " << fieldFile << std::endl;
                success = false;
            }
            else {
                // Update the correct field in OutputState based on 'view'
                // Example: appState.os->norms_vec = data;
            }
        }
        else {
            std::cerr << "File not found for " << fieldStub << " with ID: " << appState->currentFileID << std::endl;
            success = false;
        }
    }
    return success;
}

*/

void Mint3DHook::initializeFieldState() {
    if (appState->solve_params.boundary_condition == MiNT3D::BoundaryCondition::Poisson) {
        setPoissonBoundaryConditions();
    } else if (appState->solve_params.boundary_condition == MiNT3D::BoundaryCondition::SDF) {
        setSDFBoundaryConditions();
    }

    if (appState->solve_params.init_state == MiNT3D::InitState::Random) {
        setFramesToRandom(appState->solve_params.w_init_scale);
    }
}

bool Mint3DHook::simulateOneStep() {
    // int iter = 0;

    initializeLogFolder();

    if (appState->restartFromCurrentFrame == false) {
        initializeFieldState();
    }

    // do up to 200 steps of an odeco solve to initialize mint_mesh
    if (appState->solve_params.b_init_as_odeco) {
        MiNT3D::SolveParams inp_params = appState->solve_params;
        appState->solve_params.setOdecoSolve();

        // these should be the only active states.
        appState->solve_params.w_unit = inp_params.w_unit;   // .1;
        appState->solve_params.w_orthog = inp_params.w_orthog;
        appState->solve_params.w_smooth_sym = inp_params.w_smooth_sym;

        solveFromCurrentState();
        appState->solve_params = inp_params;
    }

    solveFromCurrentState();

    // while(true)
    // {
    //     std::cout << iter << std::endl;
    //     iter++;
    //     std::this_thread::sleep_for(std::chrono::milliseconds(1000));
    // }
    std::cout << "Finished Running Simulation" << std::endl;
    appState->solveStatus = "Finished Running Simulation";
    appState->keepSolving = false;
    // should do some stuff here
    return true;
}

void Mint3DHook::solveFromCurrentState() {
    Eigen::MatrixXd V = appState->V;
    Eigen::MatrixXi T = appState->T;

    CubeCover::TetMeshConnectivity mesh(T);

    MiNT3D::MiNTModel model;
    model.SetMesh(V, mesh);
    model.SetSharpFeatureEdges(appState->sharp_feature_edges);
    model.logPath = appState->logFolderPath;
    appState->solve_params.b_headless_mode = appState->headless_mode;

    // model.exp_name = appState->experiment_name;

    model.solveData->total_time = appState->total_time;
    model.solveData->assembly_time = appState->assembly_time;
    model.solveData->solve_time = appState->solve_time;
    model.solveData->line_search_time = appState->line_search_time;
    model.solveData->identity_weight = appState->identity_weight_log;
    // model.solveData->total_energy = appState->total_energy_log;

    model.solveData->energy_diff = appState->energy_diff_log;
    model.solveData->solve_residual = appState->solve_residual_log;
    model.solveData->global_scale_log = appState->global_scale_log;
    model.solveData->gradient_norm = appState->gradient_norm_log;
    model.solveData->gradient_norm_step_start = appState->gradient_norm_step_start_log;

    model.energyLog->smoothness = appState->smoothness_log;
    model.energyLog->unit_penalty = appState->unit_penalty_log;
    model.energyLog->symmetric_integrability = appState->symmetric_integrability_log;
    model.energyLog->primal_integrability = appState->primal_integrability_log;
    model.energyLog->combed_smoothness = appState->combed_smoothness_log;
    model.energyLog->aiap_combed_smoothness = appState->aiap_combed_smoothness_log;
    model.energyLog->asap_combed_smoothness = appState->asap_combed_smoothness_log;
    model.energyLog->scaled_jacobian = appState->scaled_jacobian_log;
    model.energyLog->curl_weight = appState->curl_weight_log;
    model.energyLog->total_energy_unscale = appState->total_energy_log;
    model.energyLog->global_scale_outiter = appState->global_scale_outiter_log;
    model.energyLog->fit_penalty = appState->fit_penalty_log;

    // Eigen::MatrixXd frames = model.VectorToMatrix(appState->frames, 9);
    // Eigen::MatrixXd bnd_frames = model.VectorToMatrix(appState->boundary_frames, 9);
    Eigen::MatrixXd frames = appState->frames;
    Eigen::MatrixXd bnd_frames = appState->boundary_frames;

    Eigen::MatrixXd frames_out, bnd_frames_out;

    // this is for reproduce the crashed case
    //                appState->solve_params.w_outer_step = 100;

    MiNT3D::SolveParams solve_params = appState->solve_params;

    /////////////////////////////////////////////////
    ////  Callback which updates visualization state
    /////////////////////////////////////////////////
    std::function<void()> callback = [&]() {
        render_mutex.lock();

        // if (model.solveData->total_time.size() > appState->total_time.size() + 10) {
        appState->total_time = model.solveData->total_time;
        appState->assembly_time = model.solveData->assembly_time;
        appState->solve_time = model.solveData->solve_time;
        appState->line_search_time = model.solveData->line_search_time;
        appState->identity_weight_log = model.solveData->identity_weight;
        // appState->total_energy_log = model.solveData->total_energy;
        appState->energy_diff_log = model.solveData->energy_diff;
        appState->solve_residual_log = model.solveData->solve_residual;
        appState->global_scale_log = model.solveData->global_scale_log;

        appState->gradient_norm_log = model.solveData->gradient_norm;
        appState->gradient_norm_step_start_log = model.solveData->gradient_norm_step_start;

        appState->smoothness_log = model.energyLog->smoothness;
        appState->unit_penalty_log = model.energyLog->unit_penalty;
        appState->symmetric_integrability_log = model.energyLog->symmetric_integrability;
        appState->primal_integrability_log = model.energyLog->primal_integrability;
        appState->combed_smoothness_log = model.energyLog->combed_smoothness;
        appState->asap_combed_smoothness_log = model.energyLog->asap_combed_smoothness;
        appState->aiap_combed_smoothness_log = model.energyLog->aiap_combed_smoothness;
        appState->scaled_jacobian_log = model.energyLog->scaled_jacobian;
        appState->curl_weight_log = model.energyLog->curl_weight;
        appState->total_energy_log = model.energyLog->total_energy_unscale;
        appState->global_scale_outiter_log = model.energyLog->global_scale_outiter;

        appState->fit_penalty_log = model.energyLog->fit_penalty;

        // bool reload_frames = false;
        // if (appState->currentFileID == appState->fmax) {
        //     reload_frames = true;
        // }
        bool reload_frames = true;

        fileParser = std::make_unique<FileParser>(appState->logFolderPath, appState->loadInner);

        appState->fmin = fileParser.get()->minID;
        appState->fmax = fileParser.get()->maxID;

        if (reload_frames) {
            appState->currentFileID = appState->fmax;
            loadPrimaryData();
        }

        updateRenderGeometry();
        // }

        render_mutex.unlock();
    };

    // *appState->solve_params
    std::cout << "start solve" << std::endl;
    model.SolveSmoothFrames(solve_params, frames, frames_out, &bnd_frames, &bnd_frames_out, &callback);
    std::cout << "end solve" << std::endl;

    // visualize optimization mesh

    if (!appState->headless_mode) {
        polyscope::registerTetMesh("opt", model.extended_V_, model.extended_T_);
        polyscope::getVolumeMesh("opt")->setEdgeWidth(0.6)->setTransparency(0.5);
        polyscope::getVolumeMesh("opt")->setEnabled(false);
    }
    CubeCover::exportMESH("cur_opt.mesh", model.extended_V_, model.extended_T_);

    // TODO: this is a hack, need to remove
    // Serialization::serializeMatrix(frames_out, "frames_out.bfra", 3);
    // Serialization::serializeMatrix(bnd_frames_out, "bnd_frames_out.bfra", 3);

    // Eigen::VectorXd f_flat(frames_out.size());
    // Eigen::VectorXd bnd_flat(bnd_frames_out.size());

    // for (int i = 0; i < frames_out.rows(); i++) {
    //     f_flat.segment<9>(i * 9) = frames_out.row(i);
    // }

    // for (int i = 0; i < bnd_frames_out.rows(); i++) {
    //     bnd_flat.segment<9>(i * 9) = bnd_frames_out.row(i);
    //     // std::cout << bnd_frames_out.row(i) << std::endl;
    // }

    // Eigen::VectorXd x(f_flat.size() + bnd_flat.size());
    // x << f_flat, bnd_flat;

    // reload files to browse

    // update viz state
    // std::string load_dir;
    // if (appState->loadInner) {
    //     load_dir = model.logInnerIter;
    // } else {
    //     load_dir = model.logOuterIter;
    // }

    // appState->boundary_frames = bnd_frames_out;
    // appState->frames = frames_out;

    appState->currentFileID = appState->fmax;
    callback();
    // updateRenderGeometry();
}

// bool loadMetricDrivenFrameField(std::string path);
bool Mint3DHook::loadMetricDrivenFrameField(std::string path) {
    // std::string path = appState->logFolderPath + "/metric_driven_frame_field.bfra";
    Eigen::MatrixXd frames;
    if (!Serialization::deserializeFF3FramesFromMetricDrivenFrames3D(frames, path)) {
        std::cerr << "Failed to load metric driven frame field from file: " << path << std::endl;
        return false;
    }

    appState->frames = frames;

    // set boundary frames from neighbors
    appState->boundary_frames = Eigen::MatrixXd::Zero(appState->cur_tet_mesh->nBoundaryElements(), 9);
    for (int i = 0; i < appState->cur_tet_mesh->nBoundaryElements(); i++) {
        int fid = appState->cur_tet_mesh->boundaryFace(i);

        int bound_neighbor_idx = appState->cur_tet_mesh->boundaryElementTet(i);
        appState->boundary_frames.row(i) = appState->frames.row(bound_neighbor_idx);
    }

    appState->solve_params.setSoftAlignSolve();
    appState->show_frames = true;
    appState->gui_vec_size = .01;
    // appState->solve_params.boundary_hard_constraints = MiNT3D::BoundaryHardConstraints::NoHardConstraints;

    updateRenderGeometry();
    return true;
}

// works for both .bfra and .fra files
bool Mint3DHook::loadSpecificFrameField(std::string fraFilename) {
    // std::string path = appState->logFolderPath + "/metric_driven_frame_field.bfra";
    Eigen::MatrixXd frames;
    Eigen::MatrixXi assignments;

    bool read_succeeded = CubeCover::readFrameField(fraFilename, "", appState->T, frames, assignments, true);
    if (!read_succeeded) {
        std::cerr << "Failed to load frame field: " << fraFilename << std::endl;
        return false;
    }

    // if (frames.rows() == 3 * appState->T.rows()) {
    //     appState->frames = Eigen::MatrixXd::Zero(appState->T.rows(), 9);
    //     for (int i = 0; i < appState->T.rows(); i++) {
    //         for (int j = 0; j < 3; j++) {
    //             appState->frames.row(i).segment<3>(3 * j) = frames.row(3 * i + j);
    //         }
    //     }
    // } else {
    //
    // }

    appState->frames = frames;

    std::cout << frames.row(0) << std::endl;
    std::cout << frames.row(1) << std::endl;
    std::cout << frames.row(2) << std::endl;

    CubeCover::TetMeshConnectivity mesh(appState->T);

    CubeCover::FrameField *ff = CubeCover::fromFramesAndAssignments(mesh, frames, assignments, true);
    CubeCover::writeFrameField("test.fra", "test.perm", *ff);

    // set boundary frames from neighbors
    appState->boundary_frames = Eigen::MatrixXd::Zero(appState->cur_tet_mesh->nBoundaryElements(), 9);
    for (int i = 0; i < appState->cur_tet_mesh->nBoundaryElements(); i++) {
        int fid = appState->cur_tet_mesh->boundaryFace(i);

        int bound_neighbor_idx = appState->cur_tet_mesh->boundaryElementTet(i);
        appState->boundary_frames.row(i) = appState->frames.row(bound_neighbor_idx);
    }

    appState->solve_params.setSoftAlignSolve();
    appState->show_frames = true;
    appState->gui_vec_size = .01;
    appState->solve_params.boundary_hard_constraints = MiNT3D::BoundaryHardConstraints::NoHardConstraints;

    updateRenderGeometry();
    return true;
}

// take frame field dual
void Mint3DHook::takeFrameFieldDual() {
    Eigen::MatrixXd frames = appState->frames;
    Eigen::MatrixXd boundary_frames = appState->boundary_frames;
    Eigen::MatrixXi T = appState->T;

    Eigen::MatrixXd dual_frames = Eigen::MatrixXd::Zero(T.rows(), 9);

    for (int i = 0; i < T.rows(); i++) {
        Eigen::VectorXd cur_frame = frames.row(i);
        Eigen::VectorXd dual_frame = Eigen::VectorXd::Zero(9);

        Eigen::Vector3d v0 = cur_frame.segment<3>(0);
        Eigen::Vector3d v1 = cur_frame.segment<3>(3);
        Eigen::Vector3d v2 = cur_frame.segment<3>(6);

        Eigen::Vector3d d0 = v1.cross(v0).normalized();
        Eigen::Vector3d d1 = v2.cross(v0).normalized();
        Eigen::Vector3d d2 = v2.cross(v1).normalized();

        dual_frame.segment<3>(0) = d0;
        dual_frame.segment<3>(3) = d1;
        dual_frame.segment<3>(6) = d2;

        dual_frames.row(i) = dual_frame;
    }

    appState->frames = dual_frames;

    if (appState->useBoundaryFrames) {
        Eigen::MatrixXd dual_boundary_frames = Eigen::MatrixXd::Zero(appState->cur_tet_mesh->nBoundaryElements(), 9);
        for (int i = 0; i < appState->cur_tet_mesh->nBoundaryElements(); i++) {
            Eigen::VectorXd cur_frame = boundary_frames.row(i);
            Eigen::VectorXd dual_frame = Eigen::VectorXd::Zero(9);

            Eigen::Vector3d v0 = cur_frame.segment<3>(0);
            Eigen::Vector3d v1 = cur_frame.segment<3>(3);
            Eigen::Vector3d v2 = cur_frame.segment<3>(6);

            Eigen::Vector3d d0 = v1.cross(v0).normalized();
            Eigen::Vector3d d1 = v2.cross(v0).normalized();
            Eigen::Vector3d d2 = v2.cross(v1).normalized();

            dual_frame.segment<3>(0) = d0;
            dual_frame.segment<3>(3) = d1;
            dual_frame.segment<3>(6) = d2;

            dual_boundary_frames.row(i) = dual_frame;
        }

        appState->boundary_frames = dual_boundary_frames;
    }

    updateRenderGeometry();
}

/// @brief This should only be called after the appstate has been initialized with a mesh
void Mint3DHook::resetAppState() {
    // Resetting simulation parameters to default or initial values
    // appState->currentIteration = 0;
    // appState->innerLoopIteration = 0;
    // appState->max_saved_index = 0;

    // appState->currentFileID = 0;
    // appState->maxIterations = 999999; // Default maximum iterations
    // appState->convergenceEpsilon = 1e-12;// 1e-9;
    // appState->outerLoopIteration = 0;

    appState->override_bounds.lower = 0;
    appState->override_bounds.upper = 1e-5;
    // appState->shouldReload = false;

    // appState->opt_step_start_time = std::chrono::high_resolution_clock::now();

    // Resetting mesh data
    appState->frames.setZero(appState->T.rows(), 9);                                     // TODO make this generic
    appState->boundary_frames.setZero(appState->cur_tet_mesh->nBoundaryElements(), 9);   // TODO make this generic

    appState->prev_frame_element = Field_View::Element_COUNT;

    appState->tet_centroids = Eigen::MatrixXd::Zero(appState->T.rows(), 3);
    for (int i = 0; i < appState->T.rows(); i++) {
        appState->tet_centroids.row(i) = (appState->V.row(appState->T(i, 0)) + appState->V.row(appState->T(i, 1)) +
                                          appState->V.row(appState->T(i, 2)) + appState->V.row(appState->T(i, 3))) /
                                         4.0;
    }

    // // Reinitialize the boundary conditions if needed
    // initBoundaryConditions();

    // // Initialize symmetric curl operators
    // initCurlOperators();

    if (!appState->headless_mode) {
        // compute stuff for drawing boundary
        appState->bound_centroids = Eigen::MatrixXd::Zero(appState->cur_tet_mesh->nBoundaryElements(), 3);
        appState->bound_normals = appState->bound_centroids;
        appState->bound_b1 = appState->bound_centroids;
        appState->bound_b2 = appState->bound_centroids;
        for (int i = 0; i < appState->bound_centroids.rows(); i++) {
            int boundaryFace = appState->cur_tet_mesh->boundaryFace(i);
            Eigen::Vector3d a = appState->V.row(appState->cur_tet_mesh->faceVertex(boundaryFace, 0));
            Eigen::Vector3d b = appState->V.row(appState->cur_tet_mesh->faceVertex(boundaryFace, 1));
            Eigen::Vector3d c = appState->V.row(appState->cur_tet_mesh->faceVertex(boundaryFace, 2));

            Eigen::Vector3d b1 = (b - a).normalized();
            Eigen::Vector3d b2 = (c - a).normalized();
            Eigen::Vector3d n = b1.cross(b2).normalized();

            appState->bound_centroids.row(i) = (a + b + c) / 3.0;
            appState->bound_normals.row(i) = n;
            appState->bound_b1.row(i) = b1;
            appState->bound_b2.row(i) = n.cross(b1);
        }

        // reset polyscope state

        // Optionally, re-register mesh with Polyscope if visualization needs a reset
        polyscope::removeAllStructures();
        polyscope::registerTetMesh("c", appState->V, appState->T);
        polyscope::getVolumeMesh("c")->setEdgeWidth(0.6)->setTransparency(0.5);
        polyscope::getVolumeMesh("c")->setEnabled(false);
        // polyscope::getVolumeMesh("c")->setEdgeWidth(0.6);
        polyscope::view::resetCameraToHomeView();

        polyscope::registerPointCloud("c_vecs", appState->tet_centroids)->setPointRadius(0.0);

        int num_vecs = 3;   // appState->frames.size();
        for (int v = 0; v < num_vecs; v++) {
            // Eigen::MatrixXd cur_vec = outputData->frames[v];
            Eigen::MatrixXd cur_vec = Eigen::MatrixXd::Zero(appState->frames.rows(), 3);

            double color_shift = (v + 1.) * 1.0 / num_vecs;

            auto vectorField = polyscope::getPointCloud("c_vecs")->addVectorQuantity(
                "Vector Field " + std::to_string(v), cur_vec, polyscope::VectorType::AMBIENT);

            auto vectorFieldNeg = polyscope::getPointCloud("c_vecs")->addVectorQuantity(
                "Vector Field (negative) " + std::to_string(v), (-1.) * cur_vec, polyscope::VectorType::AMBIENT);
            if (v % 3 == 0) {
                vectorField->setVectorColor(glm::vec3(color_shift, 0.1, 0.1));
                vectorFieldNeg->setVectorColor(glm::vec3(color_shift, 0.4, 0.4));
            } else if (v % 3 == 1) {
                vectorField->setVectorColor(glm::vec3(0.1, color_shift, 0.1));
                vectorFieldNeg->setVectorColor(glm::vec3(0.4, color_shift, 0.4));
            } else if (v % 3 == 2) {
                vectorField->setVectorColor(glm::vec3(0.1, 0.1, color_shift));
                vectorFieldNeg->setVectorColor(glm::vec3(0.4, 0.4, color_shift));
            }

            vectorField->setVectorRadius(0.003);
            vectorFieldNeg->setVectorRadius(0.003);
        }

        polyscope::registerPointCloud("surface_vecs", appState->bound_centroids)->setPointRadius(0.0);

        // int num_vecs = 3; // appState->frames.size();
        for (int v = 0; v < num_vecs; v++) {
            // Eigen::MatrixXd cur_vec = outputData->frames[v];
            Eigen::MatrixXd cur_vec = Eigen::MatrixXd::Zero(appState->boundary_frames.rows(), 3);

            double color_shift = .5 + (v + 1.) * 1.0 / (2 * num_vecs);

            auto vectorField =
                polyscope::getPointCloud("surface_vecs")
                    ->addVectorQuantity("Vector Field " + std::to_string(v), cur_vec, polyscope::VectorType::AMBIENT);

            auto vectorFieldNeg = polyscope::getPointCloud("surface_vecs")
                                      ->addVectorQuantity("Vector Field (negative) " + std::to_string(v),
                                                          (-1.) * cur_vec, polyscope::VectorType::AMBIENT);

            if (v % 3 == 0) {
                vectorField->setVectorColor(glm::vec3(color_shift, 0.1, 0.2));
                vectorFieldNeg->setVectorColor(glm::vec3(color_shift, 0.9, 0.8));
            } else if (v % 3 == 1) {
                vectorField->setVectorColor(glm::vec3(0.1, color_shift, 0.2));
                vectorFieldNeg->setVectorColor(glm::vec3(0.9, color_shift, 0.8));
            } else if (v % 3 == 2) {
                vectorField->setVectorColor(glm::vec3(0.1, 0.2, color_shift));
                vectorFieldNeg->setVectorColor(glm::vec3(0.9, 0.8, color_shift));
            }

            vectorField->setVectorRadius(0.005);
            vectorFieldNeg->setVectorRadius(0.005);
        }
    }
}

void Mint3DHook::initCurlOperators() {}

std::string Mint3DHook::getExperimentShortString() {
    char buffer[256];   // Adjust size as needed
    std::string result;

    if (appState->solve_params.b_smooth_combed) {
        std::ostringstream oss;
        oss << "sc_" << std::scientific << std::setprecision(2) << appState->solve_params.w_smooth_combed << "_";
        result += oss.str();
    }
    if (appState->solve_params.b_smooth_sym) {
        std::ostringstream oss;
        oss << "ss_" << std::scientific << std::setprecision(2) << appState->solve_params.w_smooth_sym << "_";
        result += oss.str();
    }
    if (appState->solve_params.b_smooth_asap_combing) {
        result += "sasap_";
    }
    if (appState->solve_params.b_int_combed) {
        std::ostringstream oss;
        oss << "icfw_" << std::scientific << std::setprecision(2) << appState->solve_params.w_int_combed_fixed_weight
            << "_";
        result += oss.str();
    }
    if (appState->solve_params.b_int_sym) {
        std::ostringstream oss;
        oss << "is_" << std::scientific << std::setprecision(2) << appState->solve_params.w_int_sym << "_";
        result += oss.str();
    }
    if (appState->solve_params.b_orthog) {
        std::ostringstream oss;
        oss << "o_" << std::scientific << std::setprecision(2) << appState->solve_params.w_orthog << "_";
        result += oss.str();
    }
    if (appState->solve_params.b_unit) {
        std::ostringstream oss;
        oss << "u_" << std::scientific << std::setprecision(2) << appState->solve_params.w_unit << "_";
        result += oss.str();
    }
    if (appState->solve_params.b_fit && !appState->solve_params.b_self_align) {
        std::ostringstream oss;
        oss << "f_" << std::scientific << std::setprecision(2) << appState->solve_params.w_fit << "_";
        result += oss.str();
    }
    if (appState->solve_params.b_self_align) {
        std::ostringstream oss;
        oss << "sa_" << std::scientific << std::setprecision(2) << appState->solve_params.w_self_align << "_";
        result += oss.str();
    }
    if (appState->solve_params.b_viscosity) {
        std::ostringstream oss;
        oss << "vsc_" << std::scientific << std::setprecision(2) << appState->solve_params.w_viscocity << "_";
        result += oss.str();
    }
    if (appState->solve_params.b_feat_align) {
        std::ostringstream oss;
        oss << "plane_1e-6_" << std::scientific << std::setprecision(2) << appState->solve_params.w_feat_align << "_";
        result += oss.str();
    }
    if (appState->solve_params.boundary_condition == MiNT3D::BoundaryCondition::Poisson) {
        std::ostringstream oss;
        oss << "bpois_";
        result += oss.str();
    }
    if (appState->solve_params.boundary_condition == MiNT3D::BoundaryCondition::SDF) {
        std::ostringstream oss;
        oss << "bsdf_";
        result += oss.str();
    }
    if (appState->solve_params.init_state == MiNT3D::InitState::Random) {
        std::ostringstream oss;
        oss << "initrand_" << std::scientific << std::setprecision(2) << appState->solve_params.w_init_scale << "_";
        result += oss.str();
    }

    // Remove the trailing space
    if (!result.empty()) {
        result.pop_back();
    }

    return result;

    // // Format the final string
    // sprintf(buffer, "Selected options: %s", result.c_str());
    // return std::string(buffer);
}

void Mint3DHook::initializeLogFolder() {
    // Init Log Folder

    auto now = std::chrono::system_clock::now();
    auto date = date::floor<date::days>(now);
    auto ymd = date::year_month_day{date};
    auto time = date::make_time(std::chrono::duration_cast<std::chrono::milliseconds>(now - date));

    // pass folder name as an argument
    // std::string exp_name = "debug";
    // std::string exp_name = "smoothness_krushkal";

    std::string exp_name = appState->exp_name;

    appState->meshName = (fs::path(appState->meshInputFilePath).stem()).string();
    // char hostname[1024];
    // gethostname(hostname, 1024);
    // std::cout << hostname << std::endl;

    // experimental setup
    // std::string exp_setup = "";
    std::string exp_setup = getExperimentShortString();

    // std::cout << exp_setup << std::endl;

    // Construct the log folder path using the mesh name and the current date-time
    std::ostringstream folderStream;
    folderStream << appState->out_dir << exp_name
                 << "/"
                 // << appState->meshName << "/"
                 //<< hostname << "_"
                 << static_cast<unsigned>(ymd.month()) << "_" << static_cast<unsigned>(ymd.day()) << "_"
                 << time.hours().count() << "_" << time.minutes().count() << "_" << appState->experiment_name << "_"
                 << appState->meshName << "_" << exp_setup;

    appState->logFolderPath = folderStream.str();
    std::filesystem::create_directories(appState->logFolderPath);

    std::string targ_mesh = appState->logFolderPath + "/" + appState->meshName + ".mesh";
    appState->meshInputFilePath = targ_mesh;

    CubeCover::exportMESH(targ_mesh, appState->V, appState->T);

    std::cout << "wrote current mesh to: " << targ_mesh << std::endl;
}

void Mint3DHook::setFramesToInvertingOnCylinder(bool isInverting, double multiplier, double noise_magnitude) {
    // // Assuming boundary faces are identified in AppState
    // Eigen::MatrixXi K;

    // appState->frames.resize(appState->T.rows(), DOFS_PER_ELEMENT);
    // appState->boundary_frames.resize(appState->cur_tet_mesh->nBoundaryElements(), DOFS_PER_ELEMENT);

    // int ntets = appState->frames.rows();
    // int nbound = appState->boundary_frames.rows();
    // int nvars = DOFS_PER_ELEMENT;
    // appState->pinned_frame_idx.resize(ntets + nbound);
    // appState->pinned_frame_idx.setZero();

    Eigen::Vector3d min = Eigen::Vector3d::Ones() * 1e10;
    Eigen::Vector3d max = Eigen::Vector3d::Ones() * -1e10;

    this->resetAppState();

    for (int i = 0; i < appState->bound_centroids.rows(); i++) {
        Eigen::VectorXd centroid = appState->bound_centroids.row(i);
        for (int j = 0; j < 3; j++) {
            if (centroid(j) < min(j)) min(j) = centroid(2);
            if (centroid(j) > max(j)) max(j) = centroid(2);
        }
    }

    double range = max(2) - min(2);
    Eigen::Vector3d center = (max - min) / 2.0;
    center(2) = 0;

    // std::cout << center << std::endl;

    int nbound = 0;

    appState->pinned_orthog_b1 = Eigen::MatrixXd::Zero(nbound, 9);
    appState->pinned_orthog_b2 = Eigen::MatrixXd::Zero(nbound, 9);

    // Initialize boundary conditions
    for (int i = 0; i < appState->bound_centroids.rows(); i++) {
        Eigen::VectorXd centroid = appState->bound_centroids.row(i);

        // Set frame orientation based on the centroid
        Eigen::Vector2d vec = Eigen::Vector2d(centroid.x() - center.x(), centroid.y() - center.y());
        double z = ((centroid(2) - min(2)) / range) * 3.1415 / 2.;   //.normalized();
        double r = vec.norm();

        vec.normalize();

        Eigen::VectorXd frame = getInvertingFrame(vec, z, isInverting, multiplier, noise_magnitude);

        appState->boundary_frames.row(i) = frame;
    }

    for (int i = 0; i < appState->tet_centroids.rows(); i++) {
        Eigen::VectorXd centroid = appState->tet_centroids.row(i);

        // Set frame orientation based on the centroid
        Eigen::Vector2d vec = Eigen::Vector2d(centroid.x() - center.x(), centroid.y() - center.y());
        double z = ((centroid(2) - min(2)) / range) * 3.1415 / 2;   //.normalized();
        double r = vec.norm();

        vec.normalize();

        Eigen::VectorXd frame = getInvertingFrame(vec, z, isInverting, multiplier, noise_magnitude);

        appState->frames.row(i) = frame;
    }
}

Eigen::VectorXd Mint3DHook::getInvertingFrame(Eigen::Vector2d vec, double z, bool isInverting, double multiplier,
                                              double noise_magnitude) {
    Eigen::Matrix3d sing;
    Eigen::Matrix3d rot;
    Eigen::Matrix3d frame_cur;

    double noise = noise_magnitude * ((double)rand() / RAND_MAX - 0.5);

    double theta = atan2(vec(1), vec(0)) * multiplier + noise;   // acos(vec(1)) * .5;
    sing << cos(theta), sin(theta), 0., -sin(theta), cos(theta), 0., 0., 0., 1.;
    rot << cos(z), 0, sin(z), 0., 1., 0., -sin(z), 0, cos(z);

    if (isInverting) {
        frame_cur = rot * sing;
    } else {
        frame_cur = sing;
    }

    // Integrable frame field on a sphere.

    Eigen::VectorXd frame = Eigen::VectorXd::Zero(9);

    frame(0) = frame_cur(1, 0);
    frame(1) = frame_cur(1, 1);
    frame(2) = frame_cur(1, 2);

    frame(3) = frame_cur(0, 0);
    frame(4) = frame_cur(0, 1);
    frame(5) = frame_cur(0, 2);

    frame(6) = frame_cur(2, 0);
    frame(7) = frame_cur(2, 1);
    frame(8) = frame_cur(2, 2);

    return frame;
}

Eigen::VectorXd Mint3DHook::solvePoissonOnTetMesh(const Eigen::MatrixXd &V, const Eigen::MatrixXi &T) {
    // Number of vertices
    const int n = V.rows();

    // Compute Laplacian and mass matrices
    Eigen::SparseMatrix<double> L, M;
    igl::cotmatrix(V, T, L);
    igl::massmatrix(V, T, igl::MASSMATRIX_TYPE_BARYCENTRIC, M);

    // Find boundary faces
    Eigen::MatrixXi F;
    igl::boundary_facets(T, F);

    // Extract unique boundary vertices
    std::set<int> boundary_vertices_set;
    for (int i = 0; i < F.rows(); ++i) {
        for (int j = 0; j < F.cols(); ++j) {
            boundary_vertices_set.insert(F(i, j));
        }
    }
    Eigen::VectorXi b(boundary_vertices_set.size());
    int idx = 0;
    for (int v : boundary_vertices_set) {
        b(idx++) = v;
    }

    // Create a map to identify interior vertices
    std::vector<int> interior_indices;
    Eigen::VectorXi index_map(n);
    for (int i = 0; i < n; ++i) {
        if (boundary_vertices_set.find(i) == boundary_vertices_set.end()) {
            interior_indices.push_back(i);
        }
        index_map(i) = -1;
    }
    for (int i = 0; i < interior_indices.size(); ++i) {
        index_map(interior_indices[i]) = i;
    }

    // Create selection matrix S
    Eigen::SparseMatrix<double> S(interior_indices.size(), n);
    std::vector<Eigen::Triplet<double>> S_triplets;
    for (int i = 0; i < interior_indices.size(); ++i) {
        S_triplets.emplace_back(i, interior_indices[i], 1.0);
    }
    S.setFromTriplets(S_triplets.begin(), S_triplets.end());

    // Adjust the system to account for Dirichlet boundary conditions
    Eigen::SparseMatrix<double> A = S * L * S.transpose();
    Eigen::VectorXd B = Eigen::VectorXd::Constant(n, 1.0);
    Eigen::VectorXd B_mod = S * M * B;

    // Use Eigen's SPQR solver
    // Eigen::SPQR<Eigen::SparseMatrix<double>> solver;
    // solver.compute(A);

    Eigen::CholmodSimplicialLDLT<Eigen::SparseMatrix<double>> solver;
    solver.analyzePattern(A);
    solver.factorize(A);

    // std::cout << A << std::endl;

    if (solver.info() != Eigen::Success) {
        std::cerr << "Decomposition failed" << std::endl;
        exit(1);
    }
    // B_mod = S * (M * B_mod);
    Eigen::VectorXd v_interior = solver.solve(B_mod);
    if (solver.info() != Eigen::Success) {
        std::cerr << "Solving failed" << std::endl;
        exit(1);
    }

    // Expand solution back to the full set of vertices
    Eigen::VectorXd v = Eigen::VectorXd::Zero(n);
    for (int i = 0; i < interior_indices.size(); ++i) {
        v(interior_indices[i]) = v_interior(i);
    }

    return v;
}

void Mint3DHook::setPoissonBoundaryConditions() {
    Eigen::MatrixXd V = appState->V;
    Eigen::MatrixXi T = appState->T;
    Eigen::MatrixXi F = appState->F;

    // Solve the Poisson equation with unit Neumann boundary conditions
    // Eigen::VectorXd S = solvePoissonEquationWithNeumannBC(V, T, F);
    Eigen::VectorXd S = solvePoissonOnTetMesh(V, T);

    CubeCover::TetMeshConnectivity mesh(T);

    Eigen::MatrixXd gradients(T.rows(), 3);   // Store gradients for each tetrahedron
    Eigen::Vector3d c1, c2;
    c1 << .1, 0, 0;
    c2 << 0, .1, 0;

    for (int i = 0; i < T.rows(); ++i) {
        Eigen::Vector4i tet = T.row(i);
        Eigen::Vector3d e0 = V.row(tet(1)) - V.row(tet(0));
        Eigen::Vector3d e1 = V.row(tet(2)) - V.row(tet(0));
        Eigen::Vector3d e2 = V.row(tet(3)) - V.row(tet(0));

        Eigen::Matrix3d J;
        J.col(0) = e0;
        J.col(1) = e1;
        J.col(2) = e2;

        Eigen::Vector3d grad_param_space;
        grad_param_space << S(tet(1)) - S(tet(0)), S(tet(2)) - S(tet(0)), S(tet(3)) - S(tet(0));

        Eigen::Vector3d grad = J.inverse().transpose() * grad_param_space;

        // hook->appState->frames.row(i).segment<3>(3) =
        //     scale * hook->appState->frames.row(i).segment<3>(3).normalized();
        // hook->appState->frames.row(i).segment<3>(6) =
        //     scale * hook->appState->frames.row(i).segment<3>(6).normalized();

        appState->frames.row(i).segment<3>(0) = -grad;

        appState->frames.row(i).segment<3>(3) = c1;
        appState->frames.row(i).segment<3>(6) = c2;
    }

    double ave_bound_norm = 0;

    for (int i = 0; i < mesh.nBoundaryElements(); i++) {
        int btet = mesh.boundaryElementTet(i);
        Eigen::Vector3d n = appState->frames.row(btet).segment<3>(0);
        appState->boundary_frames.row(i).segment<3>(0) = -n;
        appState->boundary_frames.row(i).segment<3>(3) = c1;
        appState->boundary_frames.row(i).segment<3>(6) = c2;

        ave_bound_norm += n.norm();
    }

    ave_bound_norm = ave_bound_norm / mesh.nBoundaryElements();

    appState->frames *= 1.0 / ave_bound_norm;
    appState->boundary_frames *= 1.0 / ave_bound_norm;

    for (int i = 0; i < T.rows(); i++) {
        // hook->appState->frames.row(i).segment<3>(0) *= scale;
    }
    // std::cout << S << std::endl;
    if (!appState->headless_mode) {
        polyscope::getVolumeMesh("c")->addVertexScalarQuantity("poisson problem field", -S);
    }

    // updateRenderGeometry();
}

void Mint3DHook::setSDFBoundaryConditions() {
    Eigen::MatrixXd V = appState->V;
    Eigen::MatrixXi T = appState->T;
    Eigen::MatrixXi F = appState->F;

    CubeCover::TetMeshConnectivity mesh(T);

    Eigen::MatrixXd fields(T.rows(), 9);
    fields.setZero();

    // Remove unreferenced vertices
    Eigen::MatrixXd V_new;
    Eigen::MatrixXi F_new;
    Eigen::VectorXi I;
    igl::remove_unreferenced(V, F, V_new, F_new, I);

    Eigen::VectorXd S = Eigen::MatrixXd::Zero(V.rows(), 1);   // Signed distances

    Eigen::MatrixXd C;   // Closest points
    // Eigen::VectorXi I;   // Closest simplex indices
    Eigen::VectorXd N;   // Normal at the closest points

    // SIGNED_DISTANCE_TYPE_UNSIGNED
    // SIGNED_DISTANCE_TYPE_FAST_WINDING_NUMBER
    igl::signed_distance(V, V_new, F_new, igl::SIGNED_DISTANCE_TYPE_UNSIGNED, S, I, C, N);

    // Compute the gradient of the signed distance field within each tetrahedron
    Eigen::MatrixXd gradients(T.rows(), 3);   // Store gradients for each tetrahedron
    Eigen::Vector3d c1, c2;
    // c1 << .1, 0, 0;
    // c2 << 0, .1, 0;

    for (int i = 0; i < T.rows(); ++i) {
        Eigen::Vector4i tet = T.row(i);
        Eigen::Vector3d e0 = V.row(tet(1)) - V.row(tet(0));
        Eigen::Vector3d e1 = V.row(tet(2)) - V.row(tet(0));
        Eigen::Vector3d e2 = V.row(tet(3)) - V.row(tet(0));

        Eigen::Matrix3d J;
        J.col(0) = e0;
        J.col(1) = e1;
        J.col(2) = e2;

        Eigen::Vector3d grad_param_space;
        grad_param_space << S(tet(1)) - S(tet(0)), S(tet(2)) - S(tet(0)), S(tet(3)) - S(tet(0));

        Eigen::Vector3d grad = J.inverse().transpose() * grad_param_space;

        appState->frames.row(i).segment<3>(0) = grad;
        Eigen::Vector3d r = Eigen::Vector3d::Random();
        c1 = r.cross(grad).normalized();
        c2 = grad.cross(c1).normalized();

        appState->frames.row(i).segment<3>(3) = c1 * 1e-4;
        appState->frames.row(i).segment<3>(6) = c2 * 1e-4;
        // appState->frames.row(i).segment<3>(3) = c1 ;
        // appState->frames.row(i).segment<3>(6) = c2 ;
    }

    for (int i = 0; i < mesh.nBoundaryElements(); i++) {
        int btet = mesh.boundaryElementTet(i);
        Eigen::Vector3d n = -appState->frames.row(btet).segment<3>(0);
        appState->boundary_frames.row(i).segment<3>(0) = n;

        Eigen::Vector3d r = Eigen::Vector3d::Random();
        c1 = r.cross(n).normalized();
        c2 = n.cross(c1).normalized();

        appState->boundary_frames.row(i).segment<3>(3) = c1 * 1e-4;
        appState->boundary_frames.row(i).segment<3>(6) = c2 * 1e-4;
        // appState->boundary_frames.row(i).segment<3>(3) = c1 ;
        // appState->boundary_frames.row(i).segment<3>(6) = c2 ;
    }

    if (!appState->headless_mode) {
        polyscope::getVolumeMesh("c")->addVertexScalarQuantity("signed distance field", S);
    }
}

void Mint3DHook::setFramesToRandom(double scale) {
    for (int i = 0; i < appState->cur_tet_mesh->nBoundaryElements(); i++) {
        int btet = appState->cur_tet_mesh->boundaryElementTet(i);
        Eigen::Vector3d n = appState->frames.row(btet).segment<3>(0);
        appState->boundary_frames.row(i).segment<3>(0) = -n;

        appState->boundary_frames.row(i).segment<3>(3) = scale * Eigen::Vector3d::Random().normalized();
        appState->boundary_frames.row(i).segment<3>(6) = scale * Eigen::Vector3d::Random().normalized();
    }

    for (int i = 0; i < appState->T.rows(); i++) {
        // appState->frames.row(i).segment<3>(0) = scale * Eigen::Vector3d::Random().normalized();
        appState->frames.row(i).segment<3>(0) = scale * appState->frames.row(i).segment<3>(0).normalized();
        appState->frames.row(i).segment<3>(3) = scale * Eigen::Vector3d::Random().normalized();
        appState->frames.row(i).segment<3>(6) = scale * Eigen::Vector3d::Random().normalized();
    }
}

void Mint3DHook::initializeOtherParameters() {
    // ... implementation for initializing other parameters ...
}

void Mint3DHook::initBoundaryConditions() {}

void Mint3DHook::updateOptimizationParameters() {}

void Mint3DHook::checkAndUpdateConvergence(double decrement, double energy) {}

void Mint3DHook::finalizeIteration() {}
