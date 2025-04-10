#include "CommonFunctions.h"

#include <iostream>
#include <memory>
#include <unordered_set>

#include "CubeCover/FrameField.h"
#include "CubeCover/TetMeshConnectivity.h"
#include "Surface.h"

#include "BoysColoring.h"

#include <Eigen/Dense>
#include <unsupported/Eigen/KroneckerProduct>

std::pair<Eigen::MatrixXd, Eigen::MatrixXi> GetTetSoup(const Eigen::MatrixXd& V, const Eigen::MatrixXi& T) {
    int ntet = T.rows();
    Eigen::MatrixXd soup_V(4 * ntet, 3);
    Eigen::MatrixXi soup_T(ntet, 4);
    for (int i = 0; i < ntet; i++) {
        for (int j = 0; j < 4; j++) {
            soup_V.row(4 * i + j) = V.row(T(i, j));
            soup_T(i, j) = 4 * i + j;
        }
    }
    return {soup_V, soup_T};
}

std::pair<Eigen::MatrixXd, Eigen::MatrixXi> GetTetSoup(const Eigen::MatrixXd& V,
                                                       const CubeCover::TetMeshConnectivity& mesh) {
    int ntet = mesh.nTets();
    Eigen::MatrixXd soup_V(4 * ntet, 3);
    Eigen::MatrixXi soup_T(ntet, 4);
    for (int i = 0; i < ntet; i++) {
        for (int j = 0; j < 4; j++) {
            soup_V.row(4 * i + j) = V.row(mesh.tetVertex(i, j));
            soup_T(i, j) = 4 * i + j;
        }
    }
    return {soup_V, soup_T};
}

std::pair<Eigen::MatrixXd, Eigen::MatrixXi> GetSurfaceMeshFromTetMeshWithInterior(
    const Eigen::MatrixXd& V, const CubeCover::TetMeshConnectivity& mesh) {
    Eigen::MatrixXd soup_V;
    Eigen::MatrixXi soup_T;

    std::tie(soup_V, soup_T) = GetTetSoup(V, mesh);
    CubeCover::TetMeshConnectivity soup_mesh(soup_T);

    int nfaces = 4 * soup_T.rows();
    Eigen::MatrixXi F(nfaces, 3);

    for (int i = 0; i < soup_T.rows(); i++) {
        for (int j = 0; j < 4; j++) {
            int fid = soup_mesh.tetFace(i, j);
            F.row(4 * i + j) << soup_mesh.faceVertex(fid, 0), soup_mesh.faceVertex(fid, 1),
                soup_mesh.faceVertex(fid, 2);
        }
    }
    return {soup_V, F};
}

std::pair<Eigen::MatrixXd, Eigen::MatrixXi> GetBoundarySurfaceMeshFromTetMesh(
    const Eigen::MatrixXd& V, const CubeCover::TetMeshConnectivity& mesh) {
    // make a mesh out of all of the boundary faces
    int nbdry = 0;
    int nfaces = mesh.nFaces();
    for (int i = 0; i < nfaces; i++) {
        if (mesh.isBoundaryFace(i)) nbdry++;
    }
    Eigen::MatrixXd bdry_V(V.rows(), 3);
    Eigen::MatrixXi bdry_F(nbdry, 3);

    std::unordered_map<int, int> vid_to_bryid;

    int curidx = 0;
    int curvidx = 0;
    for (int i = 0; i < nfaces; i++) {
        if (mesh.isBoundaryFace(i)) {
            for (int j = 0; j < 3; j++) {
                int vid = mesh.faceVertex(i, j);
                if (vid_to_bryid.count(vid)) {
                    bdry_F(curidx, j) = vid_to_bryid[vid];
                } else {
                    vid_to_bryid[vid] = curvidx;
                    bdry_F(curidx, j) = curvidx;
                    bdry_V.row(curvidx++) = V.row(vid);
                }
            }
            // fix triangle orientations
            int tet = mesh.faceTet(i, 1);
            if (tet == -1) {
                std::swap(bdry_F(curidx, 0), bdry_F(curidx, 1));
            }
            curidx++;
        }
    }
    bdry_V.conservativeResize(curvidx, Eigen::NoChange);
    return {bdry_V, bdry_F};
}

/**
 * @brief Finds sharp features on a triangular mesh.
 *
 * @param V Vertices of the mesh. Each row represents a vertex, with columns representing the x, y, and z coordinates.
 * @param F Faces (triangles) of the mesh. Each row represents a face, with columns containing indices into the V
 * matrix, specifying the three vertices of the triangle.
 * @param faceEdges Mapping from faces to edges. Each row corresponds to a face, and the three columns contain indices
 * into the edgeVerts matrix, representing the three edges of the face. The order of edges is implicitly defined as (v0,
 * v1), (v1, v2), (v2, v0), where v0, v1, v2 are the vertices of the face.
 * @param faceNeighbors Neighboring faces for each face. Each row corresponds to a face, and the three columns contain
 * indices into the F matrix, representing the three neighboring faces. A value of -1 indicates a boundary edge (no
 * neighbor). The order of neighbors corresponds to the order of edges in faceEdges. Specifically, neighbor i is
 * opposite to vertex i of the face.
 * @param edgeVerts Vertices for each edge. Each row represents an edge, with columns containing indices into the V
 * matrix, specifying the two vertices of the edge.
 * @param sharp_feature_threshold Angle threshold in degrees. Edges with a dihedral angle greater than this threshold
 * are considered sharp.
 * @return A pair. The first element is a vector of vectors, where each inner vector contains the sharp edges
 * (represented by their normalized direction vectors) found on a corresponding face. The second element is a pair
 * containing the sharp feature nodes and the sharp feature edges, represented by their indices.
 */

void findSharpFeatures(const Eigen::MatrixXd& V_inp, const Eigen::MatrixXi& F_inp, double sharp_feature_threshold,
                       std::vector<Eigen::Vector3d>& nodes, std::vector<Eigen::Vector2i>& edges) {
    Eigen::MatrixXi faceEdges;
    Eigen::MatrixXi faceNeighbors;
    Eigen::MatrixXi edgeVerts;

    Surface s = Surface(V_inp, F_inp);
    Eigen::MatrixXd V = s.data().V;
    Eigen::MatrixXi F = s.data().F;

    faceEdges = s.data().faceEdges;
    faceNeighbors = s.data().faceNeighbors;
    edgeVerts = s.data().edgeVerts;

    // Calculate the number of faces
    int nBoundElem = F.rows();

    // Initialize the vector to store sharp feature edges for each face
    std::vector<std::vector<Eigen::Vector3d>> sharp_feature_edges(nBoundElem);

    // Convert the angle threshold from degrees to radians
    double thresh = sharp_feature_threshold * M_PI / 180.0;

    // Preallocate memory for nodes and edges
    // std::vector<Eigen::Vector3d> nodes;
    // std::vector<Eigen::Vector2i> edges;

    nodes.reserve(2 * edgeVerts.rows());   // Maximum possible number of sharp feature nodes
    edges.reserve(faceEdges.rows());       // Maximum possible number of sharp feature edges
    nodes.clear();
    edges.clear();

    int cedgeid = 0;

    for (int i = 0; i < nBoundElem; i++) {
        Eigen::Vector3d v0 = V.row(F(i, 0));
        Eigen::Vector3d v1 = V.row(F(i, 1));
        Eigen::Vector3d v2 = V.row(F(i, 2));

        Eigen::Vector3d e1 = v1 - v0;
        Eigen::Vector3d e2 = v2 - v0;

        Eigen::Vector3d nf = e1.cross(e2).normalized();

        std::vector<Eigen::Vector3d> face_sharp_edges;

        for (int j = 0; j < 3; j++) {
            double cedge = faceEdges(i, j);
            int g = faceNeighbors(i, j);
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

            //   << std::endl;
            if (angle > thresh) {
                int faceEdge = faceEdges(i, j);

                Eigen::Vector2i e = edgeVerts.row(faceEdge);
                Eigen::Vector3d p0 = V.row(e(0));
                Eigen::Vector3d p1 = V.row(e(1));

                if (nf.dot(ng) < 0) nf = -nf;

                p0 += (nf + ng).normalized() * 0.0001;
                p1 += (nf + ng).normalized() * 0.0001;

                Eigen::Vector3d sharp_edge = (p0 - p1).normalized();

                Eigen::Vector2i ve = Eigen::Vector2i(cedgeid, cedgeid + 1);
                edges.push_back(ve);
                cedgeid += 2;

                nodes.push_back(p0);
                nodes.push_back(p1);
            }
        }

        sharp_feature_edges.push_back(face_sharp_edges);
    }

    return;
    // return {nodes, edges};
}

Eigen::MatrixXd NormalizePts(const Eigen::MatrixXd& pos, Eigen::RowVector3d& center, double& scaling_ratio) {
    Eigen::RowVector3d min_corner, max_corner;
    min_corner = pos.colwise().minCoeff().transpose();
    max_corner = pos.colwise().maxCoeff().transpose();

    center = (min_corner + max_corner) / 2;
    scaling_ratio = (max_corner - min_corner).maxCoeff();

    return NormalizePtsGivenCenterAndScalingRatio(pos, center, scaling_ratio);
}

Eigen::MatrixXd NormalizePtsGivenCenterAndScalingRatio(const Eigen::MatrixXd& pos, const Eigen::RowVector3d& center,
                                                       double scaling_ratio) {
    int npts = pos.rows();
    Eigen::MatrixXd normalize_pts = pos;
    for (int i = 0; i < npts; i++) {
        normalize_pts.row(i) = (normalize_pts.row(i) - center) / scaling_ratio;
    }
    return normalize_pts;
}

Eigen::MatrixXd VectorPts2Matrix(const std::vector<Eigen::Vector3d>& pts) {
    int n = pts.size();
    Eigen::MatrixXd mat(n, 3);
    for (int i = 0; i < n; i++) {
        mat.row(i) = pts[i].transpose();
    }
    return mat;
}

Eigen::MatrixXd ComputeGradient(const Eigen::MatrixXd& V, const CubeCover::TetMeshConnectivity& mesh,
                                const Eigen::MatrixXd& values) {
    int nvecs = values.cols();
    int ntets = mesh.nTets();

    // assume that the values are defined on the tet soup
    auto get_soup_vid = [&](int tet_id, int j) {
        return 4 * tet_id + j;   // this is corresponding to the tet(tet_id, j)
    };

    Eigen::MatrixXd grad_fields(ntets, 3 * nvecs);

    for (int vec_id = 0; vec_id < nvecs; vec_id++) {
        for (int tet_id = 0; tet_id < ntets; tet_id++) {
            Eigen::Matrix3d A;
            Eigen::Vector3d rhs;
            for (int j = 0; j < 3; j++) {
                A.row(j) = V.row(mesh.tetVertex(tet_id, j + 1)) - V.row(mesh.tetVertex(tet_id, 0));
                rhs[j] = values(get_soup_vid(tet_id, j + 1), vec_id) - values(get_soup_vid(tet_id, 0), vec_id);
            }
            if (std::abs(A.determinant()) > 1e-10) {
                Eigen::Vector3d grad = A.inverse() * rhs;
                grad_fields.row(tet_id).segment<3>(3 * vec_id) << grad[0], grad[1], grad[2];
            } else {
                grad_fields.row(tet_id).segment<3>(3 * vec_id).setZero();
            }
        }
    }
    return grad_fields;
}

std::vector<Eigen::VectorXd> GetFrameDifference(const Eigen::MatrixXd& frame0, const Eigen::MatrixXd& frame1,
                                                int nvecs) {
    int neles = frame0.rows();
    assert(frame0.cols() == frame1.cols() && frame0.rows() == neles && frame0.cols() % 3 == 0 &&
           frame0.cols() / 3 == nvecs);
    std::vector<Eigen::VectorXd> errs(nvecs, Eigen::VectorXd::Zero(neles));
    for (int vec_id = 0; vec_id < nvecs; vec_id++) {
        for (int ele_id = 0; ele_id < neles; ele_id++) {
            Eigen::RowVector3d v0 = frame0.row(ele_id).segment<3>(3 * vec_id);

            double min_err = std::numeric_limits<double>::max();
            for (int vec_id1 = 0; vec_id1 < nvecs; vec_id1++) {
                Eigen::RowVector3d v1 = frame1.row(ele_id).segment<3>(3 * vec_id1);
                double err = (v0 - v1).norm();
                min_err = std::min(min_err, err);

                err = (v0 + v1).norm();
                min_err = std::min(min_err, err);
            }

            errs[vec_id][ele_id] = min_err;
        }
    }
    return errs;
}

std::vector<Eigen::MatrixXd> GetBestMatchFrames(const std::vector<Eigen::MatrixXd>& frame_list,
                                                const Eigen::MatrixXd& frame) {
    if (frame_list.empty()) {
        return {};
    }

    int neles = frame_list[0].rows();
    int nvecs = frame_list.size();
    std::vector<Eigen::MatrixXd> best_frame_list(nvecs, Eigen::MatrixXd::Zero(neles, 3));
    std::vector<bool> used(nvecs * neles, false);

    for (int vec_id = 0; vec_id < nvecs; vec_id++) {
        for (int ele_id = 0; ele_id < neles; ele_id++) {
            Eigen::RowVector3d v0 = frame_list[vec_id].row(ele_id);

            double min_err = std::numeric_limits<double>::max();
            int best_vec_id1 = -1;
            bool best_sign = false;

            for (int vec_id1 = 0; vec_id1 < nvecs; vec_id1++) {
                if (used[ele_id * nvecs + vec_id1]) {
                    continue;
                }
                Eigen::RowVector3d v1 = frame.row(ele_id).segment<3>(3 * vec_id1);
                double err = (v0 - v1).norm();
                if (err < min_err) {
                    best_vec_id1 = vec_id1;
                    best_sign = false;
                    min_err = err;
                }
                err = (v0 + v1).norm();
                if (err < min_err) {
                    best_vec_id1 = vec_id1;
                    best_sign = true;
                    min_err = err;
                }
            }

            if (best_vec_id1 != -1) {
                used[ele_id * nvecs + best_vec_id1] = true;
                Eigen::RowVector3d best_match;
                if (best_sign) {
                    best_match = -frame.row(ele_id).segment<3>(3 * best_vec_id1);
                } else {
                    best_match = frame.row(ele_id).segment<3>(3 * best_vec_id1);
                }
                best_frame_list[vec_id].row(ele_id) = best_match;
            }
        }
    }
    return best_frame_list;
}

std::vector<Eigen::MatrixXd> ExtractFrameVectors(const Eigen::MatrixXd& frames) {
    assert(frames.cols() % 3 == 0);
    int nvecs = frames.cols() / 3;
    int ntets = frames.rows();
    std::vector<Eigen::MatrixXd> frame_vectors(nvecs, Eigen::MatrixXd::Zero(ntets, 3));
    for (int vec_id = 0; vec_id < nvecs; vec_id++) {
        for (int tet_id = 0; tet_id < ntets; tet_id++) {
            frame_vectors[vec_id].row(tet_id) = frames.row(tet_id).segment<3>(3 * vec_id);
        }
    }
    return frame_vectors;
}

void HSV2RGB(double h, double s, double v, double& r, double& g, double& b) {
    // From medit
    double f, p, q, t, hh;
    int i;
    // shift the hue to the range [0, 360] before performing calculations
    hh = ((360 + ((int)h % 360)) % 360) / 60.;
    i = (int)std::floor(hh); /* largest int <= h     */
    f = hh - i;              /* fractional part of h */
    p = v * (1.0 - s);
    q = v * (1.0 - (s * f));
    t = v * (1.0 - (s * (1.0 - f)));

    switch (i) {
        case 0:
            r = v;
            g = t;
            b = p;
            break;
        case 1:
            r = q;
            g = v;
            b = p;
            break;
        case 2:
            r = p;
            g = v;
            b = t;
            break;
        case 3:
            r = p;
            g = q;
            b = v;
            break;
        case 4:
            r = t;
            g = p;
            b = v;
            break;
        case 5:
            r = v;
            g = p;
            b = q;
            break;
    }
}

Eigen::MatrixXd PaintPhi(const Eigen::VectorXd& phi,
                         Eigen::VectorXd* brightness)   // brightness should between 0 and 1
{
    int nverts = phi.size();
    Eigen::MatrixXd color(nverts, 3);
    for (int i = 0; i < nverts; i++) {
        double r, g, b;
        //            double h = 360.0 * phi[i] / 2.0 / M_PI + 120;
        double h = 360.0 * phi[i] / 2.0 / M_PI;
        h = 360 + ((int)h % 360);   // fix for libigl bug
        double s = 1.0;
        double v = 0.5;
        if (brightness) {
            double r = (*brightness)(i);
            v = r * r / (r * r + 1);
        }
        //                v = (*brightness)(i);
        HSV2RGB(h, s, v, r, g, b);
        color(i, 0) = r;
        color(i, 1) = g;
        color(i, 2) = b;
    }
    return color;
}

void GetStreamlines(const std::vector<CubeCover::Streamline>& traces, Eigen::MatrixXd& p_start, Eigen::MatrixXd& p_end,
                    Eigen::MatrixXd& colors) {
    if (traces.empty()) {
        return;
    }
    int ntraces = traces.size();

    p_start.resize(10000, 3);
    p_end.resize(10000, 3);
    colors.resize(10000, 4);

    int id = 0;

    for (int tid = 0; tid < ntraces; tid++) {
        int nsteps = traces.at(tid).stream_pts_.size();

        for (int i = 0; i < nsteps - 1; i++) {
            Eigen::Vector3d edge =
                traces.at(tid).stream_pts_[i].start_pt_ - traces.at(tid).stream_pts_[i + 1].start_pt_;
            Eigen::Vector3d rgb_color = Boys2RGB(edge);

            if (id >= p_start.rows()) {
                p_start.conservativeResize(p_start.rows() + 10000, Eigen::NoChange);
                p_end.conservativeResize(p_end.rows() + 10000, Eigen::NoChange);
                colors.conservativeResize(colors.rows() + 10000, Eigen::NoChange);
            }

            p_start.row(id) = traces[tid].stream_pts_[i].start_pt_.transpose();
            p_end.row(id) = traces[tid].stream_pts_[i + 1].start_pt_.transpose();
            colors.row(id) << rgb_color[0], rgb_color[1], rgb_color[2], 1;
            id++;
        }
    }
    p_start.conservativeResize(id, Eigen::NoChange);
    p_end.conservativeResize(id, Eigen::NoChange);
    colors.conservativeResize(id, Eigen::NoChange);
}

void GetStreamlines(const std::vector<CubeCover::Streamline>& traces, const std::vector<Eigen::VectorXd>& errs,
                    Eigen::MatrixXd& p_start, Eigen::MatrixXd& p_end, Eigen::VectorXd& scalar_quatity,
                    Eigen::MatrixXd* colors) {
    if (traces.empty()) {
        return;
    }
    int ntraces = traces.size();

    p_start.resize(10000, 3);
    p_end.resize(10000, 3);
    scalar_quatity.resize(10000);
    if (colors) {
        colors->resize(10000, 4);
    }

    int id = 0;

    for (int tid = 0; tid < ntraces; tid++) {
        int nsteps = traces.at(tid).stream_pts_.size();

        for (int i = 0; i < nsteps - 1; i++) {
            Eigen::Vector3d edge =
                traces.at(tid).stream_pts_[i].start_pt_ - traces.at(tid).stream_pts_[i + 1].start_pt_;
            int tet_id = traces.at(tid).tet_ids_[i];
            Eigen::Vector3d rgb_color = Boys2RGB(edge);

            if (id >= p_start.rows()) {
                p_start.conservativeResize(p_start.rows() + 10000, Eigen::NoChange);
                p_end.conservativeResize(p_end.rows() + 10000, Eigen::NoChange);
                scalar_quatity.conservativeResize(scalar_quatity.size() + 10000);
                if (colors) {
                    colors->conservativeResize(colors->rows() + 10000, Eigen::NoChange);
                }
            }

            p_start.row(id) = traces[tid].stream_pts_[i].start_pt_.transpose();
            p_end.row(id) = traces[tid].stream_pts_[i + 1].start_pt_.transpose();

            if (colors) {
                colors->row(id) << rgb_color[0], rgb_color[1], rgb_color[2], 1;
            }

            double max_err = 0;
            for (int j = 0; j < errs.size(); j++) {
                max_err = std::max(max_err, errs[j][tet_id]);
            }
            scalar_quatity(id) = max_err;
            id++;
        }
    }
    p_start.conservativeResize(id, Eigen::NoChange);
    p_end.conservativeResize(id, Eigen::NoChange);
    if (colors) {
        colors->conservativeResize(id, Eigen::NoChange);
    }
    scalar_quatity.conservativeResize(id);
}

std::tuple<Eigen::MatrixXd, Eigen::MatrixXd, Eigen::MatrixXd> ConvertCurveNetWorkForRender(
    const Eigen::MatrixXd& P, const Eigen::MatrixXi& E, const Eigen::RowVector4d& const_color) {
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

double DihedralAngle(const Eigen::MatrixXd& V, const Eigen::MatrixXi& F, int f0, int f1) {
    Eigen::RowVector3d e0 = V.row(F(f0, 1)) - V.row(F(f0, 0));
    Eigen::RowVector3d e1 = V.row(F(f0, 2)) - V.row(F(f0, 0));
    Eigen::RowVector3d n0 = e0.cross(e1);

    e0 = V.row(F(f1, 1)) - V.row(F(f1, 0));
    e1 = V.row(F(f1, 2)) - V.row(F(f1, 0));
    Eigen::RowVector3d n1 = e0.cross(e1);

    double cos = n0.dot(n1) / n0.norm() / n1.norm();
    double theta = std::acos(std::clamp(cos, -1.0, 1.0));

    return M_PI - theta;
}

bool IsValidMesh(Eigen::MatrixXd& tet_pos, CubeCover::TetMeshConnectivity& tet_mesh, bool verbose) {
    bool is_orientable = true;
    for (int i = 0; i < tet_mesh.nTets(); i++) {
        Eigen::Matrix3d A;
        A.row(0) = tet_pos.row(tet_mesh.tetVertex(i, 1)) - tet_pos.row(tet_mesh.tetVertex(i, 0));
        A.row(1) = tet_pos.row(tet_mesh.tetVertex(i, 2)) - tet_pos.row(tet_mesh.tetVertex(i, 0));
        A.row(2) = tet_pos.row(tet_mesh.tetVertex(i, 3)) - tet_pos.row(tet_mesh.tetVertex(i, 0));

        if (A.determinant() < 0) {
            if (verbose) {
                std::cout << "tet " << i << " is inverted, vol: " << A.determinant() << std::endl;
            }
            is_orientable = false;
        }
    }

    if (!is_orientable) {
        std::cout << "Not orientable mesh" << std::endl;
        return false;
    }

    if (!tet_mesh.isManifold(true)) {
        std::cout << "Not manifold mesh" << std::endl;
        return false;
    }
    return true;
}

bool IsValidFrame(const CubeCover::TetMeshConnectivity& tet_mesh, const Eigen::MatrixXd& tet_frames, bool verbose) {
    using namespace CubeCover;
    Eigen::MatrixXi assiangments;
    std::unique_ptr<FrameField> frame_fields(fromFramesAndAssignments(tet_mesh, tet_frames, assiangments, false));
    frame_fields->computeLocalAssignments();

    int n_singular_edges = frame_fields->nSingularEdges();
    if (verbose) {
        std::cout << "the frames have " << n_singular_edges << " singular edges" << std::endl;
    }
    std::unordered_set<int> edge_set;
    for (int i = 0; i < n_singular_edges; i++) {
        edge_set.insert(frame_fields->singularEdge(i));
    }

    bool is_valid = true;
    for (int i = 0; i < tet_mesh.nTets(); i++) {
        int count = 0;
        for (int j = 0; j < 6; j++) {
            int eid = tet_mesh.tetEdge(i, j);
            if (edge_set.count(eid)) {
                count++;
            }
        }
        if (count >= 2) {
            if (verbose) {
                std::cout << "tet " << i << " has " << count << " singular edges" << std::endl;
            }
            is_valid = false;
        }
    }

    return is_valid;
}

Eigen::MatrixXi GetTetMatrix(const CubeCover::TetMeshConnectivity& tet_mesh) {
    Eigen::MatrixXi T(tet_mesh.nTets(), 4);
    for (int i = 0; i < tet_mesh.nTets(); i++) {
        for (int j = 0; j < 4; j++) {
            T(i, j) = tet_mesh.tetVertex(i, j);
        }
    }
    return T;
}

Eigen::MatrixXd ExpandFramePermutation(const Eigen::MatrixXd& perm, const Eigen::MatrixXd& sgn) {
    Eigen::MatrixXd P = Eigen::MatrixXd::Zero(9, 9);

    Eigen::Matrix3d I3 = Eigen::Matrix3d::Identity();
    Eigen::Matrix3d to_exp = perm * sgn;
    P = Eigen::kroneckerProduct(to_exp, I3);

    return P;
}

// Function to compute the optimal permutation matrix
double ComputeCombedL2wrtBasis(const Eigen::VectorXd& f, const Eigen::VectorXd& g, const Eigen::MatrixXd& basis) {
    Eigen::MatrixXd best_P(9, 9);
    double min_distance = std::numeric_limits<double>::max();

    // Define permutation matrices
    std::vector<Eigen::Matrix3d> perm_matrices(6);
    std::vector<Eigen::Matrix3d> sgn_matrices(8);

    Eigen::Matrix3d I3 = Eigen::Matrix3d::Identity();

    // Generate all 3! permutations
    perm_matrices[0] << 1, 0, 0, 0, 1, 0, 0, 0, 1;
    perm_matrices[1] << 0, 1, 0, 1, 0, 0, 0, 0, 1;
    perm_matrices[2] << 0, 0, 1, 1, 0, 0, 0, 1, 0;
    perm_matrices[3] << 1, 0, 0, 0, 0, 1, 0, 1, 0;
    perm_matrices[4] << 0, 1, 0, 0, 0, 1, 1, 0, 0;
    perm_matrices[5] << 0, 0, 1, 0, 1, 0, 1, 0, 0;

    // All the sign flips
    sgn_matrices[0] << 1, 0, 0, 0, 1, 0, 0, 0, 1;
    sgn_matrices[1] << -1, 0, 0, 0, 1, 0, 0, 0, 1;
    sgn_matrices[2] << 1, 0, 0, 0, -1, 0, 0, 0, 1;
    sgn_matrices[3] << -1, 0, 0, 0, -1, 0, 0, 0, 1;
    sgn_matrices[4] << 1, 0, 0, 0, 1, 0, 0, 0, -1;
    sgn_matrices[5] << -1, 0, 0, 0, 1, 0, 0, 0, -1;
    sgn_matrices[6] << 1, 0, 0, 0, -1, 0, 0, 0, -1;
    sgn_matrices[7] << -1, 0, 0, 0, -1, 0, 0, 0, -1;

    // Iterate through all permutations
    for (const auto& perm : perm_matrices) {
        for (const auto& sgn : sgn_matrices) {
            // Eigen::MatrixXd P = Eigen::kroneckerProduct(sgn, perm);
            Eigen::MatrixXd P = ExpandFramePermutation(perm, sgn);
            Eigen::MatrixXd g_perm = P * g;
            Eigen::VectorXd diff = f - g_perm;
            double distance = 0;
            for (int i = 0; i < 3; i++) {
                distance += (basis * diff.segment<3>(3 * i)).norm();
            }

            if (distance < min_distance) {
                min_distance = distance;
                // std::cout << perm << " perm \n" << sgn << " sgn " << std::endl;
                // std::cout << P << " P \n" << std::endl;
                best_P = P;
            }
        }
    }

    return min_distance;
}

void ComputePerTetFacetBasis(const Eigen::MatrixXd& V, const CubeCover::TetMeshConnectivity& mesh,
                             std::vector<std::vector<Eigen::MatrixXd>>& tet_facet_basis_,
                             std::vector<std::vector<Eigen::MatrixXd>>& tet_facet_dual_basis_,
                             std::vector<std::vector<Eigen::MatrixXd>>& tet_facet_full_basis_) {
    int ntets = mesh.nTets();
    tet_facet_basis_.resize(ntets);
    tet_facet_dual_basis_.resize(ntets);
    tet_facet_full_basis_.resize(ntets);

    for (int tetidx = 0; tetidx < ntets; tetidx++) {
        std::vector<Eigen::MatrixXd> cur_tet_bases(4);
        std::vector<Eigen::MatrixXd> cur_tet_dual_bases(4);
        std::vector<Eigen::MatrixXd> cur_tet_full_bases(4);

        for (int idx = 0; idx < 4; idx++) {
            int face_idx = mesh.tetFace(tetidx, idx);

            Eigen::Matrix3d rot_facet_to_template;
            Eigen::MatrixXd cur_facet_full_basis = Eigen::MatrixXd::Zero(3, 3);

            Eigen::Vector3d a = V.row(mesh.faceVertex(face_idx, 0));
            Eigen::Vector3d b = V.row(mesh.faceVertex(face_idx, 1));
            Eigen::Vector3d c = V.row(mesh.faceVertex(face_idx, 2));

            Eigen::Vector3d b1 = (b - a).normalized();
            Eigen::Vector3d b2 = (c - a).normalized();
            Eigen::Vector3d n = b1.cross(b2).normalized();

            cur_facet_full_basis.row(0) = b1;
            cur_facet_full_basis.row(1) = b1.cross(n);   // b2
            cur_facet_full_basis.row(2) = n;

            cur_tet_bases[idx] = cur_facet_full_basis.block<2, 3>(0, 0);
            cur_tet_dual_bases[idx] = n.transpose();
            cur_tet_full_bases[idx] = cur_facet_full_basis;
        }
        tet_facet_basis_[tetidx] = cur_tet_bases;
        tet_facet_dual_basis_[tetidx] = cur_tet_dual_bases;
        tet_facet_full_basis_[tetidx] = cur_tet_full_bases;
    }
}