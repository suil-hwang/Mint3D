#ifndef APPSTATE_H
#define APPSTATE_H

#include <Eigen/Core>
#include <Eigen/Dense>
#include <Eigen/Sparse>

#include <limits>
#include <optional>
#include <string>
#include <unordered_map>
#include <vector>

#include "FieldView.h"

#include "CubeCover/FrameField.h"
#include "CubeCover/TetMeshConnectivity.h"
#include "Surface.h"

#include "../MiNT3D/MiNTSolveParams.h"
// #include "MiNTEnergy.h"

#include <memory>

// #include <chrono>

using Field_View = Views::Field_View;

// Enum for identifying field view quantities
// enum Field_View {
//     vec_norms, delta_norms, vec_dirch, moment_dirch, primal_curl_residual, sym_curl_residual, gui_free, Element_COUNT
// };
// Views::

// enum Field_View : unsigned int;

// Struct to hold bounds for each field view quantity
struct FieldBounds {
    float upper = .9;   // std::numeric_limits<float>::max();
    float lower = .1;   // std::numeric_limits<float>::lowest();

    // std::pair<double, double> getBoundsPair() {
    //     return std::make_pair(lower, upper);
    // }
};

class OutputState {
   public:
    std::vector<Eigen::MatrixXd> frames;
    std::vector<Eigen::MatrixXd> boundary_frames;

    // derived variables
    double cur_global_objective_val;
    // Eigen::MatrixXd renderFrames;
    Eigen::VectorXd norms_vec;
    Eigen::VectorXd norms_delta;
    Eigen::VectorXd norms_moment;   // TODO IMPLEMENT THIS!
    Eigen::VectorXd thetas;         // TODO  IMPORTANT!
    Eigen::VectorXd curls_primal;

    // Current curl component of the energy, this should be a constraint that
    // is close to numerically zero at convergence.
    Eigen::VectorXd curls_sym;
    Eigen::VectorXd smoothness_primal;

    // Current smoothness energy
    Eigen::VectorXd smoothness_sym;

    // Individual components
    Eigen::VectorXd smoothness_L2;
    Eigen::VectorXd smoothness_L4;
    Eigen::VectorXd smoothness_L6;
    // Eigen::VectorXd smoothness_L2x2;

    Eigen::VectorXd E_odeco;
    Eigen::VectorXd E_fradeco;
    Eigen::VectorXd E_unit_boundary_volume;

    // Individual components
    Eigen::VectorXd curl_L2;
    Eigen::VectorXd curl_L4;
    std::vector<Eigen::VectorXd> curls_Lks;
    // Eigen::VectorXd curl_L2x2;
};

enum class VariableType { primals, moments, deltas, gammas, Element_COUNT };

// AppState holds the state of the application
class AppState {
   public:
    ~AppState() {}

    // File IO state
    std::string directoryPath;   // this is where files are loaded from
    std::string
        logFolderPath;   // this is where the output gets saved.
                         // When load from a directory logFolderPath defaults to same directory, can change this.
    std::string logOuterIterPath;

    std::string meshName;
    std::string meshInputFilePath;
    std::string exp_name;
    std::string out_dir = "../../results/";
    // std::vector<std::string> bfraFiles;
    // std::vector<std::string> bmomFiles;
    std::optional<std::string> objFilePath;
    std::optional<std::string> meshFilePath;

    int currentFileID = -1;
    bool loadLastFile = true;   // urg hacking this in.
    bool shouldReload =
        true;   // this is a dynamic var which tells updaterendergeometry to reload data from the current directory
    bool restartFromCurrentFrame = false;
    bool updateRenderGeometryNextFrameIfPaused = false;

    // Init variables
    Eigen::MatrixXd V;   // Vertex positions
    Eigen::MatrixXi T;   // Tet indices
    Eigen::MatrixXi F;   // Face indices

    bool useBoundaryFrames = true;
    int nelem;   // number of elements, i.e. T.rows() +
    int ntets;
    int nbelem;

    // virtual Eigen::VectorXd dofs_to_vars(VariableType t); // TODO: implement this.
    // This would be a map that's used inside of OptZoo functions to make things more generic.

    // Per solve metadata
    std::string solveType = "SOLVER_TYPE_NOT_SET";
    std::string solveDescription = "SOLVER_DESCRIPTION_NOT_SET";

    std::string solveStatus = "STATUS_WAS_NOT_SET";
    std::string problemFileTag = "";

    // currently unused, but maybe useful to time stuff later.
    // std::chrono::time_point<std::chrono::high_resolution_clock> opt_step_start_time;

    // Optimization variables
    std::unique_ptr<Surface>
        cur_surf;   // This initializes some more convenient data structures for building up local energies.

    std::vector<std::vector<Eigen::Vector3d>> sharp_feature_edges;
    double sharp_feature_threshold = 45; //9;   // e.g.

    std::unique_ptr<CubeCover::TetMeshConnectivity> cur_tet_mesh;
    Eigen::MatrixXd tet_centroids;     // In 3d need to use mesh data structures instead of surface.
    Eigen::MatrixXd bound_centroids;   // In 3d need to use mesh data structures instead of surface.
    Eigen::MatrixXd bound_normals;
    Eigen::MatrixXd
        bound_b1;   // a bit ugly organization, sorry.  Maybe nicer to have a std::vector of both, can change later
    Eigen::MatrixXd bound_b2;
    Eigen::VectorXi pinned_frame_idx;

    Eigen::MatrixXd pinned_orthog_b1;   // a bit ugly organization, sorry.  Maybe nicer to have a std::vector of both,
                                        // can change later
    Eigen::MatrixXd pinned_orthog_b2;

    bool headless_mode = false;
    bool keepSolving = true;

    bool subdivide_boundary_tets = true;   // false;

    ///////////////////
    /////  Solve Metrics
    ///////////////////
    int outerLoopIteration = 0;
    double cur_rel_residual = 0;
    double cur_abs_residual = 0;
    double cur_max_gradient_norm = 0;
    double cur_step_progress = 0;
    double cur_step_time = 0;
    double solve_residual = 0;
    double rhs_norm = 0;
    double solve_rel_residual = 0;
    double identity_weight = 0;
    // double obj_smoothness = 0;
    // double obj_curl = 0;
    // double obj_other_parts = 0;

    ///////////////////
    // Solver Log Stats
    ///////////////////

    std::vector<double> total_time;
    std::vector<double> assembly_time;
    std::vector<double> solve_time;
    std::vector<double> line_search_time;

    std::vector<double> identity_weight_log;

    std::vector<double> total_energy_log;
    std::vector<double> energy_diff_log;
    std::vector<double> solve_residual_log;
    std::vector<double> global_scale_log;
    std::vector<double> global_scale_outiter_log;
    std::vector<double> gradient_norm_log;
    std::vector<double> gradient_norm_step_start_log;
    std::vector<double> fit_penalty_log;

    ///////////////////

    ///////////////////
    // Energy Log Stats
    ///////////////////

    std::vector<double> smoothness_log;
    std::vector<double> unit_penalty_log;
    std::vector<double> symmetric_integrability_log;
    std::vector<double> primal_integrability_log;
    std::vector<double> combed_integrability_log;
    std::vector<double> combed_smoothness_log;

    std::vector<double> asap_combed_smoothness_log;
    std::vector<double> aiap_combed_smoothness_log;

    std::vector<double> scaled_jacobian_log;

    std::vector<double> curl_weight_log;

    ///////////////////

    MiNT3D::SolveParams solve_params;

    ///////////////////
    // frame reloading
    ///////////////////
    int fmin = 0;
    int fmax = 0;
    int fcurr = 0;

    bool loadInner = true;

    std::string experiment_name;

    ///////////////////

    std::vector<float> energy_trace;
    std::vector<float> energy_smoothness_part_trace;
    std::vector<float> energy_curl_part_trace;
    std::vector<float> energy_odeco_part_trace;
    std::vector<float> energy_fradeco_part_trace;
    std::vector<float> energy_bound_vol_part_trace;
    std::vector<float> smoothness_trace;
    std::vector<float> curl_penalty_trace;
    std::vector<float> identity_weight_trace;
    std::vector<float> timing_trace;
    void readLogFileIntoVector(const std::string& fileName, std::vector<float>& vec);
    void readAllLogFiles();

    std::vector<float> solve_rel_residual_trace;
    std::vector<float> cur_max_gradient_norm_trace;
    // std::vector<float> energy_trace;

    /// Singularities
    bool recompute_singularities = true;
    bool use_asap_combing_for_singularities = true;

    std::unique_ptr<CubeCover::FrameField> field;
    Eigen::MatrixXd Pblack;
    Eigen::MatrixXi Eblack;
    Eigen::MatrixXd Pblue;
    Eigen::MatrixXi Eblue;
    Eigen::MatrixXd Pgreen;
    Eigen::MatrixXi Egreen;

    std::vector<Eigen::Vector3d> sharp_nodes;
    std::vector<Eigen::Vector2i> sharp_edges;

    ////////////////////////
    //////  Solve State
    ////////////////////////

    // TODO: merge in the cube cover stuff.
    Eigen::MatrixXd frames;
    Eigen::MatrixXd frames_orig;
    Eigen::MatrixXd boundary_frames;
    Eigen::MatrixXd boundary_frames_orig;
    // Eigen::MatrixXd frames_viz;

    Eigen::MatrixXd moments;                        // TODO implement this!
    std::vector<Eigen::MatrixXd> frame_jacobians;   // TODO implement this!
    Eigen::MatrixXd deltas;

    // Boundary Selection Matrix
    Eigen::SparseMatrix<double> boundary_selection_matrix;   // TODO: implement this!

    // GUI state
    std::unordered_map<Field_View, FieldBounds> fieldBounds;
    // bool fieldViewActive [8] = {true, true, true, true, true, true, true, false};
    bool fieldViewActive[10] = {true, true, true, false, false, false, false, false, false, false};

    bool shouldLogData = true;
    Field_View prev_frame_element = Field_View::Element_COUNT;
    Field_View current_element;
    bool showVectorField = true;
    FieldBounds override_bounds;   //= {0.0, 1e-4};
    bool override_bounds_active = false;
    bool show_frames = false;
    bool show_frames_as_lines = true;
    bool loadedPreviousRun = false;
    int max_saved_index = 0;

    bool invert_frames_in_cylinder_test = false;
    double cylinder_multiplier = 1.;
    double noise_multiplier = 0.;

    int solve_steps = 500;   // 300

    // double L4_alpha = 0;
    float gui_vec_size = .1;
    Views::Sym_Moment_View cur_moment_view = Views::Sym_Moment_View::Total;
    Views::Sym_Curl_View cur_curl_view = Views::Sym_Curl_View::Total;

    bool LogToFile(const std::string suffix);   // Log based on fieldViewActive state
    void LogCurrentOptStats();

    // Constructor
    AppState();

    // Methods for managing AppState
    void refreshFileLists();
    void selectFile(const std::string& filename);
    void refreshData();
    void updateSliderValue(int value);
    void serializeData();
    void deserializeData();
    void zeroPassiveVars();
    void setSparseMetricFromWeights(Eigen::SparseMatrix<double>& M, const std::vector<double> weights);

    // Additional methods as needed
};

#endif   // APPSTATE_H