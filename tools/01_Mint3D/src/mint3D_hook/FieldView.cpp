#include "FieldView.h"

// Utility function to convert Field_View enums to strings
std::string Views::fieldViewToString(Field_View view) {
    switch (view) {
        case Views::Field_View::vec_norms:
            return "Vector Norms";
        case Views::Field_View::delta_norms:
            return "Delta Norms";
        case Views::Field_View::vec_dirch:
            return "Vec Dirichlet";
        case Views::Field_View::moment_dirch:
            return "Moment Dirichlet";
        case Views::Field_View::primal_curl_residual:
            return "Primal Curl Residual";
        case Views::Field_View::sym_curl_residual:
            return "Symmetric Curl Residual";
        case Views::Field_View::odeco:
            return "Orthogonality Violation";
        case Views::Field_View::fradeco:
            return "Unit Frame Violation";
        case Views::Field_View::frame_det:
            return "Violation from Unit Volume (not implemented?)";
        case Views::Field_View::gui_free:
            return "Free";
        default:
            return "Unknown";
    }
}

// Utility function to convert Field_View enums to strings
std::string Views::fieldViewToFileStub(Field_View view) {
    switch (view) {
        case Views::Field_View::vec_norms:
            return "nv_norms_vec";
        case Views::Field_View::delta_norms:
            return "nd_norms_delta";
        case Views::Field_View::vec_dirch:
            return "dv_dirch_vec";
        case Views::Field_View::moment_dirch:
            return "dm_dirch_moment";
        case Views::Field_View::primal_curl_residual:
            return "cv_curl_vec";
        case Views::Field_View::sym_curl_residual:
            return "cm_curl_mom";
        case Views::Field_View::odeco:
            return "on_odeco";
        case Views::Field_View::fradeco:
            return "on_fradeco";
        case Views::Field_View::frame_det:
            return "fd_frame_det";
        case Views::Field_View::gui_free:
            return "free";
        default:
            return "Unknown";
    }
}