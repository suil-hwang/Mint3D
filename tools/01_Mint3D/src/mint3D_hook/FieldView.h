#ifndef FIELDVIEW_H
#define FIELDVIEW_H

#include <string>

/**
 * Enumeration for the different field views.
 */

namespace Views {

enum class Field_View {
    vec_norms,
    moment_dirch,
    sym_curl_residual,
    odeco,
    fradeco,
    frame_det,
    delta_norms,
    vec_dirch,
    primal_curl_residual,
    gui_free,
    Element_COUNT
};

enum class Sym_Moment_View { Total, L2, L4, L6, Element_COUNT };

enum class Sym_Curl_View { Total, L2, L4, L6, Element_COUNT };

// enum Field_View {
//     vec_norms,
//     delta_norms,
//     vec_dirch,
//     moment_dirch,
//     primal_curl_residual,
//     curl_l2_residual,
//     curl_l4_residual,
//     gui_free,
//     Element_COUNT
// };

/**
 * Converts a field view enumeration to a string.
 *
 * @param view The field view to convert.
 * @return The name of the field view as a string.
 */
std::string fieldViewToString(Field_View view);
std::string fieldViewToFileStub(Field_View view);

};   // namespace Views

#endif   // FIELDVIEW_H
