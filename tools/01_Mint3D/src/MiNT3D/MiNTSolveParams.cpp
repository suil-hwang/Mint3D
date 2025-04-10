

// // I want to 
// #include "MiNTSolveParams.h"

// #include <Eigen/Dense>
// #include <fstream>
// #include <iostream>
// #include <nlohmann/json.hpp>

// //    bool serializeSolveParams(MiNT3D::SolveParams& params, const std::string& filepath);
//     // bool deserializeSolveParams(MiNT3D::SolveParams& params, const std::string& filepath);


// // I want to implement these  




// bool MiNTSolveParams::serializeSolveParams(MiNT3D::SolveParams& params, const std::string& filepath) {
//     try {
//         json j;

//         j["w_bound"] = params.w_bound;
//         j["w_mint"] = params.w_mint;
//         j["w_smooth"] = params.w_smooth;
//         j["w_scaled_jacobian"] = params.w_scaled_jacobian;
//         j["w_unit_norm"] = params.w_unit_norm;
//         j["w_unit_barrier"] = params.w_unit_barrier;
//         j["w_fit"] = params.w_fit;
//         j["w_global_rescale_param"] = params.getGlobalScale();

//         j["w_outer_step"] = params.w_outer_step;

//         j["inner_iter_max_steps"] = params.inner_iter_max_steps;
//         j["grad_tol"] = params.grad_tol;
//         j["xTol"] = params.xTol;
//         j["fTol"] = params.fTol;
//         j["reg"] = params.reg;

//         j["curr_total_step"] = params.total_step;
//         j["curr_outer_step"] = params.outer_step;

//         j["lambda_penalty"] = params.lambda_penalty;
//         j["use_kruskal_tensor_in_mint_term"] = params.use_kruskal_tensor_in_mint_term;


//         // data log stuff, mostly for debugging

//         j["energy_smoothness"] = params.smoothness_energy;
//         j["energy_mint"] = params.mint_energy;
//         j["energy_unit"] = params.unit_energy;
//         j["energy_primal_integrability"] = params.primal_integrability_energy;
//         j["energy_scaled_jacobian"] = params.scaled_jacobian_energy;
//         j["energy_unit_barrier"] = params.unit_barrier_energy;
//         j["energy_total"] = params.total_energy;

//         std::ofstream outFile(filepath);
//         if (!outFile.is_open()) {
//             std::cerr << "Error: Unable to open file for writing JSON: " << filepath << std::endl;
//             return false;
//         }

//         outFile << j.dump(4);   // Pretty-print with 4 spaces of indentation
//         outFile.close();
//         return true;
//     } catch (const std::exception& e) {
//         std::cerr << "Error: Unable to serialize config to JSON: " << e.what() << std::endl;
//         return false;
//     }
// }

// bool MiNTSolveParams::deserializeSolveParams(MiNT3D::SolveParams& params, const std::string& filepath) {
//     try {
//         std::ifstream inFile(filepath);
//         if (!inFile.is_open()) {
//             std::cerr << "Error: Unable to open file for deserialization: " << filepath << std::endl;
//             return false;
//         }

//         json j;
//         inFile >> j;

//         params = params.setGlobalScale(j.at("w_global_rescale_param"));
//         // the stored weights have the global scale applied

//         params.w_bound = j.at("w_bound");
//         params.w_bound_viscosity = j.at("w_bound_viscosity");
//         params.w_mint = j.at("w_mint");
//         params.w_smooth = j.at("w_smooth");
//         params.w_scaled_jacobian = j.at("w_scaled_jacobian");
//         params.w_unit_norm = j.at("w_unit_norm");
//         params.w_unit_barrier = j.at("w_unit_barrier");
//         params.w_fit = j.at("w_fit");

//         params.w_outer_step = j.at("w_outer_step");
//         params.w_attenuate = j.at("w_attenuate");

//         params.inner_iter_max_steps = j.at("inner_iter_max_steps");
//         params.grad_tol = j.at("grad_tol");
//         params.xTol = j.at("xTol");
//         params.fTol = j.at("fTol");
//         params.reg = j.at("reg");

//         params.total_step = j.at("curr_total_step");
//         params.outer_step = j.at("curr_outer_step");

//         params.lambda_penalty = j.at("lambda_penalty");
//         params.use_kruskal_tensor_in_mint_term = j.at("use_kruskal_tensor_in_mint_term");

//         return true;
//     } catch (const std::exception& e) {
//         std::cerr << "Error: Unable to write data to CSV: " << e.what() << std::endl;
//         return false;
//     }
// }