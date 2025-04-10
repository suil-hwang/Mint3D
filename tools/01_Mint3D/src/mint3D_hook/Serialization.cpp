#include "Serialization.h"
#include <Eigen/Dense>
#include <fstream>
#include <iostream>
#include <nlohmann/json.hpp>
// #include "polyscope/deps/json/include/json.hpp"

using json = nlohmann::json;

// Serialize Eigen vector to a binary file
bool Serialization::serializeVector(const Eigen::VectorXd& vec, const std::string& filepath) {
    try {
        std::ofstream outFile(filepath, std::ios::binary);
        if (!outFile.is_open()) {
            std::cerr << "Error: Unable to open file for writing: " << filepath << std::endl;
            return false;
        }

        // Write vector size
        int size = static_cast<int>(vec.size());
        outFile.write(reinterpret_cast<char*>(&size), sizeof(int));

        // Write vector data
        outFile.write(reinterpret_cast<const char*>(vec.data()), size * sizeof(double));
        outFile.close();
        return true;
    } catch (const std::exception& e) {
        std::cerr << "Error: Unable to serialize vector: " << e.what() << std::endl;
        return false;
    }
}

// Deserialize Eigen vector from a binary file
bool Serialization::deserializeVector(Eigen::VectorXd& vec, const std::string& filepath) {
    try {
        std::ifstream inFile(filepath, std::ios::binary);
        if (!inFile.is_open()) {
            std::cerr << "Error: Unable to open file for reading: " << filepath << std::endl;
            return false;
        }

        // Read vector size
        int size;
        inFile.read(reinterpret_cast<char*>(&size), sizeof(int));

        // Read vector data
        vec.resize(size);
        inFile.read(reinterpret_cast<char*>(vec.data()), size * sizeof(double));
        inFile.close();
        return true;
    } catch (const std::exception& e) {
        std::cerr << "Error: Unable to deserialize vector: " << e.what() << std::endl;
        return false;
    }
}

// Serialize Eigen matrix to a binary file
bool Serialization::serializeMatrix(const Eigen::MatrixXd& mat, const std::string& filepath, int vector_per_element) {
    try {
        std::ofstream outFile(filepath, std::ios::binary);
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

// Deserialize Eigen matrix from a binary file
bool Serialization::deserializeMatrix(Eigen::MatrixXd& mat, const std::string& filepath) {
    try {
        std::ifstream inFile(filepath, std::ios::binary);
        if (!inFile.is_open()) {
            std::cerr << "Error: Unable to open file for reading: " << filepath << std::endl;
            return false;
        }

        char* filename = new char[5];
        inFile.read(reinterpret_cast<char*>(&filename), sizeof("FRA 2"));
        // Read matrix rows and cols
        int rows, cols, vpe;
        inFile.read(reinterpret_cast<char*>(&rows), sizeof(int));
        inFile.read(reinterpret_cast<char*>(&cols), sizeof(int));
        inFile.read(reinterpret_cast<char*>(&vpe), sizeof(int));

        // Read matrix data
        mat.resize(rows, cols);
        inFile.read(reinterpret_cast<char*>(mat.data()), rows * cols * sizeof(double));
        inFile.close();
        return true;
    } catch (const std::exception& e) {
        std::cerr << "Error: Unable to deserialize matrix: " << e.what() << std::endl;
        return false;
    }
}

bool Serialization::deserializeFF3FramesFromMetricDrivenFrames3D(Eigen::MatrixXd& mat, const std::string& filepath) {
    try {
        std::ifstream inFile(filepath);
        if (!inFile.is_open()) {
            std::cerr << "Error: Unable to open file for reading: " << filepath << std::endl;
            return false;
        }

        std::string line;
        int rows, cols;

        // Read the first line to get the dimensions of the matrix
        if (std::getline(inFile, line)) {
            std::istringstream iss(line);
            if (!(iss >> cols >> rows)) {
                std::cerr << "Error reading matrix dimensions" << std::endl;
                return false;
            }
        }

        // Initialize an Eigen matrix with the given dimensions
        mat = Eigen::MatrixXd::Zero(rows, cols);
        std::cout << "rows: " << rows << ", cols: " << cols << std::endl;
        int currentRow = 0;
        while (std::getline(inFile, line)) {
            std::istringstream iss(line);
            // std::cout << line << std::endl;
            for (int col = 0; col < cols; ++col) {
                double value;
                if (!(iss >> value)) {
                    std::cerr << "Error reading matrix data at row " << currentRow << ", column " << col << std::endl;
                    return false;
                }
                mat(currentRow, col) = value;
            }
            currentRow++;
        }

        inFile.close();

        return true;
    } catch (const std::exception& e) {
        std::cerr << "Error: Unable to deserialize matrix: " << e.what() << std::endl;
        return false;
    }
}

bool Serialization::serializeSolveParams(MiNT3D::SolveParams& params, const std::string& filepath) {
    try {
        json j;

        // Serialize energy weights
        j["w_smooth_combed"] = params.w_smooth_combed;
        j["w_smooth_combed_precon"] = params.w_smooth_combed_precon;
        j["w_smooth_sym"] = params.w_smooth_sym;
        j["w_smooth_sym_precon"] = params.w_smooth_sym_precon;
        j["w_int_combed"] = params.w_int_combed;
        j["w_int_combed_fixed_weight"] = params.w_int_combed_fixed_weight;
        j["w_int_sym"] = params.w_int_sym;
        j["w_int_sym_fixed_weight"] = params.w_int_sym_fixed_weight;
        j["w_orthog"] = params.w_orthog;
        j["w_orthog_precon"] = params.w_orthog_precon;
        j["w_unit"] = params.w_unit;
        j["w_unit_precon"] = params.w_unit_precon;
        j["w_fit"] = params.w_fit;
        j["w_fit_precon"] = params.w_fit_precon;
        j["w_self_align"] = params.w_self_align;
        j["w_init_scale"] = params.w_init_scale;

        // Serialize booleans
        j["b_smooth_combed"] = params.b_smooth_combed;
        j["b_smooth_sym"] = params.b_smooth_sym;
        j["b_smooth_asap_combing"] = params.b_smooth_asap_combing;
        j["b_int_combed"] = params.b_int_combed;
        j["b_int_sym"] = params.b_int_sym;
        j["b_orthog"] = params.b_orthog;
        j["b_unit"] = params.b_unit;
        j["b_fit"] = params.b_fit;

        // Serialize boundary conditions and options
        j["boundary_condition"] = params.boundary_condition;
        j["boundary_hard_constraints"] = params.boundary_hard_constraints;
        j["init_state"] = params.init_state;
        j["connection_for_combing"] = params.connection_for_combing;
        j["fit_type"] = params.fit_type;

        j["w_bound"] = params.w_bound;
        j["w_mint"] = params.w_mint;
        j["w_smooth"] = params.w_smooth;
        j["w_scaled_jacobian"] = params.w_scaled_jacobian;
        j["w_unit_norm"] = params.w_unit_norm;
        j["w_unit_barrier"] = params.w_unit_barrier;
        j["w_fit"] = params.w_fit;
        j["w_global_rescale_param"] = params.getGlobalScale();

        j["w_outer_step"] = params.w_outer_step;

        j["inner_iter_max_steps"] = params.inner_iter_max_steps;
        j["grad_tol"] = params.grad_tol;
        j["xTol"] = params.xTol;
        j["fTol"] = params.fTol;
        j["reg"] = params.reg;

        j["curr_total_step"] = params.total_step;
        j["curr_outer_step"] = params.outer_step;

        j["lambda_penalty"] = params.lambda_penalty;
        j["b_use_kruskal_tensors_for_sym_smoothness"] = params.b_use_kruskal_tensors_for_sym_smoothness;

        // data log stuff, mostly for debugging

        j["energy_smoothness"] = params.smoothness_energy;
        j["energy_mint"] = params.mint_energy;
        j["energy_unit"] = params.unit_energy;
        j["energy_primal_integrability"] = params.primal_integrability_energy;
        j["energy_scaled_jacobian"] = params.scaled_jacobian_energy;
        j["energy_unit_barrier"] = params.unit_barrier_energy;
        j["energy_total"] = params.total_energy;

        std::ofstream outFile(filepath);
        if (!outFile.is_open()) {
            std::cerr << "Error: Unable to open file for writing JSON: " << filepath << std::endl;
            return false;
        }

        outFile << j.dump(4);   // Pretty-print with 4 spaces of indentation
        outFile.close();
        return true;
    } catch (const std::exception& e) {
        std::cerr << "Error: Unable to serialize config to JSON: " << e.what() << std::endl;
        return false;
    }
}

bool Serialization::deserializeSolveParams(MiNT3D::SolveParams& params, const std::string& filepath) {
    try {
        std::ifstream inFile(filepath);
        if (!inFile.is_open()) {
            std::cerr << "Error: Unable to open file for deserialization: " << filepath << std::endl;
            return false;
        }

        json j;
        inFile >> j;

        params = params.setGlobalScale(j.at("w_global_rescale_param"));
        // the stored weights have the global scale applied

        // Deserialize energy weights
        params.w_smooth_combed = j.at("w_smooth_combed");
        params.w_smooth_combed_precon = j.at("w_smooth_combed_precon");
        params.w_smooth_sym = j.at("w_smooth_sym");
        params.w_smooth_sym_precon = j.at("w_smooth_sym_precon");
        params.w_int_combed = j.at("w_int_combed");
        params.w_int_combed_fixed_weight = j.at("w_int_combed_fixed_weight");
        params.w_int_sym = j.at("w_int_sym");
        params.w_int_sym_fixed_weight = j.at("w_int_sym_fixed_weight");
        params.w_orthog = j.at("w_orthog");
        params.w_orthog_precon = j.at("w_orthog_precon");
        params.w_unit = j.at("w_unit");
        params.w_unit_precon = j.at("w_unit_precon");
        params.w_fit = j.at("w_fit");
        params.w_fit_precon = j.at("w_fit_precon");
        params.w_self_align = j.at("w_self_align");
        params.w_init_scale = j.at("w_init_scale");

        // Deserialize booleans
        params.b_smooth_combed = j.at("b_smooth_combed");
        params.b_smooth_sym = j.at("b_smooth_sym");
        params.b_smooth_asap_combing = j.at("b_smooth_asap_combing");
        params.b_int_combed = j.at("b_int_combed");
        params.b_int_sym = j.at("b_int_sym");
        params.b_orthog = j.at("b_orthog");
        params.b_unit = j.at("b_unit");
        params.b_fit = j.at("b_fit");

        // Deserialize boundary conditions and options
        params.boundary_condition = j.at("boundary_condition");
        params.boundary_hard_constraints = j.at("boundary_hard_constraints");
        params.init_state = j.at("init_state");
        params.connection_for_combing = j.at("connection_for_combing");
        params.fit_type = j.at("fit_type");

        params.w_bound = j.at("w_bound");
        params.w_mint = j.at("w_mint");
        params.w_smooth = j.at("w_smooth");
        params.w_scaled_jacobian = j.at("w_scaled_jacobian");
        params.w_unit_norm = j.at("w_unit_norm");
        params.w_unit_barrier = j.at("w_unit_barrier");
        params.w_fit = j.at("w_fit");

        params.w_outer_step = j.at("w_outer_step");

        params.inner_iter_max_steps = j.at("inner_iter_max_steps");
        params.grad_tol = j.at("grad_tol");
        params.xTol = j.at("xTol");
        params.fTol = j.at("fTol");
        params.reg = j.at("reg");

        params.total_step = j.at("curr_total_step");
        params.outer_step = j.at("curr_outer_step");

        params.lambda_penalty = j.at("lambda_penalty");
        params.b_use_kruskal_tensors_for_sym_smoothness = j.at("b_use_kruskal_tensors_for_sym_smoothness");

        return true;
    } catch (const std::exception& e) {
        std::cerr << "Error: Unable to write data to CSV: " << e.what() << std::endl;
        return false;
    }
}

// Serialize a vector of doubles to a CSV file
bool Serialization::writeCSV(const std::vector<double>& data, const std::string& filepath) {
    try {
        std::ofstream outFile(filepath);
        if (!outFile.is_open()) {
            std::cerr << "Error: Unable to open file for writing CSV: " << filepath << std::endl;
            return false;
        }

        for (const double& value : data) {
            outFile << value << "\n";
        }

        outFile.close();
        return true;
    } catch (const std::exception& e) {
        std::cerr << "Error: Unable to write data to CSV: " << e.what() << std::endl;
        return false;
    }
}

// Deserialize a vector of doubles from a CSV file
bool Serialization::readCSV(std::vector<double>& data, const std::string& filepath) {
    try {
        std::ifstream inFile(filepath);
        if (!inFile.is_open()) {
            std::cerr << "Error: Unable to open file for reading CSV: " << filepath << std::endl;
            return false;
        }

        data.clear();   // Clear the existing data

        double value;
        while (inFile >> value) {
            data.push_back(value);
        }

        inFile.close();
        return true;
    } catch (const std::exception& e) {
        std::cerr << "Error: Unable to read data from CSV: " << e.what() << std::endl;
        return false;
    }
}

// // Serialize JSON data to a JSON file (using nlohmann/json)
// bool Serialization::serializeJSON(const nlohmann::json& jsonData, const std::string& filepath) {
//     try {
//         std::ofstream outFile(filepath);
//         if (!outFile.is_open()) {
//             std::cerr << "Error: Unable to open file for writing JSON: " << filepath << std::endl;
//             return false;
//         }

//         outFile << jsonData.dump(4); // Pretty-print with 4 spaces of indentation
//         outFile.close();
//         return true;
//     } catch (const std::exception& e) {
//         std::cerr << "Error: Unable to serialize JSON data to file: " << e.what() << std::endl;
//         return false;
//     }
// }

// // Deserialize JSON data from a JSON file (using nlohmann/json)
// bool Serialization::deserializeJSON(nlohmann::json& jsonData, const std::string& filepath) {
//     try {
//         std::ifstream inFile(filepath);
//         if (!inFile.is_open()) {
//             std::cerr << "Error: Unable to open file for reading JSON: " << filepath << std::endl;
//             return false;
//         }

//         inFile >> jsonData;
//         inFile.close();
//         return true;
//     } catch (const std::exception& e) {
//         std::cerr << "Error: Unable to deserialize JSON data from file: " << e.what() << std::endl;
//         return false;
//     }
// }
