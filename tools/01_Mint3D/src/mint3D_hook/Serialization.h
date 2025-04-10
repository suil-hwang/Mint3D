#ifndef SERIALIZATION_H
#define SERIALIZATION_H

#include <Eigen/Dense>
#include <string>
#include <vector>
#include "../MiNT3D/MiNTSolveParams.h"

namespace Serialization {

// Eigen Serialization
bool serializeVector(const Eigen::VectorXd& vec, const std::string& filepath);
bool deserializeVector(Eigen::VectorXd& vec, const std::string& filepath);
bool serializeMatrix(const Eigen::MatrixXd& mat, const std::string& filepath, int vector_per_element = 1);
bool deserializeMatrix(Eigen::MatrixXd& mat, const std::string& filepath);

bool deserializeFF3FramesFromMetricDrivenFrames3D(Eigen::MatrixXd& mat, const std::string& filepath);

// MiNT3D Solve Data Serialization
bool serializeSolveParams(MiNT3D::SolveParams& params, const std::string& filepath);
bool deserializeSolveParams(MiNT3D::SolveParams& params, const std::string& filepath);

// CSV Serialization
bool writeCSV(const std::vector<double>& data, const std::string& filepath);
bool readCSV(std::vector<double>& data, const std::string& filepath);

}   // namespace Serialization

#endif   // SERIALIZATION_H