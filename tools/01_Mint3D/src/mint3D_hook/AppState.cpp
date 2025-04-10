#include "AppState.h"
#include "FieldView.h"
#include "FileParser.h"
#include "Serialization.h"

#include <Eigen/Sparse>
#include <iostream>

#include <fstream>

// Constructor here
AppState::AppState() { exp_name = "default"; }

void AppState::setSparseMetricFromWeights(Eigen::SparseMatrix<double> &M, const std::vector<double> weights) {}

// Implement the refreshFileLists to populate the lists of files
void AppState::refreshFileLists() { FileParser fileParser(directoryPath); }

bool AppState::LogToFile(const std::string suffix) { return true; }

void AppState::readAllLogFiles() {}
