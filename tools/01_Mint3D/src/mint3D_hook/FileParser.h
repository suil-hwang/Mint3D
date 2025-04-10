#ifndef FILE_PARSER_H
#define FILE_PARSER_H

#include <Eigen/Dense>
#include <optional>
#include <string>
#include <vector>
#include "Serialization.h"   // Include the Serialization header

#include "AppState.h"

// Define the FileType enum
enum class FileType {
    BFRA,
    BMOM,
    OBJ,
    FRA,
    MOM,
};

/**
 * The FileParser class is responsible for parsing files of different types
 * within a given directory and retrieving relevant data.
 */
class FileParser {
   public:
    // the expected directory structure is:
    // path/{name}.mesh
    // path/inner_iters
    // path/outer_iters
    explicit FileParser(const std::string& directoryPath, bool load_inner = true);

    std::string getFileWithID(const std::string& folder, const std::string& prefix, const std::string& extension,
                              int fileId);

    /*

        // Parse a file with a given ID
        bool parseFileWithID(Eigen::VectorXd& data, FileType fileType, int fileId);

        // Retrieve the largest file's data
        bool parseLargestFile(Eigen::VectorXd& data, FileType fileType);

        void parseToAppState(AppState* appState, int fileId);

        // // Load an .obj file if found
        // bool parseObjFile(Eigen::MatrixXd& V, Eigen::MatrixXi& F);

        // Update the current directory path
        void setDirectoryPath(const std::string& directoryPath);

        */

    std::string objFilePath;
    std::string meshFilePath;
    // int numIDs = 0;

    int minID = 0;
    int maxID = 0;

   private:
    std::string directoryPath;
    std::vector<std::string> conf_files;
    // std::vector<std::string> bmomFiles;

    // std::string largestBfraFile;
    // std::string largestBmomFile;

    // Helper functions to scan and sort files
    void scanDirectory();
    void findFileBounds();
};

#endif   // FILE_PARSER_H
