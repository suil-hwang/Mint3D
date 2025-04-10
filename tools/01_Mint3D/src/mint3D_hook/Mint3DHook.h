#ifndef MINT3DHOOK_H
#define MINT3DHOOK_H

#include <Eigen/Core>
#include <Eigen/Dense>
#include <string>
#include "AppState.h"   // Include AppState definition

#include "../MiNT3D/MiNTModel.h"
#include "PhysicsHook.h"

#include "FileParser.h"
#include "Surface.h"

// #include <TinyAD/ScalarFunction.hh>
// #include <TinyAD/Utils/LinearSolver.hh>

struct FileTypeInfo {
    std::string folder;
    std::string prefix;
    std::string extension;
    Eigen::MatrixXd &targetMatrix;
};

class Mint3DHook : public virtual PhysicsHook {
   public:
    Mint3DHook() : PhysicsHook() {
        appState = std::make_unique<AppState>();
        opt = std::make_unique<MiNT3D::MiNTModel>();
        //   current_element = Field_View::vec_norms;
    }
    virtual ~Mint3DHook() {}

    virtual void drawGUI();
    virtual void updateRenderGeometry();
    virtual void renderRenderGeometry();
    virtual void initSimulation();
    virtual bool simulateOneStep();
    virtual void pause();
    virtual void updateAppStateFromOptState() { return; };
    virtual void updateOptStateFromAppState() { return; };
    virtual void initConfigValues() { return; };

    // this starts the curl small and takes one step and then increases it
    // by an order of magnitude to find somewhat smooth intialization which conforms to the boundary
    // void warmup();

    void resetAppState();
    void initializeLogFolder();
    void initializeFieldState();
    void solveFromCurrentState();
    void initializeOtherParameters();
    virtual void initBoundaryConditions();
    void initCurlOperators();
    std::string getExperimentShortString();

    void setPoissonBoundaryConditions();
    void setSDFBoundaryConditions();
    void setFramesToRandom(double scale);

    void setFramesToInvertingOnCylinder(bool isInverting, double multiplier, double noise_magnitude);
    Eigen::VectorXd getInvertingFrame(Eigen::Vector2d vec, double z, bool isInverting, double multiplier,
                                      double noise_magnitude);

    // Ok so this is quite tricky, documenting here.

    // This first function will subdivide every tet which has 4 vertices on the boundary
    // This is sufficient to guarentee that the resulting mesh has the property that
    // every tet is adjacent to at most one boundary element.
    // This is necessary to ensure that the problem of solving for normal aligned boundary conditions is well posed.
    void subdivideMeshBoundary();

    // However, it is not sufficient.  In particular, the resulting mesh may still have the property that there are
    // internal facets all of whose vertices are on the boundary.
    // Such facets create issues during the parameterization step, and create artifacts when initializing with poisson
    // boundary conditions.  This function will subdivide all such facets using a different subdivision rule.

    // NOTE: This function must be called *AFTER* subdivideMeshBoundary in order to guarentee that
    //       every tet has at most one such facet.
    void subdivideLockedFacets();

    void findSharpFeatures();

    void launchFieldViewerFromCurrentState();

    bool loadPrimaryData();
    bool loadSecondaryData();

    bool loadGuiState();
    bool loadMetricDrivenFrameField(std::string path);
    bool loadSpecificFrameField(std::string fraFilename);

    void takeFrameFieldDual();

    std::unique_ptr<AppState> appState;   // Pointer to AppState instance

    std::unique_ptr<OutputState> outputData;   // threadsafe copy of the output fields.

    // Integrate the new solver
    std::unique_ptr<MiNT3D::MiNTModel> opt;

    std::unique_ptr<FileParser> fileParser;

   private:
    void updateOptimizationParameters();
    void checkAndUpdateConvergence(double decrement, double energy);
    void updateVisualizationData();
    void finalizeIteration();
    void recomputeSingularStructures();
    void updateFrameVisualization();
    void updateSingularityVisualiation();
    Eigen::VectorXd solvePoissonOnTetMesh(const Eigen::MatrixXd &V, const Eigen::MatrixXi &T);

    // Other private member variables and functions as needed
};

#endif   // MINT3DHOOK_H
