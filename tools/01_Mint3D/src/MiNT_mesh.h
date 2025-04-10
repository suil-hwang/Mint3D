#ifndef MINT_RANK1_H
#    define MINT_RANK1_H

#    include "Mint3DHook.h"
#    include "PhysicsHook.h"

#    include "Surface.h"

#    include "polyscope/polyscope.h"
#    include "polyscope/surface_mesh.h"

// #    include <TinyAD/ScalarFunction.hh>

// #    include "ADWrapper/ADFuncRunner.h"
// #    include "ADWrapper/ADFunc_TinyAD_Instance.h"

// #include <igl/on_boundary.h>

// #include <igl/writeDMAT.h>

#    ifdef _OPENMP
#        include <omp.h>
#    endif

// template<int DOFS_PER_ELEMENT>
class MiNT_mesh : public Mint3DHook {
   public:
    MiNT_mesh() : Mint3DHook() {
        appState->current_element = Field_View::vec_norms;
        appState->solveType = "mint_frame";
        appState->solveDescription = "MiNT3D solver with 3 and 3 directions per tet dofs per element";

        appState->useBoundaryFrames = true;
    }

    ~MiNT_mesh() {
        // delete _opt;
    }

    virtual void initConfigValues() {}

    virtual void drawGUI() { Mint3DHook::drawGUI(); }

    virtual void initSimulation() {
        // appState->meshName = "triangular_bipyramid";

        // appState->meshName = "tetrahedron_100";
        // appState->meshName = "octahedron_8";

        // Call Parent initialization to load mesh and initialize data structures
        // Add file parsing logic here.
        Mint3DHook::initSimulation();

        // move this inside mint2d
        appState->solveStatus = "init MiNT rank 3";
    }

    // Init state and clone it to the appstate in order to make the visualization accurate.
    void init_opt_state() {}

    virtual void initBoundaryConditions() {}

    virtual void updateRenderGeometry() { Mint3DHook::updateRenderGeometry(); }

    virtual void renderRenderGeometry() { Mint3DHook::renderRenderGeometry(); }

    virtual bool simulateOneStep() {
        std::cout << "simulateOneStep" << std::endl;
        return Mint3DHook::simulateOneStep();
    }

    // This is called after each step.
    virtual void updateAppStateFromOptState() {}

    // this is called after reloading data in order to be able to take forward steps.
    virtual void updateOptStateFromAppState() {}

   protected:
    // Read mesh and compute Tutte embedding
};

#endif   // MINT_RANK1_H

// #undef DOFS_PER_ELEMENT
