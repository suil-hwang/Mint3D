#include "polyscope/polyscope.h"

#include <thread>

#include "VizHook.h"

#include <Eigen/Dense>
#include <Eigen/SparseCholesky>
#include <Eigen/CholmodSupport>
#include <chrono>


static PhysicsHook *hook = NULL;

void toggleSimulation()
{
    if (!hook)
        return;

    if (hook->isPaused())
    {
        hook->run();
    }
    else
        hook->pause();
}

void resetSimulation()
{
    if (!hook)
        return;

    hook->reset();
}

void benchmarkCholmod()
{

        std::cout << "Test linear solver performance" << std::endl;

        std::default_random_engine gen;
        std::uniform_real_distribution<double> dist(0.0,1.0);

        int rows=10000;
        int cols=10000;

        std::vector<Eigen::Triplet<double> > tripletList;
        for(int i=0;i<rows;++i)
        {
            for(int j=0;j<cols;++j)
            {
                auto v_ij=dist(gen);
                auto v_ij_val=dist(gen);     //generate random number
                if(v_ij < 0.01)
                {
                    tripletList.push_back(Eigen::Triplet<double>(i,j,v_ij_val));      //if larger than treshold, insert it
                }
            }
        }
        Eigen::SparseMatrix<double> mat(rows,cols);
        mat.setFromTriplets(tripletList.begin(), tripletList.end()); 
        Eigen::SparseMatrix<double> id(rows,cols);
        id.setIdentity();

        

        Eigen::SparseMatrix<double> test;
        test = mat.transpose(); // + 100. * id;
        test = test * mat; //  + 1000000. * id;
        Eigen::SimplicialLLT<Eigen::SparseMatrix<double>> eigllt; 
        Eigen::CholmodSimplicialLLT<Eigen::SparseMatrix<double>> cholmodllt; 

        auto t1 = std::chrono::high_resolution_clock::now();
        eigllt.compute(test);
        auto t2 = std::chrono::high_resolution_clock::now();
        auto eigllt_ms_int = std::chrono::duration_cast<std::chrono::milliseconds>(t2 - t1);
        std::cout << "eig llt took: " << eigllt_ms_int.count() << "ms\n";
        cholmodllt.compute(test);
        auto t3 = std::chrono::high_resolution_clock::now();
        
        auto cholmodllt_ms_int = std::chrono::duration_cast<std::chrono::milliseconds>(t3 - t2);
        // std::cout << "eig llt took: " << eigllt_ms_int.count() << "ms\n";
        std::cout << "cholmod llt took: " << cholmodllt_ms_int.count() << "ms\n";


        if(eigllt.info()!=Eigen::Success) {
            std::cout << "llt failed" << std::endl;
        }
        
        if(cholmodllt.info()!=Eigen::Success) {
            std::cout << "llt failed" << std::endl;
        }
        
}

void drawGUICallback()
{
	ImGui::PushItemWidth(100); // Make ui elements 100 pixels wide,
							   // instead of full width. Must have 
							   // matching PopItemWidth() below.

    if(hook->showSimButtons() || true)
    {
        if (ImGui::CollapsingHeader("Start Simulation.", ImGuiTreeNodeFlags_DefaultOpen))
        {
            if (ImGui::Button("Run/Pause Sim"))
            {
                toggleSimulation();
            }
            if (ImGui::Button("Reset Sim"))
            {
                resetSimulation();
            }
            if (ImGui::Button("Benchmark Cholmod vs Eigen Solve"))
            {
                benchmarkCholmod();
            }
        }
    }
    hook->drawGUI();
	ImGui::PopItemWidth();
    hook->render();
}


// void callback() {

//   static int numPoints = 2000;
//   static float param = 3.14;

//   ImGui::PushItemWidth(100);

//   // Curvature
//   if (ImGui::Button("add curvature")) {
//  //   addCurvatureScalar();
//   }
  
//   // Normals 
//   if (ImGui::Button("add normals")) {
// //    computeNormals();
//   }

//   // Param
//   if (ImGui::Button("add parameterization")) {
//  //   computeParameterization();
//   }

//   // Geodesics
//   if (ImGui::Button("compute distance")) {
//  //   computeDistanceFrom();
//   }
//   ImGui::SameLine();
// //  ImGui::InputInt("source vertex", &iVertexSource);

//   ImGui::PopItemWidth();
// }

int main(int argc, char **argv) {
  // Configure the argument parser

// Comment out parser.  Should be easy to add in commandline scripting with this though.

  // args::ArgumentParser parser("A simple demo of Polyscope with libIGL.\nBy "
  //                             "Nick Sharp (nsharp@cs.cmu.edu)",
  //                             "");
  // args::Positional<std::string> inFile(parser, "mesh", "input mesh");

  // // Parse args
  // try {
  //   parser.ParseCLI(argc, argv);
  // } catch (args::Help) {
  //   std::cout << parser;
  //   return 0;
  // } catch (args::ParseError e) {
  //   std::cerr << e.what() << std::endl;

  //   std::cerr << parser;
  //   return 1;
  // }


  // Options
  polyscope::options::autocenterStructures = true;
  polyscope::view::windowWidth = 1024;
  polyscope::view::windowHeight = 1024;

  // Initialize polyscope
  polyscope::init();

  hook = new VizHook();
  hook->reset();

  polyscope::state::userCallback = drawGUICallback;

  polyscope::show();

  return 0;
}

