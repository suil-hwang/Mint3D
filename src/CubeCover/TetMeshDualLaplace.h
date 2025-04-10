

#ifndef TETMESHDUALLAPLACE_H
#define TETMESHDUALLAPLACE_H

#include <Eigen/Core>
#include <vector>
#include <Eigen/Dense>
#include <Eigen/Sparse>


/*
 * This takes in a T and V and computes the (default circumcenter) dual laplacian weights
 * specifically  If we think L_k as a matrix, then |nabla L_k|^2 dV can be discretized as:
 * 
 * |(L_k^1 - L_k^0) / dual_eLen|^2 * dual_edge_vol
 * 
 * This is heavily inspired by https://igl.ethz.ch/projects/LB3D/dualLaplace.cpp
 * 
 */

namespace CubeCover {


enum class LapWeightType {
    Combinatorial, 
    Circumcenter,
    Barycenter,
    DualEdgeLength
};



class TetMeshDualLaplace {
public:
    TetMeshDualLaplace(const Eigen::MatrixXi &T, const Eigen::MatrixXd &V);

    void barycenter(const Eigen::Matrix<double, 4, 3>& t, Eigen::Vector3d& c);
    void circumcenter(const Eigen::Matrix<double, 4, 3>& t, Eigen::Vector3d& c);
    double volume(const Eigen::Vector3d& a, const Eigen::Vector3d& b, const Eigen::Vector3d& c, const Eigen::Vector3d& d);

private: 
    std::vector<Eigen::Vector3d> tet_circumcenters; 
    std::vector<double> dual_volumes;
    std::vector<double> dual_edge_lengths;
    std::vector<double> weights; // volume over edge length squared.  

};

}


#endif
