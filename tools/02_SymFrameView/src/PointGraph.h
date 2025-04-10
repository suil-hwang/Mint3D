#pragma once
#include <Eigen/Dense>
#include <vector>
#include <queue>

struct GraphNode {
  GraphNode() = default;
  GraphNode(const Eigen::Vector3d& pos, int id = 1, double err = 0) : pos_(pos), id_(id), err_(err), time_stamp_(0) {}

  Eigen::Vector3d pos_;
  int id_ = -1;
  int time_stamp_ = 0;  // used for graph decimation
  double err_ = 0;    // used for graph decimation
  double eps_ = 1e-6;   // the threshold for same position check

  std::vector<int> pt_nei_edges_;
  std::vector<int> pt_nei_pts_;
};

struct GreaterThan {
  bool operator()(GraphNode& node0, GraphNode& node1)
  {
    return node0.err_ > node1.err_;
  };
};

struct GraphNodeHash {
  void HashCombine(std::size_t& seed, const std::size_t& hash) const {
    seed ^= hash + 0x9e3779b9 + (seed << 6) + (seed >> 2);
  }

  std::size_t operator()(const GraphNode& node) const {
    // Simple hash combining based on position and id
    std::size_t seed = 0;
    // Hashing the position vector; this is a simplified approach
    HashCombine(seed, std::hash<double>()(node.pos_.x()));
    HashCombine(seed, std::hash<double>()(node.pos_.y()));
    HashCombine(seed, std::hash<double>()(node.pos_.z()));

    return seed;
  }
};

struct GraphNodeEqual {
  bool operator()(const GraphNode& node0, const GraphNode& node1) const {
    return (node0.pos_ - node1.pos_).norm() < std::min(node0.eps_, node1.eps_);
  }
};

class PointGraph {
public:
  PointGraph() = default;
  PointGraph(const Eigen::MatrixXd& P, const Eigen::MatrixXi& E) {
    int npts = P.rows();
    pts_.reserve(npts);
    for(int i = 0; i < npts; i++) {
      GraphNode node;
      node.pos_ = P.row(i);
      node.id_ = i;
      pts_.push_back(node);
    }
    int nedges = E.rows();
    assert(E.cols() == 2);
    for(int e = 0; e < nedges; e++) {
      edges_.push_back({E(e, 0), E(e, 1)});
    }
  }
  PointGraph(const std::vector<Eigen::Vector3d>& P, const std::vector<std::vector<int>>& E) : edges_(E) {
    int npts = P.size();
    pts_.reserve(npts);
    for(int i = 0; i < npts; i++) {
      GraphNode node;
      node.pos_ = P[i];
      node.id_ = i;
      pts_.push_back(node);
    }
  }

  std::pair<Eigen::MatrixXd, Eigen::MatrixXi> GetPosEdges() {
    Eigen::MatrixXd P(pts_.size(), 3);
    for(int i = 0; i < pts_.size(); i++) {
      P.row(i) = pts_[i].pos_.transpose();
    }
    Eigen::MatrixXi E(edges_.size(), 2);
    for(int i = 0; i < edges_.size(); i++) {
      E.row(i) << edges_[i][0], edges_[i][1];
    }
    return {P, E};
  }

  bool BuildConnectivity();

  void MergeClosePoints(double eps = 1e-6);
  void MergeConsecutiveParallelSegments(double angle_eps = 0.5 * M_PI / 180);

private:
  bool InsertElementToUniqueSortedVec(std::vector<int>& vec, int elem) {
    auto insert_itr = std::lower_bound(std::begin(vec), std::end(vec), elem);
    if (insert_itr == std::end(vec) || *insert_itr != elem) {
      vec.insert(insert_itr, elem);
      return true;
    }
    return false;
  }

  bool EraseElementFromUniqueSortedVec(std::vector<int>& vec, int elem) {
    auto erase_itr = std::lower_bound(std::begin(vec), std::end(vec), elem);

    if (erase_itr != std::end(vec) && *erase_itr == elem) {
      vec.erase(erase_itr, erase_itr + 1);
      return true;
    }

    return false;
  }

  void ReindexGraph(const std::vector<bool>& valid_node_flags, const std::vector<bool>& valid_edge_flags);

public:
  std::vector<GraphNode> pts_;
  std::vector<std::vector<int>> edges_;
};