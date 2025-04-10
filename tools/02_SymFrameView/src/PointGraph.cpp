#include "PointGraph.h"

#include <iostream>
#include <unordered_map>
#include <unordered_set>

bool PointGraph::BuildConnectivity() {
  if(pts_.empty() || edges_.empty()) {
    return false;
  }

  int npts = pts_.size();

  int nedges = edges_.size();
  for(int e = 0; e < nedges; e++) {
    int v0 = edges_[e][0];
    int v1 = edges_[e][1];

    if(v0 < 0 || v0 >= npts || v1 < 0 || v1 >= npts) {
      return false;
    }

    pts_[v0].pt_nei_pts_.push_back(v1);
    pts_[v1].pt_nei_pts_.push_back(v0);

    pts_[v0].pt_nei_edges_.push_back(e);
    pts_[v1].pt_nei_edges_.push_back(e);
  }

  for(int i = 0; i < npts; i++) {
    std::sort(pts_[i].pt_nei_pts_.begin(), pts_[i].pt_nei_pts_.end());
    std::sort(pts_[i].pt_nei_edges_.begin(), pts_[i].pt_nei_edges_.end());
  }
  return true;
}

void PointGraph::MergeClosePoints(double eps) {
  for(auto& pt : pts_) {
    pt.eps_ = eps;
  }
  std::unordered_map<GraphNode, int, GraphNodeHash, GraphNodeEqual> node_map;
  std::vector<int> old_2_new_nmap(pts_.size(), -1);
  std::vector<GraphNode> merged_nodes;

  for(int i = 0; i < pts_.size(); i++) {
    if(!node_map.count(pts_[i])) {
      GraphNode new_node;
      new_node.id_ = merged_nodes.size();
      new_node.eps_ = eps;
      new_node.pos_ = pts_[i].pos_;
      merged_nodes.push_back(new_node);

      node_map[pts_[i]] = new_node.id_;
    }
    old_2_new_nmap[i] = node_map[pts_[i]];
  }

  struct pair_hash {
    std::size_t operator () (const std::pair<int, int> &pair) const {
      std::size_t seed = 0;
      auto hash_combine = [&seed](std::size_t h) {
        seed ^= h + 0x9e3779b9 + (seed<<6) + (seed>>2);
      };
      hash_combine(std::hash<int>()(pair.first));
      hash_combine(std::hash<int>()(pair.second));
      return seed;
    }
  };
  std::unordered_set<std::pair<int, int>, pair_hash> edge_set;
  std::vector<std::vector<int>> merged_edges;

  for(int e = 0; e < edges_.size(); e++) {
    int v0 = old_2_new_nmap[edges_[e][0]];
    int v1 = old_2_new_nmap[edges_[e][1]];

    if(v0 == v1) {
      continue;
    }

    if(v1 < v0) {
      std::swap(v0, v1);
    }

    if(edge_set.count(std::pair<int, int>(v0, v1))) {
      continue;
    } else {
      edge_set.insert(std::pair<int, int>(v0, v1));
      merged_edges.push_back({v0, v1});
    }
  }

  pts_ = std::move(merged_nodes);
  edges_ = std::move(merged_edges);
  BuildConnectivity();
}

void PointGraph::MergeConsecutiveParallelSegments(double eps) {
  std::priority_queue<GraphNode, std::vector<GraphNode>, GreaterThan> ns_queue;

  auto compute_node_err = [](GraphNode& node, const Eigen::Vector3d& pre_pos, const Eigen::Vector3d& next_pos) {
    Eigen::Vector3d pre_dir = node.pos_ - pre_pos, cur_dir = next_pos - node.pos_;
    double cos = pre_dir.dot(cur_dir) / pre_dir.norm() / cur_dir.norm();
    node.err_ = std::acos(std::clamp(cos, -1.0, 1.0)); // use clamp to avoid numerical issues
  };

  for(auto node : pts_) {
    if(node.pt_nei_pts_.size() != 2) {
      continue;
    }
    Eigen::Vector3d& pre_pos = pts_[node.pt_nei_pts_[0]].pos_;
    Eigen::Vector3d& next_pos = pts_[node.pt_nei_pts_[1]].pos_;
    node.time_stamp_ = 0;
    compute_node_err(node, pre_pos, next_pos);
    ns_queue.push(node);
  }

  std::vector<bool> valid_vert_flags(pts_.size(), true), valid_edge_flags(edges_.size(), true);
  int num_deleted_edges = 0;
  while(!ns_queue.empty()) {
    auto node = ns_queue.top();
    ns_queue.pop();

    if(!valid_vert_flags[node.id_]) {
      continue;
    }

    // this error is out of date
    if(node.time_stamp_ < pts_[node.id_].time_stamp_) {
      continue;
    }

    if(node.err_ > eps) {
      break;
    }

    // only handle variance-2 nodes
    assert(node.pt_nei_pts_.size() == 2);

    assert(node.time_stamp_ == pts_[node.id_].time_stamp_);

    int nid = pts_[node.id_].id_;
    valid_vert_flags[nid] = false;

    int nid0 = pts_[node.id_].pt_nei_pts_[0];
    int nid1 = pts_[node.id_].pt_nei_pts_[1];

    assert(valid_vert_flags[nid0] && valid_vert_flags[nid1]);

    // replace two edges (nid0, nid) and (nid, nid1) by a single edge (nid0, nid1)
    // we mark the edge (nid0, nid) as invalid, and update the edge (nid, nid1) to (nid0, nid1)
    int invalid_eid = pts_[node.id_].pt_nei_edges_[0];
    int valid_eid = pts_[node.id_].pt_nei_edges_[1];

    if(edges_[valid_eid][0] == nid0 || edges_[valid_eid][1] == nid0) {
      std::swap(invalid_eid, valid_eid);
    }
    valid_edge_flags[invalid_eid] = false;

    for(int i = 0; i < 2; i++) {
      if(edges_[valid_eid][i] == nid) {
        edges_[valid_eid][i] = nid0;
      }
    }

    // update pt neighboring edges
    EraseElementFromUniqueSortedVec(pts_[nid0].pt_nei_edges_, invalid_eid);
    InsertElementToUniqueSortedVec(pts_[nid0].pt_nei_edges_, valid_eid);

    // update pt neighboring points
    EraseElementFromUniqueSortedVec(pts_[nid0].pt_nei_pts_, nid);
    InsertElementToUniqueSortedVec(pts_[nid0].pt_nei_pts_, nid1);

    // update time stamps and err queue
    pts_[nid0].time_stamp_++;
    if(pts_[nid0].pt_nei_pts_.size() == 2) {
      GraphNode candidate_node = pts_[nid0];
      Eigen::Vector3d& pre_pos = pts_[nid1].pos_;
      int next_id = pts_[nid0].pt_nei_pts_[0] != nid1 ? pts_[nid0].pt_nei_pts_[0] : pts_[nid0].pt_nei_pts_[1];
      assert(valid_vert_flags[next_id]);
      Eigen::Vector3d& next_pos = pts_[next_id].pos_;
      compute_node_err(candidate_node, pre_pos, next_pos);
      ns_queue.push(candidate_node);
    }


    EraseElementFromUniqueSortedVec(pts_[nid1].pt_nei_pts_, nid);
    InsertElementToUniqueSortedVec(pts_[nid1].pt_nei_pts_, nid0);

    // update time stamps and err queue
    pts_[nid1].time_stamp_++;
    if(pts_[nid1].pt_nei_pts_.size() == 2) {
      GraphNode candidate_node = pts_[nid1];
      Eigen::Vector3d& pre_pos = pts_[nid0].pos_;
      int next_id = pts_[nid1].pt_nei_pts_[0] != nid0 ? pts_[nid1].pt_nei_pts_[0] : pts_[nid1].pt_nei_pts_[1];
      assert(valid_vert_flags[next_id]);
      Eigen::Vector3d& next_pos = pts_[next_id].pos_;
      compute_node_err(candidate_node, pre_pos, next_pos);
      ns_queue.push(candidate_node);
    }

    num_deleted_edges++;
  }

  std::cout << "num of deleted edges: " << num_deleted_edges << std::endl;
  // reindexing
  ReindexGraph(valid_vert_flags, valid_edge_flags);

}

void PointGraph::ReindexGraph(const std::vector<bool>& valid_node_flags, const std::vector<bool>& valid_edge_flags) {
  PointGraph new_graph;
  std::vector<int> old_2_new_nmap(valid_node_flags.size(), -1);
  for(int i = 0; i < valid_node_flags.size(); i++) {
    if(valid_node_flags[i]) {
      old_2_new_nmap[i] = new_graph.pts_.size();
      GraphNode node;
      node.id_ = new_graph.pts_.size();
      node.pos_ = pts_[i].pos_;
      new_graph.pts_.push_back(node);
    }
  }

  for(int i = 0; i < valid_edge_flags.size(); i++) {
    if(valid_edge_flags[i]) {
      assert(valid_node_flags[edges_[i][0]] && valid_node_flags[edges_[i][1]]);
      std::vector<int> new_edge = {old_2_new_nmap[edges_[i][0]], old_2_new_nmap[edges_[i][1]]};
      new_graph.edges_.push_back(new_edge);
    }
  }

  new_graph.BuildConnectivity();
  *this = std::move(new_graph);
}