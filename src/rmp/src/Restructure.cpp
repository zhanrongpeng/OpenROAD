// SPDX-License-Identifier: BSD-3-Clause
// Copyright (c) 2019-2025, The OpenROAD Authors

#include "rmp/Restructure.h"

#include <sys/types.h>
#include <sys/wait.h>
#include <unistd.h>

#include <cstring>
#include <ctime>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <limits>
#include <set>
#include <sstream>
#include <string>
#include <system_error>
#include <utility>
#include <vector>

#include "annealing_strategy.h"
#include "base/abc/abc.h"
#include "base/main/abcapis.h"
#include "cut/abc_init.h"
#include "cut/abc_library_factory.h"
#include "cut/blif.h"
#include "db_sta/dbNetwork.hh"
#include "db_sta/dbSta.hh"
#include "odb/db.h"
#include "odb/dbTransform.h"
#include "odb/db.h"
#include "rsz/Resizer.hh"
#include "sta/Delay.hh"
#include "sta/Graph.hh"
#include "sta/Liberty.hh"
#include "sta/Network.hh"
#include "sta/NetworkClass.hh"
#include "sta/Path.hh"
#include "sta/PathEnd.hh"
#include "sta/PathExpanded.hh"
#include "sta/PatternMatch.hh"
#include "sta/PortDirection.hh"
#include "sta/Sdc.hh"
#include "sta/Search.hh"
#include "sta/Sta.hh"
#include "utl/Logger.h"
#include "zero_slack_strategy.h"

namespace rmp {

using abc::Abc_Frame_t;
using abc::Abc_FrameGetGlobalFrame;
using abc::Abc_Start;
using abc::Abc_Stop;
using cut::Blif;
using utl::RMP;

Restructure::Restructure(utl::Logger* logger,
                         sta::dbSta* open_sta,
                         odb::dbDatabase* db,
                         rsz::Resizer* resizer,
                         est::EstimateParasitics* estimate_parasitics)
{
  logger_ = logger;
  db_ = db;
  open_sta_ = open_sta;
  resizer_ = resizer;
  estimate_parasitics_ = estimate_parasitics;

  cut::abcInit();
}

void Restructure::deleteComponents()
{
}

Restructure::~Restructure()
{
  deleteComponents();
}

void Restructure::reset()
{
  lib_file_names_.clear();
  path_insts_.clear();
}

void Restructure::resynth(sta::Corner* corner)
{
  ZeroSlackStrategy zero_slack_strategy(corner);
  zero_slack_strategy.OptimizeDesign(
      open_sta_, name_generator_, resizer_, logger_);
}

void Restructure::resynthAnnealing(sta::Corner* corner)
{
  AnnealingStrategy annealing_strategy(corner,
                                       slack_threshold_,
                                       annealing_seed_,
                                       annealing_temp_,
                                       annealing_iters_,
                                       annealing_revert_after_,
                                       annealing_init_ops_);
  annealing_strategy.OptimizeDesign(
      open_sta_, name_generator_, resizer_, logger_);
}

void Restructure::run(char* liberty_file_name,
                      float slack_threshold,
                      unsigned max_depth,
                      char* workdir_name,
                      char* abc_logfile)
{
  reset();
  block_ = db_->getChip()->getBlock();
  if (!block_) {
    return;
  }

  logfile_ = abc_logfile;
  sta::Slack worst_slack = slack_threshold;

  lib_file_names_.emplace_back(liberty_file_name);
  work_dir_name_ = workdir_name;
  work_dir_name_ = work_dir_name_ + "/";

  if (is_area_mode_) {  // Only in area mode
    removeConstCells();
  }

  getBlob(max_depth);

  if (!path_insts_.empty()) {
    runABC();

    postABC(worst_slack);
  }
}

void Restructure::selectConeByPathDelay(const ConeSelectionConfig& config) {
  open_sta_->ensureGraph();
  open_sta_->ensureLevelized();
  open_sta_->searchPreamble();
  
  auto sta_state = open_sta_->search();
  sta::VertexSet* end_points = sta_state->endpoints();
  
  logger_->report("Starting smart cone selection from {} endpoints", 
                  end_points->size());
  
  std::set<sta::Vertex*> all_selected_vertices;
  
  // 找到 worst slack 的 endpoint
  sta::Vertex* worst_endpoint = nullptr;
  sta::Slack worst_slack = std::numeric_limits<float>::max();
  
  for (auto& end_point : *end_points) {
    sta::Slack slack = open_sta_->vertexSlack(end_point, sta::MinMax::max());
    if (slack < worst_slack) {
      worst_slack = slack;
      worst_endpoint = end_point;
    }
  }
  
  if (!worst_endpoint) {
    logger_->warn(RMP, 12, "No endpoint found for cone selection");
    return;
  }
  
  logger_->report("Selecting worst endpoint with slack: {}", worst_slack);
  logger_->report("Building path tree for endpoint: {}", 
                  open_sta_->getDbNetwork()->pathName(worst_endpoint->pin()));
  
  // 构建路径树
  PathNode* root = buildPathTree(worst_endpoint, 0, config);
  
  if (root) {
    // 从树中选择节点
    std::set<sta::Vertex*> selected_vertices;
    int current_size = 0;
    selectNodesFromTree(root, config, selected_vertices, current_size);
    
    logger_->report("Selected {} vertices from this endpoint", 
                    selected_vertices.size());
    
    // 合并到总的选择集合
    all_selected_vertices.insert(selected_vertices.begin(), 
                                 selected_vertices.end());
    
    // 清理树
    cleanupPathTree(root);
  }
  
  // 收集选中的instances
  collectSelectedInstances(all_selected_vertices);
  
  logger_->report("Total selected {} instances for restructuring", 
                  path_insts_.size());
}

// 获取instance的物理坐标
std::pair<float, float> Restructure::getInstanceCoordinate(odb::dbInst* inst) {
  if (!inst) {
    return {0.0f, 0.0f};
  }
  
  const odb::dbTransform transform = inst->getTransform();
  const odb::Point origin = transform.getOffset();
  return {static_cast<float>(origin.getX()), static_cast<float>(origin.getY())};
}

Restructure::PathNode* Restructure::buildPathTree(
    sta::Vertex* vertex,
    int current_depth,
    const ConeSelectionConfig& config) {
  
  if (!vertex) {
    debugPrint(logger_, RMP, "remap", 1, "Null vertex at depth {}", current_depth);
    return nullptr;
  }
  
  if (current_depth >= config.max_cone_depth) {
    debugPrint(logger_, RMP, "remap", 1, "Reached max depth {}", current_depth);
    return nullptr;
  }
  
  PathNode* node = new PathNode();
  node->vertex = vertex;
  node->pin = vertex->pin();
  node->arrival_time = getVertexArrivalTime(vertex);
  node->slack = getVertexSlack(vertex);
  node->depth = current_depth;
  
  // 获取当前节点的物理坐标
  odb::dbITerm* term = nullptr;
  odb::dbBTerm* port = nullptr;
  odb::dbModITerm* moditerm = nullptr;
  open_sta_->getDbNetwork()->staToDb(vertex->pin(), term, port, moditerm);
  
  if (term) {
    odb::dbInst* inst = term->getInst();
    if (inst) {
      auto coords = getInstanceCoordinate(inst);
      node->x = coords.first;
      node->y = coords.second;
    }
  }
  
  debugPrint(logger_, RMP, "remap", 1, 
             "Building node at depth {}: {} (arrival: {}, slack: {}, coord: {:.2f}, {:.2f})",
             current_depth,
             open_sta_->getDbNetwork()->pathName(vertex->pin()),
             node->arrival_time,
             node->slack,
             node->x,
             node->y);
  
  // 获取fanin edges
  sta::VertexInEdgeIterator edge_iter(vertex, open_sta_->graph());
  std::vector<sta::Vertex*> fanins;
  int skipped_ports = 0;
  int skipped_seq = 0;
  
  while (edge_iter.hasNext()) {
    sta::Edge* edge = edge_iter.next();
    sta::Vertex* from_vertex = edge->from(open_sta_->graph());
    
    // 跳过端口和sequential元素
    if (open_sta_->getDbNetwork()->isTopLevelPort(from_vertex->pin())) {
      skipped_ports++;
      continue;
    }
    
    sta::LibertyCell* cell = open_sta_->getDbNetwork()->libertyCell(
        open_sta_->getDbNetwork()->instance(from_vertex->pin()));
    if (cell && cell->hasSequentials()) {
      skipped_seq++;
      continue;
    }
    
    fanins.push_back(from_vertex);
  }
  
  debugPrint(logger_, RMP, "remap", 1,
             "Found {} fanins (skipped {} ports, {} sequential)",
             fanins.size(), skipped_ports, skipped_seq);
  
  // 递归构建子节点(最多取2个最差的fanin)
  if (!fanins.empty()) {
    // 按arrival time排序,取延迟最大的
    std::sort(fanins.begin(), fanins.end(), 
              [this](sta::Vertex* a, sta::Vertex* b) {
                return getVertexArrivalTime(a) > getVertexArrivalTime(b);
              });
    
    debugPrint(logger_, RMP, "remap", 1,
               "Recursing into {} fanin(s)", std::min(2, (int)fanins.size()));
    
    // 构建左子节点(最差的fanin)
    node->left_child = buildPathTree(fanins[0], current_depth + 1, config);
    
    // 如果有第二个fanin,构建右子节点
    if (fanins.size() > 1) {
      node->right_child = buildPathTree(fanins[1], current_depth + 1, config);
    }
  } else {
    debugPrint(logger_, RMP, "remap", 1, "No valid fanins, stopping here");
  }
  
  return node;
}

float Restructure::calculateDelayGap(PathNode* left, PathNode* right) {
  if (!left || !right) {
    return std::numeric_limits<float>::max();
  }
  
  // 计算两个子节点的arrival time差异
  float gap = std::abs(left->arrival_time - right->arrival_time);
  
  // 归一化:相对于较大的arrival time
  float max_arrival = std::max(left->arrival_time, right->arrival_time);
  if (max_arrival > 0) {
    return gap / max_arrival;
  }
  
  return gap;
}

// 计算两个节点之间的物理距离（HPWL）
float Restructure::calculatePhysicalDistance(PathNode* node1, PathNode* node2) {
  if (!node1 || !node2) {
    return std::numeric_limits<float>::max();
  }
  
  // 使用曼哈顿距离作为HPWL的近似
  float dx = std::abs(node1->x - node2->x);
  float dy = std::abs(node1->y - node2->y);
  float distance = dx + dy;
  
  return distance;
}

bool Restructure::shouldSelectChild(PathNode* parent,
                                     PathNode* left,
                                     PathNode* right,
                                     const ConeSelectionConfig& config) {
  if (!parent) {
    return false;
  }
  
  // 如果只有一个子节点,直接选择
  if (!left || !right) {
    return true;
  }
  
  // 计算延迟差异比例
  float delay_gap = calculateDelayGap(left, right);
  
  // 如果差异大于阈值,需要进一步判断
  if (delay_gap > config.delay_threshold_ratio) {
    // 选择arrival time较大的（较差的）子节点
    PathNode* worse_child = (left->arrival_time > right->arrival_time) ? left : right;
    
    // 检查物理距离：如果太远就不选了
    float distance = calculatePhysicalDistance(parent, worse_child);
    debugPrint(logger_, RMP, "remap", 1,
               "Delay gap {:.2%} > threshold {:.2%}, checking distance: {:.2f} um vs threshold {:.2f} um",
               delay_gap, config.delay_threshold_ratio, 
               distance, config.distance_threshold_um);
    
    if (distance > config.distance_threshold_um) {
      debugPrint(logger_, RMP, "remap", 1,
                 "Child too far ({:.2f} um > {:.2f} um), not selecting",
                 distance, config.distance_threshold_um);
      return false;
    }
    
    // 距离在阈值内，选择较差的子节点
    return false;  // 返回false表示只选一个（worse_child）
  }
  
  // 差异小,检查物理距离后两个都选
  // 如果两个子节点都离父节点太远，可能需要谨慎选择
  float left_dist = calculatePhysicalDistance(parent, left);
  float right_dist = calculatePhysicalDistance(parent, right);
  
  debugPrint(logger_, RMP, "remap", 1,
             "Delay gap {:.2%} <= threshold {:.2%}, distances: left={:.2f} um, right={:.2f} um",
             delay_gap, config.delay_threshold_ratio, left_dist, right_dist);
  
  return true;
}

void Restructure::selectNodesFromTree(
    PathNode* node,
    const ConeSelectionConfig& config,
    std::set<sta::Vertex*>& selected_vertices,
    int& current_size) {
  
  if (!node || current_size >= config.max_cone_size) {
    debugPrint(logger_, RMP, "remap", 1,
               "Stopping: node={}, size={}/{}", 
               (void*)node, current_size, config.max_cone_size);
    return;
  }
  
  // 选择当前节点
  node->selected = true;
  selected_vertices.insert(node->vertex);
  current_size++;
  
  debugPrint(logger_, RMP, "remap", 1,
             "Selected vertex {} (total: {})",
             open_sta_->getDbNetwork()->pathName(node->pin),
             current_size);
  
  // 检查是否有子节点
  if (!node->left_child && !node->right_child) {
    debugPrint(logger_, RMP, "remap", 1, "Leaf node, stopping");
    return;  // 叶子节点,停止
  }
  
  // 判断子节点的延迟改善潜力
  float parent_slack = std::abs(node->slack);
  
  // 注释掉 slack 阈值检查,让它继续选择
  // if (parent_slack < config.min_improvement_threshold) {
  //   debugPrint(logger_, RMP, "remap", 2,
  //              "Stopping at node with slack {}, below threshold {}",
  //              parent_slack, config.min_improvement_threshold);
  //   return;
  // }
  
  // 判断是否应该选择两个子节点
  bool select_both = shouldSelectChild(node, node->left_child, 
                                       node->right_child, config);
  
  debugPrint(logger_, RMP, "remap", 1,
             "Has children: left={}, right={}, select_both={}",
             (void*)node->left_child, (void*)node->right_child, select_both);
  
  if (select_both) {
    // 两个子节点延迟差异不大,都选择
    if (node->left_child) {
      debugPrint(logger_, RMP, "remap", 2,
                 "Selecting both children, delay gap is small");
      selectNodesFromTree(node->left_child, config, 
                         selected_vertices, current_size);
    }
    if (node->right_child) {
      selectNodesFromTree(node->right_child, config, 
                         selected_vertices, current_size);
    }
  } else {
    // 差异大,只选择较差的那个继续递归
    PathNode* worse_child = nullptr;
    PathNode* better_child = nullptr;
    
    if (node->left_child && node->right_child) {
      if (node->left_child->arrival_time > node->right_child->arrival_time) {
        worse_child = node->left_child;
        better_child = node->right_child;
      } else {
        worse_child = node->right_child;
        better_child = node->left_child;
      }
      
      debugPrint(logger_, RMP, "remap", 2,
                 "Large delay gap detected, selecting worse path only. "
                 "Gap: {:.2%}", 
                 calculateDelayGap(node->left_child, node->right_child));
      
      // 较差的子节点:继续递归
      selectNodesFromTree(worse_child, config, selected_vertices, current_size);
      
      // 较好的子节点:标记为边界,选中但不递归
      better_child->selected = true;
      better_child->is_boundary = true;
      selected_vertices.insert(better_child->vertex);
      current_size++;
      
      debugPrint(logger_, RMP, "remap", 2,
                 "Better child marked as boundary, not recursing");
    } else if (node->left_child) {
      // 只有左子节点
      selectNodesFromTree(node->left_child, config, selected_vertices, current_size);
    } else if (node->right_child) {
      // 只有右子节点
      selectNodesFromTree(node->right_child, config, selected_vertices, current_size);
    }
  }
}

void Restructure::collectSelectedInstances(
    const std::set<sta::Vertex*>& vertices) {
  
  path_insts_.clear();
  
  for (sta::Vertex* vertex : vertices) {
    sta::Pin* pin = vertex->pin();
    
    odb::dbITerm* term = nullptr;
    odb::dbBTerm* port = nullptr;
    odb::dbModITerm* moditerm = nullptr;
    open_sta_->getDbNetwork()->staToDb(pin, term, port, moditerm);
    
    if (term) {
      odb::dbInst* inst = term->getInst();
      if (inst && !inst->getMaster()->isBlock()) {
        path_insts_.insert(inst);
      }
    }
  }
}

void Restructure::cleanupPathTree(PathNode* node) {
  if (!node) {
    return;
  }
  
  cleanupPathTree(node->left_child);
  cleanupPathTree(node->right_child);
  delete node;
}

float Restructure::getVertexArrivalTime(sta::Vertex* vertex) {
  sta::Arrival arrival = open_sta_->vertexArrival(vertex, sta::MinMax::max());
  return sta::delayAsFloat(arrival);
}

float Restructure::getVertexSlack(sta::Vertex* vertex) {
  sta::Slack slack = open_sta_->vertexSlack(vertex, sta::MinMax::max());
  return sta::delayAsFloat(slack);
}

void Restructure::getBlob(unsigned max_depth)
{
  open_sta_->ensureGraph();
  open_sta_->ensureLevelized();
  open_sta_->searchPreamble();

  sta::PinSet ends(open_sta_->getDbNetwork());

  getEndPoints(ends, is_area_mode_, max_depth);
  
  if (!ends.empty()) {
    if (is_area_mode_) {
      sta::PinSet boundary_points = resizer_->findFaninFanouts(ends);
      logger_->report("Found {} pins in extracted logic.",
                      boundary_points.size());
      
      for (const sta::Pin* pin : boundary_points) {
        odb::dbITerm* term = nullptr;
        odb::dbBTerm* port = nullptr;
        odb::dbModITerm* moditerm = nullptr;
        open_sta_->getDbNetwork()->staToDb(pin, term, port, moditerm);
        
        if (term && !term->getInst()->getMaster()->isBlock()) {
          path_insts_.insert(term->getInst());
        }
      }
    } else {
      ConeSelectionConfig config;
      config.max_cone_depth = max_depth;
      config.delay_threshold_ratio = 0.3;     
      config.max_cone_size = 500;              
      config.min_improvement_threshold = 0.05;
      config.distance_threshold_um = 50.0;  // 50 um threshold for physical distance

      selectConeByPathDelay(config);
    }

    logger_->report("Found {} instances for restructuring.",
                    path_insts_.size());
  }
}

void Restructure::runABC()
{
  const std::string prefix
      = work_dir_name_ + std::string(block_->getConstName());
  input_blif_file_name_ = prefix + "_crit_path.blif";
  coord_file_name_ = prefix + "_coords.txt";
  std::vector<std::string> files_to_remove;

  debugPrint(logger_,
             utl::RMP,
             "remap",
             1,
             "Constants before remap {}",
             countConsts(block_));

  Blif blif_(
      logger_, open_sta_, locell_, loport_, hicell_, hiport_, ++blif_call_id_);
  blif_.setReplaceableInstances(path_insts_);
  blif_.writeBlif(input_blif_file_name_.c_str(), !is_area_mode_);
  debugPrint(
      logger_, RMP, "remap", 1, "Writing blif file {}", input_blif_file_name_);
  files_to_remove.emplace_back(input_blif_file_name_);

  // Write instance coordinates for ABC
  writeInstanceCoordinates(coord_file_name_);
  files_to_remove.emplace_back(coord_file_name_);

  // abc optimization
  std::vector<Mode> modes;
  std::vector<pid_t> child_proc;

  if (is_area_mode_) {
    // Area Mode
    modes = {Mode::AREA_1, Mode::AREA_2, Mode::AREA_3};
  } else {
    // Delay Mode
    modes = {Mode::DELAY_1, Mode::DELAY_2, Mode::DELAY_3, Mode::DELAY_4};
  }

  child_proc.resize(modes.size(), 0);

  std::string best_blif;
  int best_inst_count = std::numeric_limits<int>::max();
  float best_delay_gain = std::numeric_limits<float>::max();

  debugPrint(
      logger_, RMP, "remap", 1, "Running ABC with {} modes.", modes.size());

  for (size_t curr_mode_idx = 0; curr_mode_idx < modes.size();
       curr_mode_idx++) {
    output_blif_file_name_
        = prefix + std::to_string(curr_mode_idx) + "_crit_path_out.blif";

    opt_mode_ = modes[curr_mode_idx];

    const std::string abc_script_file
        = prefix + std::to_string(curr_mode_idx) + "ord_abc_script.tcl";
    if (logfile_.empty()) {
      logfile_ = prefix + "abc.log";
    }

    debugPrint(logger_,
               RMP,
               "remap",
               1,
               "Writing ABC script file {}.",
               abc_script_file);

    if (writeAbcScript(abc_script_file)) {
      // call linked abc
      Abc_Start();
      Abc_Frame_t* abc_frame = Abc_FrameGetGlobalFrame();
      const std::string command = "source " + abc_script_file;
      child_proc[curr_mode_idx]
          = Cmd_CommandExecute(abc_frame, command.c_str());
      if (child_proc[curr_mode_idx]) {
        logger_->error(RMP, 6, "Error executing ABC command {}.", command);
        return;
      }
      Abc_Stop();
      // exit linked abc
      files_to_remove.emplace_back(abc_script_file);
    }
  }  // end modes

  // Inspect ABC results to choose blif with least instance count
  for (int curr_mode_idx = 0; curr_mode_idx < modes.size(); curr_mode_idx++) {
    // Skip failed ABC runs
    if (child_proc[curr_mode_idx] != 0) {
      continue;
    }

    output_blif_file_name_
        = prefix + std::to_string(curr_mode_idx) + "_crit_path_out.blif";
    const std::string abc_log_name = logfile_ + std::to_string(curr_mode_idx);

    int level_gain = 0;
    float delay = std::numeric_limits<float>::max();
    int num_instances = 0;
    bool success = readAbcLog(abc_log_name, level_gain, delay);
    if (success) {
      success
          = blif_.inspectBlif(output_blif_file_name_.c_str(), num_instances);
      logger_->report(
          "Optimized to {} instances in iteration {} with max path depth "
          "decrease of {}, delay of {}.",
          num_instances,
          curr_mode_idx,
          level_gain,
          delay);

      if (success) {
        if (is_area_mode_) {
          if (num_instances < best_inst_count) {
            best_inst_count = num_instances;
            best_blif = output_blif_file_name_;
          }
        } else {
          // Using only DELAY_4 for delay based gain since other modes not
          // showing good gains
          if (modes[curr_mode_idx] == Mode::DELAY_4) {
            best_delay_gain = delay;
            best_blif = output_blif_file_name_;
          }
        }
      }
    }
    files_to_remove.emplace_back(output_blif_file_name_);
  }

  if (best_inst_count < std::numeric_limits<int>::max()
      || best_delay_gain < std::numeric_limits<float>::max()) {
    // read back netlist
    debugPrint(logger_, RMP, "remap", 1, "Reading blif file {}.", best_blif);
    blif_.readBlif(best_blif.c_str(), block_);
    debugPrint(logger_,
               utl::RMP,
               "remap",
               1,
               "Number constants after restructure {}.",
               countConsts(block_));
  } else {
    logger_->info(
        RMP, 13, "All re-synthesis runs discarded, keeping original netlist.");
  }

  for (const auto& file_to_remove : files_to_remove) {
    if (!logger_->debugCheck(RMP, "remap", 1)) {
      std::error_code err;
      if (std::filesystem::remove(file_to_remove, err); err) {
        logger_->error(RMP, 11, "Fail to remove file {}", file_to_remove);
      }
    }
  }
}

void Restructure::postABC(float worst_slack)
{
  // Leave the parasitics up to date.
  estimate_parasitics_->estimateWireParasitics();
}
void Restructure::getEndPoints(sta::PinSet& ends,
                               bool area_mode,
                               unsigned max_depth)
{
  auto sta_state = open_sta_->search();
  sta::VertexSet* end_points = sta_state->endpoints();
  std::size_t path_found = end_points->size();
  logger_->report("Number of paths for restructure are {}", path_found);
  for (auto& end_point : *end_points) {
    if (!is_area_mode_) {
      sta::Path* path
          = open_sta_->vertexWorstSlackPath(end_point, sta::MinMax::max());
      sta::PathExpanded expanded(path, open_sta_);
      // Members in expanded include gate output and net so divide by 2
      logger_->report("Found path of depth {}", expanded.size() / 2);
      if (expanded.size() / 2 > max_depth) {
        ends.insert(end_point->pin());
        // Use only one end point to limit blob size for timing
        break;
      }
    } else {
      ends.insert(end_point->pin());
    }
  }

  // unconstrained end points
  if (is_area_mode_) {
    auto errors = open_sta_->checkTiming(false /*no_input_delay*/,
                                         false /*no_output_delay*/,
                                         false /*reg_multiple_clks*/,
                                         true /*reg_no_clks*/,
                                         true /*unconstrained_endpoints*/,
                                         false /*loops*/,
                                         false /*generated_clks*/);
    debugPrint(logger_, RMP, "remap", 1, "Size of errors = {}", errors.size());
    if (!errors.empty() && errors[0]->size() > 1) {
      sta::CheckError* error = errors[0];
      bool first = true;
      for (auto pinName : *error) {
        debugPrint(logger_, RMP, "remap", 1, "Unconstrained pin: {}", pinName);
        if (!first && open_sta_->getDbNetwork()->findPin(pinName)) {
          ends.insert(open_sta_->getDbNetwork()->findPin(pinName));
        }
        first = false;
      }
    }
    if (errors.size() > 1 && errors[1]->size() > 1) {
      sta::CheckError* error = errors[1];
      bool first = true;
      for (auto pinName : *error) {
        debugPrint(logger_, RMP, "remap", 1, "Unclocked pin: {}", pinName);
        if (!first && open_sta_->getDbNetwork()->findPin(pinName)) {
          ends.insert(open_sta_->getDbNetwork()->findPin(pinName));
        }
        first = false;
      }
    }
  }
  logger_->report("Found {} end points for restructure", ends.size());
}

int Restructure::countConsts(odb::dbBlock* top_block)
{
  int const_nets = 0;
  for (auto block_net : top_block->getNets()) {
    if (block_net->getSigType().isSupply()) {
      const_nets++;
    }
  }

  return const_nets;
}

void Restructure::removeConstCells()
{
  if (hicell_.empty() || locell_.empty()) {
    return;
  }

  odb::dbMaster* hicell_master = nullptr;
  odb::dbMTerm* hiterm = nullptr;
  odb::dbMaster* locell_master = nullptr;
  odb::dbMTerm* loterm = nullptr;

  for (auto&& lib : block_->getDb()->getLibs()) {
    hicell_master = lib->findMaster(hicell_.c_str());

    locell_master = lib->findMaster(locell_.c_str());
    if (locell_master && hicell_master) {
      break;
    }
  }
  if (!hicell_master || !locell_master) {
    return;
  }

  hiterm = hicell_master->findMTerm(hiport_.c_str());
  loterm = locell_master->findMTerm(loport_.c_str());
  if (!hiterm || !loterm) {
    return;
  }

  open_sta_->clearLogicConstants();
  open_sta_->findLogicConstants();
  std::set<odb::dbInst*> constInsts;
  int const_cnt = 1;
  for (auto inst : block_->getInsts()) {
    int outputs = 0;
    int const_outputs = 0;
    auto master = inst->getMaster();
    sta::LibertyCell* cell = open_sta_->getDbNetwork()->libertyCell(
        open_sta_->getDbNetwork()->dbToSta(master));
    if (cell->hasSequentials()) {
      continue;
    }

    for (auto&& iterm : inst->getITerms()) {
      if (iterm->getSigType() == odb::dbSigType::POWER
          || iterm->getSigType() == odb::dbSigType::GROUND) {
        continue;
      }

      if (iterm->getIoType() != odb::dbIoType::OUTPUT) {
        continue;
      }
      outputs++;
      auto pin = open_sta_->getDbNetwork()->dbToSta(iterm);
      sta::LogicValue pinVal = open_sta_->simLogicValue(pin);
      if (pinVal == sta::LogicValue::one || pinVal == sta::LogicValue::zero) {
        odb::dbNet* net = iterm->getNet();
        if (net) {
          odb::dbMaster* const_master = (pinVal == sta::LogicValue::one)
                                            ? hicell_master
                                            : locell_master;
          odb::dbMTerm* const_port
              = (pinVal == sta::LogicValue::one) ? hiterm : loterm;
          std::string inst_name = "rmp_const_" + std::to_string(const_cnt);
          debugPrint(logger_,
                     RMP,
                     "remap",
                     2,
                     "Adding cell {} inst {} for {}",
                     const_master->getName(),
                     inst_name,
                     inst->getName());
          auto new_inst
              = odb::dbInst::create(block_, const_master, inst_name.c_str());
          if (new_inst) {
            iterm->disconnect();
            new_inst->getITerm(const_port)->connect(net);
          } else {
            logger_->warn(RMP, 9, "Could not create instance {}.", inst_name);
          }
        }
        const_outputs++;
        const_cnt++;
      }
    }
    if (outputs > 0 && outputs == const_outputs) {
      constInsts.insert(inst);
    }
  }
  open_sta_->clearLogicConstants();

  debugPrint(
      logger_, RMP, "remap", 2, "Removing {} instances...", constInsts.size());

  for (auto inst : constInsts) {
    removeConstCell(inst);
  }
  logger_->report("Removed {} instances with constant outputs.",
                  constInsts.size());
}

void Restructure::removeConstCell(odb::dbInst* inst)
{
  for (auto iterm : inst->getITerms()) {
    iterm->disconnect();
  }
  odb::dbInst::destroy(inst);
}

bool Restructure::writeAbcScript(const std::string& file_name)
{
  std::ofstream script(file_name.c_str());

  if (!script.is_open()) {
    logger_->error(RMP, 3, "Cannot open file {} for writing.", file_name);
    return false;
  }

  for (const auto& lib_name : lib_file_names_) {
    // abc read_lib prints verbose by default, -v toggles to off to avoid read
    // time being printed
    std::string read_lib_str = "read_lib -v " + lib_name + "\n";
    script << read_lib_str;
  }

  script << "read_blif -n " << input_blif_file_name_ << '\n';

  if (logger_->debugCheck(RMP, "remap", 1)) {
    script << "write_verilog " << input_blif_file_name_ + std::string(".v")
           << '\n';
  }

  // Physical-aware mapping: skip rewrite/refactor, use strash + read_coords + if
  script << "strash\n";
  script << "read_coords " << coord_file_name_ << '\n';

  // Use if command with wire-aware mapping (-W flag)
  // WireDelayFactor 0.001 = 1nm per unit wirelength (scaled)
  script << "if -W 0.001\n";

  script << "write_blif " << output_blif_file_name_ << '\n';

  if (logger_->debugCheck(RMP, "remap", 1)) {
    script << "write_verilog " << output_blif_file_name_ + std::string(".v")
           << '\n';
  }

  script.close();

  return true;
}

void Restructure::writeOptCommands(std::ofstream& script)
{
  std::string choice
      = "alias choice \"fraig_store; resyn2; fraig_store; resyn2; fraig_store; "
        "fraig_restore\"";
  std::string choice2
      = "alias choice2 \"fraig_store; balance; fraig_store; resyn2; "
        "fraig_store; resyn2; fraig_store; resyn2; fraig_store; "
        "fraig_restore\"";
  script << "bdd; sop\n";

  script << "alias resyn2 \"balance; rewrite; refactor; balance; rewrite; "
            "rewrite -z; balance; refactor -z; rewrite -z; balance\""
         << '\n';
  script << choice << '\n';
  script << choice2 << '\n';

  if (opt_mode_ == Mode::AREA_3) {
    script << "choice2\n";  // << "scleanup" << std::endl;
  } else {
    script << "resyn2\n";  // << "scleanup" << std::endl;
  }

  switch (opt_mode_) {
    case Mode::DELAY_1: {
      script << "map -D 0.01 -A 0.9 -B 0.2 -M 0 -p\n";
      script << "buffer -p -c\n";
      break;
    }
    case Mode::DELAY_2: {
      script << "choice\n";
      script << "map -D 0.01 -A 0.9 -B 0.2 -M 0 -p\n";
      script << "choice\n";
      script << "map -D 0.01\n";
      script << "buffer -p -c\n"
             << "topo\n";
      break;
    }
    case Mode::DELAY_3: {
      script << "choice2\n";
      script << "map -D 0.01 -A 0.9 -B 0.2 -M 0 -p\n";
      script << "choice2\n";
      script << "map -D 0.01\n";
      script << "buffer -p -c\n"
             << "topo\n";
      break;
    }
    case Mode::DELAY_4: {
      script << "choice2\n";
      script << "amap -F 20 -A 20 -C 5000 -Q 0.1 -m\n";
      script << "choice2\n";
      script << "map -D 0.01 -A 0.9 -B 0.2 -M 0 -p\n";
      script << "buffer -p -c\n";
      break;
    }
    case Mode::AREA_2:
    case Mode::AREA_3: {
      script << "choice2\n";
      script << "amap -m -Q 0.1 -F 20 -A 20 -C 5000\n";
      script << "choice2\n";
      script << "amap -m -Q 0.1 -F 20 -A 20 -C 5000\n";
      break;
    }
    case Mode::AREA_1:
    default: {
      script << "choice2\n";
      script << "amap -m -Q 0.1 -F 20 -A 20 -C 5000\n";
      break;
    }
  }
}

void Restructure::setMode(const char* mode_name)
{
  is_area_mode_ = true;

  if (!strcmp(mode_name, "timing")) {
    is_area_mode_ = false;
    opt_mode_ = Mode::DELAY_1;
  } else if (!strcmp(mode_name, "area")) {
    opt_mode_ = Mode::AREA_1;
  } else {
    logger_->warn(RMP, 10, "Mode {} not recognized.", mode_name);
  }
}

void Restructure::setTieHiPort(sta::LibertyPort* tieHiPort)
{
  if (tieHiPort) {
    hicell_ = tieHiPort->libertyCell()->name();
    hiport_ = tieHiPort->name();
  }
}

void Restructure::setTieLoPort(sta::LibertyPort* tieLoPort)
{
  if (tieLoPort) {
    locell_ = tieLoPort->libertyCell()->name();
    loport_ = tieLoPort->name();
  }
}

void Restructure::writeInstanceCoordinates(const std::string& file_name)
{
  std::ofstream coord_file(file_name.c_str());
  if (!coord_file.is_open()) {
    logger_->error(RMP, 4, "Cannot open coordinate file {} for writing.", file_name);
    return;
  }

  // Write header comment
  coord_file << "# Instance coordinates for ABC-aware remapping\n";
  coord_file << "# Format: instance_name x y\n";

  // Write number of instances
  coord_file << path_insts_.size() << "\n";

  // Write each instance coordinate
  for (auto inst : path_insts_) {
    const odb::Point origin = inst->getOrigin();
    coord_file << inst->getName() << " " << origin.x() << " " << origin.y() << "\n";
  }

  coord_file.close();
  debugPrint(logger_, RMP, "remap", 1,
             "Wrote {} instance coordinates to {}",
             path_insts_.size(), file_name);
}

bool Restructure::readAbcLog(const std::string& abc_file_name,
                             int& level_gain,
                             float& final_delay)
{
  std::ifstream abc_file(abc_file_name);
  if (abc_file.bad()) {
    logger_->error(RMP, 2, "cannot open file {}", abc_file_name);
    return false;
  }
  debugPrint(
      logger_, utl::RMP, "remap", 1, "Reading ABC log {}.", abc_file_name);
  std::string buf;
  const char delimiter = ' ';
  bool status = true;
  std::vector<double> level;
  std::vector<float> delay;

  // read the file line by line
  while (std::getline(abc_file, buf)) {
    // convert the line in to stream:
    std::istringstream ss(buf);
    std::vector<std::string> tokens;

    // read the line, word by word
    while (std::getline(ss, buf, delimiter)) {
      tokens.push_back(buf);
    }

    if (!tokens.empty() && tokens[0] == "Error:") {
      status = false;
      logger_->warn(RMP,
                    5,
                    "ABC run failed, see log file {} for details.",
                    abc_file_name);
      break;
    }
    if (tokens.size() > 7 && tokens[tokens.size() - 3] == "lev"
        && tokens[tokens.size() - 2] == "=") {
      level.emplace_back(std::stoi(tokens[tokens.size() - 1]));
    }
    if (tokens.size() > 7) {
      std::string prev_token;
      for (std::string token : tokens) {
        if (prev_token == "delay" && token.at(0) == '=') {
          std::string delay_str = token;
          if (delay_str.size() > 1) {
            delay_str.erase(
                delay_str.begin());  // remove first char which is '='
            delay.emplace_back(std::stof(delay_str));
          }
          break;
        }
        prev_token = std::move(token);
      }
    }
  }

  if (level.size() > 1) {
    level_gain = level[0] - level[level.size() - 1];
  }
  if (!delay.empty()) {
    final_delay = delay[delay.size() - 1];  // last value in file
  }
  return status;
}
}  // namespace rmp
