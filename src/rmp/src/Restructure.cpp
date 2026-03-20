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
#include <cmath>
#include <queue>

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
#include "rsz/Resizer.hh"
#include "sta/Delay.hh"
#include "sta/Graph.hh"
#include "sta/GraphClass.hh"
#include "sta/Liberty.hh"
#include "sta/Corner.hh"
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
using abc::Abc_FrameSetWireRC;
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

  cone_config_.delay_threshold_ratio = 0.3f;
  cone_config_.max_cone_depth = 15;
  cone_config_.max_cone_size = 500;
  cone_config_.min_improvement_threshold = 0.05f;
  cone_config_.distance_threshold_um = 50.0f;
  max_depth_ = 16;

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
  max_depth_ = max_depth;

  lib_file_names_.emplace_back(liberty_file_name);
  work_dir_name_ = workdir_name;
  work_dir_name_ = work_dir_name_ + "/";

  if (is_area_mode_) {  // Only in area mode
    removeConstCells();
  }

  getBlob(max_depth_);

  if (!path_insts_.empty()) {
    runABC();

    postABC(worst_slack);
  }
}

std::pair<float, float> Restructure::getInstanceCoordinateFromPin(const sta::Pin* pin) {
  if (!pin) {
    return {0.0f, 0.0f};
  }
  
  odb::dbITerm* term = nullptr;
  odb::dbBTerm* port = nullptr;
  odb::dbModITerm* moditerm = nullptr;
  open_sta_->getDbNetwork()->staToDb(pin, term, port, moditerm);
  
  if (term) {
    odb::dbInst* inst = term->getInst();
    odb::dbBox* bbox = inst->getBBox();
    if (bbox) {
      float x = (bbox->xMin() + bbox->xMax()) / 2.0;
      float y = (bbox->yMin() + bbox->yMax()) / 2.0;
      return {x, y};
    }
  } else if (port) {
    odb::Rect rect = port->getBBox();
    float x = (rect.xMin() + rect.xMax()) / 2.0;
    float y = (rect.yMin() + rect.yMax()) / 2.0;
    return {x, y};
  }
  
  return {0.0f, 0.0f};
}

Restructure::PathNode* Restructure::buildPathTree(
    const sta::Pin* pin,
    int current_depth,
    const ConeSelectionConfig& config) {
  
  if (!pin || current_depth >= config.max_cone_depth) {
    return nullptr;
  }
  
  sta::Graph* graph = open_sta_->graph();
  sta::Vertex* vertex = graph->pinLoadVertex(pin);
  
  if (!vertex) {
    return nullptr;
  }
  
  PathNode* node = new PathNode();
  node->vertex = vertex;
  node->pin = const_cast<sta::Pin*>(pin);
  node->arrival_time = getVertexArrivalTime(vertex);
  node->slack = getVertexSlack(vertex);
  node->depth = current_depth;
  
  auto coords = getInstanceCoordinateFromPin(pin);
  node->x = coords.first;
  node->y = coords.second;
  
  node->left_child = nullptr;
  node->right_child = nullptr;
  node->selected = false;
  node->is_boundary = false;
  
  // 使用 resizer 的 findFanins 来获取 fanin pins
  sta::PinSet fanin_pins_set(open_sta_->getDbNetwork());
  fanin_pins_set.insert(pin);
  sta::PinSet fanins = resizer_->findFanins(fanin_pins_set);
  
  // 去掉输入pin本身
  for (auto fanin : fanins) {
    if (fanin != pin) {
      fanins.erase(fanin);
    }
  }
  
  // 重新获取真正的fanins
  sta::PinSet real_fanins(open_sta_->getDbNetwork());
  for (auto fanin : fanins) {
    if (fanin != pin) {
      real_fanins.insert(fanin);
    }
  }
  
  // 转为 vector 以便索引
  std::vector<const sta::Pin*> fanin_vec;
  for (auto fanin : real_fanins) {
    fanin_vec.push_back(fanin);
  }
  
  if (fanin_vec.empty()) {
    node->is_boundary = true;
    return node;
  }
  
  // 最多取2个fanin
  if (fanin_vec.size() >= 1) {
    node->left_child = buildPathTree(fanin_vec[0], current_depth + 1, config);
  }
  if (fanin_vec.size() >= 2) {
    node->right_child = buildPathTree(fanin_vec[1], current_depth + 1, config);
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
  
  // 计算延迟差异
  float delay_gap = calculateDelayGap(left, right);
  float parent_at = parent->arrival_time;
  
  // 避免除以0
  if (parent_at <= 0.0f) {
    parent_at = 1.0f;
  }
  
  float delay_gap_ratio = delay_gap / parent_at;
  
  // 计算物理距离
  float left_dist = calculatePhysicalDistance(parent, left);
  float right_dist = calculatePhysicalDistance(parent, right);
  
  debugPrint(logger_, RMP, "remap", 1,
             "Checking children: delay_gap={:.4f}, ratio={:.2%}, "
             "left_dist={:.2f}um, right_dist={:.2f}um, threshold={:.2f}um",
             delay_gap, delay_gap_ratio, left_dist, right_dist, 
             config.distance_threshold_um);
  
  // 如果两个子节点都在物理距离阈值内
  if (left_dist <= config.distance_threshold_um && 
      right_dist <= config.distance_threshold_um) {
    if (delay_gap_ratio <= config.delay_threshold_ratio) {
      // 差异小，两个都选
      debugPrint(logger_, RMP, "remap", 1,
                 "Both within distance, small delay gap - select both");
      return true;
    } else {
      // 差异大，两个都选（因为物理距离都在阈值内）
      debugPrint(logger_, RMP, "remap", 1,
                 "Both within distance, large delay gap - still select both");
      return true;
    }
  }
  
  // 至少有一个子节点超出物理距离阈值
  bool left_ok = left_dist <= config.distance_threshold_um;
  bool right_ok = right_dist <= config.distance_threshold_um;
  
  if (left_ok && !right_ok) {
    // 只有左边在阈值内
    debugPrint(logger_, RMP, "remap", 1, "Only left child within distance");
    return true;  // 选中left
  } else if (!left_ok && right_ok) {
    // 只有右边在阈值内
    debugPrint(logger_, RMP, "remap", 1, "Only right child within distance");
    return true;  // 选中right
  } else {
    // 都不在阈值内，需要选择一个更近的
    debugPrint(logger_, RMP, "remap", 1, 
               "Neither child within distance, selecting closer one");
    if (left_dist < right_dist) {
      return true;  // 选中left
    } else {
      return true;  // 选中right（如果相等）
    }
  }
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

  // Compute wire RC from estimate_parasitics_
  // wireSignalResistance() returns kohm/m (from set_layer_rc in kohm/um / 1um)
  // wireSignalCapacitance() returns pf/m (from set_layer_rc in pf/um / 1um)
  // Nangate45.rc provides: R in kohm/um, C in pf/um - already the correct per-length values.
  wire_r_per_um_ = 0.0;
  wire_c_per_um_ = 0.0;
  if (estimate_parasitics_) {
    sta::Corners* corners = open_sta_->corners();
    if (corners && corners->count() > 0) {
      sta::Corner* corner = *corners->begin();
      if (corner) {
        // kohm/m → ohm/um: multiply by 1e-3
        // pf/m → fF/um: multiply by 1e-3
        double r_kohm_per_m = estimate_parasitics_->wireSignalResistance(corner);
        double c_pf_per_m = estimate_parasitics_->wireSignalCapacitance(corner);
        wire_r_per_um_ = r_kohm_per_m * 1e-3;  // kohm/m → ohm/um
        wire_c_per_um_ = c_pf_per_m * 1e-3;    // pf/m → fF/um
        logger_->report(
            "[DEBUG] Wire RC: r_kohm/m={:.4e}, c_pf/m={:.4e} → R={:.4e} ohm/um, C={:.4e} fF/um",
            r_kohm_per_m, c_pf_per_m, wire_r_per_um_, wire_c_per_um_);
      }
    }
  }
  logger_->report(
      "[DEBUG] Wire RC from OpenROAD: R={:.4e} ohm/um, C={:.4e} fF/um",
      wire_r_per_um_, wire_c_per_um_);

  sta::PinSet ends(open_sta_->getDbNetwork());

  getEndPoints(ends, is_area_mode_, max_depth);
  
  if (!ends.empty()) {
    // Area mode 和 Delay mode 都使用 findFaninFanouts
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

      // Set wire RC from OpenROAD directly (bypass coords file)
      Abc_FrameSetWireRC(
          static_cast<float>(wire_r_per_um_),
          static_cast<float>(wire_c_per_um_));

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
  
  logger_->report("postABC called - is_area_mode: {}, path_insts size: {}", 
                  is_area_mode_, path_insts_.size());
  
  // Calculate and report wire delay, cell delay, and total delay
  if (!is_area_mode_ && !path_insts_.empty()) {
    double total_cell_delay = 0.0;
    double total_wire_delay = 0.0;
    int cell_count = 0;
    
    open_sta_->ensureGraph();
    open_sta_->searchPreamble();
    sta::Graph* graph = open_sta_->graph();
    
    // For each instance in the path, get timing information
    for (auto inst : path_insts_) {
      odb::dbInst* db_inst = inst;
      
      // Iterate through all pins of this instance
      for (odb::dbITerm* iterm : db_inst->getITerms()) {
        if (!iterm->isConnected())
          continue;
        
        const sta::Pin* pin = open_sta_->getDbNetwork()->dbToSta(iterm);
        if (!pin)
          continue;
        
        // Skip input pins - only get delays from output pins
        sta::PortDirection* dir = open_sta_->getDbNetwork()->direction(pin);
        if (!dir || !dir->isOutput())
          continue;
        
        // Get the vertex for this pin
        sta::Vertex* vertex = graph->pinLoadVertex(pin);
        if (!vertex)
          continue;
        
        // Get cell delay (from vertex arrival)
        sta::Arrival arrival = open_sta_->vertexArrival(vertex, sta::MinMax::max());
        
        // Get wire delay from the net connected to this pin
        odb::dbNet* net = iterm->getNet();
        if (!net)
          continue;
        
        // Estimate wire length from instance coordinates
        odb::dbBox* bbox = db_inst->getBBox();
        if (bbox) {
          double inst_x = (bbox->xMin() + bbox->xMax()) / 2.0;
          double inst_y = (bbox->yMin() + bbox->yMax()) / 2.0;
          
          // Find connected instances to estimate wire length
          for (odb::dbITerm* other_iterm : net->getITerms()) {
            if (other_iterm == iterm)
              continue;
            
            odb::dbInst* other_inst = other_iterm->getInst();
            if (!other_inst)
              continue;
            
            odb::dbBox* other_bbox = other_inst->getBBox();
            if (!other_bbox)
              continue;
            
            double other_x = (other_bbox->xMin() + other_bbox->xMax()) / 2.0;
            double other_y = (other_bbox->yMin() + other_bbox->yMax()) / 2.0;
            
            // Manhattan distance (approximate wire length in um)
            double wire_length_um = std::abs(inst_x - other_x) + std::abs(inst_y - other_y);
            
            // Wire delay = R * C
            // Approximate R ~ 0.1 ohm/um, C ~ 0.1 fF/um
            // wire_delay (ps) = R(ohm/um) * length(um) * C(fF/um)
            double wire_delay_ps = wire_length_um * 0.1 * 0.1; 
            
            total_wire_delay += wire_delay_ps;
          }
        }
        
        total_cell_delay += arrival * 1e9; // Convert ns to ps
        cell_count++;
      }
    }
    
    if (cell_count > 0) {
      double avg_cell_delay_ps = total_cell_delay / cell_count;
      double avg_wire_delay_ps = total_wire_delay / cell_count;
      double total_delay_ps = avg_cell_delay_ps + avg_wire_delay_ps;
      
      logger_->report("Delay Analysis:");
      logger_->report("  Average Cell Delay: {:.2f} ps", avg_cell_delay_ps);
      logger_->report("  Average Wire Delay: {:.2f} ps", avg_wire_delay_ps);
      logger_->report("  Average Total Delay: {:.2f} ps", total_delay_ps);
    }
  }
}
void Restructure::getEndPoints(sta::PinSet& ends,
                               bool area_mode,
                               unsigned max_depth)
{
  auto sta_state = open_sta_->search();
  sta::VertexSet* end_points_ptr = sta_state->endpoints();
  for (auto& end_point : *end_points_ptr) {
    if (!is_area_mode_) {
      sta::Path* path
          = open_sta_->vertexWorstSlackPath(end_point, sta::MinMax::max());
      sta::PathExpanded expanded(path, open_sta_);
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
    auto& errors = open_sta_->checkTiming(
        false /*no_input_delay*/,
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

  // Coordinates are stored in DBU (usually nanometers), convert to microns for ABC
  int dbu_per_um = block_->getDbUnitsPerMicron();

  // Write header comment with wire RC (read_coords parses "# wire_rc R C")
  // wire_r_per_um_ is already in ohm/um, wire_c_per_um_ is already in fF/um
  coord_file << "# Instance coordinates for ABC-aware remapping\n";
  coord_file << "# wire_rc " << std::scientific << wire_r_per_um_ << " " << wire_c_per_um_ << "\n";
  coord_file << "# Format: instance_name x(um) y(um)\n";

  // Write number of instances
  coord_file << path_insts_.size() << "\n";

  // Write each instance coordinate (in microns)
  for (auto inst : path_insts_) {
    const odb::Point origin = inst->getOrigin();
    coord_file << inst->getName() << " "
               << (double)origin.x() / dbu_per_um << " "
               << (double)origin.y() / dbu_per_um << "\n";
  }

  coord_file.close();
  debugPrint(logger_, RMP, "remap", 1,
             "Wrote {} instance coordinates to {} (DBU/um={}, wire RC: R={:.3e} ohm/um, C={:.3e} fF/um)",
             path_insts_.size(), file_name, dbu_per_um, wire_r_per_um_, wire_c_per_um_);
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
