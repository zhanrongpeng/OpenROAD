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
#include <map>
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

  // Priority: use estimate_parasitics_ RC (set by set_wire_rc / set_layer_rc).
  // This is the same singleton that the TCL commands set, so after
  // "set_wire_rc -layer metal3" it contains the correct per-layer RC.
  wire_r_per_um_ = 0.0;
  wire_c_per_um_ = 0.0;
  // ─── Wire RC for ABC's Elmore delay model ─────────────────────────────────
  //
  // ABC (ifTime.c) uses wire_delay = R(ohm/um) * C(fF/um) * L(um)^2 * 0.5 * 1e-3.
  //
  // We read the raw UI values (ohm/um, fF/um) that were stored by set_layer_rc
  // via the new set_wire_rc_um_cmd() API, bypassing OpenROAD's internal unit
  // conversions (which apply the library's kohm/ff scale to the resistance value).
  //
  // Elmore verification with L=30.35um, R=3.574e-3 ohm/um, C=7.516e-2 fF/um:
  //   wire_delay = 3.574e-3 * 7.516e-2 * (30.35)^2 * 0.5 * 1e-3
  //              = 0.1237 ps  ← realistic wire delay
  //
  if (estimate_parasitics_) {
    sta::Corners* corners = open_sta_->corners();
    if (corners && corners->count() > 0) {
      sta::Corner* corner = *corners->begin();
      if (corner) {
        double r_ohm_per_um = estimate_parasitics_->wireSignalResistanceUm(corner);
        double c_ff_per_um  = estimate_parasitics_->wireSignalCapacitanceUf(corner);
        if (r_ohm_per_um > 0.0 && c_ff_per_um > 0.0) {
          wire_r_per_um_ = r_ohm_per_um;  // already ohm/um
          wire_c_per_um_ = c_ff_per_um;   // already fF/um
          logger_->report(
              "[DEBUG] Wire RC from set_layer_rc (UI units): "
              "R={:.4e} ohm/um, C={:.4e} fF/um",
              wire_r_per_um_, wire_c_per_um_);
        } else {
          logger_->report(
              "[DEBUG] set_wire_rc_um not set, falling back to DB tech layer RC.");
        }
      }
    }
  }

  // Fallback: read from DB tech layer(s) directly.
  // NOTE: this path reads the RC values written by set_dblayer_wire_rc, which
  // are also contaminated by the kohm unit conversion.  Use only as last resort.
  if (wire_r_per_um_ == 0.0 && wire_c_per_um_ == 0.0) {
    odb::dbTech* tech = db_->getTech();
    odb::dbBlock* block = db_->getChip() ? db_->getChip()->getBlock() : nullptr;
    if (tech && block) {
      std::vector<odb::dbTechLayer*> signal_layers
          = estimate_parasitics_ ? estimate_parasitics_->signalLayers()
                                 : std::vector<odb::dbTechLayer*>();
      if (!signal_layers.empty()) {
        for (odb::dbTechLayer* layer : signal_layers) {
          if (!layer || layer->getType() != odb::dbTechLayerType::ROUTING)
            continue;
          const float layer_width_um
              = block->dbuToMicrons(static_cast<int>(layer->getWidth()));
          if (layer_width_um <= 0.0f)
            continue;
          const float r_per_sq   = layer->getResistance();           // ohm/sq
          const float cap_per_sq  = layer->getCapacitance();          // F/um²
          const float cap_edge    = layer->getEdgeCapacitance();      // F/um
          wire_r_per_um_ = r_per_sq / layer_width_um;                 // ohm/um
          const float c_per_um = layer_width_um * cap_per_sq
                                 + 2.0f * cap_edge;                   // F/um
          wire_c_per_um_ = c_per_um * 1e15f;                         // fF/um
          if (wire_r_per_um_ > 0.0f && wire_c_per_um_ > 0.0f) {
            logger_->report(
                "[DEBUG] Wire RC from layer \"{}\" (DB, fallback): "
                "R={:.4e} ohm/um, C={:.4e} fF/um",
                layer->getConstName(), wire_r_per_um_, wire_c_per_um_);
            break;
          }
        }
      }

      // Last resort: average over all routing layers
      if (wire_r_per_um_ == 0.0 && wire_c_per_um_ == 0.0) {
        int routing_layer_count = tech->getRoutingLayerCount();
        int n_layers = 0;
        double sum_r = 0.0;
        double sum_c = 0.0;
        for (int i = 1; i <= routing_layer_count; i++) {
          odb::dbTechLayer* layer = tech->findRoutingLayer(i);
          if (!layer || layer->getType() != odb::dbTechLayerType::ROUTING)
            continue;
          const float layer_width_um
              = block->dbuToMicrons(static_cast<int>(layer->getWidth()));
          if (layer_width_um <= 0.0f)
            continue;
          const float r_per_sq   = layer->getResistance();
          const float cap_per_sq  = layer->getCapacitance();
          const float cap_edge    = layer->getEdgeCapacitance();
          const float c_per_um = layer_width_um * cap_per_sq + 2.0f * cap_edge;
          sum_r += (r_per_sq / layer_width_um);
          sum_c += (c_per_um * 1e15f);
          n_layers++;
        }
        if (n_layers > 0) {
          wire_r_per_um_ = sum_r / n_layers;
          wire_c_per_um_ = sum_c / n_layers;
          logger_->report(
              "[DEBUG] Wire RC from DB tech ({} layer avg, fallback): "
              "R={:.4e} ohm/um, C={:.4e} fF/um",
              n_layers, wire_r_per_um_, wire_c_per_um_);
        }
      }
    }
  }

  if (wire_r_per_um_ == 0.0 || wire_c_per_um_ == 0.0) {
    logger_->warn(RMP, 99,
        "Wire RC is zero! R={:.4e} ohm/um, C={:.4e} fF/um. "
        "ABC will use fallback WireDelay coefficient. "
        "Ensure set_wire_rc or set_layer_rc was called before restructure.",
        wire_r_per_um_, wire_c_per_um_);
  }

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

  // Write net coordinates BEFORE writeBlif, because writeBlif disconnects iterms
  logger_->report("[DBG] Calling writeNetCoordinates, path_insts_ size: {}", path_insts_.size());
  writeNetCoordinates(coord_file_name_);
  logger_->report("[DBG] Done writeNetCoordinates");

  blif_.writeBlif(input_blif_file_name_.c_str(), !is_area_mode_);
  debugPrint(
      logger_, RMP, "remap", 1, "Writing blif file {}", input_blif_file_name_);
  files_to_remove.emplace_back(input_blif_file_name_);

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

  debugPrint(
      logger_, RMP, "remap", 1, "Running ABC with {} modes.", modes.size());

  for (size_t curr_mode_idx = 0; curr_mode_idx < modes.size();
       curr_mode_idx++) {
    output_blif_file_name_
        = prefix + "_" + std::to_string(curr_mode_idx) + "_crit_path_out.blif";

    opt_mode_ = modes[curr_mode_idx];

    const std::string abc_script_file
        = prefix + "_" + std::to_string(curr_mode_idx) + "_ord_abc_script.tcl";
    if (logfile_.empty()) {
      logfile_ = prefix + "abc.log";
    }

    debugPrint(logger_,
               RMP,
               "remap",
               1,
               "Writing ABC script file {}.",
               abc_script_file);

    logger_->report("[DBG] ABC script file: {}", abc_script_file);
    {
      std::ifstream sf(abc_script_file.c_str());
      if (sf.good()) {
        std::string line;
        logger_->report("[DBG] ABC script content:");
        while (std::getline(sf, line) && line.size() < 200) {
          logger_->report("  {}", line);
        }
      } else {
        logger_->report("[DBG] ABC script file NOT found!");
      }
    }

    logger_->report("[DBG] Calling writeAbcScript: {}", abc_script_file);
    if (writeAbcScript(abc_script_file)) {
      logger_->report("[DBG] writeAbcScript succeeded");
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
        logger_->error(RMP, 6, "Error executing ABC command {}. Return code: {}.",
                       command, child_proc[curr_mode_idx]);
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
        = prefix + "_" + std::to_string(curr_mode_idx) + "_crit_path_out.blif";
    const std::string abc_log_name = logfile_;  // ABC always writes to logfile_ (no suffix)
    logger_->report("[DBG] Checking ABC output BLIF: {}", output_blif_file_name_);
    {
      std::ifstream test_f(output_blif_file_name_.c_str());
      if (test_f.good()) {
        test_f.seekg(0, std::ios::end);
        long long fsize = test_f.tellg();
        logger_->report("[DBG]   File size: {} bytes", fsize);
      } else {
        logger_->report("[DBG]   File does NOT exist or cannot be opened!");
      }
    }

    int level_gain = 0;
    float delay = std::numeric_limits<float>::max();
    int num_instances = 0;
    readAbcLog(abc_log_name, level_gain, delay);

    bool blif_ok = blif_.inspectBlif(output_blif_file_name_.c_str(), num_instances);
    logger_->report(
        "Optimized to {} instances in iteration {} with max path depth "
        "decrease of {}, delay of {}.",
        num_instances,
        curr_mode_idx,
        level_gain,
        delay);

    if (blif_ok && num_instances > 0) {
      if (is_area_mode_) {
        if (num_instances < best_inst_count) {
          best_inst_count = num_instances;
          best_blif = output_blif_file_name_;
        }
      } else {
        // Using only DELAY_4 for delay based gain since other modes not
        // showing good gains
        if (modes[curr_mode_idx] == Mode::DELAY_4) {
          best_blif = output_blif_file_name_;
        }
      }
    }

    files_to_remove.emplace_back(output_blif_file_name_);
  }

  // Apply best result when we actually picked an output BLIF.
  //
  // Bug fix (timing mode): the old condition used best_delay_gain < FLT_MAX,
  // but readAbcLog() often leaves delay at FLT_MAX when the log has no
  // "delay = ..." line. DELAY_4 could still set best_blif after inspectBlif()
  // succeeded; we must not skip readBlif in that case.
  if (!best_blif.empty()) {
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

  // Calculate and report wire delay, cell delay, and total delay.
  // For each net in the cone:
  //   - Driver position: output pin center via getAvgXY()  (um)
  //   - Wirelength to each fanout: Manhattan distance      (um)
  //   - Wire delay (Elmore): R * C * L^2 / 2              (ps)
  //   - Cell delay: vertex arrival time                    (ns → ps)
  //   - Total: cell_delay + worst_wire_delay per net
  if (!is_area_mode_ && !path_insts_.empty()) {
    int dbu_per_um = block_->getDbUnitsPerMicron();

    // Use actual wire RC from OpenROAD; fall back to 0.1 ohm/um / 0.1 fF/um if unavailable
    double r_ohm_um = (wire_r_per_um_ > 0.0) ? wire_r_per_um_ : 0.1;
    double c_ff_um  = (wire_c_per_um_ > 0.0) ? wire_c_per_um_ : 0.1;

    double total_cell_delay_ps = 0.0;
    double total_wire_delay_ps = 0.0;
    int net_count = 0;
    int driver_count = 0;

    // Deduplicate by net so each net contributes once (one driver, worst fanout distance)
    std::set<odb::dbNet*> seen_nets;
    std::set<sta::Vertex*> seen_vertices;

    open_sta_->ensureGraph();
    open_sta_->searchPreamble();
    sta::Graph* graph = open_sta_->graph();

    for (auto inst : path_insts_) {
      for (odb::dbITerm* iterm : inst->getITerms()) {
        if (!iterm->isConnected())
          continue;
        if (iterm->getIoType() != odb::dbIoType::OUTPUT)
          continue;
        if (iterm->getSigType() == odb::dbSigType::POWER
            || iterm->getSigType() == odb::dbSigType::GROUND)
          continue;

        odb::dbNet* net = iterm->getNet();
        if (!net)
          continue;
        if (seen_nets.count(net))
          continue;  // already processed this net
        seen_nets.insert(net);

        const sta::Pin* pin = open_sta_->getDbNetwork()->dbToSta(iterm);
        if (!pin)
          continue;
        sta::PortDirection* dir = open_sta_->getDbNetwork()->direction(pin);
        if (!dir || !dir->isOutput())
          continue;

        // Driver position: output pin center via getAvgXY()
        int dpx, dpy;
        if (!iterm->getAvgXY(&dpx, &dpy))
          continue;
        double driver_x = static_cast<double>(dpx) / dbu_per_um;
        double driver_y = static_cast<double>(dpy) / dbu_per_um;

        // Worst fanout wirelength (Manhattan distance)
        double worst_wirelength_um = 0.0;
        for (odb::dbITerm* other_iterm : net->getITerms()) {
          if (other_iterm == iterm)
            continue;
          odb::dbInst* other_inst = other_iterm->getInst();
          if (!other_inst)
            continue;

          int opx, opy;
          if (!other_iterm->getAvgXY(&opx, &opy))
            continue;

          double ox = static_cast<double>(opx) / dbu_per_um;
          double oy = static_cast<double>(opy) / dbu_per_um;
          double dist = std::abs(driver_x - ox) + std::abs(driver_y - oy);
          if (dist > worst_wirelength_um)
            worst_wirelength_um = dist;
        }

        // Elmore wire delay: R * C * L^2 / 2  (ps)
        // R [ohm/um] * C [fF/um] * L^2 [um^2] * 0.5 * 1e-3 = fF*ohm*um^2/um^2*0.001 = pC*ohm = ps
        double wire_delay_ps = r_ohm_um * c_ff_um
                             * worst_wirelength_um * worst_wirelength_um * 0.5 * 1e-3;
        total_wire_delay_ps += wire_delay_ps;
        net_count++;

        // Cell delay: vertex arrival time (ns → ps)
        sta::Vertex* vertex = graph->pinLoadVertex(pin);
        if (vertex && !seen_vertices.count(vertex)) {
          seen_vertices.insert(vertex);
          sta::Arrival arrival = open_sta_->vertexArrival(vertex, sta::MinMax::max());
          total_cell_delay_ps += sta::delayAsFloat(arrival) * 1e9;
          driver_count++;
        }
      }
    }

    if (net_count > 0) {
      double avg_cell_delay_ps = total_cell_delay_ps / driver_count;
      double avg_wire_delay_ps = total_wire_delay_ps / net_count;
      double total_delay_ps = avg_cell_delay_ps + avg_wire_delay_ps;

      logger_->report("Delay Analysis (Elmore wire model):");
      logger_->report("  Wire RC: R={:.4e} ohm/um, C={:.4e} fF/um", r_ohm_um, c_ff_um);
      logger_->report("  Nets processed: {}", net_count);
      logger_->report("  Drivers processed: {}", driver_count);
      logger_->report("  Average Cell Delay: {:.2f} ps", avg_cell_delay_ps);
      logger_->report("  Average Wire Delay (Elmore): {:.2f} ps", avg_wire_delay_ps);
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
  
  logger_->report("[DBG] writeAbcScript: file opened successfully");

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

  // Physical-aware mapping: skip rewrite/refactor
  // Wire RC is set directly via Abc_FrameSetWireRC() in C++ before sourcing this script
  script << "strash\n";

  // Read coordinates for wire-aware mapping (sets pAbc->pNtkCoords and wire RC)
  script << "read_coords " << coord_file_name_ << '\n';

  // Use if with wire-aware mapping (-W 0.001).
  // if (FPGA mapper) outputs LUTs using cut-based evaluation with tdelay from SCL library.
  // -W 0.001: wire delay coefficient = R(ohm/um) * C(fF/um) * 0.5e-3.
  //   Combined cut delay = tdelay(cut) + WireDelay_coef * wirelength_um.
  //   Wirelength is Manhattan distance from leaf centroid to fanin.
  // NOTE: map -D 0.01 -W has a bug where -W is parsed as WireDelay argument (needs value after -W),
  //   causing WireDelay=0 and disabling wire-aware entirely. if -W 0.001 works correctly
  //   and has well-tested coordinate lookup via vNodeNameMap + vNodeDrivingPoName.
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

void Restructure::writeNetCoordinates(const std::string& file_name)
{
  std::ofstream coord_file(file_name.c_str());
  if (!coord_file.is_open()) {
    logger_->error(RMP, 4, "Cannot open coordinate file {} for writing.", file_name);
    return;
  }

  int dbu_per_um = block_->getDbUnitsPerMicron();

  // Collect all nets connected to instances in the cone, then write the
  // physical position of each net (net_name -> any connected pin position).
  //
  // Coordinate lookup in ABC (ifTime.c):
  //   - vNodeNameMap[leafId] -> blif node -> Abc_ObjName() = net name
  //   - vNodeDrivingPoName[nodeId] = PO name for CIs driving POs
  // So every net in the cone needs an entry, including:
  //   (a) OUTPUT nets  driven by a cell in path_insts_  -> getAvgXY on OUTPUT iterm
  //   (b) INPUT/PI nets whose fanin cell is outside cone -> getAvgXY on INPUT iterm
  //       (the input pin is on a cell inside path_insts_)
  //
  // Since getAvgXY() returns the pin's routing-shape center, using INPUT iterm
  // position is equally valid (both iterms of the same net are close together).
  std::set<odb::dbNet*> cone_nets;

  for (auto inst : path_insts_) {
    int iterm_count = 0;
    int iterm_connected = 0;
    int iterm_with_net = 0;
    for (odb::dbITerm* iterm : inst->getITerms()) {
      iterm_count++;
      odb::dbNet* net = iterm->getNet();
      if (!net) {
        continue;
      }
      iterm_connected++;
      if (iterm->getSigType() == odb::dbSigType::POWER
          || iterm->getSigType() == odb::dbSigType::GROUND)
        continue;
      if (net) {
        cone_nets.insert(net);
        iterm_with_net++;
      }
    }
    if (iterm_count == 0) {
      logger_->report("[DBG] Instance {} has 0 iterms", inst->getName());
    } else if (iterm_connected == 0) {
      logger_->report("[DBG] Instance {} has {} iterms, 0 connected", inst->getName(), iterm_count);
    }
  }
  logger_->report("[DBG] After loop: cone_nets={}", cone_nets.size());

  std::map<std::string, std::pair<double, double>> net_coords;
  int nets_with_position = 0;
  int nets_without_position = 0;

  for (odb::dbNet* net : cone_nets) {
    std::string netName = net->getName();

    // Try all connected iterms (input and output); use whichever has getAvgXY
    int px = 0, py = 0;
    bool found = false;
    for (odb::dbITerm* iterm : net->getITerms()) {
      if (iterm->getSigType() == odb::dbSigType::POWER
          || iterm->getSigType() == odb::dbSigType::GROUND)
        continue;
      if (iterm->getAvgXY(&px, &py)) {
        found = true;
        break;  // first hit is sufficient (all iterms of same net are close)
      }
    }

    if (found) {
      static int dbg_net = 0;
      if (dbg_net < 5) {
        dbg_net++;
        logger_->report("[DBG] writeNetCoordinates: netName=\"{}\" pin=({},{})",
            netName, px, py);
      }
      net_coords[netName] = {
        static_cast<double>(px) / dbu_per_um,
        static_cast<double>(py) / dbu_per_um
      };
      nets_with_position++;
    } else {
      nets_without_position++;
    }
  }

  logger_->report("[DBG] writeNetCoordinates: cone_nets={}, with_position={}, without_position={}",
      cone_nets.size(), nets_with_position, nets_without_position);

  coord_file << "# Net-to-output-pin coordinates for ABC wire-aware mapping\n";
  coord_file << "#wire_rc " << std::scientific << wire_r_per_um_ << " " << wire_c_per_um_ << "\n";
  coord_file << "# Format: net_name x(um) y(um)\n";
  coord_file << net_coords.size() << "\n";

  for (const auto& [net_name, coord] : net_coords) {
    coord_file << net_name << " " << coord.first << " " << coord.second << "\n";
  }

  coord_file.close();
  debugPrint(logger_, RMP, "remap", 1,
             "Wrote {} net->output_pin coordinates to {} (DBU/um={}, wire RC: R={:.3e} ohm/um, C={:.3e} fF/um)",
             net_coords.size(), file_name, dbu_per_um, wire_r_per_um_, wire_c_per_um_);
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
