// SPDX-License-Identifier: BSD-3-Clause
// Copyright (c) 2019-2025, The OpenROAD Authors

#pragma once

#include <cstdint>
#include <fstream>
#include <functional>
#include <optional>
#include <random>
#include <set>
#include <string>
#include <vector>

#include "db_sta/dbSta.hh"
#include "rsz/Resizer.hh"
#include "sta/Corner.hh"
#include "sta/Delay.hh"
#include "sta/Liberty.hh"
#include "sta/NetworkClass.hh"
#include "utl/unique_name.h"

#include "base/abc/abc.h"
#include "base/main/abcapis.h"

namespace abc {
}  // namespace abc

namespace utl {
class Logger;
}

namespace odb {
class dbDatabase;
class dbBlock;
class dbInst;
class dbNet;
class dbITerm;
}  // namespace odb

namespace est {
class EstimateParasitics;
}

namespace sta {
class dbSta;
}  // namespace sta

namespace rmp {

using utl::Logger;

enum class Mode
{
  AREA_1 = 0,
  AREA_2,
  AREA_3,
  DELAY_1,
  DELAY_2,
  DELAY_3,
  DELAY_4
};

class Restructure
{
 public:
  Restructure(utl::Logger* logger,
              sta::dbSta* open_sta,
              odb::dbDatabase* db,
              rsz::Resizer* resizer,
              est::EstimateParasitics* estimate_parasitics);
  ~Restructure();

  void reset();
  void resynth(sta::Corner* corner);
  void resynthAnnealing(sta::Corner* corner);
  void run(char* liberty_file_name,
           float slack_threshold,
           unsigned max_depth,
           char* workdir_name,
           char* abc_logfile);

  void setAnnealingSeed(std::mt19937::result_type seed)
  {
    annealing_seed_ = seed;
  }
  void setAnnealingTemp(float temp) { annealing_temp_ = temp; }
  void setAnnealingIters(unsigned iters) { annealing_iters_ = iters; }
  void setAnnealingRevertAfter(unsigned revert_after)
  {
    annealing_revert_after_ = revert_after;
  }
  void setAnnealingInitialOps(unsigned ops) { annealing_init_ops_ = ops; }
  void setSlackThreshold(sta::Slack thresh) { slack_threshold_ = thresh; }
  void setWireRC(double r_per_um, double c_per_um)
  {
    wire_r_per_um_ = r_per_um;
    wire_c_per_um_ = c_per_um;
  }
  void setMode(const char* mode_name);
  void setTieLoPort(sta::LibertyPort* loport);
  void setTieHiPort(sta::LibertyPort* hiport);

  // Smart cone selection configuration
  void setConeConfigDelayThreshold(float ratio);
  void setConeConfigMaxDepth(int depth);
  void setConeConfigMaxSize(int size);
  void setConeConfigDistanceThreshold(float distance_um);

 private:
  struct PathNode {
    sta::Vertex* vertex;
    sta::Pin* pin;
    float arrival_time;
    float slack;
    int depth;
    float x;  // X coordinate
    float y;  // Y coordinate
    PathNode* left_child;
    PathNode* right_child;
    bool selected;
    bool is_boundary;
    
    
    PathNode() : vertex(nullptr), pin(nullptr), arrival_time(0.0), 
                 slack(0.0), depth(0), x(0.0), y(0.0),
                 left_child(nullptr), 
                 right_child(nullptr), selected(false), is_boundary(false) {}
  };
  
  struct ConeSelectionConfig {
    float delay_threshold_ratio;
    int max_cone_depth;
    int max_cone_size;
    float min_improvement_threshold;
    float distance_threshold_um;

    ConeSelectionConfig()
      : delay_threshold_ratio(0.3f),
        max_cone_depth(15),
        max_cone_size(500),
        min_improvement_threshold(0.05f),
        distance_threshold_um(50.0f) {}
  };

  // Member variables
  ConeSelectionConfig cone_config_;
  unsigned max_depth_;

  // Smart cone selection methods
  // void selectConeByPathDelay(const ConeSelectionConfig& config);
  PathNode* buildPathTree(const sta::Pin* pin,
                          int current_depth,
                          const ConeSelectionConfig& config);
  bool shouldSelectChild(PathNode* parent,
                         PathNode* left,
                         PathNode* right,
                         const ConeSelectionConfig& config);
  float calculateDelayGap(PathNode* left, PathNode* right);
  float calculatePhysicalDistance(PathNode* node1, PathNode* node2);
  void selectNodesFromTree(PathNode* root, 
                           const ConeSelectionConfig& config,
                           std::set<sta::Vertex*>& selected_vertices,
                           int& current_size);
  void collectSelectedInstances(const std::set<sta::Vertex*>& vertices);
  void cleanupPathTree(PathNode* node);
  float getVertexArrivalTime(sta::Vertex* vertex);
  float getVertexSlack(sta::Vertex* vertex);
  std::pair<float, float> getInstanceCoordinateFromPin(const sta::Pin* pin);

  void deleteComponents();
  void getBlob(unsigned max_depth);
  void runABC();
  void postABC(float worst_slack);
  bool writeAbcScript(const std::string& file_name);
  void writeOptCommands(const std::string& file_name);
  // Writes the setup portion of the ABC script: read_lib + read_blif + strash + read_coords.
  // map and write_blif are executed separately in runABC() so that timing can be
  // injected via C++ API between strash/read_coords and map.
  void writeAbcScriptSetup(const std::string& file_name);
  // abc_frame is the Abc_Frame_t* returned by Abc_FrameGetGlobalFrame().
  // Declared as void* to avoid pulling in abc/main.h (which conflicts with sclLib.h).
  void injectTimingToAbc(void* abc_frame);
  // Extract CI (input) and CO (output) names from the BLIF file written by writeBlif.
  // Returns: {CI_names, CO_names} in the same order as ABC's Abc_NtkCollectCioNames.
  // This is used to build the name mapping between BLIF net names and ABC PI/CO names.
  std::pair<std::vector<std::string>, std::vector<std::string>>
  extractBlifCioNames(const std::string& blif_file);
  // Extract CI/CO names from BLIF for logging only.
  // Actual timing in or_timing_ was already populated in Restructure::run()
  // from blif_.getArrivals() / blif_.getRequireds() before writeBlif destroyed path_insts_.
  void precomputeBoundaryTimingForBlif(const std::string& blif_file);
  void writeNetCoordinates(const std::string& file_name);
  void initDB();
  void getEndPoints(sta::PinSet& ends, bool area_mode, unsigned max_depth);
  int countConsts(odb::dbBlock* top_block);
  void removeConstCells();
  void removeConstCell(odb::dbInst* inst);
  bool readAbcLog(const std::string& abc_file_name,
                  int& level_gain,
                  float& delay_val);

  Logger* logger_;
  utl::UniqueName name_generator_;
  std::string logfile_;
  std::string locell_;
  std::string loport_;
  std::string hicell_;
  std::string hiport_;
  std::string work_dir_name_;

  // db vars
  sta::dbSta* open_sta_;
  odb::dbDatabase* db_;
  rsz::Resizer* resizer_;
  est::EstimateParasitics* estimate_parasitics_;
  odb::dbBlock* block_ = nullptr;

  // Annealing
  std::optional<std::mt19937::result_type> annealing_seed_;
  std::optional<float> annealing_temp_;
  unsigned annealing_iters_ = 100;
  std::optional<unsigned> annealing_revert_after_;
  unsigned annealing_init_ops_ = 10;
  sta::Slack slack_threshold_ = 0;

  std::string input_blif_file_name_;
  std::string output_blif_file_name_;
  std::string coord_file_name_;
  std::vector<std::string> lib_file_names_;
  std::set<odb::dbInst*> path_insts_;

  // Wire RC passed directly from OpenROAD to ABC (ohm/um and fF/um)
  double wire_r_per_um_ = 0.0;
  double wire_c_per_um_ = 0.0;

  // Pre-computed timing at cone boundary nets: net_name → {arrival_ns, required_ns}
  // Populated in Restructure::run() from blif_.getArrivals() / blif_.getRequireds()
  // BEFORE writeBlif destroys path_insts_. Keys are BLIF net names (== ABC PI/CO names after strash).
  std::map<std::string, std::pair<double, double>> or_timing_;

  Mode opt_mode_{Mode::DELAY_1};
  bool is_area_mode_{false};
  int blif_call_id_{0};
};

}  // namespace rmp
