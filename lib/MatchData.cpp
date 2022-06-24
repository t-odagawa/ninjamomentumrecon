#include "MatchData.hpp"
#include <boost/log/core.hpp>
#include <boost/log/trivial.hpp>
#include <boost/log/expressions.hpp>
#include <boost/filesystem.hpp>

#include <string>
#include <vector>

#include <TFile.h>
#include <TGraphErrors.h>

const fs::path MATCH_DIRNAME("match");

const fs::path SHIFTER_EFFICIENCY_FILENAME("shifter_efficiency.root");
const fs::path TRACKER_EFFICIENCY_FILENAME("tracker_efficiency.root");

MatchData::MatchData(const std::string &file_dir_path) {
  fs::path file_dir(file_dir_path);
  
  // ReadShifterEfficiencyData(file_dir);
  ReadTrackerEfficiencyData(file_dir);

  BOOST_LOG_TRIVIAL(info) << "Matching data are initalized";
}

MatchData::MatchData(const fs::path &file_dir_path) : MatchData(file_dir_path.string()) {}

void MatchData::ReadShifterEfficiencyData(const fs::path &file_dir_path) {
  return;
}

void MatchData::ReadTrackerEfficiencyData(const fs::path &file_dir_path) {
  const fs::path filepath(file_dir_path/TRACKER_EFFICIENCY_FILENAME);

  if ( !fs::exists(filepath) )
    throw std::runtime_error("Tracker efficiency file : " + filepath.string() + " nor found");

  TFile *file = new TFile(filepath.string().c_str(), "read");
  TGraphErrors *ge = (TGraphErrors*)file->Get("ge");

  for ( int i = 0; i < ge->GetN(); i++ ) {   
    tracker_efficiency_.push_back(ge->GetPointY(i));
  }

  for ( int i = 0; i < tracker_efficiency_.size(); i++ ) {
    BOOST_LOG_TRIVIAL(trace) << "Efficiency : bin" << i
			     << ", value : " << tracker_efficiency_.at(i);
  }

  return;
}

void MatchData::GetShifterEfficiencyData(std::vector<double > &shifter_efficiency) const {
  shifter_efficiency = shifter_efficiency_;
  return;
}

void MatchData::GetTrackerEfficiencyData(std::vector<double > &tracker_efficiency) const {
  tracker_efficiency = tracker_efficiency_;
  return;
}
