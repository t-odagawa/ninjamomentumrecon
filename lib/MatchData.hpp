#ifndef MATCH_DATA_HPP
#define MATCH_DATA_HPP

#include <boost/filesystem.hpp>

#include <string>
#include <vector>

#include <TGraphErrors.h>

namespace fs = boost::filesystem;

const extern fs::path MATCH_DIRNAME;

const extern fs::path SHIFTER_EFFICIENCY_FILENAME;
const extern fs::path TRACKER_EFFICIENCY_FILENAME;

const int SHIFTER_EFFICIENCY_VECTOR_SIZE = 10;
const int TRACKER_EFFICIENCY_VECTOR_SIZE = 10;

const double BUNCH_PILEUP_EFFECT = 0.99688;
const double ST_DAQ_EFFECT = 0.9995;

class MatchData {

private:

  std::vector<double > shifter_efficiency_;
  std::vector<double > tracker_efficiency_;

  TGraphErrors *ge_;

  void ReadShifterEfficiencyData(const fs::path &file_dir_path);
  void ReadTrackerEfficiencyData(const fs::path &file_dir_path);

public:
  explicit MatchData(const std::string &file_dir_path);
  
  explicit MatchData(const fs::path &file_dir_path);

  void GetShifterEfficiencyData(std::vector<double > &shifter_efficiency) const;
  void GetTrackerEfficiencyData(std::vector<double > &tracker_efficiency) const;

};

#endif
