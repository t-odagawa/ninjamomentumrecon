#ifndef MATCH_FUNCTION_HPP
#define MATCH_FUNCTION_HPP

#include <vector>

#include "MatchData.hpp"

class MatchFunction {

private:
  MatchData match_data_;
  
  // efficiencies
  std::vector<double > shifter_efficiency_;
  std::vector<double > tracker_efficiency_;

public:
  explicit MatchFunction(const MatchData &match_data);

  double GetTrackerEfficiency(double bm_angle) const;

  double GetBmErr(double range_mom) const;

};

#endif
