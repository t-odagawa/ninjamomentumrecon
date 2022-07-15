#ifndef MATCH_FUNCTION_HPP
#define MATCH_FUNCTION_HPP

#include <vector>

#include <TVector3.h>

#include "MatchData.hpp"

class MatchFunction {

private:
  MatchData match_data_;
  
  const double x_[21] = {0., 0.1025, 1.854, 3.437, 6.812, 10.91, 31.78, 43.82, 69.50, 109.9,
			 178.7, 247.2, 512.4, 640.2, 888.5, 1248, 1825, 2383, 4509, 5532, 7524};
  const double y_[21] = {0., 47.04, 56.16, 68.02, 85.09, 100.3, 152.7, 176.4, 221.8, 286.8, 391.7,
			 494.5, 899.5, 1101, 1502, 3104, 4104, 8105, 10110, 14110};

  // efficiencies
  std::vector<double > shifter_efficiency_;
  std::vector<double > tracker_efficiency_;

public:
  explicit MatchFunction(const MatchData &match_data);

  double GetTrackerEfficiency(double bm_angle) const;

  double GetBmErr(double range_mom) const;

  int GetNumWaterPlate(int vertex_pl) const;
  int GetNumIronPlate(int vertex_pl) const;
  void ConvertFromLengthToMom(double &range_mom, double track_length) const;

  double GetDWGTrackLength(TVector3 track_direction) const;
};

#endif
