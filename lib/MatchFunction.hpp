#ifndef MATCH_FUNCTION_HPP
#define MATCH_FUNCTION_HPP

#include <vector>

#include <TVector3.h>
#include <TSpline.h>

#include "MatchData.hpp"

class MatchFunction {

private:
  MatchData match_data_;

  const double conversion_factor_h2o_CSDA_to_iron_CSDA = 1.289;
  
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
