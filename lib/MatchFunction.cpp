#include <boost/log/core.hpp>
#include <boost/log/trivial.hpp>
#include <boost/log/expressions.hpp>

#include "MatchData.hpp"
#include "MatchFunction.hpp"

MatchFunction::MatchFunction(const MatchData &match_data) : match_data_(match_data) {
  
  //match_data_.GetShifterEfficiencyData(shifter_efficiency_);
  match_data_.GetTrackerEfficiencyData(tracker_efficiency_);

  BOOST_LOG_TRIVIAL(info) << "Match functions are initialized";
}

double MatchFunction::GetTrackerEfficiency(double bm_angle) const {
  if ( bm_angle < 0. )       return 0.; // not valid
  else if ( bm_angle < 40. ) return tracker_efficiency_.at((int)(bm_angle / 5));
  else return tracker_efficiency_.at(8);
}
