#include <boost/log/core.hpp>
#include <boost/log/trivial.hpp>
#include <boost/log/expressions.hpp>

#include <TVector3.h>

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

double MatchFunction::GetBmErr(double range_mom) const {
  if ( range_mom < 300 ) return 0.079020061;
  else if ( range_mom < 400 ) return 0.079020061;
  else if ( range_mom < 500 ) return 0.080035672;
  else if ( range_mom < 600 ) return 0.086198609;
  else if ( range_mom < 700 ) return 0.082791300;
  else if ( range_mom < 800 ) return 0.071492218;
  else if ( range_mom < 900 ) return 0.074543873;
  else if ( range_mom < 1000 ) return 0.072054229;
  else if ( range_mom < 1100 ) return 0.056333879;
  else if ( range_mom < 1200 ) return 0.048178771;
  else if ( range_mom < 1300 ) return 0.037908676;
  else if ( range_mom < 1400 ) return 0.030906878;
  else return 0.036638113;
}

int MatchFunction::GetNumWaterPlate(int vertex_pl) const {
  if ( vertex_pl < 17 ) return 0;
  else return ( vertex_pl - 16 ) / 2;
}

int MatchFunction::GetNumIronPlate(int vertex_pl) const {
  if ( vertex_pl < 16 ) return (vertex_pl - 4);
  else if ( vertex_pl == 16 ) return 11;
  else return (vertex_pl - 15) / 2 + 11;
}

void MatchFunction::ConvertFromLengthToMom(double &range_mom, double track_length) const {
  for ( int i = 0; i < 20; i++ ) {
    double lefthand, righthand;
    double gapx, gapy;
    if ( track_length > x_[i] && track_length < x_[i+1] ) {
      lefthand = std::fabs(track_length - x_[i]);
      righthand = std::fabs(x_[i+1] - track_length);
      gapx = x_[i+1] - x_[i];
      gapy = y_[i+1] - y_[i];
      range_mom = y_[i] + gapy * (lefthand/gapx);
      break;
    }
  }
}

double MatchFunction::GetDWGTrackLength(TVector3 track_direction) const {
  double track_length_ = 0.;
  for ( int i = 0; i < 8; i++ ) {
    track_length_ += 60.;
  }
  track_length_ += 18.2;
  track_length_ += 24 * 7.874;

  return track_length_ * track_direction.Mag();

}
