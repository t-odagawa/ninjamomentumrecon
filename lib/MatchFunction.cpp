#include <iostream>

#include <boost/log/core.hpp>
#include <boost/log/trivial.hpp>
#include <boost/log/expressions.hpp>

#include <TVector3.h>
#include <TSpline.h>

#include <B2Const.hh>
#include <B2EmulsionSummary.hh>

#include <NTBMSummary.hh>

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

  double x_[21] = {0., 1.025, 1.854, 3.437, 6.812, 10.91, 31.78, 43.82, 69.50, 109.9,
		   178.7, 247.2, 512.4, 640.2, 888.5, 1248, 1825, 2383, 4509, 5532, 7524};
  double y_[21] = {0., 47.04, 56.16, 68.02, 85.09, 100.3, 152.7, 176.4, 221.8, 286.8, 391.7,
		   494.5, 899.5, 1101, 1502, 2103, 3104, 4104, 8105, 10110, 14110};

  TSpline3 range_spline_("range_spline_", x_, y_, 21);

  range_mom = range_spline_.Eval(track_length);
  /*
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
  */
}

double MatchFunction::GetDWGTrackLength(TVector3 track_direction) const {
  double track_length_ = 0.;
  for ( int i = 0; i < 8; i++ ) {
    track_length_ += 6.;
  }
  track_length_ += 1.82;
  track_length_ *= conversion_factor_h2o_CSDA_to_iron_CSDA;
  track_length_ += 2.4 * 7.874;

  return track_length_ * track_direction.Mag();

}

bool MatchFunction::IsMatchBMNTCandidate(int ecc, std::vector<double > position) const {
  if ( ecc != 5 ) return false;

  double center_y = 50.;
  double center_x = -50.;
  double width_y = 500.;
  double width_x = 600.;

  double y = position.at(0) - center_y;
  double x = position.at(1) - center_x;

  if ( std::fabs(y) < width_y / 2. &&
       std::fabs(x) < width_x / 2. ) return true;

  return false;
  
}

double MatchFunction::GetChisCutValue(const NTBMSummary *ntbm, int icluster) const {

  int bm_id = ntbm->GetBabyMindTrackId(icluster);
  int plane = ntbm->GetBabyMindMaximumPlane(bm_id);

  if ( plane >= 12 ) return 30.;
  else return 30. + 3. * (12 - plane);

}

double MatchFunction::GetSigmaPosX(double tangent) const {
  if ( std::fabs(tangent) < pos_angle_border )
    return std::hypot(pos_x_0, pos_x_1 * tangent);
  else
    return pos_const_x;
}

double MatchFunction::GetSigmaPosY(double tangent) const {
  if ( std::fabs(tangent) < pos_angle_border )
    return std::hypot(pos_y_0, pos_y_1 * tangent);
  else
    return pos_const_y;
}

double MatchFunction::GetSigmaAngX(double tangent) const {
  return std::hypot(ang_x_0, ang_x_1 * tangent);
}

double MatchFunction::GetSigmaAngY(double tangent) const {
  return std::hypot(ang_y_0, ang_y_1 * tangent);
}

double MatchFunction::CalculateShifterTrackerChi2(const NTBMSummary* ntbm, int icluster, const B2EmulsionSummary* emulsion) const {

  auto nt_pos = ntbm->GetNinjaPosition(icluster);
  auto nt_ang = ntbm->GetNinjaTangent(icluster);

  TVector3 bt_pos = emulsion->GetAbsolutePosition().GetValue();
  TVector3 bt_ang = emulsion->GetTangent().GetValue();

  // TSS position を tracker 座標系に変換
  bt_pos.SetX(bt_pos.X()
	      - NINJA_POS_X // NINJA box
	      - NINJA_TRACKER_POS_X // tracker box
	      );
  bt_pos.SetY(bt_pos.Y()
	      - NINJA_POS_Y
	      - NINJA_TRACKER_POS_Y
	      );

  // 外挿
  double tss_nt_distance_x = NINJA_TRACKER_POS_Z + 10. // tracker z position
    - NINJA_TSS_POS_Z  // TSS acryl center
    + 0.5 * NINJA_TSS_ATTACH_AC_THICK // half acryl
    + NINJA_ENV_THICK // envelope
    + 2 * NINJA_FILM_THICK // two films
    + NINJA_TSS_AC_THICK // SS acryl
    + NINJA_EMULSION_LAYER_THICK // emulsion gel
    + NINJA_BASE_LAYER_THICK; // base
  double tss_nt_distance_y = NINJA_TRACKER_POS_Z - 10. // tracker z position
    - NINJA_TSS_POS_Z  // TSS acryl center
    + 0.5 * NINJA_TSS_ATTACH_AC_THICK // half acryl
    + NINJA_ENV_THICK // envelope
    + 2 * NINJA_FILM_THICK // two films
    + NINJA_TSS_AC_THICK // SS acryl
    + NINJA_EMULSION_LAYER_THICK // emulsion gel
    + NINJA_BASE_LAYER_THICK; // base

  bt_pos.SetX(bt_pos.X() + bt_ang.X() * tss_nt_distance_x);
  bt_pos.SetY(bt_pos.Y() + bt_ang.Y() * tss_nt_distance_y);

  double dx = nt_pos.at(B2View::kTopView) - bt_pos.X();
  double dy = nt_pos.at(B2View::kSideView) - bt_pos.Y();
  double dax = nt_ang.at(B2View::kTopView) - bt_ang.X();
  double day = nt_ang.at(B2View::kSideView) - bt_ang.Y();

  BOOST_LOG_TRIVIAL(debug) << "Tracker : X : " << nt_pos.at(B2View::kTopView) << ", "
			   << "Y : " << nt_pos.at(B2View::kSideView) << ", "
			   << "aX : " << nt_ang.at(B2View::kTopView) << ", "
			   << "aY : " << nt_ang.at(B2View::kSideView);
  BOOST_LOG_TRIVIAL(debug) << "Shifter : X : " << bt_pos.X() << ", "
			   << "Y : " << bt_pos.Y() << ", "
			   << "aX : " << bt_ang.X() << ", "
			   << "aY : " << bt_ang.Y();
  BOOST_LOG_TRIVIAL(debug) << "dx = " << dx << ", "
			   << "dy = " << dy << ", "
			   << "dax = " << dax << ", "
			   << "day = " << day;

  double sigma_x = GetSigmaPosX(bt_ang.X());
  double sigma_y = GetSigmaPosY(bt_ang.Y());
  double sigma_ax = GetSigmaAngX(bt_ang.X());
  double sigma_ay = GetSigmaAngY(bt_ang.Y());

  double chis_dx = (dx * dx) / (sigma_x * sigma_x);
  double chis_dy = (dy * dy) / (sigma_y * sigma_y);
  double chis_dax = (dax * dax) / (sigma_ax * sigma_ax);
  double chis_day = (day * day) / (sigma_ay * sigma_ay);

  return chis_dx + chis_dy + chis_dax + chis_day;

}
