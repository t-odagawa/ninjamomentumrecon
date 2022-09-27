#ifndef MATCH_FUNCTION_HPP
#define MATCH_FUNCTION_HPP

#include <vector>
#include <cmath>

#include <TVector3.h>
#include <TSpline.h>

#include <NTBMSummary.hh>

#include "MatchData.hpp"

class MatchFunction {

private:
  MatchData match_data_;

  const double conversion_factor_h2o_CSDA_to_iron_CSDA = 1.289;

  static constexpr double pos_angle_border = 0.4;
  static constexpr double pos_x_0 = 2.5;
  static constexpr double pos_x_1 = 6.5;
  static constexpr double pos_const_x = std::hypot(pos_x_0, pos_x_1 * pos_angle_border);
  static constexpr double pos_y_0 = 2.5;
  static constexpr double pos_y_1 = 6.5;
  static constexpr double pos_const_y = std::hypot(pos_y_0, pos_y_1 * pos_angle_border);
  static constexpr double ang_x_0 = 0.040;
  static constexpr double ang_x_1 = 0.050;
  static constexpr double ang_y_0 = 0.015;
  static constexpr double ang_y_1 = 0.050;
  
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

  bool IsMatchBMNTCandidate(int ecc, std::vector<double > position) const;

  double GetChisCutValue(const NTBMSummary *ntbm, int icluster) const;

  double GetSigmaPosX(double tangent) const;
  double GetSigmaPosY(double tangent) const;
  double GetSigmaAngX(double tangent) const;
  double GetSigmaAngY(double tangent) const;

  double CalculateShifterTrackerChi2(const NTBMSummary *ntbm, int icluster, const B2EmulsionSummary* emulsion) const;

};

#endif
