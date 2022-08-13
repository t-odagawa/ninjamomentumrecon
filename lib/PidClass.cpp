#include <iostream>
#include <cmath>

#include "PidClass.hpp"

namespace Pid_data_ns
{

  DataPoint::DataPoint() {
    input_ang_min = -1.;
    input_ang_max = -1.;
    input_mom_min = -1.;
    input_mom_max = -1.;
    for ( int i = 0; i < 3; i++ ) {
      mean[i] = -1.; mean_low[i] = -1.; mean_high[i] = -1.;
      sigma[i] = -1.; sigma_low[i] = -1.; sigma_high[i] = -1.;
    }
  }

  std::ostream &operator<<(std::ostream &os, DataPoint &data) {
    os << "Angle Minimum : " << data.input_ang_min << ", "
       << "Angle Maximum : " << data.input_ang_max << ", "
       << "Momentum minimum : " << data.input_mom_min << ", "
       << "Momentum maximum : " << data.input_mom_max << ", "
       << "Electron VPH mean : " << data.mean[0]
       << data.mean_low[0]
       << "+" << data.mean_high[0] << ", "
       << "Pion/muon VPH mean : " << data.mean[1]
       << data.mean_low[1]
       << "+" << data.mean_high[1] << ", "
       << "Proton VPH mean : " << data.mean[2]
       << data.mean_low[2]
       << "+" << data.mean_high[2] << ", "
       << "Electron VPH sigma : " << data.sigma[0]
       << data.sigma_low[0]
       << "+" << data.sigma_high[0] << ", "
       << "Pion/muon VPH sigma : " << data.sigma[1]
       << data.sigma_low[1]
       << "+" << data.sigma_high[1] << ", "
       << "Proton VPH sigma : " << data.sigma[2]
       << data.sigma_low[2]
       << "+" << data.sigma_high[2];
    
    return os;
    
  }

  std::istream &operator>>(std::istream& is, DataPoint &data) {
    is >> data.input_ang_min >> data.input_ang_max
       >> data.input_mom_min >> data.input_mom_max
       >> data.mean[0] >> data.mean_low[0] >> data.mean_high[0]
       >> data.sigma[0] >> data.sigma_low[0] >> data.sigma_high[0]
       >> data.mean[1] >> data.mean_low[1] >> data.mean_high[1]
       >> data.sigma[1] >> data.sigma_low[1] >> data.sigma_high[1]
       >> data.mean[2] >> data.mean_low[2] >> data.mean_high[2]
       >> data.sigma[2] >> data.sigma_low[2] >> data.sigma_high[2];
      
    return is;
  }

  VphFuncParam::VphFuncParam() {
    input_ang_min = -1.;
    input_ang_max = -1.;
    mean_slope = NAN; mean_slope_err = NAN;
    mean_inter = NAN; mean_inter_err = NAN;
    sigma_inter = NAN; sigma_inter_err = NAN;
    sigma_scale = NAN; sigma_scale_err = NAN;
  }

  std::ostream &operator<<(std::ostream& os, VphFuncParam &param) {
    os << "Angle Minimum : " << param.input_ang_min << ", "
       << "Angle Maximum : " << param.input_ang_max << ", "
       << "Mean Slope : "<< param.mean_slope << " +/- " << param.mean_slope_err << ", "
       << "Mean Intercept : " << param.mean_inter << " +/- " << param.mean_inter_err << ", "
       << "Sigma Slope : " << param.sigma_inter << " +/- " << param.sigma_inter_err << ", "
       << "Sigma Scale : " << param.sigma_scale << " +/- " << param.sigma_scale_err;

    return os;
  }

  std::istream &operator>>(std::istream& is, VphFuncParam &param) {
    is >> param.input_ang_min >> param.input_ang_max
       >> param.mean_slope >> param.mean_slope_err
       >> param.mean_inter >> param.mean_inter_err
       >> param.sigma_inter >> param.sigma_inter_err
       >> param.sigma_scale >> param.sigma_scale_err;

    return is;
  }

  VphPionMip::VphPionMip() {
    ang_min = -1.;
    ang_max = -1.;
    pb_max = 0.;
    expect = 0.;
    bin_id = -1.;
    vph = 0.;
    entry = -1.;
    prob = -1.;
    prob_acc = -1.;
  }

  std::ostream &operator<<(std::ostream& os, VphPionMip &param) {
    os << "Angle Minimum : " << param.ang_min << ", "
       << "Angle Maximum : " << param.ang_max << ", "
       << "Pbeta Maximum : " << param.pb_max << ", "
       << "Expected value : " << param.expect << ", "
       << "Bin ID : " << param.bin_id << ", "
       << "VPH : " << param.vph << ", "
       << "Entry : " << param.entry << ", "
       << "Probability : " << param.prob << ", "
       << "Accumulated probability : " << param.prob_acc;

    return os;
  }

  std::istream &operator>>(std::istream& is, VphPionMip &param) {
    is >> param.ang_min >> param.ang_max
       >> param.pb_max >> param.expect
       >> param.bin_id >> param.vph >> param.entry
       >> param.prob >> param.prob_acc;

    return is;
  }

}
