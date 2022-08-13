#ifndef PID_CLASS_HPP
#define PID_CLASS_HPP

#include <iostream>
#include <string>
#include <vector>

namespace Pid_data_ns
{

  class DataPoint {
  public :
    double input_ang_min, input_ang_max, input_mom_min, input_mom_max;
    double mean[3], mean_low[3], mean_high[3];
    double sigma[3], sigma_low[3], sigma_high[3];
    DataPoint();

    friend std::ostream &operator<<(std::ostream& os, DataPoint &data);
    friend std::istream &operator>>(std::istream& is, DataPoint &data);
    
  };

  class VphFuncParam {
  public :
    double input_ang_min, input_ang_max;
    double mean_slope, mean_slope_err, mean_inter, mean_inter_err;
    double sigma_inter, sigma_inter_err, sigma_scale, sigma_scale_err;
    VphFuncParam();
    friend std::ostream &operator<<(std::ostream& os, VphFuncParam &param);
    friend std::istream &operator>>(std::istream& is, VphFuncParam &param);
    
  };

  class VphPionMip {
  public :
    double ang_min, ang_max, pb_max, expect, bin_id, vph, entry, prob, prob_acc;
    VphPionMip();
    friend std::ostream &operator<<(std::ostream& os, VphPionMip &param);
    friend std::istream &operator>>(std::istream& is, VphPionMip &param);
  };

}

#endif
