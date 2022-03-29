#ifndef CONNCET_CLASS_HPP
#define CONNCET_CLASS_HPP

class t2l_param {
public :
  double intercept_ax, intercept_ay, intercept_ar, intercept_al;
  double slope_ax, slope_ay, slope_ar, slope_al;
  double slope2_ax, slope2_ay, slope2_ar, slope2_al;
  double intercept_px, intercept_py, intercept_pr, intercept_pl;
  double slope_px, slope_py, slope_pr, slope_pl;
  double slope2_px, slope2_py, slope2_pr, slope2_pl;

  t2l_param& operator=(const t2l_param& obj);

};

#endif
