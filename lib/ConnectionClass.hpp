#ifndef CONNCET_CLASS_HPP
#define CONNCET_CLASS_HPP

#include <cstddef>

#include <boost/functional/hash.hpp>

class Point {
public :
  double x, y, z;
};

class FiducialArea {
public :
  int pl;
  Point p[2];
};

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

class partner_param {
public :
  double angle_accuracy_intercept_mu;
  double angle_accuracy_slope_mu;
  double angle_accuracy_intercept_partner;
  double angle_accuracy_slope_partner;
};

class Segment {
public : 
  int plate;
  unsigned int rawid;
  bool operator==(const Segment &rhs) const;
  bool operator<(const Segment &rhs) const;

  friend std::size_t hash_value(const Segment &obj) {
    std::size_t h = 0;
    boost::hash_combine(h, obj.plate);
    boost::hash_combine(h, obj.rawid);
    return h;
  }
 
};

class Group {
public :
  int start_plate;
  unsigned int start_rawid;
  std::vector<std::pair<Segment, Segment > > linklets;
};

#endif
