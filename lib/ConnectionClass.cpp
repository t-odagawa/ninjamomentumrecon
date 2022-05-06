#include "ConnectionClass.hpp"

#include <boost/unordered_map.hpp>
#include <cstddef>
#include <set>

t2l_param& t2l_param::operator=(const t2l_param& obj) {
  intercept_ax = obj.intercept_ax;
  intercept_ay = obj.intercept_ay;
  intercept_ar = obj.intercept_ar;
  intercept_al = obj.intercept_al;
  slope_ax = obj.slope_ax;
  slope_ay = obj.slope_ay;
  slope_ar = obj.slope_ar;
  slope_al = obj.slope_al;
  slope2_ax = obj.slope2_ax;
  slope2_ay = obj.slope2_ay;
  slope2_ar = obj.slope2_ar;
  slope2_al = obj.slope2_al;

  intercept_px = obj.intercept_px;
  intercept_py = obj.intercept_py;
  intercept_pr = obj.intercept_pr;
  intercept_pl = obj.intercept_pl;
  slope_px = obj.slope_px;
  slope_py = obj.slope_py;
  slope_pr = obj.slope_pr;
  slope_pl = obj.slope_pl;
  slope2_px = obj.slope2_px;
  slope2_py = obj.slope2_py;
  slope2_pr = obj.slope2_pr;
  slope2_pl = obj.slope2_pl;

  return *this;
}

bool Segment::operator==(const Segment &rhs) const {
  return plate == rhs.plate && rawid == rhs.rawid;
}

bool Segment::operator<(const Segment &rhs) const {
  if ( plate == rhs.plate )
    return rawid < rhs.rawid;
  return plate < rhs.plate;
}

bool Group::operator==(const Group &rhs) const {
  std::set<std::pair<Segment, Segment > > group_set;
  std::set<std::pair<Segment, Segment > > group_set_rhs;
  for ( auto linklet : linklets ) {
    group_set.insert(linklet);
  }
  for ( auto rhs_linklet : rhs.linklets ) {
    group_set_rhs.insert(rhs_linklet);
  }
  return group_set == group_set_rhs;

}
