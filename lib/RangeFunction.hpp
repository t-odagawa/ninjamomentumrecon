#ifndef RANGE_FUNCTION_HPP
#define RANGE_FUNCTION_HPP

#include <boost/filesystem.hpp>

#include <vector>
#include <string>

#include "RangeSpline.hpp"

namespace fs = boost::filesystem;

class RangeFunction {

private:
  RangeSpline range_spline_;
  TSpline3 iron_energy_range_spline_;
  TSpline3 iron_range_energy_spline_;
  TSpline3 water_energy_range_spline_;
  TSpline3 water_range_energy_spline_;
  TSpline3 poly_energy_range_spline_;
  TSpline3 poly_range_energy_spline_;
  TSpline3 emulsion_energy_range_spline_;
  TSpline3 emulsion_range_energy_spline_;

public:
  explicit RangeFunction(const RangeSpline &range_spline);

  double IronRangeFromEnergy(double energy) const;
  double IronEnergyFromRange(double range) const;
  double WaterRangeFromEnergy(double energy) const;
  double WaterEnergyFromRange(double range) const;
  double PolystyreneRangeFromEnergy(double energy) const;
  double PolystyreneEnergyFromRange(double range) const;
  double EmulsionRangeFromEnergy(double energy) const;
  double EmulsionEnergyFromRange(double range) const;
  
  void ModifyVectors(std::vector<double> &ax, std::vector<double> &ay, std::vector<int> &pl) const;
  
  double CalculateEnergyFromRange(std::vector<double> ax, std::vector<double> ay, std::vector<int> pl) const;
};

#endif
