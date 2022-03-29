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

  TSpline3 proton_iron_energy_range_spline_;
  TSpline3 proton_iron_range_energy_spline_;
  TSpline3 proton_water_energy_range_spline_;
  TSpline3 proton_water_range_energy_spline_;
  TSpline3 proton_poly_energy_range_spline_;
  TSpline3 proton_poly_range_energy_spline_;
  TSpline3 proton_emulsion_energy_range_spline_;
  TSpline3 proton_emulsion_range_energy_spline_;

  TSpline3 pion_iron_energy_range_spline_;
  TSpline3 pion_iron_range_energy_spline_;
  TSpline3 pion_water_energy_range_spline_;
  TSpline3 pion_water_range_energy_spline_;
  TSpline3 pion_poly_energy_range_spline_;
  TSpline3 pion_poly_range_energy_spline_;
  TSpline3 pion_emulsion_energy_range_spline_;
  TSpline3 pion_emulsion_range_energy_spline_;

public:
  explicit RangeFunction(const RangeSpline &range_spline);

  double IronRangeFromEnergy(double energy, int particle_id) const;
  double IronEnergyFromRange(double range, int particle_id) const;
  double WaterRangeFromEnergy(double energy, int particle_id) const;
  double WaterEnergyFromRange(double range, int particle_id) const;
  double PolystyreneRangeFromEnergy(double energy, int particle_id) const;
  double PolystyreneEnergyFromRange(double range, int particle_id) const;
  double EmulsionRangeFromEnergy(double energy, int paritlce_id) const;
  double EmulsionEnergyFromRange(double range, int particle_id) const;
  
  void ModifyVectors(std::vector<double> &ax, std::vector<double> &ay, std::vector<int> &pl) const;
  double CalculateEnergyFromRange(std::vector<double> ax, std::vector<double> ay, std::vector<int> pl, int particle_id, int direction) const;
  double CalculateEnergyFromRangeForward(std::vector<double> ax, std::vector<double> ay, std::vector<int> pl, int particle_id) const;
  double CalculateEnergyFromRangeBackward(std::vector<double> ax, std::vector<double> ay, std::vector<int> pl, int particle_id) const;
};

#endif
