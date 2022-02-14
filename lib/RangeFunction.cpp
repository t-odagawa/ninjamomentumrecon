#include <boost/log/core.hpp>
#include <boost/log/trivial.hpp>
#include <boost/log/expressions.hpp>
#include <boost/filesystem.hpp>

#include <iostream>
#include <string>

#include <vector>
#include <cmath>
#include <algorithm>

#include "RangeFunction.hpp"

namespace fs = boost::filesystem;

RangeFunction::RangeFunction(const RangeSpline &range_spline) : range_spline_(range_spline) {

  range_spline_.GetIronEnergyRangeSpline(iron_energy_range_spline_);
  range_spline_.GetWaterEnergyRangeSpline(water_energy_range_spline_);
  range_spline_.GetPolyEnergyRangeSpline(poly_energy_range_spline_);
  range_spline_.GetEmulsionEnergyRangeSpline(emulsion_energy_range_spline_);
  range_spline_.GetIronRangeEnergySpline(iron_range_energy_spline_);
  range_spline_.GetWaterRangeEnergySpline(water_range_energy_spline_);
  range_spline_.GetPolyRangeEnergySpline(poly_range_energy_spline_);
  range_spline_.GetEmulsionRangeEnergySpline(emulsion_range_energy_spline_);

  BOOST_LOG_TRIVIAL(info) << "Range functions are initialized";
}

double RangeFunction::IronRangeFromEnergy(double energy) const {
  return iron_energy_range_spline_.Eval(energy);
}

double RangeFunction::IronEnergyFromRange(double range) const {
  return iron_range_energy_spline_.Eval(range);
}

double RangeFunction::WaterRangeFromEnergy(double energy) const {
  return water_energy_range_spline_.Eval(energy);
}

double RangeFunction::WaterEnergyFromRange(double range) const {
  return water_range_energy_spline_.Eval(range);
}

double RangeFunction::PolystyreneRangeFromEnergy(double energy) const {
  return poly_energy_range_spline_.Eval(energy);
}

double RangeFunction::PolystyreneEnergyFromRange(double range) const {
  return poly_range_energy_spline_.Eval(range);
}

double RangeFunction::EmulsionRangeFromEnergy(double energy) const {
  return emulsion_energy_range_spline_.Eval(energy);
}

double RangeFunction::EmulsionEnergyFromRange(double range) const {
  return emulsion_range_energy_spline_.Eval(range);
}

void RangeFunction::ModifyVectors(std::vector<double> &ax, std::vector<double> &ay, std::vector<int> &pl) const {

  for ( int ipl = 0; ipl < pl.size() - 1; ipl++ ) {
    int pl_difference = pl.at(ipl + 1) - pl.at(ipl);
    if ( pl_difference > 1 ) {
      for (int jpl = 1; jpl < pl_difference; jpl++ ) {
	ax.insert(ax.begin() + ipl + jpl, ax.at(ipl));
	ay.insert(ay.begin() + ipl + jpl, ay.at(ipl));
	pl.insert(pl.begin() + ipl + jpl, pl.at(ipl));
      }
    }
  }

  std::reverse(ax.begin(), ax.end());
  std::reverse(ay.begin(), ay.end());
  std::reverse(pl.begin(), pl.end());

}

double RangeFunction::CalculateEnergyFromRange(std::vector<double> ax, std::vector<double> ay, std::vector<int> pl) const {

  if ( ax.size() != ay.size() ||
       ax.size() != pl.size() ||
       ay.size() != pl.size() ) {
    BOOST_LOG_TRIVIAL(error) << "Size of input vectors should be the same";
    std::exit(1);
  }

  double energy = 0.;
  double scale_factor = 1.;
  double range = 0.;
  double range_tmp = 0.;

  for ( int ipl = 0; ipl < pl.size(); ipl++ ) {

    scale_factor = std::sqrt(ax.at(ipl) * ax.at(ipl) + ay.at(ipl) * ay.at(ipl) + 1.);

    if ( pl.at(ipl) == 3 ) { // ISS downstream
      range  = 1. * scale_factor;
      energy = PolystyreneEnergyFromRange(range);      
      range  = 70.e-3 * scale_factor + EmulsionRangeFromEnergy(energy);
      energy = EmulsionEnergyFromRange(range);
      range  = 210.e-3 * scale_factor + PolystyreneRangeFromEnergy(energy);
      energy = PolystyreneEnergyFromRange(range);
      range  = 70.e-3 * scale_factor + EmulsionRangeFromEnergy(energy);
      energy = EmulsionEnergyFromRange(range);
    }
    else if ( pl.at(ipl) == 4 ) { // ISS upstream
      if ( ipl == 0 ) { // if stopping plate
	range  = 210.e-3 * scale_factor * 0.5;
	energy = PolystyreneEnergyFromRange(range);
	range  = 70.e-3 * scale_factor + EmulsionRangeFromEnergy(energy);
      }
      range  = 70.e-3 * scale_factor + EmulsionRangeFromEnergy(energy);
      energy = EmulsionEnergyFromRange(range);
      range  = 210.e-3 * scale_factor + PolystyreneRangeFromEnergy(energy);
      energy = PolystyreneEnergyFromRange(range);
      range  = 70.e-3 * scale_factor + EmulsionRangeFromEnergy(energy);
      energy = EmulsionEnergyFromRange(range);
    }
    else if ( pl.at(ipl) < 15 ) { // Fe ECC
      if ( ipl == 0 )  // if stopping plate
	range = 500.e-3 * scale_factor * 0.5;
      else
	range = 500.e-3 * scale_factor + IronRangeFromEnergy(energy);
      energy = IronEnergyFromRange(range);
      range  = 70.e-3 * scale_factor + EmulsionRangeFromEnergy(energy);
      energy = EmulsionEnergyFromRange(range);
      range  = 210.e-3 * scale_factor + PolystyreneRangeFromEnergy(energy);
      energy = PolystyreneEnergyFromRange(range);
      range  = 70.e-3 * scale_factor + EmulsionRangeFromEnergy(energy);
      energy = EmulsionEnergyFromRange(range);
    }
    else if ( pl.at(ipl) == 16 ) { // most downstream of water ECC
      if ( ipl == 0 ) { // if stopping plate
	range  = 210.e-3 * scale_factor * 0.5;
	energy = PolystyreneEnergyFromRange(range);
	range  = 70.e-3 * scale_factor + EmulsionRangeFromEnergy(energy);	
      }
      range  = 70.e-3 * scale_factor + EmulsionRangeFromEnergy(energy);
      energy = EmulsionEnergyFromRange(range);
      range  = 210.e-3 * scale_factor + PolystyreneRangeFromEnergy(energy);
      energy = PolystyreneEnergyFromRange(range);
      range  = 70.e-3 * scale_factor + EmulsionRangeFromEnergy(energy);
      energy = EmulsionEnergyFromRange(range);
    }
    else if ( pl.at(ipl) % 2 == 1 ) { // upstream of iron
      if ( ipl == 0 )  // if stopping plate
	range = 500.e-3 * scale_factor * 0.5;
      else
	range = 500.e-3 * scale_factor + IronRangeFromEnergy(energy);
      energy = IronEnergyFromRange(range);
      range  = 70.e-3 * scale_factor + EmulsionRangeFromEnergy(energy);
      energy = EmulsionEnergyFromRange(range);
      range  = 210.e-3 * scale_factor + PolystyreneRangeFromEnergy(energy);
      energy = PolystyreneEnergyFromRange(range);
      range  = 70.e-3 * scale_factor + EmulsionRangeFromEnergy(energy);
      energy = EmulsionEnergyFromRange(range);
      range  = 109.e-3 * scale_factor + PolystyreneRangeFromEnergy(energy);
      energy = PolystyreneEnergyFromRange(range);
    }
    else { // upstream of water
      if ( ipl == 0 )  // if stopping plate
	range = 2.3 * scale_factor * 0.5;
      else
	range = 2.3 * scale_factor + WaterRangeFromEnergy(energy);
      energy = WaterEnergyFromRange(range);
      range = 109.e-3 * scale_factor + PolystyreneRangeFromEnergy(energy);
      energy = PolystyreneEnergyFromRange(range);
      range  = 70.e-3 * scale_factor + EmulsionRangeFromEnergy(energy);
      energy = EmulsionEnergyFromRange(range);
      range  = 210.e-3 * scale_factor + PolystyreneRangeFromEnergy(energy);
      energy = PolystyreneEnergyFromRange(range);
      range  = 70.e-3 * scale_factor + EmulsionRangeFromEnergy(energy);
      energy = EmulsionEnergyFromRange(range);
    }

  }
  
  return energy;
}
