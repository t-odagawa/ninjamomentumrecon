#include <boost/log/core.hpp>
#include <boost/log/trivial.hpp>
#include <boost/log/expressions.hpp>
#include <boost/filesystem.hpp>

#include <iostream>
#include <string>

#include <vector>
#include <cmath>
#include <algorithm>

#include <B2Pdg.hh>

#include "RangeFunction.hpp"

namespace fs = boost::filesystem;

RangeFunction::RangeFunction(const RangeSpline &range_spline) : range_spline_(range_spline) {

  range_spline_.GetProtonIronEnergyRangeSpline(proton_iron_energy_range_spline_);
  range_spline_.GetProtonWaterEnergyRangeSpline(proton_water_energy_range_spline_);
  range_spline_.GetProtonPolyEnergyRangeSpline(proton_poly_energy_range_spline_);
  range_spline_.GetProtonEmulsionEnergyRangeSpline(proton_emulsion_energy_range_spline_);
  range_spline_.GetProtonIronRangeEnergySpline(proton_iron_range_energy_spline_);
  range_spline_.GetProtonWaterRangeEnergySpline(proton_water_range_energy_spline_);
  range_spline_.GetProtonPolyRangeEnergySpline(proton_poly_range_energy_spline_);
  range_spline_.GetProtonEmulsionRangeEnergySpline(proton_emulsion_range_energy_spline_);

  range_spline_.GetPionIronEnergyRangeSpline(pion_iron_energy_range_spline_);
  range_spline_.GetPionWaterEnergyRangeSpline(pion_water_energy_range_spline_);
  range_spline_.GetPionPolyEnergyRangeSpline(pion_poly_energy_range_spline_);
  range_spline_.GetPionEmulsionEnergyRangeSpline(pion_emulsion_energy_range_spline_);
  range_spline_.GetPionIronRangeEnergySpline(pion_iron_range_energy_spline_);
  range_spline_.GetPionWaterRangeEnergySpline(pion_water_range_energy_spline_);
  range_spline_.GetPionPolyRangeEnergySpline(pion_poly_range_energy_spline_);
  range_spline_.GetPionEmulsionRangeEnergySpline(pion_emulsion_range_energy_spline_);

  BOOST_LOG_TRIVIAL(info) << "Range functions are initialized";
}

double RangeFunction::IronRangeFromEnergy(double energy, int particle_id) const {
  switch (particle_id) {
  case PDG_t::kProton :
    return proton_iron_energy_range_spline_.Eval(energy);
    break;
  case PDG_t::kPiPlus :
  case PDG_t::kPiMinus :
    return pion_iron_energy_range_spline_.Eval(energy);
    break;
  default :
    throw std::invalid_argument("Particle id is not in interest");
  }
}

double RangeFunction::IronEnergyFromRange(double range, int particle_id) const {
  switch (particle_id) {
  case PDG_t::kProton :
    return proton_iron_range_energy_spline_.Eval(range);
    break;
  case PDG_t::kPiPlus :
  case PDG_t::kPiMinus :
    return pion_iron_range_energy_spline_.Eval(range);
    break;
  default :
    throw std::invalid_argument("Particle id is not in interest");
  }
}

double RangeFunction::WaterRangeFromEnergy(double energy, int particle_id) const {
  switch (particle_id) {
  case PDG_t::kProton :
    return proton_water_energy_range_spline_.Eval(energy);
    break;
  case PDG_t::kPiPlus :
  case PDG_t::kPiMinus :
    return pion_water_energy_range_spline_.Eval(energy);
    break;
  default :
    throw std::invalid_argument("Particle id is not in interest");
  }
}

double RangeFunction::WaterEnergyFromRange(double range, int particle_id) const {
  switch (particle_id) {
  case PDG_t::kProton :
    return proton_water_range_energy_spline_.Eval(range);
    break;
  case PDG_t::kPiPlus :
  case PDG_t::kPiMinus :
    return pion_water_range_energy_spline_.Eval(range);
    break;
  default :
    throw std::invalid_argument("Particle id is not in interest");
  }
}

double RangeFunction::PolystyreneRangeFromEnergy(double energy, int particle_id) const {
  switch (particle_id) {
  case PDG_t::kProton :
    return proton_poly_energy_range_spline_.Eval(energy);
    break;
  case PDG_t::kPiPlus :
  case PDG_t::kPiMinus :
    return pion_poly_energy_range_spline_.Eval(energy);
    break;
  default :
    throw std::invalid_argument("Particle id is not in interest");
  }
}

double RangeFunction::PolystyreneEnergyFromRange(double range, int particle_id) const {
  switch (particle_id) {
  case PDG_t::kProton :
    return proton_poly_range_energy_spline_.Eval(range);
    break;
  case PDG_t::kPiPlus :
  case PDG_t::kPiMinus :
    return pion_poly_range_energy_spline_.Eval(range);
    break;
  default :
    throw std::invalid_argument("Particle id is not in interest");
  }
}

double RangeFunction::EmulsionRangeFromEnergy(double energy, int particle_id) const {
  switch (particle_id) {
  case PDG_t::kProton :
    return proton_emulsion_energy_range_spline_.Eval(energy);
    break;
  case PDG_t::kPiPlus :
  case PDG_t::kPiMinus :
    return pion_emulsion_energy_range_spline_.Eval(energy);
    break;
  default :
    throw std::invalid_argument("Particle id is not in interest");
  }
}

double RangeFunction::EmulsionEnergyFromRange(double range, int particle_id) const {
  switch (particle_id) {
  case PDG_t::kProton :
    return proton_emulsion_range_energy_spline_.Eval(range);
    break;
  case PDG_t::kPiPlus :
  case PDG_t::kPiMinus :
    return pion_emulsion_range_energy_spline_.Eval(range);
    break;
  default :
    throw std::invalid_argument("Particle id is not in interest");
  }
}

void RangeFunction::ModifyVectors(std::vector<double> &ax, std::vector<double> &ay, std::vector<int> &pl) const {

  for ( int ipl = 0; ipl < pl.size() - 1; ipl++ ) {
    int pl_difference = pl.at(ipl + 1) - pl.at(ipl);
    if ( pl_difference > 1 ) {
      for (int jpl = 1; jpl < pl_difference; jpl++ ) {
	ax.insert(ax.begin() + ipl + jpl, ax.at(ipl));
	ay.insert(ay.begin() + ipl + jpl, ay.at(ipl));
	pl.insert(pl.begin() + ipl + jpl, pl.at(ipl)+jpl);
      }
    }
  }

  std::reverse(ax.begin(), ax.end());
  std::reverse(ay.begin(), ay.end());
  std::reverse(pl.begin(), pl.end());

}

double RangeFunction::CalculateEnergyFromRange(std::vector<double> ax, std::vector<double> ay, std::vector<int> pl, int particle_id, int direction) const {

  if ( ax.size() != ay.size() ||
       ax.size() != pl.size() ||
       ay.size() != pl.size() ) {
    BOOST_LOG_TRIVIAL(error) << "Size of input vectors should be the same";
    std::exit(1);
  }

  switch (direction) {
  case 1 :
    return CalculateEnergyFromRangeForward(ax, ay, pl, particle_id);
    break;
  case -1 :
    return CalculateEnergyFromRangeBackward(ax, ay, pl, particle_id);
    break;
  default :
    throw std::invalid_argument("Direction should be +/- 1");
  }

}

double RangeFunction::CalculateEnergyFromRangeForward(std::vector<double> ax, std::vector<double> ay, std::vector<int> pl, int particle_id) const {

  double energy = 0.;
  double scale_factor = 1.;
  double range = 0.;
  double range_tmp = 0.;

  for ( int ipl = 0; ipl < pl.size(); ipl++ ) {

    scale_factor = std::sqrt(ax.at(ipl) * ax.at(ipl) + ay.at(ipl) * ay.at(ipl) + 1.);

    if ( pl.at(ipl) == 3 ) { // ISS downstream
      range  = 1. * scale_factor;
      energy = PolystyreneEnergyFromRange(range, particle_id);
      range  = 70.e-3 * scale_factor + EmulsionRangeFromEnergy(energy, particle_id);
      energy = EmulsionEnergyFromRange(range, particle_id);
      range  = 210.e-3 * scale_factor + PolystyreneRangeFromEnergy(energy, particle_id);
      energy = PolystyreneEnergyFromRange(range, particle_id);
      range  = 70.e-3 * scale_factor + EmulsionRangeFromEnergy(energy, particle_id);
      energy = EmulsionEnergyFromRange(range, particle_id);
    }
    else if ( pl.at(ipl) == 4 ) { // ISS upstream
      if ( ipl == 0 ) { // if stopping plate
	range  = 210.e-3 * scale_factor * 0.5;
	energy = PolystyreneEnergyFromRange(range, particle_id);
	range  = 70.e-3 * scale_factor + EmulsionRangeFromEnergy(energy, particle_id);
      }
      range  = 70.e-3 * scale_factor + EmulsionRangeFromEnergy(energy, particle_id);
      energy = EmulsionEnergyFromRange(range, particle_id);
      range  = 210.e-3 * scale_factor + PolystyreneRangeFromEnergy(energy, particle_id);
      energy = PolystyreneEnergyFromRange(range, particle_id);
      range  = 70.e-3 * scale_factor + EmulsionRangeFromEnergy(energy, particle_id);
      energy = EmulsionEnergyFromRange(range, particle_id);
    }
    else if ( pl.at(ipl) < 15 ) { // Fe ECC
      if ( ipl == 0 )  // if stopping plate
	range = 500.e-3 * scale_factor * 0.5;
      else
	range = 500.e-3 * scale_factor + IronRangeFromEnergy(energy, particle_id);
      energy = IronEnergyFromRange(range, particle_id);
      range  = 70.e-3 * scale_factor + EmulsionRangeFromEnergy(energy, particle_id);
      energy = EmulsionEnergyFromRange(range, particle_id);
      range  = 210.e-3 * scale_factor + PolystyreneRangeFromEnergy(energy, particle_id);
      energy = PolystyreneEnergyFromRange(range, particle_id);
      range  = 70.e-3 * scale_factor + EmulsionRangeFromEnergy(energy, particle_id);
      energy = EmulsionEnergyFromRange(range, particle_id);
    }
    else if ( pl.at(ipl) == 16 ) { // most downstream of water ECC
      if ( ipl == 0 ) { // if stopping plate
	range  = 210.e-3 * scale_factor * 0.5;
	energy = PolystyreneEnergyFromRange(range, particle_id);
	range  = 70.e-3 * scale_factor + EmulsionRangeFromEnergy(energy, particle_id);
      }
      range  = 70.e-3 * scale_factor + EmulsionRangeFromEnergy(energy, particle_id);
      energy = EmulsionEnergyFromRange(range, particle_id);
      range  = 210.e-3 * scale_factor + PolystyreneRangeFromEnergy(energy, particle_id);
      energy = PolystyreneEnergyFromRange(range, particle_id);
      range  = 70.e-3 * scale_factor + EmulsionRangeFromEnergy(energy, particle_id);
      energy = EmulsionEnergyFromRange(range, particle_id);
    }
    else if ( pl.at(ipl) % 2 == 1 ) { // upstream of iron
      if ( ipl == 0 )  // if stopping plate
	range = 500.e-3 * scale_factor * 0.5;
      else
	range = 500.e-3 * scale_factor + IronRangeFromEnergy(energy, particle_id);
      energy = IronEnergyFromRange(range, particle_id);
      range  = 70.e-3 * scale_factor + EmulsionRangeFromEnergy(energy, particle_id);
      energy = EmulsionEnergyFromRange(range, particle_id);
      range  = 210.e-3 * scale_factor + PolystyreneRangeFromEnergy(energy, particle_id);
      energy = PolystyreneEnergyFromRange(range, particle_id);
      range  = 70.e-3 * scale_factor + EmulsionRangeFromEnergy(energy, particle_id);
      energy = EmulsionEnergyFromRange(range, particle_id);
      range  = 109.e-3 * scale_factor + PolystyreneRangeFromEnergy(energy, particle_id);
      energy = PolystyreneEnergyFromRange(range, particle_id);
    }
    else { // upstream of water
      if ( ipl == 0 )  // if stopping plate
	range = 2.3 * scale_factor * 0.5;
      else
	range = 2.3 * scale_factor + WaterRangeFromEnergy(energy, particle_id);
      energy = WaterEnergyFromRange(range, particle_id);
      range = 109.e-3 * scale_factor + PolystyreneRangeFromEnergy(energy, particle_id);
      energy = PolystyreneEnergyFromRange(range, particle_id);
      range  = 70.e-3 * scale_factor + EmulsionRangeFromEnergy(energy, particle_id);
      energy = EmulsionEnergyFromRange(range, particle_id);
      range  = 210.e-3 * scale_factor + PolystyreneRangeFromEnergy(energy, particle_id);
      energy = PolystyreneEnergyFromRange(range, particle_id);
      range  = 70.e-3 * scale_factor + EmulsionRangeFromEnergy(energy, particle_id);
      energy = EmulsionEnergyFromRange(range, particle_id);
    }

  }
  
  return energy;
}


double RangeFunction::CalculateEnergyFromRangeBackward(std::vector<double> ax, std::vector<double> ay, std::vector<int> pl, int particle_id) const {

  // vectors are filled from downstream to upstream
  std::reverse(ax.begin(), ax.end());
  std::reverse(ay.begin(), ay.end());
  std::reverse(pl.begin(), pl.end());
  
  double energy = 0.;
  double scale_factor = 1.;
  double range = 0.;
  double range_tmp = 0.;

  for ( int ipl = 0; ipl < pl.size(); ipl++ ) {
    
    scale_factor = std::sqrt(ax.at(ipl) * ax.at(ipl) + ay.at(ipl) * ay.at(ipl) + 1.);

    if ( pl.at(ipl) == 133 ) { // most upatream film
      range  = 109.e-3 * scale_factor;
      energy = PolystyreneEnergyFromRange(range, particle_id);
      range  = 70.e-3 * scale_factor + EmulsionRangeFromEnergy(energy, particle_id);
      energy = EmulsionEnergyFromRange(range, particle_id);
      range  = 210.e-3 * scale_factor + PolystyreneRangeFromEnergy(energy, particle_id);
      energy = PolystyreneEnergyFromRange(range, particle_id);
      range  = 70.e-3 * scale_factor + EmulsionRangeFromEnergy(energy, particle_id);
      energy = EmulsionEnergyFromRange(range, particle_id);
    }
    else if ( pl.at(ipl) >= 16 ) {
      if ( pl.at(ipl) %2 == 1 ) { // upstream of iron
	if ( ipl == 0 )  // if stopping plate
	  range = 2.3 * scale_factor * 0.5;
	else
	  range = 2.3 * scale_factor + WaterRangeFromEnergy(energy, particle_id);
	energy = WaterEnergyFromRange(range, particle_id);
	range = 109.e-3 * scale_factor + PolystyreneRangeFromEnergy(energy, particle_id);
	energy = PolystyreneEnergyFromRange(range, particle_id);
	range  = 70.e-3 * scale_factor + EmulsionRangeFromEnergy(energy, particle_id);
	energy = EmulsionEnergyFromRange(range, particle_id);
	range  = 210.e-3 * scale_factor + PolystyreneRangeFromEnergy(energy, particle_id);
	energy = PolystyreneEnergyFromRange(range, particle_id);
	range  = 70.e-3 * scale_factor + EmulsionRangeFromEnergy(energy, particle_id);
	energy = EmulsionEnergyFromRange(range, particle_id);	
      }
      else { // upstream of water
	if ( ipl == 0 )  // if stopping plate
	  range = 500.e-3 * scale_factor * 0.5;
	else
	  range = 500.e-3 * scale_factor + IronRangeFromEnergy(energy, particle_id);
	energy = IronEnergyFromRange(range, particle_id);
	range  = 70.e-3 * scale_factor + EmulsionRangeFromEnergy(energy, particle_id);
	energy = EmulsionEnergyFromRange(range, particle_id);
	range  = 210.e-3 * scale_factor + PolystyreneRangeFromEnergy(energy, particle_id);
	energy = PolystyreneEnergyFromRange(range, particle_id);
	range  = 70.e-3 * scale_factor + EmulsionRangeFromEnergy(energy, particle_id);
	energy = EmulsionEnergyFromRange(range, particle_id);
	range  = 109.e-3 * scale_factor + PolystyreneRangeFromEnergy(energy, particle_id);
	energy = PolystyreneEnergyFromRange(range, particle_id);
      }
    }
    else if ( pl.at(ipl) == 15 ) { // most upstream of Fe ECC
      if ( ipl == 0 ) { // if stopping plate
	range  = 210.e-3 * scale_factor * 0.5;
	energy = PolystyreneEnergyFromRange(range, particle_id);
	range  = 70.e-3 * scale_factor + EmulsionRangeFromEnergy(energy, particle_id);
      }
      range  = 70.e-3 * scale_factor + EmulsionRangeFromEnergy(energy, particle_id);
      energy = EmulsionEnergyFromRange(range, particle_id);
      range  = 210.e-3 * scale_factor + PolystyreneRangeFromEnergy(energy, particle_id);
      energy = PolystyreneEnergyFromRange(range, particle_id);
      range  = 70.e-3 * scale_factor + EmulsionRangeFromEnergy(energy, particle_id);
      energy = EmulsionEnergyFromRange(range, particle_id);
    }
    else if ( pl.at(ipl) >= 5 ) { // Fe ECC
      if ( ipl == 0 )  // if stopping plate
	range = 500.e-3 * scale_factor * 0.5;
      else
	range = 500.e-3 * scale_factor + IronRangeFromEnergy(energy, particle_id);
      energy = IronEnergyFromRange(range, particle_id);
      range  = 70.e-3 * scale_factor + EmulsionRangeFromEnergy(energy, particle_id);
      energy = EmulsionEnergyFromRange(range, particle_id);
      range  = 210.e-3 * scale_factor + PolystyreneRangeFromEnergy(energy, particle_id);
      energy = PolystyreneEnergyFromRange(range, particle_id);
      range  = 70.e-3 * scale_factor + EmulsionRangeFromEnergy(energy, particle_id);
      energy = EmulsionEnergyFromRange(range, particle_id);      
    }
    else {
      throw std::invalid_argument("Backward track from interaction cannot be exist in ISS");
    }
  }

  return energy;

}

double RangeFunction::CalculateProtonRangeError(double range_mom, double tangent) const {

  if ( range_mom < 200.) return range_mom * 0.03;
  else if ( range_mom < 300. ) {
    if ( tangent < 1. ) return range_mom * 0.03;
    else if ( tangent < 2. ) return range_mom * 0.04;
    else if ( tangent < 3. ) return range_mom * 0.06;
    else return range_mom * 0.02;    
  }
  else if ( range_mom < 400. ) {
    if ( tangent < 1. ) return range_mom * 0.01;
    else return range_mom * 0.02;
  }
  else if ( range_mom < 500. ) {
    if ( tangent < 1. ) return range_mom * 0.005;
    else return range_mom * 0.01;
  }
  else return range_mom * 0.01;

}
