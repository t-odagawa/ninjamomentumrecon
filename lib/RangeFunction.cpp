#include <boost/log/core.hpp>
#include <boost/log/trivial.hpp>
#include <boost/log/expressions.hpp>

#include <vector>
#include <cmath>

double IronRangeFromMomentum(double momentum) {
  return 0.;
}

double IronMomentumFromRange(double range) {
  return 0.;
}

double WaterRangeFromMomentum(double momentum) {
  return 0.;
}

double WaterMomentumFromRange(double range) {
  return 0.;
}

double PolystyreneRangeFromMomentum(double momentum) {
  return 0.;
}

double PolystyreneMomentumFromRange(double range) {
  return 0.;
}

double GelRangeFromMomentum(double momentum) {
  return 0.;
}

double GelMomentumFromRange(double range) {
  return 0.;
}

void ModifyVectors(std::vector<double> &ax, std::vector<double> &ay, std::vector<int> &pl) {

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

}

double CalculateMomentumFromRange(std::vector<double> ax, std::vector<double> ay, std::vector<int> pl) {

  if ( ax.size() != ay.size() ||
       ax.size() != pl.size() ||
       ay.size() != pl.size() ) {
    BOOST_LOG_TRIVIAL(error) << "Size of input vectors should be the same";
    std::exit(1);
  }

  ModifyVectors(ax, ay, pl);

  double momentum = 0.;
  double scale_factor = 1.;

  for ( int ipl = 0; ipl < pl.size(); ipl++ ) {

    scale_factor = std::sqrt(ax.at(ipl) * ax.at(ipl) + ay.at(ipl) * ay.at(ipl) + 1.);

    double range = 0.;
    double range_tmp = 0.;

    if ( pl.at(ipl) == 3 ) { // ISS downstream
      range += 1. * scale_factor;
      momentum = PolystyreneMomentumFromRange(range);      
      range += 70.e-3 * scale_factor + GelRangeFromMomentum(momentum);
      momentum = GelMomentumFromRange(range);
      range += 210.e-3 * scale_factor + PolystyreneRangeFromMomentum(momentum);
      momentum = PolystyreneMomentumFromRange(range);
      range += 70.e-3 * scale_factor + GelRangeFromMomentum(momentum);
      momentum = GelMomentumFromRange(range);
    }
    else if ( pl.at(ipl) == 4 ) { // ISS upstream
      if ( ipl == 0 ) { // if stopping plate
	range += 210.e-3 * scale_factor * 0.5;
	momentum = PolystyreneMomentumFromRange(range);
	range += 70.e-3 * scale_factor + GelRangeFromMomentum(momentum);
      }
      range += 70.e-3 * scale_factor + GelRangeFromMomentum(momentum);
      momentum = GelMomentumFromRange(range);
      range += 210.e-3 * scale_factor + PolystyreneRangeFromMomentum(momentum);
      momentum = PolystyreneMomentumFromRange(range);
      range += 70.e-3 * scale_factor + GelRangeFromMomentum(momentum);
      momentum = GelMomentumFromRange(range);
    }
    else if ( pl.at(ipl) < 15 ) { // Fe ECC
      if ( ipl == 0 )  // if stopping plate
	range += 500.e-3 * scale_factor * 0.5;
      else
	range += 500.e-3 * scale_factor + IronRangeFromMomentum(momentum);
      momentum = IronMomentumFromRange(range);
      range += 70.e-3 * scale_factor + GelRangeFromMomentum(momentum);
      momentum = GelMomentumFromRange(range);
      range += 210.e-3 * scale_factor + PolystyreneRangeFromMomentum(momentum);
      momentum = PolystyreneMomentumFromRange(range);
      range += 70.e-3 * scale_factor + GelRangeFromMomentum(momentum);
      momentum = GelMomentumFromRange(range);
    }
    else if ( pl.at(ipl) == 16 ) { // most downstream of water ECC
      if ( ipl == 0 ) { // if stopping plate
	range += 210.e-3 * scale_factor * 0.5;
	momentum = PolystyreneMomentumFromRange(range);
	range += 70.e-3 * scale_factor + GelRangeFromMomentum(momentum);	
      }
      range += 70.e-3 * scale_factor + GelRangeFromMomentum(momentum);
      momentum = GelMomentumFromRange(range);
      range += 210.e-3 * scale_factor + PolystyreneRangeFromMomentum(momentum);
      momentum = PolystyreneMomentumFromRange(range);
      range += 70.e-3 * scale_factor + GelRangeFromMomentum(momentum);
      momentum = GelMomentumFromRange(range);
    }
    else if ( pl.at(ipl) % 2 == 1 ) { // upstream of iron
      if ( ipl == 0 )  // if stopping plate
	range += 500.e-3 * scale_factor * 0.5;
      else
	range += 500.e-3 * scale_factor + IronRangeFromMomentum(momentum);
      momentum = IronMomentumFromRange(range);
      range += 70.e-3 * scale_factor + GelRangeFromMomentum(momentum);
      momentum = GelMomentumFromRange(range);
      range += 210.e-3 * scale_factor + PolystyreneRangeFromMomentum(momentum);
      momentum = PolystyreneMomentumFromRange(range);
      range += 70.e-3 * scale_factor + GelRangeFromMomentum(momentum);
      momentum = GelMomentumFromRange(range);      
    }
    else { // upstream of water
      if ( ipl == 0 )  // if stopping plate
	range += 2.3 * scale_factor * 0.5;
      else
	range += 2.3 * scale_factor + WaterRangeFromMomentum(momentum);
      momentum = WaterMomentumFromRange(range);
      range += 70.e-3 * scale_factor + GelRangeFromMomentum(momentum);
      momentum = GelMomentumFromRange(range);
      range += 210.e-3 * scale_factor + PolystyreneRangeFromMomentum(momentum);
      momentum = PolystyreneMomentumFromRange(range);
      range += 70.e-3 * scale_factor + GelRangeFromMomentum(momentum);
      momentum = GelMomentumFromRange(range);
    }
  }

  for ( int ipl = 0; ipl < ax.size(); ipl++ )
    momentum += std::sqrt(ax.at(ipl) * ax.at(ipl) + ay.at(ipl) * ay.at(ipl) + 1.) / 0.8;
  
  return momentum;
}
