// boost includes
#include <boost/log/core.hpp>
#include <boost/log/trivial.hpp>
#include <boost/log/expressions.hpp>

// system includes
#include <iostream>
#include <vector>
#include <algorithm>
#include <numeric>
#include <functional>

#include <TRandom.h>

// my includes
#include "McsConst.hpp"
#include "McsClass.hpp"
#include "McsFunction.hpp"

namespace logging = boost::log;

int main (int argc, char* argv[]) {

  logging::core::get()->set_filter
    (
     //logging::trivial::severity >= logging::trivial::debug
     logging::trivial::severity >= logging::trivial::info
     );

  BOOST_LOG_TRIVIAL(info) << "=========Momentum Reconstruction Start==========";

  if ( argc != 5 ) {
    BOOST_LOG_TRIVIAL(error) << "Usage : " << argv[0]
			     << " <input momch file name> <output momch file name> <material mode> <seed>";
    std::exit(1);
  }

  try {

    std::string ifilename = argv[1];
    auto ev_vec = Momentum_recon::ReadEventInformationBin(ifilename);

    Int_t num_entry = 0;

    long seed = std::atoi(argv[4]);
    if ( seed == 0 )
      gRandom->SetSeed(time(NULL));
    else
      gRandom->SetSeed(seed);

    for ( auto &ev : ev_vec ) {

      BOOST_LOG_TRIVIAL(debug) << "Entry : " << ev.groupid;
      num_entry++;
      // if ( num_entry > 10) break;
            
      if ( ev.chains.empty() ) continue;

      for ( auto &chain : ev.chains ) {

	if ( chain.base_pair.empty() ) continue;

	std::vector<Double_t > basetrack_distance = {};
	std::vector<Double_t > basetrack_distance_water = {};
	std::vector<Double_t > track_tangent = {};
	std::vector<Int_t > plate_id = {};
	std::vector<Int_t > plate_id_next = {};
	std::vector<Double_t > radial_angle_difference = {};
	std::vector<Double_t > lateral_angle_difference = {};

	double initial_pbeta = 0.;
	double radial_cut_value, lateral_cut_value;	
	double radial_cut_value_water, lateral_cut_value_water;

	Double_t radial_angle_difference_rms = 0.;
	Double_t lateral_angle_difference_rms = 0.;
	Double_t radial_angle_difference_water_rms = 0.;
	Double_t lateral_angle_difference_water_rms = 0.;
	Int_t num_angle_difference = 0;
	Int_t num_angle_difference_water = 0;

	if ( chain.direction == 1 ) { // forward

	  for ( auto itr = chain.base_pair.rbegin(); itr != chain.base_pair.rend(); itr++ ) {
	    TVector3 up_position((*itr).second.x, (*itr).second.y, (*itr).second.z);
	    TVector3 down_position((*itr).first.x, (*itr).first.y, (*itr).first.z);
	    TVector3 up_tangent((*itr).second.ax, (*itr).second.ay, 1.);
	    TVector3 down_tangent((*itr).first.ax, (*itr).first.ay, 1.);

	    track_tangent.push_back(up_tangent.Mag());

	    int plate = (*itr).second.pl;
	    int plate_next = (*itr).first.pl;
	    int plate_diff = plate - plate_next;
	    plate_id.push_back(plate);
	    plate_id_next.push_back(plate_next);

	    TVector3 basetrack_distance_pair(0.,0.,0.);
	    TVector3 basetrack_distance_water_pair(0.,0.,0.);
	    GetBasetrackDistancePair(plate, plate_next, chain.direction,
				     up_position, down_position,
				     basetrack_distance_pair,
				     basetrack_distance_water_pair,
				     kFALSE);
	    
	    basetrack_distance.push_back(basetrack_distance_pair.Mag());
	    basetrack_distance_water.push_back(basetrack_distance_water_pair.Mag());
	    
	    double radial_angle_difference_pair = RadialAngleDiffNew(up_tangent, down_tangent);
	    double lateral_angle_difference_pair = LateralAngleDiffNew(up_tangent, down_tangent);
	    
	    radial_angle_difference.push_back(radial_angle_difference_pair);
	    lateral_angle_difference.push_back(lateral_angle_difference_pair);

	    if ( plate_diff == 1 ) {
	      if ( (plate > 4 && plate < 16) ||
		   plate % 2 == 1 ) { // 鉄
		radial_angle_difference_rms += radial_angle_difference_pair * radial_angle_difference_pair;
		lateral_angle_difference_rms += lateral_angle_difference_pair * lateral_angle_difference_pair;
		num_angle_difference++;
	      }
	      else if ( plate > 17 && plate % 2 == 0 ) { // 水
		radial_angle_difference_water_rms += radial_angle_difference_pair * radial_angle_difference_pair;
		lateral_angle_difference_water_rms += lateral_angle_difference_pair * lateral_angle_difference_pair;
		num_angle_difference_water++;
	      }
	    }

	  }

	  // initial pbeta value
	  radial_angle_difference_rms /= num_angle_difference;
	  radial_angle_difference_rms = TMath::Sqrt(radial_angle_difference_rms);
	  lateral_angle_difference_rms /= num_angle_difference;
	  lateral_angle_difference_rms = TMath::Sqrt(lateral_angle_difference_rms);
	  radial_angle_difference_water_rms /= num_angle_difference_water;
	  radial_angle_difference_water_rms = TMath::Sqrt(radial_angle_difference_water_rms);
	  lateral_angle_difference_water_rms /= num_angle_difference_water;
	  lateral_angle_difference_water_rms = TMath::Sqrt(lateral_angle_difference_water_rms);

	  // fit initial parameters set
	  radial_cut_value = 3. * radial_angle_difference_rms;
	  lateral_cut_value = 3. * lateral_angle_difference_rms;
	  radial_cut_value_water = 3. * radial_angle_difference_water_rms;
	  lateral_cut_value_water = 3. * lateral_angle_difference_water_rms;
	  double radiation_length = CalcRadLength(1, track_tangent.front(), kNinjaIron);
	  initial_pbeta = MCS_SCALE_FACTOR * 13.6 / lateral_angle_difference_rms * TMath::Sqrt(radiation_length) * (1. + 0.038 * TMath::Log(radiation_length));

	} else if ( chain.direction == -1 ) { // backward
	  
	  for ( auto itr = chain.base_pair.begin(); itr != chain.base_pair.end(); itr++ ) {
	    TVector3 up_position((*itr).first.x, (*itr).first.y, (*itr).first.z);
	    TVector3 down_position((*itr).second.x, (*itr).second.y, (*itr).second.z);
	    TVector3 up_tangent((*itr).first.ax, (*itr).first.ay, 1.);
	    TVector3 down_tangent((*itr).second.ax, (*itr).second.ay, 1.);

	    track_tangent.push_back(up_tangent.Mag());

	    int plate = (*itr).first.pl;
	    int plate_next = (*itr).second.pl;
	    int plate_diff = plate_next - plate;
	    plate_id.push_back(plate);
	    plate_id_next.push_back(plate_next);

	    TVector3 basetrack_distance_pair(0.,0.,0.);
	    TVector3 basetrack_distance_water_pair(0.,0.,0.);
	    GetBasetrackDistancePair(plate, plate_next, chain.direction,
				     up_position, down_position,
				     basetrack_distance_pair,
				     basetrack_distance_water_pair,
				     kFALSE);

	    basetrack_distance.push_back(basetrack_distance_pair.Mag());
	    basetrack_distance_water.push_back(basetrack_distance_water_pair.Mag());
	    
	    double radial_angle_difference_pair = RadialAngleDiffNew(down_tangent, up_tangent);
	    double lateral_angle_difference_pair = LateralAngleDiffNew(down_tangent, up_tangent);

	    radial_angle_difference.push_back(radial_angle_difference_pair);
	    lateral_angle_difference.push_back(lateral_angle_difference_pair);

	    if ( plate_diff == 1 ) {
	      if ( (plate > 3 && plate < 15) ||
		   plate % 2 == 0 ) { //鉄
		radial_angle_difference_rms += radial_angle_difference_pair * radial_angle_difference_pair;
		lateral_angle_difference_rms += lateral_angle_difference_pair * lateral_angle_difference_pair;
		num_angle_difference++;
	      }
	      else if ( plate > 16 && plate % 2 == 1 ) { // 水
		radial_angle_difference_water_rms += radial_angle_difference_pair * radial_angle_difference_pair;
		lateral_angle_difference_water_rms += lateral_angle_difference_pair * lateral_angle_difference_pair;
		num_angle_difference_water++;
	      }
	    }

	  }

	  // initial pbeta value
	  radial_angle_difference_rms /= num_angle_difference;
	  radial_angle_difference_rms = TMath::Sqrt(radial_angle_difference_rms);
	  lateral_angle_difference_rms /= num_angle_difference;
	  lateral_angle_difference_rms = TMath::Sqrt(lateral_angle_difference_rms);
	  radial_angle_difference_water_rms /= num_angle_difference_water;
	  radial_angle_difference_water_rms = TMath::Sqrt(radial_angle_difference_water_rms);
	  lateral_angle_difference_water_rms /= num_angle_difference_water;
	  lateral_angle_difference_water_rms = TMath::Sqrt(lateral_angle_difference_water_rms);

	  // fit initial parameters set
	  radial_cut_value = 3. * radial_angle_difference_rms;
	  lateral_cut_value = 3. * lateral_angle_difference_rms;
	  radial_cut_value_water = 3. * radial_angle_difference_water_rms;
	  lateral_cut_value_water = 3. * lateral_angle_difference_water_rms;
	  double radiation_length = CalcRadLength(1, track_tangent.front(), kNinjaIron);
	  initial_pbeta = MCS_SCALE_FACTOR * 13.6 / lateral_angle_difference_rms * TMath::Sqrt(radiation_length) * (1. + 0.038 * TMath::Log(radiation_length));

	}
	else {
	  throw std::runtime_error("Direction should be +/- 1");
	  std::exit(1);
	}

	// Reconstruction
	// assume muon mass
	// assume proton mass
	// assume pion mass -> baby mind curvature にひとまず詰める？

	// Combo にするのか，鉄だけなのかは何で決める？
	int material_mode = std::atoi(argv[3]);
	// 0 : iron, 1 : water, 2 : combo
	switch ( material_mode ) {
	case kNinjaIron : 
	  BOOST_LOG_TRIVIAL(debug) << "Iron only reconstruction";
	  break;
	case kNinjaWater : 
	  BOOST_LOG_TRIVIAL(debug) << "Water only reconstruction";
	  break;
	case 2 :
	  BOOST_LOG_TRIVIAL(debug) << "Iron + water reconstruction";
	  break;
	default :
	  throw std::runtime_error("material mode not appropriate");
	}
	/*
	for ( int i = 0; i < basetrack_distance.size(); i++ ) {
	  std::cout << "basetrack distance : " << basetrack_distance.at(i) << std::endl;
	  std::cout << "basetrack distance water : " << basetrack_distance_water.at(i) << std::endl;
	  std::cout << "track tangent : " << track_tangent.at(i) << std::endl;
	  std::cout << "plate id  : " << plate_id.at(i) << std::endl;
	  std::cout << "plate id next : " << plate_id_next.at(i) << std::endl;
	  std::cout << "Radial angle difference : " << radial_angle_difference.at(i) << std::endl;
	  std::cout << "Lateral angle difference : " << lateral_angle_difference.at(i) << std::endl;
	}
	*/
	for ( int particle_id = 0; particle_id < kNumberOfNinjaMcsParticles; particle_id++ ) {

	  // lateral だけで reconstruction
	  auto results = ReconstructPBeta(initial_pbeta,
					  1, particle_id,
					  chain.direction,
					  radial_cut_value,
					  lateral_cut_value,
					  radial_cut_value_water,
					  lateral_cut_value_water,
					  kTRUE,
					  material_mode, 1,
					  basetrack_distance,
					  basetrack_distance_water,
					  track_tangent,
					  plate_id,
					  plate_id_next,
					  radial_angle_difference,
					  lateral_angle_difference);

	  double pbeta = results.at(0);
	  double pbeta_err_minus = results.at(3);
	  double pbeta_err_plus = results.at(2);
	  double momentum_minus, momentum_plus; // p +/- 1 sigma
	  // mis fitting treatment
	  if ( std::fabs(pbeta - 20.) < 1.e-4 ) {
	    double nll_min = FuncNegativeLogLikelihood(20., 1, particle_id, chain.direction,
						       radial_cut_value, lateral_cut_value,
						       radial_cut_value_water,
						       lateral_cut_value_water,
						       kTRUE,
						       material_mode, 1,
						       basetrack_distance,
						       basetrack_distance_water,
						       track_tangent,
						       plate_id,
						       plate_id_next,
						       radial_angle_difference,
						       lateral_angle_difference);
	    for ( double ipbeta = 20.; ipbeta < 10000.; ipbeta += 1. ) {
	      double nll = FuncNegativeLogLikelihood(ipbeta, 1, particle_id, chain.direction,
						     radial_cut_value, lateral_cut_value,
						     radial_cut_value_water,
						     lateral_cut_value_water,
						     kTRUE,
						     material_mode, 1,
						     basetrack_distance,
						     basetrack_distance_water,
						     track_tangent,
						     plate_id,
						     plate_id_next,
						     radial_angle_difference,
						     lateral_angle_difference);
	      if ( nll - nll_min > 1. ) {
		pbeta_err_plus = ipbeta - pbeta;
		break;
	      }
	    }
	    if ( std::fabs(pbeta_err_plus) < 1.e-4 ) {
	      pbeta_err_minus = 10000. - pbeta;
	    }
	  }
	  else if ( std::fabs(pbeta - 10000.) < 1.e-4 ) {
	    double nll_min = FuncNegativeLogLikelihood(10000., 1, particle_id, chain.direction,
						       radial_cut_value, lateral_cut_value,
						       radial_cut_value_water,
						       lateral_cut_value_water,
						       kTRUE,
						       material_mode, 1,
						       basetrack_distance,
						       basetrack_distance_water,
						       track_tangent,
						       plate_id,
						       plate_id_next,
						       radial_angle_difference,
						       lateral_angle_difference);
	    for ( double ipbeta = 10000.; ipbeta >= 20.; ipbeta -= 1. ) {
	      double nll = FuncNegativeLogLikelihood(ipbeta, 1, particle_id, chain.direction,
						     radial_cut_value, lateral_cut_value,
						     radial_cut_value_water,
						     lateral_cut_value_water,
						     kTRUE,
						     material_mode, 1,
						     basetrack_distance,
						     basetrack_distance_water,
						     track_tangent,
						     plate_id,
						     plate_id_next,
						     radial_angle_difference,
						     lateral_angle_difference);
	      if ( nll - nll_min > 1. ) {
		pbeta_err_minus = ipbeta - pbeta;
		break;
	      }
	    }
	    if ( std::fabs(pbeta_err_minus) < 1.e-4 ) {
	      pbeta_err_minus = 20. - pbeta;
	    }
	  }
	  else {
	    if ( std::fabs(pbeta_err_minus) < 1.e-4 ) {
	    pbeta_err_minus = 20. - pbeta;
	    }
	    if ( std::fabs(pbeta_err_plus) < 1.e-4 ) {
	      pbeta_err_plus = 10000. - pbeta;
	    }
	  }

	  // ある値より小さければ radial + lateral で測定し直す

	  if ( (particle_id == 0 && CalculateMomentumFromPBeta(pbeta, MCS_MUON_MASS) < 500.) || // muon
	       (particle_id == 1 && CalculateMomentumFromPBeta(pbeta, MCS_PION_MASS) < 500.) || // pion
	       (particle_id == 2 && CalculateMomentumFromPBeta(pbeta, MCS_PROTON_MASS) < 700.) ) { // proton

	    auto results = ReconstructPBeta(initial_pbeta,
					    1, particle_id,
					    chain.direction,
					    radial_cut_value,
					    lateral_cut_value,
					    radial_cut_value_water,
					    lateral_cut_value_water,
					    kTRUE,
					    material_mode, 0,
					    basetrack_distance,
					    basetrack_distance_water,
					    track_tangent,
					    plate_id,
					    plate_id_next,
					    radial_angle_difference,
					    lateral_angle_difference);

	    pbeta = results.at(0);
	    pbeta_err_minus = results.at(3);
	    pbeta_err_plus = results.at(2);

	    // mis fitting treatment
	    if ( std::fabs(pbeta - 20.) < 1.e-4 ) {
	      double nll_min = FuncNegativeLogLikelihood(20., 1, particle_id, chain.direction,
							 radial_cut_value, lateral_cut_value,
							 radial_cut_value_water,
							 lateral_cut_value_water,
							 kTRUE,
							 material_mode, 1,
							 basetrack_distance,
							 basetrack_distance_water,
							 track_tangent,
							 plate_id,
							 plate_id_next,
							 radial_angle_difference,
							 lateral_angle_difference);
	      for ( double ipbeta = 20.; ipbeta < 10000.; ipbeta += 1. ) {
		double nll = FuncNegativeLogLikelihood(ipbeta, 1, particle_id, chain.direction,
						       radial_cut_value, lateral_cut_value,
						       radial_cut_value_water,
						       lateral_cut_value_water,
						       kTRUE,
						       material_mode, 1,
						       basetrack_distance,
						       basetrack_distance_water,
						       track_tangent,
						       plate_id,
						       plate_id_next,
						       radial_angle_difference,
						       lateral_angle_difference);
		if ( nll - nll_min > 1. ) {
		  pbeta_err_plus = ipbeta - pbeta;
		  break;
		}
	      }
	      if ( std::fabs(pbeta_err_plus) < 1.e-4 ) {
		pbeta_err_plus = 10000. - pbeta;
	      }
	    }
	    else if ( std::fabs(pbeta - 10000.) < 1.e-4 ) {
	      double nll_min = FuncNegativeLogLikelihood(10000., 1, particle_id, chain.direction,
							 radial_cut_value, lateral_cut_value,
							 radial_cut_value_water,
							 lateral_cut_value_water,
							 kTRUE,
							 material_mode, 1,
							 basetrack_distance,
							 basetrack_distance_water,
							 track_tangent,
							 plate_id,
							 plate_id_next,
							 radial_angle_difference,
							 lateral_angle_difference);
	      for ( double ipbeta = 10000.; ipbeta >= 20.; ipbeta -= 1. ) {
		double nll = FuncNegativeLogLikelihood(ipbeta, 1, particle_id, chain.direction,
						       radial_cut_value, lateral_cut_value,
						       radial_cut_value_water,
						       lateral_cut_value_water,
						       kTRUE,
						       material_mode, 1,
						       basetrack_distance,
						       basetrack_distance_water,
						       track_tangent,
						       plate_id,
						       plate_id_next,
						       radial_angle_difference,
						       lateral_angle_difference);
		if ( nll - nll_min > 1. ) {
		  pbeta_err_minus = ipbeta - pbeta;
		  break;
		}
	      }
	      if ( std::fabs(pbeta_err_minus) < 1.e-4 ) {
		pbeta_err_minus = 20. - pbeta;
	      }
	    }
	    else {
	      if ( std::fabs(pbeta_err_minus) < 1.e-4 ) {
		pbeta_err_minus = 20. - pbeta;
	      }
	      if ( std::fabs(pbeta_err_plus) < 1.e-4 ) {
		pbeta_err_plus = 10000. - pbeta;
	      }
	    }	    	    
	  }

	  if ( particle_id == 0 ) { // muon
	    chain.ecc_mcs_mom[0] = CalculateMomentumFromPBeta(pbeta, MCS_MUON_MASS);
	    momentum_minus = CalculateMomentumFromPBeta(pbeta + pbeta_err_minus, MCS_MUON_MASS);
	    momentum_plus = CalculateMomentumFromPBeta(pbeta + pbeta_err_plus, MCS_MUON_MASS);
	    chain.ecc_mcs_mom_error[0][0] = chain.ecc_mcs_mom[0] - momentum_minus;
	    chain.ecc_mcs_mom_error[0][1] = momentum_plus - chain.ecc_mcs_mom[0];
	  }
	  else if ( particle_id == 1 ) { // pion
	    chain.bm_curvature_mom = CalculateMomentumFromPBeta(pbeta, MCS_PION_MASS);
	    momentum_minus = CalculateMomentumFromPBeta(pbeta + pbeta_err_minus, MCS_PION_MASS);
	    momentum_plus = CalculateMomentumFromPBeta(pbeta + pbeta_err_plus, MCS_PION_MASS);
	    chain.bm_curvature_mom_error[0] = chain.bm_curvature_mom - momentum_minus;
	    chain.bm_curvature_mom_error[1] = momentum_plus - chain.bm_curvature_mom;
	  }
	  else if ( particle_id == 2 ) { // proton
	    chain.ecc_mcs_mom[1] = CalculateMomentumFromPBeta(pbeta, MCS_PROTON_MASS);
	    momentum_minus = CalculateMomentumFromPBeta(pbeta + pbeta_err_minus, MCS_PROTON_MASS);
	    momentum_plus = CalculateMomentumFromPBeta(pbeta + pbeta_err_plus, MCS_PROTON_MASS);
	    chain.ecc_mcs_mom_error[1][0] = chain.ecc_mcs_mom[1] - momentum_minus;
	    chain.ecc_mcs_mom_error[1][1] = momentum_plus - chain.ecc_mcs_mom[1];
	  }	 
	}
      }
      
    }

    std::string ofilename = argv[2];
    Momentum_recon::WriteEventInformationBin(ofilename, ev_vec);

  } catch (const std::runtime_error &error) {
    BOOST_LOG_TRIVIAL(fatal) << "Runtime error : " << error.what();
    std::exit(1);
  } catch (const std::invalid_argument &error) {
    BOOST_LOG_TRIVIAL(fatal) << "Invalid argument error : " << error.what();
    std::exit(1);
  }

  BOOST_LOG_TRIVIAL(info) << "===========Momentum Reconstruction Finish==========";
  std::exit(0);

}
