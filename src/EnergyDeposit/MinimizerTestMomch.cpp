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

#include "McsConst.hpp"
#include "McsClass.hpp"
#include "McsFunction.hpp"

namespace logging = boost::log;

int main (int argc, char *argv[]) {

  logging::core::get()->set_filter
    (
     logging::trivial::severity >= logging::trivial::info
     //logging::trivial::severity >= logging::trivial::trace
     //logging::trivial::severity >= logging::trivial::debug
     );
  
  BOOST_LOG_TRIVIAL(info) << "==========Momentum Reconstruction Start==========";

  if ( argc != 4 ) {
    BOOST_LOG_TRIVIAL(error) << "Usage : " << argv[0]
			     << " <input momch file name> <output momch file name> <ncell = 1>";
    std::exit(1);
  }

  try {
    std::string ifilename = argv[1];
    std::ifstream ifs(ifilename, std::ios::binary);

    std::string ofilename = argv[2];
    std::ofstream ofs(ofilename, std::ios::binary);

    Momentum_recon::Event_information ev;
    Momentum_recon::Mom_chain mom_chain;
    Momentum_recon::Mom_basetrack base;
    std::pair<Momentum_recon::Mom_basetrack, Momentum_recon::Mom_basetrack > base_pair;
    int num_chain = 0;
    int num_true_chain = 0;
    int num_base = 0;
    int num_link = 0;
      
    const Int_t ncell = std::atoi(argv[3]);

    Int_t num_entry = 0;

    while ( Momentum_recon::ReadEventInformationHeader(ifs, ev, num_chain, num_true_chain) ) {
      num_entry++;

      // chains
      for ( int ichain = 0; ichain < num_chain; ichain++ ) {
	if ( Momentum_recon::ReadMomChainHeader(ifs, mom_chain, num_base, num_link) ) {
	  mom_chain.base.clear();
	  mom_chain.base_pair.clear();
	  mom_chain.base.reserve(num_base);
	  mom_chain.base_pair.reserve(num_link);
	  for ( int i = 0; i < num_base; i++ ) {
	    ifs.read((char*)& base, sizeof(Momentum_recon::Mom_basetrack));
	    mom_chain.base.push_back(base);
	  }
	  for ( int i = 0; i < num_link; i++ ) {
	    ifs.read((char*)& base_pair.first, sizeof(Momentum_recon::Mom_basetrack));
	    ifs.read((char*)& base_pair.second, sizeof(Momentum_recon::Mom_basetrack));
	    mom_chain.base_pair.push_back(base_pair);
	  }

	  if ( num_link > 0 ) {

	    // base pair vectors (backward = base_pair_vec, forward = reverse_base_pair_vec)
	    std::vector<std::pair<Momentum_recon::Mom_basetrack, Momentum_recon::Mom_basetrack> > base_pair_vec;
	    std::vector<std::pair<Momentum_recon::Mom_basetrack, Momentum_recon::Mom_basetrack> > reverse_base_pair_vec;
	    base_pair_vec.resize(num_link);
	    reverse_base_pair_vec.resize(num_link);
	    std::copy(mom_chain.base_pair.begin(), mom_chain.base_pair.end(), base_pair_vec.begin());
	    std::copy(mom_chain.base_pair.begin(), mom_chain.base_pair.end(), reverse_base_pair_vec.begin());
	    std::reverse(reverse_base_pair_vec.begin(), reverse_base_pair_vec.end());
	    
	    // input parameters for TMinuit (forward/backward)
	    std::vector<Double_t > basetrack_distance = {};
	    std::vector<Double_t > basetrack_distance_back = {};
	    std::vector<Double_t > water_basetrack_distance = {};
	    std::vector<Double_t > water_basetrack_distance_back = {};
	    std::vector<Double_t > track_tangent = {};
	    std::vector<Double_t > track_tangent_back = {};
	    std::vector<Int_t > plate_id = {};
	    std::vector<Int_t > plate_id_back = {};
	    std::vector<Double_t > radial_angle_difference = {};
	    std::vector<Double_t > radial_angle_difference_back = {};
	    std::vector<Double_t > lateral_angle_difference = {};
	    std::vector<Double_t > lateral_angle_difference_back = {};
	    Int_t iron_plate_id = -1;
	    Int_t water_plate_id = -1;
	    
	    // forward
	    for ( Int_t ipair = 0; ipair < reverse_base_pair_vec.size(); ipair++ ) {
	      TVector3 upstream_position;
	      upstream_position.SetX(reverse_base_pair_vec.at(ipair).second.x / 1.e3);
	      upstream_position.SetY(reverse_base_pair_vec.at(ipair).second.y / 1.e3);
	      upstream_position.SetZ(reverse_base_pair_vec.at(ipair).second.z / 1.e3);
	      TVector3 downstream_position;
	      downstream_position.SetX(reverse_base_pair_vec.at(ipair).first.x / 1.e3);
	      downstream_position.SetY(reverse_base_pair_vec.at(ipair).first.y / 1.e3);
	      downstream_position.SetZ(reverse_base_pair_vec.at(ipair).first.z / 1.e3);
	      TVector3 upstream_tangent;
	      upstream_tangent.SetX(reverse_base_pair_vec.at(ipair).second.ax);
	      upstream_tangent.SetY(reverse_base_pair_vec.at(ipair).second.ay);
	      upstream_tangent.SetZ(1.);
	      TVector3 downstream_tangent;
	      downstream_tangent.SetX(reverse_base_pair_vec.at(ipair).first.ax);
	      downstream_tangent.SetY(reverse_base_pair_vec.at(ipair).first.ay);
	      downstream_tangent.SetZ(1.);
	      
	      if ( reverse_base_pair_vec.at(ipair).second.pl - 1 <= 3 ) continue;
	      if ( (reverse_base_pair_vec.at(ipair).second.pl - 1) % 2 == 1 &&
		   reverse_base_pair_vec.at(ipair).second.pl - 1 >= 15 ) continue;

	      iron_plate_id = reverse_base_pair_vec.at(ipair).second.pl - 1;
	      water_plate_id = reverse_base_pair_vec.at(ipair).first.pl - 1;
	      // if ( std::abs(iron_plate_id - water_plate_id) > 1 ) continue;
	      
	      basetrack_distance.push_back((downstream_position - upstream_position).Mag() / 850.e-3);
	      track_tangent.push_back(upstream_tangent.Mag());
	      plate_id.push_back(reverse_base_pair_vec.at(ipair).second.pl - 1);
	      radial_angle_difference.push_back(RadialAngleDiffNew(upstream_tangent, downstream_tangent));
	      lateral_angle_difference.push_back(LateralAngleDiffNew(upstream_tangent, downstream_tangent));
	      if ( ipair < reverse_base_pair_vec.size() - 1 &&
		   reverse_base_pair_vec.at(ipair + 1).second.pl - 1 == water_plate_id &&
		   iron_plate_id >= 18 ) {
		TVector3 water_upstream_position;
		water_upstream_position.SetX(reverse_base_pair_vec.at(ipair + 1).second.x / 1.e3);
		water_upstream_position.SetY(reverse_base_pair_vec.at(ipair + 1).second.y / 1.e3);
		water_upstream_position.SetZ(reverse_base_pair_vec.at(ipair + 1).second.z / 1.e3);
		TVector3 water_downstream_position;
		water_downstream_position.SetX(reverse_base_pair_vec.at(ipair + 1).first.x / 1.e3);
		water_downstream_position.SetY(reverse_base_pair_vec.at(ipair + 1).first.y / 1.e3);
		water_downstream_position.SetZ(reverse_base_pair_vec.at(ipair + 1).first.z / 1.e3);
		water_basetrack_distance.push_back((water_downstream_position - water_upstream_position).Mag() / 2.868);
	      } else {
		water_basetrack_distance.push_back(0.);
	      }
	      BOOST_LOG_TRIVIAL(debug) << " Film " << reverse_base_pair_vec.at(ipair).second.pl - 1 << " and"
				       << " Film " << reverse_base_pair_vec.at(ipair).first.pl - 1 << ", "
				       << " Radial angle difference is " << radial_angle_difference.back() << ", "
				       << " Lateral angle difference is " << lateral_angle_difference.back() << ", "
				       << " Base track distance is " << basetrack_distance.back() << ", "
				       << " Water base track distance is " << water_basetrack_distance.back() << ","
				       << " Track tangent is " << track_tangent.back();
	      
	    }
	    
	    // Get forward initial pbeta value
	    Double_t vertex_ax = reverse_base_pair_vec.at(0).second.ax;
	    Double_t vertex_ay = reverse_base_pair_vec.at(0).second.ay;
	    TVector3 vertex_tangent(vertex_ax, vertex_ay, 1.);
	    
	    Double_t radial_angle_difference_rms = 0.;
	    Double_t lateral_angle_difference_rms = 0.;
	    for ( Int_t i = 0; i < radial_angle_difference.size(); i++ ) {
	      radial_angle_difference_rms += radial_angle_difference.at(i) * radial_angle_difference.at(i);
	    }
	    for ( Int_t i = 0; i < lateral_angle_difference.size(); i++ ) {
	      lateral_angle_difference_rms += lateral_angle_difference.at(i) * lateral_angle_difference.at(i);
	    }
	    radial_angle_difference_rms /= radial_angle_difference.size();
	    lateral_angle_difference_rms /= lateral_angle_difference.size();
	    radial_angle_difference_rms = TMath::Sqrt(radial_angle_difference_rms);
	    lateral_angle_difference_rms = TMath::Sqrt(lateral_angle_difference_rms);
	    Double_t radial_cut_value = 3. * radial_angle_difference_rms;
	    Double_t lateral_cut_value = 3. * lateral_angle_difference_rms;
	    Double_t radiation_length = CalcRadLength(ncell, vertex_tangent.Mag(), kNinjaIron);
	    Double_t initial_pbeta = MCS_SCALE_FACTOR * 13.6 / lateral_angle_difference_rms * TMath::Sqrt(radiation_length) * (1. * 0.038 * TMath::Log(radiation_length));
	    
	    // backward
	    Double_t radial_cut_value_back;
	    Double_t lateral_cut_value_back;
	    Double_t initial_pbeta_back;
	    
	    // Reconstruction
	    std::array<std::array<Double_t, 5>, 2> result_array = {};
	    for ( Int_t particle_id = 0; particle_id < kNumberOfNinjaMcsParticles; particle_id++ ) {
	      if ( particle_id == 1 ) continue;
	      for ( Int_t idirection = 0; idirection < kNumberOfNinjaMcsDirections; idirection++ ) {
		Int_t direction = MCS_DIRECTION[idirection];
		if ( direction == -1 ) continue; // tempraly
		switch ( direction ) {
		case kNinjaMcsForward :
		  result_array.at(idirection) = ReconstructPBeta(initial_pbeta, ncell, particle_id, direction,
								 radial_cut_value, lateral_cut_value,
								 radial_cut_value, lateral_cut_value,
								 kTRUE, 0, 0,
								 basetrack_distance,
								 water_basetrack_distance,
								 track_tangent,
								 plate_id, plate_id,
								 radial_angle_difference,
								 lateral_angle_difference);
		  break;
		case kNinjaMcsBackward :
		  result_array.at(idirection) = ReconstructPBeta(initial_pbeta_back, ncell, particle_id, direction,
								 radial_cut_value_back, lateral_cut_value_back,
								 radial_cut_value_back, lateral_cut_value_back,
								 kTRUE, 0, 0,
								 basetrack_distance_back,
								 water_basetrack_distance_back,
								 track_tangent_back,
								 plate_id_back, plate_id_back,
								 radial_angle_difference_back,
								 lateral_angle_difference_back);
		  break;
		default :
		  throw std::invalid_argument("Direction is not set properly");
		}
	      }

	      Double_t pbeta = result_array.at(0).at(0);
	      Double_t pbeta_err_minus = result_array.at(0).at(3);
	      Double_t pbeta_err_plus = result_array.at(0).at(2);
	      Double_t momentum_minus, momentum_plus; // p +/- 1sigma
	      if ( particle_id == 0 ) {
		mom_chain.ecc_mcs_mom[0] = CalculateMomentumFromPBeta(pbeta, MCS_MUON_MASS);
		momentum_minus = CalculateMomentumFromPBeta(pbeta + pbeta_err_minus, MCS_MUON_MASS);
		momentum_plus = CalculateMomentumFromPBeta(pbeta + pbeta_err_plus, MCS_MUON_MASS);
		mom_chain.ecc_mcs_mom_error[0][0] = mom_chain.ecc_mcs_mom[0] - momentum_minus;
		mom_chain.ecc_mcs_mom_error[0][1] = momentum_plus - mom_chain.ecc_mcs_mom[0];
	      }
	      else if ( particle_id == 2 ) {
		mom_chain.ecc_mcs_mom[1] = CalculateMomentumFromPBeta(pbeta, MCS_PROTON_MASS);
		momentum_minus = CalculateMomentumFromPBeta(pbeta + pbeta_err_minus, MCS_PROTON_MASS);
		momentum_plus = CalculateMomentumFromPBeta(pbeta + pbeta_err_plus, MCS_PROTON_MASS);
		mom_chain.ecc_mcs_mom_error[1][0] = mom_chain.ecc_mcs_mom[1] - momentum_minus;
		mom_chain.ecc_mcs_mom_error[1][1] = momentum_plus - mom_chain.ecc_mcs_mom[1];
	      }	      
	    }	    
	  }
	  ev.chains.push_back(mom_chain);
	}
      }


      // true chains
      for ( int itrue_chain = 0; itrue_chain < num_true_chain; itrue_chain++ ) {
	if ( Momentum_recon::ReadMomChainHeader(ifs, mom_chain, num_base, num_link) ) {
	  mom_chain.base.clear();
	  mom_chain.base_pair.clear();
	  mom_chain.base.reserve(num_base);
	  mom_chain.base_pair.reserve(num_link);
	  for ( int i = 0; i < num_base; i++ ) {
	    ifs.read((char*)& base, sizeof(Momentum_recon::Mom_basetrack));
	    mom_chain.base.push_back(base);
	  }
	  for ( int i = 0; i < num_link; i++ ) {
	    ifs.read((char*)& base_pair.first, sizeof(Momentum_recon::Mom_basetrack));
	    ifs.read((char*)& base_pair.second, sizeof(Momentum_recon::Mom_basetrack));
	    mom_chain.base_pair.push_back(base_pair);
	  }
	  ev.true_chains.push_back(mom_chain);
	}
      }


      // Write binary
      Momentum_recon::WriteEventInformationHeader(ofs, ev);
      // chains
      for ( auto ochain : ev.chains ) {
	Momentum_recon::WriteMomChainHeader(ofs, ochain);
	for ( auto obase : ochain.base ) {
	  ofs.write((char*)& obase, sizeof(Momentum_recon::Mom_basetrack));
	}
	for ( auto olink : ochain.base_pair ) {
	  ofs.write((char*)& olink.first, sizeof(Momentum_recon::Mom_basetrack));
	  ofs.write((char*)& olink.second, sizeof(Momentum_recon::Mom_basetrack));
	}
      }
      // true chains
      for ( auto otrue_chain : ev.true_chains ) {
	Momentum_recon::WriteMomChainHeader(ofs, otrue_chain);
	for ( auto obase : otrue_chain.base ) {
	  ofs.write((char*)& obase, sizeof(Momentum_recon::Mom_basetrack));
	}
	for ( auto olink : otrue_chain.base_pair ) {
	  ofs.write((char*)& olink.first, sizeof(Momentum_recon::Mom_basetrack));
	  ofs.write((char*)& olink.second, sizeof(Momentum_recon::Mom_basetrack));
	}	
      }

      ev.chains.clear();
      ev.true_chains.clear();

      // if (num_entry > 10) break;
    }
    
    ofs.close();

  } catch (const std::runtime_error &error) {
    BOOST_LOG_TRIVIAL(fatal) << "Runtime error : " << error.what();
    std::exit(1);
  } catch (const std::invalid_argument &error) {
    BOOST_LOG_TRIVIAL(fatal) << "Invalid argument : " << error.what();
    std::exit(1);
  }

  BOOST_LOG_TRIVIAL(info) << "==========Momentum Reconstruction Finish==========";
  std::exit(0);

}
