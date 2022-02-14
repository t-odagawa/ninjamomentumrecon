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

    Momentum_recon::Mom_chain mom_chain;
    Momentum_recon::Mom_basetrack base;
    std::pair<Momentum_recon::Mom_basetrack, Momentum_recon::Mom_basetrack > base_pair;
    int num_base, num_link;
      
    const Int_t ncell = std::atoi(argv[3]);
    Int_t num_entry = 0;


    while ( Momentum_recon::ReadMomChainHeader(ifs, mom_chain, num_base, num_link)) {
      num_entry++;
      // Read binary
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

      std::vector<std::pair<Momentum_recon::Mom_basetrack, Momentum_recon::Mom_basetrack> > reverse_base_pair;
      reverse_base_pair.resize(num_link);
      std::copy(mom_chain.base_pair.begin(), mom_chain.base_pair.end(), reverse_base_pair.begin());
      std::reverse(reverse_base_pair.begin(), reverse_base_pair.end());

      // input parameters for TMinuit
      std::vector<Double_t > basetrack_distance = {};
      std::vector<Double_t > water_basetrack_distance = {};
      std::vector<Double_t > track_tangent = {};
      std::vector<Int_t > plate_id = {};
      std::vector<Double_t > radial_angle_difference = {};
      std::vector<Double_t > lateral_angle_difference = {};
      Int_t iron_plate_id = -1;
      Int_t water_plate_id = -1;

      for ( Int_t ipair = 0; ipair < reverse_base_pair.size(); ipair++ ) {
	TVector3 upstream_position;
	upstream_position.SetX(reverse_base_pair.at(ipair).second.x / 1.e3);
	upstream_position.SetY(reverse_base_pair.at(ipair).second.y / 1.e3);
	upstream_position.SetZ(reverse_base_pair.at(ipair).second.z / 1.e3);
	TVector3 downstream_position;
	downstream_position.SetX(reverse_base_pair.at(ipair).first.x / 1.e3);
	downstream_position.SetY(reverse_base_pair.at(ipair).first.y / 1.e3);
	downstream_position.SetZ(reverse_base_pair.at(ipair).first.z / 1.e3);
	TVector3 upstream_tangent;
	upstream_tangent.SetX(reverse_base_pair.at(ipair).second.ax);
	upstream_tangent.SetY(reverse_base_pair.at(ipair).second.ay);
	upstream_tangent.SetZ(1.);
	TVector3 downstream_tangent;
	downstream_tangent.SetX(reverse_base_pair.at(ipair).first.ax);
	downstream_tangent.SetY(reverse_base_pair.at(ipair).first.ay);
	downstream_tangent.SetZ(1.);

	if ( reverse_base_pair.at(ipair).second.pl - 1 <= 3 ) continue;
	if ( (reverse_base_pair.at(ipair).second.pl - 1) % 2 == 1 &&
	     reverse_base_pair.at(ipair).second.pl - 1 >= 15 ) continue;
	
	basetrack_distance.push_back((downstream_position - upstream_position).Mag() / 850.e-3);
	track_tangent.push_back(upstream_tangent.Mag());
	plate_id.push_back(reverse_base_pair.at(ipair).second.pl - 1);
	radial_angle_difference.push_back(RadialAngleDiffNew(upstream_tangent, downstream_tangent));
	lateral_angle_difference.push_back(LateralAngleDiffNew(upstream_tangent, downstream_tangent));
	iron_plate_id = reverse_base_pair.at(ipair).second.pl - 1;
	water_plate_id = reverse_base_pair.at(ipair).first.pl - 1;
	if ( ipair < reverse_base_pair.size() - 1 &&
	     reverse_base_pair.at(ipair + 1).second.pl - 1 == water_plate_id &&
	     iron_plate_id >= 18 ) {
	  TVector3 water_upstream_position;
	  water_upstream_position.SetX(reverse_base_pair.at(ipair + 1).second.x / 1.e3);
	  water_upstream_position.SetY(reverse_base_pair.at(ipair + 1).second.y / 1.e3);
	  water_upstream_position.SetZ(reverse_base_pair.at(ipair + 1).second.z / 1.e3);
	  TVector3 water_downstream_position;
	  water_downstream_position.SetX(reverse_base_pair.at(ipair + 1).first.x / 1.e3);
	  water_downstream_position.SetY(reverse_base_pair.at(ipair + 1).first.y / 1.e3);
	  water_downstream_position.SetZ(reverse_base_pair.at(ipair + 1).first.z / 1.e3);
	  water_basetrack_distance.push_back((water_downstream_position - water_upstream_position).Mag() / 2.868);
	} else {
	  water_basetrack_distance.push_back(0.);
	}
	BOOST_LOG_TRIVIAL(debug) << " Film " << reverse_base_pair.at(ipair).second.pl - 1 << " and"
				 << " Film " << reverse_base_pair.at(ipair).first.pl - 1 << ", "
				 << " Radial angle difference is " << radial_angle_difference.back() << ", "
				 << " Lateral angle difference is " << lateral_angle_difference.back() << ", "
				 << " Base track distance is " << basetrack_distance.back() << ", "
				 << " Water base track distance is " << water_basetrack_distance.back() << ","
				 << " Track tangent is " << track_tangent.back();
	
      }

      // Get initial pbeta value
      Double_t vertex_ax = reverse_base_pair.at(0).second.ax;
      Double_t vertex_ay = reverse_base_pair.at(0).second.ay;
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
      Double_t radiation_length = CalcRadLength(ncell, vertex_tangent.Mag());
      Double_t initial_pbeta = MCS_SCALE_FACTOR * 13.6 / lateral_angle_difference_rms * TMath::Sqrt(radiation_length) * (1. * 0.038 * TMath::Log(radiation_length));
      
      // Reconstruction
      Int_t direction;
      std::array<std::array<Double_t, 3>, 2> result_array = {};
      std::array<Double_t, 2> log_likelihood = {};
      Double_t recon_pbeta = -1;
      for ( Int_t particle_id = 0; particle_id < kNumberOfNinjaMcsParticles; particle_id++ ) {
	for ( Int_t idirection = 0; idirection < kNumberOfNinjaMcsDirections; idirection++ ) {
	  direction = MCS_DIRECTION[idirection];
	  result_array.at(idirection) = ReconstructPBeta(initial_pbeta, ncell, particle_id, direction,
							 radial_cut_value, lateral_cut_value,
							 kTRUE,
							 basetrack_distance,
							 water_basetrack_distance,
							 track_tangent,
							 plate_id,
							 radial_angle_difference,
							 lateral_angle_difference);

	  if ( particle_id == 0) {
	    log_likelihood.at(idirection) = FuncNegativeLogLikelihood(result_array.at(idirection).at(0),
								      ncell, particle_id, direction,
								      radial_cut_value, lateral_cut_value,
								      kTRUE,
								      basetrack_distance,
								      water_basetrack_distance,
								      track_tangent,
								      plate_id,
								      radial_angle_difference,
								      lateral_angle_difference);
	  }
	}

	//if ( particle_id == 0 ) {
	if ( particle_id == 2 ) {
	  if ( log_likelihood.at(0) <= log_likelihood.at(1) ) {
	    recon_pbeta = result_array.at(0).at(0);
	  } else {
	    recon_pbeta = result_array.at(1).at(0);
	  }
	}

      }


      // Write binary
      BOOST_LOG_TRIVIAL(debug)<< "Expected momentum : " << mom_chain.ecc_mcs_mom;
      // mom_chain.ecc_mcs_mom = CalculateMomentumFromPBeta(recon_pbeta, MCS_MUON_MASS);
      mom_chain.ecc_mcs_mom = CalculateMomentumFromPBeta(recon_pbeta, MCS_PROTON_MASS);
      BOOST_LOG_TRIVIAL(debug)<< "Reconstructed momentum : " << mom_chain.ecc_mcs_mom;
      Momentum_recon::WriteMomChainHeader(ofs, mom_chain);
      for ( auto base : mom_chain.base ) {
	ofs.write((char*)& base, sizeof(Momentum_recon::Mom_basetrack));
      }
      for ( auto link : mom_chain.base_pair ) {
	ofs.write((char*)& link.first, sizeof(Momentum_recon::Mom_basetrack));
	ofs.write((char*)& link.second, sizeof(Momentum_recon::Mom_basetrack));
      }
      // if ( num_entry > 10 ) break;
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
