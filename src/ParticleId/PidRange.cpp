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

// my includes
#include "McsClass.hpp"
#include "McsConst.hpp"
#include "McsFunction.hpp"
#include "PidData.hpp"
#include "PidClass.hpp"
#include "PidFunction.hpp"
#include "RangeSpline.hpp"
#include "RangeFunction.hpp"

namespace logging = boost::log;

int main ( int argc, char* argv[]) {
  
  logging::core::get()->set_filter
    (
     logging::trivial::severity >= logging::trivial::debug
     //logging::trivial::severity >= logging::trivial::trace
     );

  BOOST_LOG_TRIVIAL(info) << "==========PID/Range measurement start==========";

  if ( argc != 5 ) {
    BOOST_LOG_TRIVIAL(error) << "Usage : " << argv[0]
			     << " <input momch file name> <output momch file name> <data dir path> <seed>";
    std::exit(1);
  }

  try {

    std::string ifilename = argv[1];
    std::vector<Momentum_recon::Event_information > ev_vec = Momentum_recon::ReadEventInformationBin(ifilename);

    const std::string data_dir_path = argv[3];
    const long seed = std::atoi(argv[4]);
    const PidData pid_data(data_dir_path);
    const PidFunction pid_function(pid_data, seed);
    const RangeSpline range_spline(data_dir_path);
    const RangeFunction range_function(range_spline);

    if ( seed == 0 )
      gRandom->SetSeed(time(NULL));
    else
      gRandom->SetSeed(seed);

    // pid_function.CheckMeanSigmaValues();

    Int_t num_entry = 0;
    
    for ( auto &ev : ev_vec ) {
      
      BOOST_LOG_TRIVIAL(debug) << "Entry : " << num_entry << ", " << ev.groupid;

      num_entry++;

      if ( ev.chains.empty() ) continue;

      for ( auto &chain : ev.chains ) {

	if ( chain.base_pair.empty() ) continue;

	// if ( chain.direction != 1 ) continue; // debug のため

	double vph = 0.;
	std::vector<double> tangent(3);
	if ( chain.direction == 1 ) {
	  tangent.at(0) = chain.base.back().ax;
	  tangent.at(1) = chain.base.back().ay;
	}
	else if ( chain.direction == -1 ) {
	  tangent.at(0) = chain.base.front().ax;
	  tangent.at(1) = chain.base.front().ay;
	}
	tangent.at(2) = 1.;

	// particle id を確認，muon はすでに終わっている
	int recon_particle_id = chain.particle_flag % 10000;
	int true_particle_id = chain.particle_flag / 10000;

	// VPH を計算 (data-driven, つまり validation ではない)
	double pb_mu_assume = CalculatePBetaFromMomentum(chain.ecc_mcs_mom[0], MCS_MUON_MASS);
	double tangent_xy = std::hypot(tangent.at(0), tangent.at(1));	
	vph = pid_function.GetVph(true_particle_id, pb_mu_assume, tangent_xy);
	// vph = pid_function.GetVph(true_particle_id, pb_mu_assume, tangent_xy, -1, 0, 1, 0); // mean plus
	// vph = pid_function.GetVph(true_particle_id, pb_mu_assume, tangent_xy, 1, 0, -1, 0); // mean minus
	// vph = pid_function.GetVph(true_particle_id, pb_mu_assume, tangent_xy, 0, 1, 0, 1); // sigma plus
	// vph = pid_function.GetVph(true_particle_id, pb_mu_assume, tangent_xy, 0, -1, 0, -1); // sigma minus

	chain.base.front().m[0].ph = vph;
	
	// p/pi likelihood を計算
	pid_function.CalcPartnerLikelihood(vph, pb_mu_assume,
					   std::hypot(tangent.at(0), tangent.at(1)),
					   chain.muon_likelihood, chain.proton_likelihood);

	if ( recon_particle_id != 0 ) continue;

	// likelihood に基づき recon particle id を決定
	if ( pb_mu_assume < 700. ) {
	  recon_particle_id = pid_function.GetReconPid(vph, pb_mu_assume, std::hypot(tangent.at(0), tangent.at(1)),
						       chain.muon_likelihood, chain.proton_likelihood);
	  chain.particle_flag += recon_particle_id;
	  BOOST_LOG_TRIVIAL(trace) << "RECON VPH : " << vph
				   << "True : " << true_particle_id << ", Recon :" << recon_particle_id
				   << "PBETA : " << pb_mu_assume << ", ANGLE : " << std::hypot(tangent.at(0), tangent.at(1))
				   << "LIKELIHOOD : " << chain.muon_likelihood << ", " << chain.proton_likelihood;
	}

	// ECC 内で partner が止まっているかを確認
	pid_function.CalculateStopFlag(chain, ev.true_chains);

	// Range momentum を計算

	std::vector<double > ax, ay;
	std::vector<int > pl;
	ax.reserve(chain.base.size());
	ay.reserve(chain.base.size());
	pl.reserve(chain.base.size());
	
	if ( chain.direction == 1 ) {
	  for ( auto pair : chain.base_pair ) {
	    double x_diff = pair.first.x - pair.second.x;
	    double y_diff = pair.first.y - pair.second.y;
	    double z_diff = pair.first.z - pair.second.z;
	    ax.push_back(x_diff / z_diff);
	    ay.push_back(y_diff / z_diff);
	    pl.push_back(pair.first.pl);
	  }
	  ax.push_back(chain.base.back().ax);
	  ay.push_back(chain.base.back().ay);
	  pl.push_back(chain.base.back().pl);
	}
	else if ( chain.direction == -1 ) {
	  ax.push_back(chain.base.front().ax);
	  ay.push_back(chain.base.front().ay);
	  pl.push_back(chain.base.front().pl);
	  for ( auto pair : chain.base_pair ) {
	    double x_diff = pair.first.x - pair.second.x;
	    double y_diff = pair.first.y - pair.second.y;
	    double z_diff = pair.first.z - pair.second.z;
	    ax.push_back(x_diff / z_diff);
	    ay.push_back(y_diff / z_diff);
	    pl.push_back(pair.second.pl);
	  }	  
	}
	
	range_function.ModifyVectors(ax, ay, pl);
	chain.ecc_range_mom[0] = range_function.CalculateEnergyFromRange(ax, ay, pl,
									 211, chain.direction);
	chain.ecc_range_mom[0] += MCS_PION_MASS;
	chain.ecc_range_mom[0] = CalculateMomentumFromEnergy(chain.ecc_range_mom[0],
							     MCS_PION_MASS);
	chain.ecc_range_mom[1] = range_function.CalculateEnergyFromRange(ax, ay, pl,
									 2212, chain.direction);
	chain.ecc_range_mom[1] += MCS_PROTON_MASS;
	chain.ecc_range_mom[1] = CalculateMomentumFromEnergy(chain.ecc_range_mom[1],
							     MCS_PROTON_MASS);
	chain.ecc_range_mom_error[1][0] = range_function.CalculateProtonRangeError(chain.ecc_range_mom[1],
										   std::hypot(tangent.at(0), tangent.at(1)));
	chain.ecc_range_mom_error[1][1] = chain.ecc_range_mom_error[1][0];
      }

    }
    
    std::string ofilename = argv[2];
    Momentum_recon::WriteEventInformationBin(ofilename, ev_vec);

  } catch ( const std::runtime_error &error ) {
    BOOST_LOG_TRIVIAL(fatal) << "Runtime error :" << error.what();
    std::exit(1);
  }

  BOOST_LOG_TRIVIAL(info) << "==========PID/Range measurement finish==========";
  std::exit(0);

}
