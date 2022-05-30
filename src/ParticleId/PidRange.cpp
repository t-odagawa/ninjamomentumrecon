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

  if ( argc != 4 ) {
    BOOST_LOG_TRIVIAL(error) << "Usage : " << argv[0]
			     << " <input momch file name> <output momch file name> <data dir path>";
    std::exit(1);
  }

  try {

    std::string ifilename = argv[1];
    std::vector<Momentum_recon::Event_information > ev_vec = Momentum_recon::ReadEventInformationBin(ifilename);

    const std::string data_dir_path = argv[3];
    const PidData pid_data(data_dir_path);
    const PidFunction pid_function(pid_data);
    const RangeSpline range_spline(data_dir_path);
    const RangeFunction range_function(range_spline);

    Int_t num_entry = 0;
    
    for ( auto &ev : ev_vec ) {
      
      BOOST_LOG_TRIVIAL(debug) << "Entry : " << num_entry << ", " << ev.groupid;

      num_entry++;

      if ( ev.chains.empty() ) continue;

      for ( auto &chain : ev.chains ) {

	if ( chain.base_pair.empty() ) continue;

	if ( chain.direction != 1 ) continue; // debug のため

	double vph = 0.;
	std::vector<double> tangent(3);
	if ( chain.direction == 1 ) {
	  tangent.at(0) = chain.base.front().ax;
	  tangent.at(1) = chain.base.front().ay;
	}
	else if ( chain.direction == -1 ) {
	  tangent.at(0) = chain.base.back().ax;
	  tangent.at(1) = chain.base.back().ay;
	}
	tangent.at(2) = 1.;

	// particle id を確認，muon はすでに終わっている
	int true_particle_id = chain.particle_flag % 10000;

	if ( true_particle_id != 2212 &&
	     std::abs(true_particle_id) != 211) continue;

	// VPH を計算 (data-driven, つまり validation ではない)
	vph = pid_function.GetVph(true_particle_id, chain.ecc_mcs_mom[0], 
				  std::hypot(tangent.at(0), tangent.at(1)));
	chain.base.front().m[0].ph = vph;
	
	// p/pi likelihood を計算
	pid_function.CalcPartnerLikelihood(vph, chain.ecc_mcs_mom[0],
					   std::hypot(tangent.at(0), tangent.at(1)),
					   chain.muon_likelihood, chain.proton_likelihood);

	// likelihood に基づき recon particle id を決定
	int recon_particle_id = pid_function.GetReconPid(chain.muon_likelihood, chain.proton_likelihood);
	chain.particle_flag += recon_particle_id * 10000;

	// ECC 内で partner が止まっているかを確認
	
	
	// Range momentum を計算

	std::vector<double > ax, ay;
	std::vector<int > pl;
	ax.reserve(chain.base.size());
	ay.reserve(chain.base.size());
	pl.reserve(chain.base.size());
	for ( auto base : chain.base ) {
	  ax.push_back(base.ax);
	  ay.push_back(base.ay);
	  pl.push_back(base.pl);
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
