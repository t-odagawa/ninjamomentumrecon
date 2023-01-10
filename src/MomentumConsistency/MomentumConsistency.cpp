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

namespace logging = boost::log;

int main (int argc, char* argv[]) {

  logging::core::get()->set_filter
    (
     logging::trivial::severity >= logging::trivial::debug
     );

  BOOST_LOG_TRIVIAL(info) << "==========Muon momentum consistency check==========";

  if ( argc != 3 ) {
    BOOST_LOG_TRIVIAL(error) << "Usage : " << argv[0]
			     << " <input momch file name> <output momch file name>";
    std::exit(1);
  }

  try {
    
    std::string ifilename = argv[1];
    auto iev_vec = Momentum_recon::ReadEventInformationBin(ifilename);

    std::vector<Momentum_recon::Event_information> oev_vec;
    oev_vec.reserve(iev_vec.size());

    for ( auto &ev : iev_vec ) {

      bool mom_consistency_flag = true;

      for ( auto &chain : ev.chains ) {
	if ( chain.particle_flag % 10000 == 13 ) {
	  double mcs_mom = chain.ecc_mcs_mom[0];
	  double range_mom = chain.bm_range_mom;

	  if ( chain.stop_flag == 0 ) {
	    if ( range_mom < mcs_mom ) break;
	    else if ( range_mom > 1000. ) break;
	    else {
	      double sigma = (range_mom - mcs_mom) / std::hypot(chain.ecc_mcs_mom_error[0][1],
								chain.bm_range_mom_error[0]);
	      if ( sigma > 2.5 ) mom_consistency_flag = false;
	      break;
	    }
	  }
	  else if ( chain.stop_flag == 1 ) {
	    if ( range_mom > 1000. ) break;
	    else if ( range_mom < mcs_mom ){
	      double sigma = (mcs_mom - range_mom) / std::hypot(chain.ecc_mcs_mom_error[0][0],
								chain.bm_range_mom_error[1]);
	      if ( sigma > 2.5 ) mom_consistency_flag = false;
	      break;
	    }
	    else {
	      double sigma = (range_mom - mcs_mom) / std::hypot(chain.ecc_mcs_mom_error[0][1],
								chain.bm_range_mom_error[0]);
	      if ( sigma > 2.5 ) mom_consistency_flag = false;
	      break;
	    }
	  }
	}
      }

      if ( mom_consistency_flag ) oev_vec.push_back(ev);

    }

    std::string ofilename = argv[2];
    Momentum_recon::WriteEventInformationBin(ofilename, oev_vec);


  } catch ( const std::runtime_error &error ) {
    BOOST_LOG_TRIVIAL(fatal) << "Runtime error :" << error.what();
    std::exit(1);
  }

  std::exit(0);

}
