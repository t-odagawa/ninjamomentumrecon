// boost includes
#include <boost/log/core.hpp>
#include <boost/log/trivial.hpp>
#include <boost/log/expressions.hpp>

// B2 includes
#include <string>
#include <fstream>
#include <iostream>
#include <vector>
#include <algorithm>
#include <utility>

// ROOT includes
#include <TFile.h>
#include <TString.h>
#include <TH1D.h>
#include <TH2D.h>

// my includes
#include "McsClass.hpp"

namespace logging = boost::log;
namespace fs = boost::filesystem;

int main ( int argc, char *argv[] ) {

  logging::core::get()->set_filter
    (
     logging::trivial::severity >= logging::trivial::debug
     );

  BOOST_LOG_TRIVIAL(info) << "==========Data momch to plots start==========";

  if ( argc != 3 ) {
    BOOST_LOG_TRIVIAL(error) << "Usage : " << argv[0]
			     << " <Input momch file> <Output plot file>";
    std::exit(1);
  }
  
  try {
    
    std::string ifilename = argv[1];
    if ( !fs::exists(ifilename) ) {
      throw std::runtime_error("File not exsits : " + ifilename);
    }
    std::ifstream ifs(ifilename);
    std::vector<Momentum_recon::Event_information > ev_vec = Momentum_recon::ReadEventInformationBin(momchfilename);
    if ( ev_vec.empty() ) {
      throw std::runtime_error("Event data empty");
    }


    

  } catch (const std::runtime_error &error) {
    BOOST_LOG_TRIVIAL(fatal) << "Runtime error : " << error.what();
    std::exit(1);
  } catch (const std::invalid_argument &error) {
    BOOST_LOG_TRIVIAL(fatal) << "Invalid argument : " << error.what();
    std::exit(1);
  }

  BOOST_LOG_TRIVIAL(info) << "==========Data momch to plots finish==========";
  std::exit(0);

}
