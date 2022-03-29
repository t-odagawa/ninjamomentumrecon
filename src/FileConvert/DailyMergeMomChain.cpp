// boost includes
#include <boost/log/core.hpp>
#include <boost/log/trivial.hpp>
#include <boost/log/expressions.hpp>
#include <boost/filesystem.hpp>
#include <boost/range/iterator_range.hpp>

// system includes
#include <fstream>
#include <iostream>
#include <sstream>
#include <vector>
#include <algorithm>

// my include
#include "McsClass.hpp"

namespace logging = boost::log;
namespace fs = boost::filesystem;

bool CompareFileName(const fs::path &lhs, const fs::path &rhs);

int main ( int argc, char *argv[] ) {

  logging::core::get()->set_filter
    (
     logging::trivial::severity >= logging::trivial::info
     //logging::trivial::severity >= logging::trivial::debug
     );

  BOOST_LOG_TRIVIAL(info) << "==========Mom chain Merge Start==========";

  if ( argc != 2 ) {
    BOOST_LOG_TRIVIAL(error) << "Usage : " << argv[0]
			     << " <Input momch file prefix>";
    std::exit(1);
  }

  try {

    std::string ifile_prefix = argv[1];
    fs::path prefix(ifile_prefix);
    BOOST_LOG_TRIVIAL(info) << "Directory : " << prefix.parent_path() << ", "
			    << "Filename prefix : " << prefix.filename();

    fs::path directory_path(prefix.parent_path());
    std::vector<fs::path> daily_filename_vector;
    for ( const auto &daily_file_path : boost::make_iterator_range( fs::directory_iterator(directory_path), {})) {
      if (!fs::is_directory(daily_file_path) && 
	  daily_file_path.path().filename() != prefix.filename())
	daily_filename_vector.push_back(daily_file_path.path().filename());
    }

    std::sort(daily_filename_vector.begin(), daily_filename_vector.end(),
	      CompareFileName);

    std::vector<Momentum_recon::Event_information> ev_vec;

    for ( auto daily_filename : daily_filename_vector ) {
      std::vector<Momentum_recon::Event_information > daily_ev_vec = Momentum_recon::ReadEventInformationBin(daily_filename.generic_string());

      ev_vec.insert(ev_vec.end(), daily_ev_vec.begin(), daily_ev_vec.end());

    }

    fs::path merge_mom_file_path(ifile_prefix);
    merge_mom_file_path += ".merge.momch";
    BOOST_LOG_TRIVIAL(info) << "Output file path : " << merge_mom_file_path.generic_string();

    Momentum_recon::WriteEventInformationBin(merge_mom_file_path.generic_string(), ev_vec);

  } catch ( const std::runtime_error &error ) {
    BOOST_LOG_TRIVIAL(fatal) << "Runtime error : " << error.what();
    std::exit(1);
  } catch ( const std::invalid_argument &error ) {
    BOOST_LOG_TRIVIAL(fatal) << "Invalid argument error : " << error.what();
    std::exit(1);
  }

  BOOST_LOG_TRIVIAL(info) << "==========Mom chain Merge Finish==========";
  std::exit(0);

}

bool CompareFileName(const fs::path &lhs, const fs::path &rhs) {
  std::string prefix = "muon_time_ecc5.momch.";
  // hard coding...
  std::string l_filename = lhs.generic_string();
  int l_year = std::atoi(l_filename.substr(prefix.length(), 4).c_str());
  int l_month; 
  int l_month_length;
  if ( l_filename.substr(prefix.length()+4+1, 2).find("_") == std::string::npos )
    l_month_length = 2;
  else
    l_month_length = 1;
  l_month = std::atoi(l_filename.substr(prefix.length()+4+1, l_month_length).c_str());
  int l_day;
  int l_day_length;
  if ( l_filename.substr(prefix.length()+4+1+l_month_length+1, 2).find("_") == std::string::npos )
    l_day_length = 2;
  else
    l_day_length = 1;
  l_day = std::atoi(l_filename.substr(prefix.length()+4+1+l_month_length+1, l_day_length).c_str());

  std::string r_filename = rhs.generic_string();
  int r_year = std::atoi(r_filename.substr(prefix.length(), 4).c_str());
  int r_month; 
  int r_month_length;
  if ( r_filename.substr(prefix.length()+4+1, 2).find("_") == std::string::npos )
    r_month_length = 2;
  else
    r_month_length = 1;
  r_month = std::atoi(r_filename.substr(prefix.length()+4+1, r_month_length).c_str());
  int r_day;
  int r_day_length;
  if ( r_filename.substr(prefix.length()+4+1+r_month_length+1, 2).find("_") == std::string::npos )
    r_day_length = 2;
  else
    r_day_length = 1;
  r_day = std::atoi(r_filename.substr(prefix.length()+4+1+r_month_length+1, r_day_length).c_str());

  if (l_year != r_year) return l_year < r_year;
  else if (l_month != r_month) return l_month < r_month;
  return l_day < r_day;
}
