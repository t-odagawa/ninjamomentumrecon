// boost includes
#include <boost/log/core.hpp>
#include <boost/log/trivial.hpp>
#include <boost/log/expressions.hpp>

// system includes
#include <fstream>
#include <iostream>
#include <iomanip>
#include <sstream>
#include <vector>
#include <map>
#include <tuple>
#include <algorithm>
#include <ctime>

// my include
#include "McsClass.hpp"

namespace logging = boost::log;

int main ( int argc, char *argv[] ) {

  logging::core::get()->set_filter
    (
     logging::trivial::severity >= logging::trivial::info
     //logging::trivial::severity >= logging::trivial::debug
     );

  BOOST_LOG_TRIVIAL(info) << "==========Mom chain Separate Start==========";

  if ( argc != 2 ) {
    BOOST_LOG_TRIVIAL(error) << "Usage : " << argv[0]
			     << " <Input momch file>";
    std::exit(1);
  }

  try {

    std::ifstream ifs(argv[1], std::ios::binary);

    // filesize取得
    ifs.seekg(0, std::ios::end);
    int64_t eofpos = ifs.tellg();
    ifs.clear();
    ifs.seekg(0, std::ios::beg);
    int64_t begpos = ifs.tellg();
    int64_t nowpos = ifs.tellg();
    int64_t size2 = eofpos - begpos;
    int64_t GB = size2 / (1000 * 1000 * 1000);
    int64_t MB = (size2 - GB * 1000 * 1000 * 1000) / (1000 * 1000);
    int64_t KB = (size2 - GB * 1000 * 1000 * 1000 - MB * 1000 * 1000) / (1000);
    if (GB > 0) {
      std::cout << "FILE size :" << GB << "." << MB << " [GB]" << std::endl;
    }
    else {
      std::cout << "FILE size :" << MB << "." << KB << " [MB]" << std::endl;
    }


    std::multimap<std::tuple<int, int, int>, Momentum_recon::Event_information> event_info_map;
    // <<year, month, day>, Event_information>

    // ファイル読み込み，multimap 生成
    std::vector<Momentum_recon::Event_information> i_ev_vec = Momentum_recon::ReadEventInformationBin((std::string)argv[1]);

    for ( auto ev : i_ev_vec ) {
      time_t unixtime = (time_t)ev.unixtime;
      tm *tm_event = localtime(&unixtime);

      int year = tm_event->tm_year + 1900;
      int month = tm_event->tm_mon + 1;
      int day = tm_event->tm_mday;
      if ( ev.entry_in_daily_file == -1 ) continue;
      BOOST_LOG_TRIVIAL(debug) << year << "/" << month << "/" << day
			       << " (" << ev.entry_in_daily_file << " )";
      event_info_map.insert(std::make_pair(std::make_tuple(year, month, day), ev));
    }

    // separate file 生成
    for ( auto itr = event_info_map.begin(); itr != event_info_map.end();
	  itr = std::next(itr, event_info_map.count(itr->first)) ) {
      int year = std::get<0>(itr->first);
      int month = std::get<1>(itr->first);
      int day = std::get<2>(itr->first);

      std::stringstream daily_filename_s;

      daily_filename_s << argv[1] << "."
		       << std::to_string(year) << "_"
		       << std::to_string(month) << "_"
		       << std::to_string(day) << ".momch";

      std::vector<Momentum_recon::Event_information> daily_ev_vec;

      auto range = event_info_map.equal_range(itr->first);

      for ( auto itr2 = range.first; itr2 != range.second; itr2++ ) {
	daily_ev_vec.push_back(itr2->second);
      }

      Momentum_recon::WriteEventInformationBin(daily_filename_s.str(), daily_ev_vec);

    }

  } catch ( const std::runtime_error &error ) {
    BOOST_LOG_TRIVIAL(fatal) << "Runtime error : " << error.what();
    std::exit(1);
  } catch ( const std::invalid_argument &error ) {
    BOOST_LOG_TRIVIAL(fatal) << "Invalid argument error : " << error.what();
    std::exit(1);
  }

  BOOST_LOG_TRIVIAL(info) << "==========Mom chain Separate Finish==========";
  std::exit(0);

}
