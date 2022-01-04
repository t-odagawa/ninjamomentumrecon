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

    Momentum_recon::Mom_chain mom_chain;
    Momentum_recon::Mom_basetrack mom_basetrack;
    std::pair<Momentum_recon::Mom_basetrack, Momentum_recon::Mom_basetrack> mom_basetrack_pair;

    std::multimap<std::tuple<int, int, int>, Momentum_recon::Mom_chain > mom_chain_map;
    // <<year, month, day>, mom_chain>

    // ファイル読み込み，multimap 生成
    int64_t num_entry = 0;

    int num_base, num_link;
    while ( Momentum_recon::ReadMomChainHeader(ifs, mom_chain, num_base, num_link) ) {

      if (num_entry % 100 == 0) {
	nowpos = ifs.tellg();
	auto size1 = nowpos - begpos;
	std::cerr << std::right << std::fixed << "\r now reading ..." << std::setw(4) << std::setprecision(1) << size1 * 100. / size2 << "%";
      }
      num_entry++;

      mom_chain.base.clear();
      mom_chain.base_pair.clear();
      mom_chain.base.reserve(num_base);
      mom_chain.base_pair.reserve(num_link);
      for ( int i = 0; i < num_base; i++ ) {
	ifs.read((char*)& mom_basetrack, sizeof(Momentum_recon::Mom_basetrack));
	mom_chain.base.push_back(mom_basetrack);
      }
      for ( int i = 0; i < num_link; i++ ) {
	ifs.read((char*)& mom_basetrack_pair.first, sizeof(Momentum_recon::Mom_basetrack));
	ifs.read((char*)& mom_basetrack_pair.second, sizeof(Momentum_recon::Mom_basetrack));
	mom_chain.base_pair.push_back(mom_basetrack_pair);
      }

      time_t unixtime = (time_t)mom_chain.unixtime;
      tm *tm_event = localtime(&unixtime);

      int year = tm_event->tm_year + 1900;
      int month = tm_event->tm_mon + 1;
      int day = tm_event->tm_mday;
      if ( mom_chain.entry_in_daily_file == -1 ) continue;
      BOOST_LOG_TRIVIAL(debug) << year << "/" << month << "/" << day
			       << " (" << mom_chain.entry_in_daily_file << " )";
      mom_chain_map.insert(std::make_pair(std::make_tuple(year, month, day), mom_chain));

    }

    auto size1 = eofpos - begpos;
    std::cerr << "\r now reading ..." << std::setw(4) << std::setprecision(1) << size1 * 100. / size2 << "%" << std::endl;

    // separate file 生成
    for ( auto itr = mom_chain_map.begin(); itr != mom_chain_map.end();
	  itr = std::next(itr, mom_chain_map.count(itr->first)) ) {
      int year = std::get<0>(itr->first);
      int month = std::get<1>(itr->first);
      int day = std::get<2>(itr->first);

      std::stringstream filename_s;

      filename_s << argv[1] << "."
		 << std::to_string(year) << "_"
		 << std::to_string(month) << "_"
		 << std::to_string(day) << ".momch";
      std::string filename = filename_s.str();
      std::ofstream daily_file(filename, std::ios::binary);

      auto range = mom_chain_map.equal_range(itr->first);

      for ( auto itr2 = range.first; itr2 != range.second; itr2++ ) {

	WriteMomChainHeader(daily_file, itr2->second);

	for ( int i = 0; i < itr2->second.base.size(); i++ ) {
	  daily_file.write((char*)& itr2->second.base.at(i), sizeof(Momentum_recon::Mom_basetrack));
	}
	for ( int i = 0; i < itr2->second.base_pair.size(); i++ ) {
	  daily_file.write((char*)& itr2->second.base_pair.at(i).first, sizeof(Momentum_recon::Mom_basetrack));
	  daily_file.write((char*)& itr2->second.base_pair.at(i).second, sizeof(Momentum_recon::Mom_basetrack));
	}
      }

      daily_file.close();

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
