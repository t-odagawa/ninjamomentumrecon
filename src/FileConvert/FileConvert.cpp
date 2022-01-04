// boost includes
#include <boost/log/core.hpp>
#include <boost/log/trivial.hpp>
#include <boost/log/expressions.hpp>

// B2 includes
#include <B2Reader.hh>
#include <B2Writer.hh>
#include <B2SpillSummary.hh>
#include <B2EventSummary.hh>
#include <B2EmulsionSummary.hh>

// ROOT includes

// system includes
#include <string>
#include <fstream>
#include <iostream>
#include <vector>
#include <algorithm>
#include <map>
#include <utility>

// my include
#include "McsClass.hpp"
#include "McsFunction.hpp"

namespace logging = boost::log;
namespace fs = boost::filesystem;


int main (int argc, char *argv[]) {
  
  logging::core::get()->set_filter
    (
     logging::trivial::severity >= logging::trivial::info
     //logging::trivial::severity >= logging::trivial::debug
     );
  
  BOOST_LOG_TRIVIAL(info) << "==========MCS Convert Start==========";
  
  if ( argc != 5 ) {
    BOOST_LOG_TRIVIAL(error) << "Usage : " << argv[0]
			     << " <Input B2 file> <Input momch file> <Output B2 file> <ECC id (0-8)> ";
    std::exit(1);
  }
  
  try {
    
    B2Reader reader((std::string)argv[1]);

    std::ifstream ifs(argv[2], std::ios::binary);

    Momentum_recon::Mom_chain mom_chain;
    Momentum_recon::Mom_basetrack mom_basetrack;
    std::pair<Momentum_recon::Mom_basetrack, Momentum_recon::Mom_basetrack> mom_basetrack_pair;

    B2Writer writer((std::string)argv[3], reader);

    int ecc_id = std::atoi(argv[4]);

    int num_entry = 0;
    int num_base, num_link;

    std::multimap<int, Momentum_recon::Mom_chain > entry_chain_map;

    while ( Momentum_recon::ReadMomChainHeader(ifs, mom_chain, num_base, num_link) ) {

      mom_chain.base.clear();
      mom_chain.base_pair.clear();
      mom_chain.base.reserve(num_base);
      mom_chain.base_pair.reserve(num_link);
      BOOST_LOG_TRIVIAL(debug) << "\nGroup ID : " << mom_chain.groupid
			       << "\nChain ID : " << mom_chain.chainid
			       << "\nUnixtime : " << mom_chain.unixtime
			       << "\nTraker track ID : " << mom_chain.tracker_track_id
			       << "\nEntry in daily file : " << mom_chain.entry_in_daily_file
			       << "\nNumber of basetracks : " << num_base
			       << "\nNumber of likelets : " << num_link;
      BOOST_LOG_TRIVIAL(debug) << "Getting basetrack information";
      for ( int i = 0; i < num_base; i++ ) {
	ifs.read((char*)& mom_basetrack, sizeof(Momentum_recon::Mom_basetrack));
	mom_chain.base.push_back(mom_basetrack);
      }
      BOOST_LOG_TRIVIAL(debug) << "Getting linklet information";
      for ( int i = 0; i < num_link; i++ ) {
	ifs.read((char*)& mom_basetrack_pair.first, sizeof(Momentum_recon::Mom_basetrack));
	ifs.read((char*)& mom_basetrack_pair.second, sizeof(Momentum_recon::Mom_basetrack));
	mom_chain.base_pair.push_back(mom_basetrack_pair);
      }      

      entry_chain_map.emplace(mom_chain.entry_in_daily_file, mom_chain);

    }

    BOOST_LOG_TRIVIAL(info) << "Multimap size : " << entry_chain_map.size();

    // Vectors for basetrack info
    std::vector<Int_t > rawid_vec;
    std::vector<Int_t > plate_vec;
    std::vector<TVector3 > film_position_vec;
    std::vector<TVector3 > tangent_vec;
    std::vector<Double_t > vph_up_vec, vph_down_vec;
    std::vector<Double_t > pixel_count_up_vec, pixel_count_down_vec;
    
    std::vector<TVector3 > absolute_position_vec;
    std::vector<TVector3 > film_position_in_down_coordinate_vec;
    std::vector<TVector3 > tangent_in_down_coordinate_vec;
    while ( reader.ReadNextSpill() > 0 ) {
      auto &input_spill_summary = reader.GetSpillSummary();
      auto &output_spill_summary = writer.GetSpillSummary();

      if ( entry_chain_map.count(reader.GetEntryNumber()) > 0 ) {

	BOOST_LOG_TRIVIAL(info) << "Entry : " << reader.GetEntryNumber() << ", "
				<< "Chain : " << entry_chain_map.count(reader.GetEntryNumber());

	auto p = entry_chain_map.equal_range(reader.GetEntryNumber());
	for ( auto itr = p.first; itr != p.second; itr++ ) {
	  auto single_mom_chain = itr->second;
	  int single_num_base = itr->second.base.size();
	  int single_num_link = itr->second.base_pair.size();
	  BOOST_LOG_TRIVIAL(debug) << "Entry in chain : " << single_mom_chain.entry_in_daily_file << ", "
				   << "Number of basetracks : " << single_num_base << ", "
				   << "Number of linklets : " << single_num_link;
	  rawid_vec.clear(); rawid_vec.resize(single_num_base);
	  plate_vec.clear(); plate_vec.resize(single_num_base);
	  film_position_vec.clear(); film_position_vec.resize(single_num_base);
	  tangent_vec.clear(); tangent_vec.resize(single_num_base);
	  vph_up_vec.clear(); vph_down_vec.clear();
	  vph_up_vec.resize(single_num_base); vph_down_vec.resize(single_num_base);
	  pixel_count_up_vec.clear(); pixel_count_down_vec.clear();
	  pixel_count_up_vec.resize(single_num_base); pixel_count_down_vec.resize(single_num_base);

	  absolute_position_vec.clear();
	  absolute_position_vec.resize(single_num_base);
	  film_position_in_down_coordinate_vec.clear();
	  film_position_in_down_coordinate_vec.resize(single_num_base);
	  tangent_in_down_coordinate_vec.clear();
	  tangent_in_down_coordinate_vec.resize(single_num_base);
      
	  TVector3 absolute_position, film_position, tangent;
	  TVector3 film_position_in_down_coordinate, tangent_in_down_coordinate;
	  Double_t vph[2];
      
	  for ( int ibase = 0; ibase < single_num_base; ibase++ ) {
	    absolute_position.SetX(single_mom_chain.base.at(ibase).x / 1e3);
	    absolute_position.SetY(single_mom_chain.base.at(ibase).y / 1e3);
	    absolute_position.SetZ(single_mom_chain.base.at(ibase).z / 1e3);
	    // offset
	    PositionAddOffset(absolute_position, ecc_id);
	    
	    tangent.SetX(single_mom_chain.base.at(ibase).ax);
	    tangent.SetY(single_mom_chain.base.at(ibase).ay);
	    tangent.SetZ(1.);
	    
	    vph[0] = (Double_t)(single_mom_chain.base.at(ibase).m[0].ph % 10000);
	    // vph[0] = (Double_t)(single_mom_chain.base.at(ibase).m[0].ph);
	    vph[1] = (Double_t)(single_mom_chain.base.at(ibase).m[1].ph % 10000);
	    // vph[1] = (Double_t)(single_mom_chain.base.at(ibase).m[1].ph);
	    rawid_vec.at(ibase) = single_mom_chain.base.at(ibase).rawid;
	    plate_vec.at(ibase) = single_mom_chain.base.at(ibase).pl - 1;
	    absolute_position_vec.at(ibase) = absolute_position;
	    tangent_vec.at(ibase) = tangent;
	    vph_up_vec.at(ibase) = vph[0];
	    vph_down_vec.at(ibase) = vph[1];
	    pixel_count_up_vec.at(ibase) = single_mom_chain.base.at(ibase).m[0].hitnum;
	    pixel_count_down_vec.at(ibase) = single_mom_chain.base.at(ibase).m[1].hitnum;	
	  }
      
	  // Linklet information write
	  for ( int ilink = 0; ilink < single_num_link; ilink++ ) {
	    film_position.SetX(single_mom_chain.base_pair.at(ilink).first.x / 1e3);
	    film_position.SetY(single_mom_chain.base_pair.at(ilink).first.y / 1e3);
	    film_position.SetZ(single_mom_chain.base_pair.at(ilink).first.z / 1e3);
	    film_position_in_down_coordinate.SetX(single_mom_chain.base_pair.at(ilink).second.x / 1e3);
	    film_position_in_down_coordinate.SetY(single_mom_chain.base_pair.at(ilink).second.y / 1e3);
	    film_position_in_down_coordinate.SetZ(single_mom_chain.base_pair.at(ilink).second.z / 1e3);
	    tangent.SetX(single_mom_chain.base_pair.at(ilink).first.ax);
	    tangent.SetY(single_mom_chain.base_pair.at(ilink).first.ay);
	    tangent.SetZ(1.);
	    tangent_in_down_coordinate.SetX(single_mom_chain.base_pair.at(ilink).second.ax);
	    tangent_in_down_coordinate.SetY(single_mom_chain.base_pair.at(ilink).second.ay);
	    tangent_in_down_coordinate.SetZ(1.);
	    auto find_itr = std::find(rawid_vec.begin(), rawid_vec.end(),
				      single_mom_chain.base_pair.at(ilink).second.rawid);
	    if ( find_itr != rawid_vec.end() ) {
	      int index = std::distance(rawid_vec.begin(), find_itr);
	      film_position_vec.at(index) = film_position;
	      film_position_in_down_coordinate_vec.at(index) = film_position_in_down_coordinate;
	      tangent_vec.at(index) = tangent;
	      tangent_in_down_coordinate_vec.at(index) = tangent_in_down_coordinate;
	    }
	  }
	  
	  std::reverse(rawid_vec.begin(), rawid_vec.end());
	  std::reverse(plate_vec.begin(), plate_vec.end());
	  std::reverse(absolute_position_vec.begin(), absolute_position_vec.end());
	  std::reverse(film_position_vec.begin(), film_position_vec.end());
	  std::reverse(tangent_vec.begin(), tangent_vec.end());
	  std::reverse(vph_up_vec.begin(), vph_up_vec.end());
	  std::reverse(vph_down_vec.begin(), vph_down_vec.end());
	  std::reverse(pixel_count_up_vec.begin(), pixel_count_up_vec.end());
	  std::reverse(pixel_count_down_vec.begin(), pixel_count_down_vec.end());
	  std::reverse(film_position_in_down_coordinate_vec.begin(),
		       film_position_in_down_coordinate_vec.end());
	  std::reverse(tangent_in_down_coordinate_vec.begin(),
		       tangent_in_down_coordinate_vec.end());
	  
	  for ( int ibase = 0; ibase < single_num_base; ibase++ ) {
	    auto &emulsion_summary = output_spill_summary.AddEmulsion();
	    emulsion_summary.SetEmulsionTrackId((UInt_t)rawid_vec.at(ibase));
	    emulsion_summary.SetParentTrackId(single_mom_chain.groupid * 10 + single_mom_chain.chainid);
	    emulsion_summary.SetAbsolutePosition(absolute_position_vec.at(ibase));
	    emulsion_summary.SetFilmPosition(film_position_vec.at(ibase));
	    emulsion_summary.SetTangent(tangent_vec.at(ibase));
	    emulsion_summary.SetFilmPositionInDownCoordinate(film_position_in_down_coordinate_vec.at(ibase));
	    emulsion_summary.SetTangentInDownCoordinate(tangent_in_down_coordinate_vec.at(ibase));
	    emulsion_summary.SetVphUp(vph_up_vec.at(ibase));
	    emulsion_summary.SetVphDown(vph_down_vec.at(ibase));
	    emulsion_summary.SetPixelCountUp(pixel_count_up_vec.at(ibase));
	    emulsion_summary.SetPixelCountDown(pixel_count_down_vec.at(ibase));
	    emulsion_summary.SetFilmType(B2EmulsionType::kECC);
	    emulsion_summary.SetEcc(ecc_id);
	    emulsion_summary.SetPlate(plate_vec.at(ibase));
	    emulsion_summary.SetMuonTrackId(single_mom_chain.tracker_track_id);
	  } // base
	} // multimap loop
      } // fi
      writer.Fill();
    }

  } catch (const std::runtime_error &error) {
    BOOST_LOG_TRIVIAL(fatal) << "Runtime error : " << error.what();
    std::exit(1);
  } catch (const std::invalid_argument &error) {
    BOOST_LOG_TRIVIAL(fatal) << "Invalid argument error : " << error.what();
    std::exit(1);
  }
    
  BOOST_LOG_TRIVIAL(info) << "==========MCS Convert Finish==========";
  std::exit(0);

}
    
