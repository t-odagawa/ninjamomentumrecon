// boost includes
#include <boost/log/core.hpp>
#include <boost/log/trivial.hpp>
#include <boost/log/expressions.hpp>

// B2 includes
#include <B2Reader.hh>
#include <B2Writer.hh>
#include <B2SpillSummary.hh>
#include <B2EmulsionSummary.hh>

// ROOT includes

// system includes
#include <fstream>
#include <iostream>
#include <vector>
#include <algorithm>
#include <ctime>

// my include
#include "FileConvert.hpp"

namespace logging = boost::log;

int main (int argc, char *argv[]) {
  
  logging::core::get()->set_filter
    (
     logging::trivial::severity >= logging::trivial::info
     //logging::trivial::severity >= logging::trivial::debug
     );
  
  BOOST_LOG_TRIVIAL(info) << "==========MCS Convert Start==========";
  
  if (argc != 4) {
    BOOST_LOG_TRIVIAL(error) << "Usage : " << argv[0]
			     << " <Input B2 file> <Input text file> <Output B2 file>";
    std::exit(1);
  }
  
  try {
    
    B2Reader reader(argv[1]);
#ifdef TEXT_MODE    
    std::ifstream ifs(argv[2]);
#else
    std::ifstream ifs(argv[2], std::ios::binary);
#endif
    Momentum_recon::Mom_chain mom_chain;

    B2Writer writer(argv[3], reader);

    int num_entry = 0;
    int num_base, num_link;
#ifdef TEXT_MODE
    while ( ifs >> mom_chain.groupid >> mom_chain.chainid
	    >> mom_chain.unixtime >> mom_chain.entry_in_daily_file
	    >> mom_chain.mom_recon
	    >> num_base >> num_link ) {
      mom_chain.base.resize(num_base);
      mom_chain.base_pair.resize(num_link);
#else
    while ( ifs.read((char*)& mom_chain, sizeof(Momentum_recon::Mom_chain)) ) {
      num_base = mom_chain.base.size();
      num_link = mom_chain.base_pair.size();
#endif
      BOOST_LOG_TRIVIAL(debug) << "Entry : " << num_entry;
      auto &spill_summary = writer.GetSpillSummary();
      
      std::vector<int> rawid_vec(num_base);
      std::vector<int> plate_vec(num_base);
      std::vector<TVector3> film_position_vec(num_base);
      std::vector<TVector3> tangent_vec(num_base);
      std::vector<Double_t> vph_up_vec(num_base), vph_down_vec(num_base);
      std::vector<Double_t> pixel_count_up_vec(num_base), pixel_count_down_vec(num_base);
      
      std::vector<TVector3> absolute_position_vec(num_base);
      std::vector<TVector3> film_position_in_down_coordinate_vec(num_base);
      std::vector<TVector3> tangent_in_down_coordinate_vec(num_base);
      
      int plate, base_id;
      double tangent_x, tangent_y;
      double position_x, position_y, position_z;
      int tmp;
      int ph[2], pixel_area[2], pixel_hit[2];
      
      TVector3 absolute_position, film_position, tangent;
      TVector3 film_position_in_down_coordinate, tangent_in_down_coordinate;
      Double_t vph[2];
      
      BOOST_LOG_TRIVIAL(debug) << "Basetrack information get";
#ifdef TEXT_MODE
      for ( int ibase = 0; ibase < num_base; ibase++ ) {
	ifs >> mom_chain.base.at(ibase).pl
	    >> mom_chain.base.at(ibase).rawid
	    >> mom_chain.base.at(ibase).ax
	    >> mom_chain.base.at(ibase).ay
	    >> mom_chain.base.at(ibase).x
	    >> mom_chain.base.at(ibase).y 
	    >> mom_chain.base.at(ibase).z
	    >> mom_chain.base.at(ibase).m[0].zone
	    >> mom_chain.base.at(ibase).m[0].view
	    >> mom_chain.base.at(ibase).m[0].imager
	    >> mom_chain.base.at(ibase).m[0].ph
	    >> mom_chain.base.at(ibase).m[0].pixelnum
	    >> mom_chain.base.at(ibase).m[0].hitnum
	    >> mom_chain.base.at(ibase).m[1].zone
	    >> mom_chain.base.at(ibase).m[1].view
	    >> mom_chain.base.at(ibase).m[1].imager
	    >> mom_chain.base.at(ibase).m[1].ph
	    >> mom_chain.base.at(ibase).m[1].pixelnum
	    >> mom_chain.base.at(ibase).m[1].hitnum;
#else
      for ( int ibase = 0; ibase < mom_chain.base.size(); ibase++ ) {
#endif

       	absolute_position.SetX(mom_chain.base.at(ibase).x / 1e3);
	absolute_position.SetY(mom_chain.base.at(ibase).y / 1e3);
	absolute_position.SetZ(mom_chain.base.at(ibase).z / 1e3);
	// offset

	tangent.SetX(mom_chain.base.at(ibase).ax);
	tangent.SetY(mom_chain.base.at(ibase).ay);
	tangent.SetZ(1.);

	vph[0] = mom_chain.base.at(ibase).m[0].ph % 10000; 
	vph[1] = mom_chain.base.at(ibase).m[1].ph % 10000;
	
	rawid_vec.at(ibase) = mom_chain.base.at(ibase).rawid;
	plate_vec.at(ibase) = mom_chain.base.at(ibase).pl - 1;
	absolute_position_vec.at(ibase) = absolute_position;
	tangent_vec.at(ibase) = tangent;
	vph_up_vec.at(ibase) = vph[0];
	vph_down_vec.at(ibase) = vph[1];
	pixel_count_up_vec.at(ibase) = mom_chain.base.at(ibase).m[0].pixelnum;
	pixel_count_down_vec.at(ibase) = mom_chain.base.at(ibase).m[1].pixelnum;
	
      }
      
      BOOST_LOG_TRIVIAL(debug) << "Linklet information get";
      int link_plate[2], link_base_id[2];
      double link_tangent_x[2], link_tangent_y[2];
      double link_position_x[2], link_position_y[2], link_position_z[2];
#ifdef TEXT_MODE
      for ( int ilink = 0; ilink < num_link; ilink++ ) {
	ifs >> mom_chain.base_pair.at(ilink).first.pl
	    >> mom_chain.base_pair.at(ilink).first.rawid
	    >> mom_chain.base_pair.at(ilink).second.pl
	    >> mom_chain.base_pair.at(ilink).second.rawid
	    >> mom_chain.base_pair.at(ilink).first.ax
	    >> mom_chain.base_pair.at(ilink).first.ay
	    >> mom_chain.base_pair.at(ilink).first.x
	    >> mom_chain.base_pair.at(ilink).first.y
	    >> mom_chain.base_pair.at(ilink).first.z
	    >> mom_chain.base_pair.at(ilink).second.ax
	    >> mom_chain.base_pair.at(ilink).second.ay
	    >> mom_chain.base_pair.at(ilink).second.x
	    >> mom_chain.base_pair.at(ilink).second.y
	    >> mom_chain.base_pair.at(ilink).second.z;
#else
      for ( int ilink = 0; ilink < mom_chain.base_pair.size(); ilink++ ) {
#endif
	film_position.SetX(mom_chain.base_pair.at(ilink).first.x / 1e3);
	film_position.SetY(mom_chain.base_pair.at(ilink).first.y / 1e3);
	film_position.SetZ(mom_chain.base_pair.at(ilink).first.z / 1e3);
	film_position_in_down_coordinate.SetX(mom_chain.base_pair.at(ilink).second.x / 1e3);
	film_position_in_down_coordinate.SetY(mom_chain.base_pair.at(ilink).second.y / 1e3);
	film_position_in_down_coordinate.SetZ(mom_chain.base_pair.at(ilink).second.z / 1e3);
	tangent.SetX(mom_chain.base_pair.at(ilink).first.ax);
	tangent.SetY(mom_chain.base_pair.at(ilink).first.ay);
	tangent.SetZ(1.);
	tangent_in_down_coordinate.SetX(mom_chain.base_pair.at(ilink).second.ax);
	tangent_in_down_coordinate.SetY(mom_chain.base_pair.at(ilink).second.ay);
	tangent_in_down_coordinate.SetZ(1.);
	auto find_itr = std::find(rawid_vec.begin(), rawid_vec.end(),
				  mom_chain.base_pair.at(ilink).second.rawid);
	int index = find_itr - rawid_vec.begin();
	film_position_in_down_coordinate_vec.at(index) = film_position_in_down_coordinate;
	tangent_vec.at(index) = tangent;
	tangent_in_down_coordinate_vec.at(index) = tangent_in_down_coordinate;
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
      
      for ( int ibase = 0; ibase < mom_chain.base.size(); ibase++ ) {
	auto &emulsion_summary = spill_summary.AddEmulsion();
	emulsion_summary.SetEmulsionTrackId((UInt_t)rawid_vec.at(ibase));
	emulsion_summary.SetParentTrackId(mom_chain.chainid);
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
	emulsion_summary.SetEcc(4);
	emulsion_summary.SetPlate(plate_vec.at(ibase));
      }

      writer.Fill();
      num_entry++;
    }

    BOOST_LOG_TRIVIAL(debug) << "Total number of entries : " << num_entry;

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
