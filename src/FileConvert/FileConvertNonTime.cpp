// boost includes
#include <boost/log/core.hpp>
#include <boost/log/trivial.hpp>
#include <boost/log/expressions.hpp>

// B2 includes
#include <B2Writer.hh>
#include <B2SpillSummary.hh>
#include <B2EventSummary.hh>
#include <B2VertexSummary.hh>
#include <B2EmulsionSummary.hh>

// system includes
#include <fstream>
#include <iostream>
#include <iomanip>
#include <vector>
#include <algorithm>
#include <chrono>

// my include
#include "McsClass.hpp"
#include "McsFunction.hpp"

namespace logging = boost::log;

int main (int argc, char *argv[]) {

  logging::core::get()->set_filter
    (
     logging::trivial::severity >= logging::trivial::info
     //logging::trivial::severity >= logging::trivial::debug
     );

  BOOST_LOG_TRIVIAL(info) << "==========MCS Convert w/o Muon ID Start==========";

  if ( argc != 4 ) {
    BOOST_LOG_TRIVIAL(error) << "Usage : " << argv[0]
			     << " <Input momch file> <Output B2 file> <ECC id (0-8)>";
    std::exit(1);
  }

  try {

    std::vector<Momentum_recon::Event_information > ev_vec = Momentum_recon::ReadEventInformationBin((std::string)argv[1]);

    B2Writer writer(argv[2]);
    int ecc_id = std::atoi(argv[3]);

    int num_entry = 0;

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

    for ( auto ev : ev_vec ) {
      auto &spill_summary = writer.GetSpillSummary();
      auto &event_summary = spill_summary.AddTrueEvent();
      event_summary.SetEventType(B2EventType::kEventIdStorage);

      auto &recon_vertex_summary = spill_summary.AddReconVertex();
      TVector3 vertex_position;
      vertex_position.SetX(ev.recon_vertex_position[0]);
      vertex_position.SetY(ev.recon_vertex_position[1]);
      vertex_position.SetZ(ev.recon_vertex_position[2]);
      recon_vertex_summary.SetRelativePosition(vertex_position);
      recon_vertex_summary.SetDetector(B2Detector::kNinja);
      recon_vertex_summary.SetIsInsideFiducialVolume(kTRUE);
      recon_vertex_summary.SetInteractionMaterial((B2Material)ev.vertex_material);
      recon_vertex_summary.SetPlane(ev.ecc_id * 1000 + ev.vertex_pl - 1);
      event_summary.SetPrimaryVertex(recon_vertex_summary);

      for ( auto chain : ev.chains ) {

	rawid_vec.clear();
	plate_vec.clear();
	film_position_vec.clear();
	tangent_vec.clear();
	vph_up_vec.clear();
	vph_down_vec.clear();
	pixel_count_up_vec.clear();
	pixel_count_down_vec.clear();
	absolute_position_vec.clear();
	film_position_in_down_coordinate_vec.clear();
	tangent_in_down_coordinate_vec.clear();

	
	
      }

    }


    Momentum_recon::Mom_chain mom_chain;
    Momentum_recon::Mom_basetrack mom_basetrack;
    std::pair<Momentum_recon::Mom_basetrack, Momentum_recon::Mom_basetrack > mom_basetrack_pair;

    int ecc_id = std::atoi(argv[3]);

    int num_entry = 0;
    int num_base, num_link;


    while ( Momentum_recon::ReadMomChainHeader(ifs, mom_chain, num_base, num_link) ) {

      // auto start = std::chrono::system_clock::now();

      if (num_entry % 100 == 0) {
	nowpos = ifs.tellg();
	auto size1 = nowpos - begpos;
	std::cerr << std::right << std::fixed << "\r now reading ..." << std::setw(4) << std::setprecision(1) << size1 * 100. / size2 << "%";
      }

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

      // auto time1 = std::chrono::system_clock::now();

      BOOST_LOG_TRIVIAL(debug) << "Entry : " << num_entry;

      auto &spill_summary = writer.GetSpillSummary();
      auto &event_summary = spill_summary.AddTrueEvent();
      event_summary.SetEventType(B2EventType::kEventIdStorage);

      rawid_vec.clear(); rawid_vec.resize(num_base);
      plate_vec.clear(); plate_vec.resize(num_base);
      film_position_vec.clear(); film_position_vec.resize(num_base);
      tangent_vec.clear(); tangent_vec.resize(num_base);
      vph_up_vec.clear(); vph_down_vec.clear();
      vph_up_vec.resize(num_base); vph_down_vec.resize(num_base);
      pixel_count_up_vec.clear(); pixel_count_down_vec.clear();
      pixel_count_up_vec.resize(num_base); pixel_count_down_vec.resize(num_base);
    
      absolute_position_vec.clear();
      absolute_position_vec.resize(num_base);
      film_position_in_down_coordinate_vec.clear();
      film_position_in_down_coordinate_vec.resize(num_base);
      tangent_in_down_coordinate_vec.clear();
      tangent_in_down_coordinate_vec.resize(num_base);

      // auto time2 = std::chrono::system_clock::now();

      TVector3 absolute_position, film_position, tangent;
      TVector3 film_position_in_down_coordinate, tangent_in_down_coordinate;
      Double_t vph[2];
      
      for ( int ibase = 0; ibase < mom_chain.base.size(); ibase++ ) {
	absolute_position.SetX(mom_chain.base.at(ibase).x / 1e3);
	absolute_position.SetY(mom_chain.base.at(ibase).y / 1e3);
	absolute_position.SetZ(mom_chain.base.at(ibase).z / 1e3);
	// offset
	PositionAddOffset(absolute_position, ecc_id);

	tangent.SetX(mom_chain.base.at(ibase).ax);
	tangent.SetY(mom_chain.base.at(ibase).ay);
	tangent.SetZ(1.);
      
	vph[0] = (Double_t)(mom_chain.base.at(ibase).m[0].ph % 10000); 
	// vph[0] = (Double_t)(mom_chain.base.at(ibase).m[0].ph); 
	vph[1] = (Double_t)(mom_chain.base.at(ibase).m[1].ph % 10000);
	// vph[1] = (Double_t)(mom_chain.base.at(ibase).m[1].ph);
	rawid_vec.at(ibase) = mom_chain.base.at(ibase).rawid;
	plate_vec.at(ibase) = mom_chain.base.at(ibase).pl - 1;
	absolute_position_vec.at(ibase) = absolute_position;
	tangent_vec.at(ibase) = tangent;
	vph_up_vec.at(ibase) = vph[0];
	vph_down_vec.at(ibase) = vph[1];
	pixel_count_up_vec.at(ibase) = mom_chain.base.at(ibase).m[0].hitnum;
	pixel_count_down_vec.at(ibase) = mom_chain.base.at(ibase).m[1].hitnum;
      }
      
      // Linklet information write
      for ( int ilink = 0; ilink < mom_chain.base_pair.size(); ilink++ ) {
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
	if ( find_itr != rawid_vec.end() ) {
	  int index = std::distance(rawid_vec.begin(), find_itr);
	  film_position_vec.at(index - 1) = film_position;
	  film_position_in_down_coordinate_vec.at(index) = film_position_in_down_coordinate;
	  tangent_vec.at(index - 1) = tangent;
	  tangent_in_down_coordinate_vec.at(index) = tangent_in_down_coordinate;
	}
      }
      
      // auto time3 = std::chrono::system_clock::now();

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
    
      // auto time4 = std::chrono::system_clock::now();

      for ( int ibase = 0; ibase < mom_chain.base.size(); ibase++ ) {
	auto &emulsion_summary = spill_summary.AddEmulsion();
	emulsion_summary.SetEmulsionTrackId((UInt_t)rawid_vec.at(ibase));
	// emulsion_summary.SetParentTrackId(mom_chain.groupid * 100000 + mom_chain.chainid);
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

      }

      writer.Fill(); 

      // auto time5 = std::chrono::system_clock::now();     
      num_entry++;
      // if (num_entry > 10) break;

      /*
      std::cout << "Binary read" << std::endl;
      std::cout << std::chrono::duration_cast<std::chrono::nanoseconds>(time1 - start).count() << std::endl;
      std::cout << "Vector setup" << std::endl;
      std::cout << std::chrono::duration_cast<std::chrono::nanoseconds>(time2 - time1).count() << std::endl;
      std::cout << "Binary to vector" << std::endl;
      std::cout << std::chrono::duration_cast<std::chrono::nanoseconds>(time3 - time2).count() << std::endl;
      std::cout << "Vector reverse" << std::endl;
      std::cout << std::chrono::duration_cast<std::chrono::nanoseconds>(time4 - time3).count() << std::endl;
      std::cout << "B2 file write" << std::endl;
      std::cout << std::chrono::duration_cast<std::chrono::nanoseconds>(time5 - time4).count() << std::endl;
      */
    }

    auto size1 = eofpos - begpos;
    std::cerr << "\r now reading ..." << std::setw(4) << std::setprecision(1) << size1 * 100. / size2 << "%" << std::endl;

    BOOST_LOG_TRIVIAL(debug) << "Total number of entries : " << num_entry;

  } catch (const std::runtime_error &error) {
    BOOST_LOG_TRIVIAL(fatal) << "Runtime error : " << error.what();
    std::exit(1);
  } catch (const std::invalid_argument & error) {
    BOOST_LOG_TRIVIAL(fatal) << "Invalid argument error : " << error.what();
    std::exit(1);
  }

  BOOST_LOG_TRIVIAL(info) << "==========MCS Convert w/o Muon ID Finish==========";
  std::exit(0);

}
