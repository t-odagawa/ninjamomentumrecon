// boost includes
#include <boost/log/core.hpp>
#include <boost/log/trivial.hpp>
#include <boost/log/expressions.hpp>

#include <boost/filesystem.hpp>
#include <boost/regex.hpp>

// B2 includes
#include <B2Reader.hh>
#include <B2SpillSummary.hh>
#include <B2EventSummary.hh>
#include <B2TrackSummary.hh>
#include <B2EmulsionSummary.hh>

// ROOT includes
#include <TStyle.h>
#include <TFile.h>
#include <TTree.h>

// system includes
#include <iostream>
#include <iomanip>
#include <fstream>
#include <vector>
#include <algorithm>

// my include
#include "/home/t2k/odagawa/NinjaMomentumRecon/src/McsCommon.cpp"


namespace logging = boost::log;
namespace fs = boost::filesystem;

class SimpleMCBaseTrackCollection{
public:
  int particle_id, plate_id, ecc_id;
  double energy_deposit_1, energy_deposit_2, ax, ay, x, y, z;
  int track_id, event_id, file_id;
  double weight;
};

int main (int argc, char *argv[]) {

  logging::core::get()->set_filter
    (
     //logging::trivial::severity >= logging::trivial::info
     logging::trivial::severity >= logging::trivial::debug
     );

  if (argc != 2) {
    BOOST_LOG_TRIVIAL(error) << "Usage : " << argv[0]
			     << " <output file name>";
    std::exit(1);
  }

  try {
  
    std::ofstream ofs(argv[1], std::ios::binary);
    SimpleMCBaseTrackCollection smbtc;
    
    const std::string target_directory_path("/home/t2k/odagawa/data/mc_data/pra1/");
    const std::string root_extension(".root");
    
    std::vector<std::string> all_matched_files;
    
    boost::filesystem::directory_iterator end_itr;
    for( boost::filesystem::directory_iterator i( target_directory_path ); i != end_itr; ++i )
      {
	
	// Skip if not a file
	if ( !boost::filesystem::is_regular_file( i->status() ) ) continue;
	
	if ( i->path().extension().string() != root_extension) continue;
	
	BOOST_LOG_TRIVIAL(debug) << i->path().filename().string();
	
	// File matches, store it
	all_matched_files.push_back( i->path().filename().string() );
      }
    
    std::sort(all_matched_files.begin(), all_matched_files.end());
    
    int fileid = -1;
    
    for (auto filename : all_matched_files) {
      
      fileid++;

      B2Reader reader(target_directory_path + filename);
      BOOST_LOG_TRIVIAL(debug) << "Input B2 filename : " << filename;
      
      int eventid = -1;
      
      while (reader.ReadNextSpill() > 0) {

	eventid++;

	auto &spill_summary = reader.GetSpillSummary();
	
	auto it_event = spill_summary.BeginTrueEvent();
	const auto *event = it_event.Next();
	double normalization = event->GetNormalization();
	double total_cross_section = event->GetTotalCrossSection() * std::pow(10, -38);
	
	
	std::vector<const B2EmulsionSummary*> emulsions;
	auto it_emulsion = spill_summary.BeginEmulsion();
	while (const auto *emulsion = it_emulsion.Next()) {
	  if (emulsion->GetParentTrackId() == 0) continue;
	  if (emulsion->GetParentTrack().GetParticlePdg() == 13) continue;
	  if (emulsion->GetFilmType() != B2EmulsionType::kECC) continue;
	  emulsions.push_back(emulsion);
	}
	
	if (emulsions.size() <= 0) continue;
	
	std::sort(emulsions.begin(), emulsions.end(), emulsion_compare);
	
	int number_of_tracks = 0;
	
	for (const auto emulsion : emulsions) {
	  smbtc.particle_id = emulsion->GetParentTrack().GetParticlePdg();
	  smbtc.plate_id = emulsion->GetPlate() + 1;
	  smbtc.ecc_id = emulsion->GetEcc() + 1;
	  smbtc.energy_deposit_1 = (emulsion->GetEdepSum() + emulsion->GetEdepDiff()) / 2.;
	  smbtc.energy_deposit_2 = (emulsion->GetEdepSum() - emulsion->GetEdepDiff()) / 2.;
	  smbtc.ax = emulsion->GetTangent().GetValue().X();
	  smbtc.ay = emulsion->GetTangent().GetValue().Y();
	  smbtc.x = emulsion->GetFilmPosition().GetValue().X() * 1000.;
	  smbtc.y = emulsion->GetFilmPosition().GetValue().Y() * 1000.;
	  smbtc.z = emulsion->GetAbsolutePosition().GetValue().Z() * 1000.;
	  smbtc.track_id = emulsion->GetParentTrackId();
	  smbtc.event_id = eventid;
	  smbtc.file_id = fileid;
	  smbtc.weight = normalization * total_cross_section * 6.02 * std::pow(10, 23) * 58.;
	  ofs.write((char*)& smbtc, sizeof(SimpleMCBaseTrackCollection));
	}       
		
      }
            
    }      

    ofs.close();

  } catch (const std::runtime_error &error){
    BOOST_LOG_TRIVIAL(fatal) << "Runtime error : " << error.what();
    std::exit(1);
  } catch (const std::invalid_argument &error) {
    BOOST_LOG_TRIVIAL(fatal) << "Invalid argument : " << error.what();
    std::exit(1);
  }
  
  std::exit(0);
  
}
