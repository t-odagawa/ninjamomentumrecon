// boost includes
#include <boost/log/core.hpp>
#include <boost/log/trivial.hpp>
#include <boost/log/expressions.hpp>

// B2 includes
#include <B2Reader.hh>
#include <B2SpillSummary.hh>
#include <B2TrackSummary.hh>
#include <B2EmulsionSummary.hh>

// ROOT includes
#include <TStyle.h>
#include <TFile.h>
#include <TCanvas.h>
#include <TH2D.h>
#include <TF1.h>

// my include
#include "/home/t2k/odagawa/NinjaMomentumRecon/src/McsCommon.cpp"

namespace logging = boost::log;

int main (int argc, char *argv[]) {

  logging::core:get()->set_filter
    (
     logging::trivial::severity >= logging::trivial::info
     // logging::trivial::severity >= logging::trivial::debug
     );

  BOOST_LOG_TRIVIAL(info) << "==========Range Plate Function Start==========";

  if (argc != 2) {
    BOOST_LOG_TRIVIAL(error) << "Usage : " << argv[0];
    std::exit(1);
  }

  try {

    gStyle->SetOptStat(0);
    gErrorIgnoreLevel = kWarning;

    B2Reader reader(argv[1]);

    while (reader.ReadNextSpill() > 0) {
      auto &spill_summary = reader.GetSpillSummary();

      std::vector<const B2EmulsionSummary*> emulsions;
      auto it_emulsion = spill_summary.BeginEmulsion();
      while (const auto *emulsion = it_emulsion.Next()) {
	if (emulsion->GetParentTrackid() == 0) continue;
	if (emulsion->GetParentTrack().GetParticlePdg() != 2212) continue;
	if (emulsion->GetFilmType() != B2EmulsionType::kECC) continue;
	emulsions.push_back(emulsion);
      }

      if (emulsions.size() <= 0) continue;

      std::sort(emulsions.begin(), emulsions.end(), emulsion_compare);

      int number_of_tracks = 0;
      std::vector<const B2EmulsionSummary*> emulsions_one_track;

      int track_id_tmp_ = emulsions.at(0)->GetParentTrackId();
      int ecc_tmp_ = emulsions.at(0)->GetEcc();

      for (const auto &emulsion : emulsions) {

	int track_id_ = emulsion->GetParentTrackId();
	int ecc_ = emulsion->GetEcc();

	if (track_id_ != track_id_tmp_ || ecc != ecc_tmp_ || emulsion == emulsions.back()) {
	  // All base tracks consisting of the same particle's track
	  if (emulsion == emulsions.back()) emulsions_one_track.push_back(emulsion);
	  number_of_tracks++;

	  int vertex_plate = emulsions_one_track.at(0)->GetPlate();


	}

      }


    }

  } catch (const std::runtime_error &error) {
    BOOST_LOG_TRIVIAL(fatal) << "Runtime error : " << error.what();
    std::exit(1);
  } catch (const std::invalid_argument &error) {
    BOOST_LOG_TRIVIAL(fatal) << "Invalid argument : " << error.what();
    std::exit(1);
  }

  BOOST_LOG_TRIVIAL(info) << "==========Range Plate Function Finish==========";
  std::exit(0);

}
