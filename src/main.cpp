// boost includes
#include <boost/log/core.hpp>
#include <boost/log/trivial.hpp>
#include <boost/log/expressions.hpp>

// B2 includes
#include <B2Reader.hh>
#include <B2SpillSummary.hh>

// ROOT includes
#include <TFile.h>
#include <TTree.h>
#include <TPDGCode.h>

#include "McsCommon.hpp"
#ifdef COORDINATE_METHOD
#include "McsCoordinateMethod.hpp"
#else
#include "McsAngleMethod.hpp"
#endif

namespace logging = boost::log;

int main (int argc, char *argv[]) {

    logging::core::get()->set_filter
    (
     //logging::trivial::severity >= logging::trivial::info
      logging::trivial::severity >= logging::trivial::debug
     );
    
  BOOST_LOG_TRIVIAL(info) << "==========NINJA Momentum Reconstruction Start==========";

  if (argc != 3) {
    BOOST_LOG_TRIVIAL(error) << "Usage : " << argv[0]
			     << " <input B2 file path> <output file path>";
    std::exit(1);
  }

#ifdef COORDINATE_METHOD
  BOOST_LOG_TRIVIAL(info) << "Coordinate method test";
#else
  BOOST_LOG_TRIVIAL(info) << "Angle method test";
#endif

  try {

    B2Reader reader(argv[1]);

    TFile *file = new TFile(argv[2], "recreate");
    TTree *tree = new TTree("tree", "ECC reconstructed momentum");
    std::vector<Double_t> true_pbeta = {};
    std::vector<Double_t> recon_pbeta = {};
    std::vector<Double_t> path_length = {};
    std::vector<Double_t> angle_difference = {};
    tree->Branch("true_pbeta", &true_pbeta);
    tree->Branch("recon_pbeta", &recon_pbeta);
    tree->Branch("path_length", &path_length);
    tree->Branch("angle_difference", &angle_difference);

    int ievent = 0;
    int particle_id = PDG_t::kMuonMinus;
    while(reader.ReadNextSpill() > 0) {
      auto &spill_summary = reader.GetSpillSummary();
      BOOST_LOG_TRIVIAL(debug) << "Entry : " << ievent;
      //true_pbeta = get_true_pbeta(spill_summary, particle_id);
#ifdef COORDINATE_METHOD
      recon_pbeta = mcs_coordinate_method(spill_summary, particle_id);
#else
      if(!mcs_angle_method(spill_summary, recon_pbeta, angle_difference, path_length, particle_id)) {
	BOOST_LOG_TRIVIAL(debug) << "pbeta is not reconstructed for entry : " << ievent;
	continue;
      }
#endif

      ievent++;

      tree->Fill();
      //if (ievent > 100) break;
      recon_pbeta.clear();
      angle_difference.clear();
      path_length.clear();
      recon_pbeta.shrink_to_fit();
      angle_difference.shrink_to_fit();
      path_length.shrink_to_fit();
    }

    file->cd();
    tree->Write();
    file->Close();

  } catch (const std::runtime_error &error) {
    BOOST_LOG_TRIVIAL(fatal) << "Runtime error : " << error.what();
    std::exit(1);
  } catch (const std::invalid_argument &error) {
    BOOST_LOG_TRIVIAL(fatal) << "Invalid argument error : " << error.what();
    std::exit(1);
  }

  BOOST_LOG_TRIVIAL(info) << "==========NINJA Momentum Reconstruction Finish==========";
  std::exit(0);

}
