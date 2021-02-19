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

  /*  logging::core::get()->set_filter
    (
     logging::trivial::severity >= logging::trivial::debug
     );
  */
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
    tree->Branch("true_pbeta", &true_pbeta);
    tree->Branch("recon_pbeta", &recon_pbeta);

    int ievent = 0;
    int particle_id = PDG_t::kMuonMinus;
    while(reader.ReadNextSpill() > 0) {
      auto &spill_summary = reader.GetSpillSummary();
      BOOST_LOG_TRIVIAL(debug) << "Entry : " << ievent;
      //true_pbeta = get_true_pbeta(spill_summary, particle_id);
#ifdef COORDINATE_METHOD
      recon_pbeta = mcs_coordinate_method(spill_summary, particle_id);
#else
      recon_pbeta = mcs_angle_method(spill_summary, particle_id);
#endif

      ievent++;
      if (recon_pbeta.size() > 0)
	tree->Fill();
      //if (ievent > 10) break;
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
