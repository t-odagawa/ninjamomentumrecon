// boost includes
#include <boost/log/core.hpp>
#include <boost/log/trivial.hpp>
#include <boost/log/expressions.hpp>

// B2 includes
#include <B2Reader.hh>
#include <B2SpillSummary.hh>
#include <B2EventSummary.hh>
#include <B2VertexSummary.hh>
#include <B2TrackSummary.hh>
#include <B2EmulsionSummary.hh>
#include <B2Enum.hh>
#include <B2Pdg.hh>

// ROOT includes
#include <TStyle.h>
#include <TFile.h>
#include <TCanvas.h>
#include <TH2D.h>
#include <TF1.h>

#include "McsConst.hpp"
#include "McsFunction.hpp"
#include "RangeFunction.hpp"
#include "RangeSpline.hpp"

namespace logging = boost::log;

int main (int argc, char *argv[]) {

  logging::core::get()->set_filter
    (
     logging::trivial::severity >= logging::trivial::info
     // logging::trivial::severity >= logging::trivial::debug
     );

  BOOST_LOG_TRIVIAL(info) << "==========Range Plate Function Start==========";

  if (argc != 5) {
    BOOST_LOG_TRIVIAL(error) << "Usage : " << argv[0]
			     << " <input B2 file name> <output root file name> <data type MC(0)/Data(1)> <spline dir path>";
    std::exit(1);
  }

  try {

    gStyle->SetOptStat(0);
    gErrorIgnoreLevel = kWarning;

    B2Reader reader(argv[1]);

    TString ofilename = argv[2];
    TFile *ofile = new TFile(ofilename, "recreate");
    TTree *otree = new TTree("tree", "tree");
    Int_t entry_in_daily_file;
    Double_t weight;
    Int_t true_particle_id;
    Int_t npl;
    Int_t direction;
    Int_t stop_flag;
    Double_t true_momentum;
    Double_t recon_momentum;
    std::vector<Int_t > pl;
    std::vector<Double_t > ax, ay;
    std::vector<Double_t > energy_deposit, vph, pixel_count;
    otree->Branch("entry_in_daily_file", &entry_in_daily_file, "entry_in_daily_file/I");
    otree->Branch("weight", &weight, "weight/D");
    otree->Branch("true_particle_id", &true_particle_id, "true_particle_id/I");
    otree->Branch("npl", &npl, "npl/I");
    otree->Branch("direction", &direction, "direction/I");
    otree->Branch("stop_flag", &stop_flag, "stop_flag/I");
    otree->Branch("true_momentum", &true_momentum, "true_momentum/D");
    otree->Branch("recon_momentum", &recon_momentum, "recon_momentum/D");
    otree->Branch("pl", &pl);
    otree->Branch("ax", &ax);
    otree->Branch("ay", &ay);
    otree->Branch("energy_deposit", &energy_deposit);
    otree->Branch("vph", &vph);
    otree->Branch("pixel_count", &pixel_count);

    Int_t datatype = std::atoi(argv[3]);
    Int_t ecc_id = 4;

    const std::string spline_dir_path = argv[4];
    const RangeSpline range_spline(spline_dir_path);
    const RangeFunction range_function(range_spline);

    int num_entry = 0;

    while (reader.ReadNextSpill() > 0) {
      auto &spill_summary = reader.GetSpillSummary();

      // muon id      

      std::vector<const B2EmulsionSummary*> emulsions;
      auto it_emulsion = spill_summary.BeginEmulsion();
      while ( const auto *emulsion = it_emulsion.Next() ) {
	if ( datatype == B2DataType::kMonteCarlo &&
	     emulsion->GetParentTrackId() == 0 ) continue; // track id should be assigned when MC
	if ( emulsion->GetFilmType() != B2EmulsionType::kECC ) continue;
	if ( emulsion->GetEcc() != ecc_id ) continue;
	emulsions.push_back(emulsion);
      }

      if ( emulsions.empty() ) continue;

      entry_in_daily_file = reader.GetEntryNumber();

      std::sort(emulsions.begin(), emulsions.end(), EmulsionCompare);

      auto it_event = spill_summary.BeginTrueEvent();
      const auto *event = it_event.Next();
      double norm = event->GetNormalization();
      auto &primary_vertex_summary = event->GetPrimaryVertex();
      double total_cross_section = primary_vertex_summary.GetTotalCrossSection();
      weight = norm * total_cross_section * 1.e-38 * 6.02e23;

      // Separate emulsions into chains
      std::vector<std::vector<const B2EmulsionSummary* > > emulsion_single_chains;
      Int_t track_id_tmp_ = emulsions.at(0)->GetParentTrackId();
      std::vector<const B2EmulsionSummary* > emulsion_single_chain;
      for ( Int_t iemulsion = 0; iemulsion < emulsions.size(); iemulsion++ ) {
	if ( emulsions.at(iemulsion)->GetParentTrackId() == track_id_tmp_ ) {
	  emulsion_single_chain.push_back(emulsions.at(iemulsion));
	} else {
	  emulsion_single_chains.push_back(emulsion_single_chain);
	  emulsion_single_chain.clear();
	  track_id_tmp_ = emulsions.at(iemulsion)->GetParentTrackId();
	  emulsion_single_chain.push_back(emulsions.at(iemulsion));
	}
      }
      emulsion_single_chains.push_back(emulsion_single_chain);
      
      // loop for each chain
      for ( auto chain : emulsion_single_chains ) {

	npl = chain.size();

	// Get true information
	Int_t particle_pdg = chain.at(0)->GetParentTrack().GetParticlePdg();
	if ( datatype == B2DataType::kMonteCarlo ) {
	  true_particle_id = particle_pdg;
	  Double_t direction_tmp_ = chain.at(0)->GetTangent().GetValue().Z();
	  if ( direction_tmp_ > 0 ) direction = 1;
	  else direction = -1;
	  Int_t momentum_ = chain.at(0)->GetParentTrack().GetInitialAbsoluteMomentum().GetValue();
	  if ( B2Pdg::IsChargedPion(particle_pdg) )
	    true_momentum = momentum_;
	  else if ( particle_pdg == PDG_t::kProton )
	    true_momentum = momentum_;
	  else
	    BOOST_LOG_TRIVIAL(trace) << "Particle is not in interest";
	}

	if ( particle_pdg != PDG_t::kProton &&
	     !B2Pdg::IsChargedPion(particle_pdg) ) continue;

	stop_flag = (Int_t)IsStopInEccFiducial(chain, direction);

	for ( const auto emulsion : chain ) {
	  TVector3 tangent = emulsion->GetTangent().GetValue();
	  tangent = (1. / tangent.Z()) * tangent;
	  SmearTangentVector(tangent);
	  if (std::fabs(tangent.X()) > 4 || std::fabs(tangent.Y()) > 4 ) continue;
	  pl.push_back(emulsion->GetPlate() + 1);
	  ax.push_back(tangent.X());
	  ay.push_back(tangent.Y());
	  energy_deposit.push_back(emulsion->GetEdepSum());
	  vph.push_back(emulsion->GetVphUp() + emulsion->GetVphDown());
	  pixel_count.push_back(emulsion->GetPixelCountUp() + emulsion->GetPixelCountDown());
	}

	if ( pl.empty() ) continue;
	
	range_function.ModifyVectors(ax, ay, pl);
	recon_momentum = range_function.CalculateEnergyFromRange(ax, ay, pl, particle_pdg, direction);
	switch ( particle_pdg ) {
	case PDG_t::kProton :
	  recon_momentum += MCS_PROTON_MASS;
	  recon_momentum = CalculateMomentumFromEnergy(recon_momentum, MCS_PROTON_MASS);
	  break;
	case PDG_t::kPiPlus :
	case PDG_t::kPiMinus :
	  recon_momentum += MCS_PION_MASS;
	  recon_momentum = CalculateMomentumFromEnergy(recon_momentum, MCS_PION_MASS);
	  break;
	}

	otree->Fill();
	pl.clear();
	ax.clear(); ay.clear();
	vph.clear(); pixel_count.clear();
      }

      // if ( num_entry > 10) break;
      num_entry++;

    }

    ofile->cd();
    otree->Write();
    ofile->Close();

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
