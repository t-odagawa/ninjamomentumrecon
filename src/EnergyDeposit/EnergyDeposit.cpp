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

  logging::core::get()->set_filter
    (
     logging::trivial::severity >= logging::trivial::info
     // logging::trivial::severity >= logging::trivial::debug
     );

  BOOST_LOG_TRIVIAL(info) << "==========Energy Deposit Study Start==========";

  if (argc != 1) {
    BOOST_LOG_TRIVIAL(error) << "Usage : " << argv[0];
    std::exit(1);
  }

  try {

    gStyle->SetOptStat(0);
    gErrorIgnoreLevel = kWarning;

    TString ofilename = "momentum_drop.root";
    TFile *ofile = new TFile(ofilename, "recreate");
    TTree *otree = new TTree("tree", "tree");
    Double_t momentum_drop, momentum_drop_err;
    Double_t momentum, sideview, topview;
    otree->Branch("momentum_drop", &momentum_drop, "momentum_drop/D");
    otree->Branch("momentum_drop_err", &momentum_drop_err, "momentum_drop_err/D");
    otree->Branch("momentum", &momentum, "momentum/D");
    otree->Branch("sideview", &sideview, "sideview/D");
    otree->Branch("topview", &topview, "topview/D");

    TF1 *linear = new TF1("linear","[0] + [1] * x");

    TCanvas *c = new TCanvas("c", "c");
    TString pdfname = "momentum_drop.pdf";
    c->Print(pdfname + "[", "pdf");

    for (Int_t imomentum = 0; imomentum < 13; imomentum++) {
      momentum = 300 + imomentum * 100;
      for (Int_t isideview = 0; isideview < 6; isideview++){
	sideview = isideview * 15;
	for (Int_t itopview = 0; itopview < 6; itopview++) {
	  topview = itopview * 15;

	  BOOST_LOG_TRIVIAL(debug) << "Momentum = " << (Int_t)momentum << " [MeV/c], "
				   << "Side angle = " << (Int_t)sideview << " [deg], "
				   << "Top angle = " << (Int_t)topview << " [deg]";

	  TString filename = Form("/home/t2k/odagawa/data/mc_data/particlegun/particlegun_muon_%dMeV_m350_0_m230_%d_%d.root",
				  (Int_t)momentum, (Int_t)sideview, (Int_t)topview);
	  B2Reader reader(filename);

	  TH2D *hist_mom_pl = new TH2D("hist_mom_pl","Momentum drop inside the ECC;pl;Momentum [MeV/c]", 133, 0, 133, 100, momentum - 200., momentum);

	  while (reader.ReadNextSpill() > 0) {
	    auto &spill_summary = reader.GetSpillSummary();
	    
	    // Get emulsion tracks
	    std::vector<const B2EmulsionSummary*> emulsions;
	    auto it_emulsion = spill_summary.BeginEmulsion();
	    while (const auto *emulsion = it_emulsion.Next()) {
	      if (emulsion->GetParentTrackId() == 0) continue;
	      if (emulsion->GetParentTrack().GetParticlePdg() != 13) continue;
	      if (emulsion->GetFilmType() != B2EmulsionType::kECC) continue;
	      emulsions.push_back(emulsion);
	    }
	    
	    if (emulsions.size() <= 0) continue;
	    
	    std::sort(emulsions.begin(), emulsions.end(), emulsion_compare);
	    
	    int track_id_tmp_ = emulsions.at(0)->GetParentTrackId();
	    int ecc_tmp_ = emulsions.at(0)->GetEcc();
	    
	    for (const auto &emulsion : emulsions) {
	      
	      int track_id_ = emulsion->GetParentTrackId();
	      int ecc_ = emulsion->GetEcc();
	      
	      // Only consider muon? track in the same ECC
	      if (track_id_ != track_id_tmp_ || ecc_ != ecc_tmp_) continue;
	      hist_mom_pl->Fill(emulsion->GetPlate(), emulsion->GetMomentum().GetValue().Mag());
	    }
            
	  }

	  c->cd();
	  hist_mom_pl->SetTitle(Form("Momentum drop inside the ECC (p_{0} = %d MeV/c, #theta_{x} = %d deg., #theta_{y} = %d deg.);pl;Momentum [MeV/c]",
				     (Int_t)momentum, (Int_t)sideview, (Int_t)topview));
	  hist_mom_pl->Draw("colz");
	  c->Print(pdfname);

	  // Remove small entry bins
	  for (int ihistx = 1; ihistx <= hist_mom_pl->GetNbinsX(); ihistx++) {
	    for (int ihisty = 1; ihisty <= hist_mom_pl->GetNbinsY(); ihisty++) {
	      if (hist_mom_pl->GetBinContent(ihistx, ihisty)
		  < 0.1 * hist_mom_pl->GetBinContent(hist_mom_pl->GetMaximumBin()))
		hist_mom_pl->SetBinContent(ihistx, ihisty, 0.);
	    }
	  }

	  c->cd();
	  hist_mom_pl->Draw("colz");
	  linear->SetParameter(1, 0);
	  linear->SetParError(1, 0);
	  hist_mom_pl->Fit(linear,"Q");
	  c->Print(pdfname);

	  BOOST_LOG_TRIVIAL(info) << "Momentum drop = " << linear->GetParameter(1) << " +/- "
				  << linear->GetParError(1) << " [MeV/c/plate]";
	  momentum_drop = linear->GetParameter(1);
	  momentum_drop_err = linear->GetParError(1);

	  otree->Fill();

	}
      }
    }
    
    c->Print(pdfname + "]", "pdf");
    ofile->cd();
    otree->Write();
    ofile->Close();
    
  } catch (const std::runtime_error &error) {
    BOOST_LOG_TRIVIAL(fatal) << "Runtime error : " << error.what();
    std::exit(1);
  } catch (const std::invalid_argument &error) {
    BOOST_LOG_TRIVIAL(fatal) << "Invalid argument error : " << error.what();
    std::exit(1);
  }

  BOOST_LOG_TRIVIAL(info) << "==========Energy Deposit Study Finish==========";
  std::exit(0);

}
