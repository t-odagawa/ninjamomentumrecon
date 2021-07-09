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
#include <TCanvas.h>
#include <TFile.h>
#include <TH2D.h>
#include <TF1.h>

// my include 
#include "/home/t2k/odagawa/NinjaMomentumRecon/src/McsCommon.cpp"

namespace logging = boost::log;

int main (int argc, char *argv[]) {

  logging::core::get()->set_filter
    (
     logging::trivial::severity >= logging::trivial::info
     // loggin::trivial::severity >= logging::trivial::debug
     );
  
  BOOST_LOG_TRIVIAL(info) << "==========Bethe Bloch Fit Start==========";
  
  if (argc != 1) {
    BOOST_LOG_TRIVIAL(error) << "Usage : " << argv[0];
    std::exit(1);
  }
  
  try {
    
    gStyle->SetOptStat(0);
    gErrorIgnoreLevel = kWarning;
    const Double_t muon_mass = 105.658 * MeV;
    
    
    TCanvas *c = new TCanvas("c", "c");
    //TString pdfname = "bethe_bloch_fit.pdf";
    TString pdfname = "bethe_bloch_fit_new.pdf";
    c->Print(pdfname + "[", "pdf");

    Double_t momentum, sideview, topview;
    Double_t x, y, z;
    z = -183.08;

    TFile *ofile = new TFile("bethe_bloch_fit.root","recreate");
    TTree *otree = new TTree("tree", "tree");
    Int_t plate;
    Double_t momentum_first, energy_first, momentum_next, energy_next;
    Double_t beta;
    Double_t dz;
    Double_t enedep;

    otree->Branch("plate", &plate, "plate/I");
    otree->Branch("momentum_first", &momentum_first, "momentum_first/D");
    otree->Branch("energy_first", &energy_first, "energy_first/D");
    otree->Branch("momentum_next", &momentum_next, "momentum_next/D");
    otree->Branch("energy_next", &energy_next, "energy_next/D");
    otree->Branch("beta", &beta, "beta/D");
    otree->Branch("dz", &dz, "dz/D");
    otree->Branch("enedep", &enedep, "enedep/D");

    TH2D *hist_bethe_bloch_iron  = new TH2D("hist_bethe_bloch_iron",  ";beta;Energy deposit [MeV/unit];", 1000, 0, 1, 1000, 0, 10);
    TH2D *hist_bethe_bloch_water = new TH2D("hist_bethe_bloch_water", ";beta;Energy deposit [MeV/unit];", 1000, 0, 1, 1000, 0, 10);

    TString filename = "/home/t2k/odagawa/data/mc_data/particlegun/particlegun_for_bethebloch/particlegun_muon_merge.root";
    B2Reader reader(filename);
    
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
      
      for (Int_t iemulsion = 0; iemulsion < emulsions.size(); iemulsion++) {
	
	const auto emulsion = emulsions.at(iemulsion);
	
	int track_id_ = emulsion->GetParentTrackId();
	int ecc_ = emulsion->GetEcc();
	
	if (track_id_ != track_id_tmp_ || ecc_ != ecc_tmp_) continue;
	// Use only energy deposit inside the tracking unit
	if (emulsion == emulsions.back()) continue;
	if (emulsion->GetPlate() < 16) continue;
	
	plate = emulsion->GetPlate();
	momentum_first = emulsion->GetMomentum().GetValue().Mag();
	energy_first = std::hypot(momentum_first, muon_mass);
	momentum_next = emulsions.at(iemulsion+1)->GetMomentum().GetValue().Mag();
	energy_next = std::hypot(momentum_next, muon_mass);
	beta = momentum_first / energy_first;

	dz = emulsion->GetTangent().GetValue().Mag();
	//enedep = emulsion->GetEdepSum();
	enedep = energy_first - energy_next;
	enedep /= dz;

	if (beta > 0.9) {
	  BOOST_LOG_TRIVIAL(debug) << "Momentum difference exceeds 100 MeV/c : "
				   << "Entry : " << reader.GetEntryNumber();
	  //std::exit(1);
	  continue;
	}
	
	if (emulsion->GetPlate()%2==0)
	  hist_bethe_bloch_iron->Fill(beta, enedep);
	else
	  hist_bethe_bloch_water->Fill(beta, enedep);
	
	otree->Fill();
      }      
      
    }
    
    TF1 *f_bethe_bloch = new TF1("f_bethe_bloch", "[0] / x / x * (log (x / (1 - x * x)) + [1]) + [2]");
    c->cd();
    hist_bethe_bloch_iron->Draw("colz");
    hist_bethe_bloch_iron->Fit(f_bethe_bloch, "", "");
    //hist_bethe_bloch->Fit(f_bethe_bloch, "", "", 0.4, 0.8);
    f_bethe_bloch->Draw("same");
    c->SaveAs(pdfname, "pdf");

    hist_bethe_bloch_water->Draw("colz");
    hist_bethe_bloch_water->Fit(f_bethe_bloch, "", "");
    //hist_bethe_bloch->Fit(f_bethe_bloch, "", "", 0.4, 0.8);
    f_bethe_bloch->Draw("same");
    c->SaveAs(pdfname, "pdf");

    c->SaveAs(pdfname + "]" , "pdf");
    
    ofile->cd();
    hist_bethe_bloch_iron->Write();
    hist_bethe_bloch_water->Write();
    otree->Write();
    ofile->Close();
    
  } catch (const std::runtime_error &error) {
    BOOST_LOG_TRIVIAL(fatal) << "Runtime error : " << error.what();
    std::exit(1);
  } catch (const std::invalid_argument &error) {
    BOOST_LOG_TRIVIAL(fatal) << "Invalid argument : " << error.what();
    std::exit(1);
  }

  BOOST_LOG_TRIVIAL(info) << "==========Bethe Bloch Fit Finish==========";
  std::exit(0);

}
