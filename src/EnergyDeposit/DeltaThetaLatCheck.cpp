#include <iostream>

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
#include "/home/t2k/odagawa/NinjaMomentumRecon/src/McsAngleMethod.cpp"

namespace logging = boost::log;

int main (int argc, char *argv[]) {

  logging::core::get()->set_filter
    (
     //logging::trivial::severity >= logging::trivial::info
     logging::trivial::severity >= logging::trivial::debug
     );

  BOOST_LOG_TRIVIAL(info) << "==========Angle Difference Fit Start==========";

  if (argc != 1) {
    BOOST_LOG_TRIVIAL(error) << "Usage :" << argv[0];
    std::exit(1);
  }
  
  try {
    
    gStyle->SetOptStat(0);
    gErrorIgnoreLevel = kWarning;
    
    const Double_t muon_mass = 105.658;//MeV/c2
    
    TCanvas *c = new TCanvas("c", "c");
    TString pdfname = "angle_difference_fit_radcheck_ang.pdf";
    c->Print(pdfname + "[", "pdf");
    
    Double_t momentum, sideview, topview;
    Double_t x, y, z;
    z = -183.08;
    Double_t angle, unit_path_length;
    Double_t x0, x1, y0, y1, z0, z1;
    Double_t center_lat = 0;
    Double_t center_lat_err = 0;
    Double_t sigma_lat = 0;
    Double_t sigma_lat_err = 0;
    Double_t pbeta, beta;
    
    TFile *ofile = new TFile("angle_difference_fit_radcheck_ang.root", "recreate");
    TTree *otree = new TTree("tree", "tree");
    otree->Branch("momentum", &momentum, "momentum/D");
    otree->Branch("pbeta", &pbeta, "pbeta/D");
    otree->Branch("beta", &beta, "beta/D");
    otree->Branch("sideview", &sideview, "sideview/D");
    otree->Branch("topview", &topview, "topview/D");
    otree->Branch("angle", &angle, "angle/D");
    /*
    otree->Branch("x0", &x0, "x0/D");
    otree->Branch("y0", &y0, "y0/D");
    otree->Branch("z0", &z0, "z0/D");
    otree->Branch("x1", &x1, "x1/D");
    otree->Branch("y1", &y1, "y1/D");
    otree->Branch("z1", &z1, "z1/D");
    */
    otree->Branch("unit_path_length", &unit_path_length, "unit_path_length/D");
    otree->Branch("center_lat", &center_lat, "center_lat/D");
    otree->Branch("center_lat_err", &center_lat_err, "center_lat_err/D");
    otree->Branch("sigma_lat", &sigma_lat, "sigma_lat/D");
    otree->Branch("sigma_lat_err", &sigma_lat_err, "sigma_lat_err/D");

    TString idirname = "/home/t2k/odagawa/data/mc_data/particlegun/particlegun_for_deltatheta";

    for (Int_t imomentum = 1; imomentum < 30; imomentum++) { // 50 MeV is not good to fit by Gauss
      for (Int_t isideview = 0; isideview < 10; isideview++) {
	for (Int_t itopview = 0; itopview < 10; itopview++) {
	  momentum = (imomentum + 1) * 50;
	  sideview = isideview * 5;
	  topview = itopview * 5;
	  
	  angle = TMath::Hypot(TMath::Tan(sideview * TMath::DegToRad()), TMath::Tan(topview * TMath::DegToRad()));
	  unit_path_length = TMath::Hypot(angle, 1.);
	  angle = TMath::ATan(angle);
	  
	  if (sideview == 0) {
	    y = -12.4;
	  } else {
	    y = -132.4;
	  }
	  
	  if (topview == 0) {
	    x = -350;
	  } else {
	    x = -470;
	  }
	  
	  Double_t energy = TMath::Hypot(momentum, muon_mass);
	  beta = momentum / energy;
	  pbeta = momentum * momentum / energy;
	  
	  TString ifilename = idirname + "/";
	  ifilename += Form("particlegun_muon_%dMeV_%d_%.1f_%.2f_%d_%d.root",
			    (Int_t)momentum, (Int_t)x, y, z, (Int_t)sideview, (Int_t)topview);
	  B2Reader reader(ifilename);
	  
	  Double_t maxbin;
	  if (imomentum == 0) maxbin = 0.5;
	  else if (imomentum < 3) maxbin = 0.5;
	  else if (imomentum < 9) maxbin = 0.1;
	  else if (imomentum < 25) maxbin = 0.05;
	  else maxbin = 0.02;
    
	  TH1D *hist_angle_difference = new TH1D("hist", "", 100, -maxbin, maxbin);    
	  TF1 *gaus = new TF1("gaus","gaus");
	  
	  while (reader.ReadNextSpill() > 0) {
	    auto &input_spill_summary = reader.GetSpillSummary();
	    
	    // Get the most upstream emulsion tracks
	    std::vector<const B2EmulsionSummary*> emulsions;
	    auto it_emulsion = input_spill_summary.BeginEmulsion();
	    while (const auto *emulsion = it_emulsion.Next()) {
	      if (emulsion->GetParentTrackId() == 0) continue;
	      if (emulsion->GetParentTrack().GetParticlePdg() != 13) continue;
	      if (emulsion->GetFilmType() != B2EmulsionType::kECC) continue;
	      if (emulsion->GetEcc() != 4 ) continue;
	      if (emulsion->GetPlate() != 130 && emulsion->GetPlate() != 129) continue;
	      emulsions.push_back(emulsion);
	    }
	    
	    if (emulsions.size() != 2) continue;
	    
	    std::sort(emulsions.begin(), emulsions.end(), emulsion_compare);
	    
	    if (emulsions.at(0)->GetPlate() != 130 ||
		emulsions.at(1)->GetPlate() != 129) continue;
	    
	    const auto emulsion_up = emulsions.at(0);
	    const auto emulsion_down = emulsions.at(1);
	    
	    TVector3 tangent_up = emulsion_up->GetTangent().GetValue();
	    TVector3 tangent_down = emulsion_down->GetTangent().GetValue();
	    if (topview == 0)
	      //hist_angle_difference->Fill(TMath::ATan(tangent_up.X()) - TMath::ATan(tangent_down.X()));
	      hist_angle_difference->Fill(TMath::ATan(tangent_up.Y()) - TMath::ATan(tangent_down.Y()));
	    else if (sideview == 0)
	      //hist_angle_difference->Fill(TMath::ATan(tangent_up.Y()) - TMath::ATan(tangent_down.Y()));
	      hist_angle_difference->Fill(TMath::ATan(tangent_up.X()) - TMath::ATan(tangent_down.X()));
	    else {
	      //hist_angle_difference->Fill(get_angle_difference_lateral(tangent_up, tangent_down, tangent_up));
	      //hist_angle_difference->Fill(TMath::ATan(get_tangent_difference_lateral(tangent_up, tangent_down, tangent_up)));
	      //hist_angle_difference->Fill(TMath::ATan(get_tangent_difference_radial(tangent_up, tangent_down)));
	      hist_angle_difference->Fill(get_angle_difference_radial(tangent_up, tangent_down));
	    }
	    emulsions.clear();
	    emulsions.shrink_to_fit();

	  }

	  c->cd();
	  //hist_angle_difference->SetTitle(Form("Lateral angle difference (momentum = %d MeV/c, #theta_{x} = %d deg., #theta_{y} = %d deg.);#Delta#theta_{lat};Entries",
	  //				       (Int_t)momentum, (Int_t)topview, (Int_t)sideview));
	  hist_angle_difference->SetTitle(Form("Radial angle difference (momentum = %d MeV/c, #theta_{x} = %d deg., #theta_{y} = %d deg.);#Delta#theta_{rad};Entries",
	  				       (Int_t)momentum, (Int_t)topview, (Int_t)sideview));      
	  hist_angle_difference->Draw();
	  Double_t stddev = hist_angle_difference->GetStdDev();
	  hist_angle_difference->Fit(gaus, "", "", -stddev, stddev);
	  gaus->Draw("same");
	  c->Print(pdfname, "pdf");

	  center_lat = gaus->GetParameter(1);
	  center_lat_err = gaus->GetParError(1);
	  sigma_lat = gaus->GetParameter(2);
	  sigma_lat_err = gaus->GetParError(2);

	  otree->Fill();

	  delete hist_angle_difference;
	  delete gaus;
	  
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
  }

  std::exit(0);

}
