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
#include <B2Pdg.hh>

// ROOT includes
#include <TStyle.h>
#include <TCanvas.h>
#include <TFile.h>
#include <TH2D.h>
#include <TF1.h>

// my include
#include "McsConst.hpp"
#include "McsFunction.hpp"

namespace logging = boost::log;

int main (int argc, char *argv[]) {

  logging::core::get()->set_filter
    (
     //logging::trivial::severity >= logging::trivial::info
     logging::trivial::severity >= logging::trivial::debug
     );

  BOOST_LOG_TRIVIAL(info) << "==========Angle Difference Fit Start==========";

  if ( argc != 6 ) {
    BOOST_LOG_TRIVIAL(error) << "Usage :" << argv[0]
			     << " <output rootfile name> <output pdf name> <mode (lateral(0)/ax(1)/ay(2)/radial angle(3))> <certain topview angle> <certain sideview angle>";
    std::exit(1);
  }

  const Int_t analysis_mode = std::atoi(argv[3]);
  if (analysis_mode != 0 && analysis_mode != 1 && analysis_mode != 2) {
    BOOST_LOG_TRIVIAL(error) << "Analysis mode should be 0 - 2";
    std::exit(1);
  }

  try {

    gStyle->SetOptStat(0);
    gErrorIgnoreLevel = kWarning;

    TCanvas *c = new TCanvas("c", "c");
    TString pdfname = argv[2];

    c->Print(pdfname + "[", "pdf");

    Double_t momentum, sideview, topview;
    Double_t x, y, z;
    z = -183.08;
    Double_t angle, unit_path_length;
    Double_t center_lat = 0;
    Double_t center_lat_err = 0;
    Double_t sigma_lat = 0;
    Double_t sigma_lat_err = 0;
    Double_t pbeta, beta;

    TString ofilename = argv[1];
    TFile *ofile = new TFile(ofilename, "recreate");
    TTree *otree = new TTree("tree", "tree");
    otree->Branch("momentum", &momentum, "momentum/D");
    otree->Branch("pbeta", &pbeta, "pbeta/D");
    otree->Branch("beta", &beta, "beta/D");
    otree->Branch("sideview", &sideview, "sideview/D");
    otree->Branch("topview", &topview, "topview/D");
    otree->Branch("angle", &angle, "angle/D");
    otree->Branch("unit_path_length", &unit_path_length, "unit_path_length/D");
    otree->Branch("center_lat", &center_lat, "center_lat/D");
    otree->Branch("center_lat_err", &center_lat_err, "center_lat_err/D");
    otree->Branch("sigma_lat", &sigma_lat, "sigma_lat/D");
    otree->Branch("sigma_lat_err", &sigma_lat_err, "sigma_lat_err/D");

    TString idirname = "/home/t2k/odagawa/data/mc_data/particlegun/particlegun_for_deltatheta";

    for ( Int_t imomentum = 1; imomentum < 30; imomentum++ ) { // 50 MeV is not good to fit by Gauss
      for ( Int_t isideview = 0; isideview < 10; isideview++ ) {
	for ( Int_t itopview = 0; itopview < 10; itopview++ ) {
	  
	  if ( std::atoi(argv[4]) >= 0 ) {
	    if ( itopview != std::atoi(argv[4]) ) continue;
	  }
	  
	  if ( std::atoi(argv[5]) >= 0 ) {
	    if ( isideview != std::atoi(argv[5]) ) continue;
	  }

	  momentum = (imomentum + 1) * 50;
	  sideview = isideview * 5;
	  topview = itopview * 5;

	  angle = TMath::Hypot(TMath::Tan(sideview * TMath::RadToDeg()), TMath::Tan(topview * TMath::RadToDeg()));
	  unit_path_length = TMath::Hypot(angle, 1.);
	  angle = TMath::ATan(angle);

	  if ( sideview == 0 ) {
	    y = -12.4;
	  } else {
	    y = -132.4;
	  }

	  if ( topview == 0 ) {
	    x = -350;
	  } else {
	    x = -470;
	  }

	  Double_t energy = TMath::Hypot(momentum, MCS_MUON_MASS);
	  beta = momentum / energy;
	  pbeta = momentum * momentum / energy;

	  TString ifilename = idirname + "/";
	  ifilename += Form("particlegun_muon_%dMeV_%d_%.1f_%.2f_%d_%d.root",
			    (Int_t)momentum, (Int_t)x, y, z, (Int_t)sideview, (Int_t)topview);
	  B2Reader reader(ifilename);

	  Double_t maxbin;
	  if ( imomentum == 0 ) maxbin = 0.2;
	  else if ( imomentum < 3 ) maxbin = 0.1;
	  else if ( imomentum < 9 ) maxbin = 0.05;
	  else if ( imomentum < 25 ) maxbin = 0.02;
	  else maxbin = 0.01;
    
	  TH1D *hist_angle_difference = new TH1D("hist", "", 100, -maxbin, maxbin);    
	  TF1 *gaus = new TF1("gaus","gaus");

	  while ( reader.ReadNextSpill() > 0 ) {
	    auto &input_spill_summary = reader.GetSpillSummary();

	    // Get the most upstream emulsion tracks
	    std::vector<const B2EmulsionSummary* > emulsions;
	    auto it_emulsion = input_spill_summary.BeginEmulsion();
	    while ( const auto *emulsion = it_emulsion.Next() ) {
	      if ( emulsion->GetParentTrackId() == 0 ) continue;
	      if ( emulsion->GetParentTrack().GetParticlePdg() != PDG_t::kMuonMinus) continue;
	      if ( emulsion->GetFilmType() != B2EmulsionType::kECC) continue;
	      if ( emulsion->GetEcc() != 4 ) continue;
	      if ( emulsion->GetPlate() != 130 && emulsion->GetPlate() != 129 ) continue;
	      emulsions.push_back(emulsion);
	    }

	    if ( emulsions.size() != 2 ) continue;

	    std::sort(emulsions.begin(), emulsions.end(), EmulsionCompare);

	    if ( emulsions.at(0)->GetPlate() != 130 ||
		 emulsions.at(1)->GetPlate() != 129 ) continue;

	    const auto emulsion_up = emulsions.at(0);
	    const auto emulsion_down = emulsions.at(1);

	    TVector3 tangent_up = emulsion_up->GetTangent().GetValue();
	    TVector3 tangent_down = emulsion_down->GetTangent().GetValue();
	    if ( sideview == 0 && topview == 0 ) {
	      if ( analysis_mode != 2 )
		hist_angle_difference->Fill(TMath::ATan(tangent_up.X()) - TMath::ATan(tangent_down.X()));
	      else if ( analysis_mode == 2 )
		hist_angle_difference->Fill(TMath::ATan(tangent_up.Y()) - TMath::ATan(tangent_down.Y()));
	    } else {
	      if ( analysis_mode == 0 ) {
		hist_angle_difference->Fill(RadialAngleDiffNew(tangent_up, tangent_down));
	      }
	      else if ( analysis_mode == 1 )
		hist_angle_difference->Fill(TMath::ATan(tangent_up.X()) - TMath::ATan(tangent_down.X()));
	      else if ( analysis_mode == 2 )
	      hist_angle_difference->Fill(TMath::ATan(tangent_up.Y()) - TMath::ATan(tangent_down.Y()));
	    }

	    emulsions.clear();
	    emulsions.shrink_to_fit();

	  }

	  c->cd();
	  if ( analysis_mode == 0 ) {
	    hist_angle_difference->SetTitle(Form("Lateral angle difference (momentum = %d MeV/c, #theta_{x} = %d deg., #theta_{y} = %d deg.);#Delta#theta_{lat};Entries",
						 (Int_t)momentum, (Int_t)topview, (Int_t)sideview));
	  }
	  else if ( analysis_mode == 1 ) {
	    hist_angle_difference->SetTitle(Form("Angle difference (momentum = %d MeV/c, #theta_{x} = %d deg., #theta_{y} = %d deg.);#Delta#theta_{x};Entries",
						 (Int_t)momentum, (Int_t)topview, (Int_t)sideview));
	  }
	  else if ( analysis_mode == 2 ) {
	    hist_angle_difference->SetTitle(Form("Angle difference (momentum = %d MeV/c, #theta_{x} = %d deg., #theta_{y} = %d deg.);#Delta#theta_{y};Entries",
						 (Int_t)momentum, (Int_t)topview, (Int_t)sideview));
	  }
	  hist_angle_difference->Draw();
	  Double_t stddev = hist_angle_difference->GetStdDev();
	  hist_angle_difference->Fit(gaus, "", "", -0.9 * stddev, 0.9 * stddev);
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
