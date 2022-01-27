// boost includes
#include <boost/log/core.hpp>
#include <boost/log/trivial.hpp>
#include <boost/log/expressions.hpp>

#include <boost/filesystem.hpp>
#include <boost/regex.hpp>

// B2 includes
#include <B2Pdg.hh>
#include <B2Reader.hh>
#include <B2SpillSummary.hh>
#include <B2TrackSummary.hh>
#include <B2EmulsionSummary.hh>

// ROOT includes
#include <TStyle.h>
#include <TPad.h>
#include <TCanvas.h>
#include <TLine.h>
#include <TFile.h>
#include <TH2D.h>
#include <TF1.h>
#include <TGraphErrors.h>

// system include
#include <iostream>

// my include 
#include "McsConst.hpp"
#include "McsFunction.hpp"

namespace logging = boost::log;
namespace fs = boost::filesystem;

int main (int argc, char *argv[]) {

  logging::core::get()->set_filter
    (
     logging::trivial::severity >= logging::trivial::info
     //logging::trivial::severity >= logging::trivial::debug
     );
  
  BOOST_LOG_TRIVIAL(info) << "==========Bethe Bloch Fit Start==========";
  
  if ( argc != 4 ) {
    BOOST_LOG_TRIVIAL(error) << "Usage : " << argv[0]
			     << " <output root file name> <output pdf file name> <particle id>";
    std::exit(1);
  }
  
  try {
    
    gStyle->SetOptStat(0);
    gErrorIgnoreLevel = kWarning;

    Int_t particle_id = std::atoi(argv[3]);
    std::string particle_name;
    Double_t particle_mass;
    if ( B2Pdg::IsMuonPlusOrMinus(particle_id) ) {
      particle_name = "muon";
      particle_mass = MCS_MUON_MASS;
    }
    else if ( B2Pdg::IsChargedPion(particle_id) ) {
      particle_name = "pion";
      particle_mass = MCS_PION_MASS;
    }
    else if ( particle_id == PDG_t::kProton ) {
      particle_name = "proton";
      particle_mass = MCS_PROTON_MASS;
    }
    else
      throw std::invalid_argument("Particle ID is not in interest");

    // Input B2 files    
    const std::string target_directory_path("/home/t2k/odagawa/data/mc_data/particlegun/particlegun_for_bethebloch_new/" + particle_name + "/");
    const std::string root_extension(".root");

    std::vector<std::string> all_matched_files;

    fs::directory_iterator end_itr;
    for( fs::directory_iterator i( target_directory_path ); i != end_itr; ++i ) {
      // Skip if not a file
      if ( !fs::is_regular_file( i->status() ) ) continue;

      if ( i->path().extension().string() != root_extension ) continue;

      BOOST_LOG_TRIVIAL(debug) << i->path().filename().string();

      // File matches, store it
      all_matched_files.push_back( i->path().filename().string() );
    }

    std::sort(all_matched_files.begin(), all_matched_files.end());

    // Output root file
    TFile *ofile = new TFile(argv[1],"recreate");
    TTree *otree = new TTree("o_tree", "o_tree");
    TTree *paramtree = new TTree("p_tree", "p_tree");
    Double_t iron_edep;
    Double_t iron_edep_err;
    Double_t water_edep;
    Double_t water_edep_err;
    Double_t beta;
    otree->Branch("iron_edep", &iron_edep, "iron_edep/D");
    otree->Branch("iron_edep_err", &iron_edep_err, "iron_edep_err/D");
    otree->Branch("water_edep", &water_edep, "water_edep/D");
    otree->Branch("water_edep_err", &water_edep_err, "water_edep_err/D");
    otree->Branch("beta", &beta, "beta/D");


    std::vector<Double_t > iron_param = {};
    std::vector<Double_t > water_param = {};
    paramtree->Branch("iron_param", &iron_param);
    paramtree->Branch("water_param", &water_param);

    // Output pdf
    TCanvas *c = new TCanvas("c", "c");
    TString pdfname = argv[2];
    c->Print(pdfname + "[", "pdf");


    std::vector<Double_t > iron_edep_vec, iron_edep_err_vec;
    std::vector<Double_t > water_edep_vec, water_edep_err_vec;
    std::vector<Double_t > beta_vec, beta_err_vec;

    TF1 *f_landau = new TF1("f_landau", "landau");
    TLine *line_mean = new TLine();
    line_mean->SetLineColor(kRed);
    line_mean->SetLineWidth(3);

    // Loop for input files    
    for ( auto filename : all_matched_files ) {

      B2Reader reader(target_directory_path + filename);
      
      TH1D *hist_enedep_iron  = new TH1D("hist_enedep_iron",
					 "Energy Deposit Iron "  + (TString)filename + ";Energy deposit [MeV];Entries",
					 2000, 0., 20.);
      TH1D *hist_enedep_water = new TH1D("hist_enedep_water",
					 "Energy Deposit Water " + (TString)filename + ";Energy deposit [MeV];Entries",
					 2000, 0., 20.);
      
      beta = 0;
      int nbeta = 0;

      while ( reader.ReadNextSpill() > 0 ) {

	auto &spill_summary = reader.GetSpillSummary();

	// Get emulsion tracks
	std::vector<const B2EmulsionSummary* > emulsions;
	auto it_emulsion = spill_summary.BeginEmulsion();
	while ( const auto *emulsion = it_emulsion.Next() ) {
	  if ( emulsion->GetParentTrackId() == 0 ) continue;
	  if ( emulsion->GetParentTrack().GetParticlePdg() != particle_id ) continue;
	  if ( emulsion->GetFilmType() != B2EmulsionType::kECC ) continue;
	  if ( emulsion->GetEcc() != 4 ) continue;
	  if ( emulsion->GetPlate() < 127 ) continue;
	  emulsions.push_back(emulsion);
	}

	if ( emulsions.size() < 4 ) continue;
	
	std::sort(emulsions.begin(), emulsions.end(), EmulsionCompare);

	for ( Int_t iemulsion = 0; iemulsion < emulsions.size(); iemulsion++ ) {

	  const auto emulsion = emulsions.at(iemulsion);
	  // Use only energy deposit inside the tracking unit
	  if ( emulsion == emulsions.back() ) continue;
	  
	  Double_t momentum_first = emulsion->GetMomentum().GetValue().Mag();
	  Double_t energy_first = std::hypot(momentum_first, particle_mass);
	  Double_t beta_first = momentum_first / energy_first;
	  Double_t momentum_next = emulsions.at(iemulsion+1)->GetMomentum().GetValue().Mag();
	  Double_t energy_next = std::hypot(momentum_next, particle_mass);
	  Double_t beta_next = momentum_next / energy_next;
	  // beta += beta_first;
	  beta += beta_next;
	  // beta += (beta_first + beta_next) / 2.
	  nbeta++;
	  
	  Double_t dz = emulsion->GetTangent().GetValue().Mag();
	  //enedep = emulsion->GetEdepSum();
	  Double_t enedep = energy_first - energy_next;
	  enedep /= dz;

	  BOOST_LOG_TRIVIAL(debug) << "Upstream plate : " << emulsion->GetPlate() << ", "
				   << "Downstream pate : " << emulsions.at(iemulsion+1)->GetPlate() << ", "
				   << "Energy deposit : " << enedep << " MeV/unit";

	  if (emulsion->GetPlate() == 130 &&
	      emulsions.at(iemulsion+1)->GetPlate() == 129) { // Across one iron plate
	    hist_enedep_iron->Fill(enedep);
	  } else if (emulsion->GetPlate() == 129 &&
		     emulsions.at(iemulsion+1)->GetPlate() == 128) { // Across one water layer
	    hist_enedep_water->Fill(enedep);
	  } else continue;
	  	 
	}

      }

      c->cd();

      if ( hist_enedep_iron->GetEntries() < 9000. ||
	   hist_enedep_water->GetEntries() < 9000. ) continue;

      // Iron energy deposit
      hist_enedep_iron->Draw("");
      iron_edep = hist_enedep_iron->GetMean();
      iron_edep_err = hist_enedep_iron->GetMeanError();

      f_landau->SetParameter(1, 0.5);
      hist_enedep_iron->Fit(f_landau);
      gPad->Update(); // Call this function before GetUymin()/GetUymax()
      line_mean->DrawLine(iron_edep, gPad->GetUymin(),
			  iron_edep, gPad->GetUymax());
      c->SaveAs(pdfname, "pdf");
      
      iron_edep_vec.push_back(iron_edep);
      iron_edep_err_vec.push_back(iron_edep_err);

      // Water energy deposit
      hist_enedep_water->Draw("");
      water_edep = hist_enedep_water->GetMean();
      water_edep_err = hist_enedep_water->GetMeanError();

      f_landau->SetParameter(1, 0.5);
      hist_enedep_water->Fit(f_landau);
      gPad->Update();
      line_mean->DrawLine(water_edep, gPad->GetUymin(),
			  water_edep, gPad->GetUymax());
      c->SaveAs(pdfname, "pdf");

      water_edep_vec.push_back(water_edep);
      water_edep_err_vec.push_back(water_edep_err);

      // Beta
      beta /= (Double_t)nbeta;
      beta_vec.push_back(beta);
      beta_err_vec.push_back(0.);

      //if (filename == all_matched_files.at(1))
      //break;

      otree->Fill();

    }

    TF1 *f_bethe_bloch = new TF1("f_bethe_bloch", "[0] * ( (log(x*x / (1-x*x)) + [1]) / x / x - 1 )");
    const Int_t number_of_files = beta_vec.size();
    TGraphErrors *ge_iron = new TGraphErrors(number_of_files,
					     &beta_vec[0], &iron_edep_vec[0],
					     &beta_err_vec[0], &iron_edep_err_vec[0]);
    ge_iron->SetTitle("Bethe Bloch function across one iron plate;#beta;Energy deposit [MeV/unit]");
    ge_iron->GetYaxis()->SetRangeUser(0.5, 10.);
    TGraphErrors *ge_water = new TGraphErrors(number_of_files,
					      &beta_vec[0], &water_edep_vec[0],
					      &beta_err_vec[0], &water_edep_err_vec[0]);
    ge_water->SetTitle("Bethe Bloch fuction across one water layer;#beta;Energy deposit [MeV/unit]");
    ge_water->GetYaxis()->SetRangeUser(0.5, 10.);

    c->cd();

    // Iron energy deposit fitting
    ge_iron->Draw("AP");
    c->SaveAs(pdfname, "pdf");

    f_bethe_bloch->SetParameter(0, 6.e-2);
    f_bethe_bloch->SetParameter(1, 8.);
    ge_iron->Fit(f_bethe_bloch, "", "", 0.3, 0.9);
    f_bethe_bloch->Draw("SAME");
    for (Int_t i = 0; i < 2; i++)
      iron_param.push_back(f_bethe_bloch->GetParameter(i));
    c->SaveAs(pdfname, "pdf");

    // Water energy deposit fitting
    ge_water->Draw("AP");
    c->SaveAs(pdfname, "pdf");
    ge_water->Fit(f_bethe_bloch, "", "", 0.3, 0.9);
    f_bethe_bloch->Draw("SAME");
    for (Int_t i = 0; i < 2; i++)
      water_param.push_back(f_bethe_bloch->GetParameter(i));
    c->SaveAs(pdfname, "pdf");

    c->SaveAs(pdfname + "]" , "pdf");
    
    ofile->cd();
    ge_iron->SetName("ge_iron");
    ge_iron->Write();
    ge_water->SetName("ge_water");
    ge_water->Write();
    paramtree->Fill();
    paramtree->Write();
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
