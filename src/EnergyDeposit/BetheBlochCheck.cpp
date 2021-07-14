// boost includes
#include <boost/log/core.hpp>
#include <boost/log/trivial.hpp>
#include <boost/log/expressions.hpp>

#include <boost/filesystem.hpp>
#include <boost/regex.hpp>

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
#include <TGraphErrors.h>

// system include
#include <iostream>

// my include 
#include "/home/t2k/odagawa/NinjaMomentumRecon/src/McsCommon.cpp"

namespace logging = boost::log;
namespace fs = boost::filesystem;

int main (int argc, char *argv[]) {

  logging::core::get()->set_filter
    (
     // logging::trivial::severity >= logging::trivial::info
     logging::trivial::severity >= logging::trivial::debug
     );
  
  BOOST_LOG_TRIVIAL(info) << "==========Bethe Bloch Fit Start==========";
  
  if (argc != 3) {
    BOOST_LOG_TRIVIAL(error) << "Usage : " << argv[0]
			     << " <output root file name> <output pdf file name>";
    std::exit(1);
  }
  
  try {
    
    gStyle->SetOptStat(0);
    gErrorIgnoreLevel = kWarning;
    const Double_t muon_mass = 105.658 * MeV;

    // Input B2 files    
    const std::string target_directory_path("/home/t2k/odagawa/data/mc_data/particlegun/particlegun_for_bethebloch_new/");
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
    TTree *otree = new TTree("tree", "tree");
    Double_t momentum_first, energy_first, momentum_next, energy_next;
    Double_t beta;
    Double_t dz;
    Double_t enedep;
    std::vector<Double_t> iron_param = {};
    std::vector<Double_t> water_param = {};
    otree->Branch("iron_param", &iron_param);
    otree->Branch("water_param", &water_param);

    // Output pdf
    TCanvas *c = new TCanvas("c", "c");
    TString pdfname = argv[2];
    c->Print(pdfname + "[", "pdf");


    std::vector<Double_t> enedep_iron, enedep_err_iron;
    std::vector<Double_t> enedep_water, enedep_err_water;
    std::vector<Double_t> beta_vec, beta_err_vec;

    TF1 *f_landau = new TF1("f_landau", "landau");

    // Loop for input files    
    for (auto filename : all_matched_files) {

      B2Reader reader(target_directory_path + filename);

      TH1D *hist_enedep_iron  = new TH1D("hist_enedep_iron",  "Energy Deposit Iron " + (TString)filename,  300, 0., 3.);
      TH1D *hist_enedep_water = new TH1D("hist_enedep_water", "Energy Deposit Water " + (TString)filename, 300, 0., 3.);

      double beta = 0;
      int nbeta = 0;

      while (reader.ReadNextSpill() > 0) {

	auto &spill_summary = reader.GetSpillSummary();

	// Get emulsion tracks
	std::vector<const B2EmulsionSummary*> emulsions;
	auto it_emulsion = spill_summary.BeginEmulsion();
	while (const auto *emulsion = it_emulsion.Next()) {
	  if (emulsion->GetParentTrackId() == 0) continue;
	  if (emulsion->GetParentTrack().GetParticlePdg() != 13) continue;
	  if (emulsion->GetFilmType() != B2EmulsionType::kECC) continue;
	  if (emulsion->GetEcc() != 4) continue;
	  if (emulsion->GetPlate() < 128) continue;
	  emulsions.push_back(emulsion);
	}

	if (emulsions.size() < 3) continue;
	
	std::sort(emulsions.begin(), emulsions.end(), emulsion_compare);

	for (Int_t iemulsion = 0; iemulsion < emulsions.size(); iemulsion++) {

	  const auto emulsion = emulsions.at(iemulsion);
	  // Use only energy deposit inside the tracking unit
	  if (emulsion == emulsions.back()) continue;
	  
	  momentum_first = emulsion->GetMomentum().GetValue().Mag();
	  energy_first = std::hypot(momentum_first, muon_mass);
	  momentum_next = emulsions.at(iemulsion+1)->GetMomentum().GetValue().Mag();
	  energy_next = std::hypot(momentum_next, muon_mass);
	  beta += momentum_first / energy_first;
	  nbeta++;
	  
	  dz = emulsion->GetTangent().GetValue().Mag();
	  //enedep = emulsion->GetEdepSum();
	  enedep = energy_first - energy_next;
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
      hist_enedep_iron->Draw("");
      hist_enedep_iron->Fit(f_landau);
      f_landau->Draw("SAME");
      c->SaveAs(pdfname, "pdf");

      enedep_iron.push_back(f_landau->GetParameter(1));
      enedep_err_iron.push_back(f_landau->GetParError(1));

      hist_enedep_water->Draw("");
      hist_enedep_water->Fit(f_landau);
      f_landau->Draw("SAME");
      c->SaveAs(pdfname, "pdf");

      enedep_water.push_back(f_landau->GetParameter(1));
      enedep_err_water.push_back(f_landau->GetParError(1));

      beta_vec.push_back(beta / nbeta);
      beta_err_vec.push_back(0.);

      //if (filename == all_matched_files.at(1))
      //break;

    }

    TF1 *f_bethe_bloch = new TF1("f_bethe_bloch", "[0] / x / x * (log (x / (1 - x * x)) + [1]) + [2]");
    f_bethe_bloch->SetParLimits(0, 0., 100.);
    const Int_t number_of_files = all_matched_files.size();
    TGraphErrors *ge_iron = new TGraphErrors(number_of_files,
					     &beta_vec[0], &enedep_iron[0],
					     &beta_err_vec[0], &enedep_err_iron[0]);
    ge_iron->SetTitle("Bethe Bloch function across one iron plate;#beta;Energy deposit [MeV/unit]");
    ge_iron->GetYaxis()->SetRangeUser(0., 5.);
    TGraphErrors *ge_water = new TGraphErrors(number_of_files,
					      &beta_vec[0], &enedep_water[0],
					      &beta_err_vec[0], &enedep_err_water[0]);
    ge_water->SetTitle("Bethe Bloch fuction across one water layer;#beta;Energy deposit [MeV/unit]");
    ge_water->GetYaxis()->SetRangeUser(0., 5.);

    c->cd();
    ge_iron->Draw("AP");
    c->SaveAs(pdfname, "pdf");
    ge_iron->Fit(f_bethe_bloch, "", "", 0.7, 1.);
    f_bethe_bloch->Draw("SAME");
    for (Int_t i = 0; i < 3; i++)
      iron_param.push_back(f_bethe_bloch->GetParameter(i));
    c->SaveAs(pdfname, "pdf");

    ge_water->Draw("AP");
    c->SaveAs(pdfname, "pdf");
    ge_water->Fit(f_bethe_bloch, "", "", 0.7, 1.);
    f_bethe_bloch->Draw("SAME");
    for (Int_t i = 0; i < 3; i++)
      water_param.push_back(f_bethe_bloch->GetParameter(i));
    c->SaveAs(pdfname, "pdf");

    c->SaveAs(pdfname + "]" , "pdf");
    
    ofile->cd();
    ge_iron->Write();
    ge_water->Write();
    otree->Fill();
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
