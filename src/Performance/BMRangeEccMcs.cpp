// boost includes
#include <boost/log/core.hpp>
#include <boost/log/trivial.hpp>
#include <boost/log/expressions.hpp>

// ROOT includes
#include <TString.h>
#include <TFile.h>
#include <TTree.h>
#include <TMath.h>

// system includes
#include <vector>
#include <algorithm>
#include <iostream>

// my includes
#include <NTBMSummary.hh>

namespace logging = boost::log;

int main (int argc, char* argv[]) {

  logging::core::get()->set_filter
    (
     logging::trivial::severity >= logging::trivial::info
     // logging::trivial::severity >= logging::trivial::debug
     );

  BOOST_LOG_TRIVIAL(info) << "==========Range vs. MCS Check Start==========";

  if (argc != 5) {
    BOOST_LOG_TRIVIAL(error) << "Usage : " << argv[0]
			     << " <input NTBM file name> <input MCS file name> <output root file name> <MC(0)/Data(1)>";
    std::exit(1);
  }

  try {

    // Input NTBM file
    TString ntbmfilename = argv[1];
    TFile *ntbmfile = new TFile(ntbmfilename, "read");
    TTree *ntbmtree = (TTree*)ntbmfile->Get("tree");
    NTBMSummary *ntbm = nullptr;
    ntbmtree->SetBranchAddress("NTBMSummary", &ntbm);

    // Input MCS file
    TString mcsfilename = argv[2];
    TFile *mcsfile = new TFile(mcsfilename, "read");
    TTree *mcstree = (TTree*)mcsfile->Get("tree");
    Int_t entry_in_daily_file;
    Int_t npl;
    Int_t muon_track_id = 0;
    Double_t true_pbeta;
    Double_t initial_pbeta;
    Double_t recon_pbeta;
    Double_t recon_pbeta_err;
    Double_t log_likelihood[3][2] = {}; // 3 particles (mu, pi, p), 2 directions (up, down)
    Double_t recon_pbeta_candidate[3][2] = {};
    Double_t recon_pbeta_err_candidate[3][2] = {};
    Int_t fit_status[3][2] = {};
    mcstree->SetBranchAddress("entry_in_daily_file", &entry_in_daily_file);
    mcstree->SetBranchAddress("npl", &npl);
    //mcstree->SetBranchAddress("muon_track_id", &muon_track_id);
    mcstree->SetBranchAddress("true_pbeta",  &true_pbeta);
    mcstree->SetBranchAddress("initial_pbeta", &initial_pbeta);
    mcstree->SetBranchAddress("recon_pbeta", &recon_pbeta);
    mcstree->SetBranchAddress("recon_pbeta_err", &recon_pbeta_err);
    mcstree->SetBranchAddress("log_likelihood", log_likelihood);
    mcstree->SetBranchAddress("recon_pbeta_candidate", recon_pbeta_candidate);
    mcstree->SetBranchAddress("recon_pbeta_err_candidate", recon_pbeta_err_candidate);
    mcstree->SetBranchAddress("fit_status", fit_status);

    TString ofilename = argv[3];
    TFile *ofile = new TFile(ofilename, "recreate");
    TTree *otree = new TTree("tree", "tree");
    Double_t mcs_pbeta, mcs_momentum, best_log_likelihood;
    Double_t range_pbeta, range_momentum;
    otree->Branch("entry_in_daily_file", &entry_in_daily_file, "entry_in_daily_file/I");
    otree->Branch("npl", &npl, "npl/I");
    otree->Branch("mcs_pbeta", &mcs_pbeta, "mcs_pbeta/D");
    otree->Branch("mcs_momentum", &mcs_momentum, "mcs_momentum/D");
    otree->Branch("best_log_likelihood", &best_log_likelihood, "best_log_likelihood/D");
    otree->Branch("range_pbeta", &range_pbeta, "range_pbeta/D");
    otree->Branch("range_momentum", &range_momentum, "range_momentum/D");


    for ( Int_t imcsentry = 0; imcsentry < mcstree->GetEntries(); imcsentry++ ) {
      mcstree->GetEntry(imcsentry);
      ntbmtree->GetEntry(entry_in_daily_file - 1);
      if (ntbm->GetEntryInDailyFile() != entry_in_daily_file)
	BOOST_LOG_TRIVIAL(debug) << "Entry in daily file is not the same!";

      mcs_pbeta = recon_pbeta;
      Double_t mass = 105.658;
      Double_t mcs_energy = 0.5 * (mcs_pbeta + TMath::Hypot(mcs_pbeta, 2. * mass));
      mcs_momentum = TMath::Sqrt(mcs_energy * mcs_energy - mass * mass);

      best_log_likelihood = log_likelihood[0][0];

      if ( ntbm->GetMomentumType(muon_track_id) != 0 ) continue;
      range_momentum = ntbm->GetMomentum(muon_track_id);
      Double_t range_energy = TMath::Hypot(range_momentum, mass);
      range_pbeta = range_momentum * range_momentum / range_energy;

      otree->Fill();

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

  BOOST_LOG_TRIVIAL(info) << "==========Range vs MCS Check Finish==========";
  std::exit(0);


}
