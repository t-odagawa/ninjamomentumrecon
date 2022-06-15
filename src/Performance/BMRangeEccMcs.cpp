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
#include "McsConst.hpp"
#include "McsFunction.hpp"

namespace logging = boost::log;

int main (int argc, char* argv[]) {

  logging::core::get()->set_filter
    (
     //logging::trivial::severity >= logging::trivial::info
      logging::trivial::severity >= logging::trivial::debug
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
    Double_t mcs_pbeta_err;
    Double_t range_pbeta, range_momentum;
    otree->Branch("entry_in_daily_file", &entry_in_daily_file, "entry_in_daily_file/I");
    otree->Branch("npl", &npl, "npl/I");
    otree->Branch("mcs_pbeta", &mcs_pbeta, "mcs_pbeta/D");
    otree->Branch("mcs_momentum", &mcs_momentum, "mcs_momentum/D");
    otree->Branch("mcs_pbeta_err", &mcs_pbeta_err, "mcs_pbeta_err/D");
    otree->Branch("best_log_likelihood", &best_log_likelihood, "best_log_likelihood/D");
    otree->Branch("range_pbeta", &range_pbeta, "range_pbeta/D");
    otree->Branch("range_momentum", &range_momentum, "range_momentum/D");

    int tmp_bm_entry = -1;
    int entry_in_daily_file_previous;
    int entry_in_daily_file_next;

    for ( Int_t imcsentry = 0; imcsentry < mcstree->GetEntries(); imcsentry++ ) {
      // MCS file 側の multi を除く
      if ( imcsentry != 0 ) {
	mcstree->GetEntry(imcsentry - 1);
	entry_in_daily_file_previous = entry_in_daily_file;
      }
      if ( imcsentry != mcstree->GetEntries() - 1 ) {
	mcstree->GetEntry(imcsentry + 1);
	entry_in_daily_file_next = entry_in_daily_file;
      }
      
      mcstree->GetEntry(imcsentry);

      if ( imcsentry != 0 &&
	   entry_in_daily_file == entry_in_daily_file_previous ) continue;
      if ( imcsentry != mcstree->GetEntries() - 1 &&
	   entry_in_daily_file == entry_in_daily_file_next ) continue;

      ntbmtree->GetEntry(entry_in_daily_file - 1);

      mcs_pbeta = recon_pbeta;
      mcs_pbeta_err = recon_pbeta_err;
      mcs_momentum = CalculateMomentumFromPBeta(mcs_pbeta, MCS_MUON_MASS);

      BOOST_LOG_TRIVIAL(debug) << "Target entry in daily file : " << entry_in_daily_file;
      if ( ntbm->GetEntryInDailyFile() != entry_in_daily_file ) {
	BOOST_LOG_TRIVIAL(debug) << "Entry in daily file is not the same!";
	BOOST_LOG_TRIVIAL(debug) << "Entry in daily file in MCS file : " << entry_in_daily_file;
	BOOST_LOG_TRIVIAL(debug) << "Entry in daily file in NTBM file : " << ntbm->GetEntryInDailyFile();
      }

      best_log_likelihood = log_likelihood[0][0];
      
      if ( ntbm->GetNumberOfTracks() == 0 ) continue;      

      int n_3d_clusters = 0;
      int cluster_id = -1;

      // NTBM 側の multi を削除
      for ( int icluster = 0; icluster < ntbm->GetNumberOfNinjaClusters(); icluster ++ ) {
	if ( ntbm->GetNumberOfHits(icluster).at(0) != 0 &&
	     ntbm->GetNumberOfHits(icluster).at(1) != 0 ) {
	  n_3d_clusters++;
	  cluster_id = icluster;
	}	  
      }
      if ( n_3d_clusters != 1 ) continue;

      muon_track_id = ntbm->GetBabyMindTrackId(cluster_id);
      if ( ntbm->GetMomentumType(muon_track_id) != 0 ) continue;

      range_momentum = ntbm->GetMomentum(muon_track_id);
      range_pbeta = CalculatePBetaFromMomentum(range_momentum, MCS_MUON_MASS);
      if ( mcs_momentum > range_momentum + 500. &&
	   range_momentum < 1000. ) {
	std::cout << "Entry : " << entry_in_daily_file << std::endl;
	std::cout << "MCS : " << mcs_momentum << ", " << "Range : " << range_momentum << std::endl;
      }
      

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
