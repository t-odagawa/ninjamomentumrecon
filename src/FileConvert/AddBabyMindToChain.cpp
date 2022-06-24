//system includes
#include <string>
#include <fstream>
#include <iostream>
#include <vector>
#include <algorithm>

// boost includes
#include <boost/log/core.hpp>
#include <boost/log/trivial.hpp>
#include <boost/log/expressions.hpp>
#include <boost/filesystem.hpp>

// root includes
#include <TRandom.h>

// B2 includes
#include <B2Reader.hh>
#include <B2SpillSummary.hh>
#include <B2TrackSummary.hh>

// NTBM includes
#include <NTBMConst.hh>
#include <NTBMSummary.hh>

// my includes
#include "McsClass.hpp"
#include "McsFunction.hpp"
#include "MatchData.hpp"
#include "MatchFunction.hpp"

namespace logging = boost::log;
namespace fs = boost::filesystem;

int main ( int argc, char* argv[] ) {

  logging::core::get()->set_filter
    (
     logging::trivial::severity >= logging::trivial::info
     );

  BOOST_LOG_TRIVIAL(info) << "==========Add Baby MIND information to MC chain information start==========";

  if ( argc != 6 ) {
    BOOST_LOG_TRIVIAL(error) << "Usage : " << argv[0]
			     << " <Input B2 file> <Input NTBM file> <Input momch file> <Output momch file> <data dir path>";
  }

  try {

    gRandom->SetSeed(time(NULL));

    B2Reader reader((std::string)argv[1]);
    std::string ntbmfilename = argv[2];
    /*
    TFile *ntbmfile = new TFile(ntbmfilename, "read");
    TTree *ntbmtree = (TTree*)ntbmfile->Get("tree");
    NTBMSummary *ntbm = nullptr;
    ntbmtree->SetBranchAddress("NTBMSummary", &ntbm);
    */
    auto ev_vec = Momentum_recon::ReadEventInformationBin((std::string)argv[3]);

    const std::string data_dir_path = argv[5];
    const MatchData match_data(data_dir_path);
    const MatchFunction match_function(match_data);
    // Shifter efficiency data?

    for ( auto &ev : ev_vec ) {

      if ( ev.chains.empty() ) continue;

      reader.ReadSpill(ev.groupid);
      auto &spill_summary = reader.GetSpillSummary();

      bool fill_flag = true;

      // Baby MIND の飛跡の中で一番長いものを muon とする

      // ECC 中の muon id された飛跡を見つける
      for ( auto &chain : ev.chains ) {

	if ( chain.base.front().pl != 3 &&
	     chain.base.front().pl != 4 ) continue;

	



	

	auto it_recon_track = spill_summary.BeginReconTrack();
	while ( auto *track = it_recon_track.Next() ) {

	  


	  // 対応する Baby MIND range momentum などを追加
	  if ( track->GetIsStopping() ) {
	    chain.stop_flag = 1;
	  }
	  else {
	    chain.stop_flag = 0;
	  }
	  chain.bm_range_mom = track->GetFinalAbsoluteMomentum().GetValue();
	  chain.bm_range_mom_error[0] = track->GetFinalAbsoluteMomentum().GetError();
	  chain.bm_range_mom_error[1] = track->GetFinalAbsoluteMomentum().GetError();
	}
	
	

	
	// Detection/connection efficiency を計算
	// 乱数によってイベントを見つけるかを確認する
	double efficiency = match_function.GetTrackerEfficiency(15.);
	double prob = gRandom->Uniform();
	if ( prob < 1. - efficiency ) fill_flag = false;	

	break;

      }

      if ( !fill_flag ) {
	for ( auto &chain : ev.chains ) {
	  chain.base.clear();
	  chain.base_pair.clear();
	}
	ev.chains.clear();
      }
      
    }

    Momentum_recon::WriteEventInformationBin(argv[2], ev_vec);

    
  } catch ( const std::runtime_error &error ) {
    BOOST_LOG_TRIVIAL(fatal) << "Runtime error : " << error.what();
    std::exit(1);
  }

  BOOST_LOG_TRIVIAL(info) << "==========Finish==========";
  std::exit(0);

}
