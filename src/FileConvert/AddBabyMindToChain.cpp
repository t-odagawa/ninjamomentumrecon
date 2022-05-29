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

// B2 includes
#include <B2Reader.hh>
#include <B2SpillSummary.hh>

// NTBM includes
#include <NTBMConst.hh>
#include <NTBMSummary.hh>

// my includes
#include "McsClass.hpp"
#include "McsFunction.hpp"

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

    B2Reader reader((std::string)argv[1]);
    std::string ntbmfilename = argv[2];
    TFile *ntbmfile = new TFile(ntbmfilename, "read");
    TTree *ntbmtree = (TTree*)ntbmfile->Get("tree");
    NTBMSummary *ntbm = nullptr;
    // ntbmtree->SetBranchAddress("NTBMSummary", &ntbm);

    std::vector<Momentum_recon::Event_information > ev_vec = Momentum_recon::ReadEventInformationBin((std::string)argv[3]);
    std::vector<Momentum_recon::Event_information > ev_vec_out = {};
    ev_vec_out.reserve(ev_vec.size());

    const std::string data_dir_path = argv[5];
    // Shifter efficiency data?

    for ( auto ev : ev_vec ) {

      if ( ev.chains.empty() ) continue;

      reader.ReadSpill(ev.groupid);
      auto &spill_summary = reader.GetSpillSummary();

      // ECC 中の muon id された飛跡を見つける
      for ( auto &chain : ev.chains ) {

	if ( chain.particle_flag % 10000 != 13 ) continue;



	// 対応するトラックが Baby MIND にあるかを確認する

	

	const auto it_recon_track = spill_summary.BeginReconTrack();
	while ( auto *track != it_recon_track.Next() ) {

	  


	  // 対応する Baby MIND range momentum などを追加
	  chain.bm_range_mom = track->GetFinalAbsoluteMomentum().GetValue();
	  chain.bm_range_mom_err[0] = track->GetFinalAbsoluteMomentum().GetError();
	  chain.bm_range_mom_err[1] = track->GetFinalAbsoluteMomentum().GetError();
	}
	
	

	
	// Detection/connection efficiency を計算
	// 乱数によってイベントを見つけるかを確認する
	double efficiency = GetEfficiency();
	double prob = gRandom->Uniform();
	if ( prob < efficiency ) fill_flag = false;	

      }

      if ( fill_flag ) ev_vec_out.push_back(ev);

    }

    Momentum_recon::WriteEventInformationBin((std::string)argv[4], ev_vec_out);
    
  } catch ( const std::runtime_error &error ) {
    BOOST_LOG_TRIVIAL(fatal) << "Runtime error : " << error.what();
    std::exit(1);
  }

  BOOST_LOG_TRIVIAL(info) << "==========Finish==========";
  std::exit(0);

}
