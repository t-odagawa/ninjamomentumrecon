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
#include <TMath.h>
#include <TVector3.h>

// B2 includes
#include <B2Const.hh>
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
     // logging::trivial::severity >= logging::trivial::debug
     //logging::trivial::severity >= logging::trivial::trace
     );

  BOOST_LOG_TRIVIAL(info) << "==========Add Baby MIND information to MC chain information start==========";

  if ( argc != 6 ) {
    BOOST_LOG_TRIVIAL(error) << "Usage : " << argv[0]
			     << " <Input B2 file> <Input NTBM file> <Input momch file> <Output momch file> <data dir path>";
    std::exit(1);
  }

  try {

    gRandom->SetSeed(time(NULL));

    B2Reader reader((std::string)argv[1]);

    TFile *ntbmfile = new TFile(argv[2], "read");
    TTree *ntbmtree = (TTree*)ntbmfile->Get("tree");
    NTBMSummary *ntbm = nullptr;
    ntbmtree->SetBranchAddress("NTBMSummary", &ntbm);

    auto ev_vec = Momentum_recon::ReadEventInformationBin((std::string)argv[3]);
    std::vector<Momentum_recon::Event_information> ev_vec_out = {};

    const std::string data_dir_path = argv[5];
    const MatchData match_data(data_dir_path);
    const MatchFunction match_function(match_data);
    // Shifter efficiency data?

    const double bm_nt_distance = (BABYMIND_POS_Z + BM_SECOND_LAYER_POS) 
      - (NINJA_POS_Z + NINJA_ECC_POS_Z + NINJA_DESIC_DEPTH / 2. - NINJA_DESIC_THICK
	 - NINJA_ENV_THICK - 3. * NINJA_FILM_THICK - NINJA_SS_AC_THICK);
    BOOST_LOG_TRIVIAL(debug) << "Baby MIND-tracker distance = " << bm_nt_distance << " mm";


    for ( auto &ev : ev_vec ) {

      if ( ev.chains.empty() ) continue;

      BOOST_LOG_TRIVIAL(debug) << "Entry : " << ev.groupid;
      // if ( ev.groupid > 64 ) break;
      reader.ReadSpill(ev.groupid);
      auto &spill_summary = reader.GetSpillSummary();

      ntbmtree->GetEntry(ev.groupid);

      bool fill_flag = true;

      int bm_stop_flag = -1;
      int bm_charge = 0;
      double range_mom = -1;
      double muon_angle = -1.;
      std::vector<double > muon_position;
      std::vector<double > muon_tangent;
      double track_length = -1.;

      // Baby MIND の飛跡の中で一番長いものを muon とする
      for ( int itrack = 0; itrack < ntbm->GetNumberOfTracks(); itrack++ ) {
	if ( ntbm->GetTrackLengthTotal(itrack) > track_length ) {
	  BOOST_LOG_TRIVIAL(trace) << "Baby MIND muon track update by longer length";
	  bm_stop_flag = ntbm->GetMomentumType(itrack);
	  bm_charge = ntbm->GetCharge(itrack);
	  range_mom = ntbm->GetMomentum(itrack);
	  auto tangent = ntbm->GetBabyMindTangent(itrack);
	  muon_angle = std::atan(std::hypot(tangent.at(0),tangent.at(1))) * TMath::RadToDeg();
	  muon_position = ntbm->GetBabyMindPosition(itrack);
	  muon_position.at(0) += BABYMIND_POS_Y;
	  muon_position.at(1) += BABYMIND_POS_X;
	  muon_tangent = ntbm->GetBabyMindTangent(itrack);
	}
      }

      // ECC 中の muon id された飛跡を見つける
      // pl3 or 4 に basetrack を持つ chain の中から Baby MIND との位置ズレが
      // 最も小さいもの？
      // 角度で適当にカットしないとおかしな chance coincidence がありえる?
      // pl3 or 4 に basetrack が一つもない場合はすでに消去されている
      if ( ntbm->GetNumberOfTracks() < 1 ) {
	fill_flag = false;
      }
      else {
	double dx = 100000;
	double dy = 100000;

	int muon_chain_id = -1;	
	int ichain = 0;
	for ( auto chain : ev.chains ) {

	  ichain++;
	  if ( chain.base.front().pl != 3 &&
	       chain.base.front().pl != 4 ) continue;
	  
	  double x = chain.base.front().x / 1000.;
	  double y = chain.base.front().y / 1000.; // um -> mm
	  x = x - 125. + NINJA_POS_X + NINJA_ECC_POS_X;
	  y = y - 125. + NINJA_POS_Y + NINJA_FV_IRON_POS_Y;
	  double ax = chain.base.front().ax;
	  double ay = chain.base.front().ay;

	  double ex_dx = x + ax * bm_nt_distance - muon_position.at(1);
	  double ex_dy = y + ay * bm_nt_distance - muon_position.at(0);
	  
	  if ( ex_dx * ex_dx + ex_dy * ex_dy < dx * dx + dy * dy ) {
	    BOOST_LOG_TRIVIAL(trace) << "ECC muon track update by smaller position difference";
	    muon_chain_id = ichain - 1;
	  }

	}

	if ( muon_chain_id >= 0 ) {
	
	  auto &muon_chain = ev.chains.at(muon_chain_id);
	  if ( bm_stop_flag == 0 ) {
	    muon_chain.stop_flag = 1;
	  }
	  else {
	    muon_chain.stop_flag = 0;
	  }
	  muon_chain.particle_flag += 13; // true * 10000 + recon
	  muon_chain.charge_sign = bm_charge;

	  int nwater = match_function.GetNumIronPlate(ev.vertex_pl);
	  int niron = match_function.GetNumWaterPlate(ev.vertex_pl);
	  TVector3 track_direction(muon_chain.base.back().ax, muon_chain.base.back().ay, 1.);
	  track_length += track_direction.Mag() * (2.3 * nwater + 0.5 * niron * 8.03);
	  track_length += match_function.GetDWGTrackLength(track_direction);
	  match_function.ConvertFromLengthToMom(range_mom, track_length);

	  muon_chain.bm_range_mom = range_mom;
	  muon_chain.bm_range_mom_error[0] = range_mom * match_function.GetBmErr(range_mom);
	  muon_chain.bm_range_mom_error[1] = range_mom * match_function.GetBmErr(range_mom);
	  // Detection/connection efficiency を計算
	  // weight に efficiency をかける
	  double efficiency = match_function.GetTrackerEfficiency(muon_angle);
	  efficiency *= 0.93; // preliminary shifter efficiency
	  ev.weight *= efficiency;
	}
	else { fill_flag = false; }

      }
      
      if ( !fill_flag ) {
	for ( auto &chain : ev.chains ) {
	  chain.base.clear();
	  chain.base.shrink_to_fit();
	  chain.base_pair.clear();
	  chain.base_pair.shrink_to_fit();
	}
	ev.chains.clear();
	ev.chains.shrink_to_fit();
      }

      ev_vec_out.push_back(ev);

    }

    Momentum_recon::WriteEventInformationBin(argv[4], ev_vec_out);

  } catch ( const std::runtime_error &error ) {
    BOOST_LOG_TRIVIAL(fatal) << "Runtime error : " << error.what();
    std::exit(1);
  }

  BOOST_LOG_TRIVIAL(info) << "==========Finish==========";
  std::exit(0);

}
