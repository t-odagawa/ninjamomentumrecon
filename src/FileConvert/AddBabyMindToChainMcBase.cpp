// system includes
#include <string>
#include <fstream>
#include <iostream>
#include <vector>
#include <map>
#include <algorithm>

// boost includes
#include <boost/log/core.hpp>
#include <boost/log/trivial.hpp>
#include <boost/log/expressions.hpp>

// root includes
#include <TFile.h>
#include <TTree.h>
#include <TMath.h>
#include <TVector3.h>

// B2 includes
#include <B2Const.hh>
#include <B2Enum.hh>
#include <B2Reader.hh>
#include <B2SpillSummary.hh>
#include <B2TrackSummary.hh>
#include <B2EmulsionSummary.hh>

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

int main (int argc, char* argv[]) {

  logging::core::get()->set_filter
    (
     logging::trivial::severity >= logging::trivial::info
     // logging::trivial::severity >= logging::trivial::debug
     );

  BOOST_LOG_TRIVIAL(info) << "==========Connect BM tracks to MC chains==========";

  if ( argc != 6 ) {
    BOOST_LOG_TRIVIAL(error) << "Usage : " << argv[0]
			     << " <Input B2 file> <Input NTBM file> <Input momch file> <Outout momch file> <data dir path>";
    std::exit(1);
  }

  try {

    B2Reader reader(argv[1]);

    TFile *ntbmfile = new TFile(argv[2], "read");
    auto ntbmtree = (TTree*)ntbmfile->Get("tree");
    NTBMSummary *ntbm = nullptr;
    ntbmtree->SetBranchAddress("NTBMSummary", &ntbm);

    auto ev_vec = Momentum_recon::ReadEventInformationBin(static_cast<std::string>(argv[3]));
    std::vector<Momentum_recon::Event_information > ev_vec_out = {};
    ev_vec_out.reserve(ev_vec.size());

    const std::string data_dir_path = argv[5];
    const MatchData match_data(data_dir_path);
    const MatchFunction match_function(match_data);

    for ( auto &ev : ev_vec ) {

      if ( ev.chains.empty() ) continue;      

      bool fill_flag = true;

      // BOOST_LOG_TRIVIAL(info) << "Entry : " << ev.groupid;
      BOOST_LOG_TRIVIAL(debug) << "Entry : " << ev.groupid;

      reader.ReadSpill(ev.groupid);
      auto &spill_summary = reader.GetSpillSummary();
      std::map<int, const B2EmulsionSummary* > tss_emulsions; // chain id, emulsion summary
      auto it_emulsion = spill_summary.BeginEmulsion();
      while ( const auto emulsion = it_emulsion.Next() ) {
	if ( emulsion->GetParentTrackId() == 0 ) continue;
	if ( emulsion->GetFilmType() != B2EmulsionType::kShifter ) continue;
	if ( emulsion->GetPlate() != 14 ) continue;
	tss_emulsions.insert(std::make_pair(emulsion->GetParentTrackId(), emulsion));
      }

      ntbmtree->GetEntry(ev.groupid);

      // 2d cluster を持ってくる
      std::vector<int > bm_2d_cluster_vec = {};
      for ( int icluster = 0; icluster < ntbm->GetNumberOfNinjaClusters(); icluster++ ) {
	if ( ntbm->GetBabyMindTrackId(icluster) < 0 ) continue;
	auto position = ntbm->GetNinjaPosition(icluster);

	// position が ECC5 の範囲外であれば matching はされない
	if ( match_function.IsMatchBMNTCandidate(5, position) ) {
	  BOOST_LOG_TRIVIAL(debug) << "X : " << position.at(B2View::kTopView) << ", "
				   << "Y : " << position.at(B2View::kSideView);
	  bm_2d_cluster_vec.push_back(icluster);
	}

      }

      if ( bm_2d_cluster_vec.empty() ) {
	for ( auto &chain : ev.chains ) {
	  chain.base.clear();
	  chain.base.shrink_to_fit();
	  chain.base_pair.clear();
	  chain.base_pair.shrink_to_fit();
	}
	ev.chains.clear();
	ev.chains.shrink_to_fit();
	ev_vec_out.push_back(ev);
	continue;
      }

      // Tracker SS までつながっている basetrack chain id を mom chain から持ってくる
      std::vector<std::pair<int, const B2EmulsionSummary* > > chain_tss_emulsions; // chain id, emulsion summary
      for ( auto chain : ev.chains ) {
	if ( chain.base.front().pl != 3 &&
	     chain.base.front().pl != 4 ) continue;
	
	int true_track_id = chain.chainid;
	if ( tss_emulsions.find(true_track_id) != tss_emulsions.end() ) {
	  const auto emulsion = tss_emulsions.at(true_track_id);
	  chain_tss_emulsions.push_back(std::make_pair(true_track_id, emulsion));
	}
      }

      // 各 2d cluster に対して chi2 が最も小さいもの，かつ，ある値より小さい場合に matching する
      std::multimap<int, std::pair<int, double > > match_bt_nt_map; // chain id, cluster id, chi2
      for ( auto cluster_id : bm_2d_cluster_vec ) {

	bool match_flag = false;
	int chain_tss_id = -1;
	double chi2_max = match_function.GetChisCutValue(ntbm, cluster_id); // chi2 がある値以上は matching とは判断しない
	double chi2 = chi2_max;
	for ( auto chain_tss : chain_tss_emulsions ) {
	  const auto emulsion = chain_tss.second;
	  double chi2_tmp = match_function.CalculateShifterTrackerChi2(ntbm, cluster_id, emulsion);
	  if ( chi2_tmp < chi2 ) {
	    BOOST_LOG_TRIVIAL(debug) << "Chi2 update : " << chi2 << " -> " << chi2_tmp;
	    chi2 = chi2_tmp;
	    chain_tss_id = chain_tss.first;
	  }
	}

	if ( chi2 < chi2_max ) match_flag = true;

	if ( match_flag ) {
	  BOOST_LOG_TRIVIAL(debug) << "Cluster : " << cluster_id << " and "
				   << "Basetrack : " << chain_tss_id <<  " are matched";
	  auto pair = std::make_pair(cluster_id, chi2);
	  match_bt_nt_map.insert(std::make_pair(chain_tss_id, pair));
	}
	
      }

      // 各 basetrack に対して複数の 2d cluster が接続されていた場合は最良のものだけを残す
      std::map<int, int> match_bt_nt_unique_map; // chain id, baby mind id
      for ( auto itr = match_bt_nt_map.begin();
	    itr != match_bt_nt_map.end();
	    itr++ ) {

	auto range = match_bt_nt_map.equal_range(itr->first);

	auto pair = range.first->second;
	int bm_track_id = ntbm->GetBabyMindTrackId(pair.first);
	double chi2 = pair.second;

	for ( auto itr1 = range.first; itr1 != range.second; itr1++ ) {
	  auto pair = itr1->second;
	  double bm_track_id_tmp = ntbm->GetBabyMindTrackId(pair.first);
	  double chi2_tmp = pair.second;
	  if ( chi2_tmp < chi2 ) {
	    bm_track_id = bm_track_id_tmp;
	    chi2 = chi2_tmp;
	  }
	}

	match_bt_nt_unique_map.insert(std::make_pair(itr->first, bm_track_id));

      }

      // BM track length が最も長いものを muon として判断する
      int muon_chain_id = -1;
      int muon_track_id = -1;
      double track_length = -1.;
      for ( auto itr = match_bt_nt_unique_map.begin();
	    itr != match_bt_nt_unique_map.end();
	    itr++ ) {
	int track_id_tmp = itr->second;

	if ( track_id_tmp < 0 ) continue;

	double track_length_tmp = ntbm->GetTrackLengthTotal(track_id_tmp);

	if ( track_length_tmp > track_length ) {
	  muon_chain_id = itr->first;
	  muon_track_id = track_id_tmp;
	  track_length = track_length_tmp;
	}

      }

      BOOST_LOG_TRIVIAL(debug) << "Chain id : " << muon_chain_id << ", BM id : " << muon_track_id;

      if ( muon_chain_id < 0 || muon_track_id < 0 ) {
	for ( auto &chain : ev.chains ) {
	  chain.base.clear();
	  chain.base.shrink_to_fit();
	  chain.base_pair.clear();
	  chain.base_pair.shrink_to_fit();
	}
	ev.chains.clear();
	ev.chains.shrink_to_fit();
	ev_vec_out.push_back(ev);
	continue;
      }

      int id = -1;
      for ( auto chain : ev.chains ) {
	id++;
	if ( chain.chainid == muon_chain_id ) break;
      }

      auto &muon_chain = ev.chains.at(id);
      
      int bm_stop_flag = ntbm->GetMomentumType(muon_track_id);

      if ( bm_stop_flag == 0 ) {
	muon_chain.stop_flag = 1;
      }
      else {
	muon_chain.stop_flag = 0;
      }
      muon_chain.particle_flag += 13;
      muon_chain.charge_sign = ntbm->GetCharge(muon_track_id);

      double muon_track_length = ntbm->GetTrackLengthTotal(muon_track_id);
      double range_mom = -1.;
      int nwater = match_function.GetNumWaterPlate(ev.vertex_pl % 1000);
      int niron = match_function.GetNumIronPlate(ev.vertex_pl % 1000);
      TVector3 track_direction(muon_chain.base.front().ax, muon_chain.base.front().ay, 1.);
      muon_track_length += track_direction.Mag() * (0.23 * nwater * 1.289 + 0.05 * niron * 8.03);
      muon_track_length += match_function.GetDWGTrackLength(track_direction);
      match_function.ConvertFromLengthToMom(range_mom, muon_track_length);

      muon_chain.bm_range_mom = range_mom;
      muon_chain.bm_range_mom_error[0] = range_mom * match_function.GetBmErr(range_mom);
      muon_chain.bm_range_mom_error[1] = range_mom * match_function.GetBmErr(range_mom);
      
      // ECC-shifter の efficiency は weight として実装
      ev.weight *= 0.99;
     
      ev_vec_out.push_back(ev);
 
    }

    Momentum_recon::WriteEventInformationBin(argv[4], ev_vec_out);


  } catch (const std::runtime_error &error) {
    BOOST_LOG_TRIVIAL(fatal) << "Runtime error : " << error.what();
    std::exit(1);
  }

  std::exit(0);

}
