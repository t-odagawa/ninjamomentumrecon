// boost includes
#include <boost/log/core.hpp>
#include <boost/log/trivial.hpp>
#include <boost/log/expressions.hpp>

// B2 includes
#include <B2Reader.hh>
#include <B2Writer.hh>
#include <B2SpillSummary.hh>
#include <B2EventSummary.hh>
#include <B2VertexSummary.hh>
#include <B2EmulsionSummary.hh>

// ROOT includes

// system includes
#include <string>
#include <fstream>
#include <iostream> 
#include <vector>
#include <algorithm>
#include <map>
#include <utility>

// my includes
#include "McsClass.hpp"
#include "McsFunction.hpp"
#include "ConnectionData.hpp"
#include "ConnectionClass.hpp"
#include "ConnectionFunction.hpp"

namespace logging = boost::log;
namespace fs = boost::filesystem;

bool EmulsionCompareDownToUp(const B2EmulsionSummary *lhs,
			     const B2EmulsionSummary *rhs) {
  if ( lhs->GetEcc() != rhs->GetEcc() ) 
    return lhs->GetEcc() < rhs->GetEcc();
  else if ( lhs->GetParentTrackId() != rhs->GetParentTrackId() )
    return lhs->GetParentTrackId() < rhs->GetParentTrackId();
  else
    return lhs->GetPlate() < rhs->GetPlate();
}


int main ( int argc, char* argv[] ) {
  
  logging::core::get()->set_filter
    (
     logging::trivial::severity >= logging::trivial::debug
     // logging::trivial::severity >= logging::trivial::trace
     );

  BOOST_LOG_TRIVIAL(info) << "==========MC output to chain like data conversion start==========";

  if ( argc != 5 ) {
    BOOST_LOG_TRIVIAL(error) << "Usage : " << argv[0]
			     << " <Input B2 file> <Output event info file> <ECC id> <data dir path>";
    std::exit(1);
  }

  try {

    B2Reader reader((std::string)argv[1]);
    std::string ofilename = argv[2];

    std::vector<Momentum_recon::Event_information > ev_vec = {};

    const Int_t ecc_id = std::atoi(argv[3]);


    const std::string data_dir_path = argv[4];
    const ConnectionData connection_data(data_dir_path);
    const ConnectionFunction connection_function(connection_data);

    while ( reader.ReadNextSpill() > 0 ) {

      if ( reader.GetEntryNumber() != 100 ) continue;
      
      Momentum_recon::Event_information ev;

      auto &spill_summary = reader.GetSpillSummary();

      // Event information
      ev.groupid = reader.GetEntryNumber();
      ev.entry_in_daily_file = reader.GetEntryNumber();
      ev.ecc_id = ecc_id;

      // MC weight & beam info
      auto it_event = spill_summary.BeginTrueEvent();
      const auto *event = it_event.Next();
      auto &primary_vertex_summary = event->GetPrimaryVertex();
      ev.weight = primary_vertex_summary.GetMcWeight();
      ev.nu_energy = event->GetPrimaryParticleEnergy().GetValue();
      TVector3 nu_mom = event->GetPrimaryParticleMomentum().GetValue();
      ev.nu_ax = nu_mom.X() / nu_mom.Z();
      ev.nu_ay = nu_mom.Y() / nu_mom.Z();
      TVector3 vertex_position = primary_vertex_summary.GetAbsolutePosition().GetValue();
      connection_function.CalcPosInEccCoordinate(vertex_position, ecc_id);
      ev.true_vertex_position[0] = vertex_position.X();
      ev.true_vertex_position[1] = vertex_position.Y();
      ev.true_vertex_position[2] = vertex_position.Z();

      // Get true emulsion tracks
      std::vector<const B2EmulsionSummary* > emulsions = {};
      connection_function.GetTrueEmulsionTracks(emulsions, spill_summary, ecc_id);

      if ( emulsions.empty() ) continue;
      std::sort(emulsions.begin(), emulsions.end(), EmulsionCompareDownToUp);

      std::cout << "Entry : " << reader.GetEntryNumber() << std::endl;

      // Get true chains
      std::vector<std::vector<const B2EmulsionSummary* > > chains = {};
      connection_function.GetTrueEmulsionChains(chains, emulsions);

      // Add true information
      connection_function.AddTrueChainsToEventInfo(ev, chains, ecc_id);

      // Smear
      std::vector<B2EmulsionSummary* > emulsions_smeared;
      connection_function.SmearEmulsions(emulsions_smeared, emulsions);
      

      // Apply detection efficiency
      std::vector<B2EmulsionSummary* > emulsions_detected;
      connection_function.ApplyDetectionEfficiency(emulsions_detected, emulsions_smeared, ecc_id);

      // Fiducial volume cut
      std::vector<B2EmulsionSummary* > emulsions_detected_in_fv;
      connection_function.ApplyFVCut(emulsions_detected_in_fv, emulsions_detected, ecc_id);     
      // connection_data.DrawFiducialAreaData(ecc_id, 12, emulsions_detected, reader.GetEntryNumber());
      
      // Linklet
      // Black についてはVPH <-> PID が同じものだけを対象に loose につなぐ
      // Multi 消しはしなくてよい
      std::vector<std::pair<B2EmulsionSummary*, B2EmulsionSummary* > > linklets;
      connection_function.GenerateLinklet(linklets, emulsions_detected_in_fv);

      // Group reconstruction
      std::vector<std::pair<B2EmulsionSummary*, std::vector<std::pair<B2EmulsionSummary*, B2EmulsionSummary* > > > > groups;
      connection_function.GenerateGroup(groups, linklets, emulsions_detected_in_fv);

      // 連結成分を取り出す -> それ以上上流に行けないtrack, 下流に行けないtrack を引っ張ってくる
      // 引っ張ってきたtrack をつかって再接続
      std::vector<std::pair<B2EmulsionSummary*, std::vector<std::pair<B2EmulsionSummary*, B2EmulsionSummary* > > > > groups_reconnected;
      connection_function.ReconnectGroups(groups_reconnected, groups, emulsions_detected_in_fv);
      for ( auto group : groups_reconnected ) {
	std::cout << "start_seg : PL" << group.first->GetPlate() + 1 << ", " << group.first->GetEmulsionTrackId() << std::endl;
	for ( auto linklet : group.second ) {
	  std::cout << "(PL" << linklet.first->GetPlate() + 1 << ", " << linklet.first->GetEmulsionTrackId() << ", "
		    << "PL" << linklet.second->GetPlate() + 1 << ", " << linklet.second->GetEmulsionTrackId() << ")" << std::endl; 
	}
      }

      // group をほどくアルゴリズムに突っ込む (if necessary)

      // vertex plate を確定させる
      std::pair<B2EmulsionSummary*, std::vector<std::pair<B2EmulsionSummary*, B2EmulsionSummary*> > > muon_group;
      /*
      if ( connection_function.SelectMuonGroup(groups_reconnected, muon_group) ) {

	// vertex plate に attach する basetrack を探す
	std::vector<B2EmulsionSummary* > emulsions_partner;
	//connection_function.PartnerSearch(vertex_track, emulsions_partner, emulsions);

	// vertex plate に attach する basetrack に対して group を作る
	// 4pl さきまで見るかどうかが先程との違い
	// connection_function.GenerateGroupPartner(groups_partner, linklets, emulsions_partner);

	// attach した group の再接続を行う?
	// connection_function.ReconnectGroupsPartner(groups_reconnected_partner, groups_partner);

	// group を event information に変換
	// connection_function.AddGroupsToEventInfo(ev, group);
	
      }
*/
      ev_vec.push_back(ev);

    }

    Momentum_recon::WriteEventInformationBin(ofilename, ev_vec);

  } catch (const std::runtime_error &error) {
    BOOST_LOG_TRIVIAL(fatal) << "Runtime error : " << error.what();
    std::exit(1);
  } catch (const std::invalid_argument &error) {
    BOOST_LOG_TRIVIAL(fatal) << "Invalid argument error : " << error.what();
    std::exit(1);
  }

  BOOST_LOG_TRIVIAL(info) << "==========MC output to chain like data conversion finish==========";
  std::exit(0);

}
