// boost includes
#include <boost/log/core.hpp>
#include <boost/log/trivial.hpp>
#include <boost/log/expressions.hpp>

// B2 includes
#include <B2Reader.hh>
#include <B2Writer.hh>
#include <B2SpillSummary.hh>
#include <B2EventSummary.hh>
#include <B2TrackSummary.hh>
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
     logging::trivial::severity >= logging::trivial::info
     // logging::trivial::severity >= logging::trivial::debug
     // logging::trivial::severity >= logging::trivial::trace
     );

  BOOST_LOG_TRIVIAL(info) << "==========MC output to chain like data conversion start==========";

  if ( argc != 6 ) {
    BOOST_LOG_TRIVIAL(error) << "Usage : " << argv[0]
			     << " <Input B2 file> <Output event info file> <ECC id> <material id (0:water/2:iron)> <data dir path>";
    std::exit(1);
  }

  try {

    B2Reader reader((std::string)argv[1]);
    std::string ofilename = argv[2];

    std::vector<Momentum_recon::Event_information > ev_vec = {};

    const Int_t ecc_id = std::atoi(argv[3]);

    const Int_t material_id = std::atoi(argv[4]);
    if ( material_id != B2Material::kWater &&
	 material_id != B2Material::kIron )
      throw std::invalid_argument("Water or iron is only accepted material");

    const std::string data_dir_path = argv[5];
    const ConnectionData connection_data(data_dir_path);
    const ConnectionFunction connection_function(connection_data);

    while ( reader.ReadNextSpill() > 0 ) {

      // if ( reader.GetEntryNumber() > 1000 ) continue;
      // if ( reader.GetEntryNumber() != 1489 ) continue;

      Momentum_recon::Event_information ev;

      auto &spill_summary = reader.GetSpillSummary();

      // Event information
      ev.groupid = reader.GetEntryNumber();
      ev.entry_in_daily_file = reader.GetEntryNumber();

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
      ev.vertex_pl = (primary_vertex_summary.GetPlane() % 1000 + 1) * 1000;
      ev.ecc_id = (primary_vertex_summary.GetPlane() / 1000 + 1) * 10;

      // Get true emulsion tracks
      std::vector<const B2EmulsionSummary* > emulsions = {};
      connection_function.GetTrueEmulsionTracks(emulsions, spill_summary, ecc_id);

      if ( emulsions.empty() ) continue;
      std::sort(emulsions.begin(), emulsions.end(), EmulsionCompareDownToUp);
      // std::cout << "Entry:" << reader.GetEntryNumber() << std::endl;
      // Get true chains
      std::vector<std::vector<const B2EmulsionSummary* > > chains = {};
      connection_function.GetTrueEmulsionChains(chains, emulsions);

      // Add true information
      connection_function.AddTrueChainsToEventInfo(ev, chains, ecc_id, material_id);

      // Smear
      std::vector<B2EmulsionSummary* > emulsions_smeared;
      connection_function.SmearEmulsions(emulsions_smeared, emulsions);
      

      // Apply detection efficiency
      std::vector<B2EmulsionSummary* > emulsions_detected;
      connection_function.ApplyDetectionEfficiency(emulsions_detected, emulsions_smeared, ecc_id);

      // Fiducial volume cut
      std::vector<B2EmulsionSummary* > emulsions_detected_in_fv;
      connection_function.ApplyFVCut(emulsions_detected_in_fv, emulsions_detected, ecc_id);     
      // if ( reader.GetEntryNumber() == 22 )
      //   connection_data.DrawFiducialAreaData(ecc_id, 37, emulsions_detected, reader.GetEntryNumber());
      
      // Linklet
      // Black については VPH <-> PID が同じものだけを対象に loose につなぐ
      // Multi 消しは background に効くだけなはずなので省略する
      std::vector<std::pair<B2EmulsionSummary*, B2EmulsionSummary* > > linklets;
      connection_function.GenerateLinklet(linklets, emulsions_detected_in_fv);

      // Group reconstruction
      std::vector<std::pair<B2EmulsionSummary*, std::vector<std::pair<B2EmulsionSummary*, B2EmulsionSummary* > > > > groups;
      connection_function.GenerateGroup(groups, linklets, emulsions_detected_in_fv);

      // 連結成分を取り出す -> それ以上上流に行けない track, 下流に行けない track を引っ張ってくる
      // 引っ張ってきた track をつかって再接続
      std::vector<std::pair<B2EmulsionSummary*, std::vector<std::pair<B2EmulsionSummary*, B2EmulsionSummary* > > > > groups_reconnected;
      connection_function.ReconnectGroups(groups_reconnected, groups, emulsions_detected_in_fv);

      // group をほどくアルゴリズムに突っ込む部分は S/N を良くするために行っているが
      // event-base の MC では省略する
      std::vector<std::pair<B2EmulsionSummary*, std::vector<std::pair<B2EmulsionSummary*, B2EmulsionSummary* > > > > groups_modified;
      connection_function.ModifyGroups(groups_modified, groups_reconnected, emulsions_detected_in_fv);
      /*
      for ( auto group : groups_modified ) {
	std::cout << "start track : PL" << group.first->GetPlate() + 1 << ", " << group.first->GetEmulsionTrackId() << std::endl;
	for ( auto linklet : group.second ) {
	  std::cout << "( PL" << linklet.first->GetPlate() + 1  << ", " << linklet.first->GetEmulsionTrackId() << ", "
		    << "  PL" << linklet.second->GetPlate() + 1 << ", " << linklet.second->GetEmulsionTrackId() << ")" << std::endl;
	}
      }
      */
      // vertex plate を確定させる
      std::pair<B2EmulsionSummary*, std::vector<std::pair<B2EmulsionSummary*, B2EmulsionSummary*> > > muon_group;
      std::vector<std::pair<B2EmulsionSummary*, std::vector<std::pair<B2EmulsionSummary*, B2EmulsionSummary* > > > > groups_partner;
      B2EmulsionSummary* vertex_track;

      if ( connection_function.SelectMuonGroup(groups_modified, muon_group, vertex_track, emulsions_detected_in_fv) ) {

	int stop_flag = -1;

	if ( vertex_track->GetPlate() + 1 > 131 ) {
	  stop_flag = 0;
	  ev.vertex_material = -2;
	} // penetrate
	else if ( !connection_function.JudgeEdgeOut(vertex_track, ecc_id) ) { // stop check
	  stop_flag = 1; // ECC stop
	  // vertex plate に attach する basetrack を探す
	  std::vector<B2EmulsionSummary* > emulsions_partner;
	  TVector3 recon_vertex(0., 0., 0.);
	  connection_function.PartnerSearch(vertex_track, groups_modified,
	  				    emulsions_partner, emulsions_detected_in_fv,
	  				    ecc_id,
	  				    recon_vertex);

	  ev.recon_vertex_position[0] = recon_vertex.X();
	  ev.recon_vertex_position[1] = recon_vertex.Y();
	  ev.recon_vertex_position[2] = recon_vertex.Z();

	  ev.vertex_material = connection_function.GetVertexMaterial(recon_vertex);

	  // partner の group は chain の再接続の linklet を加えて 4pl 先まで行っているが
	  // background に効くだけなはずなので省略する
	  
	  // vertex plate に attach する basetrack の group をもってくる
	  connection_function.SelectPartnerGroups(groups_modified, groups_partner, emulsions_partner, emulsions_detected_in_fv);

	}
	else { // side escape
	  stop_flag = 2;
	  ev.vertex_material = -3;
	}

	// group を event information に変換
	connection_function.AddGroupsToEventInfo(ev, vertex_track, muon_group, groups_partner, emulsions_detected_in_fv, ecc_id);
	ev.vertex_pl += vertex_track->GetPlate() + 1;
	ev.ecc_id += ecc_id + 1;
      }

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
