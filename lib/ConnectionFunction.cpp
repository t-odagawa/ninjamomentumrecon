#include <boost/log/trivial.hpp>
#include <boost/log/core.hpp>
#include <boost/log/expressions.hpp>
#include <boost/unordered_set.hpp>
#include <boost/unordered_map.hpp>

#include <vector>
#include <set>
#include <unordered_set>
#include <map>
#include <unordered_map>

#include <TRandom.h>
#include <TVector3.h>

#include <B2Pdg.hh>
#include <B2Enum.hh>
#include <B2Const.hh>
#include <B2TrackSummary.hh>
#include <B2EmulsionSummary.hh>

#include "McsClass.hpp"
#include "ConnectionFunction.hpp"
#include "ConnectionClass.hpp"
#include "ConnectionData.hpp"

bool CompareParticlePdg(std::vector<const B2EmulsionSummary* > lhs,
					    std::vector<const B2EmulsionSummary* > rhs) {
  if ( std::abs(lhs.at(0)->GetParentTrack().GetParticlePdg()) != std::abs(rhs.at(0)->GetParentTrack().GetParticlePdg()))
    return std::abs(lhs.at(0)->GetParentTrack().GetParticlePdg()) < std::abs(rhs.at(0)->GetParentTrack().GetParticlePdg());
  return false;
}

ConnectionFunction::ConnectionFunction(const ConnectionData &connection_data) : connection_data_(connection_data) {

  // Get allowances
  connection_data_.GetETECCConnectData(et_ecc_param_);
  connection_data_.GetETECCFeConnectData(et_ecc_fe_param_);
  connection_data_.GetETECCFeFeConnectData(et_ecc_fe_fe_param_);
  connection_data_.GetBlackFeConnectData(black_fe_param_);
  connection_data_.GetBlackReFeConnectData(black_re_fe_param_);
  connection_data_.GetBlackReFeFeConnectData(black_re_fe_fe_param_);
  connection_data_.GetBlackReFeWaterConnectData(black_re_fe_water_param_);
  connection_data_.GetBlackWaterConnectData(black_water_param_);
  connection_data_.GetFeConnectData(fe_param_);
  connection_data_.GetFeFeConnectData(fe_fe_param_);
  connection_data_.GetFeFeFeConnectData(fe_fe_fe_param_);
  connection_data_.GetFeWaterConnectData(fe_water_param_);
  connection_data_.GetFeWaterFeConnectData(fe_water_fe_param_);
  connection_data_.GetReETECCFeFeFeConnectData(re_et_ecc_fe_fe_fe_param_);
  connection_data_.GetReETECCFeFeFeFeConnectData(re_et_ecc_fe_fe_fe_fe_param_);
  connection_data_.GetReETECCFeFeFeWaterConnectData(re_et_ecc_fe_fe_fe_water_param_);
  connection_data_.GetReETECCFeFeWaterConnectData(re_et_ecc_fe_fe_water_param_);
  connection_data_.GetReETECCFeWaterConnectData(re_et_ecc_fe_water_param_);
  connection_data_.GetReFeConnectData(re_fe_param_);
  connection_data_.GetReFeFeConnectData(re_fe_fe_param_);
  connection_data_.GetReFeFeFeConnectData(re_fe_fe_fe_param_);
  connection_data_.GetReFeFeFeFeConnectData(re_fe_fe_fe_fe_param_);
  connection_data_.GetReFeFeFeFeFeConnectData(re_fe_fe_fe_fe_fe_param_);
  connection_data_.GetReFeWaterConnectData(re_fe_water_param_);
  connection_data_.GetReFeWaterFeConnectData(re_fe_water_fe_param_);
  connection_data_.GetReFeWaterFeWaterConnectData(re_fe_water_fe_water_param_);
  connection_data_.GetReFeWaterFeWaterFeConnectData(re_fe_water_fe_water_fe_param_);
  connection_data_.GetReWaterFeWaterConnectData(re_water_fe_water_param_);
  connection_data_.GetWaterConnectData(water_param_);
  
  for ( int ecc = 0; ecc < 9; ecc++ ) {
    connection_data_.GetFiducialAreaData(ecc, ecc_fiducial_[ecc]);
  }

  BOOST_LOG_TRIVIAL(info) << "Connection functions are initilized";

}

void ConnectionFunction::GetTrueEmulsionTracks(std::vector<const B2EmulsionSummary* > &emulsions,
					       B2SpillSummary &spill_summary, int ecc_id) const {
  auto it_emulsion = spill_summary.BeginEmulsion();
  while ( const auto emulsion = it_emulsion.Next() ) {
    if ( emulsion->GetParentTrackId() == 0 ) continue;
    if ( emulsion->GetFilmType() != B2EmulsionType::kECC ) continue;
    if ( emulsion->GetEcc() != ecc_id ) continue;
    if ( emulsion->GetPlate() < 2 ) continue;
    emulsions.push_back(emulsion);
  }
  return;
}

void ConnectionFunction::GetTrueEmulsionChains(std::vector<std::vector<const B2EmulsionSummary* > > &chains,
					       std::vector<const B2EmulsionSummary* > &emulsions) const {
  if ( emulsions.empty() )
    throw std::invalid_argument("emulsions should have non zero size : " + emulsions.size());

  int track_id_tmp = emulsions.front()->GetParentTrackId();
  std::vector<const B2EmulsionSummary* > chain;
  chain.reserve(emulsions.size());
  for ( auto & emulsion : emulsions ) {
    if ( emulsion->GetParentTrackId() == track_id_tmp )
      chain.push_back(emulsion);
    else {
      chains.push_back(chain);
      chain.clear();
      track_id_tmp = emulsion->GetParentTrackId();
      chain.push_back(emulsion);
    } 
  }
  chains.push_back(chain);
  std::sort(chains.begin(), chains.end(), CompareParticlePdg);
  return;

}

void ConnectionFunction::AddTrueChainsToEventInfo(Momentum_recon::Event_information &ev,
						  std::vector<std::vector<const B2EmulsionSummary* > > &chains,
						  int ecc_id) const {

  Momentum_recon::Mom_chain mom_chain;
  Momentum_recon::Mom_basetrack mom_basetrack;
  std::pair<Momentum_recon::Mom_basetrack, Momentum_recon::Mom_basetrack> mom_basetrack_pair;

  for ( auto chain : chains ) {
    int num_base = chain.size();
    int num_link = num_base - 1;
    mom_chain.base.clear();
    mom_chain.base_pair.clear();
    mom_chain.base.reserve(num_base);
    mom_chain.base_pair.reserve(num_link);

    mom_chain.chainid = chain.front()->GetParentTrackId();
    mom_chain.stop_flag = 0;
    mom_chain.particle_flag = chain.front()->GetParentTrack().GetParticlePdg();
    if ( chain.front()->GetTangent().GetValue().Z() > 0 )
      mom_chain.direction = 1;
    else
      mom_chain.direction = -1;
    mom_chain.charge_sign = 0;
    mom_chain.bm_range_mom = chain.front()->GetParentTrack().GetInitialAbsoluteMomentum().GetValue();

    double downstream_position_z = 0.;

    for ( auto emulsion : chain ) {

      TVector3 position = emulsion->GetAbsolutePosition().GetValue();
      CalcPosInEccCoordinate(position, ecc_id);
      TVector3 tangent = emulsion->GetTangent().GetValue();
      tangent = (1./tangent.Z()) * tangent;

      mom_basetrack.pl = emulsion->GetPlate() + 1;
      mom_basetrack.rawid = emulsion->GetEmulsionTrackId();
      mom_basetrack.x = position.X();
      mom_basetrack.y = position.Y();
      mom_basetrack.z = position.Z();
      mom_basetrack.ax = tangent.X();
      mom_basetrack.ay = tangent.Y();
      mom_basetrack.m[0].zone = 0;
      mom_basetrack.m[0].view = 0;
      mom_basetrack.m[0].imager = 0;
      mom_basetrack.m[0].ph = std::min((int)((emulsion->GetEdepSum() + emulsion->GetEdepDiff()) * 1000 / 2), 9999);
      mom_basetrack.m[0].pixelnum = 0;
      mom_basetrack.m[0].hitnum = 0;
      mom_basetrack.m[1].zone = 0;
      mom_basetrack.m[1].view = 0;
      mom_basetrack.m[1].imager = 0;
      mom_basetrack.m[1].ph = std::min((int)((emulsion->GetEdepSum() - emulsion->GetEdepDiff()) * 1000 / 2), 9999);
      mom_basetrack.m[1].pixelnum = 0;
      mom_basetrack.m[1].hitnum = 0;

      mom_chain.base.push_back(mom_basetrack);

      mom_basetrack_pair.first = mom_basetrack_pair.second;
      mom_basetrack_pair.first.z = 0.;
      mom_basetrack_pair.second.pl = emulsion->GetPlate() + 1;
      mom_basetrack_pair.second.rawid = emulsion->GetEmulsionTrackId();
      mom_basetrack_pair.second.x = position.X();
      mom_basetrack_pair.second.y = position.Y();
      mom_basetrack_pair.second.z = position.Z() - downstream_position_z;
      mom_basetrack_pair.second.ax = tangent.X();
      mom_basetrack_pair.second.ay = tangent.Y();
      if ( emulsion != chain.front() )
	mom_chain.base_pair.push_back(mom_basetrack_pair);

      downstream_position_z = position.Z();

    } // emulsion
    ev.true_chains.push_back(mom_chain);
  } // chain

  return;
}

void ConnectionFunction::CalcPosInEccCoordinate(TVector3 &position, int ecc_id) const {

  // center of ECC5 dessicator
  position.SetX(position.X() - NINJA_POS_X - NINJA_ECC_POS_X);
  position.SetY(position.Y() - NINJA_POS_Y - NINJA_ECC_POS_Y);
  position.SetZ(position.Z() - NINJA_POS_Z - NINJA_ECC_POS_Z);

  // film coordinate
  position.SetX(position.X()
		+ 0.5 * NINJA_ECC_FILM_XY);
  position.SetY(position.Y()
		+ 0.5 * NINJA_DESIC_HEIGHT
		- NINJA_DESIC_THICK
		- NINJA_ENV_THICK);
  position.SetZ(position.Z()
		- 0.5 * NINJA_DESIC_DEPTH
		+ NINJA_DESIC_THICK
		+ NINJA_ENV_THICK
		+ NINJA_EMULSION_LAYER_THICK
		+ NINJA_BASE_LAYER_THICK);

  // move to each ECC
  position.SetX(position.X() 
		+ NINJA_ECC_GAP_X * (1 - ecc_id % 3));
  position.SetY(position.Y()
		+ NINJA_ECC_GAP_Y * (ecc_id / 3 - 1));

  // mm to um
  position.SetX(position.X() * 1000.);
  position.SetY(position.Y() * 1000.);
  position.SetZ(position.Z() * 1000.);

}

void ConnectionFunction::SmearEmulsions(std::vector<B2EmulsionSummary* > &emulsions_smeared,
					std::vector<const B2EmulsionSummary* > &emulsions) const {

  for ( auto emulsion : emulsions ) {
    B2EmulsionSummary *emulsion_tmp = new B2EmulsionSummary();

    int ecc_id = emulsion->GetEcc();

    double z_diff_smear = 0.;

    TVector3 position = emulsion->GetAbsolutePosition().GetValue();
    position = 1000. * position; // mm -> um
    z_diff_smear = -position.Z();
    SmearPosition(position);
    z_diff_smear += position.Z(); // z difference before & after smearing ( smeared - true )
    z_diff_smear *= 1.e-3; // um -> mm
    TVector3 tangent  = emulsion->GetTangent().GetValue();
    double sign = tangent.Z();
    tangent = sign * tangent;
    SmearTangentVector(tangent);
    emulsion_tmp->SetEmulsionTrackId(emulsion->GetEmulsionTrackId());
    emulsion_tmp->SetParentTrack(emulsion->GetParentTrack());
    // emulsion_tmp->SetParentTrackId(emulsion->GetParentTrackId());
    emulsion_tmp->SetMomentum(emulsion->GetMomentum());
    emulsion_tmp->SetAbsolutePosition(1.e-3 * position); // um -> mm
    position.SetX(position.X() + 1000. * (0.5 * NINJA_ECC_FILM_XY
					  - NINJA_POS_X - NINJA_ECC_POS_X
					  - NINJA_ECC_GAP_X * (std::floor(ecc_id % 3) - 1)));
    position.SetY(position.Y() + 1000. * (0.5 * NINJA_RACK_HEIGHT
					  - NINJA_POS_Y - (1030. + NINJA_DESIC_THICK + NINJA_ENV_THICK)
					  - NINJA_ECC_GAP_Y * (1 - std::floor(ecc_id / 3.))));
    position.SetZ(0.);
    emulsion_tmp->SetFilmPosition(1.e-3 * position); // um -> mm
    double filmpos_in_down_z;
    if ( emulsion->GetPlate() == 2 ) filmpos_in_down_z = z_diff_smear;
    else if ( emulsion->GetPlate() == 3 ) filmpos_in_down_z = - NINJA_FILM_THICK + z_diff_smear;
    else if ( emulsion->GetPlate() < 15 ) filmpos_in_down_z = - NINJA_FILM_THICK - NINJA_IRON_LAYER_THICK + z_diff_smear;
    else if ( emulsion->GetPlate() == 15 ) filmpos_in_down_z = - NINJA_FILM_THICK - 2 * NINJA_ENV_THICK + z_diff_smear;
    else if ( emulsion->GetPlate() % 2 == 0 ) filmpos_in_down_z = - NINJA_FILM_THICK - NINJA_IRON_LAYER_THICK + z_diff_smear;
    else filmpos_in_down_z = - NINJA_FILM_THICK - NINJA_WATER_LAYER_THICK - 2 * NINJA_ENV_THICK + z_diff_smear;
    position.SetZ(filmpos_in_down_z * 1000.);
    emulsion_tmp->SetFilmPositionInDownCoordinate(1.e-3 * position); // um -> mm
    emulsion_tmp->SetTangent((1./sign) * tangent);
    emulsion_tmp->SetTangentInDownCoordinate((1./sign) * tangent);
    emulsion_tmp->SetEdepSum(emulsion_tmp->GetEdepSum());
    emulsion_tmp->SetEdepDiff(emulsion_tmp->GetEdepDiff());
    emulsion_tmp->SetFilmType(B2EmulsionType::kECC);
    emulsion_tmp->SetEcc(ecc_id);
    emulsion_tmp->SetPlate(emulsion->GetPlate());
    emulsions_smeared.push_back(emulsion_tmp);
  }
  return;
}

void ConnectionFunction::SmearPosition(TVector3 &position /*um*/) const {
  position.SetX(gRandom->Gaus(position.X(), 0.2));
  position.SetY(gRandom->Gaus(position.Y(), 0.2));
  position.SetZ(gRandom->Gaus(position.Z(), 2.));
}

void ConnectionFunction::ApplyDetectionEfficiency(std::vector<B2EmulsionSummary* > &emulsions_detected,
						  std::vector<B2EmulsionSummary* > &emulsions_smeared,
						  int ecc_id) const {
  if ( ecc_id < 0 || ecc_id > 8 )
    throw std::invalid_argument("ECC id not valid : " + ecc_id);

  for ( auto emulsion : emulsions_smeared ) {
    TVector3 tangent = emulsion->GetTangent().GetValue();
    if ( std::fabs(tangent.X()) < 4.0 &&
	 std::fabs(tangent.Y()) < 4.0 ) {
      emulsions_detected.push_back(emulsion);
      BOOST_LOG_TRIVIAL(trace) << "Emulsion : " << emulsion->GetEmulsionTrackId()
			       << " is detected : "
			       << " PL" << emulsion->GetPlate() + 1;
    }
  }

  return;

}

void ConnectionFunction::ApplyFVCut(std::vector<B2EmulsionSummary* > &emulsions_detected_in_fv,
				    std::vector<B2EmulsionSummary* > &emulsions_detected,
				    int ecc_id) const {
  if ( ecc_id < 0 || ecc_id > 8 )
    throw std::invalid_argument("ECC id not valid : " + ecc_id);

  for ( auto emulsion : emulsions_detected ) {
    std::vector<FiducialArea> fa = ecc_fiducial_[ecc_id].at(emulsion->GetPlate() + 1);
    if ( JudgeFiducialArea(fa, emulsion) ) {
      emulsions_detected_in_fv.push_back(emulsion);
      BOOST_LOG_TRIVIAL(trace) << "Emulsion : " << emulsion->GetEmulsionTrackId()
			       << " is in fiducial volume : "
			       << " PL" << emulsion->GetPlate() + 1
			       << " ("
			       << emulsion->GetFilmPosition().GetValue().X() << ", "
			       << emulsion->GetFilmPosition().GetValue().Y() << ")";
    }
  }
  
  return;
    
}

void ConnectionFunction::GenerateLinklet(std::vector<std::pair<B2EmulsionSummary*, B2EmulsionSummary* > > &linklet,
					 std::vector<B2EmulsionSummary* > &emulsions) const {

  for ( auto emulsion : emulsions ) {
    // Black 判定
    // Black だった場合は同じ track 由来に限りつなぐ
    bool is_black_ = false;
    int parent_track_id_ = emulsion->GetParentTrackId();
    if ( emulsion->GetParentTrack().GetParticlePdg() == PDG_t::kProton ||
	 ( B2Pdg::IsChargedPion(emulsion->GetParentTrack().GetParticlePdg()) &&
	   emulsion->GetParentTrack().GetFinalAbsoluteMomentum().GetValue() < 100. ) ) {
      is_black_ = true;
    }
    // plate を確認
    int plate = emulsion->GetPlate();
    for ( auto jemulsion : emulsions) {
      if ( plate < 12 ) { // 鉄ECC
	if ( jemulsion->GetPlate() == plate + 1 ) {
	  if ( JudgeConnect(emulsion, jemulsion, fe_param_) ) {
	    linklet.push_back(std::make_pair(emulsion, jemulsion));
	  }
	  else if ( is_black_ && jemulsion->GetParentTrackId() == parent_track_id_ &&
		    JudgeConnect(emulsion, jemulsion, black_fe_param_) ) {
	    linklet.push_back(std::make_pair(emulsion, jemulsion));
	  }
	}
	else if ( jemulsion->GetPlate() == plate + 2 ) {
	  if ( JudgeConnect(emulsion, jemulsion, fe_fe_param_) ) {
	    linklet.push_back(std::make_pair(emulsion, jemulsion));
	  }
	}
	else if ( jemulsion->GetPlate() == plate + 3 ) {
	  if ( JudgeConnect(emulsion, jemulsion, fe_fe_fe_param_) ) {
	    linklet.push_back(std::make_pair(emulsion, jemulsion));
	  }
	}
	else continue;
      }
      else if ( plate == 12 ) { // 鉄ECC-水ECCまたぎ
	if ( jemulsion->GetPlate() == 13 ) {
	  if ( JudgeConnect(emulsion, jemulsion, fe_param_) ) {
	    linklet.push_back(std::make_pair(emulsion, jemulsion));
	  }
	  else if ( is_black_ && jemulsion->GetParentTrackId() == parent_track_id_ &&
		    JudgeConnect(emulsion, jemulsion, black_fe_param_) ) {
	    linklet.push_back(std::make_pair(emulsion, jemulsion));
	  }
	}
	else if ( jemulsion->GetPlate() == 14 ) {
	  if ( JudgeConnect(emulsion, jemulsion, fe_fe_param_) ) {
	    linklet.push_back(std::make_pair(emulsion, jemulsion));
	  }	  
	}
	else if ( jemulsion->GetPlate( ) == 15 ) {
	  if ( JudgeConnect(emulsion, jemulsion, et_ecc_fe_fe_param_) ) {
	    linklet.push_back(std::make_pair(emulsion, jemulsion));
	  }	  
	}
	else continue;
      }
      else if ( plate == 13 ) {
	if ( jemulsion->GetPlate() == 14 ) {
	  if ( JudgeConnect(emulsion, jemulsion, fe_param_) ) {
	    linklet.push_back(std::make_pair(emulsion, jemulsion));
	  }
	  else if ( is_black_ && jemulsion->GetParentTrackId() == parent_track_id_ &&
		    JudgeConnect(emulsion, jemulsion, black_fe_param_) ) {
	    linklet.push_back(std::make_pair(emulsion, jemulsion));
	  }
	}
	else if ( jemulsion->GetPlate() == 15 ) {
	  if ( JudgeConnect(emulsion, jemulsion, et_ecc_fe_param_) ) {
	    linklet.push_back(std::make_pair(emulsion, jemulsion));
	  }
	}
	else if ( jemulsion->GetPlate() == 16 ) {
	  if ( JudgeConnect(emulsion, jemulsion, et_ecc_fe_fe_param_) ) {
	    linklet.push_back(std::make_pair(emulsion, jemulsion));
	  }
	}
	else continue;
      }    
      else if ( plate == 14 ) {
	if ( jemulsion->GetPlate() == 15 ) {
	  if ( JudgeConnect(emulsion, jemulsion, et_ecc_param_) ) {
	    linklet.push_back(std::make_pair(emulsion, jemulsion));
	  }
	}
	else if ( jemulsion->GetPlate() == 16 ) {
	  if ( JudgeConnect(emulsion, jemulsion, et_ecc_fe_param_) ) {
	    linklet.push_back(std::make_pair(emulsion, jemulsion));;
	  }
	}
	else continue;
      }
      else if ( plate == 130 ) { // 最上流例外処理
	if ( jemulsion->GetPlate() == 131 ) {
	  if ( JudgeConnect(emulsion, jemulsion, water_param_) ) {
	    linklet.push_back(std::make_pair(emulsion, jemulsion));
	  }
	  else if ( is_black_ && jemulsion->GetParentTrackId() == parent_track_id_ &&
		    JudgeConnect(emulsion, jemulsion, black_water_param_) ) {
	    linklet.push_back(std::make_pair(emulsion, jemulsion));
	  }
	}
	else if ( jemulsion->GetPlate() == 132 ) {
	  if ( JudgeConnect(emulsion, jemulsion, fe_water_param_) ) {
	    linklet.push_back(std::make_pair(emulsion, jemulsion));
	  }
	}
	else continue;
      }
      else if ( plate == 131 ) {
	if ( jemulsion->GetPlate() == 132 ) {
	  if ( JudgeConnect(emulsion, jemulsion, fe_param_) ) {
	    linklet.push_back(std::make_pair(emulsion, jemulsion));
	  }
	  else if ( is_black_ && jemulsion->GetParentTrackId() == parent_track_id_ &&
		    JudgeConnect(emulsion, jemulsion, black_fe_param_) ) {
	    linklet.push_back(std::make_pair(emulsion, jemulsion));
	  }
	}
	else continue;
      }      
      else if ( plate % 2 == 1 ) { // 鉄の下流
	if ( jemulsion->GetPlate() == plate + 1 ) {
	  if ( JudgeConnect(emulsion, jemulsion, fe_param_) ) {
	    linklet.push_back(std::make_pair(emulsion, jemulsion));
	  }
	  else if ( is_black_ && jemulsion->GetParentTrackId() == parent_track_id_ &&
		    JudgeConnect(emulsion, jemulsion, black_fe_param_) ) {
	    linklet.push_back(std::make_pair(emulsion, jemulsion));
	  }
	}
	else if ( jemulsion->GetPlate() == plate + 2 ) {
	  if ( JudgeConnect(emulsion, jemulsion, fe_water_param_) ) {
	    linklet.push_back(std::make_pair(emulsion, jemulsion));
	  }
	}
	else if ( jemulsion->GetPlate() == plate + 3 ) {
	  if ( JudgeConnect(emulsion, jemulsion, fe_water_fe_param_) ) {
	    linklet.push_back(std::make_pair(emulsion, jemulsion));
	  }	  
	}
	else continue;
      }    
      else if ( plate % 2 == 0 ) { // 水の下流
	if ( jemulsion->GetPlate() == plate + 1 ) {
	  if ( JudgeConnect(emulsion, jemulsion, water_param_) ) {
	    linklet.push_back(std::make_pair(emulsion, jemulsion));
	  }
	  else if ( is_black_ && jemulsion->GetParentTrackId() == parent_track_id_ &&
		    JudgeConnect(emulsion, jemulsion, black_water_param_) ) {
	    linklet.push_back(std::make_pair(emulsion, jemulsion));
	  }
	}
	else if ( jemulsion->GetPlate() == plate + 2 ) {
	  if ( JudgeConnect(emulsion, jemulsion, fe_water_param_) ) {
	    linklet.push_back(std::make_pair(emulsion, jemulsion));
	  }
	}
	else continue;
      }
    }
  }

  return;

}

void ConnectionFunction::GenerateGroup(std::vector<std::pair<B2EmulsionSummary*, std::vector<std::pair<B2EmulsionSummary*, B2EmulsionSummary* > > > > &groups,
				       std::vector<std::pair<B2EmulsionSummary*, B2EmulsionSummary* > > &linklets,
				       std::vector<B2EmulsionSummary* > &emulsions) const {
  boost::unordered_multimap<Segment, Segment > link_next;
  boost::unordered_multimap<Segment, Segment > link_prev;
  for ( auto linklet : linklets ) {
    std::pair<Segment, Segment > linklet_seg;
    if ( linklet.first->GetPlate() > linklet.second->GetPlate() ) 
      std::swap(linklet.first, linklet.second);
    LinkletConvertB2ToSeg(linklet, linklet_seg);
    link_next.insert(std::make_pair(linklet_seg.first, linklet_seg.second));
    link_next.insert(std::make_pair(linklet_seg.second, linklet_seg.first));
  }

  std::vector<Group > groups_seg;
  for ( auto emulsion : emulsions ) {
    std::vector<std::pair<Segment, Segment > > group_next, group_prev;
    Group group;
    group.start_plate = emulsion->GetPlate() + 1;
    group.start_rawid = emulsion->GetEmulsionTrackId();
    GenerateGroupNext(emulsion, link_next, link_prev, group_next);
    GenerateGroupPrev(emulsion, link_next, link_prev, group_prev);
    GroupMerge(emulsion, group_next, group_prev, group);
    groups_seg.push_back(group);
  }

  for ( auto group_seg : groups_seg ) {
    std::vector<std::pair<B2EmulsionSummary*, B2EmulsionSummary* > > linklets;
    for ( auto link_seg : group_seg.linklets ) {
      std::pair<B2EmulsionSummary*, B2EmulsionSummary* > linklet;
      LinkletConvertSegToB2(link_seg, linklet, emulsions);
      linklets.push_back(linklet);
    }
    B2EmulsionSummary* start_emulsion;
    ConvertRawIdToB2(group_seg.start_rawid, start_emulsion, emulsions);
    groups.push_back(std::make_pair(start_emulsion, linklets));
  }

  return;

}

void ConnectionFunction::ConvertRawIdToB2(unsigned int rawid, B2EmulsionSummary* &emulsion,
					  std::vector<B2EmulsionSummary*> &emulsions ) const {
  for ( auto iemulsion : emulsions ) {
    if ( iemulsion->GetEmulsionTrackId() == rawid ) 
      emulsion = iemulsion;
  }

  return;
}

void ConnectionFunction::LinkletConvertB2ToSeg(std::pair<B2EmulsionSummary*, B2EmulsionSummary* > &linklet_b2,
					       std::pair<Segment, Segment > &linklet_seg) const {
  Segment seg0, seg1;
  seg0.plate = linklet_b2.first->GetPlate() + 1;
  seg0.rawid = linklet_b2.first->GetEmulsionTrackId();
  seg1.plate = linklet_b2.second->GetPlate() + 1;
  seg1.rawid = linklet_b2.second->GetEmulsionTrackId();
  linklet_seg = std::make_pair(seg0, seg1);
  return;
}

void ConnectionFunction::LinkletConvertSegToB2(std::pair<Segment, Segment > &linklet_seg,
					       std::pair<B2EmulsionSummary*, B2EmulsionSummary* > &linklet_b2,
					       std::vector<B2EmulsionSummary* > &emulsions) const {
  Segment seg0, seg1;
  seg0 = linklet_seg.first;
  seg1 = linklet_seg.second;
  for ( auto emulsion : emulsions ) {
    if ( seg0.plate == emulsion->GetPlate() + 1 &&
	 seg0.rawid == emulsion->GetEmulsionTrackId() ) { linklet_b2.first = emulsion; }
    else if ( seg1.plate == emulsion->GetPlate() + 1 &&
	      seg1.rawid == emulsion->GetEmulsionTrackId() ) { linklet_b2.second = emulsion; }    
  }
  return;
}

void ConnectionFunction::GroupConvertB2ToSeg(std::pair<B2EmulsionSummary*, std::vector<std::pair<B2EmulsionSummary*, B2EmulsionSummary* > > > &group_b2,
					     Group &group_seg) const {
  B2EmulsionSummary* start_b2 = group_b2.first;
  group_seg.start_plate = start_b2->GetPlate() + 1;
  group_seg.start_rawid = start_b2->GetEmulsionTrackId();
  for ( auto linklet : group_b2.second ) {
    std::pair<Segment, Segment > linklet_seg;
    LinkletConvertB2ToSeg(linklet, linklet_seg);
    group_seg.linklets.push_back(linklet_seg);
  }
  return;
}
 
void ConnectionFunction::GroupConvertSegToB2(Group &group_seg,
					     std::pair<B2EmulsionSummary*, std::vector<std::pair<B2EmulsionSummary*, B2EmulsionSummary*> > > &group_b2,
					     std::vector<B2EmulsionSummary* > &emulsions) const {
  B2EmulsionSummary* start_b2;
  ConvertRawIdToB2(group_seg.start_rawid, start_b2, emulsions);
  group_b2.first = start_b2;
  std::vector<std::pair<B2EmulsionSummary*, B2EmulsionSummary* > > linklets_b2;
  for ( auto linklet : group_seg.linklets ) {
    std::pair<B2EmulsionSummary*, B2EmulsionSummary* > linklet_b2;
    LinkletConvertSegToB2(linklet, linklet_b2, emulsions);
    linklets_b2.push_back(linklet_b2);
  }
  group_b2.second = linklets_b2;
  return;
}

void ConnectionFunction::GenerateGroupNext(B2EmulsionSummary* emulsion,
					   boost::unordered_multimap<Segment, Segment > &link_next,
					   boost::unordered_multimap<Segment, Segment > &link_prev,
					   std::vector<std::pair<Segment, Segment > > &group_next) const {

  std::set<std::pair<Segment, Segment > > all_link;
  std::set<Segment > seen;
  std::set<Segment > add, now;

  Segment seg;
  seg.plate = emulsion->GetPlate() + 1;
  seg.rawid = emulsion->GetEmulsionTrackId();

  add.insert(seg);
  while ( !add.empty() ) {
    now = add;
    add.clear();

    for ( auto itr = now.begin(); itr != now.end(); itr++ ) {
      seen.insert(*itr);
    }

    for ( auto itr = now.begin(); itr != now.end(); itr++ ) {
      if ( link_next.find(*itr) == link_next.end() ) continue;
      auto range = link_next.equal_range(*itr);
      for ( auto res = range.first; res != range.second; res++ ) {
	if ( seen.count(res->second) == 0 )
	  add.insert(res->second);
	if ( res->first.plate < res->second.plate ) {
	  all_link.insert(std::make_pair(res->first, res->second));
	} else {
	  all_link.insert(std::make_pair(res->second, res->first));
	}
      }
    }
  }

  seen.clear();

  for ( auto itr = all_link.begin(); itr != all_link.end(); itr++ ) {
    add.insert(itr->first);
    add.insert(itr->second);
  }

  while ( !add.empty() ) {
    now = add;
    add.clear();
    for ( auto itr = now.begin(); itr != now.end(); itr++ ) {
      seen.insert(*itr);
    }
    for ( auto itr = now.begin(); itr != now.end(); itr++ ) {
      if ( link_prev.find(*itr) == link_prev.end() ) continue;
      auto range = link_prev.equal_range(*itr);
      for ( auto res = range.first; res != range.second; res++ ) {
	if ( res->second.plate <= emulsion->GetPlate() + 1 ) continue;
	
	if ( seen.count(res->second) == 0 ) {
	  add.insert(res->second);
	}
	if ( res->first.plate < res->second.plate ) {
	  all_link.insert(std::make_pair(res->first, res->second));
	} else {
	  all_link.insert(std::make_pair(res->second, res->first));
	}
      }
    }
  }

  if ( all_link.size() ) {
    for ( auto itr = all_link.begin(); itr != all_link.end(); itr++ ) {
      group_next.push_back(std::make_pair(itr->first, itr->second));
    }
  }

  return;
  
}

void ConnectionFunction::GenerateGroupPrev(B2EmulsionSummary* emulsion,
					   boost::unordered_multimap<Segment, Segment > &link_next,
					   boost::unordered_multimap<Segment, Segment > &link_prev,
					   std::vector<std::pair<Segment, Segment > > &group_prev) const {

  std::set<std::pair<Segment, Segment > > all_link; // group_prev に含むべき全linklet
  std::set<Segment > seen; // すでにall_link に入るか確認したbasetrack
  std::set<Segment > add, now; // add : now :

  Segment seg;
  seg.plate = emulsion->GetPlate() + 1;
  seg.rawid = emulsion->GetEmulsionTrackId();

  add.insert(seg);
  while ( !add.empty() ) {
    now = add;
    add.clear();

    for ( auto itr = now.begin(); itr != now.end(); itr++ ) {
      seen.insert(*itr);
    }

    for ( auto itr = now.begin(); itr != now.end(); itr++ ) {
      if ( link_prev.find(*itr) == link_prev.end() ) continue;
      auto range = link_prev.equal_range(*itr);
      for ( auto res = range.first; res != range.second; res++ ) {
	if ( seen.count(res->second) == 0 )
	  add.insert(res->second);
	if ( res->first.plate < res->second.plate ) {
	  all_link.insert(std::make_pair(res->first, res->second));
	} else {
	  all_link.insert(std::make_pair(res->second, res->first));
	}
      }
    }
  }

  seen.clear();

  for ( auto itr = all_link.begin(); itr != all_link.end(); itr++ ) {
    add.insert(itr->first);
    add.insert(itr->second);
  }

  while ( !add.empty() ) {
    now = add;
    add.clear();
    for ( auto itr = now.begin(); itr != now.end(); itr++ ) {
      seen.insert(*itr);
    }
    for ( auto itr = now.begin(); itr != now.end(); itr++ ) {
      if (link_next.find(*itr) == link_next.end() ) continue;
      auto range = link_next.equal_range(*itr);
      for ( auto res = range.first; res != range.second; res++ ) {
	if ( res->second.plate >= emulsion->GetPlate() + 1 ) continue;
	
	if ( seen.count(res->second) == 0 ) {
	  add.insert(res->second);
	}
	if ( res->first.plate < res->second.plate ) {
	  all_link.insert(std::make_pair(res->first, res->second));
	} else {
	  all_link.insert(std::make_pair(res->second, res->first));
	}
      }
    }
  }

  if ( all_link.size() ) {
    for ( auto itr = all_link.begin(); itr != all_link.end(); itr++ ) {
      group_prev.push_back(std::make_pair(itr->first, itr->second));
    }
  }

  return;

}

void ConnectionFunction::GroupMerge(B2EmulsionSummary* emulsion,
				    std::vector<std::pair<Segment, Segment > > &group_next,
				    std::vector<std::pair<Segment, Segment > > &group_prev,
				    Group &group) const {
  std::set<std::pair<Segment, Segment > > linklet_all;
  for ( int i = 0; i < group_next.size(); i++ ) {
    linklet_all.insert(group_next.at(i));
  }
  for ( int i = 0; i < group_prev.size(); i++ ) {
    linklet_all.insert(group_prev.at(i));
  }

  if ( linklet_all.empty() ) {
    group.linklets.clear();
    group.linklets.shrink_to_fit();
    return;
  }

  group.start_plate = emulsion->GetPlate() + 1;
  group.start_rawid = emulsion->GetEmulsionTrackId();
  for ( auto itr = linklet_all.begin(); itr != linklet_all.end(); itr++ ) {
    group.linklets.push_back(*itr);
  }
  return;
}

void ConnectionFunction::ReconnectGroups(std::vector<std::pair<B2EmulsionSummary*, std::vector<std::pair<B2EmulsionSummary*, B2EmulsionSummary* > > > > &groups_reconnected,
					 std::vector<std::pair<B2EmulsionSummary*, std::vector<std::pair<B2EmulsionSummary*, B2EmulsionSummary* > > > > &groups,
					 std::vector<B2EmulsionSummary* > &emulsions) const {

  // group を Group に変換
  std::vector<Segment > start_vec;
  std::map<Segment, Group > groups_seg_map;

  for ( auto group : groups ) {
    Group group_seg;
    GroupConvertB2ToSeg(group, group_seg);
    Segment start_seg;
    start_seg.plate = group_seg.start_plate;
    start_seg.rawid = group_seg.start_rawid;
    start_vec.push_back(start_seg);
    groups_seg_map.insert(std::make_pair(start_seg, group_seg));
  }

  // reconnection するための track candidate を探す
  boost::unordered_multimap<Segment, Segment > connect_cand_ups; // 最上流候補
  boost::unordered_multimap<Segment, Segment > connect_cand_downs; // 最下流候補
  boost::unordered_multimap<Segment, Segment > connect_cand_ups_start; // 最上流候補から group を検索
  boost::unordered_multimap<Segment, Segment > connect_cand_downs_start; // 最下流候補から group を検索

  // 各 group について最上流・最下流の basetrack を探索
  // 見つけたら unordered_multimap に詰める
  for ( auto mapitr = groups_seg_map.begin(); mapitr != groups_seg_map.end(); mapitr++ ) {
    Segment start_seg = mapitr->first;
    Group group_seg = mapitr->second;
    CollectEdgeTracks(connect_cand_ups, connect_cand_downs,
		      connect_cand_ups_start, connect_cand_downs_start,
		      group_seg, start_seg);
  }

  // 追加するべき linklet の集合
  std::set<std::pair<Segment, Segment > > linklets_appended;

  // 各 connect_cand_up について connect_cand_down からつながる飛跡があるかを確認
  // あれば set につめる
  for ( auto mapitr = groups_seg_map.begin(); mapitr != groups_seg_map.end(); mapitr++ ) {
    Segment seg = mapitr->first;
    Group target_group = mapitr->second;

    auto range = connect_cand_ups.equal_range(seg);

    for ( auto itr = range.first; itr != range.second; itr++ ) {
      
      Segment connect_cand_up_seg = itr->second;
      B2EmulsionSummary* connect_cand_up;
      ConvertRawIdToB2(connect_cand_up_seg.rawid,
		       connect_cand_up, emulsions);

      // black の判定
      bool is_black_ = false;
      int parent_track_id_ = connect_cand_up->GetParentTrackId();
      if ( connect_cand_up->GetParentTrack().GetParticlePdg() == PDG_t::kProton ||
	   ( B2Pdg::IsChargedPion(connect_cand_up->GetParentTrack().GetParticlePdg()) &&
	    connect_cand_up->GetParentTrack().GetFinalAbsoluteMomentum().GetValue() < 100. ) ) {
	is_black_ = true;
      }
      // plate を確認
      int plate = connect_cand_up->GetPlate() + 1;

      // 接続探索候補は connect_cand_downs の中で start_seg が違うものすべて
      for ( auto jmapitr = groups_seg_map.begin(); jmapitr != groups_seg_map.end(); jmapitr++ ) {
	Segment jseg = jmapitr->first;
	Group connect_group = jmapitr->second;
	if ( jseg == seg ) continue;

	auto jrange = connect_cand_downs.equal_range(jseg);
	
	for ( auto jitr = jrange.first; jitr != jrange.second; jitr++ ) {
	  Segment connect_cand_down_seg = jitr->second;
	  B2EmulsionSummary* connect_cand_down;
	  ConvertRawIdToB2(connect_cand_down_seg.rawid,
			   connect_cand_down, emulsions);
	  if ( plate == 3 ) { continue; } // 2pl 以上の chain のみが対象
	  else if ( plate < 11 ) { // 鉄ECC
	    if ( connect_cand_down->GetPlate() == plate + 1 ) {
	      if ( JudgeConnect(connect_cand_up, connect_cand_down, re_fe_param_) ) {
		linklets_appended.insert(std::make_pair(connect_cand_up_seg, connect_cand_down_seg));
	      }
	      else if ( is_black_ && connect_cand_down->GetParentTrackId() == parent_track_id_ &&
			JudgeConnect(connect_cand_up, connect_cand_down, black_re_fe_param_) ) {
		linklets_appended.insert(std::make_pair(connect_cand_up_seg, connect_cand_down_seg));
	      }
	    }
	    else if ( connect_cand_down->GetPlate() == plate + 2 ) {
	      if ( JudgeConnect(connect_cand_up, connect_cand_down, re_fe_fe_param_) ) {
		linklets_appended.insert(std::make_pair(connect_cand_up_seg, connect_cand_down_seg));
	      }
	      else if ( is_black_ && connect_cand_down->GetParentTrackId() == parent_track_id_ &&
			JudgeConnect(connect_cand_up, connect_cand_down, black_re_fe_fe_param_) ) {
		linklets_appended.insert(std::make_pair(connect_cand_up_seg, connect_cand_down_seg));
	      }
	    }
	    else if ( connect_cand_down->GetPlate() == plate + 3 ) {
	      if ( JudgeConnect(connect_cand_up, connect_cand_down, re_fe_fe_fe_param_) ) {
		linklets_appended.insert(std::make_pair(connect_cand_up_seg, connect_cand_down_seg));
	      }
	    }
	    else if ( connect_cand_down->GetPlate() == plate + 4 ) {
	      if ( JudgeConnect(connect_cand_up, connect_cand_down, re_fe_fe_fe_fe_param_) ) {
		linklets_appended.insert(std::make_pair(connect_cand_up_seg, connect_cand_down_seg));
	      }
	    }
	    else if ( connect_cand_down->GetPlate() == plate + 5 ) {
	      if ( JudgeConnect(connect_cand_up, connect_cand_down, re_fe_fe_fe_fe_fe_param_) ) {
		linklets_appended.insert(std::make_pair(connect_cand_up_seg, connect_cand_down_seg));
	      }
	    }
	    else continue;
	  }
	  else if ( plate == 11 ) { // 鉄ECC-水ECCまたぎ
	    if ( connect_cand_down->GetPlate() == plate + 1 ) {
	      if ( JudgeConnect(connect_cand_up, connect_cand_down, re_fe_param_) ) {
		linklets_appended.insert(std::make_pair(connect_cand_up_seg, connect_cand_down_seg));
	      }
	      else if ( is_black_ && connect_cand_down->GetParentTrackId() == parent_track_id_ &&
			JudgeConnect(connect_cand_up, connect_cand_down, black_re_fe_param_) ) {
		linklets_appended.insert(std::make_pair(connect_cand_up_seg, connect_cand_down_seg));
	      }
	    }
	    else if ( connect_cand_down->GetPlate() == plate + 2 ) {
	      if ( JudgeConnect(connect_cand_up, connect_cand_down, re_fe_fe_param_) ) {
		linklets_appended.insert(std::make_pair(connect_cand_up_seg, connect_cand_down_seg));
	      }
	      else if ( is_black_ && connect_cand_down->GetParentTrackId() == parent_track_id_ &&
			JudgeConnect(connect_cand_up, connect_cand_down, black_re_fe_fe_param_) ) {
		linklets_appended.insert(std::make_pair(connect_cand_up_seg, connect_cand_down_seg));
	      }
	    }
	    else if ( connect_cand_down->GetPlate() == plate + 3 ) {
	      if ( JudgeConnect(connect_cand_up, connect_cand_down, re_fe_fe_fe_param_) ) {
		linklets_appended.insert(std::make_pair(connect_cand_up_seg, connect_cand_down_seg));
	      }
	    }
	    else if ( connect_cand_down->GetPlate() == plate + 4 ) {
	      if ( JudgeConnect(connect_cand_up, connect_cand_down, re_fe_fe_fe_fe_param_) ) {
		linklets_appended.insert(std::make_pair(connect_cand_up_seg, connect_cand_down_seg));
	      }
	    }
	    else if ( connect_cand_down->GetPlate() == plate + 5 ) {
	      if ( JudgeConnect(connect_cand_up, connect_cand_down, re_et_ecc_fe_fe_fe_fe_param_) ) {
		linklets_appended.insert(std::make_pair(connect_cand_up_seg, connect_cand_down_seg));
	      }
	    }
	    else continue;
	  }
	  else if ( plate == 12 ) {
	    if ( connect_cand_down->GetPlate() == plate + 1 ) {
	      if ( JudgeConnect(connect_cand_up, connect_cand_down, re_fe_param_) ) {
		linklets_appended.insert(std::make_pair(connect_cand_up_seg, connect_cand_down_seg));
	      }
	      else if ( is_black_ && connect_cand_down->GetParentTrackId() == parent_track_id_ &&
			JudgeConnect(connect_cand_up, connect_cand_down, black_re_fe_param_) ) {
		linklets_appended.insert(std::make_pair(connect_cand_up_seg, connect_cand_down_seg));
	      }
	    }
	    else if ( connect_cand_down->GetPlate() == plate + 2 ) {
	      if ( JudgeConnect(connect_cand_up, connect_cand_down, re_fe_fe_param_) ) {
		linklets_appended.insert(std::make_pair(connect_cand_up_seg, connect_cand_down_seg));
	      }
	      else if ( is_black_ && connect_cand_down->GetParentTrackId() == parent_track_id_ &&
			JudgeConnect(connect_cand_up, connect_cand_down, black_re_fe_fe_param_) ) {
		linklets_appended.insert(std::make_pair(connect_cand_up_seg, connect_cand_down_seg));
	      }
	    }
	    else if ( connect_cand_down->GetPlate() == plate + 3 ) {
	      if ( JudgeConnect(connect_cand_up, connect_cand_down, re_fe_fe_fe_param_) ) {
		linklets_appended.insert(std::make_pair(connect_cand_up_seg, connect_cand_down_seg));
	      }
	    }
	    else if ( connect_cand_down->GetPlate() == plate + 4 ) {
	      if ( JudgeConnect(connect_cand_up, connect_cand_down, re_et_ecc_fe_fe_fe_param_) ) {
		linklets_appended.insert(std::make_pair(connect_cand_up_seg, connect_cand_down_seg));
	      }
	    }
	    else if ( connect_cand_down->GetPlate() == plate + 5 ) {
	      if ( JudgeConnect(connect_cand_up, connect_cand_down, re_et_ecc_fe_fe_fe_fe_param_) ) {
		linklets_appended.insert(std::make_pair(connect_cand_up_seg, connect_cand_down_seg));
	      }
	    }
	    else continue;
	  }
	  else if ( plate == 13 ) {
	    if ( connect_cand_down->GetPlate() == plate + 1 ) {
	      if ( JudgeConnect(connect_cand_up, connect_cand_down, re_fe_param_) ) {
		linklets_appended.insert(std::make_pair(connect_cand_up_seg, connect_cand_down_seg));
	      }
	      else if ( is_black_ && connect_cand_down->GetParentTrackId() == parent_track_id_ &&
			JudgeConnect(connect_cand_up, connect_cand_down, black_re_fe_param_) ) {
		linklets_appended.insert(std::make_pair(connect_cand_up_seg, connect_cand_down_seg));
	      }
	    }
	    else if ( connect_cand_down->GetPlate() == plate + 2 ) {
	      if ( JudgeConnect(connect_cand_up, connect_cand_down, re_fe_fe_param_) ) {
		linklets_appended.insert(std::make_pair(connect_cand_up_seg, connect_cand_down_seg));
	      }
	      else if ( is_black_ && connect_cand_down->GetParentTrackId() == parent_track_id_ &&
			JudgeConnect(connect_cand_up, connect_cand_down, black_re_fe_fe_param_) ) {
		linklets_appended.insert(std::make_pair(connect_cand_up_seg, connect_cand_down_seg));
	      }
	    }
	    else if ( connect_cand_down->GetPlate() == plate + 4 ) {
	      if ( JudgeConnect(connect_cand_up, connect_cand_down, re_et_ecc_fe_fe_fe_fe_param_) ) {
		linklets_appended.insert(std::make_pair(connect_cand_up_seg, connect_cand_down_seg));
	      }
	    }
	    else if ( connect_cand_down->GetPlate() == plate + 5 ) {
	      if ( JudgeConnect(connect_cand_up, connect_cand_down, re_et_ecc_fe_fe_fe_water_param_) ) {
		linklets_appended.insert(std::make_pair(connect_cand_up_seg, connect_cand_down_seg));
	      }
	    }
	    else continue;
	  }
	  else if ( plate == 14 ) {
	    if ( connect_cand_down->GetPlate() == plate + 1 ) {
	      if ( JudgeConnect(connect_cand_up, connect_cand_down, re_fe_param_) ) {
		linklets_appended.insert(std::make_pair(connect_cand_up_seg, connect_cand_down_seg));
	      }
	      else if ( is_black_ && connect_cand_down->GetParentTrackId() == parent_track_id_ &&
			JudgeConnect(connect_cand_up, connect_cand_down, black_re_fe_param_) ) {
		linklets_appended.insert(std::make_pair(connect_cand_up_seg, connect_cand_down_seg));
	      }
	    }
	    else if ( connect_cand_down->GetPlate() == plate + 2 ) {
	      if ( is_black_ && connect_cand_down->GetParentTrackId() == parent_track_id_ &&
			JudgeConnect(connect_cand_up, connect_cand_down, black_re_fe_param_) ) {
		linklets_appended.insert(std::make_pair(connect_cand_up_seg, connect_cand_down_seg));
	      }
	    }
	    else if ( connect_cand_down->GetPlate() == plate + 4 ) {
	      if ( JudgeConnect(connect_cand_up, connect_cand_down, re_et_ecc_fe_fe_water_param_) ) {
		linklets_appended.insert(std::make_pair(connect_cand_up_seg, connect_cand_down_seg));
	      }
	    }
	    else if ( connect_cand_down->GetPlate() == plate + 5 ) {
	      if ( JudgeConnect(connect_cand_up, connect_cand_down, re_et_ecc_fe_fe_fe_water_param_) ) {
		linklets_appended.insert(std::make_pair(connect_cand_up_seg, connect_cand_down_seg));
	      }
	    }
	    else continue;
	  }
	  else if ( plate == 15 ) {
	    if ( connect_cand_down->GetPlate() == plate + 2 ) {
	      if ( is_black_ && connect_cand_down->GetParentTrackId() == parent_track_id_ &&
			JudgeConnect(connect_cand_up, connect_cand_down, black_re_fe_param_) ) {
		linklets_appended.insert(std::make_pair(connect_cand_up_seg, connect_cand_down_seg));
	      }
	    }
	    else if ( connect_cand_down->GetPlate() == plate + 3 ) {
	      if ( JudgeConnect(connect_cand_up, connect_cand_down, re_et_ecc_fe_water_param_) ) {
		linklets_appended.insert(std::make_pair(connect_cand_up_seg, connect_cand_down_seg));
	      }
	    }
	    else if ( connect_cand_down->GetPlate() == plate + 4 ) {
	      if ( JudgeConnect(connect_cand_up, connect_cand_down, re_et_ecc_fe_fe_water_param_) ) {
		linklets_appended.insert(std::make_pair(connect_cand_up_seg, connect_cand_down_seg));
	      }
	    }
	    else continue;
	  }
	  else if ( plate == 128 ) { // 最上流例外処理
	    if ( connect_cand_down->GetPlate() == plate + 1 ) {
	      if ( JudgeConnect(connect_cand_up, connect_cand_down, re_fe_param_) ) {
		linklets_appended.insert(std::make_pair(connect_cand_up_seg, connect_cand_down_seg));
	      }
	      else if ( is_black_ && connect_cand_down->GetParentTrackId() == parent_track_id_ &&
			JudgeConnect(connect_cand_up, connect_cand_down, black_re_fe_param_) ) {
		linklets_appended.insert(std::make_pair(connect_cand_up_seg, connect_cand_down_seg));
	      }
	    }
	    else if ( connect_cand_down->GetPlate() == plate + 2 ) {
	      if ( JudgeConnect(connect_cand_up, connect_cand_down, re_fe_water_param_) ) {
		linklets_appended.insert(std::make_pair(connect_cand_up_seg, connect_cand_down_seg));
	      }
	      else if ( is_black_ && connect_cand_down->GetParentTrackId() == parent_track_id_ &&
			JudgeConnect(connect_cand_up, connect_cand_down, black_re_fe_water_param_) ) {
		linklets_appended.insert(std::make_pair(connect_cand_up_seg, connect_cand_down_seg));
	      }	  
	    }
	    else if ( connect_cand_down->GetPlate() == plate + 3 ) {
	      if ( JudgeConnect(connect_cand_up, connect_cand_down, re_fe_water_fe_param_) ) {
		linklets_appended.insert(std::make_pair(connect_cand_up_seg, connect_cand_down_seg));
	      }
	    }
	    else if ( connect_cand_down->GetPlate() == plate + 4 ) {
	      if ( JudgeConnect(connect_cand_up, connect_cand_down, re_fe_water_fe_water_param_) ) {
		linklets_appended.insert(std::make_pair(connect_cand_up_seg, connect_cand_down_seg));
	      }	  
	    }
	    else continue;
	  }
	  else if ( plate == 129 ) {
	    if ( connect_cand_down->GetPlate() == plate + 2 ) {
	      if ( JudgeConnect(connect_cand_up, connect_cand_down, re_fe_water_param_) ) {
		linklets_appended.insert(std::make_pair(connect_cand_up_seg, connect_cand_down_seg));
	      }
	      else if ( is_black_ && connect_cand_down->GetParentTrackId() == parent_track_id_ &&
			JudgeConnect(connect_cand_up, connect_cand_down, black_re_fe_water_param_) ) {
		linklets_appended.insert(std::make_pair(connect_cand_up_seg, connect_cand_down_seg));
	      }
	    }
	    else if ( connect_cand_down->GetPlate() == plate + 3 ) {
	      if ( JudgeConnect(connect_cand_up, connect_cand_down, re_water_fe_water_param_) ) {
		linklets_appended.insert(std::make_pair(connect_cand_up_seg, connect_cand_down_seg));
	      }
	    } 
	    else continue;
	  }
	  else if ( plate == 130 ) {
	    if ( connect_cand_down->GetPlate() == plate + 1 ) {
	      if ( JudgeConnect(connect_cand_up, connect_cand_down, re_fe_param_) ) {
		linklets_appended.insert(std::make_pair(connect_cand_up_seg, connect_cand_down_seg));
	      }
	      else if ( is_black_ && connect_cand_down->GetParentTrackId() == parent_track_id_ &&
			JudgeConnect(connect_cand_up, connect_cand_down, black_re_fe_param_) ) {
		linklets_appended.insert(std::make_pair(connect_cand_up_seg, connect_cand_down_seg));
	      }
	    }
	    else if ( connect_cand_down->GetPlate() == plate + 2 ) {
	      if ( JudgeConnect(connect_cand_up, connect_cand_down, re_fe_water_param_) ) {
		linklets_appended.insert(std::make_pair(connect_cand_up_seg, connect_cand_down_seg));
	      }
	      else if ( is_black_ && connect_cand_down->GetParentTrackId() == parent_track_id_ &&
			JudgeConnect(connect_cand_up, connect_cand_down, black_re_fe_water_param_) ) {
		linklets_appended.insert(std::make_pair(connect_cand_up_seg, connect_cand_down_seg));
	      }
	    }
	    else continue;
	  }
	  else if ( plate == 131 ) { continue; }
	  else if ( plate == 132 ) { continue; }
	  else if ( plate % 2 == 0 ) { // 鉄下流 plate
	    if ( connect_cand_down->GetPlate() == plate + 1 ) {
	      if ( JudgeConnect(connect_cand_up, connect_cand_down, re_fe_param_) ) {
		linklets_appended.insert(std::make_pair(connect_cand_up_seg, connect_cand_down_seg));		
	      }
	      else if ( is_black_ && connect_cand_down->GetParentTrackId() == parent_track_id_ &&
			JudgeConnect(connect_cand_up, connect_cand_down, black_re_fe_param_) ) {
		linklets_appended.insert(std::make_pair(connect_cand_up_seg, connect_cand_down_seg));
	      }
	    }
	    else if ( connect_cand_down->GetPlate() == plate + 2 ) {
	      if ( JudgeConnect(connect_cand_up, connect_cand_down, re_fe_water_param_) ) {
		linklets_appended.insert(std::make_pair(connect_cand_up_seg, connect_cand_down_seg));
	      }
	      else if ( is_black_ && connect_cand_down->GetParentTrackId() == parent_track_id_ &&
			JudgeConnect(connect_cand_up, connect_cand_down, black_re_fe_water_param_) ) {
		linklets_appended.insert(std::make_pair(connect_cand_up_seg, connect_cand_down_seg));
	      }	  
	    }
	    else if ( connect_cand_down->GetPlate() == plate + 3 ) {
	      if ( JudgeConnect(connect_cand_up, connect_cand_down, re_fe_water_fe_param_) ) {
		linklets_appended.insert(std::make_pair(connect_cand_up_seg, connect_cand_down_seg));
	      }
	    }
	    else if ( connect_cand_down->GetPlate() == plate + 4 ) {
	      if ( JudgeConnect(connect_cand_up, connect_cand_down, re_fe_water_fe_water_param_) ) {
		linklets_appended.insert(std::make_pair(connect_cand_up_seg, connect_cand_down_seg));
	      }	  
	    }
	    else if ( connect_cand_down->GetPlate() == plate + 5 ) {
	      if ( JudgeConnect(connect_cand_up, connect_cand_down, re_fe_water_fe_water_fe_param_) ) {
		linklets_appended.insert(std::make_pair(connect_cand_up_seg, connect_cand_down_seg));
	      }	  
	    }
	    else continue;
	  }
	  else if ( plate % 2 == 1 ) {
	    if ( connect_cand_down->GetPlate() == plate + 2 ) {
	      if ( JudgeConnect(connect_cand_up, connect_cand_down, re_fe_water_param_) ) {
		linklets_appended.insert(std::make_pair(connect_cand_up_seg, connect_cand_down_seg));
	      }
	      else if ( is_black_ && connect_cand_down->GetParentTrackId() == parent_track_id_ &&
			JudgeConnect(connect_cand_up, connect_cand_down, black_re_fe_water_param_) ) {
		linklets_appended.insert(std::make_pair(connect_cand_up_seg, connect_cand_down_seg));
	      }	  
	    }
	    else if ( connect_cand_down->GetPlate() == plate + 3 ) {
	      if ( JudgeConnect(connect_cand_up, connect_cand_down, re_water_fe_water_param_) ) {
		linklets_appended.insert(std::make_pair(connect_cand_up_seg, connect_cand_down_seg));
	      }
	    }
	    else if ( connect_cand_down->GetPlate() == plate + 4 ) {
	      if ( JudgeConnect(connect_cand_up, connect_cand_down, re_fe_water_fe_water_param_) ) {
		linklets_appended.insert(std::make_pair(connect_cand_up_seg, connect_cand_down_seg));
	      }	  
	    }
	    else continue;
	  }
	} // jitr

      } // jmapitr

    } // itr

  } // mapitr


  // 各 start_seg を持つ group について linklets_appended の要素を持つか確認して
  // マージするべき set を作成
  std::map<Segment, std::set<std::pair<Segment, Segment > > > groups_reconnected_seg;
  for ( auto linklet : linklets_appended ) {

    // 対象の linklet を edge に持つ group を検索
    auto range_down = connect_cand_downs_start.equal_range(linklet.second);

    for ( auto itr_down = range_down.first; itr_down != range_down.second; itr_down++ ) {
      auto group_up = groups_seg_map.at(itr_down->second); // 最下流に linklet を持つ group は上流側
      Segment start_seg_up;
      start_seg_up.plate = group_up.start_plate;
      start_seg_up.rawid = group_up.start_rawid;
      auto range_up = connect_cand_ups_start.equal_range(linklet.first);
      for ( auto itr_up = range_up.first; itr_up != range_up.second; itr_up++ ) {
	auto group_down = groups_seg_map.at(itr_up->second);
	Segment start_seg_down;
	start_seg_down.plate = group_down.start_plate;
	start_seg_down.rawid = group_down.start_rawid;

	// 下流側 group
	std::set<std::pair<Segment, Segment > > linklets_down;
	if ( groups_reconnected_seg.count(start_seg_down) != 0 ) {
	  linklets_down = groups_reconnected_seg.at(start_seg_down);
	}
	linklets_down.insert(std::make_pair(linklet.first, linklet.second));
	for ( auto ilinklet : group_down.linklets ) {
	  linklets_down.insert(std::make_pair(ilinklet.first, ilinklet.second));
	}
	for ( auto ilinklet : group_up.linklets ) {
	  linklets_down.insert(std::make_pair(ilinklet.first, ilinklet.second));
	}
	if ( groups_reconnected_seg.count(start_seg_down) == 0)
	  groups_reconnected_seg.insert(std::make_pair(start_seg_down, linklets_down));
	else
	  groups_reconnected_seg.at(start_seg_down) = linklets_down;
	
	// 上流側 group
	std::set<std::pair<Segment, Segment > > linklets_up;
	if ( groups_reconnected_seg.count(start_seg_up) != 0 ) {
	  linklets_up = groups_reconnected_seg.at(start_seg_up);
	}
	linklets_up.insert(std::make_pair(linklet.first, linklet.second));
	for ( auto ilinklet : group_down.linklets ) {
	  linklets_up.insert(std::make_pair(ilinklet.first, ilinklet.second));
	}
	for ( auto ilinklet : group_up.linklets ) {
	  linklets_up.insert(std::make_pair(ilinklet.first, ilinklet.second));
	}
	if ( groups_reconnected_seg.count(start_seg_up) == 0)
	  groups_reconnected_seg.insert(std::make_pair(start_seg_up, linklets_down));
	else
	  groups_reconnected_seg.at(start_seg_up) = linklets_up;  
      }
    }

  }

  // B2 に変換
  for ( auto start : start_vec ) {
    B2EmulsionSummary* start_b2;
    ConvertRawIdToB2(start.rawid, start_b2, emulsions);
    auto group_seg = groups_seg_map.at(start);
    if ( groups_reconnected_seg.count(start) == 0 ) {
      std::vector<std::pair<B2EmulsionSummary*, B2EmulsionSummary* > > linklets_b2;
      for ( auto linklet_seg : group_seg.linklets ) {
	std::pair<B2EmulsionSummary*, B2EmulsionSummary* > linklet_b2;    
	LinkletConvertSegToB2(linklet_seg, linklet_b2, emulsions);
	linklets_b2.push_back(linklet_b2);
      }
      groups_reconnected.push_back(std::make_pair(start_b2, linklets_b2));
    }
    else {
      auto linklet_set = groups_reconnected_seg.at(start);
      std::vector<std::pair<B2EmulsionSummary*, B2EmulsionSummary* > > linklets_b2;
      for ( auto linklet_seg : linklet_set ) {
	std::pair<B2EmulsionSummary*, B2EmulsionSummary* > linklet_b2;    
	LinkletConvertSegToB2(linklet_seg, linklet_b2, emulsions);
	linklets_b2.push_back(linklet_b2);
      }
      groups_reconnected.push_back(std::make_pair(start_b2, linklets_b2));
    }
  }

  return;

}

void ConnectionFunction::CollectEdgeTracks(boost::unordered_multimap<Segment, Segment > &upstream_tracks,
					   boost::unordered_multimap<Segment, Segment > &downstream_tracks,
					   boost::unordered_multimap<Segment, Segment > &upstream_tracks_start,
					   boost::unordered_multimap<Segment, Segment > &downstream_tracks_start,
					   Group &group, Segment &start_seg) const {

  for ( auto linklet : group.linklets ) {
    bool is_most_up = true;
    bool is_most_down = true;
    for ( auto jlinklet : group.linklets ) {
      if ( linklet.first == jlinklet.second ) {
	is_most_down = false;
      }
      if ( linklet.second == jlinklet.first ) {
	is_most_up = false;
      }
      else continue;
    }

    if ( is_most_down ) {
      auto range = downstream_tracks.equal_range(start_seg);
      bool is_already_listed = false;
      for ( auto itr = range.first; itr != range.second; itr++ ) {
	if ( itr->second == linklet.first ) is_already_listed = true;
      }
      if ( !is_already_listed ) {
	BOOST_LOG_TRIVIAL(trace) << "Downstream : PL" << linklet.first.plate << ", " << linklet.first.rawid;
	downstream_tracks.insert(std::make_pair(start_seg, linklet.first));
	downstream_tracks_start.insert(std::make_pair(linklet.first, start_seg));
      }
    }
    if ( is_most_up ) {
      auto range = upstream_tracks.equal_range(start_seg);
      bool is_already_listed = false;
      for ( auto itr = range.first; itr != range.second; itr++ ) {
	if ( itr->second == linklet.second ) is_already_listed = true;
      }
      if ( !is_already_listed ) {
	BOOST_LOG_TRIVIAL(trace) << "Upstream : PL" << linklet.second.plate << ", " << linklet.second.rawid;
	upstream_tracks.insert(std::make_pair(start_seg, linklet.second));
	upstream_tracks_start.insert(std::make_pair(linklet.second, start_seg));
      }
    }

  }

  return;

}

bool ConnectionFunction::SelectMuonGroup(std::vector<std::pair<B2EmulsionSummary*, std::vector<std::pair<B2EmulsionSummary*, B2EmulsionSummary* > > > > &groups,
					 std::pair<B2EmulsionSummary*, std::vector<std::pair<B2EmulsionSummary*, B2EmulsionSummary* > > > &muon_group) const {

  std::vector<Group > groups_seg;
  for ( auto group : groups ) {
    Group group_seg;
    GroupConvertB2ToSeg(group, group_seg);
    groups_seg.push_back(group_seg);
  }

  boost::unordered_multimap<Segment, Segment > ups;
  boost::unordered_multimap<Segment, Segment > downs;
  boost::unordered_multimap<Segment, Segment > ups_start;
  boost::unordered_multimap<Segment, Segment > downs_start;
  for ( auto group_seg : groups_seg ) {
    Segment start_seg;
    start_seg.plate = group_seg.start_plate;
    start_seg.rawid = group_seg.start_rawid;
    CollectEdgeTracks(ups, downs, ups_start, downs_start,
		      group_seg, start_seg);
  }

  return false;  

}

/*
void ConnectionFunction::GenrerateGroupNextPartner(B2EmulsionSummary* emulsion,
						   boost::unordered_multimap<Segment, Segment > &link_next,
						   boost::unordered_multimap<Segment, Segment > &link_prev,
						   std::vector<std::pair<Segment, Segment > > &group_next_partner) const {
  std::set<std::pair<Segment, Segment > > all_link;
  std::set<Segment > seen;
  std::set<Segment > add, now;

  Segment seg;
  seg.plate = emulsion->GetPlate();
  seg.rawid = emulsion->GetEmulsionTrackId();

  add.insert(seg);
  while ( !add.empty() ) {
    now = add;
    add.clear();

    for ( auto itr = now.begin(); itr != now.end(); itr++ ) {
      seen.insert(*itr);
    }

    for ( auto itr = now.begin(); itr != now.end(); itr++ ) {
      if ( link_next.find(*itr) == link_next.end() ) continue;
      auto range = link_next.equal_range(*itr);
      for ( auto res = range.first; res != range.second; res++ ) {
	if ( seen.count(res->second) == 0 )
	  add.insert(res->second);
	if ( res->first.plate < res->second.plate ) {
	  all_link.insert(std::make_pair(res->first, res->second));
	} else {
	  all_link.insert(std::make_pair(res->second, res->first));
	}
      }
    }
  }

  seen.clear();

  for ( auto itr = all_link.begin(); itr != all_link.end(); itr++ ) {
    add.insert(itr->first);
    add.insert(itr->second);
  }

  while ( !add.empty() ) {
    now = add;
    add.clear();
    for ( auto itr = now.begin(); itr != now.end(); itr++ ) {
      seen.insert(*itr);
    }
    for ( auto itr = now.begin(); itr != now.end(); itr++ ) {
      if ( link_prev.find(*itr) == link_prev.end() ) continue;
      auto range = link_prev.equal_range(*itr);
      for ( auto res = range.first; res != range.second; res++ ) {
	if ( res->second.plate < emulsion->GetPlate() - 4 ) continue;
	
	if ( seen.count(res->second) == 0 ) {
	  add.insert(res->second);
	}
	if ( res->first.plate < res->second.plate ) {
	  all_link.insert(std::make_pair(res->first, res->second));
	} else {
	  all_link.insert(std::make_pair(res->second, res->first));
	}
      }
    }
  }

  if ( all_link.size() ) {
    for ( auto itr = all_link.begin(); itr != all_link.end(); itr++ ) {
      group_next.push_back(std::make_pair(itr->first, itr->second));
    }
  }

  return;

}


void ConnectionFunction::GenerateGroupPrevPartner(B2EmulsionSummary* emulsion,
						  boost::unordered_multimap<Segment, Segment > &link_next,
						  boost::unordered_multimap<Segment, Segment > &link_prev,
						  std::vector<std::pair<Segment, Segment > > &group_prev_partner) const {
  std::set<std::pair<Segment, Segment > > all_link; // group_prev に含むべき全linklet
  std::set<Segment > seen; // すでにall_link に入るか確認したbasetrack
  std::set<Segment > add, now; // add : now :

  Segment seg;
  seg.plate = emulsion->GetPlate();
  seg.rawid = emulsion->GetEmulsionTrackId();

  add.insert(seg);
  while ( !add.empty() ) {
    now = add;
    add.clear();

    for ( auto itr = now.begin(); itr != now.end(); itr++ ) {
      seen.insert(*itr);
    }

    for ( auto itr = now.begin(); itr != now.end(); itr++ ) {
      if ( link_prev.find(*itr) == link_prev.end() ) continue;
      auto range = link_prev.equal_range(*itr);
      for ( auto res = range.first; res != range.second; res++ ) {
	if ( seen.count(res->second) == 0 )
	  add.insert(res->second);
	if ( res->first.plate < res->second.plate ) {
	  all_link.insert(std::make_pair(res->first, res->second));
	} else {
	  all_link.insert(std::make_pair(res->second, res->first));
	}
      }
    }
  }

  seen.clear();

  for ( auto itr = all_link.begin(); itr != all_link.end(); itr++ ) {
    add.insert(itr->first);
    add.insert(itr->second);
  }

  while ( !add.empty() ) {
    now = add;
    add.clear();
    for ( auto itr = now.begin(); itr != now.end(); itr++ ) {
      seen.insert(*itr);
    }
    for ( auto itr = now.begin(); itr != now.end(); itr++ ) {
      if (link_next.find(*itr) == link_next.end() ) continue;
      auto range = link_next.equal_range(*itr);
      for ( auto res = range.first; res != range.second; res++ ) {
	if ( res->second.plate > emulsion->GetPlate() + 4 ) continue;
	
	if ( seen.count(res->second) == 0 ) {
	  add.insert(res->second);
	}
	if ( res->first.plate < res->second.plate ) {
	  all_link.insert(std::make_pair(res->first, res->second));
	} else {
	  all_link.insert(std::make_pair(res->second, res->first));
	}
      }
    }
  }

  if ( all_link.size() ) {
    for ( auto itr = all_link.begin(); itr != all_link.end(); itr++ ) {
      group_prev.push_back(std::make_pair(itr->first, itr->second));
    }
  }

  return;

}

*/
void ConnectionFunction::AddGroupsToEventInfo(Momentum_recon::Event_information &ev,
					      std::vector<std::vector<B2EmulsionSummary* > > &groups) const {
  Momentum_recon::Mom_chain mom_chain;
  Momentum_recon::Mom_basetrack mom_basetrack;
  std::pair<Momentum_recon::Mom_basetrack, Momentum_recon::Mom_basetrack > mom_basetrack_pair;

  for ( auto group : groups ) {
    int num_base = group.size();
    int num_link = num_base - 1;
    mom_chain.base.clear();
    mom_chain.base_pair.clear();
    mom_chain.base.reserve(num_base);
    mom_chain.base_pair.reserve(num_link);

    mom_chain.chainid = group.front()->GetParentTrackId();
    mom_chain.stop_flag = 0;
    mom_chain.particle_flag = group.front()->GetParentTrack().GetParticlePdg();
    if ( group.front()->GetPlate() == ev.vertex_pl )
      mom_chain.direction = 1;
    else if ( group.front()->GetPlate() == ev.vertex_pl + 1 )
      mom_chain.direction = -1;
    mom_chain.charge_sign = 0;
    
    double downstream_position_z = 0.;

    for ( auto emulsion : group ) {
      TVector3 position = emulsion->GetAbsolutePosition().GetValue();
      CalcPosInEccCoordinate(position, 4);
      TVector3 tangent = emulsion->GetTangent().GetValue();
      tangent = (1./tangent.Z()) * tangent;

      mom_basetrack.pl = emulsion->GetPlate() + 1;
      mom_basetrack.rawid = emulsion->GetEmulsionTrackId();
      mom_basetrack.x = position.X();
      mom_basetrack.y = position.Y();
      mom_basetrack.z = position.Z();
      mom_basetrack.ax = tangent.X();
      mom_basetrack.ay = tangent.Y();
      mom_basetrack.m[0].zone = 0;
      mom_basetrack.m[0].view = 0;
      mom_basetrack.m[0].imager = 0;
      mom_basetrack.m[0].ph = std::min((int)((emulsion->GetEdepSum() + emulsion->GetEdepDiff()) * 1000 / 2), 9999);
      mom_basetrack.m[0].pixelnum = 0;
      mom_basetrack.m[0].hitnum = 0;
      mom_basetrack.m[1].zone = 0;
      mom_basetrack.m[1].view = 0;
      mom_basetrack.m[1].imager = 0;
      mom_basetrack.m[1].ph = std::min((int)((emulsion->GetEdepSum() - emulsion->GetEdepDiff()) * 1000 / 2), 9999);
      mom_basetrack.m[1].pixelnum = 0;
      mom_basetrack.m[1].hitnum = 0;

      mom_chain.base.push_back(mom_basetrack);

      mom_basetrack_pair.first = mom_basetrack_pair.second;
      mom_basetrack_pair.first.z = 0.;
      mom_basetrack_pair.second.pl = emulsion->GetPlate() + 1;
      mom_basetrack_pair.second.rawid = emulsion->GetEmulsionTrackId();
      mom_basetrack_pair.second.x = position.X();
      mom_basetrack_pair.second.y = position.Y();
      mom_basetrack_pair.second.z = position.Z() - downstream_position_z;
      mom_basetrack_pair.second.ax = tangent.X();
      mom_basetrack_pair.second.ay = tangent.Y();
      if ( emulsion != group.front() )
	mom_chain.base_pair.push_back(mom_basetrack_pair);

      downstream_position_z = position.Z();
      
    }

    ev.chains.push_back(mom_chain);
  } // chains

  return;
}


bool ConnectionFunction::JudgeConnect(B2EmulsionSummary* down, B2EmulsionSummary* up, t2l_param param) const {
  return JudgeConnectXY(down, up, param) && JudgeConnectRL(down, up, param);
}

// down : alignementの基準, up : alignementで補正されうる
bool ConnectionFunction::JudgeConnectXY(B2EmulsionSummary* down, B2EmulsionSummary* up, t2l_param param) const {

  double d_pos_x, d_pos_y, d_ang_x, d_ang_y;
  double all_pos_x, all_pos_y, all_ang_x, all_ang_y;

  TVector3 tangent_down  = down->GetTangent().GetValue();
  TVector3 tangent_up    = up->GetTangent().GetValue();
  tangent_down = (1. / tangent_down.Z()) * tangent_down;
  tangent_up = (1./ tangent_up.Z()) * tangent_up;
  TVector3 position_down = down->GetFilmPosition().GetValue();
  TVector3 position_up   = up->GetFilmPositionInDownCoordinate().GetValue();
  position_down.SetZ(down->GetAbsolutePosition().GetValue().Z());
  position_up.SetZ(up->GetAbsolutePosition().GetValue().Z());
  position_down = 1000. * position_down; // mm -> um
  position_up = 1000. * position_up;

  all_ang_x = param.intercept_ax
    + param.slope_ax * std::fabs(tangent_down.X())
    + param.slope2_ax * std::pow(tangent_down.X(), 2.);
  all_ang_y = param.intercept_ay
    + param.slope_ay * std::fabs(tangent_down.Y())
    + param.slope2_ay * std::pow(tangent_down.Y(), 2.);
  all_pos_x = param.intercept_px
    + param.slope_px * std::fabs(tangent_down.X())
    + param.slope2_px * std::pow(tangent_down.X(), 2.);
  all_pos_y = param.intercept_py
    + param.slope_py * std::fabs(tangent_down.Y())
    + param.slope2_py * std::pow(tangent_down.Y(), 2.);

  d_pos_x = position_up.X()
    - position_down.X()
    - (tangent_down.X() + tangent_up.X()) / 2. * (position_up.Z() - position_down.Z());
  d_pos_y= position_up.Y()
    - position_down.Y()
    - (tangent_down.Y() + tangent_up.Y()) / 2. * (position_up.Z() - position_down.Z());
  d_ang_x = tangent_up.X() - tangent_down.X();
  d_ang_y = tangent_up.Y() - tangent_down.Y();
  if ( std::fabs(d_ang_x) > std::fabs(all_ang_x) ) return false;
  if ( std::fabs(d_ang_y) > std::fabs(all_ang_y) ) return false;

  if ( std::fabs(d_pos_x) > std::fabs(all_pos_x) ) return false;
  if ( std::fabs(d_pos_y) > std::fabs(all_pos_y) ) return false;

  return true;

}

bool ConnectionFunction::JudgeConnectRL(B2EmulsionSummary* down, B2EmulsionSummary* up, t2l_param param) const {

  double angle, d_pos_r, d_pos_l, d_ang_r, d_ang_l;
  double all_pos_r, all_pos_l, all_ang_r, all_ang_l;

  TVector3 tangent_down  = down->GetTangent().GetValue();
  TVector3 tangent_up    = up->GetTangentInDownCoordinate().GetValue();
  tangent_down = (1. / tangent_down.Z()) * tangent_down;
  tangent_up = (1. / tangent_up.Z()) * tangent_up;
  TVector3 position_down = down->GetFilmPosition().GetValue();
  TVector3 position_up   = up->GetFilmPositionInDownCoordinate().GetValue();
  position_down.SetZ(down->GetAbsolutePosition().GetValue().Z());
  position_up.SetZ(up->GetAbsolutePosition().GetValue().Z());
  position_down = 1000. * position_down; // mm -> um
  position_up = 1000. * position_up;
  
  angle = std::sqrt(tangent_down.X() * tangent_down.X() + tangent_down.Y() * tangent_down.Y());
  if ( angle < 0.01 ) return true;

  all_ang_r = param.intercept_ar
    + param.slope_ar * angle
    + param.slope2_ar * angle * angle;
  all_ang_l = param.intercept_al
    + param.slope_al * angle
    + param.slope2_al * angle * angle;
  all_pos_r = param.intercept_pr
    + param.slope_pr * angle
    + param.slope2_pr * angle * angle;
  all_pos_l = param.intercept_pl
    + param.slope_pl * angle
    + param.slope2_pl * angle * angle;

  d_ang_r = (tangent_up.X() - tangent_down.X()) * tangent_down.X()
    + (tangent_up.Y() - tangent_down.Y()) * tangent_down.Y();
  d_ang_l = (tangent_up.X() - tangent_down.X()) * tangent_down.Y()
    - (tangent_up.Y() - tangent_down.Y()) * tangent_down.X();

  if ( std::fabs(d_ang_r) > std::fabs(all_ang_r) ) return false;
  if ( std::fabs(d_ang_l) > std::fabs(all_ang_l) ) return false;

  CalculatePositionDifference(down, up, d_pos_r, d_pos_l);

  if ( std::fabs(d_pos_r) > std::fabs(all_pos_r) ) return false;
  if ( std::fabs(d_pos_l) > std::fabs(all_pos_l) ) return false;

  return true;

}

void ConnectionFunction::CalculatePositionDifference(B2EmulsionSummary* down, B2EmulsionSummary* up,
						     double &dr, double &dl) const {

  TVector3 pos0 = down->GetFilmPosition().GetValue();
  TVector3 pos1 = up->GetFilmPositionInDownCoordinate().GetValue();
  pos0.SetZ(down->GetAbsolutePosition().GetValue().Z());
  pos1.SetZ(up->GetAbsolutePosition().GetValue().Z());
  TVector3 dir0 = down->GetTangent().GetValue();
  TVector3 dir1 = up->GetTangent().GetValue();
  pos0 = 1000. * pos0; // mm -> um
  pos1 = 1000. * pos1;
  dir0 = (1. / dir0.Z()) * dir0;
  dir1 = (1. / dir1.Z()) * dir1;

  // 外挿基準点を1:1に内分した点に設定
  TVector3 base_point = 0.5 * (pos0 + pos1);
  TVector3 difference = pos1 - pos0;

  Double_t ratio0 = -1 * ((pos0 - base_point) * difference) / (dir0 * difference);
  Double_t ratio1 = -1 * ((pos1 - base_point) * difference) / (dir1 * difference);
  
  TVector3 extrapolate0 = pos0 + ratio0 * dir0;
  TVector3 extrapolate1 = pos1 + ratio1 * dir1;

  TVector3 unit_r, unit_l;
  unit_l.SetXYZ(-1 + difference.Y(), difference.X(), 0.);
  unit_r.SetXYZ(-1. * difference.X() * difference.Z(),
		-1. * difference.Y() * difference.Z(),
		difference.X() * difference.X() + difference.Y() * difference.Y());
  unit_l = unit_l.Unit();
  unit_r = unit_r.Unit();

  dr = (extrapolate1 - extrapolate0) * unit_r;
  dl = (extrapolate1 - extrapolate0) * unit_l;
}

bool ConnectionFunction::JudgeFiducialArea(std::vector<FiducialArea> &area, B2EmulsionSummary *emulsion) const {
  
  std::map<double, Point> point_map;
  TVector3 b_pos = emulsion->GetFilmPosition().GetValue();
  b_pos = 1000. * b_pos; // mm -> um
  TVector3 b_ang = emulsion->GetTangent().GetValue();
  b_ang = (1. / b_ang.Z()) * b_ang;
  /*
  double ex_x, ex_y, distance;
  for ( auto itr = area.begin(); itr != area.end(); itr++ ) {
    ex_x = b_pos.X() + b_ang.X() * (itr->p[0].z - b_pos.Z());
    ex_y = b_pos.Y() + b_ang.Y() * (itr->p[0].z - b_pos.Z());
    distance = std::pow(ex_x - itr->p[0].x, 2.) + std::pow(ex_y - itr->p[0].y, 2.);
    point_map.insert(std::make_pair(distance, itr->p[0]));
  }
  // 外挿先から距離が一番近い点のz座標を使用
  double z = point_map.begin()->second.z;
  double x = b_pos.X() + b_ang.X() * (z - b_pos.Z());
  double y = b_pos.Y() + b_ang.Y() * (z - b_pos.Z());
  */
  double x = b_pos.X();
  double y = b_pos.Y();
  // (x, y) からx軸正の方向に直線を引き，その直線と多角形の辺が何回交わるか
  // 下から上に交わったとき+1, 上から下に交わったとき-1
  int wn = 0;
  double vt;
  for ( auto itr = area.begin(); itr != area.end(); itr++ ) {
    // 上向きの辺，点Pがy軸方向について始点と終点の間にある，ただし終点は含まない
    if ( itr->p[0].y <= y && itr->p[1].y > y ) {
      // 辺は点Pより右側にある，ただし重ならない
      // 辺が点Pと同じ高さになる位置を特定し，そのときのxの値を点Pのx座標と比較する
      vt = (y - itr->p[0].y) / (itr->p[1].y - itr->p[0].y);
      if ( x < itr->p[0].x + vt * (itr->p[1].x - itr->p[0].x)) {
	++wn; // 上向きの辺と交差
      }
    }
    // 下向きの辺，点Pがy軸方向について始点と終点の間にある，ただし始点は含まない
    else if ( itr->p[0].y > y && itr->p[1].y <= y ) {
      // 辺は点Pよりも右側にある．ただし重ならない
      // 辺が点Pと同じ高さになる位置を特定し，そのときのxの値を点Pのx座標と比較する      
      vt = (y - itr->p[0].y) / (itr->p[1].y - itr->p[0].y);
      if ( x < itr->p[0].x + vt * (itr->p[1].x - itr->p[0].x)) {
	--wn; // 下向きの辺と交差
      }
    }
  }

  if ( wn >= 1 ) return true;
  return false;

}
/*
void ConnectionFunction::PartnerSearch(B2EmulsionSummary* target,
				       std::vector<B2EmulsionSummary* > &emulsions_partner,
				       std::vector<B2EmulsionSummary* > &emulsions,
				       partner_param &param) {
  for ( auto emulsion : emulsions ) {
    if ( emulsion->GetEmulsionTrackId() == target->GetEmulsionTrackId() ) continue;
    if ( emulsion->GetPlate() != target->GetPlate() &&
	 emulsion->GetPlate() != target->GetPlate() +1 ) continue;
    // MD & OA が入るか
    
    emulsions_partner.push_back(emulsion);
  }

  return;
}
*/
