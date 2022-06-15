#include <boost/log/trivial.hpp>
#include <boost/log/core.hpp>
#include <boost/log/expressions.hpp>
#include <boost/unordered_set.hpp>
#include <boost/unordered_map.hpp>

#include <vector>
#include <array>
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
    connection_data_.GetEfficiencyData(ecc, ecc_efficiency_[ecc]);
  }

  gRandom->SetSeed(time(NULL));

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
						  int ecc_id,
						  int material_id) const {

  Momentum_recon::Mom_chain mom_chain;
  Momentum_recon::Mom_basetrack mom_basetrack;
  std::pair<Momentum_recon::Mom_basetrack, Momentum_recon::Mom_basetrack> mom_basetrack_pair;

  for ( auto chain : chains ) {
    int num_base = chain.size();
    int num_link = num_base - 1;

    if ( material_id == B2Material::kWater &&
	 num_base < 2 ) continue;
    else if ( material_id == B2Material::kIron &&
	      num_base < 3 ) continue;

    mom_chain.base.clear();
    mom_chain.base_pair.clear();
    mom_chain.base.reserve(num_base);
    mom_chain.base_pair.reserve(num_link);

    mom_chain.chainid = chain.front()->GetParentTrackId();

    mom_chain.particle_flag = chain.front()->GetParentTrack().GetParticlePdg();
    if ( chain.front()->GetTangent().GetValue().Z() > 0 )
      mom_chain.direction = 1;
    else
      mom_chain.direction = -1;

    // stop flag
    TVector3 stop_position;
    double stop_mom = 1000.;
    if ( mom_chain.direction == 1 ) {
      stop_position = chain.front()->GetParentTrack().GetFinalPosition().GetValue();
      stop_mom = chain.front()->GetMomentum().GetValue().Mag();
    }
    else if ( mom_chain.direction == -1 ) {
      stop_position = chain.back()->GetParentTrack().GetFinalPosition().GetValue();
      stop_mom = chain.back()->GetMomentum().GetValue().Mag();
    }

    CalcPosInEccCoordinate(stop_position, ecc_id);
    if ( 0. < stop_position.X() &&
	 stop_position.X() < 250.e3 &&
	 0. < stop_position.Y() &&
	 stop_position.Y() < 250.e3 &&
	 -228.e3 < stop_position.Z() &&
	 stop_position.Z() < 0. ) {
      mom_chain.stop_flag = 2;
    }
    else
      mom_chain.stop_flag = 0;
    
    if ( mom_chain.particle_flag == 2212 ) {
      if ( stop_mom > 200. )
	mom_chain.stop_flag = 0;
    }
    else if ( stop_mom > 50.) {
      mom_chain.stop_flag = 0;
    }

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
		+ NINJA_ENV_THICK);

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
    TVector3 tangent_vec = emulsion->GetTangent().GetValue();
    double tangent = std::hypot(tangent_vec.X(), tangent_vec.Y());
    if ( tangent >= 4.0 ) continue;

    int bin_id = GetAngleBinId(tangent);
    double efficiency = ecc_efficiency_[ecc_id].at(emulsion->GetPlate() + 1).at(bin_id).efficiency;

    if ( emulsion->GetParentTrack().GetParticlePdg() == PDG_t::kProton ||
	 (B2Pdg::IsChargedPion(emulsion->GetParentTrack().GetParticlePdg()) &&
	  emulsion->GetParentTrack().GetFinalAbsoluteMomentum().GetValue() < 100.) ) {
      emulsions_detected.push_back(emulsion);
    }
    else if ( gRandom->Uniform() < efficiency ) {
      emulsions_detected.push_back(emulsion);
    }
    else {
      BOOST_LOG_TRIVIAL(trace) << "Emulsion : " << emulsion->GetEmulsionTrackId()
			       << " is not detected : "
			       << " PL" << emulsion->GetPlate() + 1;
    }
  }

  return;

}

int ConnectionFunction::GetAngleBinId(double tangent) const {

  tangent = std::fabs(tangent);

  if ( tangent < 0.1 ) return 0;
  else if ( tangent < 0.3 ) return 1;
  else if ( tangent < 0.5 ) return 2;
  else if ( tangent < 0.7 ) return 3;
  else if ( tangent < 0.9 ) return 4;
  else if ( tangent < 1.1 ) return 5;
  else if ( tangent < 1.3 ) return 6;
  else if ( tangent < 1.5 ) return 7;
  else if ( tangent < 2.0 ) return 8;
  else if ( tangent < 2.5 ) return 9;
  else if ( tangent < 3.0 ) return 10;
  else if ( tangent < 3.5 ) return 11;
  else if ( tangent < 4.0 ) return 12;
  else throw std::runtime_error("Large angle efficiency cannot be calculated");

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


void ConnectionFunction::ModifyGroups(std::vector<std::pair<B2EmulsionSummary*, std::vector<std::pair<B2EmulsionSummary*, B2EmulsionSummary* > > > > &groups_modified,
				      std::vector<std::pair<B2EmulsionSummary*, std::vector<std::pair<B2EmulsionSummary*, B2EmulsionSummary* > > > > &groups_reconnected,
				      std::vector<B2EmulsionSummary* > &emulsions) const {
  for ( auto group : groups_reconnected ) {
    std::pair<B2EmulsionSummary*, std::vector<std::pair<B2EmulsionSummary*, B2EmulsionSummary* > > > group_modified;
    ModifyGroup(group_modified, group, emulsions);
    if ( IsValidGroup(group_modified) )
      groups_modified.push_back(group_modified);
  }

  return;
}

void ConnectionFunction::ModifyGroup(std::pair<B2EmulsionSummary*, std::vector<std::pair<B2EmulsionSummary*, B2EmulsionSummary* > > > &group_modified,
				     std::pair<B2EmulsionSummary*, std::vector<std::pair<B2EmulsionSummary*, B2EmulsionSummary* > > > &group,
				     std::vector<B2EmulsionSummary* > &emulsions) const {

  group_modified.first = group.first;

  // basetrack の set にする
  std::set<Segment > segment_set;
  for ( auto linklet_b2 : group.second ) {
    std::pair<Segment, Segment > linklet_seg;
    LinkletConvertB2ToSeg(linklet_b2, linklet_seg);
    segment_set.insert(linklet_seg.first);
    segment_set.insert(linklet_seg.second);
  }

  // plate ごとに multimap につめる
  std::multimap<int, B2EmulsionSummary* > b2_multimap;
  double ax_ave = 0.;
  double ay_ave = 0.;
  double ax_dev = 0.;
  double ay_dev = 0.;
  for ( auto itr = segment_set.begin(); itr != segment_set.end(); itr++ ) {
    B2EmulsionSummary* tmp;
    ConvertRawIdToB2((*itr).rawid, tmp, emulsions);    
    b2_multimap.insert(std::make_pair((*itr).plate, tmp));
    TVector3 tangent = tmp->GetTangent().GetValue();
    tangent = (1. / tangent.Z()) * tangent;
    ax_ave += tangent.X(); ay_ave += tangent.Z();
    ax_dev += tangent.X() * tangent.X();
    ay_dev += tangent.Y() * tangent.Y();
  }

  ax_ave /= segment_set.size(); ay_ave /= segment_set.size();
  ax_dev /= segment_set.size(); ay_dev /= segment_set.size();
  ax_dev = std::sqrt(ax_dev - ax_ave * ax_ave);
  ay_dev = std::sqrt(ay_dev - ay_ave * ay_ave);

  // 各 plate について 1 basetrack までにする
  // 取り除かれた basetrack の list を生成
  std::vector<int> removed_tracks;
  for ( int iplate = 0; iplate < 134; iplate++ ) {
    if ( b2_multimap.count(iplate) < 2 ) continue;
    auto range = b2_multimap.equal_range(iplate);
    int ncount = 0;
    double smallest_angle_diff = 0.;
    int smallest_rawid;
    for ( auto itr = range.first; itr != range.second; itr++ ) {
      B2EmulsionSummary* tmp = (*itr).second;
      TVector3 tangent = tmp->GetTangent().GetValue();
      tangent = (1. / tangent.Z()) * tangent;
      double angle_diff_tmp = std::hypot((tangent.X() - ax_ave) / ax_dev,
					 (tangent.Y() - ay_ave) / ay_dev);
      if ( ncount == 0 ) {
	smallest_rawid = tmp->GetEmulsionTrackId();
	smallest_angle_diff = angle_diff_tmp;
      }
      else {
	if ( angle_diff_tmp < smallest_angle_diff ) {
	  removed_tracks.push_back(smallest_rawid);
	  smallest_rawid = tmp->GetEmulsionTrackId();
	  smallest_angle_diff = angle_diff_tmp;
	}
	else {
	  removed_tracks.push_back(tmp->GetEmulsionTrackId());
	}
      }
      ncount++;
    }
  }

  // 取り除くべき basetrack を含む linklet を取り除きながら，group_modified に push back
  for ( auto linklet_b2 : group.second ) {
    bool is_removed = false;
    for ( auto removed_track : removed_tracks ) {      
      if ( linklet_b2.first->GetEmulsionTrackId() == removed_track ||
	   linklet_b2.second->GetEmulsionTrackId() == removed_track ) {
	is_removed = true;
      }
    }
    if ( is_removed ) continue;
    group_modified.second.push_back(linklet_b2);
  }

  return;

}

bool ConnectionFunction::IsValidGroup(std::pair<B2EmulsionSummary*, std::vector<std::pair<B2EmulsionSummary*, B2EmulsionSummary* > > > &group) const {
  for ( auto linklet : group.second ) {
    if ( linklet.first->GetEmulsionTrackId() == group.first->GetEmulsionTrackId() ||
	 linklet.second->GetEmulsionTrackId() == group.first->GetEmulsionTrackId() ) return true;
  }
  return false;
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
					 std::pair<B2EmulsionSummary*, std::vector<std::pair<B2EmulsionSummary*, B2EmulsionSummary* > > > &muon_group,
					 B2EmulsionSummary* &vertex_track,
					 std::vector<B2EmulsionSummary* > &emulsions) const {
  bool find_flag = false;
  int vertex_plate = -1;

  std::vector<Group > groups_seg;

  // まずは muon を探す
  for ( auto group : groups ) {
    if ( B2Pdg::IsMuonPlusOrMinus(group.first->GetParentTrack().GetParticlePdg()) &&
	 !group.second.empty() ) {
      Group group_seg;
      GroupConvertB2ToSeg(group, group_seg);
      if ( group_seg.start_plate == 3 ||
	   group_seg.start_plate == 4 )
	groups_seg.push_back(group_seg);
    }
  }

  // muon が見つかっていなかった際は pion を探しておく
  if ( groups_seg.empty() ) {
    for ( auto group : groups ) {
      if ( B2Pdg::IsChargedPion(group.first->GetParentTrack().GetParticlePdg()) &&
	   !group.second.empty()) {
	Group group_seg;
	GroupConvertB2ToSeg(group, group_seg);
	if ( group_seg.start_plate == 3 ||
	     group_seg.start_plate == 4 )
	  groups_seg.push_back(group_seg);
      }
    }
  }

  std::vector<Segment > mu_start_seg_vec;
  std::vector<Group > mu_group_seg_vec;
  boost::unordered_multimap<Segment, Segment > ups;
  boost::unordered_multimap<Segment, Segment > downs;
  boost::unordered_multimap<Segment, Segment > ups_start;
  boost::unordered_multimap<Segment, Segment > downs_start;
  for ( auto group_seg : groups_seg ) {
    Segment start_seg;
    start_seg.plate = group_seg.start_plate;
    start_seg.rawid = group_seg.start_rawid;
    mu_start_seg_vec.push_back(start_seg);
    mu_group_seg_vec.push_back(group_seg);
    CollectEdgeTracks(ups, downs, ups_start, downs_start,
		      group_seg, start_seg);
  }

  if ( mu_start_seg_vec.size() == 1 ||
       mu_start_seg_vec.size() == 2 ) {
    BOOST_LOG_TRIVIAL(debug) << "Muon is detected : size : " << mu_start_seg_vec.size();
    GroupConvertSegToB2(mu_group_seg_vec.at(0), muon_group, emulsions);
    if ( mu_start_seg_vec.size() == 2 ) {
      if ( mu_group_seg_vec.at(0) != mu_group_seg_vec.at(1)) {    
	BOOST_LOG_TRIVIAL(warning) << "Groups not identical";
      }
    }

    auto range = ups.equal_range(mu_start_seg_vec.front());
    for ( auto itr = range.first; itr != range.second; itr++ ) {     
      Segment upstream_seg = itr->second;
      if ( upstream_seg.plate > vertex_plate ) { // 一番深いところを vertex_plate とする
	ConvertRawIdToB2(upstream_seg.rawid, vertex_track, emulsions);
	vertex_plate = vertex_track->GetPlate() + 1;
      }
      find_flag = true;
    }
  }

  return find_flag;

}

void ConnectionFunction::PartnerSearch(B2EmulsionSummary* vertex_track,
				       std::vector<std::pair<B2EmulsionSummary*, std::vector<std::pair<B2EmulsionSummary*, B2EmulsionSummary* > > > > &groups,
				       std::vector<B2EmulsionSummary* > &emulsions_partner,
				       std::vector<B2EmulsionSummary* > &emulsions,
				       int ecc_id,
				       TVector3 &recon_vertex) const {
  int num_partner = 0;
  TVector3 recon_vertex_tmp(0., 0., 0.);

  int vertex_plate = vertex_track->GetPlate() + 1;
  for ( auto emulsion : emulsions ) {
    if ( emulsion == vertex_track ) continue;
    if ( emulsion->GetPlate() + 1 == vertex_plate ||
	 emulsion->GetPlate() + 1 == vertex_plate + 1 ) {
      if ( JudgePartnerTrack(vertex_track, emulsion, groups, ecc_id, recon_vertex_tmp) ) {
	num_partner++;
	recon_vertex += recon_vertex_tmp;
	emulsions_partner.push_back(emulsion);
      }      
    }
  }
  
  if ( num_partner == 0 ) { // one prong vertex
    TVector3 position = vertex_track->GetAbsolutePosition().GetValue();
    CalcPosInEccCoordinate(position, ecc_id);
    TVector3 tangent = vertex_track->GetTangent().GetValue();
    if ( vertex_plate < 14 ) {
      recon_vertex = position - 1000. * (NINJA_EMULSION_LAYER_THICK + 0.5 * NINJA_IRON_LAYER_THICK) * tangent;
    } else if ( vertex_plate == 14 ) {
      recon_vertex = position - 1000. * (NINJA_EMULSION_LAYER_THICK + NINJA_ENV_THICK) * tangent;
    } else if ( vertex_plate%2 == 0 ) {
      recon_vertex = position - 1000. * (NINJA_EMULSION_LAYER_THICK + NINJA_ENV_THICK + 0.5 * NINJA_WATER_LAYER_THICK) * tangent;
    } else if ( vertex_plate%2 == 1 ) {
      recon_vertex = position - 1000. * (NINJA_EMULSION_LAYER_THICK + NINJA_IRON_LAYER_THICK) * tangent;
    }
  }
  else {
    recon_vertex = (1. / num_partner) * recon_vertex;
  }

  return;

}

void ConnectionFunction::SelectPartnerGroups(std::vector<std::pair<B2EmulsionSummary*, std::vector<std::pair<B2EmulsionSummary*, B2EmulsionSummary* > > > > &groups,
					     std::vector<std::pair<B2EmulsionSummary*, std::vector<std::pair<B2EmulsionSummary*, B2EmulsionSummary* > > > > &groups_partner,
					     std::vector<B2EmulsionSummary* > &emulsions_partner,
					     std::vector<B2EmulsionSummary* > &emulsions) const {

  // groups の中から emulsions_partner に対応する start segment をもつ Group をみつけて
  // groups_partner に push_back するだけ
  std::vector<Segment > partner_seg_vec;
  for ( auto emulsion : emulsions_partner ) {
    Segment tmp;
    tmp.plate = emulsion->GetPlate() + 1;
    tmp.rawid = emulsion->GetEmulsionTrackId();
    partner_seg_vec.push_back(tmp);
  }

  std::map<Segment, Group > groups_seg_map;
  for ( auto group : groups ) {
    Group group_seg;
    GroupConvertB2ToSeg(group, group_seg);
    Segment start_seg;
    start_seg.plate = group_seg.start_plate;
    start_seg.rawid = group_seg.start_rawid;
    groups_seg_map.insert(std::make_pair(start_seg, group_seg));
  }

  for ( auto partner_seg : partner_seg_vec ) {
    if ( groups_seg_map.find(partner_seg) == groups_seg_map.end() ) continue;
    Group group_partner_seg = groups_seg_map.at(partner_seg);
    std::pair<B2EmulsionSummary*, std::vector<std::pair<B2EmulsionSummary*, B2EmulsionSummary* > > > group_partner;
    GroupConvertSegToB2(group_partner_seg, group_partner, emulsions);
    groups_partner.push_back(group_partner);
  }

  return;

}

void ConnectionFunction::AddGroupsToEventInfo(Momentum_recon::Event_information &ev,
					      B2EmulsionSummary* vertex_track,
					      std::pair<B2EmulsionSummary*, std::vector<std::pair<B2EmulsionSummary*, B2EmulsionSummary* > > > &muon_group,
					      std::vector<std::pair<B2EmulsionSummary*, std::vector<std::pair<B2EmulsionSummary*, B2EmulsionSummary* > > > > &groups_partner,
					      std::vector<B2EmulsionSummary* > &emulsions,
					      int ecc_id) const {
  AddGroupToEventInfo(ev, vertex_track, muon_group, muon_group, emulsions, ecc_id); // muon
  
  for ( auto group_partner : groups_partner ) {
    AddGroupToEventInfo(ev, vertex_track, muon_group, group_partner, emulsions, ecc_id);
  }
  return;
  
}

void ConnectionFunction::AddGroupToEventInfo(Momentum_recon::Event_information &ev,
					     B2EmulsionSummary* vertex_track,
					     std::pair<B2EmulsionSummary*, std::vector<std::pair<B2EmulsionSummary*, B2EmulsionSummary* > > > &muon_group,
					     std::pair<B2EmulsionSummary*, std::vector<std::pair<B2EmulsionSummary*, B2EmulsionSummary* > > > &group,
					     std::vector<B2EmulsionSummary* > &emulsions,
					     int ecc_id) const {

  Momentum_recon::Mom_chain mom_chain;
  Momentum_recon::Mom_basetrack mom_basetrack;
  std::pair<Momentum_recon::Mom_basetrack, Momentum_recon::Mom_basetrack > mom_basetrack_pair;

  std::set<Segment > segs_set;
  for ( auto linklet : group.second ) {
    std::pair<Segment, Segment > linklet_seg;
    LinkletConvertB2ToSeg(linklet, linklet_seg);
    segs_set.insert(linklet_seg.first);
    segs_set.insert(linklet_seg.second);
  }

  std::vector<B2EmulsionSummary* > chain;
  for ( auto itr = segs_set.begin(); itr != segs_set.end(); itr++ ) {
    B2EmulsionSummary* tmp;
    ConvertRawIdToB2((*itr).rawid, tmp, emulsions);
    chain.push_back(tmp);
  }

  int num_base = chain.size();
  int num_link = num_base - 1;
  if ( num_base < 2 ) return;
  mom_chain.base.clear();
  mom_chain.base_pair.clear();
  mom_chain.base.reserve(num_base);
  mom_chain.base_pair.reserve(num_link);

  mom_chain.chainid = chain.front()->GetParentTrackId();
  mom_chain.stop_flag = -1;
  mom_chain.particle_flag = chain.front()->GetParentTrack().GetParticlePdg();
  mom_chain.particle_flag = std::fabs(mom_chain.particle_flag) * 10000;

  if ( muon_group.first->GetEmulsionTrackId() == group.first->GetEmulsionTrackId() ) { // muon
    mom_chain.direction = 1;
    mom_chain.particle_flag += 13;
  }
  else if ( vertex_track->GetPlate() == group.first->GetPlate() ) { // forward partner
    mom_chain.direction = 1;
    if ( group.first->GetPlate() < 16 ) { // iron ECC
      if ( group.second.back().second->GetPlate() - group.second.front().first->GetPlate() < 1 ) return;
    }
    else if ( group.first->GetPlate() % 2 == 0 ) { // water
      if ( group.second.back().second->GetPlate() - group.second.front().first->GetPlate() < 1 ) return;
    }
    else { // iron
      if ( group.second.back().second->GetPlate() - group.second.front().first->GetPlate() < 2 ) return;
    }
  }
  else if ( vertex_track->GetPlate() + 1 == group.first->GetPlate() ) { // backward
    mom_chain.direction = -1;
    if ( group.first->GetPlate() < 16 ) { // iron ECC
      if ( group.second.back().second->GetPlate() - group.second.front().first->GetPlate() < 1 ) return;
    }
    else if ( group.first->GetPlate() % 2 == 1 ) { // water
      if ( group.second.back().second->GetPlate() - group.second.front().first->GetPlate() < 1 ) return;
    }
    else { // iron
      if ( group.second.back().second->GetPlate() - group.second.front().first->GetPlate() < 2 ) return;
    }

  }
  mom_chain.charge_sign = 0;

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
    mom_basetrack.m[0].ph = 0;
    mom_basetrack.m[0].pixelnum = 0;
    mom_basetrack.m[0].hitnum = 0;
    mom_basetrack.m[1].zone = 0;
    mom_basetrack.m[1].view = 0;
    mom_basetrack.m[1].imager = 0;
    mom_basetrack.m[1].ph = 0;
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

  } // muon chain
  
  ev.chains.push_back(mom_chain);

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

bool ConnectionFunction::JudgeFiducialArea(const std::vector<FiducialArea> &area, B2EmulsionSummary *emulsion) const {
  
  std::map<double, Point> point_map;
  TVector3 b_pos = emulsion->GetFilmPosition().GetValue();
  b_pos = 1000. * b_pos; // mm -> um

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

bool ConnectionFunction::JudgeEdgeOut(B2EmulsionSummary* vertex_track, int ecc_id) const {
  int vertex_plate = vertex_track->GetPlate() + 1;
  TVector3 position = vertex_track->GetFilmPosition().GetValue();
  TVector3 tangent = vertex_track->GetTangent().GetValue();
  tangent = (1./tangent.Z()) * tangent;

  bool edge_out_flag = false;

  for ( int ipl = 1; ipl < 5; ipl++ ) {
    double z_difference = 0.;
    if ( vertex_plate < 12 ) { // 鉄ECC
      z_difference = -1. * (NINJA_FILM_THICK + NINJA_IRON_LAYER_THICK) * ipl;
    } else if ( vertex_plate == 12 ) {
      if ( ipl < 4 )
	z_difference = -1. * (NINJA_FILM_THICK + NINJA_IRON_LAYER_THICK) * ipl;
      else
	z_difference = -1. * (NINJA_FILM_THICK + NINJA_IRON_LAYER_THICK) * 3
	  - 2. * NINJA_ENV_THICK - NINJA_FILM_THICK;
    } else if ( vertex_plate == 13 ) { // 鉄-水ECC 接続部
      if ( ipl < 3 ) 
	z_difference = -1. * (NINJA_FILM_THICK + NINJA_IRON_LAYER_THICK) * ipl;
      else if ( ipl == 3 )
	z_difference = -1. * (NINJA_FILM_THICK + NINJA_IRON_LAYER_THICK) * 2
	  - 2. * NINJA_ENV_THICK - NINJA_FILM_THICK;
      else
	z_difference = -1. * (NINJA_FILM_THICK + NINJA_IRON_LAYER_THICK) * 2
	  - 2. * NINJA_ENV_THICK - NINJA_FILM_THICK - NINJA_IRON_LAYER_THICK - NINJA_FILM_THICK;
    } else if ( vertex_plate == 14 ) {
      if ( ipl == 1 )
	z_difference = -1. * (NINJA_FILM_THICK + NINJA_IRON_LAYER_THICK);
      else if ( ipl == 2 )
	z_difference = -1. * (NINJA_FILM_THICK + NINJA_IRON_LAYER_THICK)
	  - 2. * NINJA_ENV_THICK - NINJA_FILM_THICK;
      else if ( ipl == 3 )
	z_difference = -1. * (NINJA_FILM_THICK + NINJA_IRON_LAYER_THICK)
	  - 2. * NINJA_ENV_THICK - NINJA_FILM_THICK - NINJA_IRON_LAYER_THICK - NINJA_FILM_THICK;
      else
	z_difference= -1. * (NINJA_FILM_THICK + NINJA_IRON_LAYER_THICK)
	  - 2. * NINJA_ENV_THICK - NINJA_FILM_THICK
	  - 1. * (NINJA_FILM_THICK + NINJA_IRON_LAYER_THICK)
	  - 1. * (NINJA_FILM_THICK + NINJA_WATER_LAYER_THICK + 2. * NINJA_ENV_THICK);
    } else if ( vertex_plate == 15 ) {
      z_difference = -2. * NINJA_ENV_THICK - NINJA_FILM_THICK
	- 1. * (NINJA_FILM_THICK + NINJA_IRON_LAYER_THICK) * (ipl / 2)
	- 1. * (NINJA_FILM_THICK + NINJA_WATER_LAYER_THICK + 2. * NINJA_ENV_THICK) * ((ipl-1) / 2);
    } else if ( vertex_plate == 130 ) { // 最上流例外処理
      if ( ipl < 4 )
	z_difference =	- 1. * (NINJA_FILM_THICK + NINJA_IRON_LAYER_THICK) * ((ipl+1) / 2)
	  - 1. * (NINJA_FILM_THICK + NINJA_WATER_LAYER_THICK + 2. * NINJA_ENV_THICK) * (ipl / 2);
      else continue;       	
    } else if ( vertex_plate == 131 ) {
      if ( ipl < 3 )
	z_difference = - 1. * (NINJA_FILM_THICK + NINJA_IRON_LAYER_THICK) * (ipl / 2)
	  - 1. * (NINJA_FILM_THICK + NINJA_WATER_LAYER_THICK + 2. * NINJA_ENV_THICK) * ((ipl+1) / 2);
      else continue;
    } else if ( vertex_plate % 2 == 0 ) { // 鉄下流
      z_difference = - 1. * (NINJA_FILM_THICK + NINJA_IRON_LAYER_THICK) * ((ipl+1) / 2)
	- 1. * (NINJA_FILM_THICK + NINJA_WATER_LAYER_THICK + 2. * NINJA_ENV_THICK) * (ipl / 2);
    } else { // 水下流
      z_difference = - 1. * (NINJA_FILM_THICK + NINJA_IRON_LAYER_THICK) * (ipl / 2)
	- 1. * (NINJA_FILM_THICK + NINJA_WATER_LAYER_THICK + 2. * NINJA_ENV_THICK) * ((ipl+1) / 2);
    }

    double x = position.X() + tangent.X() * z_difference;
    double y = position.Y() + tangent.Y() * z_difference;
    B2EmulsionSummary* tmp = new B2EmulsionSummary();
    tmp->SetFilmPosition(TVector3(x, y, 0.));
    if ( !JudgeFiducialArea(ecc_fiducial_[ecc_id].at(vertex_plate + ipl), tmp) ) {
      BOOST_LOG_TRIVIAL(debug) << "Edge out : PL" << vertex_plate + ipl;
      edge_out_flag = true;
    }
  }
  
  return edge_out_flag;
  
}

bool ConnectionFunction::JudgePartnerTrack(B2EmulsionSummary* vertex_track,
					   B2EmulsionSummary* partner_track,
					   std::vector<std::pair<B2EmulsionSummary*, std::vector<std::pair<B2EmulsionSummary*, B2EmulsionSummary* > > > > &groups,
					   int ecc_id,
					   TVector3 &recon_vertex) const {
  TVector3 parent_pos = vertex_track->GetAbsolutePosition().GetValue();
  parent_pos = 1000. * parent_pos; // mm -> um
  TVector3 daughter_pos = partner_track->GetAbsolutePosition().GetValue();
  daughter_pos = 1000. * daughter_pos;
  TVector3 parent_dir = vertex_track->GetTangent().GetValue();
  TVector3 daughter_dir = partner_track->GetTangent().GetValue();

  int vertex_plate = 0;

  std::array<Double_t, 2 > z_range;

  int vertex_material = GetVertexTargetMaterial(vertex_track->GetPlate() + 1);
  if ( vertex_material == B2Material::kWater ) {
    z_range.at(0) = (- NINJA_WATER_LAYER_THICK - NINJA_FILM_THICK - 2 * NINJA_ENV_THICK) * 1.e3 - 1000.;
    z_range.at(1) = (+ NINJA_BASE_LAYER_THICK + NINJA_EMULSION_LAYER_THICK) * 1.e3 + 30.;
  }
  else if ( vertex_material == B2Material::kIron ) {
    z_range.at(0) = (- NINJA_IRON_LAYER_THICK - NINJA_FILM_THICK) * 1.e3 - 50.;
    z_range.at(1) = (+ NINJA_BASE_LAYER_THICK + NINJA_EMULSION_LAYER_THICK) * 1.e3 + 30.;
  }

  z_range.at(0) += parent_pos.Z();
  z_range.at(1) += parent_pos.Z();

  std::array<Double_t, 2 > extrapolate_z;

  Double_t minimum_distance = GetMinimumDistance(parent_pos, daughter_pos,
						 parent_dir, daughter_dir,
						 z_range, extrapolate_z,
						 ecc_id, recon_vertex);
  Double_t tangent = std::hypot(parent_dir.X(), parent_dir.Y());
  Double_t tangent_partner = std::hypot(daughter_dir.X(), daughter_dir.Y());
  if ( minimum_distance < std::sqrt(std::pow(extrapolate_z[0] * (0.04 * tangent + 0.04) + 5, 2) +
				    std::pow(extrapolate_z[1] * (0.04 * tangent_partner + 0.04) + 5, 2)) ) {

    // 対応する group の端を探して，start_seg が端かを確認する
    Segment start_seg;
    start_seg.plate = partner_track->GetPlate() + 1;
    start_seg.rawid = partner_track->GetEmulsionTrackId();
    Group group_seg;
    for ( auto group : groups ) {
      Segment tmp;
      tmp.plate = group.first->GetPlate() + 1;
      tmp.rawid = group.first->GetEmulsionTrackId();
      if ( tmp != start_seg ) continue;
      GroupConvertB2ToSeg(group, group_seg);
      break;
    }

    boost::unordered_multimap<Segment, Segment > ups;
    boost::unordered_multimap<Segment, Segment > downs;
    boost::unordered_multimap<Segment, Segment > ups_start;
    boost::unordered_multimap<Segment, Segment > downs_start;
    CollectEdgeTracks(ups, downs, ups_start, downs_start,
		      group_seg, start_seg);

    if ( vertex_track->GetPlate() == partner_track->GetPlate() ) { // 最上流確認
      for ( auto itr = ups.begin(); itr != ups.end(); itr++ ) {
	if ( (*itr).second.plate > start_seg.plate ) return false;
      }
    }
    else if ( vertex_track->GetPlate() + 1 == partner_track->GetPlate() ) { // 最下流確認
      for ( auto itr = downs.begin(); itr != downs.end(); itr++ ) {
	if ( (*itr).second.plate < start_seg.plate ) return false;
      }
    }    
    return true;
  }
  else {
    return false;
  }

}

double ConnectionFunction::GetMinimumDistance(TVector3 parent_pos, TVector3 daughter_pos,
					      TVector3 parent_dir, TVector3 daughter_dir,
					      std::array<Double_t,2> z_range, std::array<Double_t,2> &extrapolate_z,
					      int ecc_id, TVector3 &recon_vertex) const {

  std::array<Double_t,2> extrapolate_distance;
  TVector3 position_difference = daughter_pos - parent_pos;
  // Almost parallel
  if ( TMath::ACos((parent_dir * daughter_dir) / (parent_dir.Mag() * daughter_dir.Mag())) < 1.e-4) {
    extrapolate_distance.at(0) = (parent_pos.Z() + daughter_pos.Z()) / 2. - parent_pos.Z();
    extrapolate_distance.at(1) = (parent_pos.Z() + daughter_pos.Z()) / 2. - daughter_pos.Z();
  }
  else {
    Double_t delta = parent_dir.Mag2() * daughter_dir.Mag2() - (parent_dir * daughter_dir) * (parent_dir * daughter_dir);
    extrapolate_distance.at(0) = ( 1 * (position_difference * parent_dir) * daughter_dir.Mag2()
				   - (parent_dir * daughter_dir) * (position_difference * daughter_dir) ) / delta;
    extrapolate_distance.at(1) = (-1 * (position_difference * daughter_dir) * parent_dir.Mag2()
				   + (parent_dir * daughter_dir) * (position_difference * parent_dir) ) / delta;
  }
  // z_range.at(0) : small, z_range.at(1) : large
  if ( z_range.at(0) > z_range.at(1) ) {
    std::swap(z_range.at(0), z_range.at(1));
  }

  if ( parent_pos.Z() + extrapolate_distance.at(0) < z_range.at(0) ||
       daughter_pos.Z() + extrapolate_distance.at(1) * daughter_dir.Z() < z_range.at(0)) {
    extrapolate_distance.at(0) = z_range.at(0) - parent_pos.Z();
    extrapolate_distance.at(1) = z_range.at(0) - daughter_pos.Z();
  }
  else if ( parent_pos.Z() + extrapolate_distance.at(0) > z_range.at(1) ||
	    daughter_pos.Z() + extrapolate_distance.at(1) * daughter_dir.Z() > z_range.at(1)) {
    extrapolate_distance.at(0) = z_range.at(1) - parent_pos.Z();
    extrapolate_distance.at(1) = z_range.at(1) - daughter_pos.Z();
  }

  extrapolate_z.at(0) = extrapolate_distance.at(0);
  extrapolate_z.at(1) = extrapolate_distance.at(1);
  TVector3 calculate_parent_position = parent_pos + extrapolate_distance.at(0) * parent_dir;
  TVector3 calculate_daughter_position = daughter_pos + extrapolate_distance.at(1) * daughter_dir;
  TVector3 distance_vec = calculate_parent_position - calculate_daughter_position;
  recon_vertex = calculate_daughter_position + 0.5 * distance_vec;
  recon_vertex = 1.e-3 * recon_vertex;
  CalcPosInEccCoordinate(recon_vertex, ecc_id);
  return distance_vec.Mag();

}

int ConnectionFunction::GetVertexMaterial(TVector3 recon_vertex) const {
  double z_pos = recon_vertex.Z();
  z_pos = 1.e-3 * z_pos; // um -> mm

  if ( z_pos >= -15 * NINJA_FILM_THICK
       - 11 * NINJA_IRON_LAYER_THICK
       - NINJA_SS_AC_THICK
       - NINJA_ENV_THICK ) { // Iron ECC
    z_pos  = z_pos + 4 * NINJA_FILM_THICK + NINJA_SS_AC_THICK; // most downstream iron -> origin
    if ( z_pos >= NINJA_EMULSION_LAYER_THICK + NINJA_BASE_LAYER_THICK ) return 5; // emulsion
    else if ( z_pos >= NINJA_EMULSION_LAYER_THICK ) return B2Material::kCarbon; // base
    else if ( z_pos >= 0 ) return 5; // emulsion
    else {
      int film_id = (int)(-z_pos / (NINJA_FILM_THICK + NINJA_IRON_LAYER_THICK));
      z_pos = z_pos + film_id * (NINJA_FILM_THICK + NINJA_IRON_LAYER_THICK); // certain iron downstream -> origin
      if ( z_pos >= -NINJA_IRON_LAYER_THICK ) {
	return B2Material::kIron;
      }
      else if ( z_pos >= -NINJA_IRON_LAYER_THICK - NINJA_EMULSION_LAYER_THICK ) {
	return 5; // emulsion
      }
      else if ( z_pos >= -NINJA_IRON_LAYER_THICK - NINJA_EMULSION_LAYER_THICK - NINJA_BASE_LAYER_THICK ) {
	return B2Material::kCarbon; // base
      }
      else return 5; // emulsion
    }
  }
  else if ( z_pos >= -133 * NINJA_FILM_THICK
	    - 58 * NINJA_WATER_LAYER_THICK
	    - (59 * 2 + 1) * NINJA_ENV_THICK
	    - 70 * NINJA_IRON_LAYER_THICK
	    - NINJA_SS_AC_THICK ) { // Water ECC
    z_pos = z_pos + 16 * NINJA_FILM_THICK + NINJA_SS_AC_THICK
      + 11 * NINJA_IRON_LAYER_THICK + 2 * NINJA_ENV_THICK; // most downstream iron in Water ECC -> origin
    if ( z_pos >= NINJA_EMULSION_LAYER_THICK + NINJA_BASE_LAYER_THICK ) return 5; // emulsion
    else if ( z_pos >= NINJA_EMULSION_LAYER_THICK ) return B2Material::kCarbon; // base
    else if ( z_pos >= 0. ) return 5; // emulsion
    else {
      int unit_id = (int)(-z_pos / (2 * NINJA_FILM_THICK + NINJA_IRON_LAYER_THICK
				    + NINJA_WATER_LAYER_THICK + 2 * NINJA_ENV_THICK));
      z_pos = z_pos + unit_id * (2 * NINJA_FILM_THICK + NINJA_IRON_LAYER_THICK
				 + NINJA_WATER_LAYER_THICK + 2 * NINJA_ENV_THICK); // certain iron downstream -> origin
      if ( z_pos >= -NINJA_IRON_LAYER_THICK ) return B2Material::kIron;
      else if ( z_pos >= -NINJA_IRON_LAYER_THICK - NINJA_EMULSION_LAYER_THICK ) return 5; // emulsion
      else if ( z_pos >= -NINJA_IRON_LAYER_THICK - NINJA_EMULSION_LAYER_THICK
		- NINJA_BASE_LAYER_THICK ) return B2Material::kCarbon; // base
      else if ( z_pos >= -NINJA_IRON_LAYER_THICK - NINJA_FILM_THICK ) return 5; // emulsion
      else if ( z_pos >= -NINJA_IRON_LAYER_THICK - NINJA_FILM_THICK - NINJA_ENV_THICK ) return 6; // envelope
      else if ( z_pos >= -NINJA_IRON_LAYER_THICK - NINJA_FILM_THICK - NINJA_ENV_THICK
		- NINJA_WATER_LAYER_THICK ) return B2Material::kWater;
      else if ( z_pos >= NINJA_IRON_LAYER_THICK - NINJA_FILM_THICK - NINJA_ENV_THICK
		- NINJA_WATER_LAYER_THICK - NINJA_ENV_THICK ) return 6; // envelope
      else if ( z_pos >= NINJA_IRON_LAYER_THICK - NINJA_FILM_THICK - NINJA_WATER_LAYER_THICK
		- 2 * NINJA_ENV_THICK - NINJA_EMULSION_LAYER_THICK ) return 5; // emulsion
      else if ( z_pos >= NINJA_IRON_LAYER_THICK - NINJA_FILM_THICK - NINJA_WATER_LAYER_THICK
		- 2 * NINJA_ENV_THICK - NINJA_EMULSION_LAYER_THICK - NINJA_BASE_LAYER_THICK ) return B2Material::kCarbon; // base
      else return 5; // emulsion
    }
  }
  return -1;
}

int ConnectionFunction::GetVertexTargetMaterial(int vertex_plate) const {
  if ( vertex_plate < 17 )
    return B2Material::kIron;
  else if ( vertex_plate % 2 == 1 )
    return B2Material::kWater;
  else
    return B2Material::kIron;
}
