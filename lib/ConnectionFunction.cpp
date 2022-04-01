#include <boost/log/core.hpp>
#include <boost/log/trivial.hpp>
#include <boost/log/expressions.hpp>

#include <vector>

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
  
  BOOST_LOG_TRIVIAL(info) << "Connection functions are initilized";

}

void ConnectionFunction::GetTrueEmulsionTracks(std::vector<const B2EmulsionSummary* > &emulsions,
					       B2SpillSummary &spill_summary, int ecc_id) const {
  auto it_emulsion = spill_summary.BeginEmulsion();
  while ( const auto emulsion = it_emulsion.Next() ) {
    if ( emulsion->GetParentTrackId() == 0 ) continue;
    if ( emulsion->GetFilmType() != B2EmulsionType::kECC ) continue;
    if ( emulsion->GetEcc() != ecc_id ) continue;
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

  // mm -> um
  position.SetX(position.X() * 1.e3);
  position.SetY(position.Y() * 1.e3);
  position.SetZ(position.Z() * 1.e3);

}

void ConnectionFunction::SmearEmulsions(std::vector<B2EmulsionSummary* > &emulsions_smeared,
					std::vector<const B2EmulsionSummary* > &emulsions) const {
  for ( auto emulsion : emulsions ) {
    B2EmulsionSummary *emulsion_tmp = new B2EmulsionSummary();
    TVector3 position = emulsion->GetAbsolutePosition().GetValue();
    SmearPosition(position);

    TVector3 tangent  = emulsion->GetTangent().GetValue();
    double sign = tangent.Z();
    tangent = sign * tangent;
    SmearTangentVector(tangent);

    emulsion_tmp->SetEmulsionTrackId(emulsion->GetEmulsionTrackId());
    emulsion_tmp->SetParentTrack(emulsion->GetParentTrack());
    // emulsion_tmp->SetParentTrackId(emulsion->GetParentTrackId());
    emulsion_tmp->SetMomentum(emulsion->GetMomentum());
    emulsion_tmp->SetAbsolutePosition(position);
    emulsion_tmp->SetTangent((1./sign) * tangent);
    emulsion_tmp->SetTangentInDownCoordinate((1./sign) * tangent);
    emulsion_tmp->SetEdepSum(emulsion_tmp->GetEdepSum());
    emulsion_tmp->SetEdepDiff(emulsion_tmp->GetEdepDiff());
    emulsion_tmp->SetFilmType(B2EmulsionType::kECC);
    emulsion_tmp->SetEcc(emulsion->GetEcc());
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

void ConnectionFunction::GenerateLinklet(std::vector<std::pair<B2EmulsionSummary*, B2EmulsionSummary* > > &linklet,
					 std::vector<B2EmulsionSummary* > &emulsions) const {

  for ( auto emulsion : emulsions ) {
    // Black 判定
    // Black だった場合は同じ track 由来に限りつなぐ
    // std::cout << emulsion->GetParentTrack().GetParticlePdg() << std::endl;
    if ( B2Pdg::IsMuonPlusOrMinus(emulsion->GetParentTrack().GetParticlePdg()) )
      {}
    // plate を確認
    int plate = emulsion->GetPlate();

    if ( plate > 16 ) {

    }
    else if {

    }



  }

}

// down : alignementの基準, up : alignementで補正されうる
bool ConnectionFunction::JudgeConnectXY(B2EmulsionSummary* down, B2EmulsionSummary* up, t2l_param &param) {

  double d_pos_x, d_pos_y, d_ang_x, d_ang_y;
  double all_pos_x, all_pos_y, all_ang_x, all_ang_y;

  TVector3 tangent_down  = down->GetTangent().GetValue();
  TVector3 tangent_up    = up->GetTangent().GetValue();
  TVector3 position_down = down->GetFilmPosition().GetValue();
  TVector3 position_up   = up->GetFilmPositionInDownCoordinate().GetValue();

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
    + param.slope_px * std::fabs(tangent_down.Y())
    + param.slope2_px * std::pow(tangent_down.Y(), 2.);

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

bool ConnectionFunction::JudgeConnectRL(B2EmulsionSummary* down, B2EmulsionSummary* up, t2l_param &param) {

  double angle, d_pos_r, d_pos_l, d_ang_r, d_ang_l;
  double all_pos_r, all_pos_l, all_ang_r, all_ang_l;

  TVector3 tangent_down  = down->GetTangent().GetValue();
  TVector3 tangent_up    = up->GetTangentInDownCoordinate().GetValue();
  tangent_down = (1. / tangent_down.Z()) * tangent_down;
  tangent_up = (1. / tangent_up.Z()) * tangent_up;
  TVector3 position_down = down->GetFilmPosition().GetValue();
  TVector3 position_up   = up->GetFilmPositionInDownCoordinate().GetValue();
  
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
						     double &dr, double &dl) {

  TVector3 pos0 = down->GetFilmPosition().GetValue();
  TVector3 pos1 = up->GetFilmPositionInDownCoordinate().GetValue();
  TVector3 dir0 = down->GetTangent().GetValue();
  TVector3 dir1 = up->GetTangent().GetValue();
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

bool ConnectionFunction::JudgeFiducialArea(std::vector<FiducialArea> &area, B2EmulsionSummary *emulsion) {
  
  std::map<double, Point> point_map;
  TVector3 b_pos = emulsion->GetFilmPosition().GetValue();
  TVector3 b_ang = emulsion->GetTangent().GetValue();

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
