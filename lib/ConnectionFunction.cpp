std::vector<const B2EmulsionSummary* > ConnectionFunction::GetTrueEmulsionTracks(B2SpillSummary &spill_summary, Int_t ecc_id) {
  std::vector<const B2EmulsionSummary* > emulsions;
  auto it_emulsion = spill_summary.BeginEmulsion();
  while ( const auto emulsion = it_emulsion.Next() ) {
    if ( emulsion->GetParentTrackId() == 0 ) continue;
    if ( emulsion->GetFilmType() != B2FilmType::kECC ) continue;
    if ( emulsion->GetEcc() != ecc_id ) continue;
    emulsions.push_back(emulsion);
  }
  return emulsions;
}

std::vector<std::vector<const B2EmulsionSummary* > > ConnectionFunction::GetTrueEmulsionChains(std::vector<const B2EmulsionSummary* > emulsions) {
  std::vector<std::vector<const B2EmulsionSummary* > > chains;
  std::vector<const B2EmulsionSummary* > chain;
  Int_t track_id_tmp = emulsions.front()->GetParentTrackId();
  for ( auto &emulsion : emulsions ) {
    if ( emulsion->GetParentTrackId() == track_id_tmp )
      chain.push_back(chain);
    else {
      chains.push_back(chain);
      chain.clear();
      track_id_tmp = emulsion.GetParentTrackId();
      chain.push_back(emulsion);
    }
  }
  chains.push_back(chain);
  return chains;
}

void ConnectionFunction::ReadFiducialVolumeData(std::string path) {

}

void ConnectionFunction::ReadAllowanceData(std::string path) {
  
}

std::vector<B2EmusionSummary* > ConnectionFunction::ModifyChain(std::vector<B2EmulsionSummary* > chain) {
  
}

// down : alignementの基準, up : alignementで補正されうる
bool ConnectionFunction::JudgeConnectXY(B2EmulsionSummary* down, B2EmulsionSummary* up, t2l_param &param) {

  double d_pos_x, d_pos_y, d_ang_x, d_ang_y;
  double all_pos_x, all_pos_y, all_ang_x, all_ang_y;

  TVector3 tangent_down  = down.GetTangent().GetValue();
  TVector3 tangent_up    = up.GetTangent().GetValue();
  TVector3 position_down = down.GetFilmPosition().GetValue();
  TVector3 position_up   = up.GetFilmPositionInDownCoordinate().GetValue();

  all_ang_x = param.intercept_ax
    + param.slope_ax * std::fabs(tangent_down.X())
    + param.slope2_ax * std::pow(tangent_down.X(). 2.);
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

  TVector3 tangent_down  = down.GetTangent().GetValue();
  TVector3 tangent_up    = up.GetTangentInDownCoordinate().GetValue();
  tangent_down = (1. / tangent_down.Z()) * tangent_down;
  tangent_up = (1. / tangent_up.Z()) * tangent_up;
  TVector3 position_down = down.GetFilmPosition().GetValue();
  TVector3 position_up   = up.GetFilmPositionInDownCoordinate().GetValue();
  
  angle = std::sqrt(tangent_down.X() * tangent_down.X() + tangent_down.Y() * tangetn_down.Y());
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

  TVector3 pos0 = down.GetFilmPosition().GetValue();
  TVector3 pos1 = up.GetFilmPositionInDownCoordinate().GetValue();
  TVector3 dir0 = down.GetTangent().GetValue();
  TVector3 dir1 = up.GetTangent().GetValue();
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
  
  std::map<double, point> point_map;
  TVector3 b_pos = emulsion.GetFilmPosition().GetValue();
  TVector3 b_ang = emulsion.GetTangent().GetValue();

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
