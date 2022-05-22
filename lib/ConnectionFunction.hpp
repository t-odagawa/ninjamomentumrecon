#ifndef CONNECTION_FUNCTION_HPP
#define CONNECTION_FUNCTION_HPP

#include <boost/filesystem.hpp>
#include <boost/unordered_set.hpp>
#include <boost/unordered_map.hpp>

#include <vector>
#include <array>
#include <string>
#include <set>
#include <unordered_set>
#include <map>
#include <unordered_map>

#include <B2SpillSummary.hh>
#include <B2EmulsionSummary.hh>

#include <TVector3.h>

#include "McsClass.hpp"
#include "McsFunction.hpp"
#include "ConnectionClass.hpp"
#include "ConnectionData.hpp"

namespace fs = boost::filesystem;

bool CompareParticlePdg(std::vector<const B2EmulsionSummary* > lhs,
			std::vector<const B2EmulsionSummary* > rhs);

class ConnectionFunction {

private:
  ConnectionData connection_data_;

  // Allowances
  t2l_param et_ecc_param_;
  t2l_param et_ecc_fe_param_;
  t2l_param et_ecc_fe_fe_param_;
  t2l_param black_fe_param_;
  t2l_param black_re_fe_param_;
  t2l_param black_re_fe_fe_param_;
  t2l_param black_re_fe_water_param_;
  t2l_param black_water_param_;
  t2l_param fe_param_;
  t2l_param fe_fe_param_;
  t2l_param fe_fe_fe_param_;
  t2l_param fe_water_param_;
  t2l_param fe_water_fe_param_;
  t2l_param re_et_ecc_fe_fe_fe_param_;
  t2l_param re_et_ecc_fe_fe_fe_fe_param_;
  t2l_param re_et_ecc_fe_fe_fe_water_param_;
  t2l_param re_et_ecc_fe_fe_water_param_;
  t2l_param re_et_ecc_fe_water_param_;
  t2l_param re_fe_param_;
  t2l_param re_fe_fe_param_;
  t2l_param re_fe_fe_fe_param_;
  t2l_param re_fe_fe_fe_fe_param_;
  t2l_param re_fe_fe_fe_fe_fe_param_;
  t2l_param re_fe_water_param_;
  t2l_param re_fe_water_fe_param_;
  t2l_param re_fe_water_fe_water_param_;
  t2l_param re_fe_water_fe_water_fe_param_;
  t2l_param re_water_fe_water_param_;
  t2l_param water_param_;

  std::map<int, std::vector<FiducialArea > > ecc_fiducial_[9];
  std::map<int, std::vector<Efficiency > > ecc_efficiency_[9];

public:
  explicit ConnectionFunction(const ConnectionData &connection_data);

  void GetTrueEmulsionTracks(std::vector<const B2EmulsionSummary* > &emulsions,
			     B2SpillSummary &spill_summary, int ecc_id) const;
  void GetTrueEmulsionChains(std::vector<std::vector<const B2EmulsionSummary* > > &chains,
			     std::vector<const B2EmulsionSummary* > &emulsions) const;
  void AddTrueChainsToEventInfo(Momentum_recon::Event_information &ev,
				std::vector<std::vector<const B2EmulsionSummary* > > &chains,
				int ecc_id, int material_id) const;
  void CalcPosInEccCoordinate(TVector3 &position, int ecc_id) const;
  void SmearEmulsions(std::vector<B2EmulsionSummary* > &emulsions_smeared,
		      std::vector<const B2EmulsionSummary*> &emulsions) const;
  void SmearPosition(TVector3 &position) const;

  void ApplyDetectionEfficiency(std::vector<B2EmulsionSummary* > &emulsions_detected,
				std::vector<B2EmulsionSummary* > &emulsions_smeared,
				int ecc_id) const;

  int GetAngleBinId(double tangent) const;

  void ApplyFVCut(std::vector<B2EmulsionSummary* > &emulsions_detected_in_fv,
		  std::vector<B2EmulsionSummary* > &emulsions_detected,
		  int ecc_id) const;

  void GenerateLinklet(std::vector<std::pair<B2EmulsionSummary*, B2EmulsionSummary* > > &linklet,
		       std::vector<B2EmulsionSummary* > &emulsions) const;

  void GenerateGroup(std::vector<std::pair<B2EmulsionSummary*, std::vector<std::pair<B2EmulsionSummary*, B2EmulsionSummary* > > > > &groups,
		     std::vector<std::pair<B2EmulsionSummary*, B2EmulsionSummary* > > &linklets,
		     std::vector<B2EmulsionSummary* > &emulsions) const;

  void ConvertRawIdToB2(unsigned int rawid, B2EmulsionSummary* &emulsion,
			std::vector<B2EmulsionSummary* > &emulsions ) const;

  void LinkletConvertB2ToSeg(std::pair<B2EmulsionSummary*, B2EmulsionSummary* > &linklet_b2,
			     std::pair<Segment, Segment > &linklet_seg) const;

  void LinkletConvertSegToB2(std::pair<Segment, Segment > &linklet_seg,
			     std::pair<B2EmulsionSummary*, B2EmulsionSummary* > &linklet_b2,
			     std::vector<B2EmulsionSummary* > &emulsions) const;

  void GroupConvertB2ToSeg(std::pair<B2EmulsionSummary*, std::vector<std::pair<B2EmulsionSummary*, B2EmulsionSummary* > > > &group_b2, Group &group_seg) const;

  void GroupConvertSegToB2(Group &group_seg,
			   std::pair<B2EmulsionSummary*, std::vector<std::pair<B2EmulsionSummary*, B2EmulsionSummary* > > > &group_b2,
			   std::vector<B2EmulsionSummary* > &emulsions) const;
 
  void GenerateGroupNext(B2EmulsionSummary* emulsion,
			 boost::unordered_multimap<Segment, Segment > &link_next,
			 boost::unordered_multimap<Segment, Segment > &link_prev,
			 std::vector<std::pair<Segment, Segment > > &group_next) const;
  void GenerateGroupPrev(B2EmulsionSummary* emulsion,
			 boost::unordered_multimap<Segment, Segment > &link_next,
			 boost::unordered_multimap<Segment, Segment > &link_prev,
			 std::vector<std::pair<Segment, Segment > > &group_prev) const;
  void GroupMerge(B2EmulsionSummary* emulsion,
		  std::vector<std::pair<Segment, Segment > > &group_next,
		  std::vector<std::pair<Segment, Segment > > &group_prev,
		  Group &group) const;

  void GenerateGroupNextPartner(B2EmulsionSummary* emulsion,
				boost::unordered_multimap<Segment, Segment > &link_next,
				boost::unordered_multimap<Segment, Segment > &link_prev,
				std::vector<std::pair<Segment, Segment > > &group_next_partner) const;
  void GenerateGroupPrevPartner(B2EmulsionSummary* emulsion,
				boost::unordered_multimap<Segment, Segment > &link_next,
				boost::unordered_multimap<Segment, Segment > &link_prev,
				std::vector<std::pair<Segment, Segment > > &group_prev_partner) const;
  
  void ReconnectGroups(std::vector<std::pair<B2EmulsionSummary*, std::vector<std::pair<B2EmulsionSummary*, B2EmulsionSummary* > > > > &groups_reconnected,
		       std::vector<std::pair<B2EmulsionSummary*, std::vector<std::pair<B2EmulsionSummary*, B2EmulsionSummary* > > > > &groups,
		       std::vector<B2EmulsionSummary* > &emulsions) const;

  void ModifyGroups(std::vector<std::pair<B2EmulsionSummary*, std::vector<std::pair<B2EmulsionSummary*, B2EmulsionSummary* > > > > &groups_modified,
					std::vector<std::pair<B2EmulsionSummary*, std::vector<std::pair<B2EmulsionSummary*, B2EmulsionSummary* > > > > &groups_reconnected,
					std::vector<B2EmulsionSummary* > &emulsions) const;

  void ModifyGroup(std::pair<B2EmulsionSummary*, std::vector<std::pair<B2EmulsionSummary*, B2EmulsionSummary* > > > &group_modified,
				       std::pair<B2EmulsionSummary*, std::vector<std::pair<B2EmulsionSummary*, B2EmulsionSummary* > > > &group,
				       std::vector<B2EmulsionSummary* > &emulsions) const;

  bool IsValidGroup(std::pair<B2EmulsionSummary*, std::vector<std::pair<B2EmulsionSummary*, B2EmulsionSummary* > > > &group) const;

  bool SelectMuonGroup(std::vector<std::pair<B2EmulsionSummary*, std::vector<std::pair<B2EmulsionSummary*, B2EmulsionSummary* > > > > &groups,
		       std::pair<B2EmulsionSummary*, std::vector<std::pair<B2EmulsionSummary*, B2EmulsionSummary* > > > &muon_group,
		       B2EmulsionSummary* &vertex_track,
		       std::vector<B2EmulsionSummary* > &emulsions) const;

  void CollectEdgeTracks(boost::unordered_multimap<Segment, Segment > &upstream_tracks,
			 boost::unordered_multimap<Segment, Segment > &downstream_tracks,
			 boost::unordered_multimap<Segment, Segment > &upstream_tracks_start,
			 boost::unordered_multimap<Segment, Segment > &downstream_tracks_start,
			 Group &group, Segment &start_seg) const;

  void PartnerSearch(B2EmulsionSummary* vertex_track,
		     std::vector<std::pair<B2EmulsionSummary*, std::vector<std::pair<B2EmulsionSummary*, B2EmulsionSummary* > > > > &groups,
		     std::vector<B2EmulsionSummary*> &emulsions_partner,
		     std::vector<B2EmulsionSummary*> &emulsions,
		     int ecc_id,
		     TVector3 &recon_vertex) const;
  void SelectPartnerGroups(std::vector<std::pair<B2EmulsionSummary*, std::vector<std::pair<B2EmulsionSummary*, B2EmulsionSummary* > > > > &groups,
			   std::vector<std::pair<B2EmulsionSummary*, std::vector<std::pair<B2EmulsionSummary*, B2EmulsionSummary* > > > > &groups_partner,
			   std::vector<B2EmulsionSummary* > &emulsions_partner,
			   std::vector<B2EmulsionSummary* > &emulsions) const;

  void AddGroupsToEventInfo(Momentum_recon::Event_information &ev,
			    B2EmulsionSummary* vertex_track,
			    std::pair<B2EmulsionSummary*, std::vector<std::pair<B2EmulsionSummary*, B2EmulsionSummary* > > > &muon_group,
			    std::vector<std::pair<B2EmulsionSummary*, std::vector<std::pair<B2EmulsionSummary*, B2EmulsionSummary* > > > > &groups_partner,
			    std::vector<B2EmulsionSummary* > &emulsions,
			    int ecc_id) const;

  void AddGroupToEventInfo(Momentum_recon::Event_information &ev,
			   B2EmulsionSummary* vertex_track,
			   std::pair<B2EmulsionSummary*, std::vector<std::pair<B2EmulsionSummary*, B2EmulsionSummary* > > > &muon_group,
			   std::pair<B2EmulsionSummary*, std::vector<std::pair<B2EmulsionSummary*, B2EmulsionSummary* > > > &group,
			   std::vector<B2EmulsionSummary* > &emulsions,
			   int ecc_id) const;

  bool JudgeConnect(B2EmulsionSummary* down, B2EmulsionSummary* up, t2l_param param) const;
  bool JudgeConnectXY(B2EmulsionSummary* down, B2EmulsionSummary* up, t2l_param param) const;
  bool JudgeConnectRL(B2EmulsionSummary* down, B2EmulsionSummary* up, t2l_param param) const;
  void CalculatePositionDifference(B2EmulsionSummary* down, B2EmulsionSummary* up, double &dr, double &dl) const;
  bool JudgeFiducialArea(const std::vector<FiducialArea> &area, B2EmulsionSummary *emulsion) const;
  bool JudgeEdgeOut(B2EmulsionSummary* vertex_track, int ecc_id) const;
  bool JudgePartnerTrack(B2EmulsionSummary* vertex_track, B2EmulsionSummary* partner_track,\
			 std::vector<std::pair<B2EmulsionSummary*, std::vector<std::pair<B2EmulsionSummary*, B2EmulsionSummary* > > > > &groups,
			 int ecc_id, TVector3 &recon_vertex) const;
  double GetMinimumDistance(TVector3 parent_pos, TVector3 daughter_pos,
			    TVector3 parent_dir, TVector3 daughter_dir,
			    std::array<Double_t,2> z_range,
			    std::array<Double_t,2> &extrapolate_z,
			    int ecc_id,
			    TVector3 &recon_vertex) const;
  int GetVertexMaterial(TVector3 recon_vertex) const;
  int GetVertexTargetMaterial(int vertex_plate) const;
};

#endif
