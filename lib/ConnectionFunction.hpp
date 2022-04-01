#ifndef CONNECTION_FUNCTION_HPP
#define CONNECTION_FUNCTION_HPP

#include <boost/filesystem.hpp>

#include <vector>
#include <string>

#include <B2SpillSummary.hh>
#include <B2EmulsionSummary.hh>

#include "McsClass.hpp"
#include "McsFunction.hpp"
#include "ConnectionClass.hpp"
#include "ConnectionData.hpp"

namespace fs = boost::filesystem;



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


public:
  explicit ConnectionFunction(const ConnectionData &connection_data);

  void GetTrueEmulsionTracks(std::vector<const B2EmulsionSummary* > &emulsions,
			     B2SpillSummary &spill_summary, int ecc_id) const;
  void GetTrueEmulsionChains(std::vector<std::vector<const B2EmulsionSummary* > > &chains,
			     std::vector<const B2EmulsionSummary* > &emulsions) const;
  void AddTrueChainsToEventInfo(Momentum_recon::Event_information &ev,
				std::vector<std::vector<const B2EmulsionSummary* > > &chains,
				int ecc_id) const;
  void CalcPosInEccCoordinate(TVector3 &position, int ecc_id) const;
  void SmearEmulsions(std::vector<B2EmulsionSummary* > &emulsions_smeared,
		      std::vector<const B2EmulsionSummary*> &emulsions) const;
  void SmearPosition(TVector3 &position) const;

  void GenerateLinklet(std::vector<std::pair<B2EmulsionSummary*, B2EmulsionSummary* > > &linklet,
		       std::vector<B2EmulsionSummary* > &emulsions) const;

  bool JudgeConnectXY(B2EmulsionSummary* down, B2EmulsionSummary* up, t2l_param &param);
  bool JudgeConnectRL(B2EmulsionSummary* down, B2EmulsionSummary* up, t2l_param &param);
  void CalculatePositionDifference(B2EmulsionSummary* down, B2EmulsionSummary* up, double &dr, double &dl);
  bool JudgeFiducialArea(std::vector<FiducialArea> &area, B2EmulsionSummary *emulsion);
};

#endif
