#ifndef PID_FUNCTION_HPP
#define PID_FUNCTION_HPP

#include <vector>
#include <map>

#include <TRandom3.h>
#include <TH1D.h>

#include "McsClass.hpp"
#include "PidClass.hpp"
#include "PidData.hpp"

class PidFunction {

private:
  PidData pid_data_;
  TRandom3 r_;

  // Pid likelihood parameters
  std::map<int, std::map<double, Pid_data_ns::DataPoint > > likelihood_param_map_;
  std::map<int, double > vph_cross_pbeta_map_;
  std::map<int, double > vph_sigma_thr_map_;

  // VPH function parameters
  std::map<int, Pid_data_ns::VphFuncParam > vph_func_param_map_;

  // Pion MIP parameters
  std::map<int, std::map<double, Pid_data_ns::VphPionMip > > vph_pion_mip_map_;
  std::map<int, double > vph_pion_mip_mean_map_;
  std::map<int, double > vph_pion_mip_thr_map_;
  std::map<int, TH1D* > vph_pion_mip_hist_map_;

  void GenerateVphMeanCrossPBeta();
  void GenerateSigmaThr();

public:
  explicit PidFunction(const PidData &pid_data);

  int GetPidAngBinId(double tangent) const;

  double GetVph(int true_particle_id,
		double ecc_mcs_mom,
		double tangent) const;

  void CalcPartnerLikelihood(double vph,
			     double ecc_mcs_mom,
			     double tangent,
			     double &muon_likelihood,
			     double &proton_likelihood) const;

  int GetReconPid(double vph, double pbeta, double tangent,
		  double muon_likelihood, double proton_likelihood) const;

  double CalcVphPionMipProb(double vph, double tangent) const;

  double CalcVphMuon(double pbeta, double tangent) const;
  double CalcVphSigmaMuon(double pbeta, double tangent) const;
  double CalcVphProton(double pbeta, double tangent) const;
  double CalcVphSigmaProton(double pbeta, double tangent) const;

  std::map<double, Pid_data_ns::DataPoint > GetParamMapByTangent(double tangent) const;
  double GetPBetaCrossPointByTangent(double tangent) const;
  double GetSigmaThrByTangent(double tangent) const;
  Pid_data_ns::VphFuncParam GetFuncParamByTangent(double tangent) const;

  double CalcMomentumVphEmulsionFit(double pbeta, double tangent) const;

  void CalculateStopFlag(Momentum_recon::Mom_chain &chain,
			 std::vector<Momentum_recon::Mom_chain> true_chains) const;
  
};

#endif
