#ifndef PID_FUNCTION_HPP
#define PID_FUNCTION_HPP

#include <vector>
#include <map>

#include <TRandom3.h>

#include "McsClass.hpp"
#include "PidClass.hpp"
#include "PidData.hpp"

class PidFunction {

private:
  PidData pid_data_;
  TRandom3 r_;

  // Pid likelihood bin
  std::vector<std::pair<double, double > > pid_ang_bin_vec_;
  std::map<std::pair<double, double >, std::vector<std::pair<double, double > > > pid_mom_bin_map_;  
  // Pid likelihood parameters
  std::vector<Pid_data_ns::DataPoint > likelihood_param_;
  std::map<std::pair<int, int >, Pid_data_ns::DataPoint > likelihood_param_bin_map_;

  std::map<std::pair<double, double >, Pid_data_ns::VphFuncParam > vph_func_param_map_;

  std::map<std::pair<double, double >, double > mean_p_thr_map_;
  std::map<std::pair<double, double >, double > sigma_p_thr_map_;

  // data/fit VPH threshold map
  std::map<std::pair<double, double >, double > vph_mean_cross_map_;
  std::map<std::pair<double, double >, double > vph_sigma_cross_map_;

  void GenerateMeanPThrMap();
  void GenerateLikelihoodParamBinMap();
  void GenerateVphMeanCrossPoint();
  void GenerateVphSigmaCrossPoint();

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

  double CalcMomentumVphEmulsionFit(double pbeta, std::pair<double, double > ang_bin) const;

  double CalcVphMuon(double pbeta, double tangent) const;
  double CalcVphSigmaMuon(double pbeta, double tangent) const;
  double CalcVphProton(double pbeta, double tangent) const;
  double CalcVphSigmaProton(double pbeta, double tangent) const;

  int GetReconPid(double muon_likelihood, double proton_likelihood) const;
  void CalculateStopFlag(Momentum_recon::Mom_chain &chain,
			 std::vector<Momentum_recon::Mom_chain> true_chains) const;

  void CheckMeanSigmaValues() const;
  
};

#endif
