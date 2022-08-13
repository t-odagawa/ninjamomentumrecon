#include <boost/log/core.hpp>
#include <boost/log/trivial.hpp>
#include <boost/log/expressions.hpp>

#include <TMath.h>
#include <TRandom3.h>

#include "PidClass.hpp"
#include "PidData.hpp"
#include "PidFunction.hpp"
#include "AbsorberMaterial_Emulsion_NINJA_Run6.h"

PidFunction::PidFunction(const PidData &pid_data) : pid_data_(pid_data) {

  gRandom->SetSeed(time(NULL));

  pid_data_.GetLikelihoodParam(likelihood_param_map_);
  pid_data_.GetVphFuncParamMapData(vph_func_param_map_);

  pid_data_.GetPionMipParamMap(vph_pion_mip_map_);
  pid_data_.GetPionMipMeanMap(vph_pion_mip_mean_map_);
  pid_data_.GetPionMipThrMap(vph_pion_mip_thr_map_);
  pid_data_.GetPionPdfHistograms(vph_pion_mip_hist_map_);

  GenerateVphMeanCrossPBeta();
  GenerateSigmaThr();

  BOOST_LOG_TRIVIAL(info) << "Pid functions are initialized";

}

void PidFunction::GenerateVphMeanCrossPBeta() {

  for ( auto itr = likelihood_param_map_.begin();
	itr != likelihood_param_map_.end();
	itr++ ) { // angle loop

    double pbeta_cross = -1.;
    double vph_bethe[2] = {};
    double vph_data[2] = {};

    auto datapoints = itr->second;
    auto tmp = datapoints.begin()->second;
    double angle = (tmp.input_ang_min + tmp.input_ang_max) / 2.;

    for ( double pbeta = 1.; pbeta <= 1200.; pbeta += 1. ) {
      vph_bethe[0] = CalcMomentumVphEmulsionFit(pbeta / 1000., angle);
      vph_bethe[1] = CalcMomentumVphEmulsionFit((pbeta + 1.) / 1000., angle);
      Pid_data_ns::DataPoint data[2];
      data[0] = datapoints.lower_bound(pbeta)->second;
      if ( datapoints.lower_bound(pbeta) == datapoints.begin() ) {
	data[1] = std::next(datapoints.lower_bound(pbeta), 1)->second;
      }
      else {
	data[1] = std::next(datapoints.lower_bound(pbeta), -1)->second;
      }

      double pbeta_mean[2];
      pbeta_mean[0] = (data[0].input_mom_min + data[0].input_mom_max) / 2.;
      pbeta_mean[1] = (data[1].input_mom_min + data[1].input_mom_max) / 2.;

      vph_data[0] = (data[0].mean[2] - data[1].mean[2]) / (pbeta_mean[0] - pbeta_mean[1]) * (pbeta - pbeta_mean[0])
	+ data[0].mean[2];
      vph_data[1] = (data[0].mean[2] - data[1].mean[2]) / (pbeta_mean[0] - pbeta_mean[1]) * (pbeta + 1. - pbeta_mean[0])
	+ data[0].mean[2];

      if ( vph_bethe[0] >= vph_data[0] &&
	   vph_bethe[1] < vph_data[1] ) {
	pbeta_cross = pbeta + std::min(1., 1. * (vph_bethe[0] - vph_data[0])
				       / ((vph_data[1] - vph_bethe[1]) + (vph_bethe[0] - vph_data[0])));
	vph_cross_pbeta_map_.insert(std::make_pair(itr->first, pbeta_cross));

	BOOST_LOG_TRIVIAL(trace) << "Angle : (" << data[0].input_ang_min << ", " << data[0].input_ang_max
				 << ") -> VPH Threshold : " << (vph_data[0] + vph_data[1]) / 2.
				 << " at pbeta = " << pbeta_cross;

	break;
      }
    } // pbeta loop

  }

  return;

}

void PidFunction::GenerateSigmaThr() {

  for ( auto itr = likelihood_param_map_.begin();
	itr != likelihood_param_map_.end();
	itr++ ) { // angle loop
    auto datapoints = itr->second;

    for ( auto itr1 = datapoints.begin(); itr1 != datapoints.end(); itr1++ ) {
      auto data = itr1->second;
      if ( 1200. < (data.input_mom_min + data.input_mom_max) / 2. &&
	   (data.input_mom_min + data.input_mom_max) / 2. < 1300. ) {
	vph_sigma_thr_map_.insert(std::make_pair(itr->first, data.sigma[1]));

	BOOST_LOG_TRIVIAL(trace) << "Angle : (" << data.input_ang_min << ", " << data.input_ang_max
				 << ") -> Sigma Threshold : " << data.sigma[1]
				 << " at pbeta = " << (data.input_mom_min + data.input_mom_max) / 2.;

	break;
      }
    }
  }

  return;

}

// 与えられた particle id, pbeta, tangent をもとに
// data driven な mean/sigma を使って gauss で VPH を計算
double PidFunction::GetVph(int true_particle_id,
			   double pbeta,
			   double tangent) const {

  double vph = -1;
  
  if ( true_particle_id == 2212 ) {
    double mean = CalcVphProton(pbeta, tangent);
    double sigma = CalcVphSigmaProton(pbeta, tangent);
    vph = gRandom->Gaus(mean, sigma);
  }
  else {
    double pbeta_thr;
    if ( vph_pion_mip_thr_map_.upper_bound((int)(tangent * 10)) == vph_pion_mip_thr_map_.end() ) {
      pbeta_thr = vph_pion_mip_thr_map_.rbegin()->second;
    }
    else {
      pbeta_thr = vph_pion_mip_thr_map_.upper_bound((int)(tangent * 10))->second;
    }
    if ( pbeta > pbeta_thr + 100. ) { // MIP
      TH1D *hist;
      if ( vph_pion_mip_hist_map_.upper_bound((int)(tangent * 10)) == vph_pion_mip_hist_map_.end() ) {
	hist = vph_pion_mip_hist_map_.rbegin()->second;
      }
      else {
	hist = vph_pion_mip_hist_map_.upper_bound((int)(tangent * 10))->second;
      }
      vph = hist->GetRandom();
    }
    else {
      double mean = CalcVphMuon(pbeta, tangent);
      double sigma = CalcVphSigmaMuon(pbeta, tangent);
      vph = gRandom->Gaus(mean, sigma);
    }
  }

  if ( vph < 0 ) vph = 0.;
  return vph;

}

void PidFunction::CalcPartnerLikelihood(double vph, double pbeta, double tangent,
					double &muon_likelihood,
					double &proton_likelihood) const {

  double proton_mean = CalcVphProton(pbeta, tangent);
  double proton_sigma = CalcVphSigmaProton(pbeta, tangent);
  proton_likelihood = (1 + std::erf((vph - proton_mean) / (std::sqrt(2) * proton_sigma))) / 2.;

  double pbeta_thr;
  if ( vph_pion_mip_thr_map_.upper_bound((int)(tangent * 10)) == vph_pion_mip_thr_map_.end() ) {
    pbeta_thr = vph_pion_mip_thr_map_.rbegin()->second;
  }
  else {
    pbeta_thr = vph_pion_mip_thr_map_.upper_bound((int)(tangent * 10))->second;
  }

  if ( pbeta < pbeta_thr + 100. ) {
    double muon_mean = CalcVphMuon(pbeta, tangent);
    double muon_sigma = CalcVphSigmaMuon(pbeta, tangent);
    muon_likelihood = 1 - (1 + std::erf((vph - muon_mean) / (std::sqrt(2) * muon_sigma))) / 2.;
  }
  else {
    muon_likelihood = CalcVphPionMipProb(vph, tangent);
  }
  
  return;

}

int PidFunction::GetReconPid(double vph, double pbeta, double tangent,
			     double muon_likelihood, double proton_likelihood) const {

  double pbeta_thr;
  if ( vph_pion_mip_thr_map_.upper_bound((int)(tangent * 10)) == vph_pion_mip_thr_map_.end() ) {
    pbeta_thr = vph_pion_mip_thr_map_.rbegin()->second;
  }
  else {
    pbeta_thr = vph_pion_mip_thr_map_.upper_bound((int)(tangent * 10))->second;
  }

  double likelihood_ratio = muon_likelihood / (muon_likelihood + proton_likelihood);
  double proton_mean = CalcVphProton(pbeta, tangent);
  double pion_mean;

  if ( pbeta < pbeta_thr + 100. ) {
    if ( vph_pion_mip_mean_map_.upper_bound((int)(tangent * 10)) == vph_pion_mip_mean_map_.end() ) {
      pion_mean = vph_pion_mip_mean_map_.rbegin()->second;
    }
    else {
      pion_mean = vph_pion_mip_mean_map_.upper_bound((int)(tangent * 10))->second;
    }
  }
  else {
    pion_mean = CalcVphMuon(pbeta, tangent);
  }

  if ( vph > proton_mean ) {
    return 2212;
  }
  else if ( vph < pion_mean ) {
    return 211;
  }
  else if ( likelihood_ratio <= 0.5 ) {
    return 2212;
  }
  else 
    return 211;

}
			     
double PidFunction::CalcVphPionMipProb(double vph, double tangent) const {

  double sum = 0.;
  Pid_data_ns::VphPionMip point[2];

  std::map<double, Pid_data_ns::VphPionMip > param;
  if ( vph_pion_mip_map_.upper_bound((int)(tangent * 10)) == vph_pion_mip_map_.end() ) {
    param = vph_pion_mip_map_.rbegin()->second;
  }
  else {
    param = vph_pion_mip_map_.upper_bound((int)(tangent * 10))->second;
  }
  

  for ( auto itr = param.rbegin(); itr != std::next(param.rend(),-1); itr++ ) {
    point[0] = std::next(itr, 1)->second;
    point[1] = itr->second;
    if ( point[1].vph > vph ) {
      sum += (point[0].prob + point[1].prob) / 2.;
    }
    else if ( point[0].vph > vph ) {
      double prob_middle;
      prob_middle =  point[0].prob + (point[1].prob - point[0].prob) / (point[1].vph - point[0].vph)
	* (vph - point[0].vph);
      double prob_middle_area = (prob_middle + point[1].prob) * (point[1].vph - vph) / 2.;
      double prob_middle_scale = prob_middle_area / ((point[1].prob + point[0].prob) * (point[1].vph - point[0].vph) / 2.)
	* ((point[0].prob + point[1].prob) / 2.);

      sum += prob_middle_scale;
      break;
    }
  }

  return sum;

}

double PidFunction::CalcVphMuon(double pbeta, double tangent) const {

  double mean = 0.;
  
  // 角度ビンを決めて datapoint map を選ぶ
  auto datapoints = GetParamMapByTangent(tangent);

  Pid_data_ns::DataPoint data[2];
  if ( datapoints.rbegin()->second.input_mom_min < pbeta ) {
    data[0] = datapoints.rbegin()->second;
    data[1] = std::next(datapoints.rbegin(), 1)->second;
  }
  else {
    data[0] = datapoints.lower_bound(pbeta)->second;
    if ( datapoints.lower_bound(pbeta) == datapoints.begin() ) {
      data[1] = std::next(datapoints.lower_bound(pbeta), 1)->second;
    }
    else {
      data[1] = std::next(datapoints.lower_bound(pbeta), -1)->second;
    }
  }

  double pbeta_mean[2];
  pbeta_mean[0] = (data[0].input_mom_min + data[0].input_mom_max) / 2.;
  pbeta_mean[1] = (data[1].input_mom_min + data[1].input_mom_max) / 2.;
  mean = (data[0].mean[1] - data[1].mean[1]) / (pbeta_mean[0] - pbeta_mean[1]) * (pbeta - pbeta_mean[0])
    + data[0].mean[1];
  return mean;

}

double PidFunction::CalcVphSigmaMuon(double pbeta, double tangent) const {

  double sigma = 0.;

  // 角度ビンを決めて datapoint map を選ぶ
  auto datapoints = GetParamMapByTangent(tangent);

  Pid_data_ns::DataPoint data[2];
  if ( datapoints.rbegin()->second.input_mom_min < pbeta ) {
    data[0] = datapoints.rbegin()->second;
    data[1] = std::next(datapoints.rbegin(), 1)->second;
  }
  else { // pbeta が最小 bin より小さい場合
    data[0] = datapoints.lower_bound(pbeta)->second;
    if ( datapoints.lower_bound(pbeta) == datapoints.begin() ) {
      data[1] = std::next(datapoints.lower_bound(pbeta), 1)->second;
    }
    else {
      data[1] = std::next(datapoints.lower_bound(pbeta), -1)->second;
    }    
  }

  double pbeta_mean[2];
  pbeta_mean[0] = (data[0].input_mom_min + data[0].input_mom_max) / 2.;
  pbeta_mean[1] = (data[1].input_mom_min + data[1].input_mom_max) / 2.;
  sigma = (data[0].sigma[1] - data[1].sigma[1]) / (pbeta_mean[0] - pbeta_mean[1]) * (pbeta - pbeta_mean[0])
    + data[0].sigma[1];
  return sigma;

}

double PidFunction::CalcVphProton(double pbeta, double tangent) const {

  double mean = 0.;

  double pbeta_cross_point = GetPBetaCrossPointByTangent(tangent);

  if ( pbeta < pbeta_cross_point ) {

    auto datapoints = GetParamMapByTangent(tangent);

    Pid_data_ns::DataPoint data[2];
    data[0] = datapoints.lower_bound(pbeta)->second;
    if ( datapoints.lower_bound(pbeta) == datapoints.begin() ) {
      data[1] = std::next(datapoints.lower_bound(pbeta), 1)->second;
    }
    else {
      data[1] = std::next(datapoints.lower_bound(pbeta), -1)->second;
    }

    double pbeta_mean[2];
    pbeta_mean[0] = (data[0].input_mom_min + data[0].input_mom_max) / 2.;
    pbeta_mean[1] = (data[1].input_mom_min + data[1].input_mom_max) / 2.;
    mean = (data[0].mean[2] - data[1].mean[2]) / (pbeta_mean[0] - pbeta_mean[1]) * (pbeta - pbeta_mean[0])
      + data[0].mean[2];
  }
  else {
    mean = CalcMomentumVphEmulsionFit(pbeta / 1000., tangent);
  }

  return mean;

}

double PidFunction::CalcVphSigmaProton(double pbeta, double tangent) const {

  double sigma = 0.;

  double pbeta_cross_point = GetPBetaCrossPointByTangent(tangent);

  double vph = CalcMomentumVphEmulsionFit(pbeta / 1000., tangent);
  double vph_thr = CalcMomentumVphEmulsionFit(pbeta / 1000., tangent);

  if ( pbeta < 200. ) {
    auto datapoints = GetParamMapByTangent(tangent);

    Pid_data_ns::DataPoint data[2];
    if ( datapoints.rbegin()->second.input_mom_min < pbeta ) {
      data[0] = datapoints.rbegin()->second;
      data[1] = std::next(datapoints.rbegin(), 1)->second;
    }
    else {
      data[0] = datapoints.lower_bound(pbeta)->second;
      if ( datapoints.lower_bound(pbeta) == datapoints.begin() ) {
	data[1] = std::next(datapoints.lower_bound(pbeta), 1)->second;
      }
      else {
	data[1] = std::next(datapoints.lower_bound(pbeta), -1)->second;
      }
    }

    double pbeta_mean[2];
    pbeta_mean[0] = (data[0].input_mom_min + data[0].input_mom_max) / 2.;
    pbeta_mean[1] = (data[1].input_mom_min + data[1].input_mom_max) / 2.;
    sigma = (data[0].sigma[2] - data[1].sigma[2]) / (pbeta_mean[0] - pbeta_mean[1]) * (pbeta - pbeta_mean[0])
      + data[0].sigma[2];
    return sigma;
  }

  auto func_param = GetFuncParamByTangent(tangent);
  sigma = func_param.sigma_inter + func_param.sigma_scale * std::sqrt(vph);
  double sigma_thr = GetSigmaThrByTangent(tangent);
  return std::max(sigma, sigma_thr);

}

std::map<double, Pid_data_ns::DataPoint > PidFunction::GetParamMapByTangent(double tangent) const {

  int iang = tangent * 10;

  if ( likelihood_param_map_.upper_bound(iang) == likelihood_param_map_.end() ) {
    return likelihood_param_map_.rbegin()->second;
  }
  else {
    return likelihood_param_map_.upper_bound(iang)->second;
  }

}

double PidFunction::GetPBetaCrossPointByTangent(double tangent) const {

  int iang = tangent * 10;
  
  if ( vph_cross_pbeta_map_.upper_bound(iang) == vph_cross_pbeta_map_.end() ) {
    return vph_cross_pbeta_map_.rbegin()->second;
  }
  else {
    return vph_cross_pbeta_map_.upper_bound(iang)->second;
  }

}

double PidFunction::GetSigmaThrByTangent(double tangent) const {

  int iang = tangent * 10;
  
  if ( vph_sigma_thr_map_.upper_bound(iang) == vph_sigma_thr_map_.end() ) {
    return vph_sigma_thr_map_.rbegin()->second;
  }
  else {
    return vph_sigma_thr_map_.upper_bound(iang)->second;
  }

}

Pid_data_ns::VphFuncParam PidFunction::GetFuncParamByTangent(double tangent) const {

  int iang = tangent * 10;

  if ( vph_func_param_map_.upper_bound(iang) == vph_func_param_map_.end() ) {
    return vph_func_param_map_.rbegin()->second;
  }
  else {
    return vph_func_param_map_.upper_bound(iang)->second;
  }

}

double PidFunction::CalcMomentumVphEmulsionFit(double pbeta, double tangent) const {

  auto func_param_ = GetFuncParamByTangent(tangent);

  double dedx = 0.;
  
  double mass = 938.e-3; // proton mass
  double z = 1.; // charge
  
  // information constants
  double electron_mass = 0.511e-3; // electron mass
  double K = 0.307075; // 4 pi N A re^2 me c^2 [MeV mol^-1 cm^2] 

  for ( int i = 0; i < 12; i++ ) {

    switch (i) {
    case 0  : Absorber = Ag; break;
    case 1  : Absorber = Br; break;
    case 2  : Absorber = I;  break;
    case 3  : Absorber = C;  break;
    case 4  : Absorber = N;  break;
    case 5  : Absorber = O;  break;
    case 6  : Absorber = H;  break;
    case 7  : Absorber = S;  break;
    case 8  : Absorber = Na; break;
    case 9  : Absorber = Fe; break;
    case 10 : Absorber = Au; break;
    case 11 : Absorber = Cl; break;
    default : break;
    }

    // information of absorber
    double Z = mtrl_Z[Absorber];
    double A = mtrl_A[Absorber];
    double I = mtrl_I[Absorber] * 1e-9; // GeV
    double density = mtrl_Density[Absorber]; // g cm-3

    // Calculation of dE/dx
    double E = ( pbeta + TMath::Sqrt(pbeta * pbeta + 4 * mass * mass)) / 2.; // Energy [GeV]
    double p = TMath::Sqrt(pbeta * E); // Momentum [GeV/c]
    double beta = p / E;
    double gamma = 1. / TMath::Sqrt(1. - beta * beta);

    double W = (2. * electron_mass * beta * beta * gamma * gamma)
      / (1. + 2 * gamma * electron_mass / mass + (electron_mass / mass) * (electron_mass / mass));

    // Density effect
    double  C_bar = mtrl_C_bar[Absorber];
    double      a = mtrl_a[Absorber];
    double      m = mtrl_m[Absorber];
    double     x1 = mtrl_x1[Absorber];
    double     x0 = mtrl_x0[Absorber];
    double delta0 = mtrl_delta0[Absorber];

    double x = TMath::Log10(p / mass);
    double density_effect = 0.;
    
    if ( x >= x1 )
      density_effect = (2. * TMath::Log(10) * x  - C_bar) / 2.;
    else if ( x >= x0 )
      density_effect = (2. * TMath::Log(10) * x  - C_bar + a * TMath::Power(x1-x,m)) / 2.;
    else
      density_effect = (delta0 * TMath::Power(10., 2.*(x-x0))) / 2.;

    double wratio = mtrl_wratio[Absorber];

    dedx += wratio * (K * z * z * (Z / A) * (1. / beta / beta)
		      * (0.5 * TMath::Log(2. * electron_mass * beta * beta * gamma * gamma * W / (I * I)) - beta * beta - density_effect));
    
  }

  double slope = func_param_.mean_slope;
  double intercept = func_param_.mean_inter;
  
  return slope * dedx + intercept;

}

void PidFunction::CalculateStopFlag(Momentum_recon::Mom_chain &chain,
				    std::vector<Momentum_recon::Mom_chain> true_chains) const {

  int recon_particle_flag = chain.particle_flag % 10000;

  if ( recon_particle_flag != 2212 &&
       std::abs(recon_particle_flag) != 211 ) return;
  if ( chain.stop_flag != -1 ) return;

  for ( auto true_chain : true_chains ) {
    if ( chain.chainid == true_chain.chainid ) {
      
      // 最下端が true と一致しているか，かつ
      // true が ECC 中で stop しているか
      bool edge_identity_flag = false;
      // true の direction を確認
      int direction = true_chain.direction;
      if ( direction == 1 && 
	   chain.base.front().rawid == true_chain.base.front().rawid ) {
	  edge_identity_flag = true;
      }
      else if ( direction == -1 &&
		chain.base.back().rawid == true_chain.base.back().rawid ) {
	edge_identity_flag = true;
      }

      if ( true_chain.stop_flag == 2 &&
	   edge_identity_flag ) {
	chain.stop_flag = 2;
      }
      else {
	chain.stop_flag = 0;
      }
      break;
    }
  }

  return;

}
