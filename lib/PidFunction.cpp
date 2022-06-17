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

  pid_data_.GetAngBinVectorData(pid_ang_bin_vec_);
  pid_data_.GetMomBinMapData(pid_mom_bin_map_);
  pid_data_.GetLikelihoodParam(likelihood_param_);
  pid_data_.GetVphFuncParamMapData(vph_func_param_map_);

  GenerateLikelihoodParamBinMap();
  GenerateVphMeanCrossPoint();
  GenerateVphSigmaCrossPoint();

  BOOST_LOG_TRIVIAL(info) << "Pid functions are initialized";

}

void PidFunction::GenerateLikelihoodParamBinMap() {

  for ( int iang = 0; iang < pid_ang_bin_vec_.size(); iang++ ) {
    auto ang_bin = pid_ang_bin_vec_.at(iang);
    auto pid_mom_bin_vec_ = pid_mom_bin_map_.at(ang_bin);
    for ( int imom = 0; imom < pid_mom_bin_vec_.size(); imom++ ) {
      auto mom_bin = pid_mom_bin_vec_.at(imom);
      auto bin = std::make_pair(iang, imom);
      for ( auto param : likelihood_param_ ) {
	if ( param.input_ang_min == ang_bin.first &&
	     param.input_ang_max == ang_bin.second &&
	     param.input_mom_min == mom_bin.first &&
	     param.input_mom_max == mom_bin.second ) {
	  likelihood_param_bin_map_.insert(std::make_pair(bin, param));
	}
      }
    }
  }

  return;
  
}

// data と fit の境界を決める pbeta に対応する VPH mean の値を求める
void PidFunction::GenerateVphMeanCrossPoint() {

  // 角度ごとに mean cross point を作る
  double mom_data_fit_thr = 200.;

  for ( auto ang_bin : pid_ang_bin_vec_ ) {
    // data と fit の境界を決める pbeta に対応する VPH
    double vph = CalcMomentumVphEmulsionFit(mom_data_fit_thr / 1000., ang_bin);
    vph_mean_cross_map_.insert(std::make_pair(ang_bin, vph));

    BOOST_LOG_TRIVIAL(trace) << "Angle : (" << ang_bin.first << ", " << ang_bin.second << ")"
			     << " -> Threshold VPH : " << vph;
			     
  }

  return;

}

// pion と proton の sigma が一緒になる境界を決める
void PidFunction::GenerateVphSigmaCrossPoint() {


  // 角度ごとに sigma cross point を作る
  double mom_pi_p_thr_min = 1200.;
  double mom_pi_p_thr_max = 1300.;

  for ( int i_ang_bin = 0; i_ang_bin < pid_ang_bin_vec_.size(); i_ang_bin++ ) {

    auto ang_bin = pid_ang_bin_vec_.at(i_ang_bin);

    // pion と proton の sigma が同じとみなされるときの　VPH

    auto pid_mom_bin_vec = pid_mom_bin_map_.at(ang_bin);
    auto param = vph_func_param_map_.at(ang_bin);

    for ( int i_mom_bin = 0; i_mom_bin < pid_mom_bin_vec.size(); i_mom_bin++) {

      auto mom_bin = pid_mom_bin_vec.at(i_mom_bin);
      auto data = likelihood_param_bin_map_.at(std::make_pair(i_ang_bin, i_mom_bin));

      double mom_bin_mean = (mom_bin.first + mom_bin.second) / 2.;
      if ( mom_pi_p_thr_min < mom_bin_mean &&
	   mom_bin_mean < mom_pi_p_thr_max ) { // 1200--1300 MeV/c の間で cross するはず
	for ( auto pbeta = mom_bin.first; pbeta <= mom_bin.second; pbeta += 1.0 ) {
	  double vph = CalcMomentumVphEmulsionFit(pbeta / 1000., ang_bin);
	  double sigma = param.sigma_inter + param.sigma_scale * TMath::Sqrt(vph);
	  // bin の中での vph から計算した sigma が bin の sigma と一致したところを threshold とする
	  if ( sigma - data.sigma[1] < 1. ) {
	    vph_sigma_cross_map_.insert(std::make_pair(ang_bin, vph));
	    BOOST_LOG_TRIVIAL(trace) << "Angle : (" << ang_bin.first << ", " << ang_bin.second << ")"
				     << " -> Threshold sigma VPH : " << vph;
	  }
	  break;
	} // pbeta
	break;
      } // fi mom
    } // imom
  } // iang

  return;

}


// pid_ang_bin_vec_ で何番目の bin に入るかを確認
// 最大値より大きい場合は pid_ang_bin_vec_.size() とする
int PidFunction::GetPidAngBinId(double tangent) const {
  
  for ( int i = 0; i < pid_ang_bin_vec_.size(); i++ ) {
    auto pid_ang_bin_ = pid_ang_bin_vec_.at(i);
    if ( tangent >= pid_ang_bin_.first &&
	 tangent < pid_ang_bin_.second ) {
      return i;
    }
  }

  return pid_ang_bin_vec_.size();
  
}

// 与えられた particle id, pbeta, tangent をもとに
// data driven な mean/sigma を使って gauss で VPH を計算
double PidFunction::GetVph(int true_particle_id,
			   double ecc_mcs_mom,
			   double tangent) const {

  // tangent と pbeta に従って VPH の mean と sigma を求める
  double mean = -1;
  double sigma = 0.;
  if ( true_particle_id == 2212 ) {
    mean = CalcVphProton(ecc_mcs_mom, tangent);
    sigma = CalcVphSigmaProton(ecc_mcs_mom, tangent);
  }
  else {
    mean = CalcVphMuon(ecc_mcs_mom, tangent);
    sigma = CalcVphSigmaMuon(ecc_mcs_mom, tangent);
  }

  double vph = gRandom->Gaus(mean, sigma);
  if ( vph < 0 ) vph = 0.;
  return vph;

}

// 与えられた vph, pbeta, tangent をもとに
// likelihood を計算 (何 sigma 離れているか)
void PidFunction::CalcPartnerLikelihood(double vph, double ecc_mcs_mom, double tangent,
					double &muon_likelihood,
					double &proton_likelihood) const {
  
  // Reconstructed pbeta に対応する 各粒子の mean/sigma を求める
  double muon_mean = CalcVphMuon(ecc_mcs_mom, tangent);
  double proton_mean = CalcVphProton(ecc_mcs_mom, tangent);
  double muon_sigma = CalcVphSigmaMuon(ecc_mcs_mom, tangent);
  double proton_sigma = CalcVphSigmaProton(ecc_mcs_mom, tangent);

  // 実際の VPH をもとにlikelihood を計算
  muon_likelihood = (vph - muon_mean) / muon_sigma;
  proton_likelihood = (proton_mean - vph) / proton_sigma;
  /*
  muon_likelihood = TMath::Gaus(vph, muon_mean,
				muon_sigma, kTRUE);
  proton_likelihood = TMath::Gaus(vph, proton_mean,
				  proton_sigma, kTRUE);
  */
  return;
}

// pbeta, tangent から muon vph を計算する
double PidFunction::CalcVphMuon(double pbeta, double tangent) const {

  int ang_bin_id = GetPidAngBinId(tangent);

  // 最大値より大きな角度の場合は最大の bin の値をそのまま使う
  if ( ang_bin_id == pid_ang_bin_vec_.size() ) {
    ang_bin_id = ang_bin_id - 1;
  }

  // ang_bin_id に対応する ang_bin に対応する mom_bin の vector
  auto pid_mom_bin = pid_mom_bin_map_.at(pid_ang_bin_vec_.at(ang_bin_id));
  // 最小 bin の平均より小さい場合は最小の bin の値をそのまま使う
  if ( pbeta < (pid_mom_bin.front().first + pid_mom_bin.front().second) / 2. ) {
    auto param = likelihood_param_bin_map_.at(std::make_pair(ang_bin_id, 0));
    return param.mean[1];
  }

  // mom bin の 平均同士を比べて内挿に使う値を決定する
  double vph_min = -1.;
  double vph_max = -1.;
  double pb_min = -1.;
  double pb_max = -1.;

  for ( int i = 0; i < pid_mom_bin.size(); i++ ) {
    if ( i == 0 ) continue; // すでに return されているが一応
    auto param = likelihood_param_bin_map_.at(std::make_pair(ang_bin_id, i));
    double input_mom_mean = (param.input_mom_min + param.input_mom_max) / 2.;
    if ( pbeta < input_mom_mean ) { // 初めて bin の平均より小さくなったときが正しい bin
      auto param_prev = likelihood_param_bin_map_.at(std::make_pair(ang_bin_id, i-1));
      vph_min = param_prev.mean[1];
      vph_max = param.mean[1];
      pb_min = (param_prev.input_mom_min + param_prev.input_mom_max) / 2.;
      pb_max = input_mom_mean;
      break;
    }
  }

  if ( vph_min < 0. ) { // 見つからなかったら最後の bin の値をそのまま使う
    int mom_bin_id = pid_mom_bin.size() - 1;
    return likelihood_param_bin_map_.at(std::make_pair(ang_bin_id, mom_bin_id)).mean[1];
  }

  // 直線で内挿して決定
  return vph_min + (vph_max - vph_min) / (pb_max - pb_min) * (pbeta - pb_min);

}


// pbeta tangent から muon vph sigma を計算する
double PidFunction::CalcVphSigmaMuon(double pbeta, double tangent) const {
  int ang_bin_id = GetPidAngBinId(tangent);

  // 最大値より大きな角度の場合は最大の bin の値をそのまま使う
  if ( ang_bin_id == pid_ang_bin_vec_.size() ) {
    ang_bin_id = ang_bin_id - 1;
  }

  // ang_bin_id に対応する ang_bin に対応する mom_bin の vector
  auto pid_mom_bin = pid_mom_bin_map_.at(pid_ang_bin_vec_.at(ang_bin_id));

  // bin の中で sigma は一定
  for ( int i = 0; i < pid_mom_bin.size(); i++ ) {
    auto param = likelihood_param_bin_map_.at(std::make_pair(ang_bin_id, i));
    if ( pbeta < param.input_mom_max ) { // 初めて mom_max より小さくなったところが正しい bin
      return param.sigma[1];
    }
  }

  // 最後の bin の値をそのまま使う
  int mom_bin_id = pid_mom_bin.size() - 1;
  auto param = likelihood_param_bin_map_.at(std::make_pair(ang_bin_id, mom_bin_id));  
  return param.sigma[1];
  
}


// pbeta tangent から proton vph sigma を計算する
double PidFunction::CalcVphProton(double pbeta, double tangent) const {

  int ang_bin_id = GetPidAngBinId(tangent);

  // 最大値より大きな角度の場合は最大の bin の値をそのまま使う
  if ( ang_bin_id == pid_ang_bin_vec_.size() ) {
    ang_bin_id = ang_bin_id - 1;
  }

  auto ang_bin = pid_ang_bin_vec_.at(ang_bin_id);
  double vph = CalcMomentumVphEmulsionFit(pbeta / 1000., ang_bin);
  // 200 MeV/c に対応する点より値が小さかったら fit した値を使う
  if ( vph < vph_mean_cross_map_.at(ang_bin))
    return vph;

  auto pid_mom_bin = pid_mom_bin_map_.at(ang_bin);
  // 最小 bin の平均より小さい場合は最小の bin の値をそのまま使う
  if ( pbeta < (pid_mom_bin.front().first + pid_mom_bin.front().second) / 2. ) {
    auto param = likelihood_param_bin_map_.at(std::make_pair(ang_bin_id, 0));
    return param.mean[2];
  }

  // mom bin の平均同士を比べて内挿に使う値を決定する
  double vph_min = -1.;
  double vph_max = -1.;
  double pb_min = -1.;
  double pb_max = -1.;

  for ( int i = 0; i < pid_mom_bin.size(); i++ ) {
    if ( i == 0 ) continue; // すでに return されているが一応
    auto param = likelihood_param_bin_map_.at(std::make_pair(ang_bin_id, i));
    double input_mom_mean = (param.input_mom_min + param.input_mom_max) / 2.;
    if ( pbeta < input_mom_mean ) { // 初めて bin の平均より小さくなったときが正しい bin
      auto param_prev = likelihood_param_bin_map_.at(std::make_pair(ang_bin_id, i-1));
      vph_min = param_prev.mean[2];
      vph_max = param.mean[2];
      pb_min = (param_prev.input_mom_min + param_prev.input_mom_max) / 2.;
      pb_max = input_mom_mean;

      // pb max が 200 MeV/c より大きくなったら範囲を変更
      if ( pb_max > 200. ) {
	pb_max = 200.; vph_max = vph_mean_cross_map_.at(ang_bin);
      }

      break;
    }
  }


  // 直線で内挿して決定
  return vph_min + (vph_max - vph_min) / (pb_max - pb_min) * (pbeta - pb_min);

}

double PidFunction::CalcVphSigmaProton(double pbeta, double tangent) const {

  int ang_bin_id = GetPidAngBinId(tangent);

  // 最大値より大きな角度の場合は最大の bin の値をそのまま使う
  if ( ang_bin_id == pid_ang_bin_vec_.size() ) {
    ang_bin_id = ang_bin_id - 1;
  }

  auto ang_bin = pid_ang_bin_vec_.at(ang_bin_id);
  auto func_param = vph_func_param_map_.at(ang_bin);
  
  double vph = CalcMomentumVphEmulsionFit(pbeta / 1000., ang_bin);

  // low momentum はデータをそのまま使う
  if ( vph > vph_mean_cross_map_.at(ang_bin) ) {
    // ang_bin_id に対応する ang_bin に対応する mom_bin の vector
    auto pid_mom_bin = pid_mom_bin_map_.at(pid_ang_bin_vec_.at(ang_bin_id));
    
    // bin の中で sigma は一定
    for ( int i = 0; i < pid_mom_bin.size(); i++ ) {
      auto param = likelihood_param_bin_map_.at(std::make_pair(ang_bin_id, i));
      if ( pbeta < param.input_mom_max ) { // 初めて mom_max より小さくなったところが正しい bin
	return param.sigma[2];
      }
    }    
  }
  // pi/p が分けられない点では pi をそのまま延長
  else if ( vph < vph_sigma_cross_map_.at(ang_bin) ) {
    vph = vph_sigma_cross_map_.at(ang_bin);
  }

  return func_param.sigma_inter + func_param.sigma_scale * TMath::Sqrt(vph);

}




double PidFunction::CalcMomentumVphEmulsionFit(double pbeta, std::pair<double, double > ang_bin) const {

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

  auto param = vph_func_param_map_.at(ang_bin);
  double slope = param.mean_slope;
  double intercept = param.mean_inter;
  
  return slope * dedx + intercept;

}


int PidFunction::GetReconPid(double muon_likelihood, double proton_likelihood) const {
  
  if ( muon_likelihood < 0. ) {
    return 211;
  }
  else if ( proton_likelihood < 0. ) {
    return 2212;
  }
  else {
    double likelihood_ratio = muon_likelihood / (muon_likelihood + proton_likelihood);
    if ( likelihood_ratio < 0.5 ) {
      return 211;
    }
    else {
      return 2212;
    }
  }
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

void PidFunction::CheckMeanSigmaValues() const {

  for ( double pbeta = 1.; pbeta < 1500.; pbeta++ ) {
    double proton_mean = CalcVphProton(pbeta, 0.);
    double proton_sigma = CalcVphSigmaProton(pbeta, 0.);
    std::cout << "Momentum : " << pbeta << ", "
	      << "Mean : " << proton_mean << ", "
	      << "Sigma : " << proton_sigma << std::endl;
  }

  std::exit(1);

}
