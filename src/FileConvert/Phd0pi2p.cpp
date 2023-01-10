#include <boost/log/core.hpp>
#include <boost/log/trivial.hpp>
#include <boost/log/expressions.hpp>
#include <boost/filesystem.hpp>

#include <iostream>
#include <cmath>
#include <string>

#include <TFile.h>
#include <TH1D.h>
#include <TMath.h>

#include "McsClass.hpp"
#include "McsFunction.hpp"
#include "/home/t2k/odagawa/NinjaMCStudy/lib/HistogramStyle.hpp"
#include "/home/t2k/odagawa/NinjaMCStudy/lib/DrawConst.hpp"

namespace fs = boost::filesystem;
namespace logging = boost::log;

int main(int argc, char* argv[]) {

  if ( argc != 3 ) {
    BOOST_LOG_TRIVIAL(error) << "Usage : " << argv[0]
			     << " <input momch filename> <output histogram root filename>";
    return 1;
  }

  const double neutrino_beam_ax = -2.325e-2;
  const double neutrino_beam_ay = -8.075e-2;
  const double neutrino_beam_thetax = std::atan(neutrino_beam_ax);
  const double neutrino_beam_thetay = std::atan(neutrino_beam_ay);

  std::string filename = argv[1];
  if ( !fs::exists(filename) )
    throw std::runtime_error("File not found : " + filename);

  auto ev_vec = Momentum_recon::ReadEventInformationBin(filename);

  TFile *ofile = new TFile(argv[2], "recreate");
  TH1D *hist_muon_mom = new TH1D("hist_muon_mom", "", muon_mom_bin_size-1, muon_mom_bins);
  TH1D *hist_muon_mom_mcs = new TH1D("hist_muon_mom_mcs", "", muon_mom_bin_size-1, muon_mom_bins);
  TH1D *hist_muon_mom_range = new TH1D("hist_muon_mom_range", "", muon_mom_bin_size-1, muon_mom_bins);
  TH1D *hist_proton_mom = new TH1D("hist_proton_mom", "", hadron_mom_bin_size-1, hadron_mom_bins);
  TH1D *hist_proton_mom_low = new TH1D("hist_proton_mom_low", "", hadron_mom_bin_size-1, hadron_mom_bins);
  TH1D *hist_proton_mom_high = new TH1D("hist_proton_mom_high", "", hadron_mom_bin_size-1, hadron_mom_bins);
  TH1D *hist_proton_mom_mcs = new TH1D("hist_proton_mom_mcs", "", hadron_mom_bin_size-1, hadron_mom_bins);
  TH1D *hist_proton_mom_low_mcs = new TH1D("hist_proton_mom_low_mcs", "", hadron_mom_bin_size-1, hadron_mom_bins);
  TH1D *hist_proton_mom_high_mcs = new TH1D("hist_proton_mom_high_mcs", "", hadron_mom_bin_size-1, hadron_mom_bins);
  TH1D *hist_proton_mom_range = new TH1D("hist_proton_mom_range", "", hadron_mom_bin_size-1, hadron_mom_bins);
  TH1D *hist_proton_mom_low_range = new TH1D("hist_proton_mom_low_range", "", hadron_mom_bin_size-1, hadron_mom_bins);
  TH1D *hist_proton_mom_high_range = new TH1D("hist_proton_mom_high_range", "", hadron_mom_bin_size-1, hadron_mom_bins);
  TH1D *hist_muon_ang_deg = new TH1D("hist_muon_ang_deg", "", muon_deg_bin_size-1, muon_deg_bins);
  TH1D *hist_muon_ang_cos = new TH1D("hist_muon_ang_cos", "", muon_cos_bin_size-1, muon_cos_bins);
  TH1D *hist_proton_ang_deg = new TH1D("hist_proton_ang_deg", "", hadron_deg_bin_size-1, hadron_deg_bins);
  TH1D *hist_proton_ang_deg_low = new TH1D("hist_proton_ang_deg_low", "", hadron_deg_bin_size-1, hadron_deg_bins);
  TH1D *hist_proton_ang_deg_high = new TH1D("hist_proton_ang_deg_high", "", hadron_deg_bin_size-1, hadron_deg_bins);
  TH1D *hist_proton_ang_cos = new TH1D("hist_proton_ang_cos", "", hadron_cos_bin_size-1, hadron_cos_bins);
  TH1D *hist_proton_ang_cos_low = new TH1D("hist_proton_ang_cos_low", "", hadron_cos_bin_size-1, hadron_cos_bins);
  TH1D *hist_proton_ang_cos_high = new TH1D("hist_proton_ang_cos_high", "", hadron_cos_bin_size-1, hadron_cos_bins);
  
  // Proton kinematics
  TH1D *hist_open_ang_deg = new TH1D("hist_open_ang_deg", "", open_deg_bin_size-1, open_deg_bins);
  TH1D *hist_open_ang_cos = new TH1D("hist_open_ang_cos", "", open_cos_bin_size-1, open_cos_bins);
  TH1D *hist_mom_vecsum = new TH1D("hist_mom_vecsum", "", mom_vecsum_bin_size-1, mom_vecsum_bins);
  TH1D *hist_mom_scasum = new TH1D("hist_mom_scasum", "", mom_scasum_bin_size-1, mom_scasum_bins);
  TH1D *hist_mom_ratio = new TH1D("hist_mom_ratio", "", mom_ratio_bin_size-1, mom_ratio_bins);

  // TKI
  TH1D *hist_dptt = new TH1D("hist_dptt", "", dptt_bin_size-1, dptt_bins);
  TH1D *hist_dptt_range = new TH1D("hist_dptt_range", "", dptt_bin_size-1, dptt_bins);
  TH1D *hist_dpt = new TH1D("hist_dpt", "", dpt_2p_bin_size-1, dpt_2p_bins);
  TH1D *hist_dpt_range = new TH1D("hist_dpt_range", "", dpt_2p_bin_size-1, dpt_2p_bins);
  TH1D *hist_pn = new TH1D("hist_pn", "", pn_bin_size-1, pn_bins);
  TH1D *hist_pn_range = new TH1D("hist_pn_range", "", pn_bin_size-1, pn_bins);
  TH1D *hist_dalphat = new TH1D("hist_dalphat", "", dalphat_2p_bin_size-1, dalphat_2p_bins);
  TH1D *hist_dalphat_range = new TH1D("hist_dalphat_range", "", dalphat_2p_bin_size-1, dalphat_2p_bins);
  TH1D *hist_cosdat = new TH1D("hist_cosdat", "", cosdat_2p_bin_size-1, cosdat_2p_bins);
  TH1D *hist_cosdat_range = new TH1D("hist_cosdat_range", "", cosdat_2p_bin_size-1, cosdat_2p_bins);

  for ( auto ev : ev_vec ) {

    int num_muon = 0;
    int num_proton = 0;
    int num_nonproton = 0;

    bool sample_flag = false;
    for ( auto chain : ev.chains ) {
      double pb_mu_assume = CalculatePBetaFromMomentum(chain.ecc_mcs_mom[0], MCS_MUON_MASS);
      
      double vph = 0.;
      int num_vph_base = 0;
      for ( auto base : chain.base ) {
	if ( base.m[0].zone % 10000 > 0 && base.m[1].zone % 10000 > 0 ) {
	  vph += (base.m[0].zone % 10000 + base.m[1].zone % 10000);
	  num_vph_base += 2;
	}
      }
      vph /= (double)num_vph_base;

      if ( chain.particle_flag == 13 ) {
	num_muon++;
      }
      else if ( chain.particle_flag == 2212 &&
		(pb_mu_assume < 700. || vph > 125.) ) {
	num_proton++;
      }
      else {
	num_nonproton++;
      }
    }

    std::cout << "muon : " <<  num_muon << ", proton : " << num_proton << ", other : " << num_nonproton << std::endl;
    if ( num_muon == 1 && num_proton == 2 && num_nonproton == 0 ) {
      sample_flag = true;
    }


    if ( !sample_flag ) continue;
    
    bool momentum_consistency_flag = true;
    
    for ( auto chain : ev.chains ) {
      if ( chain.particle_flag == 13 ) {
	double mcs_mom = chain.ecc_mcs_mom[0];
	double range_mom = chain.bm_range_mom;
	
	if ( range_mom < 0 ) {
	  momentum_consistency_flag = false;
	  break;
	}
	
	if ( chain.stop_flag == 0 ) {
	  if ( range_mom < mcs_mom ) break;
	  else if ( range_mom > 1000. ) break;
	  else {
	    double sigma = (range_mom - mcs_mom) / std::hypot(chain.bm_range_mom_error[0],
							      chain.ecc_mcs_mom_error[0][1]);
	    if ( sigma > 2.5 ) {
	      momentum_consistency_flag = false;		
	    }
	    break;
	  }
	}
	else {
	  if ( range_mom > 1000. ) break;
	  else if ( range_mom < mcs_mom ) {
	    double sigma = (mcs_mom - range_mom) / std::hypot(chain.bm_range_mom_error[1],
							      chain.ecc_mcs_mom_error[0][0]);
	    if ( sigma > 2.5 ) {
	      momentum_consistency_flag = false;
	    }
	    break;
	  }
	  else {
	    double sigma = (range_mom - mcs_mom) / std::hypot(chain.bm_range_mom_error[0],
							      chain.ecc_mcs_mom_error[0][1]);
	    if ( sigma > 2.5 ) {
	      momentum_consistency_flag = false;
	    }
	    break;
	  }
	}
      }
    }
    
    if ( !momentum_consistency_flag ) {
      continue;
    }

    std::cout << "Group : " << ev.groupid << std::endl;
    std::cout << "Entry : " << ev.entry_in_daily_file << std::endl;
    std::cout << "Timestamp : " << ev.unixtime << std::endl;
    
    num_proton = 0;
    bool muon_stop_flag = false;
    bool proton_low_stop_flag = false;
    bool proton_high_stop_flag = false;
    TVector3 muon_tangent;
    TVector3 proton_low_tangent;
    TVector3 proton_high_tangent;
    int proton_low_direction = 0;
    int proton_high_direction = 0;
    double muon_momentum = -1.;
    double proton_low_momentum = -1.;
    double proton_high_momentum = -1.;
    double muon_theta_deg = 0.;
    double proton_low_theta_deg = 0.;
    double proton_high_theta_deg = 0.;
    double muon_cosine = 0.;
    double proton_low_cosine = 0.;
    double proton_high_cosine = 0.;
    TVector3 muon_momentum_vec;
    TVector3 proton_low_momentum_vec;
    TVector3 proton_high_momentum_vec;
    
    for ( auto chain : ev.chains ) {
      
      double ax = 0.;
      double ay = 0.;
      if ( chain.direction == 1 ) {
	ax = chain.base.back().ax;
	ay = chain.base.back().ay;
      }
      else if ( chain.direction == -1 ) {
	ax = chain.base.front().ax;
	ay = chain.base.front().ay;
      }
      
      double thetax = std::atan(ax);
      double thetay = std::atan(ay);
      
      thetax -= neutrino_beam_thetax;
      thetay -= neutrino_beam_thetay;
      
      ax = std::tan(thetax);
      ay = std::tan(thetay);
      
      double tangent = chain.direction * std::hypot(ax, ay);
      double theta = std::atan(tangent);
      double theta_deg = theta * TMath::RadToDeg();
      if ( chain.direction == -1 ) theta_deg += 180.;
      double cosine = std::cos(theta_deg * TMath::DegToRad());
      
      if ( chain.particle_flag == 13 ) {
	muon_tangent.SetXYZ(ax, ay, 1.);
	muon_theta_deg = theta_deg;
	muon_cosine = cosine;
	if ( chain.stop_flag == 1 ) {
	  muon_stop_flag = true;
	  muon_momentum = chain.bm_range_mom;
	}
	else {
	  muon_stop_flag = false;
	  muon_momentum = chain.ecc_mcs_mom[0];
	}
	muon_momentum_vec = (muon_momentum / muon_tangent.Mag()) * muon_tangent;
      }
      else {
	num_proton++;
	if ( num_proton == 1 ) {
	  proton_low_direction = chain.direction;
	  proton_low_tangent.SetXYZ(ax, ay, 1.);
	  proton_low_theta_deg = theta_deg;
	  proton_low_cosine = cosine;
	  if ( chain.stop_flag == 1 ) {
	    proton_low_stop_flag = true;
	    proton_low_momentum = chain.ecc_range_mom[1];
	  }
	  else {
	    proton_low_stop_flag = false;
	    proton_low_momentum = chain.ecc_mcs_mom[1];
	  }
	  proton_low_momentum_vec = (proton_low_momentum * proton_low_direction / proton_low_tangent.Mag()) * proton_low_tangent;
	}
	else if ( num_proton == 2 ) {
	  proton_high_direction = chain.direction;
	  proton_high_tangent.SetXYZ(ax, ay, 1.);
	  proton_high_theta_deg = theta_deg;
	  proton_high_cosine = cosine;
	  if ( chain.stop_flag == 1 ) {
	    proton_high_stop_flag = true;
	    proton_high_momentum = chain.ecc_range_mom[1];
	  }
	  else {
	    proton_high_stop_flag = false;
	    proton_high_momentum = chain.ecc_mcs_mom[1];
	  }
	  proton_high_momentum_vec = (proton_high_momentum * proton_high_direction / proton_high_tangent.Mag()) * proton_high_tangent;
	}
      }
    }
    
    if ( proton_low_momentum > proton_high_momentum ) {
      std::swap(proton_low_direction, proton_high_direction);
      std::swap(proton_low_tangent, proton_high_tangent);
      std::swap(proton_low_theta_deg, proton_high_theta_deg);
      std::swap(proton_low_cosine, proton_high_cosine);
      std::swap(proton_low_stop_flag, proton_high_stop_flag);
      std::swap(proton_low_momentum, proton_high_momentum);
      std::swap(proton_low_momentum_vec, proton_high_momentum_vec);
    }
    
    // Calculate
    
    double open_cos = (proton_low_momentum_vec * proton_high_momentum_vec)
      / proton_low_momentum_vec.Mag() / proton_high_momentum_vec.Mag();
    double open_ang = std::acos(open_cos) * TMath::RadToDeg();
    double mom_vecsum = (proton_low_momentum_vec + proton_high_momentum_vec).Mag();
    double mom_scasum = proton_low_momentum_vec.Mag() + proton_high_momentum_vec.Mag();
    double mom_ratio = proton_low_momentum_vec.Mag() / proton_high_momentum_vec.Mag();
    
    TVector3 ztt(-1. * muon_momentum_vec.Y(), muon_momentum_vec.X(), 0.);
    ztt = (1. / std::hypot(muon_momentum_vec.X(), muon_momentum_vec.Y())) * ztt;
    double dptt = ztt * proton_high_momentum_vec + ztt * proton_low_momentum_vec;
    
    TVector2 muon_mom_vec_2d(muon_momentum_vec.X(), muon_momentum_vec.Y());
    TVector2 proton_mom_low_vec_2d(proton_low_momentum_vec.X(), proton_low_momentum_vec.Y());
    TVector2 proton_mom_high_vec_2d(proton_high_momentum_vec.X(), proton_high_momentum_vec.Y());
    TVector2 dpt_vec = muon_mom_vec_2d + proton_mom_low_vec_2d + proton_mom_high_vec_2d;
    
    double muon_energy = std::sqrt(muon_mass * muon_mass + muon_momentum * muon_momentum);
    double proton_high_energy = std::sqrt(proton_mass * proton_mass + proton_high_momentum * proton_high_momentum);
    double proton_low_energy = std::sqrt(proton_mass * proton_mass + proton_low_momentum * proton_low_momentum);
    double residual_mass = 23.; // MeV
    
    double pl = 0.5 * (proton_mass + muon_momentum_vec.Z() + proton_high_momentum_vec.Z() + proton_low_momentum_vec.Z()
		       - muon_energy - proton_high_energy - proton_low_energy)
      - 0.5 * (dpt_vec.Mod2() + residual_mass * residual_mass)
      / (proton_mass + muon_momentum_vec.Z() + proton_high_momentum_vec.Z() + proton_low_momentum_vec.Z()
	 - muon_energy - proton_high_energy - proton_low_energy);
    double pn = std::sqrt(pl * pl + dpt_vec.Mod2());
    
    double dalphat = std::acos((-1. * muon_mom_vec_2d * dpt_vec) / muon_mom_vec_2d.Mod() / dpt_vec.Mod()) * TMath::RadToDeg();
    double cosdat = (-1. * muon_mom_vec_2d * dpt_vec) / muon_mom_vec_2d.Mod() / dpt_vec.Mod();
    
    std::cout << "Muon momentum : " << muon_momentum << " [MeV/c], " << muon_stop_flag  << ", " << muon_theta_deg << " [deg], "
	      << "Proton momentum : " << proton_high_momentum << " [MeV/c], " << proton_high_stop_flag << ", " << proton_high_theta_deg << " [deg], "
	      << "Proton momentum : " << proton_low_momentum << " [MeV/c], " << proton_low_stop_flag << ", " << proton_low_theta_deg << " [deg]" << std::endl;

    // Fill histograms
    hist_muon_mom->Fill(muon_momentum);
    hist_proton_mom->Fill(proton_low_momentum);
    hist_proton_mom->Fill(proton_high_momentum);
    hist_proton_mom_low->Fill(proton_low_momentum);
    hist_proton_mom_high->Fill(proton_high_momentum);
    hist_muon_ang_deg->Fill(muon_theta_deg);
    hist_muon_ang_cos->Fill(muon_cosine);
    hist_proton_ang_deg->Fill(proton_low_theta_deg);
    hist_proton_ang_deg->Fill(proton_high_theta_deg);
    hist_proton_ang_deg_low->Fill(proton_low_theta_deg);
    hist_proton_ang_deg_high->Fill(proton_high_theta_deg);
    hist_proton_ang_cos->Fill(proton_low_cosine);
    hist_proton_ang_cos->Fill(proton_high_cosine);
    hist_proton_ang_cos_low->Fill(proton_low_cosine);
    hist_proton_ang_cos_high->Fill(proton_high_cosine);
    
    hist_open_ang_deg->Fill(open_ang);
    hist_open_ang_cos->Fill(open_cos);
    hist_mom_vecsum->Fill(mom_vecsum);
    hist_mom_scasum->Fill(mom_scasum);
    hist_mom_ratio->Fill(mom_ratio);
    
    hist_dptt->Fill(dptt);
    hist_dpt->Fill(dpt_vec.Mod());
    hist_pn->Fill(pn);
    hist_dalphat->Fill(dalphat);
    hist_cosdat->Fill(cosdat);
    
    if ( muon_stop_flag ) {
      hist_muon_mom_range->Fill(muon_momentum);
      hist_dptt_range->Fill(dptt);
      hist_dpt_range->Fill(dpt_vec.Mod());
      hist_pn_range->Fill(pn);
      hist_dalphat_range->Fill(dalphat);
      hist_cosdat_range->Fill(cosdat);      
    }
    else {
      hist_muon_mom_mcs->Fill(muon_momentum);
    }
    
    if ( proton_low_stop_flag ) {
      hist_proton_mom_low_range->Fill(proton_low_momentum);      
    }
    else {
      hist_proton_mom_low_mcs->Fill(proton_low_momentum);
    }
    
    if ( proton_high_stop_flag ) {
      hist_proton_mom_high_range->Fill(proton_high_momentum);      
    }
    else {
      hist_proton_mom_high_mcs->Fill(proton_high_momentum);
    }
    
  }
  
  ofile->cd();
  hist_muon_mom->Write();
  hist_muon_mom_range->Write();
  hist_muon_mom_mcs->Write();
  hist_proton_mom->Write();
  hist_proton_mom_range->Write();
  hist_proton_mom_mcs->Write();
  hist_proton_mom_low->Write();
  hist_proton_mom_low_range->Write();
  hist_proton_mom_low_mcs->Write();
  hist_proton_mom_high->Write();
  hist_proton_mom_high_range->Write();
  hist_proton_mom_high_mcs->Write();
  hist_muon_ang_deg->Write();
  hist_muon_ang_cos->Write();
  hist_proton_ang_deg->Write();
  hist_proton_ang_cos->Write();
  hist_proton_ang_deg_low->Write();
  hist_proton_ang_deg_high->Write();
  hist_proton_ang_cos_low->Write();
  hist_proton_ang_cos_high->Write();
  hist_open_ang_deg->Write();
  hist_open_ang_cos->Write();
  hist_mom_vecsum->Write();
  hist_mom_scasum->Write();
  hist_mom_ratio->Write();
  hist_dptt->Write();
  hist_dptt_range->Write();
  hist_dpt->Write();
  hist_dpt_range->Write();
  hist_pn->Write();
  hist_pn_range->Write();
  hist_dalphat->Write();
  hist_dalphat_range->Write();
  hist_cosdat->Write();
  hist_cosdat_range->Write();
  ofile->Close();
  
  return 0;
  
}
