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
    throw std::runtime_error("File not fund : " + filename);
  
  auto ev_vec = Momentum_recon::ReadEventInformationBin(filename);

  TFile *ofile = new TFile(argv[2], "recreate");
  TH1D *hist_muon_mom = new TH1D("hist_muon_mom", "", muon_mom_bin_size-1, muon_mom_bins);
  TH1D *hist_muon_mom_mcs = new TH1D("hist_muon_mom_mcs", "", muon_mom_bin_size-1, muon_mom_bins);
  TH1D *hist_muon_mom_range = new TH1D("hist_muon_mom_range", "", muon_mom_bin_size-1, muon_mom_bins);
  TH1D *hist_proton_mom = new TH1D("hist_proton_mom", "", hadron_mom_bin_size-1, hadron_mom_bins);
  TH1D *hist_proton_mom_mcs = new TH1D("hist_proton_mom_mcs", "", hadron_mom_bin_size-1, hadron_mom_bins);
  TH1D *hist_proton_mom_range = new TH1D("hist_proton_mom_range", "", hadron_mom_bin_size-1, hadron_mom_bins);
  TH1D *hist_muon_ang_deg = new TH1D("hist_muon_ang_deg", "", muon_deg_bin_size-1, muon_deg_bins);
  TH1D *hist_muon_ang_cos = new TH1D("hist_muon_ang_cos", "", muon_cos_bin_size-1, muon_cos_bins);
  TH1D *hist_proton_ang_deg = new TH1D("hist_proton_ang_deg", "", hadron_deg_bin_size-1, hadron_deg_bins);
  TH1D *hist_proton_ang_cos = new TH1D("hist_proton_ang_cos", "", hadron_cos_bin_size-1, hadron_cos_bins);

  TH1D *hist_q2 = new TH1D("hist_q2","", q2_bin_size-1, q2_bins);
  TH1D *hist_q2_range = new TH1D("hist_q2_range", "", q2_bin_size-1, q2_bins);
  TH1D *hist_nu_ene_recon = new TH1D("hist_nu_ene_recon", "", nu_ene_recon_bin_size-1, nu_ene_recon_bins);
  TH1D *hist_nu_ene_recon_range = new TH1D("hist_nu_ene_recon_range", "", nu_ene_recon_bin_size-1, nu_ene_recon_bins);
  
  TH1D *hist_dpt = new TH1D("hist_dpt", "", dpt_bin_size-1, dpt_bins);
  TH1D *hist_dpt_range = new TH1D("hist_dpt_range", "", dpt_bin_size-1, dpt_bins);
  TH1D *hist_dalphat = new TH1D("hist_dalphat", "", dalphat_bin_size-1, dalphat_bins);
  TH1D *hist_dalphat_range = new TH1D("hist_dalphat_range", "", dalphat_bin_size-1, dalphat_bins);
  TH1D *hist_cosdat = new TH1D("hist_cosdat", "", cosdat_bin_size-1, cosdat_bins);
  TH1D *hist_cosdat_range = new TH1D("hist_cosdat_range", "", cosdat_bin_size-1, cosdat_bins);
  TH1D *hist_dphit = new TH1D("hist_dphit", "", dphit_bin_size-1, dphit_bins);
  TH1D *hist_dphit_range = new TH1D("hist_dphit_range", "", dphit_bin_size-1, dphit_bins);
  TH1D *hist_cosdphit = new TH1D("hist_cosdphit", "", cosdphit_bin_size-1, cosdphit_bins);
  TH1D *hist_cosdphit_range = new TH1D("hist_cosdphit_range", "", cosdphit_bin_size-1, cosdphit_bins);
  TH1D *hist_dptx = new TH1D("hist_dptx", "", dptx_bin_size-1, dptx_bins);
  TH1D *hist_dptx_range = new TH1D("hist_dptx_range", "", dptx_bin_size-1, dptx_bins);
  TH1D *hist_dpty = new TH1D("hist_dpty", "", dpty_bin_size-1, dpty_bins);
  TH1D *hist_dpty_range = new TH1D("hist_dpty_range", "", dpty_bin_size-1, dpty_bins);

  for ( auto ev : ev_vec ) {
    
    int num_muon = 0;
    int num_proton = 0;
    int num_nonproton = 0;

    std::cout << "Group : " << ev.groupid << std::endl;
    std::cout << "Entry : " << ev.entry_in_daily_file << std::endl;
    std::cout << "Timestamp : " << ev.unixtime << std::endl;

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
	   (pb_mu_assume < 700. || vph > 125) ) {
	num_proton++;
      }
      else {
	num_nonproton++;
      }
    }

    if ( num_muon == 1 && num_proton == 1 && num_nonproton == 0 ) {
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

    bool muon_stop_flag = false;
    bool proton_stop_flag = false;
    TVector3 muon_tangent;
    TVector3 proton_tangent;
    int proton_direction = 0;
    double muon_momentum = -1.;
    double proton_momentum = -1.;
    double muon_theta_deg = 0.;
    double proton_theta_deg = 0.;
    double muon_cosine = 0.;
    double proton_cosine = 0.;
    TVector3 muon_momentum_vec;
    TVector3 proton_momentum_vec;

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
	else if ( chain.stop_flag == 0 ) {
	  muon_stop_flag = false;
	  muon_momentum = chain.ecc_mcs_mom[0];
	}
	muon_momentum_vec = (muon_momentum / muon_tangent.Mag()) * muon_tangent;
      }
      else {
	proton_direction = chain.direction;
	proton_tangent.SetXYZ(ax, ay, 1.);
	proton_theta_deg = theta_deg;
	proton_cosine = cosine;
	if ( chain.stop_flag == 1 ) {
	  proton_stop_flag = true;
	  proton_momentum = chain.ecc_range_mom[1];
	}
	else {
	  proton_stop_flag = false;
	  proton_momentum = chain.ecc_mcs_mom[1];
	}
	proton_momentum_vec = (proton_momentum * proton_direction / proton_tangent.Mag()) * proton_tangent;
      }

    }

    // Calculate
    double muon_energy = std::sqrt(muon_mass * muon_mass + muon_momentum * muon_momentum);
    double recon_nu_ene = (neutron_mass * muon_energy - muon_mass * muon_mass / 2.
			   + (proton_mass * proton_mass - neutron_mass * neutron_mass) / 2.)
      / (neutron_mass - muon_energy + muon_momentum * muon_cosine);
    double recon_q2 = - muon_mass * muon_mass + 2 * recon_nu_ene * (muon_energy - muon_momentum * muon_cosine);
    
    TVector2 mumom_vec_2d(muon_momentum_vec.X(), muon_momentum_vec.Y());
    TVector2 promom_vec_2d(proton_momentum_vec.X(), proton_momentum_vec.Y());
    
    TVector2 dpt_vec = mumom_vec_2d + promom_vec_2d;
    double dalphat = std::acos((-1. * mumom_vec_2d * dpt_vec) / mumom_vec_2d.Mod() / dpt_vec.Mod()) * TMath::RadToDeg();
    double cosdat = (-1. * mumom_vec_2d * dpt_vec) / mumom_vec_2d.Mod() / dpt_vec.Mod();
    double dphit = std::acos((-1. * mumom_vec_2d * promom_vec_2d) / mumom_vec_2d.Mod() / promom_vec_2d.Mod()) * TMath::RadToDeg();
    double cosdphit = (-1. * mumom_vec_2d * promom_vec_2d) / mumom_vec_2d.Mod() / promom_vec_2d.Mod();
    double dptx = TMath::Sign(1., mumom_vec_2d.X()) * dpt_vec.Mod() * std::sin(dalphat);
    double dpty = dpt_vec.Mod() * cosdat;
    
    // Fill histograms
    hist_muon_mom->Fill(muon_momentum);
    hist_proton_mom->Fill(proton_momentum);
    hist_muon_ang_deg->Fill(muon_theta_deg);
    hist_muon_ang_cos->Fill(muon_cosine);
    hist_proton_ang_deg->Fill(proton_theta_deg);
    hist_proton_ang_cos->Fill(proton_cosine);
    hist_q2->Fill(recon_q2);
    hist_nu_ene_recon->Fill(recon_nu_ene);
    hist_dpt->Fill(dpt_vec.Mod());
    hist_dalphat->Fill(dalphat);
    hist_cosdat->Fill(cosdat);
    hist_dphit->Fill(dphit);
    hist_cosdphit->Fill(cosdphit);
    hist_dptx->Fill(dptx);
    hist_dpty->Fill(dpty);
    
    if ( muon_stop_flag ) {
      hist_muon_mom_range->Fill(muon_momentum);
      hist_q2_range->Fill(recon_q2);
      hist_nu_ene_recon_range->Fill(recon_nu_ene);
      hist_dpt_range->Fill(dpt_vec.Mod());
      hist_dalphat_range->Fill(dalphat);
      hist_cosdat_range->Fill(cosdat);
      hist_dphit_range->Fill(dphit);
      hist_cosdphit_range->Fill(cosdphit);
      hist_dptx_range->Fill(dptx);
      hist_dpty_range->Fill(dpty);
    }
    else {
      hist_muon_mom_mcs->Fill(muon_momentum);
    }
    
    if ( proton_stop_flag ) {
      hist_proton_mom_range->Fill(proton_momentum);       
    }
    else {
      hist_proton_mom_mcs->Fill(proton_momentum);
    }
       
  }

  ofile->cd();

  hist_muon_mom->Write();
  hist_muon_mom_range->Write();
  hist_muon_mom_mcs->Write();
  hist_muon_ang_deg->Write();
  hist_muon_ang_cos->Write();
  hist_proton_mom->Write();
  hist_proton_mom_range->Write();
  hist_proton_mom_mcs->Write();
  hist_proton_ang_deg->Write();
  hist_proton_ang_cos->Write();
  hist_q2->Write();
  hist_q2_range->Write();
  hist_nu_ene_recon->Write();
  hist_nu_ene_recon_range->Write();
  hist_dpt->Write();
  hist_dpt_range->Write();
  hist_dalphat->Write();
  hist_dalphat_range->Write();
  hist_cosdat->Write();
  hist_cosdat_range->Write();
  hist_dphit->Write();
  hist_dphit_range->Write();
  hist_cosdphit->Write();
  hist_cosdphit_range->Write();
  hist_dptx->Write();
  hist_dptx_range->Write();
  hist_dpty->Write();
  hist_dpty_range->Write();
  ofile->Close();

  return 0;

}
