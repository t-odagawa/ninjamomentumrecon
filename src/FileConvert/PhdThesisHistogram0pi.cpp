// NuFact2022 のためにヒストグラムを作る

#include <boost/log/core.hpp>
#include <boost/log/trivial.hpp>
#include <boost/log/expressions.hpp>
#include <boost/filesystem.hpp>

#include <iostream>
#include <cmath>
#include <string>

#include <TFile.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TGraphAsymmErrors.h>
#include <TMath.h>

#include "McsClass.hpp"
#include "McsFunction.hpp"
#include "/home/t2k/odagawa/NinjaMCStudy/lib/HistogramStyle.hpp"

namespace fs = boost::filesystem;

int main (int argc, char* argv[]) {

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
  TH1D *hist_multi = new TH1D("hist_multi", "", total_multi_bin_size-1, total_multi_bins);
  TH1D *hist_multi_p = new TH1D("hist_multi_p", "", hadron_multi_bin_size-1, hadron_multi_bins);
  TH1D *hist_multi_pi = new TH1D("hist_multi_pi", "", hadron_multi_bin_size-1, hadron_multi_bins);
  TH1D *hist_multi_other = new TH1D("hist_multi_other", "", hadron_multi_bin_size-1, hadron_multi_bins);
  TH1D *hist_muon_mom = new TH1D("hist_muon_mom", "", muon_mom_bin_size-1, muon_mom_bins);
  TH1D *hist_muon_mom_mcs = new TH1D("hist_muon_mom_mcs", "", muon_mom_bin_size-1, muon_mom_bins);
  TH1D *hist_muon_mom_range = new TH1D("hist_muon_mom_range", "", muon_mom_bin_size-1, muon_mom_bins);
  TH1D *hist_pion_mom = new TH1D("hist_pion_mom", "", hadron_mom_bin_size-1, hadron_mom_bins);
  TH1D *hist_pion_mom_mcs = new TH1D("hist_pion_mom_mcs", "", hadron_mom_bin_size-1, hadron_mom_bins);
  TH1D *hist_pion_mom_range = new TH1D("hist_pion_mom_range", "", hadron_mom_bin_size-1, hadron_mom_bins);
  TH1D *hist_proton_mom = new TH1D("hist_proton_mom", "", hadron_mom_bin_size-1, hadron_mom_bins);
  TH1D *hist_proton_mom_mcs = new TH1D("hist_proton_mom_mcs", "", hadron_mom_bin_size-1, hadron_mom_bins);
  TH1D *hist_proton_mom_range = new TH1D("hist_proton_mom_range", "", hadron_mom_bin_size-1, hadron_mom_bins);
  TH1D *hist_muon_ang_deg = new TH1D("hist_muon_ang_deg", "", muon_deg_bin_size-1, muon_deg_bins);
  TH1D *hist_muon_ang_cos = new TH1D("hist_muon_ang_cos", "", muon_cos_bin_size-1, muon_cos_bins);
  TH1D *hist_pion_ang_deg = new TH1D("hist_pion_ang_deg", "", hadron_deg_bin_size-1, hadron_deg_bins);
  TH1D *hist_pion_ang_cos = new TH1D("hist_pion_ang_cos", "", hadron_cos_bin_size-1, hadron_cos_bins);
  TH1D *hist_proton_ang_deg = new TH1D("hist_proton_ang_deg", "", hadron_deg_bin_size-1, hadron_deg_bins);
  TH1D *hist_proton_ang_cos = new TH1D("hist_proton_ang_cos", "", hadron_cos_bin_size-1, hadron_cos_bins);

  TH1D *hist_muon_mom_mcs_all = new TH1D("hist_muon_mom_mcs_all", "", 10, 0, 2000.);
  TH1D *hist_muon_mom_mcs_fine = new TH1D("hist_muon_mom_mcs_fine", "", 10, 0, 2000.);

  TH2D *hist_muon_mom_range_mcs = new TH2D("hist_muon_mom_range_mcs","",20,0.,2000.,20,0.,2000.);
  TGraphAsymmErrors *g_muon_mom_range_mcs = new TGraphAsymmErrors();
  g_muon_mom_range_mcs->SetName("g_muon_mom_range_mcs");
  int ipoint = 0;

  TGraph *g_mu = new TGraph();
  TGraph *g_mu_mcs = new TGraph();
  TGraph *g_mu_range = new TGraph();
  TGraph *g_mu_mcs_all = new TGraph();
  g_mu->SetName("g_mu");
  g_mu_mcs->SetName("g_mu_mcs");
  g_mu_range->SetName("g_mu_range");
  g_mu_mcs_all->SetName("g_mu_mcs_all");
  g_mu->SetMarkerStyle(20);
  g_mu_mcs->SetMarkerStyle(20);
  g_mu_range->SetMarkerStyle(20);
  g_mu_mcs_all->SetMarkerStyle(20);
  g_mu->SetMarkerColor(kGreen);
  g_mu_mcs->SetMarkerColor(kMagenta);
  g_mu_range->SetMarkerColor(kRed);
  g_mu_mcs_all->SetMarkerColor(kRed);
  int imu = 0;
  int imu_mcs = 0;
  int imu_range = 0;

  TGraph *g_p = new TGraph();
  TGraph *g_p_mcs = new TGraph();
  TGraph *g_p_range = new TGraph();
  g_p->SetName("g_p");
  g_p_mcs->SetName("g_p_mcs");
  g_p_range->SetName("g_p_range");
  g_p->SetMarkerStyle(20);
  g_p_mcs->SetMarkerStyle(20);
  g_p_range->SetMarkerStyle(20);
  g_p->SetMarkerColor(kGreen);
  g_p_mcs->SetMarkerColor(kMagenta);
  g_p_range->SetMarkerColor(kRed);
  int ip = 0;
  int ip_mcs = 0;
  int ip_range = 0;

  TGraph *g_pi = new TGraph();
  TGraph *g_pi_mcs = new TGraph();
  TGraph *g_pi_range = new TGraph();
  g_pi->SetName("g_pi");
  g_pi_mcs->SetName("g_pi_mcs");
  g_pi_range->SetName("g_pi_range");
  g_pi->SetMarkerStyle(20);
  g_pi_mcs->SetMarkerStyle(20);
  g_pi_range->SetMarkerStyle(20);
  g_pi->SetMarkerColor(kGreen);
  g_pi_mcs->SetMarkerColor(kMagenta);
  g_pi_range->SetMarkerColor(kRed);
  int ipi = 0;
  int ipi_mcs = 0;
  int ipi_range = 0;
  
  for ( auto ev : ev_vec ) {
   
    int multiplicity = 0;
    int num_muon = 0;
    int num_proton = 0;
    int num_pion = 0;
    int num_other = 0;

    // select 0pi events
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
		//		chain.ecc_mcs_mom[0] < 700.) {
		(pb_mu_assume < 700. || vph > 125) ) {
	num_proton++; 
      }
      else if ( chain.particle_flag == 211 &&
		//		chain.ecc_mcs_mom[0] < 700. ) {
		(pb_mu_assume < 700. && vph <= 125) ) {
	num_pion++;
      } 
      else {
	num_other++;
      }
    }

    bool sample_flag = false;
    if ( num_muon == 1 && num_pion == 0 && num_other == 0 ) sample_flag =true;

    // if ( num_pion != 0 || num_other != 0 || num_proton != 2 ) continue;
    // if ( num_pion != 0 || num_other != 0 ) continue;
    if ( !sample_flag ) continue;
    num_proton = 0;
    num_pion = 0;
    num_other = 0;

    std::cout << "Group : " << ev.groupid << std::endl;
    std::cout << "Entry : " << ev.entry_in_daily_file << std::endl;
    std::cout << "Timestamp : " << ev.unixtime << std::endl;

    // Momentum consistency (from 2022.10.31)
    bool momentum_consistency_flag = true;
    for ( auto chain : ev.chains ) {
      if ( chain.particle_flag == 13 ) {
	double mcs_mom = chain.ecc_mcs_mom[0];
	double range_mom = chain.bm_range_mom;

	if ( range_mom < 0 ) {
	  momentum_consistency_flag = false;
	  break;
	}

	if ( chain.stop_flag == 0 ) { // BM non-stop のとき
	  if ( range_mom < mcs_mom ) break; // Range < MCS のときは momentum consistency はわからない
	  else if ( range_mom > 1000. ) break; // Range > 1 GeV/c のときは MCS saturation で consistency は難しい
	  else {
	    double sigma = (range_mom - mcs_mom) / std::hypot(chain.bm_range_mom_error[0],
							      chain.ecc_mcs_mom_error[0][1]);
	    if ( sigma > 2.5 ) {
	      std::cout << "================================" << std::endl;
	      std::cout << "Momentum consistency cut : " << ev.groupid << std::endl;
	      std::cout << "Not-Stop : MCS : " << mcs_mom << ", " << "Range : " << range_mom << std::endl;
	      std::cout << "Sigma : " << sigma << std::endl;
	      momentum_consistency_flag = false;
	    }
	    break;
	  }
	}
	else { // BM stop のとき
	  if ( range_mom > 1000. ) break; // Range > 1 GeV/c のときは MCS saturation で consistency は難しい
	  else if ( range_mom < mcs_mom ) {
	    double sigma = (mcs_mom - range_mom) / std::hypot(chain.bm_range_mom_error[1],
							      chain.ecc_mcs_mom_error[0][0]);
	    if ( sigma > 2.5 ) {
	      std::cout << "================================" << std::endl;
	      std::cout << "Momentum consistency cut : " << ev.groupid << std::endl;
	      std::cout << "Stop : MCS : " << mcs_mom << ", " << "Range : " << range_mom << std::endl;
	      std::cout << "Sigma : " << sigma << std::endl;
	      momentum_consistency_flag = false;
	    }

	    break;
	  }
	  else {
	    double sigma = (range_mom - mcs_mom) / std::hypot(chain.bm_range_mom_error[0],
							      chain.ecc_mcs_mom_error[0][1]);
	    if ( sigma > 2.5 ) {
	      std::cout << "================================" << std::endl;
	      std::cout << "Momentum consistency cut : " << ev.groupid << std::endl;
	      std::cout << "Stop : MCS : " << mcs_mom << ", " << "Range : " << range_mom << std::endl;
	      std::cout << "Sigma : " << sigma << std::endl;
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

    for ( auto chain : ev.chains ) {

      int pl_diff = chain.base.front().pl - chain.base.back().pl;
      if ( chain.base.size() == 2 && pl_diff == 1 ) continue;

      multiplicity++;

      double ax, ay;
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

      double tangent = chain.direction * std::hypot(std::tan(thetax), std::tan(thetay));
      double theta = std::atan(tangent);
      double theta_deg = theta * TMath::RadToDeg();
      if ( chain.direction == -1 ) theta_deg += 180.;
      double cosine = std::cos(theta_deg * TMath::DegToRad());

      double momentum = -1.;

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
	if ( chain.stop_flag == 1 ) {
	  momentum = chain.bm_range_mom;
	  hist_muon_mom_range->Fill(momentum);
	  hist_muon_mom_range_mcs->Fill(chain.bm_range_mom, chain.ecc_mcs_mom[0]);
	  g_muon_mom_range_mcs->SetPoint(ipoint, chain.bm_range_mom, chain.ecc_mcs_mom[0]);
	  g_muon_mom_range_mcs->SetPointError(ipoint,
					      chain.bm_range_mom_error[0], chain.bm_range_mom_error[1],
					      chain.ecc_mcs_mom_error[0][0], chain.ecc_mcs_mom_error[0][1]);
	  ipoint++;
	  g_mu_range->SetPoint(imu_range, momentum, theta_deg);
	  imu_range++;
	}
	else if ( chain.stop_flag == 0 ) {
	  momentum = chain.ecc_mcs_mom[0];
	  hist_muon_mom_mcs->Fill(momentum);
	  hist_muon_mom_mcs_fine->Fill(momentum);
	  g_mu_mcs->SetPoint(imu_mcs, momentum, theta_deg);
	  imu_mcs++;
	}
	
	std::cout << "Muon momentum : " << chain.ecc_mcs_mom[0] << ", " << chain.bm_range_mom << std::endl;
	std::cout << "Muon momentum : " << momentum << " MeV/c, " << "Stop flag : " << chain.stop_flag << std::endl;

	hist_muon_mom->Fill(momentum);
	hist_muon_ang_deg->Fill(theta_deg);
	hist_muon_ang_cos->Fill(cosine);

	hist_muon_mom_mcs_all->Fill(chain.ecc_mcs_mom[0]);

	g_mu->SetPoint(imu, momentum, theta_deg);
	g_mu_mcs_all->SetPoint(imu, chain.ecc_mcs_mom[0], theta_deg);
	imu++;

      }
      else if ( chain.particle_flag == 2212 &&
		//		chain.ecc_mcs_mom[0] < 700.) {
		(pb_mu_assume < 700. || vph > 125) ) {
	num_proton++; 
	if ( chain.stop_flag == 1 ) {
	  momentum = chain.ecc_range_mom[1];
	  hist_proton_mom_range->Fill(momentum);
	  g_p_range->SetPoint(ip_range, momentum, theta_deg);
	  ip_range++;
	}
	else {
	  momentum = chain.ecc_mcs_mom[1];
	  hist_proton_mom_mcs->Fill(momentum);
	  g_p_mcs->SetPoint(ip_mcs, momentum, theta_deg);
	  ip_mcs++;
	}

	std::cout << "Proton momentum : " << momentum << " MeV/c, " << "Stop flag : " << chain.stop_flag << std::endl;

	hist_proton_mom->Fill(momentum);
	hist_proton_ang_deg->Fill(theta_deg);
	hist_proton_ang_cos->Fill(cosine);

	g_p->SetPoint(ip, momentum, theta_deg);
	ip++;

      }
      else if ( chain.particle_flag == 211 &&
		//		chain.ecc_mcs_mom[0] < 700. ) {
		(pb_mu_assume < 700. && vph <= 125) ) {
	num_pion++;
	momentum = chain.bm_curvature_mom;
	hist_pion_mom_mcs->Fill(momentum);
	g_pi_mcs->SetPoint(ipi_mcs, momentum, theta_deg);
	ipi_mcs++;

	std::cout << "Pion momentum : " << momentum << " MeV/c, " << "Stop flag : " << chain.stop_flag << std::endl;
	if ( 200 < momentum && momentum < 400 && 100 < theta_deg && theta_deg< 120 )
	  std::cout << "Pion excess contribution ? : NSEG : " << chain.base.size() << ", Plate : " << ev.vertex_pl << std::endl;
	hist_pion_mom->Fill(momentum);
	hist_pion_ang_deg->Fill(theta_deg);
	hist_pion_ang_cos->Fill(cosine);
	
	g_pi->SetPoint(ipi, momentum, theta_deg);
	ipi++;
	
      }   
    } 

    std::cout << "Multiplicity : " << multiplicity
	      << " = Muon : 1 "
	      << " + Proton : " << num_proton
	      << " + Pion : " << num_pion 
	      << " + Other : " << multiplicity - num_proton - num_pion - 1 << std::endl;
    if  ( num_proton == 2 )  std::cout << "**********Group : " << ev.groupid << std::endl;
    hist_multi->Fill(multiplicity);
    hist_multi_p->Fill(num_proton);
    hist_multi_pi->Fill(num_pion);
    hist_multi_other->Fill(multiplicity - num_proton - num_pion - 1);
    
  }

  ofile->cd();
  hist_multi->Write();
  hist_multi_p->Write();
  hist_multi_pi->Write();
  hist_multi_other->Write();
  hist_muon_mom->Write();
  hist_muon_mom_mcs->Write();
  hist_muon_mom_range->Write();
  hist_proton_mom->Write();
  hist_proton_mom_mcs->Write();
  hist_proton_mom_range->Write();
  hist_pion_mom->Write();
  hist_pion_mom_mcs->Write();
  hist_pion_mom_range->Write();
  hist_muon_ang_deg->Write();
  hist_muon_ang_cos->Write();
  hist_proton_ang_deg->Write();
  hist_proton_ang_cos->Write();
  hist_pion_ang_deg->Write();
  hist_pion_ang_cos->Write();

  hist_muon_mom_mcs_fine->Write();
  hist_muon_mom_mcs_all->Write();

  hist_muon_mom_range_mcs->Write();
  g_muon_mom_range_mcs->Write();

  g_mu->Write();
  g_mu_mcs->Write();
  g_mu_range->Write();
  g_mu_mcs_all->Write();
  g_p->Write();
  g_p_mcs->Write();
  g_p_range->Write();
  g_pi->Write();
  g_pi_mcs->Write();
  g_pi_range->Write();

  ofile->Close();

  std::exit(0);

}
