// PhD thesis のためにヒストグラムを作る

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

namespace fs = boost::filesystem;

int main(int argc, char* argv[]) {

  if ( argc != 3 ) {
    BOOST_LOG_TRIVIAL(error) << "Usage : " << argv[0]
			     << " <input momch filename> <output graph root filename>";
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
  // PID systematics
  TGraphAsymmErrors *g_proton_pid = new TGraphAsymmErrors();
  TGraphAsymmErrors *g_pion_pid = new TGraphAsymmErrors();
  TGraphAsymmErrors *g_other_pid = new TGraphAsymmErrors();
  TGraphAsymmErrors *g_proton_pid_nocut = new TGraphAsymmErrors();
  TGraphAsymmErrors *g_pion_pid_nocut = new TGraphAsymmErrors();
  g_proton_pid->SetMarkerStyle(20);
  g_pion_pid->SetMarkerStyle(20);
  g_other_pid->SetMarkerStyle(20);
  g_proton_pid_nocut->SetMarkerStyle(25);
  g_pion_pid_nocut->SetMarkerStyle(25);
  g_proton_pid->SetMarkerColor(kRed);
  g_pion_pid->SetMarkerColor(kBlue);
  g_other_pid->SetMarkerColor(kGreen);
  g_proton_pid_nocut->SetMarkerColor(kMagenta);
  g_pion_pid_nocut->SetMarkerColor(kCyan);
  int ip = 0;
  int ipi = 0;
  int iother = 0;
  int ip_nocut = 0;
  int ipi_nocut = 0;
  // Momentum consistency

  TGraphAsymmErrors *g_muon_mom_consistency = new TGraphAsymmErrors();
  TGraphAsymmErrors *g_muon_mom_consistency_stop = new TGraphAsymmErrors();
  TGraphAsymmErrors *g_muon_mom_consistency_nonstop = new TGraphAsymmErrors();
  TGraphAsymmErrors *g_muon_mom_consistency_remain = new TGraphAsymmErrors();
  TGraphAsymmErrors *g_muon_mom_consistency_cut = new TGraphAsymmErrors();
  TGraphAsymmErrors *g_muon_mom_consistency_stop_remain = new TGraphAsymmErrors();
  TGraphAsymmErrors *g_muon_mom_consistency_stop_cut = new TGraphAsymmErrors();
  TGraphAsymmErrors *g_muon_mom_consistency_nonstop_remain = new TGraphAsymmErrors();
  TGraphAsymmErrors *g_muon_mom_consistency_nonstop_cut = new TGraphAsymmErrors();
  g_muon_mom_consistency->SetMarkerStyle(20);
  g_muon_mom_consistency_stop->SetMarkerStyle(20);
  g_muon_mom_consistency_nonstop->SetMarkerStyle(20);
  g_muon_mom_consistency_remain->SetMarkerStyle(24);
  g_muon_mom_consistency_cut->SetMarkerStyle(20);
  g_muon_mom_consistency_stop_remain->SetMarkerStyle(24);
  g_muon_mom_consistency_stop_cut->SetMarkerStyle(20);
  g_muon_mom_consistency_nonstop_remain->SetMarkerStyle(24);
  g_muon_mom_consistency_nonstop_cut->SetMarkerStyle(20);
  g_muon_mom_consistency->SetMarkerColor(kRed);
  g_muon_mom_consistency_stop->SetMarkerColor(kBlue);
  g_muon_mom_consistency_nonstop->SetMarkerColor(kGreen);
  g_muon_mom_consistency_remain->SetMarkerColor(kRed);
  g_muon_mom_consistency_cut->SetMarkerColor(kRed);
  g_muon_mom_consistency_stop_remain->SetMarkerColor(kBlue);
  g_muon_mom_consistency_stop_cut->SetMarkerColor(kBlue);
  g_muon_mom_consistency_nonstop_remain->SetMarkerColor(kGreen);
  g_muon_mom_consistency_nonstop_cut->SetMarkerColor(kGreen);
  int icons = 0;
  int icons_stop = 0;
  int icons_nonstop = 0;
  int icons_remain = 0;
  int icons_cut = 0;
  int icons_stop_remain = 0;
  int icons_stop_cut = 0;
  int icons_nonstop_remain = 0;
  int icons_nonstop_cut = 0;

  // vtxz
  TH1D *hist_vtx_pl = new TH1D("hist_vtx_pl", ";Vertex Plate;Entries",
			       70,0,70);
  
  
  for ( auto ev : ev_vec ) {
    std::cout << ev.groupid << " : " << ev.unixtime << ", " << ev.entry_in_daily_file << std::endl;
    hist_vtx_pl->Fill(ev.vertex_pl);

    bool momentum_consistency_flag = false;

    for ( auto chain : ev.chains ) {

      if ( chain.particle_flag == 13 ) {

	// Momentum consistency graphs

	double sigma = -9999.;
	if ( chain.ecc_mcs_mom[0] > chain.bm_range_mom )
	  sigma = (chain.ecc_mcs_mom[0] - chain.bm_range_mom) / std::hypot(chain.ecc_mcs_mom_error[0][0], chain.bm_range_mom_error[1]);
	else
	  sigma = (chain.bm_range_mom - chain.ecc_mcs_mom[0]) / std::hypot(chain.ecc_mcs_mom_error[0][1], chain.bm_range_mom_error[0]);

	if ( chain.stop_flag == 0 ) { // non stop
	  g_muon_mom_consistency_nonstop->SetPoint(icons_nonstop,
						   chain.bm_range_mom, chain.ecc_mcs_mom[0]);
	  g_muon_mom_consistency_nonstop->SetPointError(icons_nonstop,
							chain.bm_range_mom_error[0], chain.bm_range_mom_error[1],
							chain.ecc_mcs_mom_error[0][0], chain.ecc_mcs_mom_error[0][1]);
	  icons_nonstop++;

	  if ( chain.bm_range_mom < chain.ecc_mcs_mom[0] ) momentum_consistency_flag = true;
	  else if ( chain.bm_range_mom > 1.e3 ) momentum_consistency_flag = true;
	  else if ( sigma > 2.5 ) momentum_consistency_flag = false;
	  else momentum_consistency_flag = true;

	  if ( momentum_consistency_flag ) {
	    g_muon_mom_consistency_nonstop_remain->SetPoint(icons_nonstop_remain,
							    chain.bm_range_mom, chain.ecc_mcs_mom[0]);
	    g_muon_mom_consistency_nonstop_remain->SetPointError(icons_nonstop_remain,
								 chain.bm_range_mom_error[0], chain.bm_range_mom_error[1],
								 chain.ecc_mcs_mom_error[0][0], chain.ecc_mcs_mom_error[0][1]);
	    icons_nonstop_remain++;
	  }
	  else {
	    g_muon_mom_consistency_nonstop_cut->SetPoint(icons_nonstop_cut,
							 chain.bm_range_mom, chain.ecc_mcs_mom[0]);
	    g_muon_mom_consistency_nonstop_cut->SetPointError(icons_nonstop_cut,
							      chain.bm_range_mom_error[0], chain.bm_range_mom_error[1],
							      chain.ecc_mcs_mom_error[0][0], chain.ecc_mcs_mom_error[0][1]);
	    icons_nonstop_cut++;
	  }
	}
	else if ( chain.stop_flag == 1 ) { // stop
	  g_muon_mom_consistency_stop->SetPoint(icons_stop,
						chain.bm_range_mom, chain.ecc_mcs_mom[0]);
	  g_muon_mom_consistency_stop->SetPointError(icons_stop,
						     chain.bm_range_mom_error[0], chain.bm_range_mom_error[1],
						     chain.ecc_mcs_mom_error[0][0], chain.ecc_mcs_mom_error[0][1]);
	  icons_stop++;
	  
	  if ( chain.bm_range_mom > 1.e3 ) momentum_consistency_flag = true;
	  else if ( sigma > 2.5 ) momentum_consistency_flag = false;
	  else momentum_consistency_flag = true;
	  
	  if ( momentum_consistency_flag ) {
	    g_muon_mom_consistency_stop_remain->SetPoint(icons_stop_remain,
							 chain.bm_range_mom, chain.ecc_mcs_mom[0]);
	    g_muon_mom_consistency_stop_remain->SetPointError(icons_stop_remain,
							      chain.bm_range_mom_error[0], chain.bm_range_mom_error[1],
							      chain.ecc_mcs_mom_error[0][0], chain.ecc_mcs_mom_error[0][1]);
	    icons_stop_remain++;
	  }
	  else {
	    g_muon_mom_consistency_stop_cut->SetPoint(icons_stop_cut,
						      chain.bm_range_mom, chain.ecc_mcs_mom[0]);
	    g_muon_mom_consistency_stop_cut->SetPointError(icons_stop_cut,
							   chain.bm_range_mom_error[0], chain.bm_range_mom_error[1],
							   chain.ecc_mcs_mom_error[0][0], chain.ecc_mcs_mom_error[0][1]);
	    icons_stop_cut++;
	  }
	}

	g_muon_mom_consistency->SetPoint(icons, chain.bm_range_mom, chain.ecc_mcs_mom[0]);
	g_muon_mom_consistency->SetPointError(icons,
					      chain.bm_range_mom_error[0], chain.bm_range_mom_error[1],
					      chain.ecc_mcs_mom_error[0][0], chain.ecc_mcs_mom_error[0][1]);
	icons++;

	if ( momentum_consistency_flag ) {
	  g_muon_mom_consistency_remain->SetPoint(icons_remain,
						  chain.bm_range_mom, chain.ecc_mcs_mom[0]);
	  g_muon_mom_consistency_remain->SetPointError(icons_remain,
						       chain.bm_range_mom_error[0], chain.bm_range_mom_error[1],
						       chain.ecc_mcs_mom_error[0][0], chain.ecc_mcs_mom_error[0][1]);
	  icons_remain++;
	}
	else {
	  g_muon_mom_consistency_cut->SetPoint(icons_cut,
					       chain.bm_range_mom, chain.ecc_mcs_mom[0]);
	  g_muon_mom_consistency_cut->SetPointError(icons_cut,
						    chain.bm_range_mom_error[0], chain.bm_range_mom_error[1],
						    chain.ecc_mcs_mom_error[0][0], chain.ecc_mcs_mom_error[0][1]);
	  icons_cut++;
	}
	
	std::cout << "Stop flag : " << chain.stop_flag << ", "
		  << "Momentum by MCS : " << chain.ecc_mcs_mom[0] << " [MeV/c], "
		  << "Momentum by range : " << chain.bm_range_mom << " [MeV/c], "
		  << "Sigma : " << sigma << std::endl;
	continue;
      }

      if ( chain.base.size() < 2 ) continue;
      if ( !momentum_consistency_flag ) continue;

      double ax = -1.;
      double ay = -1.;
      if ( chain.direction == 1 ) {
	ax = chain.base.back().ax;
	ay = chain.base.back().ay;
      }
      else {
	ax = chain.base.front().ax;
	ay = chain.base.front().ay;
      }
      double tangent = std::hypot(ax, ay);
      // if ( tangent < 0.4 || tangent > 0.5 ) continue;
      
      double vph = 0.;
      double vph_dev = 0.;
      int num_vph_base = 0;
      for ( auto base : chain.base ) {
	if ( base.m[0].zone % 10000 > 0 && base.m[1].zone % 10000 > 0 ) {
	  vph += (base.m[0].zone % 10000 + base.m[1].zone % 10000);
	  vph_dev += (base.m[0].zone % 10000) * (base.m[0].zone % 10000)
	    + (base.m[1].zone % 10000) * (base.m[1].zone % 10000);
	  num_vph_base++;
	  num_vph_base++;
	}
      }
      vph /= (double)num_vph_base;
      vph_dev = vph_dev / (double)num_vph_base - vph * vph;
      vph_dev = std::sqrt(vph_dev);

      double pbeta = chain.ecc_mcs_mom[0];
      double pbeta_minus = chain.ecc_mcs_mom[0] - chain.ecc_mcs_mom_error[0][0];
      double pbeta_plus = chain.ecc_mcs_mom[0] + chain.ecc_mcs_mom_error[0][1];
      pbeta = CalculatePBetaFromMomentum(pbeta, MCS_MUON_MASS);
      pbeta_minus = CalculatePBetaFromMomentum(pbeta_minus, MCS_MUON_MASS);
      pbeta_plus = CalculatePBetaFromMomentum(pbeta_plus, MCS_MUON_MASS);

      if ( chain.particle_flag == 2212 ) {
	g_proton_pid_nocut->SetPoint(ip_nocut, pbeta, vph);
	g_proton_pid_nocut->SetPointError(ip_nocut, pbeta - pbeta_minus, pbeta_plus - pbeta, vph_dev, vph_dev);
	ip_nocut++;
	if ( pbeta < 700. || vph > 125 ) {
	  g_proton_pid->SetPoint(ip, pbeta, vph);
	  g_proton_pid->SetPointError(ip, pbeta - pbeta_minus, pbeta_plus - pbeta, vph_dev, vph_dev);
	  ip++;
	}
      }
      else if ( chain.particle_flag == 211 ) {
	g_pion_pid_nocut->SetPoint(ipi_nocut, pbeta, vph);
	g_pion_pid_nocut->SetPointError(ipi_nocut, pbeta - pbeta_minus, pbeta_plus - pbeta, vph_dev, vph_dev);
	ipi_nocut++;
	if ( pbeta < 700. && vph <= 125 ) {
	  g_pion_pid->SetPoint(ipi, pbeta, vph);
	  g_pion_pid->SetPointError(ipi, pbeta - pbeta_minus, pbeta_plus - pbeta, vph_dev, vph_dev);
	  ipi++;
	}
      }
      
      // non pid
      if ( pbeta > 700. && vph <= 125 ) {
	g_other_pid->SetPoint(iother, pbeta, vph);
	g_other_pid->SetPointError(iother, pbeta - pbeta_minus, pbeta_plus - pbeta, vph_dev, vph_dev);
	iother++;
      }

      std::cout << "PID : " << chain.particle_flag << ", "
		<< "Pbeta (mu) : " << pbeta << " - " << pbeta - pbeta_minus << " + " << pbeta_plus - pbeta << " [MeV/c] ,"
		<< "VPH : " << vph << " +/- " << vph_dev << ", "
		<< "Segment : " << chain.base.size() << " ,"
		<< "Stop flag : " << chain.stop_flag << std::endl;

    } // chain
  }

  ofile->cd();
  hist_vtx_pl->Write();
  g_proton_pid->Write("g_proton_pid");
  g_pion_pid->Write("g_pion_pid");
  g_other_pid->Write("g_other_pid");
  g_proton_pid_nocut->Write("g_proton_pid_nocut");
  g_pion_pid_nocut->Write("g_pion_pid_nocut");

  g_muon_mom_consistency->Write("g_muon_mom_consistency");
  g_muon_mom_consistency_stop->Write("g_muon_mom_consistency_stop");
  g_muon_mom_consistency_nonstop->Write("g_muon_mom_consistency_nonstop");
  g_muon_mom_consistency_remain->Write("g_muon_mom_consistency_remain");
  g_muon_mom_consistency_cut->Write("g_muon_mom_consistency_cut");
  g_muon_mom_consistency_stop_remain->Write("g_muon_mom_consistency_stop_remain");
  g_muon_mom_consistency_stop_cut->Write("g_muon_mom_consistency_stop_cut");
  g_muon_mom_consistency_nonstop_remain->Write("g_muon_mom_consistency_nonstop_remain");
  g_muon_mom_consistency_nonstop_cut->Write("g_muon_mom_consistency_nonstop_cut");


  ofile->Close();

  return 0;

}
