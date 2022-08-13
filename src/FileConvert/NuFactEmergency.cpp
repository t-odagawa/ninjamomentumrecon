#include <boost/log/core.hpp>
#include <boost/log/trivial.hpp>
#include <boost/log/expressions.hpp>
#include <boost/filesystem.hpp>

#include <string>
#include <ctime>
#include <sstream>
#include <iostream>

#include <TFile.h>
#include <TTree.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TVector3.h>

#include <B2Const.hh>
#include <B2Enum.hh>

#include <NTBMSummary.hh>
#include <NTBMConst.hh>

#include "McsClass.hpp"
#include "McsFunction.hpp"
#include "MatchData.hpp"
#include "MatchFunction.hpp"
#include "RangeSpline.hpp"
#include "RangeFunction.hpp"

namespace logging = boost::log;
namespace fs = boost::filesystem;

int main ( int argc, char* argv[] ) {

  logging::core::get()->set_filter
    (
     //logging::trivial::severity >= logging::trivial::info
     logging::trivial::severity >= logging::trivial::debug
     );

  if ( argc != 7 ) {
    BOOST_LOG_TRIVIAL(error) << "Usage : " << argv[0]
			     << " <Input momch file> <NTBM prefix> <Old NTBM prefix> <Output momch file> <output root file> <Data dir path>";
    std::exit(1);
  }

  if ( !fs::exists((std::string)argv[1]) ) {
    throw std::runtime_error("File not found : " + (std::string)argv[1]);
  }
  auto ev_vec = Momentum_recon::ReadEventInformationBin((std::string)argv[1]);

  const std::string data_dir_path = argv[6];
  const MatchData match_data(data_dir_path);
  const MatchFunction match_function(match_data);
  const RangeSpline range_spline(data_dir_path);
  const RangeFunction range_function(range_spline);

  const double ecc_tracker_distance_x = 80.;
  const double ecc_tracker_distance_y = 60.;

  TH1D *hist_npl = new TH1D("hist_npl","", 18,0,18);
  TH1D *hist_range = new TH1D("hist_range", "", 20,0,2000.);
  TH2D *hist_npl_range = new TH2D("hist_npl_range", "", 18, 0, 18, 20, 0., 2000.);

  for ( auto &ev : ev_vec ) {

    BOOST_LOG_TRIVIAL(info) << "Entry : " << ev.entry_in_daily_file;
    BOOST_LOG_TRIVIAL(info) << "Timestamp : " << ev.unixtime;
    
    for ( auto &chain : ev.chains ) {

      if ( chain.particle_flag == 13 ) {
      
	time_t unixtime = (time_t)ev.unixtime;
	tm *tm_event = localtime(&unixtime);
	
	int year  = tm_event->tm_year + 1900;
	int month = tm_event->tm_mon + 1;
	int day   = tm_event->tm_mday;
	BOOST_LOG_TRIVIAL(debug) << year << "/" << month << "/" << day;

	std::stringstream filename_ss;
	filename_ss << argv[2] << "_" << year << "_" << month << "_" << day << ".root";
	if ( !fs::exists(filename_ss.str()) ) continue;
	TFile *file = new TFile(filename_ss.str().c_str(), "read");
	TTree *tree = (TTree*)file->Get("tree");
	NTBMSummary *ntbm = nullptr;
	tree->SetBranchAddress("NTBMSummary", &ntbm);

	std::stringstream old_filename_ss;
	old_filename_ss << argv[3] << "_" << year << "_" << month << "_" << day << ".root";
	if ( !fs::exists(old_filename_ss.str()) ) continue;
	TFile *old_file = new TFile(old_filename_ss.str().c_str(), "read");
	TTree *old_tree = (TTree*)old_file->Get("tree");
	NTBMSummary *old_ntbm = nullptr;
	old_tree->SetBranchAddress("NTBMSummary", &old_ntbm);

	tree->GetEntry(ev.entry_in_daily_file-1);
	old_tree->GetEntry(ev.entry_in_daily_file-1);

	int old_cluster_id = -1;
	for ( int i = 0; i < old_ntbm->GetNumberOfNinjaClusters(); i++ ) {
	  if ( old_ntbm->GetBabyMindTrackId(i) == ev.tracker_track_id ) {
	    old_cluster_id = i;
	    break;
	  }
	}
	auto old_position = old_ntbm->GetNinjaPosition(old_cluster_id);
	std::cout << old_position.at(0) << ", " << old_position.at(1) << std::endl;
	
	double dx = 300.;
	double dy = 200.;
	double dpos = std::hypot(dx, dy);	

	int muon_track_id = -1;

	for ( int i = 0; i < ntbm->GetNumberOfNinjaClusters(); i++ ) {
	  
	  if ( ntbm->GetBabyMindTrackId(i) < 0 ) continue;

	  auto position = ntbm->GetNinjaPosition(i);

	  double dx_tmp = position.at(B2View::kTopView) - old_position.at(B2View::kTopView);
	  double dy_tmp = position.at(B2View::kSideView) - old_position.at(B2View::kSideView);

	  if ( std::hypot(dx_tmp, dy_tmp) < dpos ) {
	    dx = dx_tmp, dy = dy_tmp;
	    dpos = std::hypot(dx, dy);
	    muon_track_id = ntbm->GetBabyMindTrackId(i);
	  }

	}

	std::cout << "Reconstructed position distance : (" << dx << ", " << dy << ")" << std::endl;

	double momentum = -1;
	double momentum_err = -1;

	if ( muon_track_id < 0 ) { // Not in new file...
	  chain.charge_sign = old_ntbm->GetCharge(ev.tracker_track_id);
	  momentum = old_ntbm->GetMomentum(ev.tracker_track_id);
	  if ( momentum > 0 ) chain.stop_flag = 1;
	  else chain.stop_flag = 0;
	}
	else {
	  chain.charge_sign = ntbm->GetCharge(muon_track_id);
	  if ( ntbm->GetMomentumType(muon_track_id) == 0 ) { // range
	    chain.stop_flag = 1;
	    
	    double range_mom = -1.;
	    double track_length = ntbm->GetTrackLengthTotal(muon_track_id); // g/cm2
	    TVector3 track_direction(chain.base.front().ax,
				     chain.base.front().ay, 1.);
	    int nwater = match_function.GetNumWaterPlate(ev.vertex_pl);
	    int niron = match_function.GetNumIronPlate(ev.vertex_pl);
	    track_length += match_function.GetDWGTrackLength(track_direction);
	    track_length += track_direction.Mag() * (.23 * nwater * 1.289 + 0.05 * niron * 8.03);
	    match_function.ConvertFromLengthToMom(range_mom, track_length);
	    momentum = range_mom;
	    momentum_err = momentum * match_function.GetBmErr(momentum);
	    hist_npl->Fill(ntbm->GetBabyMindMaximumPlane(muon_track_id));
	    hist_npl_range->Fill(ntbm->GetBabyMindMaximumPlane(muon_track_id),
				 momentum);
	  }
	  else {
	    chain.stop_flag = 0;
	  }
	}

	chain.bm_range_mom = momentum;
	chain.bm_range_mom_error[0] = momentum_err;
	chain.bm_range_mom_error[1] = momentum_err;
	hist_range->Fill(momentum);
	if ( chain.stop_flag == 1 ) 
	  std::cout << "Range : " << chain.bm_range_mom << ", "
		    << "MCS : " << chain.ecc_mcs_mom[0] << std::endl;
      }
      else { // partner

	std::vector<double > tangent(3);
	if ( chain.direction == 1 ) {
	  tangent.at(0) = chain.base.back().ax;
	  tangent.at(1) = chain.base.back().ay;
	}
	else if ( chain.direction == -1 ) {
	  tangent.at(0) = chain.base.front().ax;
	  tangent.at(1) = chain.base.front().ay;
	}
	tangent.at(2) = 1.;

	std::vector<double > ax, ay;
	std::vector<int > pl;
	ax.reserve(chain.base.size());
	ay.reserve(chain.base.size());
	pl.reserve(chain.base.size());
	if ( chain.direction == 1 ) {
	  for ( auto pair : chain.base_pair ) {
	    double x_diff = pair.first.x - pair.second.x;
	    double y_diff = pair.first.y - pair.second.y;
	    double z_diff = pair.first.z - pair.second.z;
	    ax.push_back(x_diff / z_diff);
	    ay.push_back(y_diff / z_diff);
	    pl.push_back(pair.first.pl);
	  }
	  ax.push_back(chain.base.back().ax);
	  ay.push_back(chain.base.back().ay);
	  pl.push_back(chain.base.back().pl);
	}
	else if ( chain.direction == -1 ) {
	  ax.push_back(chain.base.front().ax);
	  ay.push_back(chain.base.front().ay);
	  pl.push_back(chain.base.front().pl);
	  for ( auto pair : chain.base_pair ) {
	    double x_diff = pair.first.x - pair.second.x;
	    double y_diff = pair.first.y - pair.second.y;
	    double z_diff = pair.first.z - pair.second.z;
	    ax.push_back(x_diff / z_diff);
	    ay.push_back(y_diff / z_diff);
	    pl.push_back(pair.second.pl);
	  }
	}

	range_function.ModifyVectors(ax, ay, pl);
	double pion_range_ene = range_function.CalculateEnergyFromRange(ax, ay, pl, 211, chain.direction);
	pion_range_ene += MCS_PION_MASS;
	chain.ecc_range_mom[0] = CalculateMomentumFromEnergy(pion_range_ene, MCS_PION_MASS);

	double proton_range_ene = range_function.CalculateEnergyFromRange(ax, ay, pl, 2212, chain.direction);
	proton_range_ene += MCS_PROTON_MASS;
	chain.ecc_range_mom[1] = CalculateMomentumFromEnergy(proton_range_ene, MCS_PROTON_MASS);
	chain.ecc_range_mom_error[1][0] = range_function.CalculateProtonRangeError(chain.ecc_range_mom[1],
										   std::hypot(tangent.at(0), tangent.at(1)));
	chain.ecc_range_mom_error[1][1] = chain.ecc_range_mom_error[1][0];

      }
      
    }

  }

  Momentum_recon::WriteEventInformationBin((std::string)argv[4], ev_vec);
  
  TFile *ofile = new TFile(argv[5], "recreate");
  hist_npl->Write();
  hist_range->Write();
  hist_npl_range->Write();
  ofile->Close();

  std::exit(0);

}


