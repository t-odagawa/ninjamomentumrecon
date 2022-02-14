// boost includes
#include <boost/log/core.hpp>
#include <boost/log/trivial.hpp>
#include <boost/log/expressions.hpp>

// ROOT includes
#include <TFile.h>
#include <TTree.h>
#include <TString.h>
#include <TVector3.h>
#include <TMath.h>
#include <TRandom3.h>
#include <TGraph.h>
#include <TCanvas.h>
#include <TAxis.h>
#include <TLine.h>
#include <TSpline.h>

// system includes
#include <vector>
#include <array>

// B2 includes
#include <B2Reader.hh>
#include <B2SpillSummary.hh>
#include <B2EventSummary.hh>
#include <B2VertexSummary.hh>
#include <B2TrackSummary.hh>
#include <B2EmulsionSummary.hh>
#include <B2Pdg.hh>

// MCS includes
#include "McsConst.hpp"
#include "McsFunction.hpp"

namespace logging = boost::log;

int main (int argc, char *argv[]) {

  logging::core::get()->set_filter
    (
     logging::trivial::severity >= logging::trivial::info
     //logging::trivial::severity >= logging::trivial::trace
     //logging::trivial::severity >= logging::trivial::debug
     );
  
  BOOST_LOG_TRIVIAL(info) << "==========Momentum Reconstruction Start==========";

  if ( argc != 8 ) {
    BOOST_LOG_TRIVIAL(error) << "Usage : " << argv[0]
			     << " <input particle gun MC file name> <output root file name> <ncell> <initial pbeta bias> <true(0)/smear(1)> <MC(0)/Data(1)> <Draw LL option>";
    std::exit(1);
  }
    
  try {

    gErrorIgnoreLevel = kWarning;

    // Input file
    B2Reader reader(argv[1]);

    // Output file
    TString ofilename = argv[2];
    TFile *ofile = new TFile(ofilename, "recreate");
    TTree *otree = new TTree("tree", "tree");

    const Int_t ncell = std::atoi(argv[3]); // Ncell will be sweeped in the next step

    Int_t entry_in_daily_file;
    Double_t weight = 1.;
    Double_t true_pbeta = -1;
    Int_t true_particle_id = -1;
    Int_t muon_track_id;
    Int_t true_direction;
    Int_t npl;
    std::vector<Int_t> pl;
    std::vector<Double_t> ax, ay;
    std::vector<Double_t> vph, pixel_count;
    Double_t initial_pbeta;
    Double_t recon_pbeta;
    Double_t recon_pbeta_err;
    Double_t log_likelihood[3][2] = {}; // 3 particles (mu, pi, p), 2 directions (up, down)
    Double_t recon_pbeta_candidate[3][2] = {};
    Double_t recon_pbeta_err_candidate[3][2] = {};
    Int_t fit_status[3][2] = {};
    // any other variables?
    otree->Branch("entry_in_daily_file", &entry_in_daily_file, "entry_in_daily_file/I");
    otree->Branch("weight", &weight, "weight/D");
    otree->Branch("true_pbeta",  &true_pbeta,  "true_pbeta/D");
    otree->Branch("true_particle_id", &true_particle_id, "true_particle_id/I");
    otree->Branch("muon_track_id", &muon_track_id, "muon_track_id/I");
    otree->Branch("true_direction", &true_direction, "true_direction/I");
    otree->Branch("npl", &npl, "npl/I");
    otree->Branch("pl", &pl);
    otree->Branch("ax", &ax);
    otree->Branch("ay", &ay);
    otree->Branch("vph", &vph);
    otree->Branch("pixel_count", &pixel_count);
    otree->Branch("initial_pbeta", &initial_pbeta, "initial_pbeta/D");
    otree->Branch("recon_pbeta", &recon_pbeta, "recon_pbeta/D");
    otree->Branch("recon_pbeta_err", &recon_pbeta_err, "recon_pbeta_err/D");
    otree->Branch("log_likelihood", log_likelihood, "log_likelihood[3][2]/D");
    otree->Branch("recon_pbeta_candidate", recon_pbeta_candidate, "recon_pbeta_candidate[3][2]/D");
    otree->Branch("recon_pbeta_err_candidate", recon_pbeta_err_candidate,
		  "recon_pbeta_err_candidate[3][2]/D");
    otree->Branch("fit_status", fit_status, "fit_status[3][2]/I");
    
    Bool_t smear_flag = atoi(argv[5]);
    Int_t datatype = atoi(argv[6]);

    if ( smear_flag ) {
      if ( datatype == B2DataType::kMonteCarlo ) 
	BOOST_LOG_TRIVIAL(info) << "========== Smear flag : TRUE ==========";
      else if ( datatype == B2DataType::kRealData ) {
	BOOST_LOG_TRIVIAL(info) << "========== Smear flag : FALSE (Real data are never smeared)==========";
	smear_flag = false;
      }
    }
    else 
      BOOST_LOG_TRIVIAL(info) << "========== Smear flag : FALSE ==========";

    Bool_t scan_flag = atoi(argv[7]);

    TCanvas *c = new TCanvas("c", "c");
    TLine *l = new TLine();
    l->SetLineColor(kRed);
    std::vector<Double_t> x, y;
    std::vector<Double_t> x_plus, y_plus;
    std::vector<Double_t> x_minus, y_minus;
    
    if ( scan_flag ) {
      BOOST_LOG_TRIVIAL(info) << "========== Scan flag : TRUE ==========";
      c->Print(ofilename + ".pdf" + "[", "pdf");
    }
    else {
      BOOST_LOG_TRIVIAL(info) << "========== Scan flag : FALSE ==========";
    }

    const Int_t ecc_id = 4;

    // Start pbeta reconstruction process
    while ( reader.ReadNextSpill() > 0 ) {
      // if ( reader.GetEntryNumber() > 10 ) break;
      auto &spill_summary = reader.GetSpillSummary();

      if ( datatype == B2DataType::kMonteCarlo ) {
	auto it_event = spill_summary.BeginTrueEvent();
	const auto *event = it_event.Next();
	Double_t norm = event->GetNormalization();
	
	auto &primary_vertex_summary = event->GetPrimaryVertex();
	Double_t total_cross_section = primary_vertex_summary.GetTotalCrossSection();
	weight = norm * total_cross_section * 1e-38 * 6.02e23;
      }

      // Get emulsion tracks
      std::vector<const B2EmulsionSummary*> emulsions;
      auto it_emulsion = spill_summary.BeginEmulsion();
      while ( const auto *emulsion = it_emulsion.Next() ) {
	if ( datatype == B2DataType::kMonteCarlo &&
	    emulsion->GetParentTrackId() == 0 ) continue; // track id should be assigned when MC
	if ( emulsion->GetFilmType() != B2EmulsionType::kECC ) continue; // only ECC films
	if ( emulsion->GetEcc() != ecc_id ) continue;
	emulsions.push_back(emulsion);
      }

      if ( emulsions.empty() ) continue;

      std::sort(emulsions.begin(), emulsions.end(), EmulsionCompare);
      entry_in_daily_file = reader.GetEntryNumber();

      // Separate emulsions into chains
      std::vector<std::vector<const B2EmulsionSummary* > > emulsion_single_chains;
      Int_t track_id_tmp_ = emulsions.at(0)->GetParentTrackId();
      std::vector<const B2EmulsionSummary* > emulsion_single_chain;
      for ( Int_t iemulsion = 0; iemulsion < emulsions.size(); iemulsion++ ) {
	if ( emulsions.at(iemulsion)->GetParentTrackId() == track_id_tmp_ ) {
	  emulsion_single_chain.push_back(emulsions.at(iemulsion));
	} else {
	  emulsion_single_chains.push_back(emulsion_single_chain);
	  emulsion_single_chain.clear();
	  track_id_tmp_ = emulsions.at(iemulsion)->GetParentTrackId();
	  emulsion_single_chain.push_back(emulsions.at(iemulsion));
	}
      }
      emulsion_single_chains.push_back(emulsion_single_chain);

      // loop for each chain
      for ( auto chain : emulsion_single_chains ) {
	npl = chain.size();
	if ( chain.at(0)->GetTangent().GetValue().Z() )
	  true_direction = 1;
	else true_direction = -1;

	// Get true information
	if ( datatype == B2DataType::kMonteCarlo ) {
	  Int_t particle_pdg_ = chain.at(0)->GetParentTrack().GetParticlePdg();
	  true_particle_id = particle_pdg_;
	  Int_t momentum_ = chain.at(0)->GetParentTrack().GetInitialAbsoluteMomentum().GetValue();
	  if ( B2Pdg::IsMuonPlusOrMinus(particle_pdg_) )
	    true_pbeta = CalculatePBetaFromMomentum(momentum_, MCS_MUON_MASS);
	  else if ( B2Pdg::IsChargedPion(particle_pdg_) )
	    true_pbeta = CalculatePBetaFromMomentum(momentum_, MCS_PION_MASS);	  
	  else if ( particle_pdg_ == PDG_t::kProton )
	    true_pbeta = CalculatePBetaFromMomentum(momentum_, MCS_PROTON_MASS);
	  else 
	    BOOST_LOG_TRIVIAL(warning) << "Particle is not in interest";
	}
	
	for ( const auto emulsion : chain ) {
	  pl.push_back(emulsion->GetPlate() + 1);
	  ax.push_back(emulsion->GetTangent().GetValue().X());
	  ay.push_back(emulsion->GetTangent().GetValue().Y());
	  vph.push_back(emulsion->GetVphUp() + emulsion->GetVphDown());
	  pixel_count.push_back(emulsion->GetPixelCountUp() + emulsion->GetPixelCountDown());
	}

	// input parameter vectors for TMinuit
	std::vector<Double_t > basetrack_distance = {};
	std::vector<Double_t > water_basetrack_distance = {};
	std::vector<Double_t > track_tangent = {};
	std::vector<Int_t > plate_id = {};
	std::vector<Double_t > radial_angle_difference = {};
	std::vector<Double_t > lateral_angle_difference = {};

	for ( Int_t iemulsion_up = 0; iemulsion_up < chain.size(); iemulsion_up++ ) {
	  const auto emulsion_up = chain.at(iemulsion_up);
	  // Only use films across one iron plate
	  if ( emulsion_up->GetPlate() <= 3 ) continue;
	  if ( emulsion_up->GetPlate()%2 == 1 &&
	       emulsion_up->GetPlate() >= 15 ) continue;

	  TVector3 position_up = emulsion_up->GetFilmPositionInDownCoordinate().GetValue();
	  TVector3 tangent_up = emulsion_up->GetTangentInDownCoordinate().GetValue();
	  tangent_up = (1. / tangent_up.Z()) * tangent_up;
	  for ( Int_t iemulsion_down = iemulsion_up + 1; iemulsion_down < chain.size(); iemulsion_down++ ) {
	    const auto emulsion_down = chain.at(iemulsion_down);
	    
	    Int_t plate_difference = emulsion_up->GetPlate() - emulsion_down->GetPlate();
	    if ( plate_difference != 2 * ncell - 1 ) continue;

	    TVector3 position_down = emulsion_down->GetFilmPosition().GetValue();
	    TVector3 tangent_down = emulsion_down->GetTangent().GetValue();
	    tangent_down = (1. / tangent_down.Z()) * tangent_down;
	    
	    TVector3 displacement = position_down - position_up;
	    if ( smear_flag ) SmearDistanceVector(displacement, kNinjaIron);
	    basetrack_distance.push_back(displacement.Mag() / 850.e-3); // 500 + 350 um
	    // vector push back should be done after the corresponding downstream film found
	    track_tangent.push_back(tangent_up.Mag());
	    plate_id.push_back(emulsion_up->GetPlate());

	    // Water distance calculation
	    if ( emulsion_up->GetPlate() >= 18 && iemulsion_down < chain.size() - 1 ) {
	      const auto emulsion_2pl_down = chain.at(iemulsion_down + 1);
	      TVector3 position_down_for_water = emulsion_down->GetFilmPositionInDownCoordinate().GetValue();
	      TVector3 position_2pl_down = emulsion_2pl_down->GetFilmPosition().GetValue();
	      TVector3 water_displacement = position_2pl_down - position_down_for_water;
	      if ( smear_flag ) SmearDistanceVector(water_displacement, kNinjaWater);
	      water_basetrack_distance.push_back(water_displacement.Mag() / 2.868);
	    } else {
	      water_basetrack_distance.push_back(0.);
	    }

	    Double_t angle_difference_radial_new, angle_difference_lateral_new;
	    Double_t radial_angle_precision_new = 0.;
	    Double_t lateral_angle_precision_new = 0.;
	    Double_t tangent_accuracy_radial_new = 0.;
	    Double_t tangent_accuracy_lateral_new = 0.;
	    angle_difference_radial_new = RadialAngleDiffNew(tangent_up, tangent_down);
	    angle_difference_lateral_new = LateralAngleDiffNew(tangent_up, tangent_down);
	    if ( smear_flag ) {
	      radial_angle_precision_new = RadialAnglePrecision(TMath::Hypot(tangent_up.X(), tangent_up.Y()));
	      lateral_angle_precision_new = LateralAnglePrecision(TMath::Hypot(tangent_up.X(), tangent_up.Y()));
	    }
	    radial_angle_difference.push_back(gRandom->Gaus(angle_difference_radial_new, radial_angle_precision_new));
	    lateral_angle_difference.push_back(gRandom->Gaus(angle_difference_lateral_new, lateral_angle_precision_new));
	    
	    BOOST_LOG_TRIVIAL(debug) << " Film " << emulsion_up->GetPlate() << " and"
				     << " Film " << emulsion_down->GetPlate()
				     << " Radial angle difference is " << radial_angle_difference.back() << ", "
				     << " Lateral angle difference is " << lateral_angle_difference.back() << ", "
				     << " Base track distance is " << basetrack_distance.back() << ", "
				     << " Water base track distance is " << water_basetrack_distance.back() << ","
				     << " Track tangnet is " << track_tangent.back();
	    break;

	  }

	}

	// Get initial pbeta value
	TVector3 vertex_tangent = chain.at(0)->GetTangent().GetValue();
	if ( smear_flag ) SmearTangentVector(vertex_tangent);

	Double_t radial_angle_difference_rms = 0.;
	Double_t lateral_angle_difference_rms = 0.;
	for ( Int_t ipair = 0; ipair < radial_angle_difference.size(); ipair++ ) {
	  radial_angle_difference_rms += radial_angle_difference.at(ipair) * radial_angle_difference.at(ipair);
	  lateral_angle_difference_rms += lateral_angle_difference.at(ipair) * lateral_angle_difference.at(ipair);
	}

	if ( smear_flag || datatype == B2DataType::kRealData ) {
	  radial_angle_difference_rms /= radial_angle_difference.size();
	  lateral_angle_difference_rms /= lateral_angle_difference.size();      
	} else {
	  radial_angle_difference_rms = (radial_angle_difference_rms + lateral_angle_difference_rms) 
	    / (radial_angle_difference.size() + lateral_angle_difference.size());
	}
	radial_angle_difference_rms = TMath::Sqrt(radial_angle_difference_rms);
	lateral_angle_difference_rms = TMath::Sqrt(lateral_angle_difference_rms);
	if ( !smear_flag )
	  lateral_angle_difference_rms = radial_angle_difference_rms; // angle difference are equivalent in true-level

	Double_t radial_cut_value, lateral_cut_value;
	radial_cut_value = 3. * radial_angle_difference_rms;
	lateral_cut_value = 3. * lateral_angle_difference_rms;
	// It is better if cut value is too small, accept all angle differences
	Double_t radiation_length = CalcRadLength(ncell, vertex_tangent.Mag());
	initial_pbeta = MCS_SCALE_FACTOR * 13.6 / lateral_angle_difference_rms * TMath::Sqrt(radiation_length) * (1. + 0.038 * TMath::Log(radiation_length));
	// only used in the performance study
	const Double_t initial_pbeta_bias = std::atof(argv[4]);
	initial_pbeta *= initial_pbeta_bias;

	Int_t direction;
	std::array<std::array<Double_t, 3>, 2> result_array = {};
	for ( Int_t particle_id = 0; particle_id < kNumberOfNinjaMcsParticles; particle_id++ ) {
	  for ( Int_t idirection = 0; idirection < kNumberOfNinjaMcsDirections; idirection++ ) {
	    
	    direction = MCS_DIRECTION[idirection];
	    
	    result_array.at(idirection) = ReconstructPBeta(initial_pbeta, ncell, particle_id, direction,
							   radial_cut_value, lateral_cut_value,
							   smear_flag || (datatype == B2DataType::kRealData),
							   basetrack_distance,
							   water_basetrack_distance,
							   track_tangent,
							   plate_id,
							   radial_angle_difference,
							   lateral_angle_difference);
	    recon_pbeta_candidate[particle_id][idirection] = result_array.at(idirection).at(0);
	    recon_pbeta_err_candidate[particle_id][idirection] = result_array.at(idirection).at(1);
	    fit_status[particle_id][idirection] = (Int_t)result_array.at(idirection).at(2);
	    log_likelihood[particle_id][idirection] = FuncNegativeLogLikelihood(result_array.at(idirection).at(0),
										ncell, particle_id, direction,
										radial_cut_value, lateral_cut_value,
										smear_flag || (datatype == B2DataType::kRealData),
										basetrack_distance,
										water_basetrack_distance,
										track_tangent,
										plate_id,
										radial_angle_difference,
										lateral_angle_difference);

	    if ( scan_flag && direction == 1 && particle_id == 0) {
	      int i = 1;
	      do {
		x_plus.push_back(result_array.at(0).at(0) + 1. * i);
		y_plus.push_back(FuncNegativeLogLikelihood(result_array.at(0).at(0) + 1. * i,
							   ncell, 0, 1,
							   radial_cut_value, lateral_cut_value,
							   smear_flag || (datatype == B2DataType::kRealData),
							   basetrack_distance,
							   water_basetrack_distance,
							   track_tangent,
							   plate_id,
							   radial_angle_difference,
							   lateral_angle_difference) - log_likelihood[0][0]);
		i++;
	      } while ( y_plus.back() <= 10. && i < 1000);
	      i = 1;
	      do {
		x_minus.push_back(result_array.at(0).at(0) - 1. * i);
		y_minus.push_back(FuncNegativeLogLikelihood(result_array.at(0).at(0) - 1. * i,
							    ncell, 0, 1,
							    radial_cut_value, lateral_cut_value,
							    smear_flag || (datatype == B2DataType::kRealData),
							    basetrack_distance,
							    water_basetrack_distance,
							    track_tangent,
							    plate_id,
							    radial_angle_difference,
							    lateral_angle_difference) - log_likelihood[0][0]);
		i++;
	      } while ( (y_minus.back() <= 10. || 
			 x_minus.back() > 0) &&
			i < 1000 );
	      std::reverse(x_minus.begin(), x_minus.end());
	      std::reverse(y_minus.begin(), y_minus.end());
	      x.insert(x.end(), x_minus.begin(), x_minus.end());
	      y.insert(y.end(), y_minus.begin(), y_minus.end());
	      x.push_back(result_array.at(0).at(0));
	      y.push_back(0.);
	      x.insert(x.end(), x_plus.begin(), x_plus.end());
	      y.insert(y.end(), y_plus.begin(), y_plus.end());
	      
	      TGraph *g = new TGraph(x.size(), &x[0], &y[0]);
	      g->SetMarkerStyle(4);
	      g->SetMarkerSize(1);
	      g->SetTitle(Form("Negative Log Likelihood Scan (Entry %lld);Reconstructed p#beta [MeV/c];-2ln(L) + 2ln(L_{max})",
			       reader.GetEntryNumber()));
	      g->GetYaxis()->SetRangeUser(-.5, 20.);
	      g->Draw("AP");
	      l->DrawLine(true_pbeta, -0.5, true_pbeta, 20.);
	      c->Print(ofilename + ".pdf", "pdf");
	      x.clear(); x_plus.clear(); x_minus.clear();
	      y.clear(); y_plus.clear(); y_minus.clear();
	    } // scan flag	    

	  } // direction loop
	  
	  if ( particle_id == 0 ) {
	    if ( log_likelihood[0][0] <= log_likelihood[0][1] ) {
	      recon_pbeta = result_array.at(0).at(0);
	      recon_pbeta_err = result_array.at(0).at(1);
	    } else {
	      recon_pbeta = result_array.at(1).at(0);
	      recon_pbeta_err = result_array.at(1).at(1);
	    }
	  }
	} // particle id loop
	
	// Fill 
	otree->Fill();
	ax.clear(); ay.clear();
	vph.clear(); pixel_count.clear();

      } // chain loop

    }
    
    ofile->cd();
    otree->Write();
    ofile->Close();

    if ( scan_flag )
      c->Print(ofilename + ".pdf" + "]", "pdf");
      
  } catch (const std::runtime_error &error) {
    BOOST_LOG_TRIVIAL(fatal) << "Runtime error : " << error.what();
    std::exit(1);
  } catch (const std::invalid_argument &error) {
    BOOST_LOG_TRIVIAL(fatal) << "Invalid argument : " << error.what();
    std::exit(1);
  }

  BOOST_LOG_TRIVIAL(info) << "==========Momentum Reconstruction Finish==========";
  std::exit(0);

}
