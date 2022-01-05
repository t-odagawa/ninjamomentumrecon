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

// system includes
#include <vector>
#include <array>

// B2 includes
#include <B2Reader.hh>
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

  if ( argc != 7 ) {
    BOOST_LOG_TRIVIAL(error) << "Usage : " << argv[0]
			     << " <input particle gun MC file name> <output root file name> <ncell> <initial pbeta bias> <true(0)/smear(1)> <MC(0)/Data(1)>";
    std::exit(1);
  }
    
  try {

    // Input file
    B2Reader reader(argv[1]);

    // Output file
    TString ofilename = argv[2];
    TFile *ofile = new TFile(ofilename, "recreate");
    TTree *otree = new TTree("tree", "tree");

    Int_t entry_in_daily_file;
    Int_t npl;
    std::vector<Double_t> ax, ay;
    std::vector<Double_t> vph, pixel_count;
    Double_t true_pbeta;
    Double_t initial_pbeta;
    Double_t recon_pbeta;
    Double_t recon_pbeta_err;
    Double_t log_likelihood[3][2] = {}; // 3 particles (mu, pi, p), 2 directions (up, down)
    Double_t recon_pbeta_candidate[3][2] = {};
    Double_t recon_pbeta_err_candidate[3][2] = {};
    Int_t fit_status[3][2] = {};
    // any other variables?
    otree->Branch("entry_in_daily_file", &entry_in_daily_file, "entry_in_daily_file/I");
    otree->Branch("npl", &npl, "npl/I");
    otree->Branch("ax", &ax);
    otree->Branch("ay", &ay);
    otree->Branch("vph", &vph);
    otree->Branch("pixel_count", &pixel_count);
    otree->Branch("true_pbeta",  &true_pbeta,  "true_pbeta/D");
    otree->Branch("initial_pbeta", &initial_pbeta, "initial_pbeta/D");
    otree->Branch("recon_pbeta", &recon_pbeta, "recon_pbeta/D");
    otree->Branch("recon_pbeta_err", &recon_pbeta_err, "recon_pbeta_err/D");
    otree->Branch("log_likelihood", log_likelihood, "log_likelihood[3][2]/D");
    otree->Branch("recon_pbeta_candidate", recon_pbeta_candidate, "recon_pbeta_candidate[3][2]/D");
    otree->Branch("recon_pbeta_err_candidate", recon_pbeta_err_candidate,
		  "recon_pbeta_err_candidate[3][2]/D");
    otree->Branch("fit_status", fit_status, "fit_status[3][2]/I");
    
    Bool_t smear_flag = atoi(argv[5]);
    if ( smear_flag )
      BOOST_LOG_TRIVIAL(info) << "========== Smear flag : TRUE ==========";
    else 
      BOOST_LOG_TRIVIAL(info) << "========== Smear flag : FALSE ==========";

    Int_t datatype = atoi(argv[6]);

    // Start pbeta reconstruction process
    while ( reader.ReadNextSpill() > 0 ) {
      // if ( reader.GetEntryNumber() > 10 ) break;
      auto &spill_summary = reader.GetSpillSummary();

      // Get emulsion tracks
      std::vector<const B2EmulsionSummary*> emulsions;
      auto it_emulsion = spill_summary.BeginEmulsion();
      while ( const auto *emulsion = it_emulsion.Next() ) {
	if ( datatype == B2DataType::kMonteCarlo &&
	    emulsion->GetParentTrackId() == 0 ) continue; // track id should be assigned when MC
	if ( datatype == B2DataType::kMonteCarlo &&
	     emulsion->GetParentTrack().GetParticlePdg() != PDG_t::kMuonMinus ) continue; // only muon when MC
	if ( emulsion->GetFilmType() != B2EmulsionType::kECC ) continue; // only ECC films
	if ( emulsion->GetEcc() != 4 ) continue; // only ECC5
	emulsions.push_back(emulsion);
      }

      if ( emulsions.empty() ) continue;

      entry_in_daily_file = reader.GetEntryNumber();
      npl = emulsions.size();

      std::sort(emulsions.begin(), emulsions.end(), EmulsionCompare);

      for ( const auto emulsion : emulsions ) {
	ax.push_back(emulsion->GetTangent().GetValue().X());
	ay.push_back(emulsion->GetTangent().GetValue().Y());
	vph.push_back(emulsion->GetVphUp() + emulsion->GetVphDown());
	pixel_count.push_back(emulsion->GetPixelCountUp() + emulsion->GetPixelCountDown());
      }

      // Get true information
      if ( datatype == B2DataType::kMonteCarlo ) {
	Int_t track_id_tmp_ = emulsions.at(0)->GetParentTrackId();
	Int_t particle_pdg_tmp_ = emulsions.at(0)->GetParentTrack().GetParticlePdg();
	Int_t ecc_tmp_ = emulsions.at(0)->GetEcc();
	if ( B2Pdg::IsMuonPlusOrMinus(particle_pdg_tmp_) )
	  true_pbeta = CalculatePBetaFromMomentum(emulsions.at(0)->GetMomentum().GetValue().Mag(), MCS_MUON_MASS);
	else if ( particle_pdg_tmp_ == PDG_t::kProton )
	  true_pbeta = CalculatePBetaFromMomentum(emulsions.at(0)->GetMomentum().GetValue().Mag(), MCS_PROTON_MASS);
	else if ( B2Pdg::IsChargedPion(particle_pdg_tmp_) )
	  true_pbeta = CalculatePBetaFromMomentum(emulsions.at(0)->GetMomentum().GetValue().Mag(), MCS_PION_MASS);
      }

      // input parameters for TMinuit
      const Int_t ncell = std::atoi(argv[3]); // Ncell will be sweeped in the next step
      std::vector<Double_t > basetrack_distance = {};
      std::vector<Double_t > water_basetrack_distance = {};
      std::vector<Double_t > track_tangent = {};
      std::vector<Int_t > plate_id = {};
      std::vector<Double_t > radial_angle_difference = {};
      std::vector<Double_t > lateral_angle_difference = {};

      for ( Int_t iemulsion_up = 0; iemulsion_up < emulsions.size(); iemulsion_up++ ) {

	const auto emulsion_up = emulsions.at(iemulsion_up);
	// Only use films across one iron plate
	 if ( emulsion_up->GetPlate() <= 3 ) continue;
	 if ( emulsion_up->GetPlate()%2 == 1 &&
	      emulsion_up->GetPlate() >= 15) continue; 

	TVector3 tangent_up = emulsion_up->GetTangentInDownCoordinate().GetValue();
	tangent_up = (1. / tangent_up.Z()) * tangent_up;
	TVector3 position_up = emulsion_up->GetFilmPositionInDownCoordinate().GetValue();
	for ( Int_t iemulsion_down = iemulsion_up + 1; iemulsion_down < emulsions.size(); iemulsion_down++ ) {
	  const auto emulsion_down = emulsions.at(iemulsion_down);
	  
	  Int_t plate_difference = emulsion_up->GetPlate() - emulsion_down->GetPlate();
	  
	  if ( plate_difference == 2 * ncell - 1 ) {
	    TVector3 tangent_down = emulsion_down->GetTangent().GetValue();
	    tangent_down = (1. / tangent_down.Z()) * tangent_down;
	    
	    TVector3 position_down = emulsion_down->GetFilmPosition().GetValue();
	    TVector3 displacement = position_down - position_up;
	    if ( datatype == B2DataType::kMonteCarlo && smear_flag ) { // smear
	      displacement = SmearDistanceVector(displacement, kNinjaIron);
	    }
	    switch (datatype) {
	    case B2DataType::kMonteCarlo :
	      basetrack_distance.push_back(displacement.Mag() / 850.e-3); // Nominal distance = 500 + 350 um
	      break;
	    case B2DataType::kRealData :
	      basetrack_distance.push_back(displacement.Mag() / 850.e-3);
	      break;
	    }
	    // Water distance
	    if ( emulsion_up->GetPlate() >= 18 && iemulsion_down < emulsions.size() - 1 ) {
	      const auto emulsion_2pl_down = emulsions.at(iemulsion_down + 1);
	      TVector3 position_down_for_water = emulsion_down->GetFilmPositionInDownCoordinate().GetValue();
	      TVector3 position_2pl_down = emulsion_2pl_down->GetFilmPosition().GetValue();
	      TVector3 water_displacement = position_2pl_down - position_down_for_water;
	      if ( datatype == B2DataType::kMonteCarlo && smear_flag ) { // smear
		water_displacement = SmearDistanceVector(water_displacement, kNinjaWater);
	      }
	      switch (datatype) {
	      case B2DataType::kMonteCarlo :
		water_basetrack_distance.push_back(water_displacement.Mag() / 2.868); // Nominal distance = 2300 + 109 * 2 + 350 um
		break;
	      case B2DataType::kRealData :
		water_basetrack_distance.push_back(water_displacement.Mag() / 2.868);
		break;
	      }
	    } else {
	      water_basetrack_distance.push_back(0.);
	    }

	    track_tangent.push_back(tangent_up.Mag());
	    plate_id.push_back(emulsion_up->GetPlate());

	    Double_t angle_difference_radial_new, angle_difference_lateral_new;
	    Double_t radial_angle_precision_new = 0.;
	    Double_t lateral_angle_precision_new = 0.;
	    Double_t tangent_accuracy_radial_new = 0.;
	    Double_t tangent_accuracy_lateral_new = 0.;
	    angle_difference_radial_new = RadialAngleDiffNew(tangent_up, tangent_down);
	    angle_difference_lateral_new = LateralAngleDiffNew(tangent_up, tangent_down);
	    if ( datatype == B2DataType::kMonteCarlo && smear_flag ) { // smear
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
      }
      
      // Get initial pbeta value
      
      TVector3 vertex_tangent = emulsions.at(0)->GetTangent().GetValue();
      if ( datatype == B2DataType::kMonteCarlo && smear_flag ) { // smear 
	vertex_tangent = SmearTangentVector(vertex_tangent);
      }

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
      if ( datatype == B2DataType::kMonteCarlo && !smear_flag )
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
	}

	if ( particle_id == 0 ) {
	  if ( log_likelihood[0][0] <= log_likelihood[0][1] ) {
	    recon_pbeta = result_array.at(0).at(0);
	    recon_pbeta_err = result_array.at(0).at(1);
	  } else {
	    recon_pbeta = result_array.at(1).at(0);
	    recon_pbeta_err = result_array.at(1).at(1);
	  }
	}
	
      }

      if ( recon_pbeta > 4999 ){
	//BOOST_LOG_TRIVIAL(info) << reader.GetEntryNumber();
	//std::exit(1);
      }
      otree->Fill();
      ax.clear(); ay.clear();
      vph.clear(); pixel_count.clear();

    }
    
    ofile->cd();
    otree->Write();
    ofile->Close();
    
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
