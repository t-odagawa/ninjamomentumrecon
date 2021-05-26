// boost includes
#include <boost/log/core.hpp>
#include <boost/log/trivial.hpp>
#include <boost/log/expressions.hpp>

// ROOT includes
#include <TFile.h>
#include <TMinuit.h>
#include <TMath.h>

// system includes
#include <vector>
#include <array>

// B2 includes
#include <B2Reader.hh>

#include "/home/t2k/odagawa/NinjaMomentumRecon/src/McsCommon.cpp"
#include "/home/t2k/odagawa/NinjaMomentumRecon/src/McsAngleMethod.cpp"
#include "MinimizerTest.hpp"

namespace logging = boost::log;

// Function to be minimized
// par[0]  : reconstructed pbeta obtained from minimum log likelihood
// par[1]  : Ncell (number of skipped film between the pair)
// par[2]  : Initial angle (true information or smeared if necessary) dz
// par[3]  : direction (1 or -1)
// par[4]  : number of pairs of the basetracks
// par[5-] : angle differences between basetracks
void LogLikelihood(Int_t &npar, Double_t *gin, Double_t &f, Double_t *par, Int_t iflag) {
  // Log likelihood
  Double_t ll = 0.; // constant (= n * log(2pi) / 2) is ignored

  // Momentum
  Double_t pbeta = par[0];
  
  for (Int_t ifilm = 0; ifilm < (UInt_t)par[4]; ifilm++) {

    Double_t sigma = SigmaAtIfilm(pbeta, (UInt_t)par[1], par[2]);

    if (TMath::Abs(par[5 + ifilm]) > 2.5 * sigma) continue;

    ll += 2 * TMath::Log(sigma) + par[5 + ifilm] * par[5 + ifilm] / sigma / sigma;

    Double_t energy = CalculateEnergyFromPBeta(pbeta, muon_mass);
    Double_t beta = CalculateBetaFromPBeta(pbeta, muon_mass);
    energy = energy - par[3] * par[2] * (BetheBlochIron(beta) + BetheBlochWater(beta));
    pbeta = CalculatePBetaFromEnergy(energy, muon_mass);
  }
  
  f = ll;
  
}

Double_t CalculateBetaFromPBeta(Double_t pbeta, Double_t mass) {
  Double_t energy = CalculateEnergyFromPBeta(pbeta, mass);
  Double_t momentum = CalculateMomentumFromPBeta(pbeta, mass);
  return momentum / energy;
}

Double_t CalculateEnergyFromMomentum(Double_t momentum, Double_t mass) {
  return TMath::Hypot(momentum, mass);
}

Double_t CalculateEnergyFromPBeta(Double_t pbeta, Double_t mass) {
  return 0.5 * (pbeta + TMath::Hypot(pbeta, 2. * mass));
}

Double_t CalculatePBetaFromEnergy(Double_t energy, Double_t mass) {
  Double_t momentum = CalculateMomentumFromEnergy(energy, mass);
  return momentum * momentum / energy;
}

Double_t CalculatePBetaFromMomentum(Double_t momentum, Double_t mass) {
  Double_t energy = CalculateEnergyFromMomentum(momentum, mass);
  return CalculatePBetaFromEnergy(energy, mass);
}

Double_t CalculateMomentumFromPBeta(Double_t pbeta, Double_t mass) {
  Double_t energy = CalculateEnergyFromPBeta(pbeta, mass);
  return CalculateMomentumFromEnergy(energy, mass);
}

Double_t CalculateMomentumFromEnergy(Double_t energy, Double_t mass) {
  return TMath::Sqrt(energy * energy - mass * mass);
}

Double_t BetheBlochIron(Double_t beta) {
  std::array<Double_t, 3> funcpar = CalculateBetheBlochIronParameters();
  return funcpar.at(0) / beta / beta * (TMath::Log(beta / TMath::Sqrt(1 - beta * beta)) + funcpar.at(1)) + funcpar.at(2);
}

std::array<Double_t, 3> CalculateBetheBlochIronParameters() {
  std::array<Double_t, 3> return_array = {0.0170, 26.96, 0.0890};
  return return_array;
}

Double_t BetheBlochWater(Double_t beta) {
  std::array<Double_t, 3> funcpar = CalculateBetheBlochWaterParameters();
  return funcpar.at(0) / beta / beta * (TMath::Log(beta / TMath::Sqrt(1 - beta * beta)) + funcpar.at(1)) + funcpar.at(2);
}

std::array<Double_t, 3> CalculateBetheBlochWaterParameters() {
  std::array<Double_t, 3> return_array = {0.0119, 35.39, 0.0724};
  return return_array;
}

Double_t SigmaAtIfilm(Double_t pbeta, UInt_t ncell, Double_t dz) {
  Double_t beta = 1.;
  Double_t radiation_length = calculate_radiation_length(ncell, dz);
  return 13.6 / pbeta * TMath::Sqrt(radiation_length) * (1. + 0.038 * TMath::Log(radiation_length / beta));
}
  

std::array<Double_t, 2> ReconstructPBeta(UInt_t ncell, Double_t dz, std::vector<Double_t> angle_difference) {


  const UInt_t num_of_param = 5 + angle_difference.size();

  TMinuit *min = new TMinuit(num_of_param); // TMinuit(n), n = number of parameters
  min->SetPrintLevel(-1);
  min->SetFCN(LogLikelihood); // minimize LogLikelihood function
  int ierflg = 0; // output level (0 = normal)
  

  // Parameter Names
  TString parname[num_of_param];
  parname[0] = "Reconstructed momentum";
  parname[1] = "Ncell";
  parname[2] = "Initial angle (dz)";
  parname[3] = "Track direction";
  parname[4] = "Number of angle differences";
  for (Int_t ipar = 5; ipar < num_of_param; ipar++) {
    parname[ipar] = Form("Angle difference %d", ipar - 5);
  }

  // Initial values
  Double_t vstart[num_of_param];
  vstart[0] = 1000;
  vstart[1] = ncell;
  vstart[2] = dz;
  vstart[3] = 1;
  vstart[4] = angle_difference.size();
  for (Int_t ipar = 5; ipar < num_of_param; ipar++) {
    vstart[ipar] = angle_difference.at(ipar - 5);
  }
  
  // Step
  Double_t step[num_of_param];
  step[0] = 1.;
  for (Int_t ipar = 1; ipar < num_of_param; ipar++) {
    step[ipar] = 0.;
  }

  // Parameter setting
  // Bethe-Bloch is applicable to pbeta > ~20 (beta > 0.4)
  min->mnparm(0, parname[0], vstart[0], step[0], 20, 2000, ierflg);
  for (Int_t ipar = 1; ipar < num_of_param; ipar++) {
    min->mnparm(ipar, parname[ipar], vstart[ipar], step[ipar], 0, 0, ierflg);
  }
  
  
  // Fix parameters
  for (Int_t ipar = 1; ipar < num_of_param; ipar++) {
    min->FixParameter(ipar);
  }
  
  // array to be used for several steps
  Double_t arglist[2] = {};
  
  // Set delta log-likelihood corresponds to 1 sigma
  arglist[0] = 1.;
  min->mnexcm("SET ERR", arglist, 1, ierflg);
  
  // Execute minimizer (migrad method)
  arglist[0] = 1000; // maximum number of calls
  arglist[1] = 1; // tolerance
  min->mnexcm("MIGRAD", arglist, 2, ierflg);
  
  Double_t rec_mom, rec_mom_err;
  min->GetParameter(0, rec_mom, rec_mom_err);
  std::array<Double_t,2> return_array = {rec_mom, rec_mom_err};

  BOOST_LOG_TRIVIAL(info) << "Reconstructed Momentum = "
			   << rec_mom << " +/- "
			   << rec_mom_err << " [MeV/c]";
  
  delete min;

  return return_array;
 
}


int main (int argc, char *argv[]) {

  logging::core::get()->set_filter
    (
     logging::trivial::severity >= logging::trivial::info
     //logging::trivial::severity >= logging::trivial::debug
     );
  
  BOOST_LOG_TRIVIAL(info) << "==========Momentum Reconstruction Start==========";

  if (argc != 4) {
    BOOST_LOG_TRIVIAL(error) << "Usage : " << argv[0]
			     << " <input particle gun MC file name> <output root file name> <ncell>";
    std::exit(1);
  }
    
  try {

    B2Reader reader(argv[1]);

    TString ofilename = argv[2];
    TFile *ofile = new TFile(ofilename, "recreate");
    TTree *otree = new TTree("tree", "tree");

    Double_t true_pbeta, recon_pbeta;
    Double_t recon_pbeta_err;
    //Double_t log_likelihood; // any other variables?
    otree->Branch("true_pbeta",  &true_pbeta,  "true_pbeta/D");
    otree->Branch("recon_pbeta", &recon_pbeta, "recon_pbeta/D");
    otree->Branch("recon_pbeta_err", &recon_pbeta_err, "recon_pbeta_err/D");
    
    while (reader.ReadNextSpill() > 0) {
      auto &spill_summary = reader.GetSpillSummary();

      // Get emulsion tracks
      std::vector<const B2EmulsionSummary*> emulsions;
      auto it_emulsion = spill_summary.BeginEmulsion();
      while (const auto *emulsion = it_emulsion.Next()) {
	if (emulsion->GetParentTrackId() == 0) continue; // track id should be assigned
	if (emulsion->GetParentTrack().GetParticlePdg() != 13) continue; // only muon
	if (emulsion->GetFilmType() != B2EmulsionType::kECC) continue; // only ECC films
	if (emulsion->GetEcc() != 4) continue; // only ECC5
	if (emulsion->GetPlate () < 15) continue; // no Iron ECC films
	emulsions.push_back(emulsion);
      }

      if (emulsions.size() <= 0) continue;

      std::sort(emulsions.begin(), emulsions.end(), emulsion_compare);

      Int_t track_id_tmp_ = emulsions.at(0)->GetParentTrackId();
      Int_t ecc_tmp_ = emulsions.at(0)->GetEcc();
      true_pbeta = emulsions.at(0)->GetMomentum().GetValue().Mag();
      true_pbeta = true_pbeta * true_pbeta / TMath::Sqrt(true_pbeta * true_pbeta + muon_mass * muon_mass);
      TVector3 vertex_tangent = emulsions.at(0)->GetTangent().GetValue();
      Double_t momentum_difference_max = 0.;

      // input parameters for TMinuit
      const Int_t ncell = std::atoi(argv[3]); // Ncell will be sweeped in the next step
      Double_t dz = vertex_tangent.Mag(); // dz is calculated from the initial tangent
      std::vector<Double_t> angle_difference = {};

      for (Int_t iemulsion_up = 0; iemulsion_up < emulsions.size(); iemulsion_up++) {
	// lateral angle differences are pushed back to the vector
	// always use the upstream film in a tracking unit as a start film
	const auto emulsion_up = emulsions.at(iemulsion_up);
	if (emulsion_up->GetPlate()%2 == 1) continue; 
	TVector3 tangent_up = emulsion_up->GetTangent().GetValue();
	
	for (Int_t iemulsion_down = iemulsion_up + 1; iemulsion_down < emulsions.size(); iemulsion_down++) {
	  const auto emulsion_down = emulsions.at(iemulsion_down);
	  Int_t plate_difference = emulsion_up->GetPlate() - emulsion_down->GetPlate();
	  Int_t momentum_difference = emulsion_up->GetMomentum().GetValue().Mag()
	    - emulsion_down->GetMomentum().GetValue().Mag();
	  if (momentum_difference > momentum_difference_max) momentum_difference_max = momentum_difference;

	  if (plate_difference == 2 * ncell - 1) {
	    TVector3 tangent_down = emulsion_down->GetTangent().GetValue();
	    angle_difference.push_back(TMath::ATan(get_angle_difference_lateral(tangent_up,
										tangent_down,
										vertex_tangent)));
	    
	    BOOST_LOG_TRIVIAL(debug) << "Angle difference between"
				     << " film " << emulsion_up->GetPlate() << " and"
				     << " film " << emulsion_down->GetPlate()
				     << " is " << angle_difference.back();
	    break;
	  }
	}
      }


      std::array<Double_t,2> result_array = ReconstructPBeta(ncell, dz, angle_difference);
      recon_pbeta = result_array.at(0);
      recon_pbeta_err = result_array.at(1);
      otree->Fill();
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

// These functions are not used anymore

Double_t MomentumDrop(Double_t momentum, UInt_t ncell, Double_t dz) {
  std::array<Double_t, 3> funcpar = CalculateFunctionParameters(dz);
  return funcpar.at(0) / (momentum - funcpar.at(1)) + funcpar.at(2);
}

std::array<Double_t, 3> CalculateFunctionParameters(Double_t dz) {
  std::array<Double_t, 3> return_array = {4.19, 256.1, 0.626 * dz};
  return return_array;
}

Double_t CalculateRadLenUnit(UInt_t ncell) {
  if (ncell > 10)
    throw std::out_of_range("Ncell should be less than 10");

  Double_t rad_len_unit = 0;

  for (Int_t imaterial = 0; imaterial < kNumberOfNinjaMaterials; imaterial++) {
    int num_of_layers = 0;
    switch (imaterial) {
    case kNinjaIron : 
      num_of_layers = ncell;
      break;
    case kNinjaWater :
      num_of_layers = ncell - 1;
      break;
    case kNinjaGel :
      num_of_layers = 4 * ncell - 2;
      break;
    case kNinjaBase :
    case kNinjaPacking :
      num_of_layers = 2 * ncell - 2;
      break;
    }
    rad_len_unit += num_of_layers * MATERIAL_THICK[imaterial] / RAD_LENGTH[imaterial];
  }

  return rad_len_unit;

}
