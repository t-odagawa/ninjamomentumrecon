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
// par[0]           : reconstructed pbeta obtained from minimum log likelihood
// par[1]           : Ncell (number of skipped film between the pair)
// par[2]           : particle id (muon, charged pion, proton)
// par[3]           : direction (1 or -1)
// par[4]           : radial angle difference cut value
// par[5]           : lateral angle difference cut value
// par[6]           : true(0)/smear(1) flag
// par[7]           : number of pairs of the basetracks ( = N )
// par[   8 -  N+7] : basetrack distance (used for energy deposition and radiation length calculation)
// par[ N+8 - 2N+7] : downstream water basetrack distance (used for energy deposition)
// par[2N+8 - 3N+7] : upstream track tangent
// par[3N+8 - 4N+7] : upstream plate id
// par[4N+8 - 5N+7] : radial angle differences between basetracks
// par[5N+8 - 6N+7] : lateral angle differences between basetracks
void LogLikelihood(Int_t &npar, Double_t *gin, Double_t &f, Double_t *par, Int_t iflag) {

  // Initialzation of parameters
  // Momentum
  Double_t pbeta = par[0];
  // Ncell
  UInt_t ncell = (UInt_t)par[1];
  // particle id
  Int_t particle_id = (Int_t)par[2];
  // direction
  Int_t direction = (Int_t)par[3];
  // radial angle difference cut value
  Double_t radial_cut_value = par[4];
  // lateral angle difference cut value
  Double_t lateral_cut_value = par[5];
  // Smear flag
  Bool_t smear_flag = (Bool_t) par[6];
  // number of pairs
  const UInt_t number_of_pairs = (UInt_t)par[7];
  // basetrack distance
  std::vector<Double_t> basetrack_distance = {};
  basetrack_distance.resize(number_of_pairs);
  // downstream water basetrack distance
  std::vector<Double_t> water_basetrack_distance = {};
  water_basetrack_distance.resize(number_of_pairs);
  // upstream track tangent
  std::vector<Double_t> track_tangent = {};
  track_tangent.resize(number_of_pairs);
  // upstream plate id
  std::vector<Int_t> plate_id = {};
  plate_id.resize(number_of_pairs);
  // radial angle differences
  std::vector<Double_t> radial_angle_difference = {};
  radial_angle_difference.resize(number_of_pairs);
  // lateral angle differences
  std::vector<Double_t> lateral_angle_difference = {};
  lateral_angle_difference.resize(number_of_pairs);
  for(Int_t ipairs = 0; ipairs < number_of_pairs; ipairs++) {
    basetrack_distance.at(ipairs)       =        par[8                       + ipairs];
    water_basetrack_distance.at(ipairs) =        par[8 + 1 * number_of_pairs + ipairs];
    track_tangent.at(ipairs)            =        par[8 + 2 * number_of_pairs + ipairs];
    plate_id.at(ipairs)                 = (Int_t)par[8 + 3 * number_of_pairs + ipairs];
    radial_angle_difference.at(ipairs)  =        par[8 + 4 * number_of_pairs + ipairs];
    lateral_angle_difference.at(ipairs) =        par[8 + 5 * number_of_pairs + ipairs];
  }
    
  f = FuncLogLikelihood(pbeta, ncell, particle_id, direction,
			radial_cut_value, lateral_cut_value,
			smear_flag,
			basetrack_distance,
			water_basetrack_distance,
			track_tangent,
			plate_id,
			radial_angle_difference,
			lateral_angle_difference);
  
}

Double_t FuncLogLikelihood(Double_t pbeta,
			   UInt_t ncell,
			   Int_t particle_id,
			   Int_t direction,
			   Double_t radial_cut_value,
			   Double_t lateral_cut_value,
			   Bool_t smear_flag,
			   std::vector<Double_t> basetrack_distance,
			   std::vector<Double_t> water_basetrack_distance,
			   std::vector<Double_t> track_tangent,
			   std::vector<Int_t> plate_id,
			   std::vector<Double_t> radial_angle_difference,
			   std::vector<Double_t> lateral_angle_difference) {

  // vector size check
  std::vector<std::size_t> par_vect_sizes = {};
  par_vect_sizes.push_back(basetrack_distance.size());
  par_vect_sizes.push_back(water_basetrack_distance.size());
  par_vect_sizes.push_back(track_tangent.size());
  par_vect_sizes.push_back(plate_id.size());
  par_vect_sizes.push_back(radial_angle_difference.size());
  par_vect_sizes.push_back(lateral_angle_difference.size());

  for (Int_t ivect = 0; ivect < par_vect_sizes.size() - 1; ivect++) {
    if (par_vect_sizes.at(ivect) != par_vect_sizes.at(ivect + 1)) {
    BOOST_LOG_TRIVIAL(error) << "Vector size different";
    BOOST_LOG_TRIVIAL(error) << "Basetrack distance size : " << basetrack_distance.size() << ", "
			     << "Water basetrack distance size : " << water_basetrack_distance.size() << ", "
			     << "Tangent size : " << track_tangent.size() << ", "
			     << "Plate size : " << plate_id.size() << ", "
			     << "Radial angle difference size : " << radial_angle_difference.size() << ", "
			     << "Lateral angle difference size : " << lateral_angle_difference.size();
    std::exit(1);
    }
  }

  const Int_t number_of_pairs = basetrack_distance.size();

  // Calculate negative log likelihood
  Double_t ll = 0; // constant (= n * log(2pi) / 2) is ignored

  for (Int_t ipairs = 0; ipairs < number_of_pairs; ipairs++) {

    Double_t radial_sigma, lateral_sigma;

    if (smear_flag){
      radial_sigma = RadialSigmaAtIfilm(pbeta, ncell, basetrack_distance.at(ipairs),
					track_tangent.at(ipairs));
      lateral_sigma = LateralSigmaAtIfilm(pbeta, ncell, basetrack_distance.at(ipairs),
					  track_tangent.at(ipairs));
    } else {
      radial_sigma = HighlandSigmaAtIfilm(pbeta, ncell, basetrack_distance.at(ipairs));
      lateral_sigma = HighlandSigmaAtIfilm(pbeta, ncell, basetrack_distance.at(ipairs));
    }

    if (TMath::Abs(radial_angle_difference.at(ipairs)) < radial_cut_value &&
    	TMath::Abs(lateral_angle_difference.at(ipairs)) < lateral_cut_value ) {
      ll += 2 * TMath::Log(radial_sigma)
	+ radial_angle_difference.at(ipairs) * radial_angle_difference.at(ipairs) / radial_sigma / radial_sigma;
      ll += 2 * TMath::Log(lateral_sigma)
	+ lateral_angle_difference.at(ipairs) * lateral_angle_difference.at(ipairs) / lateral_sigma / lateral_sigma;
    }
    // else break;

    Double_t energy = CalculateEnergyFromPBeta(pbeta, PARTICLE_MASS[particle_id]);
    Double_t beta = CalculateBetaFromPBeta(pbeta, PARTICLE_MASS[particle_id]);

    // Consider energy deposit (upstream -> downstream)
    if (plate_id.at(ipairs) > 16) {
      energy -= direction * ncell * basetrack_distance.at(ipairs) * EnergyDepositIron(beta);
      // Water smearing is considered separately
      energy -= direction * ncell * water_basetrack_distance.at(ipairs) * EnergyDepositWater(beta);
    // } if else (plate_id.at(ipairs) > 15 && plate_id.at(ipairs)%2 == 1) {
    } else {
      energy -= direction * ncell * basetrack_distance.at(ipairs) * EnergyDepositIron(beta);
    }

    if (energy <= PARTICLE_MASS[particle_id]) break;
    pbeta = CalculatePBetaFromEnergy(energy, PARTICLE_MASS[particle_id]);

  }

  return ll;

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

Double_t EnergyDepositIron(Double_t beta) {
  std::vector<Double_t> funcpar = CalculateEnergyDepositIronParameters();
  Double_t enedep = 0;
  for (Int_t i = 0; i < funcpar.size(); i++)
    enedep += funcpar.at(i) * TMath::Power(beta, i);
  return enedep;
}

std::vector<Double_t> CalculateEnergyDepositIronParameters() {
  std::vector<Double_t> return_vec;
   Double_t fit_result[5] = {11.3616,
			    -41.0421,
			    63.9158,
			    -47.3527,
			    13.7393};
  for(Int_t i = 0; i < 5; i++)
    return_vec.push_back(fit_result[i]);
  return return_vec;
}

Double_t EnergyDepositWater(Double_t beta) {
  std::vector<Double_t> funcpar = CalculateEnergyDepositWaterParameters();
  Double_t enedep = 0;
  for (Int_t i = 0; i < funcpar.size(); i++)
    enedep += funcpar.at(i) * TMath::Power(beta, i);
  return enedep;
}

std::vector<Double_t> CalculateEnergyDepositWaterParameters() {
  std::vector<Double_t> return_vec;
  Double_t fit_result[5] = {16.677,
			    -68.1774,
			    114.435,
			    -88.6891,
			    26.3072};
  for(Int_t i = 0; i < 5; i++)
    return_vec.push_back(fit_result[i]);
  return return_vec;
}

Double_t HighlandSigmaAtIfilm(Double_t pbeta, UInt_t ncell, Double_t dz) {
  Double_t beta = 1.;
  Double_t radiation_length = calculate_radiation_length(ncell, dz);
  return scale_factor * 13.6 / pbeta * TMath::Sqrt(radiation_length) * (1. + 0.038 * TMath::Log(radiation_length / beta));
}

Double_t RadialSigmaAtIfilm(Double_t pbeta, UInt_t ncell, Double_t dz, Double_t tangent) {
  Double_t sigma_highland = HighlandSigmaAtIfilm(pbeta, ncell, dz);
  Double_t radial_accuracy = new_radial_tangent_accuracy(tangent, kNinjaIron);
  return std::hypot(sigma_highland, radial_accuracy);
}


Double_t LateralSigmaAtIfilm(Double_t pbeta, UInt_t ncell, Double_t dz, Double_t tangent) {
  Double_t sigma_highland = HighlandSigmaAtIfilm(pbeta, ncell, dz);
  Double_t lateral_accuracy = new_lateral_tangent_accuracy(tangent, kNinjaIron);
  return std::hypot(sigma_highland, lateral_accuracy);
}

std::array<Double_t, 3> ReconstructPBeta(Double_t initial_pbeta,
					 UInt_t ncell,
					 Int_t particle_id,
					 Int_t direction,
					 Double_t radial_cut_value,
					 Double_t lateral_cut_value,
					 Bool_t smear_flag,
					 std::vector<Double_t> basetrack_distance,
					 std::vector<Double_t> water_basetrack_distance,
					 std::vector<Double_t> track_tangent,
					 std::vector<Int_t> plate_id,
					 std::vector<Double_t> radial_angle_difference,
					 std::vector<Double_t> lateral_angle_difference) {

  const UInt_t num_of_param = 8 + 6 * basetrack_distance.size();

  TMinuit *min = new TMinuit(num_of_param); // TMinuit(n), n = number of parameters
  min->SetPrintLevel(-1);
  min->SetFCN(LogLikelihood); // minimize LogLikelihood function
  int ierflg = 0; // output level (0 = normal)
  
  // Parameter Names
  TString parname[num_of_param];
  parname[0] = "Reconstructed pbeta";
  parname[1] = "Ncell";
  parname[2] = "Particle id";
  parname[3] = "Track direction";
  parname[4] = "Radial angle difference cut value";
  parname[5] = "Lateral angle difference cut value";
  parname[6] = "Smear flag";
  parname[7] = "Number of track pairs";
  for (Int_t ipar = 0; ipar < basetrack_distance.size(); ipar++) {
    parname[8 +                                 ipar] = Form("Basetrack distance %d",       ipar);
    parname[8 + 1 * basetrack_distance.size() + ipar] = Form("Water basetrack distance %d", ipar);
    parname[8 + 2 * basetrack_distance.size() + ipar] = Form("Track tangent %d",            ipar);
    parname[8 + 3 * basetrack_distance.size() + ipar] = Form("Plate id %d",                 ipar);
    parname[8 + 4 * basetrack_distance.size() + ipar] = Form("Radial angle difference %d",  ipar);
    parname[8 + 5 * basetrack_distance.size() + ipar] = Form("Lateral angle difference %d", ipar);
  }

  // Initial values
  Double_t vstart[num_of_param];
  vstart[0] = initial_pbeta;
  vstart[1] = ncell;
  vstart[2] = particle_id;
  vstart[3] = direction;
  vstart[4] = radial_cut_value;
  vstart[5] = lateral_cut_value;
  vstart[6] = smear_flag;
  vstart[7] = basetrack_distance.size();
  for (Int_t ipar = 0; ipar < basetrack_distance.size(); ipar++) {
    vstart[8 +                                 ipar] = basetrack_distance.at(ipar);
    vstart[8 + 1 * basetrack_distance.size() + ipar] = water_basetrack_distance.at(ipar);
    vstart[8 + 2 * basetrack_distance.size() + ipar] = track_tangent.at(ipar);
    vstart[8 + 3 * basetrack_distance.size() + ipar] = plate_id.at(ipar);
    vstart[8 + 4 * basetrack_distance.size() + ipar] = radial_angle_difference.at(ipar);
    vstart[8 + 5 * basetrack_distance.size() + ipar] = lateral_angle_difference.at(ipar);
  }
  
  // Step
  Double_t step[num_of_param];
  step[0] = 1.;
  for (Int_t ipar = 1; ipar < num_of_param; ipar++) {
    step[ipar] = 0.;
  }

  // Parameter setting
  // Bethe-Bloch is applicable to pbeta > ~20 (beta > 0.4)
  min->mnparm(0, parname[0], vstart[0], step[0], 20, 5000, ierflg);
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
  arglist[1] = .1; // tolerance
  min->mnexcm("MIGRAD", arglist, 2, ierflg);
  
  Double_t rec_mom, rec_mom_err;
  Int_t fit_status;
  min->GetParameter(0, rec_mom, rec_mom_err);
  fit_status = min->GetStatus();
  std::array<Double_t,3> return_array = {rec_mom, rec_mom_err, (Double_t)fit_status};

  BOOST_LOG_TRIVIAL(trace) << "Reconstructed Momentum = "
			  << rec_mom << " +/- "
			  << rec_mom_err << " [MeV/c]";
  
  delete min;

  return return_array;
 
}


int main (int argc, char *argv[]) {

  logging::core::get()->set_filter
    (
     logging::trivial::severity >= logging::trivial::info
     //logging::trivial::severity >= logging::trivial::trace
     //logging::trivial::severity >= logging::trivial::debug
     );
  
  BOOST_LOG_TRIVIAL(info) << "==========Momentum Reconstruction Start==========";

  if (argc != 6) {
    BOOST_LOG_TRIVIAL(error) << "Usage : " << argv[0]
			     << " <input particle gun MC file name> <output root file name> <ncell> <initial pbeta bias> <true(0)/smear(1)>";
    std::exit(1);
  }
    
  try {

    // Input file
    B2Reader reader(argv[1]);

    // Output file
    TString ofilename = argv[2];
    TFile *ofile = new TFile(ofilename, "recreate");
    TTree *otree = new TTree("tree", "tree");

    Double_t true_pbeta;
    Double_t initial_pbeta;
    Double_t recon_pbeta;
    Double_t recon_pbeta_err;
    Double_t log_likelihood[3][2] = {}; // 3 particles (mu, pi, p), 2 directions (up, down)
    Int_t fit_status[3][2] = {};
    // any other variables?
    otree->Branch("true_pbeta",  &true_pbeta,  "true_pbeta/D");
    otree->Branch("initial_pbeta", &initial_pbeta, "initial_pbeta/D");
    otree->Branch("recon_pbeta", &recon_pbeta, "recon_pbeta/D");
    otree->Branch("recon_pbeta_err", &recon_pbeta_err, "recon_pbeta_err/D");
    otree->Branch("log_likelihood", log_likelihood, "log_likelihood[3][2]/D");
    otree->Branch("fit_status", fit_status, "fit_status[3][2]/I");
    
    Bool_t smear_flag = atoi(argv[5]);
    if (smear_flag)
      BOOST_LOG_TRIVIAL(info) << "========== Smear flag : TRUE ==========";
    else 
      BOOST_LOG_TRIVIAL(info) << "========== Smear flag : FALSE ==========";

    // Start pbeta reconstruction process
    while (reader.ReadNextSpill() > 0) {
      // if (reader.GetEntryNumber() != 6594) continue;
      auto &spill_summary = reader.GetSpillSummary();

      // Get emulsion tracks
      std::vector<const B2EmulsionSummary*> emulsions;
      auto it_emulsion = spill_summary.BeginEmulsion();
      while (const auto *emulsion = it_emulsion.Next()) {
	if (emulsion->GetParentTrackId() == 0) continue; // track id should be assigned
	if (emulsion->GetParentTrack().GetParticlePdg() != 13) continue; // only muon
	if (emulsion->GetFilmType() != B2EmulsionType::kECC) continue; // only ECC films
	if (emulsion->GetEcc() != 4) continue; // only ECC5
	//if (emulsion->GetPlate() <= 14) continue; // only ECC5 w/o Iron ECC
	emulsions.push_back(emulsion);
      }

      if (emulsions.size() <= 0) continue;

      std::sort(emulsions.begin(), emulsions.end(), emulsion_compare);

      // Get true information
      Int_t track_id_tmp_ = emulsions.at(0)->GetParentTrackId();
      Int_t particle_pdg_tmp_ = emulsions.at(0)->GetParentTrack().GetParticlePdg();
      Int_t ecc_tmp_ = emulsions.at(0)->GetEcc();
      if (particle_pdg_tmp_ == 13)
	true_pbeta = CalculatePBetaFromMomentum(emulsions.at(0)->GetMomentum().GetValue().Mag(), MCS_MUON_MASS);
      else if (particle_pdg_tmp_ == 2212)
	true_pbeta = CalculatePBetaFromMomentum(emulsions.at(0)->GetMomentum().GetValue().Mag(), MCS_PROTON_MASS);
      else if (std::abs(particle_pdg_tmp_) == 211)
	true_pbeta = CalculatePBetaFromMomentum(emulsions.at(0)->GetMomentum().GetValue().Mag(), MCS_PION_MASS);

      // When there are in-flight interactions,
      // the event is not suitable for the evaluation of the reconstruction
      Double_t momentum_difference_max = 0.;

      // input parameters for TMinuit
      const Int_t ncell = std::atoi(argv[3]); // Ncell will be sweeped in the next step
      std::vector<Double_t> basetrack_distance = {};
      std::vector<Double_t> water_basetrack_distance = {};
      std::vector<Double_t> track_tangent = {};
      std::vector<Int_t> plate_id = {};
      std::vector<Double_t> radial_angle_difference = {};
      std::vector<Double_t> lateral_angle_difference = {};

      for (Int_t iemulsion_up = 0; iemulsion_up < emulsions.size(); iemulsion_up++) {

	const auto emulsion_up = emulsions.at(iemulsion_up);
	// Only use films across one iron plate
	if (emulsion_up->GetPlate() <= 3) continue;
	if (emulsion_up->GetPlate()%2 == 1 &&
	    emulsion_up->GetPlate() >= 15) continue; 

	TVector3 tangent_up = emulsion_up->GetTangent().GetValue();
	tangent_up = (1. / tangent_up.Z()) * tangent_up;
	TVector3 position_up = emulsion_up->GetAbsolutePosition().GetValue();

	for (Int_t iemulsion_down = iemulsion_up + 1; iemulsion_down < emulsions.size(); iemulsion_down++) {
	  const auto emulsion_down = emulsions.at(iemulsion_down);

	  Int_t plate_difference = emulsion_up->GetPlate() - emulsion_down->GetPlate();

	  if (plate_difference == 2 * ncell - 1) {
	    TVector3 tangent_down = emulsion_down->GetTangent().GetValue();
	    tangent_down = (1. / tangent_down.Z()) * tangent_down;
	    // smear
	    // tangent_up = smear_tangent_vector(tangent_up, kNinjaIron);
	    // tangent_down = smear_tangent_vector(tangent_down, kNinjaIron);

	    TVector3 position_down = emulsion_down->GetAbsolutePosition().GetValue();
	    TVector3 displacement = position_down - position_up;
	    if (smear_flag) { // smear
	      // displacement = smear_distance_vector(displacement, kNinjaIron);
	    }
	    basetrack_distance.push_back(displacement.Mag() / 850.e-3); // Nominal distance = 500 + 350 um

	    if (emulsion_up->GetPlate() >= 18 && iemulsion_down < emulsions.size() - 1) {
	      const auto emulsion_2pl_down = emulsions.at(iemulsion_down + 1);
	      TVector3 position_2pl_down = emulsion_2pl_down->GetAbsolutePosition().GetValue();
	      TVector3 water_displacement = position_2pl_down - position_down;
	      if (smear_flag) { // smear
		// water_displacement = smear_distance_vector(water_displacement, kNinjaWater);
	      }
	      water_basetrack_distance.push_back(water_displacement.Mag() / 2.868); // Nominal distance = 2300 + 109 * 2 + 350 um
	    } else {
	      water_basetrack_distance.push_back(0.);
	    }

	    track_tangent.push_back(tangent_up.Mag());
	    plate_id.push_back(emulsion_up->GetPlate());

	    Double_t angle_difference_radial_new, angle_difference_lateral_new;
	    Double_t tangent_accuracy_radial_new = 0.;
	    Double_t tangent_accuracy_lateral_new = 0.;
	    angle_difference_radial_new = get_angle_difference_radial_new(tangent_up, tangent_down);
	    angle_difference_lateral_new = get_angle_difference_lateral_new(tangent_up, tangent_down);
	    if (smear_flag) { // smear
	      tangent_accuracy_radial_new = new_radial_tangent_accuracy(tangent_up.Mag(), kNinjaIron);
	      tangent_accuracy_lateral_new = new_lateral_tangent_accuracy(tangent_up.Mag(), kNinjaIron);
	    }
	    radial_angle_difference.push_back(gRandom->Gaus(angle_difference_radial_new, tangent_accuracy_radial_new));
	    lateral_angle_difference.push_back(gRandom->Gaus(angle_difference_lateral_new, tangent_accuracy_lateral_new));

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
      if (smear_flag) { // smear 
	vertex_tangent = smear_tangent_vector(vertex_tangent, kNinjaIron);
      }

      Double_t radial_angle_difference_rms = 0.;
      Double_t lateral_angle_difference_rms = 0.;
      for (Int_t ipair = 0; ipair < radial_angle_difference.size(); ipair++) {
	radial_angle_difference_rms += radial_angle_difference.at(ipair) * radial_angle_difference.at(ipair);
	lateral_angle_difference_rms += lateral_angle_difference.at(ipair) * lateral_angle_difference.at(ipair);
      }

      if (smear_flag) {
	radial_angle_difference_rms /= radial_angle_difference.size();
	lateral_angle_difference_rms /= lateral_angle_difference.size();      
      } else {
	radial_angle_difference_rms = (radial_angle_difference_rms + lateral_angle_difference_rms) 
	  / (radial_angle_difference.size() + lateral_angle_difference.size());
      }
      radial_angle_difference_rms = TMath::Sqrt(radial_angle_difference_rms);
      lateral_angle_difference_rms = TMath::Sqrt(lateral_angle_difference_rms);
      if (!smear_flag) lateral_angle_difference_rms = radial_angle_difference_rms; // angle difference are equivalent in true-level

      Double_t radial_cut_value, lateral_cut_value;
      radial_cut_value = 3. * radial_angle_difference_rms;
      lateral_cut_value = 3. * lateral_angle_difference_rms;

      Double_t radiation_length = calculate_radiation_length(ncell, vertex_tangent.Mag());
      initial_pbeta = scale_factor * 13.6 / radial_angle_difference_rms * TMath::Sqrt(radiation_length) * (1. + 0.038 * TMath::Log(radiation_length));
      const Double_t initial_pbeta_bias = std::atof(argv[4]);
      initial_pbeta *= initial_pbeta_bias;

      Int_t particle_id = 0;
      Int_t direction;
      Int_t fit_status_tmp = 0;
      std::array<std::array<Double_t, 3>, 2> result_array = {};
      for (Int_t idirection = 0; idirection < kNumberOfNinjaMcsDirections; idirection++) {

	if (idirection != 0) break;
	direction = MCS_DIRECTION[idirection];

	result_array.at(idirection) = ReconstructPBeta(initial_pbeta, ncell, particle_id, direction,
						       radial_cut_value, lateral_cut_value,
						       smear_flag,
						       basetrack_distance,
						       water_basetrack_distance,
						       track_tangent,
						       plate_id,
						       radial_angle_difference,
						       lateral_angle_difference);
	fit_status[particle_id][idirection] = (Int_t)result_array.at(idirection).at(2);
	log_likelihood[particle_id][idirection] = FuncLogLikelihood(result_array.at(idirection).at(0),
								    ncell, particle_id, direction,
								    radial_cut_value, lateral_cut_value,
								    smear_flag,
								    basetrack_distance,
								    water_basetrack_distance,
								    track_tangent,
								    plate_id,
								    radial_angle_difference,
								    lateral_angle_difference);
      }

      if (log_likelihood[particle_id][0] <= log_likelihood[particle_id][1]) {
	recon_pbeta = result_array.at(0).at(0);
	recon_pbeta_err = result_array.at(0).at(1);
      } else {
	recon_pbeta = result_array.at(1).at(0);
	recon_pbeta_err = result_array.at(1).at(1);
      }

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
