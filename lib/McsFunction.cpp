#include <boost/log/core.hpp>
#include <boost/log/trivial.hpp>
#include <boost/log/expressions.hpp>

#include <iostream>
#include <vector>
#include <array>

#include <B2Const.hh>

#include <TMath.h>
#include <TMinuit.h>
#include <TVector3.h>
#include <TRandom3.h>

#include "McsConst.hpp"
#include "McsFunction.hpp"

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
void NegativeLogLikelihood(Int_t &npar, Double_t *gin, Double_t &f, Double_t *par, Int_t iflag) {
  
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
  std::vector<Double_t > basetrack_distance = {};
  basetrack_distance.resize(number_of_pairs);
  // downstream water basetrack distance
  std::vector<Double_t > water_basetrack_distance = {};
  water_basetrack_distance.resize(number_of_pairs);
  // upstream track tangent
  std::vector<Double_t > track_tangent = {};
  track_tangent.resize(number_of_pairs);
  // upstream plate id
  std::vector<Int_t > plate_id = {};
  plate_id.resize(number_of_pairs);
  // radial angle differences
  std::vector<Double_t > radial_angle_difference = {};
  radial_angle_difference.resize(number_of_pairs);
  // lateral angle differences
  std::vector<Double_t > lateral_angle_difference = {};
  lateral_angle_difference.resize(number_of_pairs);
  for(Int_t ipairs = 0; ipairs < number_of_pairs; ipairs++) {
    basetrack_distance.at(ipairs)       =        par[8                       + ipairs];
    water_basetrack_distance.at(ipairs) =        par[8 + 1 * number_of_pairs + ipairs];
    track_tangent.at(ipairs)            =        par[8 + 2 * number_of_pairs + ipairs];
    plate_id.at(ipairs)                 = (Int_t)par[8 + 3 * number_of_pairs + ipairs];
    radial_angle_difference.at(ipairs)  =        par[8 + 4 * number_of_pairs + ipairs];
    lateral_angle_difference.at(ipairs) =        par[8 + 5 * number_of_pairs + ipairs];
  }
  
  f = FuncNegativeLogLikelihood(pbeta, ncell, particle_id, direction,
				radial_cut_value, lateral_cut_value,
				smear_flag,
				basetrack_distance,
				water_basetrack_distance,
				track_tangent,
				plate_id,
				radial_angle_difference,
				lateral_angle_difference);
  
}

// NegativeLogLikelihoodWater()

// NegativeLogLikelihoodBoth()

Double_t FuncNegativeLogLikelihood(Double_t pbeta,
				   UInt_t ncell,
				   Int_t particle_id,
				   Int_t direction,
				   Double_t radial_cut_value,
				   Double_t lateral_cut_value,
				   Bool_t smear_flag,
				   std::vector<Double_t > basetrack_distance,
				   std::vector<Double_t > water_basetrack_distance,
				   std::vector<Double_t > track_tangent,
				   std::vector<Int_t > plate_id,
				   std::vector<Double_t > radial_angle_difference,
				   std::vector<Double_t > lateral_angle_difference) {
  
  // vector size check
  std::vector<std::size_t > par_vect_sizes = {};
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
  Double_t nll = 0; // constant (= n * log(2pi) / 2) is ignored
  Double_t radial_sigma, lateral_sigma;
  for (Int_t ipairs = 0; ipairs < number_of_pairs; ipairs++) {
    
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
      nll += 2 * TMath::Log(radial_sigma)
      	+ radial_angle_difference.at(ipairs) * radial_angle_difference.at(ipairs) / radial_sigma / radial_sigma;
      nll += 2 * TMath::Log(lateral_sigma)
      	+ lateral_angle_difference.at(ipairs) * lateral_angle_difference.at(ipairs) / lateral_sigma / lateral_sigma;
    }
    
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
  
  return nll;
  
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
  if ( beta > 0.8 )
    return EnergyDepositIronPoly(beta);
  else
    return EnergyDepositIronBetheBloch(beta);
}

Double_t EnergyDepositIronPoly(Double_t beta) {
  Double_t enedep = 0;
  Int_t i = 0;
  for ( auto parameter : IRON_ENEDEP_FUNC_PAR ) {
    enedep += parameter * TMath::Power(beta, i);
    i++;
  }
  return enedep;
}

Double_t EnergyDepositIronBetheBloch(Double_t beta) {
  return IRON_BB_FUNC_PAR.at(0) * ( (TMath::Log(beta * beta / (1 - beta * beta)) + IRON_BB_FUNC_PAR.at(1)) / beta / beta - 1 );
}

Double_t EnergyDepositWater(Double_t beta) {
  if ( beta > 0.8 ) 
    return EnergyDepositWaterPoly(beta);
  else
    return EnergyDepositWaterBetheBloch(beta);
}

Double_t EnergyDepositWaterPoly(Double_t beta) {
  Double_t enedep = 0;
  Int_t i = 0;
  for ( auto parameter : WATER_ENEDEP_FUNC_PAR ) {
    enedep += parameter * TMath::Power(beta, i);
    i++;
  }
  return enedep;
}

Double_t EnergyDepositWaterBetheBloch(Double_t beta) {
  return WATER_BB_FUNC_PAR.at(0) * ( (TMath::Log(beta * beta / (1 - beta * beta)) + WATER_BB_FUNC_PAR.at(1)) / beta / beta - 1 );
}

Double_t HighlandSigmaAtIfilm(Double_t pbeta, UInt_t ncell, Double_t dz) {
  Double_t beta = 1.;
  Double_t radiation_length = CalcRadLength(ncell, dz);
  return MCS_SCALE_FACTOR * 13.6 / pbeta * TMath::Sqrt(radiation_length) * (1. + 0.038 * TMath::Log(radiation_length / beta));
}

Double_t RadialSigmaAtIfilm(Double_t pbeta, UInt_t ncell, Double_t dz, Double_t tangent) {
  Double_t sigma_highland = HighlandSigmaAtIfilm(pbeta, ncell, dz);
  Double_t tangent_xy = TMath::Sqrt(tangent * tangent - 1);
  Double_t radial_precision = RadialAnglePrecision(tangent_xy);
  return TMath::Hypot(sigma_highland, radial_precision);
}

Double_t LateralSigmaAtIfilm(Double_t pbeta, UInt_t ncell, Double_t dz, Double_t tangent) {
  Double_t sigma_highland = HighlandSigmaAtIfilm(pbeta, ncell, dz);
  Double_t tangent_xy = TMath::Sqrt(tangent * tangent - 1);
  Double_t lateral_precision = LateralAnglePrecision(tangent_xy);
  return TMath::Hypot(sigma_highland, lateral_precision);
}

std::pair<Double_t, Double_t > GetMinMaxPlot(Double_t tangent) {

  std::pair<Double_t, Double_t> min_max_pair;
  if (tangent < 0.05) {min_max_pair.first = 0.; min_max_pair.second = 0.05;}
  else if (tangent < 0.15) {min_max_pair.first = 0.05; min_max_pair.second = 0.15;}
  else if (tangent < 0.25) {min_max_pair.first = 0.15; min_max_pair.second = 0.25;}
  else if (tangent < 0.35) {min_max_pair.first = 0.25; min_max_pair.second = 0.35;}
  else if (tangent < 0.45) {min_max_pair.first = 0.35; min_max_pair.second = 0.45;}
  else if (tangent < 0.55) {min_max_pair.first = 0.45; min_max_pair.second = 0.55;}
  else if (tangent < 0.65) {min_max_pair.first = 0.55; min_max_pair.second = 0.65;}
  else if (tangent < 0.75) {min_max_pair.first = 0.65; min_max_pair.second = 0.75;}
  else if (tangent < 0.85) {min_max_pair.first = 0.75; min_max_pair.second = 0.85;}
  else if (tangent < 0.95) {min_max_pair.first = 0.85; min_max_pair.second = 0.95;}
  else if (tangent < 1.05) {min_max_pair.first = 0.95; min_max_pair.second = 1.05;}
  else if (tangent < 1.15) {min_max_pair.first = 1.05; min_max_pair.second = 1.15;}
  else if (tangent < 1.25) {min_max_pair.first = 1.15; min_max_pair.second = 1.25;}
  else if (tangent < 1.35) {min_max_pair.first = 1.25; min_max_pair.second = 1.35;}
  else if (tangent < 1.55) {min_max_pair.first = 1.35; min_max_pair.second = 1.55;}
  else if (tangent < 1.85) {min_max_pair.first = 1.55; min_max_pair.second = 1.85;}
  else if (tangent < 2.15) {min_max_pair.first = 1.85; min_max_pair.second = 2.15;}
  else if (tangent < 2.45) {min_max_pair.first = 2.15; min_max_pair.second = 2.45;}
  else if (tangent < 2.75) {min_max_pair.first = 2.45; min_max_pair.second = 2.75;}
  else if (tangent < 3.05) {min_max_pair.first = 2.75; min_max_pair.second = 3.05;}
  else if (tangent < 3.35) {min_max_pair.first = 3.05; min_max_pair.second = 3.35;}
  else if (tangent < 3.65) {min_max_pair.first = 3.35; min_max_pair.second = 3.65;}
  else if (tangent < 3.95) {min_max_pair.first = 3.65; min_max_pair.second = 3.95;}
  else if (tangent < 4.25) {min_max_pair.first = 3.95; min_max_pair.second = 4.25;}
  else if (tangent < 4.55) {min_max_pair.first = 4.25; min_max_pair.second = 4.55;}
  else if (tangent < 4.85) {min_max_pair.first = 4.55; min_max_pair.second = 4.85;}
  else if (tangent < 5.)   {min_max_pair.first = 4.85; min_max_pair.second = 5.;}

  return min_max_pair;
}

Double_t GetRadialPlot(Double_t tangent) {
  if (tangent < 0.1) return 0.00171;
  else if (tangent < 0.2) return 0.00234;
  else if (tangent < 0.3) return 0.00335;
  else if (tangent < 0.4) return 0.00424;
  else if (tangent < 0.5) return 0.00505;
  else if (tangent < 0.6) return 0.00550;
  else if (tangent < 0.7) return 0.00596;
  else if (tangent < 0.8) return 0.00623;
  else if (tangent < 0.9) return 0.00628;
  else if (tangent < 1.0) return 0.00643;
  else if (tangent < 1.1) return 0.00654;
  else if (tangent < 1.2) return 0.00642;
  else if (tangent < 1.3) return 0.00642;
  else if (tangent < 1.4) return 0.00628;
  else if (tangent < 1.7) return 0.00604;
  else if (tangent < 2.0) return 0.00557;
  else if (tangent < 2.3) return 0.00501;
  else if (tangent < 2.6) return 0.00457;
  else if (tangent < 2.9) return 0.00437;
  else if (tangent < 3.2) return 0.00431;
  else if (tangent < 3.5) return 0.00418;
  else if (tangent < 3.8) return 0.00409;
  else if (tangent < 4.1) return 0.00399;
  else if (tangent < 4.4) return 0.00387;
  else if (tangent < 4.7) return 0.00380;
  else if (tangent < 5.0) return 0.00368;
  else return 0.00368;
}

Double_t GetLateralPlot(Double_t tangent) {
  if (tangent < 0.1) return 0.00131;
  else if (tangent < 0.2) return 0.00130;
  else if (tangent < 0.3) return 0.00131;
  else if (tangent < 0.4) return 0.00129;
  else if (tangent < 0.5) return 0.00129;
  else if (tangent < 0.6) return 0.00119;
  else if (tangent < 0.7) return 0.00110;
  else if (tangent < 0.8) return 0.00102;
  else if (tangent < 0.9) return 0.000952;
  else if (tangent < 1.0) return 0.000889;
  else if (tangent < 1.1) return 0.000822;
  else if (tangent < 1.2) return 0.000765;
  else if (tangent < 1.3) return 0.000707;
  else if (tangent < 1.4) return 0.000663;
  else if (tangent < 1.7) return 0.000587;
  else if (tangent < 2.0) return 0.000496;
  else return 0.000496;
}

Double_t RadialAnglePrecision(Double_t tangent) {
  auto min_max_pair = GetMinMaxPlot(tangent);
  if (min_max_pair.first < 0.05)
    return GetRadialPlot(min_max_pair.first);
  else if (min_max_pair.second > 4.85) 
    return GetRadialPlot(min_max_pair.second);
  else 
    return GetRadialPlot(min_max_pair.first) + 
      (GetRadialPlot(min_max_pair.second) - GetRadialPlot(min_max_pair.first)) / 
      (min_max_pair.second - min_max_pair.first) * (tangent - min_max_pair.first);
}

Double_t LateralAnglePrecision(Double_t tangent) {
  auto min_max_pair = GetMinMaxPlot(tangent);
  if (min_max_pair.first < 0.05)
    return GetLateralPlot(min_max_pair.first);
  else if (min_max_pair.second > 1.85) 
    return GetLateralPlot(min_max_pair.second);
  else 
    return GetLateralPlot(min_max_pair.first) + 
      (GetLateralPlot(min_max_pair.second) - GetLateralPlot(min_max_pair.first)) / 
      (min_max_pair.second - min_max_pair.first) * (tangent - min_max_pair.first);
}

std::array<Double_t, 5> ReconstructPBeta(Double_t initial_pbeta,
					 UInt_t ncell,
					 Int_t particle_id,
					 Int_t direction,
					 Double_t radial_cut_value,
					 Double_t lateral_cut_value,
					 Bool_t smear_flag,
					 std::vector<Double_t > basetrack_distance,
					 std::vector<Double_t > water_basetrack_distance,
					 std::vector<Double_t > track_tangent,
					 std::vector<Int_t > plate_id,
					 std::vector<Double_t > radial_angle_difference,
					 std::vector<Double_t > lateral_angle_difference) {

  const UInt_t num_of_param = 8 + 6 * basetrack_distance.size();

  TMinuit *min = new TMinuit(num_of_param); // TMinuit(n), n = number of parameters
  min->SetPrintLevel(-1);
  min->SetFCN(NegativeLogLikelihood); // minimize NegativeLogLikelihood function
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
  // Pbeta range may be better to set depending on initial_pbeta
  min->mnparm(0, parname[0], vstart[0], step[0], 20, 10000, ierflg);
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
  min->mnexcm("MINOS", arglist, 2, ierflg);
  
  Double_t rec_mom;
  Double_t rec_mom_err, rec_mom_err_plus, rec_mom_err_minus;
  Double_t tmp;
  Int_t fit_status;
  min->GetParameter(0, rec_mom, tmp);
  min->mnerrs(0, rec_mom_err_plus, rec_mom_err_minus, rec_mom_err, tmp);
  fit_status = min->GetStatus();
  std::array<Double_t, 5 > return_array = {rec_mom,
					   rec_mom_err,
					   rec_mom_err_plus,
					   rec_mom_err_minus,
					   (Double_t)fit_status};
  
  delete min;

  return return_array;

}

void PositionAddOffset(TVector3 &absolute_position, int ecc_id) {
  // move to each ECC
  absolute_position.SetX(absolute_position.X()
			 - NINJA_ECC_GAP_X * (1 - ecc_id % 3));
  absolute_position.SetY(absolute_position.Y()
			 - NINJA_ECC_GAP_Y * (ecc_id / 3 - 1));

  // film coordinate to ECC coordinate
  absolute_position.SetX(absolute_position.X()
			 - 0.5 * NINJA_ECC_FILM_XY);
  absolute_position.SetY(absolute_position.Y()
			 + NINJA_ENV_THICK
			 + NINJA_DESIC_THICK
			 - 0.5 * NINJA_DESIC_HEIGHT);
  absolute_position.SetZ(absolute_position.Z()
			 - NINJA_BASE_LAYER_THICK
			 - NINJA_EMULSION_LAYER_THICK
			 - NINJA_ENV_THICK
			 - NINJA_DESIC_THICK
			 + 0.5 * NINJA_DESIC_DEPTH);

  // ECC coordinate to global coordinate
  absolute_position.SetX(absolute_position.X() + NINJA_POS_X + NINJA_ECC_POS_X);
  absolute_position.SetY(absolute_position.Y() + NINJA_POS_Y + NINJA_ECC_POS_Y);
  absolute_position.SetZ(absolute_position.Z() + NINJA_POS_Z + NINJA_ECC_POS_Z);

}

Double_t CalcRadLength(Int_t skip, Double_t dz) {

  if (skip >= MAX_NUM_SKIP)
    throw std::out_of_range("skip should be less than MAX_NUM_SKIP");

  Double_t rad_length = 0.;

  for ( Int_t material = 0; material < kNumberOfNinjaMaterials; material++) {

    int num_layers = 0;

    switch (material) {
    case kNinjaIron : 
      num_layers = skip;
      break;
    case kNinjaWater :
      num_layers = skip - 1;
      break;
    case kNinjaGel :
      num_layers = 4 * skip - 2;
      break;
    case kNinjaBase :
      num_layers = 2 * skip - 1;
      break;
    case kNinjaPacking :
      num_layers = 2 * skip - 2;
      break;
    }
    rad_length += num_layers * dz * MATERIAL_THICK[material] / RAD_LENGTH[material];
  }

  return rad_length;

}

bool EmulsionCompare(const B2EmulsionSummary *lhs, const B2EmulsionSummary *rhs) {
  if ( lhs->GetEcc() != rhs->GetEcc() ) 
    return lhs->GetEcc() < rhs->GetEcc();
  else if ( lhs->GetParentTrackId() != rhs->GetParentTrackId() )
    return lhs->GetParentTrackId() < rhs->GetParentTrackId();
  else
    return lhs->GetPlate() > rhs->GetPlate();
}

void SmearDistanceVector(TVector3 &distance, Int_t material) {

  if ( distance.Mag() < 1e-10 ) return;

  Double_t delta_x = gRandom->Gaus(0., XY_ALIGN_ACCURACY[material]);
  Double_t delta_y = gRandom->Gaus(0., XY_ALIGN_ACCURACY[material]);
  Double_t delta_z = gRandom->Gaus(0., Z_ALIGN_ACCURACY[material]);

  distance.SetX(distance.X() + delta_x);
  distance.SetY(distance.Y() + delta_y);
  distance.SetZ(distance.Z() + delta_z);

}

void SmearTangentVector(TVector3 &tangent) {

  Double_t delta_ay, delta_ax;

  if (std::fabs(tangent.X()) < LateralTangentAccuracy()) {
    delta_ax = - gRandom->Gaus(0., LateralTangentAccuracy());
    delta_ay =   gRandom->Gaus(0., RadialTangentAccuracy(tangent.Y()));
  } else {
    
    Double_t tan_theta = std::hypot(tangent.X(), tangent.Y());
    
    Double_t tan_phi = tangent.Y() / tangent.X();
    Double_t cos_phi = tangent.X() / std::fabs(tangent.X())
      * std::sqrt(1. / (1 + tan_phi * tan_phi));
    Double_t sin_phi = tan_phi * cos_phi;
    
    Double_t delta_radial = gRandom->Gaus(0., RadialTangentAccuracy(tan_theta));
    Double_t delta_lateral = gRandom->Gaus(0., LateralTangentAccuracy());
    
    delta_ax = delta_radial * cos_phi - delta_lateral * sin_phi;
    delta_ay = delta_radial * sin_phi + delta_lateral * cos_phi;    
  }
  
  
  tangent.SetX(tangent.X() + delta_ax);
  tangent.SetY(tangent.Y() + delta_ay);
  
}


Double_t RadialTangentAccuracy(Double_t tangent) {

  if ( tangent < 2.5 ) {
    return std::sqrt(2) * std::sqrt(XY_POSITION_ACCURACY * XY_POSITION_ACCURACY
				     + tangent * tangent * Z_POSITION_ACCURACY * Z_POSITION_ACCURACY)/ 210.e-3;
  } else {
    return 1.64e-2 * (tangent - 2.5) + 2.37e-2;
  }

}

Double_t LateralTangentAccuracy() {
  return std::sqrt(2) / 210.e-3 * XY_POSITION_ACCURACY;
}


Double_t RadialAngleDiffNew(TVector3 tangent_up, TVector3 tangent_down) {

  if (std::fabs(tangent_up.X()) < 1e-2 &&
      std::fabs(tangent_up.Y()) < 1e-2)
    return std::atan(tangent_down.Y()) - std::atan(tangent_up.Y());
  
  double a = ( - tangent_up.X() * tangent_down.X() - tangent_up.Y() * tangent_down.Y()
	       + tangent_up.X() * tangent_up.X()   + tangent_up.Y() * tangent_up.Y())
				  / std::sqrt( ( tangent_up.X() * tangent_up.X() + tangent_up.Y() * tangent_up.Y() ) 
					       * tangent_up.Mag2() );
  double b = (tangent_up.X() * tangent_down.X() + tangent_up.Y() * tangent_down.Y() + 1.)
				  / tangent_up.Mag();
  return std::atan(a / b);
}

Double_t LateralAngleDiffNew(TVector3 tangent_up, TVector3 tangent_down) {
  
  if (std::fabs(tangent_up.X()) < 1e-2 &&
      std::fabs(tangent_up.Y()) < 1e-2)
    return std::atan(tangent_down.X()) - std::atan(tangent_up.X());
  
  double a = ( - tangent_up.Y() * tangent_down.X() + tangent_up.X() * tangent_down.Y() )
				  / std::sqrt(tangent_up.X() * tangent_up.X() + tangent_up.Y() * tangent_up.Y());
  double b = ( tangent_up.X() * tangent_down.X() + tangent_up.Y() * tangent_down.Y() + 1.)
				  / tangent_up.Mag();
  return std::atan(a / b);
}


Int_t GetChainDirection(std::vector<const B2EmulsionSummary*> emulsions, Int_t vertex_pl) {
  return 1;
}

Bool_t IsStopInEccFiducial(std::vector<const B2EmulsionSummary*> emulsions, Int_t direction) {

  switch (direction) {
  case 1 : 
    if ( emulsions.back()->GetFilmPosition().GetValue().Y() > 5. &&
	 emulsions.back()->GetFilmPosition().GetValue().Y() < NINJA_ECC_FILM_XY - 5. &&
	 emulsions.back()->GetFilmPosition().GetValue().X() > 5. &&
	 emulsions.back()->GetFilmPosition().GetValue().X() < NINJA_ECC_FILM_XY - 5. &&\
	 emulsions.back()->GetPlate() > 3)
      return true;
    else 
      return false;
    break;
  case -1 : 
    if ( emulsions.front()->GetFilmPosition().GetValue().Y() > 5. &&
	 emulsions.front()->GetFilmPosition().GetValue().Y() < NINJA_ECC_FILM_XY - 5. &&
	 emulsions.front()->GetFilmPosition().GetValue().X() > 5. &&
	 emulsions.front()->GetFilmPosition().GetValue().X() < NINJA_ECC_FILM_XY - 5. &&
	 emulsions.front()->GetPlate() < 132)
      return true;
    else 
      return false;
    break;
  default : 
    throw std::invalid_argument("Direction should be +/- 1");
  }
}
