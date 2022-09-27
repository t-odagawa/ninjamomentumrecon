#ifndef MCS_FUNCTION_HPP
#define MCS_FUNCTION_HPP

#include <vector>
#include <array>

#include <TVector3.h>

#include <B2EmulsionSummary.hh>

#include "McsConst.hpp"

///> NLL wrapper for MINUIT SetFCN function
// par[0]             : reconstructed pbeta obtained from minimum log likelihood
// par[1]             : Ncell (number of skipped film between the pair)
// par[2]             : particle id (muon, charged pion, proton)
// par[3]             : direction (1 or -1)
// par[4]             : radial angle difference cut value
// par[5]             : lateral angle difference cut value
// par[6]             : radial angle difference water cut value
// par[7]             : lateral angle difference water cut value
// par[8]             : true(0)/smear(1) flag
// par[9]             : reconstruction material mode
// par[10]            : radial lateral mode
// par[11]            : number of pairs of the basetracks ( = N )
// par[   12 -  N+11] : basetrack distance (used for energy deposition and radiation length calculation)
// par[ N+12 - 2N+11] : downstream water basetrack distance (used for energy deposition)
// par[2N+12 - 3N+11] : upstream track tangent
// par[3N+12 - 4N+11] : upstream plate id
// par[4N+12 - 5N+11] : downstream plate id
// par[5N+12 - 6N+11] : radial angle differences between basetracks
// par[6N+12 - 7N+11] : lateral angle differences between basetracks
void NegativeLogLikelihood(Int_t &npar, Double_t *gin, Double_t &f, Double_t *par, Int_t iflag);

///> NLL function value for given information
Double_t FuncNegativeLogLikelihood(Double_t pbeta,
				   UInt_t ncell,
				   Int_t particle_id,
				   Int_t direction,
				   Double_t radial_cut_value,
				   Double_t lateral_cut_value,
				   Double_t radial_cut_value_water,
				   Double_t lateral_cut_value_water,
				   Bool_t smear_flag,
				   Int_t material_mode,
				   Int_t radlat_mode,
				   std::vector<Double_t > basetrack_distance,
				   std::vector<Double_t > basetrack_distance_water,
				   std::vector<Double_t > track_tangent,
				   std::vector<Int_t > plate_id,
				   std::vector<Int_t > plate_id_next,
				   std::vector<Double_t > radial_angle_difference,
				   std::vector<Double_t > laterl_angle_difference);

void GetSigmaAtEachStep(Double_t pbeta,
			UInt_t ncell,
			Bool_t smear_flag,
			Int_t material,
			Double_t basetrack_distance,
			Double_t track_tangent,
			Double_t &radial_sigma,
			Double_t &lateral_sigma);

Double_t FuncNegativeLogLikelihoodWater(Double_t pbeta, UInt_t ncell,
					Int_t particle_id, Int_t direction,
					Double_t radial_cut_value_water,
					Double_t lateral_cut_value_water,
					Bool_t smear_flag,
					Int_t radlat_mode,
					std::vector<Double_t > basetrack_distance,
					std::vector<Double_t > basetrack_distance_water,
					std::vector<Double_t > track_tangent,
					std::vector<Int_t > plate_id,
					std::vector<Int_t > plate_id_next,
					std::vector<Double_t > radial_angle_difference,
					std::vector<Double_t > lateral_angle_difference);

Double_t FuncNegativeLogLikelihoodIron(Double_t pbeta, UInt_t ncell,
				       Int_t particle_id, Int_t direction,
				       Double_t radial_cut_value, Double_t lateral_cut_value,
				       Bool_t smear_flag,
				       Int_t radlat_mode,
				       std::vector<Double_t > basetrack_distance,
				       std::vector<Double_t > basetrack_distance_water,
				       std::vector<Double_t > track_tangent,
				       std::vector<Int_t > plate_id,
				       std::vector<Int_t > plate_id_next,
				       std::vector<Double_t > radial_angle_difference,
				       std::vector<Double_t > lateral_angle_difference);

Double_t FuncNegativeLogLikelihoodCombo(Double_t pbeta, UInt_t ncell,
					Int_t particle_id, Int_t direction,
					Double_t radial_cut_value, Double_t lateral_cut_value,
					Double_t radial_cut_value_water,
					Double_t lateral_cut_value_water,
					Bool_t smear_flag,
					Int_t radlat_mode,
					std::vector<Double_t > basetrack_distance,
					std::vector<Double_t > basetrack_distance_water,
					std::vector<Double_t > track_tangent,
					std::vector<Int_t > plate_id,
					std::vector<Int_t > plate_id_next,
					std::vector<Double_t > radial_angle_difference,
					std::vector<Double_t > lateral_angle_difference);

void GetBasetrackDistancePair(int plate, int plate_next, int direction,
			      TVector3 up_position, TVector3 down_position,
			      TVector3 &basetrack_distance_pair,
			      TVector3 &basetrack_distance_water_pair,
			      Bool_t smear_flag);

///> pbeta -> beta
Double_t CalculateBetaFromPBeta(Double_t pbeta, Double_t mass);
///> p -> E
Double_t CalculateEnergyFromMomentum(Double_t momentum, Double_t mass);
///> pbeta -> E
Double_t CalculateEnergyFromPBeta(Double_t pbeta, Double_t mass);
///> E -> pbeta
Double_t CalculatePBetaFromEnergy(Double_t pbeta, Double_t mass);
///> p -> pbeta
Double_t CalculatePBetaFromMomentum(Double_t momentum, Double_t mass);
///> pbeta -> p
Double_t CalculateMomentumFromPBeta(Double_t pbeta, Double_t mass);
///> E -> p
Double_t CalculateMomentumFromEnergy(Double_t energy, Double_t mass);


///> Iron energy deposition calculation
Double_t EnergyDepositIron(Double_t beta);
Double_t EnergyDepositIronPoly(Double_t beta);
Double_t EnergyDepositIronBetheBloch(Double_t beta);
std::vector<Double_t > CalculateEnergyDepositIronParameters();

///> Water energy deposition calculation
Double_t EnergyDepositWater(Double_t beta);
Double_t EnergyDepositWaterPoly(Double_t beta);
Double_t EnergyDepositWaterBetheBloch(Double_t beta);
std::vector<Double_t > CalculateEnergyDepositWaterParameters();

///> Sigma of scattering angles due to naive MCS
Double_t HighlandSigmaAtIfilm(Double_t pbeta, UInt_t ncell, Double_t dz, Int_t material);
///> Sigma of radial scattering angles
Double_t RadialSigmaAtIfilm(Double_t pbeta, UInt_t ncell, Double_t dz,
			    Int_t material, Double_t tangent);
///> Sigma of lateral scattering angles
Double_t LateralSigmaAtIfilm(Double_t pbeta, UInt_t ncell, Double_t dz,
			     Int_t material, Double_t tangent);

///> Calculate angle precisions for radial and lateral directions
std::pair<Double_t, Double_t > GetMinMaxPlot(Double_t tangent);
Double_t GetRadialPlot(Double_t tangent);
Double_t GetRadialErrPlot(Double_t tangent);
Double_t GetLateralPlot(Double_t tangent);
Double_t GetLateralErrPlot(Double_t tangent);
Double_t RadialAnglePrecision(Double_t tangent);
Double_t LateralAnglePrecision(Double_t tangent);

///> Pbeta reconstruction using MINUIT minimizer
std::array<Double_t, 5> ReconstructPBeta(Double_t initial_pbeta,
					 UInt_t ncell,
					 Int_t particle_id,
					 Int_t direction,
					 Double_t radial_cut_value,
					 Double_t lateral_cut_value,
					 Double_t radial_cut_value_water,
					 Double_t lateral_cut_value_water,
					 Bool_t smear_flag,
					 Int_t material_mode,
					 Int_t radlat_mode,
					 std::vector<Double_t > basetrack_distance,
					 std::vector<Double_t > basetrack_distance_water,
					 std::vector<Double_t > track_tangent,
					 std::vector<Int_t > plate_id,
					 std::vector<Int_t > plate_id_next,
					 std::vector<Double_t > radial_angle_difference,
					 std::vector<Double_t > lateral_angle_differnece);

///> Position offset for ECC to global coordinate
void PositionAddOffset(TVector3 &absolute_position, int ecc_id);

///> Calculate unit radiation length
Double_t CalcRadLength(Int_t skip, Double_t dz, Int_t material);
Double_t CalcRadLengthIron(Int_t skip, Double_t dz);
Double_t CalcRadLengthWater(Int_t skip, Double_t dz);

///> Emulsion summary comparator
bool EmulsionCompare(const B2EmulsionSummary *lhs, const B2EmulsionSummary *rhs);

///> Smearing functions
void SmearDistanceVector(TVector3 &distance, Int_t material);
void SmearTangentVector(TVector3 &tangent);

///> Tangent accuracies
Double_t RadialTangentAccuracy(Double_t tangent);
Double_t LateralTangentAccuracy();

///> New angle differences
Double_t RadialAngleDiffNew(TVector3 tangent_up, TVector3 tangent_down);
Double_t LateralAngleDiffNew(TVector3 tangent_up, TVector3 tangent_down);


///> MC chain functions
Int_t GetChainDirection(std::vector<const B2EmulsionSummary*> emulsions, Int_t vertex_pl);
Bool_t IsStopInEccFiducial(std::vector<const B2EmulsionSummary*> emulsions, Int_t direction);

#endif
