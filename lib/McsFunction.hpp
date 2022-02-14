#ifndef MCS_FUNCTION_HPP
#define MCS_FUNCTION_HPP

#include <vector>
#include <array>

#include <TVector3.h>

#include <B2EmulsionSummary.hh>

#include "McsConst.hpp"

///> NLL wrapper for MINUIT SetFCN function
void NegativeLogLikelihood(Int_t &npar, Double_t *gin, Double_t &f, Double_t *par, Int_t iflag);

///> NLL function value for given information
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
				   std::vector<Double_t > laterl_angle_difference);

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
Double_t HighlandSigmaAtIfilm(Double_t pbeta, UInt_t ncell, Double_t dz);
///> Sigma of radial scattering angles
Double_t RadialSigmaAtIfilm(Double_t pbeta, UInt_t ncell, Double_t dz, Double_t tangent);
///> Sigma of lateral scattering angles
Double_t LateralSigmaAtIfilm(Double_t pbeta, UInt_t ncell, Double_t dz, Double_t tangent);

///> Calculate angle precisions for radial and lateral directions
std::pair<Double_t, Double_t > GetMinMaxPlot(Double_t tangent);
Double_t GetRadialPlot(Double_t tangent);
Double_t GetLateralPlot(Double_t tangent);
Double_t RadialAnglePrecision(Double_t tangent);
Double_t LateralAnglePrecision(Double_t tangent);

///> Pbeta reconstruction using MINUIT minimizer
std::array<Double_t, 3> ReconstructPBeta(Double_t initial_pbeta,
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
					 std::vector<Double_t > lateral_angle_differnece);

///> Position offset for ECC to global coordinate
void PositionAddOffset(TVector3 &absolute_position, int ecc_id);

///> Calculate unit radiation length
Double_t CalcRadLength(Int_t skip, Double_t dz);

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


///> MC chain booleans
Bool_t IsStopInEccFiducial(std::vector<const B2EmulsionSummary*> emulsions);

#endif
