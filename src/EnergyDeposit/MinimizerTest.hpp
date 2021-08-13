const static Double_t muon_mass = 105.658;// MeV/c2
const static Double_t scale_factor = 1.0348;

void LogLikelihood(Int_t &npar, Double_t *gin, Double_t &f, Double_t *par, Int_t iflag);

Double_t FuncLogLikelihood(Double_t pbeta,
			   UInt_t ncell,
			   Int_t particle_id,
			   Int_t direction,
			   Double_t radial_cut_value,
			   Double_t laterla_cut_value,
			   Bool_t smear_flag,
			   std::vector<Double_t> basetrack_distance,
			   std::vector<Double_t> water_basetrack_distance,
			   std::vector<Double_t> track_tangent,
			   std::vector<Int_t> plate_id,
			   std::vector<Double_t> radial_angle_difference,
			   std::vector<Double_t> lateral_angle_difference);

Double_t CalculateBetaFromPBeta(Double_t pbeta, Double_t mass);

Double_t CalculateEnergyFromMomentum(Double_t momentum, Double_t mass);

Double_t CalculateEnergyFromPBeta(Double_t pbeta, Double_t mass);

Double_t CalculatePBetaFromEnergy(Double_t energy, Double_t mass);

Double_t CalculatePBetaFromMomentum(Double_t momentum, Double_t mass);

Double_t CalculateMomentumFromPBeta(Double_t pbeta, Double_t mass);

Double_t CalculateMomentumFromEnergy(Double_t energy, Double_t mass);

Double_t EnergyDepositIron(Double_t beta);

std::vector<Double_t> CalculateEnergyDepositIronParameters();

Double_t EnergyDepositWater(Double_t beta);

std::vector<Double_t> CalculateEnergyDepositWaterParameters();

Double_t HighlandSigmaAtIfilm(Double_t pbeta, UInt_t ncell, Double_t dz);
Double_t RadialSigmaAtIfilm(Double_t pbeta, UInt_t ncell, Double_t dz, Double_t tangent);
Double_t LateralSigmaAtIfilm(Double_t pbeta, UInt_t ncell, Double_t dz, Double_t tangent);


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
					 std::vector<Double_t> lateral_angle_difference);
