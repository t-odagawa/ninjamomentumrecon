const static Double_t muon_mass = 105.658;// MeV/c2

void LogLikelihood(Int_t &npar, Double_t *gin, Double_t &f, Double_t *par, Int_t iflag);

Double_t CalculateBetaFromPBeta(Double_t pbeta, Double_t mass);

Double_t CalculateEnergyFromMomentum(Double_t momentum, Double_t mass);

Double_t CalculateEnergyFromPBeta(Double_t pbeta, Double_t mass);

Double_t CalculatePBetaFromEnergy(Double_t energy, Double_t mass);

Double_t CalculatePBetaFromMomentum(Double_t momentum, Double_t mass);

Double_t CalculateMomentumFromPBeta(Double_t pbeta, Double_t mass);

Double_t CalculateMomentumFromEnergy(Double_t energy, Double_t mass);

Double_t BetheBlochIron(Double_t beta);

std::array<Double_t, 3> CalculateBetheBlochIronParameters();

Double_t BetheBlochWater(Double_t beta);

std::array<Double_t, 3> CalculateBetheBlochWaterParameters();

Double_t SigmaAtIfilm(Double_t momentum, UInt_t ncell, Double_t dz);

std::array<Double_t, 2> ReconstructPBeta(UInt_t ncell, Double_t dz, std::vector<Double_t> angle_difference);

// Not used function
Double_t CalculateRadLenUnit(UInt_t ncell);
Double_t MomentumDrop(Double_t momentum, UInt_t ncell, Double_t dz);
std::array<Double_t, 3> CalculateFunctionParameters(Double_t dz);

