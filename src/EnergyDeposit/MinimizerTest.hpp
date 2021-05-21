const static Double_t muon_mass = 105.658;// MeV/c2

void LogLikelihood(Int_t &npar, Double_t *gin, Double_t &f, Double_t *par, Int_t iflag);

Double_t MomentumDrop(Double_t momentum, UInt_t ncell, Double_t dz);

std::array<Double_t, 3> CalculateFunctionParameters(Double_t dz);

Double_t SigmaAtIfilm(Double_t momentum, UInt_t ncell, Double_t dz);

std::array<Double_t, 2> ReconstructMomentum(UInt_t ncell, Double_t dz, std::vector<Double_t> angle_difference);

// Not used function
Double_t CalculateRadLenUnit(UInt_t ncell);
