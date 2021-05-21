// boost includes
#include <boost/log/core.hpp>
#include <boost/log/trivial.hpp>
#include <boost/log/expressions.hpp>

// ROOT includes
#include <TFile.h>
#include <TMinuit.h>
#include <TMath.h>

// B2 includes
#include <B2Reader.hh>

//#include "/home/t2k/odagawa/NinjaMomentumRecon/src/McsCommon.hpp"
#include "/home/t2k/odagawa/NinjaMomentumRecon/src/McsCommon.cpp"
#include "/home/t2k/odagawa/NinjaMomentumRecon/src/McsAngleMethod.cpp"
#include "MinimizerTest.hpp"

namespace logging = boost::log;

// Function to be minimized
// par[0] : reconstructed momentum obtained from minimum log likelihood
// par[1] : number of pairs of the basetracks
// par[2] : Initial angle (true information or smeared if necessary)
// par[3] : Ncell (number of skipped film between the pair)
// par[4] : direction (1 or -1)
void LogLikelihood(Int_t &npar, Double_t *gin, Double_t &f, Double_t *par, Int_t iflag) {
  // Log likelihood
  Double_t ll = 0.; // constant (= n * log(2pi) / 2) is ignored

  // Momentum
  Double_t momentum = par[0];
  
  for (Int_t ifilm = 0; ifilm < (Int_t)par[1]; ifilm++) {
    momentum = momentum + par[4] * MomentumDrop(par[0], par[2], (Int_t)par[3]);
    Double_t sigma = SigmaAtIfilm(momentum, par[2], (Int_t)par[3]);
    //ll += TMath::Log(sigma) + 0.5 * delta_angle[ifilm] * delta_angle[ifilm] / sigma / sigma;
  }
  
  f = ll;
  
}

Double_t MomentumDrop(Double_t momentum, Double_t angle, Int_t ncell) {
  return 0.;
}

Double_t SigmaAtIfilm(Double_t momentum, Double_t angle, Int_t ncell) {
  Double_t beta = 1.;
  Double_t dz = TMath::Sqrt(TMath::Tan(angle * TMath::DegToRad()) * TMath::Tan(angle * TMath::DegToRad()) + 1.);
  Double_t radiation_length = calculate_radiation_length(ncell, dz);
  return 13.6 / momentum * TMath::Sqrt(radiation_length) * (1. + 0.038 * TMath::Log(radiation_length / beta));
}
  

std::array<Double_t, 2> ReconstructMomentum() {

  TMinuit *min = new TMinuit(1); // TMinuit(n), n = number of parameters
  min->SetPrintLevel(1);
  min->SetFCN(LogLikelihood); // minimize LogLikelihood function
  int ierflg = 0; // output level (0 = normal)
  
  // Set parameters and constants    
  // Initial values
  Double_t vstart[4];
  //vstart[0] = std::atof(argv[0]);
  //vstart[1] = std::atof(argv[1]);
  
  // Step
  Double_t step[4];
  step[0] = 1.;
  step[1] = 0.;
  
  min->mnparm(0, "Reconstructed Momentum", vstart[0], step[0], 300, 2000, ierflg);
  min->mnparm(1, "", vstart[1], step[1], 0, 0, ierflg);
  
  
  // Fix parameters
  min->FixParameter(1);
  min->FixParameter(2);
  min->FixParameter(3);
  
  // array to be used for several steps
  Double_t arglist[10] = {};
  
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

  BOOST_LOG_TRIVIAL(debug) << "Reconstructed Momentum = "
			   << rec_mom << " +/- "
			   << rec_mom_err << " [MeV/c]";
  
  delete min;

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


int main (int argc, char *argv[]) {

  logging::core::get()->set_filter
    (
     logging::trivial::severity >= logging::trivial::debug
     );
  
  BOOST_LOG_TRIVIAL(info) << "==========Momentum Reconstruction Start==========";

  if (argc != 2) {
    BOOST_LOG_TRIVIAL(error) << "Usage : " << argv[0]
			     << " <input particle gun MC file name>";
    std::exit(1);
  }
    
  try {

    B2Reader reader(argv[1]);
    
    while (reader.ReadNextSpill() > 0) {
      
    }
    
    
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
