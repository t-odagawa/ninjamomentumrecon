#include <fstream>
#include <iostream>
#include <string>
#include <vector>

#include <TFile.h>
#include <TGraphErrors.h>

#include "McsClass.hpp"

int main (int argc, char* argv[]) {
  
  if ( argc != 4 ) {
    std::cerr << "Usage : " << argv[0]
	      << " <input momch file1> <input momch file2> <output root file>" << std::endl;
    return 1;
  }

  std::string ifilename1 = argv[1];
  std::string ifilename2 = argv[2];

  auto ev_vec1 = Momentum_recon::ReadEventInformationBin(ifilename1);
  auto ev_vec2 = Momentum_recon::ReadEventInformationBin(ifilename2);
  
  if ( ev_vec1.size() != ev_vec2.size() ) {
    std::cerr << "Vector size is not consistent" << std::endl;
    return 1;
  }

  TFile *ofile = new TFile(argv[3], "recreate");
  TGraphErrors *g_muon_mom = new TGraphErrors();
  g_muon_mom->SetMarkerColor(kRed);
  g_muon_mom->SetMarkerStyle(20);

  for ( int i = 0; i < ev_vec1.size(); i++ ) {
    auto ev1 = ev_vec1.at(i);
    auto ev2 = ev_vec2.at(i);

    double mom1, mom2;

    for ( auto chain : ev1.chains ) {
      if ( chain.particle_flag != 13 ) continue;
      mom1 = chain.ecc_mcs_mom[0];
    }
    for ( auto chain : ev2.chains ) {
      if ( chain.particle_flag != 13 ) continue;
      mom2 = chain.ecc_mcs_mom[0];
    }

    g_muon_mom->SetPoint(i, mom1, mom2);

  }

  ofile->cd();
  g_muon_mom->Write("g_muon_mom");
  ofile->Close();

  return 0;

}
