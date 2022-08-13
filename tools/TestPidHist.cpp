#include <map>
#include <string>

#include <TFile.h>
#include <TH1D.h>

#include "PidData.hpp"

int main (int argc, char* argv[]) {

  TFile *ofile = new TFile(argv[1], "recreate");

  const std::string data_file_dir = argv[2];
  
  PidData pid_data_(data_file_dir);

  std::map<int, TH1D* > hist_map;
  pid_data_.GetPionPdfHistograms(hist_map);

  ofile->cd();
  for ( auto itr = hist_map.begin(); itr != hist_map.end(); itr++ ) {
    (*itr).second->Write();
  }
  ofile->Close();

  return 0;
}
