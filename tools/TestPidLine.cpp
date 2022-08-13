#include <map>
#include <string>

#include <TString.h>
#include <TCanvas.h>
#include <TGraph.h>

#include "PidData.hpp"
#include "PidFunction.hpp"

int main (int argc, char* argv[]) {

  TString canvasname = argv[1];
  TCanvas *c = new TCanvas("c", "c");
  c->Print(canvasname + "[", "pdf");

  const std::string data_file_dir = argv[2];
  PidData pid_data_(data_file_dir);
  PidFunction pid_function_(pid_data_);

  std::map<int, std::map<double, Pid_data_ns::VphPionMip > > map;
  pid_data_.GetPionMipParamMap(map);

  for ( auto itr = map.begin(); itr != map.end(); itr++ ) {

    int i_pi = 0;
    int i_p = 0;
    TGraph *g_pi = new TGraph();
    g_pi->SetMarkerStyle(20);
    g_pi->SetMarkerSize(0.5);
    g_pi->SetMarkerColor(kRed);
    TGraph *g_p = new TGraph();
    g_p->SetMarkerStyle(20);
    g_p->SetMarkerSize(0.5);
    g_p->SetMarkerColor(kBlue);

    auto each_map = itr->second;
    auto param = each_map.begin()->second;
    double tangent = (param.ang_min + param.ang_max) / 2.;

    for ( double pb = 10; pb < 1500; pb = pb + 10 ) {
      for ( double vph = 0; vph < 350; vph = vph + 5  ) {
	double proton_likelihood, pion_likelihood;
	pid_function_.CalcPartnerLikelihood(vph, pb, tangent,
					    pion_likelihood, proton_likelihood);
	
	int recon_particle_id = pid_function_.GetReconPid(vph, pb, tangent,
							  pion_likelihood, proton_likelihood);

	if ( recon_particle_id == 211 ) {
	  g_pi->SetPoint(i_pi, pb, vph);
	  i_pi++;
	}
	else if ( recon_particle_id == 2212 ) {
	  g_p->SetPoint(i_p, pb, vph);
	  i_p++;
	}
      }
    }

    c->DrawFrame(0,0,1500,350);
    g_pi->Draw("P");
    g_p->Draw("PSAME");
    c->Print(canvasname, "pdf");
  }
  
  c->Print(canvasname + "]", "pdf");

  return 0;

}
