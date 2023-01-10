void CheckAngleResAngle() {

  //string inputdir = "/home/t2k/odagawa/data/mc_data/particlegun/pbeta_recon/muon/true_level";
  string inputdir = "/home/t2k/odagawa/data/mc_data/particlegun/pbeta_recon/muon/smear";

  const double muon_mass = 105.658;

  vector<double> bias;
  vector<double> width;
  vector<double> tangent;
  vector<double> bias_err;
  vector<double> width_err;
  vector<double> tangent_err;

  //TFile *ofile = new TFile("~/muon_wo_angres_ang_graph.root","recreate");
  // TFile *ofile = new TFile("~/muon_w_angres_ang_graph.root","recreate");
  TFile *ofile = new TFile("~/muon_w_angres_ang_graph_1GeV.root","recreate");

  for ( int i = 0; i < 8; i++ ) {
    for (int j = 0; j < 8; j++ ) {

      stringstream inputfile_ss;
      // inputfile_ss << inputdir << "/test_500_"
      inputfile_ss << inputdir << "/test_1000_"
		   << (int)5* i << "_" << (int)5*j
		   << "_new_rad_lat_iron.root";
      TFile *inputfile = new TFile(inputfile_ss.str().c_str(), "read");
      TTree *inputtree = (TTree*)inputfile->Get("tree");

      TH1D *hist = new TH1D(Form("hist_%d_%d", i, j), "", 200, -1, 1);

      inputtree->Draw(Form("(1/recon_pbeta_candidate[0][0] - 1/true_pbeta) * true_pbeta >> hist_%d_%d", i, j));

      TF1 *gaus = new TF1("gaus", "gaus", -1, 1);
      
      hist->Fit(gaus, "Q", "N", -1 ,1);
      
      ofile->cd();
      hist->Write();

      double angle_x = 5. * i;
      double angle_y = 5. * j;

      bias.push_back(gaus->GetParameter(1));
      width.push_back(gaus->GetParameter(2));
      tangent.push_back(hypot(tan(angle_x * TMath::DegToRad()), tan(angle_y * TMath::DegToRad())));
      bias_err.push_back(gaus->GetParError(1));
      width_err.push_back(gaus->GetParError(2));
      tangent_err.push_back(0.);      

    }
  }

  TGraphErrors *ge_bias = new TGraphErrors(bias.size(),
					   &tangent[0], &bias[0],
					   &tangent_err[0], &bias_err[0]);
  ge_bias->SetName("ge_bias");

  TGraphErrors *ge_width = new TGraphErrors(width.size(),
					    &tangent[0], &width[0],
					    &tangent_err[0], &width_err[0]);
  ge_width->SetName("ge_width");

  ofile->cd();
  ge_bias->Write();
  ge_width->Write();
  ofile->Close();

}
