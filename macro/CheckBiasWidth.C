void CheckBiasWidth() {

  gStyle->SetOptStat(0);

  TString dirname = "/home/t2k/odagawa/data/mc_data/particlegun/pbeta_recon/muon/true_level/";
  Int_t mode = 0; // 0 : momentum, 1 : angle

  Double_t true_pbeta;
  Double_t recon_pbeta;
  Double_t initial_pbeta;

  TF1 *gaus = new TF1("gaus", "gaus", -1., 1.);

  std::vector<Double_t > pbeta_vec = {};
  std::vector<Double_t > bias_vec = {};
  std::vector<Double_t > width_vec = {};
  std::vector<Double_t > pbeta_err_vec = {};
  std::vector<Double_t > bias_err_vec = {};
  std::vector<Double_t > width_err_vec = {};

  TString ofilename = "/home/t2k/odagawa/mcs_bias_width_momentum.root";
  TFile *ofile = new TFile(ofilename, "recreate"); 

  if ( mode == 0 ) {
    
    for ( Int_t imom = 0; imom < 40; imom++ ) {
      Int_t momentum = 50 * (imom + 1);
      if ( momentum < 300 ) continue;

      TString filename = dirname + Form("test_%d_0_0_new_rad_lat_iron.root", momentum);
      TFile *file = new TFile(filename, "read");
      TTree *tree = (TTree*)file->Get("tree");
      
      tree->SetBranchAddress("true_pbeta", &true_pbeta);
      tree->SetBranchAddress("recon_pbeta", &recon_pbeta);

      tree->GetEntry(0);

      TH1D *hist = new TH1D("hist", Form("1/p#beta residual (%d MeV/c w/ S tuning);(1/p#beta_{recon}-1/p#beta_{true})p#beta_{true};Events", momentum), 200, -1.0, 1.0);

      tree->Draw("(1 / recon_pbeta - 1 / true_pbeta) * true_pbeta >> hist", "", "");
      hist->Fit(gaus, "", "", -0.5, 0.5);
      pbeta_vec.push_back(true_pbeta);
      pbeta_err_vec.push_back(0.);
      bias_vec.push_back(gaus->GetParameter(1));
      bias_err_vec.push_back(gaus->GetParError(1));
      width_vec.push_back(gaus->GetParameter(2));
      width_err_vec.push_back(gaus->GetParError(2));

    }

    TGraphErrors *ge_bias = new TGraphErrors(pbeta_vec.size(),
					     &pbeta_vec[0], &bias_vec[0],
					     &pbeta_err_vec[0], &bias_err_vec[0]);
    ge_bias->SetMarkerStyle(20);
    ge_bias->SetTitle("Bias (w/ S tune);p#beta [MeV/c];Bias");
    // TF1 *line = new TF1("line", "[0]", 300., 2000.);
    // ge->Fit(line, "", "", 400., 1900.);    

    TGraphErrors *ge_width = new TGraphErrors(pbeta_vec.size(),
					      &pbeta_vec[0], &width_vec[0],
					      &pbeta_err_vec[0], &width_err_vec[0]);
    ge_width->SetMarkerStyle(20);
    ge_width->SetTitle("Width (w/ S tune);p#beta [MeV/c];Width");

    ge_bias->SetName("ge_bias");
    ge_width->SetName("ge_width");
    ofile->cd();
    ge_bias->Write();
    ge_width->Write();
    ofile->Close();

  }

}
