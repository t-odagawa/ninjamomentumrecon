void CheckInitialBias() {

  gStyle->SetOptStat(0);
  //TString dirname = "/home/t2k/odagawa/data/mc_data/particlegun/pbeta_recon/test20210726/";
  //TString dirname = "/home/t2k/odagawa/data/mc_data/particlegun/pbeta_recon/test20210908/before_s_tune/";
  TString dirname = "/home/t2k/odagawa/data/mc_data/particlegun/pbeta_recon/test20210908/true_level/";
  //TString dirname = "/home/t2k/odagawa/data/mc_data/particlegun/pbeta_recon/test20210908/before_s_tune_larger_initial/";

  Double_t true_pbeta;
  Double_t recon_pbeta;
  Double_t initial_pbeta;

  TF1 *gaus = new TF1("gaus", "gaus", -1., 1.);

  std::vector<Double_t> pbeta_vec = {};
  std::vector<Double_t> bias_vec = {};
  std::vector<Double_t> pbeta_err_vec = {};
  std::vector<Double_t> bias_err_vec = {};
 
  //TString pdfname = "~/mcs_bias_plot.pdf";
  TString pdfname = "~/mcs_bias_plot_after_tune.pdf";
  //TString pdfname = "~/mcs_bias_plot_larger_initial.pdf";
  TCanvas *c = new TCanvas("c","c");
  c->SaveAs(pdfname + "[", "pdf");

  for (Int_t imom = 0; imom < 40; imom++) {
    Int_t momentum = 50 * (imom + 1);

    if (momentum < 300) continue;

    // TString filename = Form("test_%d_0_0_new_rad_1.0.root", pbeta);
    TString filename = Form("test_%d_0_0_new_rad_lat_iron_true_after_tune.root", momentum);
    //TString filename = Form("test_%d_0_0_new_rad_lat_iron_true_before_tune.root", momentum);
    //TString filename = Form("test_%d_0_0_new_rad_lat_iron_true_before_tune_1.1.root", momentum);
    filename = dirname + filename;
    TFile *file = new TFile(filename, "read");
    TTree *tree = (TTree*)file->Get("tree");

    tree->SetBranchAddress("true_pbeta", &true_pbeta);
    tree->SetBranchAddress("recon_pbeta", &recon_pbeta);
    tree->SetBranchAddress("initial_pbeta", &initial_pbeta);

    tree->GetEntry(0);

    // TH1D *hist = new TH1D("hist", Form("1/p#beta residual (%d MeV/c w/o S tuning);(1/p#beta_{recon}-1/p#beta_{true})p#beta_{true};Events", momentum), 200, -1.0, 1.0);
    TH1D *hist = new TH1D("hist", Form("1/p#beta residual (%d MeV/c w/ S tuning);(1/p#beta_{recon}-1/p#beta_{true})p#beta_{true};Events", momentum), 200, -1.0, 1.0);

    c->cd();
    tree->Draw("(1 / recon_pbeta - 1 / true_pbeta) * true_pbeta >> hist", "", "");
    hist->Fit(gaus, "", "", -0.5, 0.5);
    c->SaveAs(pdfname, "pdf");
    pbeta_vec.push_back(true_pbeta);
    pbeta_err_vec.push_back(0.);
    bias_vec.push_back(gaus->GetParameter(1));
    bias_err_vec.push_back(gaus->GetParError(1));

  }

  TGraphErrors *ge = new TGraphErrors(pbeta_vec.size(),
				      &pbeta_vec[0], &bias_vec[0],
				      &pbeta_err_vec[0], &bias_err_vec[0]);
  ge->SetMarkerStyle(20);
  //ge->SetTitle("Bias dependence on initial value (w/o S tune);p#beta [MeV/c];Bias");
  ge->SetTitle("Bias dependence on initial value (w/ S tune);p#beta [MeV/c];Bias");
  ge->Draw("AP");

  TF1 *line = new TF1("line", "[0]", 300., 2000.);
  ge->Fit(line, "", "", 400., 1900.);

  c->SaveAs(pdfname, "pdf");
  c->SaveAs(pdfname + "]", "pdf");

}
