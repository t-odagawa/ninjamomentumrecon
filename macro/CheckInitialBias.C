void CheckInitialBias() {

  gStyle->SetOptStat(0);
  TString dirname = "/home/t2k/odagawa/data/mc_data/particlegun/pbeta_recon/test20210726/";

  Double_t true_pbeta;
  Double_t recon_pbeta;
  Double_t initial_pbeta;

  TF1 *gaus = new TF1("gaus", "gaus", -1., 1.);

  std::vector<Double_t> pbeta_vec = {};
  std::vector<Double_t> bias_vec = {};
  std::vector<Double_t> pbeta_err_vec = {};
  std::vector<Double_t> bias_err_vec = {};
 
  TString pdfname = "~/mcs_bias_plot.pdf";
  TCanvas *c = new TCanvas("c","c");
  c->SaveAs(pdfname + "[", "pdf");

  for (Int_t ipbeta = 0; ipbeta < 30; ipbeta++) {
    Int_t pbeta = 50 * (ipbeta + 1);

    if (pbeta < 300) continue;

    TString filename = Form("test_%d_0_0_new_rad_1.0.root", pbeta);
    filename = dirname + filename;
    TFile *file = new TFile(filename, "read");
    TTree *tree = (TTree*)file->Get("tree");
    tree->SetBranchAddress("true_pbeta", &true_pbeta);
    tree->SetBranchAddress("recon_pbeta", &recon_pbeta);
    tree->SetBranchAddress("initial_pbeta", &initial_pbeta);

    tree->GetEntry(0);

    TH1D *hist = new TH1D("hist", "hist", 200, -1.0, 1.0);

    c->cd();
    tree->Draw("(1 / recon_pbeta - 1 / true_pbeta) * true_pbeta >> hist", "", "");
    hist->Fit(gaus);
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
  ge->SetTitle("Bias dependence on initial value;p#beta [MeV/c];Bias");
  ge->Draw("AP");

  TF1 *line = new TF1("line", "[0]", 0., 1600.);
  ge->Fit(line,"","",500., 1500.);

  c->SaveAs(pdfname, "pdf");
  c->SaveAs(pdfname + "]", "pdf");

}
