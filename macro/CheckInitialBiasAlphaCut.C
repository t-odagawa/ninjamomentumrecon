void CheckInitialBiasAlphaCut() {

  gStyle->SetOptStat(0);
  TString dirname = "/home/t2k/odagawa/data/mc_data/particlegun/pbeta_recon/test20210726/";

  Double_t true_pbeta_1, true_pbeta_2;
  Double_t recon_pbeta_1, recon_pbeta_2;
  Double_t initial_pbeta_1, initial_pbeta_2;

  TF1 *gaus = new TF1("gaus", "gaus", -1., 1.);

  std::vector<Double_t> pbeta_vec = {};
  std::vector<Double_t> bias_vec = {};
  std::vector<Double_t> pbeta_err_vec = {};
  std::vector<Double_t> bias_err_vec = {};
 
  TString pdfname = "~/mcs_bias_plot_alpha_cut.pdf";
  TCanvas *c = new TCanvas("c","c");
  c->SaveAs(pdfname + "[", "pdf");

  for (Int_t ipbeta = 0; ipbeta < 30; ipbeta++) {
    Int_t pbeta = 50 * (ipbeta + 1);

    if (pbeta < 300) continue;

    TString filename1 = Form("test_%d_0_0_new_rad_1.0.root", pbeta);
    filename1 = dirname + filename1;
    TFile *file1 = new TFile(filename1, "read");
    TTree *tree1 = (TTree*)file1->Get("tree");
    tree1->SetBranchAddress("true_pbeta", &true_pbeta_1);
    tree1->SetBranchAddress("recon_pbeta", &recon_pbeta_1);
    tree1->SetBranchAddress("initial_pbeta", &initial_pbeta_1);

    TString filename2 = Form("test_%d_0_0_new_rad_1.1.root", pbeta);
    filename2 = dirname + filename2;
    TFile *file2 = new TFile(filename2, "read");
    TTree *tree2 = (TTree*)file2->Get("tree");
    tree2->SetBranchAddress("true_pbeta", &true_pbeta_2);
    tree2->SetBranchAddress("recon_pbeta", &recon_pbeta_2);
    tree2->SetBranchAddress("initial_pbeta", &initial_pbeta_2);


    TH1D *hist = new TH1D("hist", "hist", 200, -1.0, 1.0);

    for (Int_t ientry = 0; ientry < tree1->GetEntries(); ientry++) {
      tree1->GetEntry(ientry);
      tree2->GetEntry(ientry);

      if (std::fabs(recon_pbeta_1 / recon_pbeta_2 - 1.) < 1e-3) {
	hist->Fill((1 / recon_pbeta_1 - 1 / true_pbeta_1) * true_pbeta_1);
      }

    }

    c->cd();
    hist->Draw("");
    hist->Fit(gaus);
    c->SaveAs(pdfname, "pdf");
    pbeta_vec.push_back(true_pbeta_1);
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
