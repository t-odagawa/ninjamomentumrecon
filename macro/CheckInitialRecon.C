void CheckInitialRecon() {

  gStyle->SetOptStat(0);

  TString filename = "/home/t2k/odagawa/data/plots/mcs_test/test_500_0_0_new_rad.root";
  TFile *file = new TFile(filename, "read");
  TTree *tree = (TTree*)file->Get("tree");

  Double_t true_pbeta;
  Double_t recon_pbeta;
  Double_t initial_pbeta;

  tree->SetBranchAddress("true_pbeta", &true_pbeta);
  tree->SetBranchAddress("recon_pbeta", &recon_pbeta);
  tree->SetBranchAddress("initial_pbeta", &initial_pbeta);

  Double_t x1, y1, x2, y2;
  x1 = 50; y1 = 50; x2 = 750; y2 = 750;

  TH2D *hist_initial_recon = new TH2D("hist_initial_recon", "Initial vs reconstructed pbeta;Initial p#beta [MeV/c];Recon. p#beta",
				      100, x1, x2, 100, y1, y2);

  tree->Draw("recon_pbeta : initial_pbeta >> hist_initial_recon", "", "colz");

  TLine *line = new TLine(x1, y1, x2, y2);
  line->SetLineColor(kRed);
  line->SetLineWidth(2);
  line->Draw();


}
