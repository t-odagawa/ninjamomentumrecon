void CheckAlphaDifference() {

  gStyle->SetOptStat(0);

  TString filename1 = "/home/t2k/odagawa/data/plots/mcs_test/test_500_0_0_new_rad_1.0.root";
  TFile *file1 = new TFile(filename1, "read");
  TTree *tree1 = (TTree*)file1->Get("tree");
  Double_t recon_pbeta_1;
  Double_t true_pbeta_1;
  tree1->SetBranchAddress("recon_pbeta", &recon_pbeta_1);
  tree1->SetBranchAddress("true_pbeta", &true_pbeta_1);

  TString filename2 = "/home/t2k/odagawa/data/plots/mcs_test/test_500_0_0_new_rad_1.1.root";
  TFile *file2 = new TFile(filename2, "read");
  TTree *tree2 = (TTree*)file2->Get("tree");
  Double_t recon_pbeta_2;
  tree2->SetBranchAddress("recon_pbeta", &recon_pbeta_2);

  TString filename3 = "/home/t2k/odagawa/data/plots/mcs_test/test_500_0_0_new_rad_0.9.root";
  TFile *file3 = new TFile(filename3, "read");
  TTree *tree3 = (TTree*)file3->Get("tree");
  Double_t recon_pbeta_3;
  tree3->SetBranchAddress("recon_pbeta", &recon_pbeta_3);

  const Int_t nentry = tree1->GetEntries();

  const Double_t max_pbeta = 700.;

  TH2D *hist = new TH2D("hist", "Initial value difference;p#beta_{recon} (#alpha = 0);p#beta_{recon} (#alpha = 0.1)", 350, 0, max_pbeta, 350, 0, max_pbeta);
  //TH2D *hist = new TH2D("hist", "Initial value difference;p#beta_{recon} (#alpha = 0);p#beta_{recon} (#alpha = -0.1)", 350, 0, max_pbeta, 350, 0, max_pbeta);
  TH1D *hist_pbeta = new TH1D("hist_pbeta", "Reconstructed p#beta;(1/p#beta_{recon} - 1/p#beta_{true})p#beta_{true};Entries", 100, -1., 1.);
  //TH1D *hist_pbeta = new TH1D("hist_pbeta", "Reconstructed p#beta;p#beta_{recon};Entries", 100, 0., 1000.);

  TF1 *gaus = new TF1("gaus", "gaus", -1., 1.);

  for(Int_t ientry = 0; ientry < nentry; ientry++) {
    tree1->GetEntry(ientry);
    tree2->GetEntry(ientry);
    tree3->GetEntry(ientry);
    hist->Fill(recon_pbeta_1, recon_pbeta_2);
    //hist->Fill(recon_pbeta_1, recon_pbeta_3);

    if (fabs(recon_pbeta_1 / recon_pbeta_2 - 1.) < 1.e-3) {
    //if (fabs(recon_pbeta_1 / recon_pbeta_3 - 1.) < 1.e-3) {
      hist_pbeta->Fill((1/recon_pbeta_1 - 1/true_pbeta_1) * true_pbeta_1);
    }
  }

  TCanvas *c = new TCanvas("c", "c");
  c->cd();
  hist->Draw("colz");

  //c->SaveAs("~/mcs_alpha_compariton.pdf");

  TLine *line = new TLine();
  line->SetLineColor(kRed);
  line->SetLineWidth(2);

  // line->DrawLine(0, 0, max_pbeta, max_pbeta);
  c->Update();

  hist_pbeta->Draw("hist");
  hist_pbeta->Fit(gaus,"","",-0.5, 0.5);

}
