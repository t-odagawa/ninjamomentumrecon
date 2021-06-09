void CheckHighland() {

  TString filename = "/home/t2k/odagawa/NinjaMomentumRecon/build/angle_difference_fit.root";
  TFile *file = new TFile(filename, "read");
  TTree *tree = (TTree*)file->Get("tree");

  Double_t beta, pbeta, momentum;
  Double_t sigma_lat, sigma_lat_err;

  tree->SetBranchAddress("beta", &beta);
  tree->SetBranchAddress("pbeta", &pbeta);
  tree->SetBranchAddress("momentum", &momentum);

  tree->SetBranchAddress("sigma_lat", &sigma_lat);
  tree->SetBranchAddress("sigma_lat_err", &sigma_lat_err);

  TString pdfname = "./check_highland.pdf";
  TCanvas *c = new TCanvas("c", "c");
  c->Print(pdfname + "[", "pdf");

  const Int_t nentry = tree->GetEntries();
  std::vector<Double_t> pbeta_sigma_prod = {};
  std::vector<Double_t> pbeta_vec = {};
  std::vector<Double_t> pbeta_sigma_prod_err = {};
  std::vector<Double_t> pbeta_err = {};

  for (Int_t ientry = 0; ientry < nentry; ientry++) {
    tree->GetEntry(ientry);
    pbeta_vec.push_back(pbeta);
    pbeta_err.push_back(0.);
    pbeta_sigma_prod.push_back(pbeta * sigma_lat);
    pbeta_sigma_prod_err.push_back(pbeta * sigma_lat_err);
  }

  TGraphErrors *ge = new TGraphErrors(pbeta_sigma_prod.size(),
				      &pbeta_vec[0], &pbeta_sigma_prod[0],
				      &pbeta_err[0], &pbeta_sigma_prod_err[0]);
  ge->SetMarkerStyle(10);
  ge->SetTitle("#sigma_{lat} x p#betac;p#beta [MeV/c];#sigma_{lat} x p#betac [MeV]");

  TF1 *constant = new TF1("constant", "[0]");

  c->cd();
  ge->Draw("AP");
  ge->Fit(constant, "", "", 100, 1500);
  c->Print(pdfname, "pdf");

  c->Print(pdfname + "]", "pdf");

}
