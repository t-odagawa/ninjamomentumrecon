void CheckHighland() {

  TString filename = "/home/t2k/odagawa/NinjaMomentumRecon/build/angle_difference_fit.root";
  TFile *file = new TFile(filename, "read");
  TTree *tree = (TTree*)file->Get("tree");

  Double_t beta, pbeta, momentum;
  Double_t sigma_lat, sigma_lat_err;
  Double_t sideview, topview;

  tree->SetBranchAddress("beta", &beta);
  tree->SetBranchAddress("pbeta", &pbeta);
  tree->SetBranchAddress("momentum", &momentum);

  tree->SetBranchAddress("sideview", &sideview);
  tree->SetBranchAddress("topview", &topview);

  tree->SetBranchAddress("sigma_lat", &sigma_lat);
  tree->SetBranchAddress("sigma_lat_err", &sigma_lat_err);

  TString pdfname = "~/check_highland.pdf";
  TCanvas *c = new TCanvas("c", "c");
  c->Print(pdfname + "[", "pdf");

  const Int_t nentry = tree->GetEntries();

  Double_t pbeta_sigma_prod[29][10][10] = {};
  Double_t pbeta_sigma_prod_err[29][10][10] = {};
  Double_t pbeta_array[29][10][10] = {};
  
  for (Int_t ientry = 0; ientry < nentry; ientry++) {

    tree->GetEntry(ientry);

    if (momentum == 50) continue;
    Int_t imomentum = momentum / 50 - 2;
    Int_t isideview = sideview / 5;
    Int_t itopview = topview / 5;

    pbeta_sigma_prod[imomentum][isideview][itopview] = pbeta * sigma_lat;
    pbeta_sigma_prod_err[imomentum][isideview][itopview] = pbeta * sigma_lat_err;
    pbeta_array[imomentum][isideview][itopview] = pbeta;
    /*
    std::cout << "pbeta_sigma_prod : " << pbeta_sigma_prod[imomentum][isideview][itopview] << ", "
	      << "pbeta_sigma_prod_err : " << pbeta_sigma_prod_err[imomentum][isideview][itopview] << ", "
	      << "pbeta_array : " << pbeta_array[imomentum][isideview][itopview] << std::endl;
    */
  }

  std::vector<Double_t> unit_path_length_vec = {};
  std::vector<Double_t> pbeta_sigma_ave_vec = {};
  std::vector<Double_t> unit_path_length_err_vec = {};
  std::vector<Double_t> pbeta_sigma_ave_err_vec = {};

  for (Int_t isideview = 0; isideview < 10; isideview++) {
    for (Int_t itopview = 0; itopview < 10; itopview++) {

      Double_t tangent = TMath::Hypot(TMath::Tan(isideview * 5. * TMath::DegToRad()),
				      TMath::Tan(itopview  * 5. * TMath::DegToRad()));
      Double_t unit_path_length = TMath::Hypot(tangent, 1.);
      unit_path_length_vec.push_back(unit_path_length);
      unit_path_length_err_vec.push_back(0.);

      Double_t x[29], y[29], xe[29], ye[29];

      for (Int_t imomentum = 0; imomentum < 29; imomentum++) {
	x[imomentum] = pbeta_array[imomentum][isideview][itopview];
	y[imomentum] = pbeta_sigma_prod[imomentum][isideview][itopview];
	xe[imomentum] = 0.;
	ye[imomentum] = pbeta_sigma_prod_err[imomentum][isideview][itopview];
      }

      TGraphErrors *ge = new TGraphErrors(29, x, y, xe, ye);
      ge->SetMarkerStyle(10);
      ge->SetTitle(Form("#sigma_{lat} x p#betac (#theta_{x} = %d, #theta_{y} = %d);p#beta [MeV/c];#sigma_{lat} x p#betac [MeV]",
			isideview * 5, itopview * 5));
      TF1 *constant = new TF1("constant", "[0]");
      c->cd();
      ge->Draw("AP");
      ge->Fit(constant, "", "", 100, 1500);
      c->Print(pdfname, "pdf");
      
      pbeta_sigma_ave_vec.push_back(constant->GetParameter(0));
      pbeta_sigma_ave_err_vec.push_back(constant->GetParError(0));

    }
  }

  TGraphErrors *highland_check_ge = new TGraphErrors(unit_path_length_vec.size(),
						     &unit_path_length_vec[0],
						     &pbeta_sigma_ave_vec[0],
						     &unit_path_length_err_vec[0],
						     &pbeta_sigma_ave_err_vec[0]);
  highland_check_ge->SetMarkerStyle(10);
  //gPad->SetLogy();
  //gPad->SetLogx();
  highland_check_ge->SetTitle(";#sqrt{1 + tan^{2}#theta};#sigma_{lat} x p#betac [MeV]");
  //TF1 *highland_func = new TF1("highland_func", "[0] * sqrt([1] * x) * (1 + 0.038 * log([1] * x))", 1.1, 1.6);
  TF1 *highland_func = new TF1("highland_func", "13.6 * pow([0] * x, [1])", 1., 1.8);
  highland_func->SetParameter(0, 0.03);
  highland_func->SetParameter(1, 1.);
  c->cd();
  highland_check_ge->Draw("AP");
  c->Print(pdfname, "pdf");
  highland_check_ge->Fit(highland_func, "", "", 1.1, 1.7);
  highland_func->Draw("same");
  c->Print(pdfname, "pdf");

  c->Print(pdfname + "]", "pdf");

}
