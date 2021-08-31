//#define STRAIGHT 1
#define ANGLEDEF 1
#define COLHIST 1

void CheckSmearEffect() {

  gStyle->SetOptStat(0);

  TString inputdirname = "/home/t2k/odagawa/data/mc_data/particlegun/pbeta_recon/test20210818/smear_all";
  TString outputdirname = "/home/t2k/odagawa/data/plots/mcs_test/20210818smear";
#ifdef STRAIGHT
  TString outputpdfname = "/check_smear_effect_mom.pdf";
#elif ANGLEDEF
  //TString outputpdfname = "/check_smear_effect_angle.pdf";
  TString outputpdfname = "/check_smear_effect_angle_1000mev.pdf";
#else
  TString outputpdfname = "/check_smear_effect.pdf";
#endif
  outputpdfname = outputdirname + outputpdfname;


  TCanvas *c = new TCanvas("c", "c");
  c->Print(outputpdfname + "[", "pdf");

  TF1 *gaus = new TF1("gaus", "gaus", -1, 1);
  Double_t bias[40][15][15] = {};
  Double_t sigma[40][15][15] = {};
  Double_t bias_err[40][15][15] = {};
  Double_t sigma_err[40][15][15] = {};

#ifdef ANGLEDEF
  for (Int_t imomentum = 20; imomentum < 21; imomentum++) { // only one momentum
#else
  for (Int_t imomentum = 1; imomentum < 41; imomentum++) {
#endif
    Double_t momentum = imomentum * 50.;

#ifdef COLHIST
#ifdef ANGLEDEF
    TH2D *hist_bias = new TH2D("hist_bias", "", 4, -2.5, 17.5, 10, -2.5, 47.5);
    TH2D *hist_sigma = new TH2D("hist_sigma", "", 4, -2.5, 17.5, 10, -2.5, 47.5);
#else
    TH2D *hist_bias = new TH2D("hist_bias", "", 16, -2.5, 77.5, 16, -2.5, 77.5);
    TH2D *hist_sigma = new TH2D("hist_sigma", "", 16, -2.5, 77.5, 16, -2.5, 77.5);
#endif
    hist_bias->SetTitle(Form("Bias (%d MeV/c);#theta_{x} [deg];#theta_{y} [deg]", (Int_t)momentum));
    hist_sigma->SetTitle(Form("Sigma (%d MeV/c);#theta_{x} [deg];#theta_{y} [deg]", (Int_t)momentum));
#endif

#ifdef STRAIGHT
    for (Int_t isideview = 0; isideview < 1; isideview++) { // only straight track
#elif ANGLEDEF
    for (Int_t isideview = 0; isideview < 4; isideview++) { // use max # of angles
#else
    for (Int_t isideview = 0; isideview < 16; isideview++) {
#endif
      Int_t sideview = isideview * 5;

#ifdef STRAIGHT
      for (Int_t itopview = 0; itopview < 1; itopview++) { // only straight track
#elif ANGLEDEF
      for (Int_t itopview = 0; itopview < 10; itopview++) { // use max # of angles
#else
      for (Int_t itopview = 0; itopview < 16; itopview++) {
#endif
      	Int_t topview = itopview * 5;
	cout << "Momentum : " << momentum << endl;
	cout << "Sideview : " << sideview << endl;
	cout << "Topview : " << topview << endl;
	TString inputfilename = Form("/test_%d_%d_%d_new_rad_lat_iron_smear.root",
				     (Int_t)momentum, sideview, topview);
	inputfilename = inputdirname + inputfilename;
	TFile *inputfile = new TFile(inputfilename, "read");
	TTree *inputtree = (TTree*)inputfile->Get("tree");
	Double_t true_pbeta;
	Double_t recon_pbeta;
	inputtree->SetBranchAddress("true_pbeta", &true_pbeta);
	inputtree->SetBranchAddress("recon_pbeta", &recon_pbeta);

	TH1D *hist_residual = new TH1D("hist_residual",
				       Form("1/p#beta residual (%d MeV/c, %d, %d)",
					    (Int_t)momentum, sideview, topview),
				       100 ,-1, 1);
	for (Int_t ientry = 0; ientry < inputtree->GetEntries(); ientry++) {
	  inputtree->GetEntry(ientry);
	  if (recon_pbeta >= 3000) continue; // reconstruction saturated
	  hist_residual->Fill((1/recon_pbeta - 1/true_pbeta) * true_pbeta);
	}

	c->cd();
	hist_residual->Draw("hist");

	hist_residual->Fit(gaus, "", "", -1, 1);
	bias[imomentum-1][isideview][itopview] = gaus->GetParameter(1);
	sigma[imomentum-1][isideview][itopview] = gaus->GetParameter(2);
	bias_err[imomentum-1][isideview][itopview] = gaus->GetParError(1);
	sigma_err[imomentum-1][isideview][itopview] = gaus->GetParError(2);
	// cout << "Bias : " << bias[imomentum-1][isideview][itopview] << endl;
	// cout << "Sigma : " << sigma[imomentum-1][isideview][itopview] << endl;
	gaus->Draw("same");
	c->Update();
	c->Print(outputpdfname,"pdf");

#ifdef COLHIST
	hist_bias->Fill(sideview, topview, bias[imomentum-1][isideview][itopview]);
	hist_sigma->Fill(sideview, topview, sigma[imomentum-1][isideview][itopview]);
#endif

      }
    }

#ifdef COLHIST
    c->cd();
    hist_bias->Draw("colz text");
    c->Update();
    c->Print(outputpdfname, "pdf");
    hist_sigma->Draw("colz text");
    c->Update();
    c->Print(outputpdfname, "pdf");
#endif

    delete hist_bias;
    delete hist_sigma;

  }

  std::vector<Double_t> x, y, xe, ye;
  TLine *line = new TLine();
  line->SetLineStyle(2);
  line->SetLineColor(kRed);

#ifdef STRAIGHT
  // straight penetrating muon momentum dependence graph
  for (Int_t i = 1; i < 40; i++) {
    x.push_back((i+1) * 50.);
    y.push_back(bias[i][0][0]);
    xe.push_back(0.);
    ye.push_back(sigma[i][0][0]);
  }

  TGraphErrors *g_0_0 = new TGraphErrors(39,
					 &x[0],  &y[0],
					 &xe[0], &ye[0]);
  g_0_0->SetMarkerStyle(21);
  g_0_0->SetTitle("ECC penetrating muon (0,0);Muon momentum [MeV/c];Bias and #sigma");
  c->cd();
  g_0_0->Draw("AP");
  gPad->Update(); // Call this function before GetUymin()/GetUymax()
  line->DrawLine(gPad->GetUxmin(), 0.,
		 gPad->GetUxmax(), 0.);
  c->Print(outputpdfname, "pdf");
  x.clear(); x.shrink_to_fit();
  y.clear(); y.shrink_to_fit();
  xe.clear(); xe.shrink_to_fit();
  ye.clear(); ye.shrink_to_fit();

  // straight penetrating muon momentum dependence graph bias only
  for (Int_t i = 1; i < 40; i++) {
    x.push_back((i+1) * 50.);
    y.push_back(bias[i][0][0]);
    xe.push_back(0.);
    ye.push_back(bias_err[i][0][0]);
  }

  TGraphErrors *g_0_0_bias = new TGraphErrors(39,
					      &x[0],  &y[0],
					      &xe[0], &ye[0]);
  g_0_0_bias->SetMarkerStyle(21);
  g_0_0_bias->SetTitle("ECC penetrating muon (0,0) Bias;Muon momentum [MeV/c];Bias");
  c->cd();
  g_0_0_bias->Draw("AP");
  gPad->Update(); // Call this function before GetUymin()/GetUymax()
  line->DrawLine(gPad->GetUxmin(), 0.,
		 gPad->GetUxmax(), 0.);
  c->Print(outputpdfname, "pdf");

  x.clear(); x.shrink_to_fit();
  y.clear(); y.shrink_to_fit();
  xe.clear(); xe.shrink_to_fit();
  ye.clear(); ye.shrink_to_fit();

  // straight penetrating muon momentum dependence graph sigma only

  for (Int_t i = 1; i < 40; i++) {
    x.push_back((i+1) * 50.);
    y.push_back(sigma[i][0][0]);
    xe.push_back(0.);
    ye.push_back(sigma_err[i][0][0]);
  }

  TGraphErrors *g_0_0_sigma = new TGraphErrors(39,
					       &x[0],  &y[0],
					       &xe[0], &ye[0]);
  g_0_0_sigma->SetMarkerStyle(21);
  g_0_0_sigma->SetTitle("ECC penetrating muon (0,0) Sigma;Muon momentum [MeV/c];#sigma");
  c->cd();
  g_0_0_sigma->Draw("AP");
  c->Print(outputpdfname, "pdf");

#endif

#ifdef ANGLEDEF
  // 500 MeV/c muon angle dependence graph 9
  // 1 GeV/c muon angle dependence graph 19
  const Int_t momentum_index = 19;
  Int_t momentum = 50 * (momentum_index + 1);
  for (Int_t i = 0; i < 4; i++) {
    Double_t sideview = i * 5. * TMath::DegToRad(); // rad
    for (Int_t j = 0; j < 10; j++) {
      Double_t topview = j * 5. * TMath::DegToRad(); // rad
      Double_t tangent = TMath::Hypot(TMath::Tan(sideview), TMath::Tan(topview));
      Double_t unit_length = TMath::Sqrt(1 + tangent * tangent);
      x.push_back(unit_length);
      y.push_back(bias[momentum_index][i][j]);
      xe.push_back(0.);
      ye.push_back(sigma[momentum_index][i][j]);
    }
  }

  TGraphErrors *g_500mev = new TGraphErrors(40,
					    &x[0],  &y[0],
					    &xe[0], &ye[0]);
  g_500mev->SetMarkerStyle(21);
  g_500mev->SetTitle(Form("ECC penetrating muon (%d MeV/c);Initial angle (#sqrt{1+tan^{2}#theta});Bias and #sigma", momentum));
  c->cd();
  g_500mev->Draw("AP");
  gPad->Update(); // Call this function before GetUymin()/GetUymax()
  line->DrawLine(gPad->GetUxmin(), 0.,
		 gPad->GetUxmax(), 0.);
  c->Print(outputpdfname, "pdf");
  x.clear(); x.shrink_to_fit();
  y.clear(); y.shrink_to_fit();
  xe.clear(); xe.shrink_to_fit();
  ye.clear(); ye.shrink_to_fit();

  // 500 MeV/c muon angle dependence graph bias only

  for (Int_t i = 0; i < 4; i++) {
    Double_t sideview = i * 5. * TMath::DegToRad(); // rad
    for (Int_t j = 0; j < 10; j++) {
      Double_t topview = j * 5. * TMath::DegToRad(); // rad
      Double_t tangent = TMath::Hypot(TMath::Tan(sideview), TMath::Tan(topview));
      Double_t unit_length = TMath::Sqrt(1 + tangent * tangent);
      x.push_back(unit_length);
      y.push_back(bias[momentum_index][i][j]);
      xe.push_back(0.);
      ye.push_back(bias_err[momentum_index][i][j]);
    }
  }

  TGraphErrors *g_500mev_bias = new TGraphErrors(40,
						 &x[0],  &y[0],
						 &xe[0], &ye[0]);
  g_500mev_bias->SetMarkerStyle(21);
  g_500mev_bias->SetTitle(Form("ECC penetrating muon (%d MeV/c) Bias;Initial angle (#sqrt{1+tan^{2}#theta});Bias", momentum));
  c->cd();
  g_500mev_bias->Draw("AP");
  gPad->Update(); // Call this function before GetUymin()/GetUymax()
  line->DrawLine(gPad->GetUxmin(), 0.,
		 gPad->GetUxmax(), 0.);
  c->Print(outputpdfname, "pdf");

  x.clear(); x.shrink_to_fit();
  y.clear(); y.shrink_to_fit();
  xe.clear(); xe.shrink_to_fit();
  ye.clear(); ye.shrink_to_fit();

  // 500 MeV/c muon angle dependence graph sigma only

  for (Int_t i = 0; i < 4; i++) {
    Double_t sideview = i * 5. * TMath::DegToRad(); // rad
    for (Int_t j = 0; j < 10; j++) {
      Double_t topview = j * 5. * TMath::DegToRad(); // rad
      Double_t tangent = TMath::Hypot(TMath::Tan(sideview), TMath::Tan(topview));
      Double_t unit_length = TMath::Sqrt(1 + tangent * tangent);
      x.push_back(unit_length);
      y.push_back(sigma[momentum_index][i][j]);
      xe.push_back(0.);
      ye.push_back(sigma_err[momentum_index][i][j]);
    }
  }

  TGraphErrors *g_500mev_sigma = new TGraphErrors(40,
						  &x[0],  &y[0],
						  &xe[0], &ye[0]);
  g_500mev_sigma->SetMarkerStyle(21);
  g_500mev_sigma->SetTitle(Form("ECC penetrating muon (%d MeV/c) Sigma;Initial angle (#sqrt{1+tan^{2}#theta});#sigma", momentum));
  c->cd();
  g_500mev_sigma->Draw("AP");
  c->Print(outputpdfname, "pdf");

#endif

  c->Print(outputpdfname + "]", "pdf");

}

