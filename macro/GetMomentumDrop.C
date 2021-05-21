void GetMomentumDrop() {

  gErrorIgnoreLevel = kWarning;
  gStyle->SetOptFit(0);
  gStyle->SetOptStat(0);
  
  TString filename = "momentum_drop.root";
  TFile *file = new TFile(filename, "open");
  TTree *tree = (TTree*)file->Get("tree");

  Double_t momentum_drop, momentum_drop_err;
  Double_t momentum, sideview, topview;
  tree->SetBranchAddress("momentum_drop", &momentum_drop);
  tree->SetBranchAddress("momentum_drop_err", &momentum_drop_err);
  tree->SetBranchAddress("momentum", &momentum);
  tree->SetBranchAddress("sideview", &sideview);
  tree->SetBranchAddress("topview", &topview);

  TCanvas *canvas = new TCanvas("canvas", "canvas");
  TString pdfname = "momentum_drop_function.pdf";
  canvas->SaveAs(pdfname + "[", "pdf");

  const Int_t nmomentum = 13;
  const Int_t nsideview = 6;
  const Int_t ntopview  = 6;
  
  Double_t x[nmomentum][nsideview][ntopview];
  Double_t y[nmomentum][nsideview][ntopview];
  Double_t xe[nmomentum][nsideview][ntopview];
  Double_t ye[nmomentum][nsideview][ntopview];
  
  Int_t imomentum[nsideview][ntopview] = {{}};
  
  for (Int_t i = 0; i < tree->GetEntries(); i++) {
    tree->GetEntry(i);
    Int_t isideview = (Int_t)sideview/15;
    Int_t itopview  = (Int_t)topview/15;

    x[imomentum[isideview][itopview]][isideview][itopview] = momentum;
    y[imomentum[isideview][itopview]][isideview][itopview] = momentum_drop;
    xe[imomentum[isideview][itopview]][isideview][itopview] = 0;
    ye[imomentum[isideview][itopview]][isideview][itopview] = momentum_drop_err;
    imomentum[isideview][itopview]++;
    cout <<"Momentum : " << momentum << ", "
	 << "Sideview : " << isideview * 15 << ",  "
	 << "Topview : " << itopview * 15 << ", "
	 << "imomentum :" << imomentum[isideview][itopview] << endl; 
  }

  // Parameters dependent to the initial angle
  Double_t a[nsideview][ntopview], b[nsideview][ntopview], c[nsideview][ntopview];
  Double_t ae[nsideview][ntopview], be[nsideview][ntopview], ce[nsideview][ntopview];
  
  for (Int_t isideview = 0; isideview < nsideview; isideview++) {
    for (Int_t itopview = 0; itopview < ntopview; itopview++) {

      // Make graph with error bars
      Double_t gx[nmomentum], gy[nmomentum], gxe[nmomentum], gye[nmomentum];
      Int_t igmomentum = 0;

      while (igmomentum < nmomentum) {
	gx[igmomentum] = x[igmomentum][isideview][itopview];
	gy[igmomentum] = y[igmomentum][isideview][itopview];
	gxe[igmomentum] = xe[igmomentum][isideview][itopview];
	gye[igmomentum] = ye[igmomentum][isideview][itopview];
	cout << gx[igmomentum] << " " << gy[igmomentum] << " " << gxe[igmomentum] << " " << gye[igmomentum] << endl;
	igmomentum++;
      }
      

      TGraphErrors *ge = new TGraphErrors(nmomentum, gx, gy, gxe, gye);
      ge->SetTitle(Form("Momentum drop (#theta_{x} = %d deg., #theta_{y} = %d deg.);Initial momentum [MeV/c];Momentum drop (dp/dpl)",
			isideview * 15, itopview * 15));
      ge->GetYaxis()->SetTitleOffset(1.1);
      ge->SetMarkerStyle(20);
      canvas->cd();
      ge->Draw("AP");

      // Fitting
      TF1 *f = new TF1("f", "[0] / (x - [1]) + [2]");
      //TF1 *f = new TF1("f", "-[0] * log(x-[1]) + [2] * x + [3]");
      //TF1 *f = new TF1("f", "[0] * exp(-[1] * x) + [3]");
      ge->Fit(f, "", "", 200, 1500);

      a[isideview][itopview] = f->GetParameter(0);
      b[isideview][itopview] = f->GetParameter(1);
      c[isideview][itopview] = f->GetParameter(2);
      ae[isideview][itopview] = f->GetParError(0);
      be[isideview][itopview] = f->GetParError(1);
      ce[isideview][itopview] = f->GetParError(2);
      
      canvas->SaveAs(pdfname, "pdf");
  
    }
  }

  // Only angles less than 45 degree is used for now
  vector<Double_t> va, vb, vc;
  vector<Double_t> vae, vbe, vce;
  vector<Double_t> vangle, vanglee;
  for (Int_t isideview = 0; isideview < 4; isideview++) {
    for (Int_t itopview = 0; itopview < 4; itopview++) {
      Double_t angle = TMath::ATan(TMath::Hypot(TMath::Tan(isideview * 15 * TMath::DegToRad()), TMath::Tan(itopview * 15 * TMath::DegToRad()))) * TMath::RadToDeg();
      if (ae[isideview][itopview] > 100) continue; // ignore unsuccessful fitting
      cout << "Angle : " << angle << endl;
      va.push_back(a[isideview][itopview]);
      vae.push_back(ae[isideview][itopview]);
      vb.push_back(b[isideview][itopview]);
      vbe.push_back(be[isideview][itopview]);
      vc.push_back(c[isideview][itopview]);
      vce.push_back(ce[isideview][itopview]);
      vangle.push_back(angle);
      vanglee.push_back(0);
    }
  }

  canvas->SaveAs(pdfname + "]", "pdf");
  
  TString parpdfname = "momentum_drop_function_parameter.pdf";
  canvas->SaveAs(parpdfname + "[", "pdf");

  TF1 *constant = new TF1("constant","[0]");
  TF1 *fc = new TF1("fc", "[0] * sqrt(1 + tan(0.01745329251 * x) * tan(0.01745329251 * x))");
  fc->SetParameter(0, 0.6);
  
  TGraphErrors *agraph = new TGraphErrors(va.size(), &vangle[0], &va[0], &vanglee[0], &vae[0]);
  agraph->SetMarkerStyle(20);
  agraph->SetTitle(";Angle [deg];a");
  agraph->Draw("AP");
  agraph->Fit(constant);
  std::cout << "A(theta) = "
	    << constant->GetParameter(0) << " +/- "
	    << constant->GetParError(0) << std::endl;
  canvas->SaveAs(parpdfname, "pdf");

  TGraphErrors *bgraph = new TGraphErrors(vb.size(), &vangle[0], &vb[0], &vanglee[0], &vbe[0]);
  bgraph->SetMarkerStyle(20);
  bgraph->SetTitle(";Angle [deg];b");
  bgraph->Draw("AP");
  bgraph->Fit(constant);
  std::cout << "B(theta) = "
	    << constant->GetParameter(0) << " +/- "
	    << constant->GetParError(0) << std::endl;
  canvas->SaveAs(parpdfname, "pdf");

  TGraphErrors *cgraph = new TGraphErrors(vc.size(), &vangle[0], &vc[0], &vanglee[0], &vce[0]);
  cgraph->SetMarkerStyle(20);
  cgraph->SetTitle(";Angle [deg];c");
  cgraph->Draw("AP");
  cgraph->Fit(fc);
  canvas->SaveAs(parpdfname, "pdf");

  canvas->SaveAs(parpdfname + "]", "pdf");

}
