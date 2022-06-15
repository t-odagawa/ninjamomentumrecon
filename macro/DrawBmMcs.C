double GetBmErr(double range) {
  if ( range < 300 ) return 0.079020061;
  else if ( range < 400 ) return 0.079020061;
  else if ( range < 500 ) return 0.080035672;
  else if ( range < 600 ) return 0.086198609;
  else if ( range < 700 ) return 0.082791300;
  else if ( range < 800 ) return 0.071492218;
  else if ( range < 900 ) return 0.074543873;
  else if ( range < 1000 ) return 0.072054229;
  else if ( range < 1100 ) return 0.056333879;
  else if ( range < 1200 ) return 0.048178771;
  else if ( range < 1300 ) return 0.037908676;
  else if ( range < 1400 ) return 0.030906878;
  else return 0.036638113;
}

void DrawBmMcs() {

  gRandom->SetSeed(0);

  const double muon_mass = 105.658;

  TString filename = "/home/t2k/odagawa/data/mcs-recon/ecc5/bm-compare/mcs_recon_bm_compare.root";
  TFile *file = new TFile(filename, "read");
  TTree *tree = (TTree*)file->Get("tree");

  Double_t mcs_pbeta, range_pbeta;
  Double_t range_momentum;
  Double_t mcs_pbeta_err;

  tree->SetBranchAddress("mcs_pbeta", &mcs_pbeta);
  tree->SetBranchAddress("mcs_pbeta_err", &mcs_pbeta_err);
  tree->SetBranchAddress("range_pbeta", &range_pbeta);
  tree->SetBranchAddress("range_momentum", &range_momentum);

  TH2D *hist = new TH2D("hist",";Range-base p#beta [MeV/c];MCS p#beta [MeV/c]",
			50, 0, 1500, 50, 0, 1500);

  std::vector<double> mcs_ge = {};
  std::vector<double> range_ge = {};
  std::vector<double> mcs_err_min_ge  = {};
  std::vector<double> mcs_err_max_ge  = {};
  std::vector<double> range_err_min_ge = {};
  std::vector<double> range_err_max_ge = {};

  for ( int i = 0; i < tree->GetEntries(); i++ ) {
    
    tree->GetEntry(i);

    hist->Fill(range_pbeta, mcs_pbeta);
    
    if ( gRandom->Uniform() > 0.9 ) {
      range_ge.push_back(range_pbeta);
      mcs_ge.push_back(mcs_pbeta);
      double range_err_frac = GetBmErr(range_momentum);
      double range_mom_min = range_momentum * ( 1 - range_err_frac);
      double range_mom_max = range_momentum * ( 1 + range_err_frac);
      double range_pbeta_min = range_mom_min * range_mom_min / sqrt(range_mom_min * range_mom_min + muon_mass * muon_mass);
      double range_pbeta_max = range_mom_max * range_mom_max / sqrt(range_mom_max * range_mom_max + muon_mass * muon_mass);
      range_err_min_ge.push_back(range_pbeta - range_pbeta_min);
      range_err_max_ge.push_back(range_pbeta_max - range_pbeta);
      mcs_err_min_ge.push_back(mcs_pbeta_err);
      mcs_err_max_ge.push_back(mcs_pbeta_err);
    }
    
  }

  TGraphAsymmErrors *ge = new TGraphAsymmErrors(range_ge.size(),
						&range_ge[0], &mcs_ge[0],
						&range_err_min_ge[0], &range_err_max_ge[0], 
						&mcs_err_min_ge[0], &mcs_err_max_ge[0]);

  ge->SetName("ge");
  ge->SetMarkerStyle(4);
  ge->SetLineWidth(2);

  TCanvas *c1 = new TCanvas("c1", "c1");
  gStyle->SetOptStat(0);
  hist->Draw("colz");
  hist->GetYaxis()->CenterTitle();
  hist->GetXaxis()->CenterTitle();
  ge->Draw("SAME P");

  TLine *l = new TLine();
  l->SetLineColor(kRed);
  l->SetLineWidth(1);
  l->SetLineStyle(kDashed);
  l->DrawLine(0,0,1500,1500);

  gPad->SetGrid();
   

}
