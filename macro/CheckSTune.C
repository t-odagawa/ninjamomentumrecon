void CheckSTune() {

  //string inputdir = "/home/t2k/odagawa/data/mc_data/particlegun/pbeta_recon/test20220523";
  string inputdir = "/home/t2k/odagawa/data/mc_data/particlegun/pbeta_recon/muon/true_level";

  const double muon_mass = 105.658;

  vector<double> bias;
  vector<double> pbeta;
  vector<double> bias_err;
  vector<double> pbeta_err;

  // TFile *ofile = new TFile("~/muon_wo_s_tune_graph.root", "recreate");  
  TFile *ofile = new TFile("~/muon_w_s_tune_graph.root", "recreate");  

  for ( int i = 6; i < 41; i++ ) {
    int momentum = i * 50;
    double energy = sqrt(momentum * momentum + muon_mass * muon_mass);
    double ipbeta = momentum * momentum / energy;
    stringstream inputfile_ss;
    //inputfile_ss << inputdir << "/test_" << momentum << "_0_0_new_rad_wo_s_tune.root";
    inputfile_ss << inputdir << "/test_" << momentum << "_0_0_new_rad_lat_iron.root";
    TFile *inputfile = new TFile(inputfile_ss.str().c_str(), "read");
    TTree *inputtree = (TTree*)inputfile->Get("tree");
    
    TH1D *hist = new TH1D(Form("hist_%d", momentum), "", 200, -1, 1);
    
    inputtree->Draw(Form("(1/recon_pbeta_candidate[0][0] - 1/true_pbeta) * true_pbeta >> hist_%d", momentum));

    TF1 *gaus = new TF1("gaus", "gaus", -1, 1);

    hist->Fit(gaus, "", "", -1, 1);

    ofile->cd();
    hist->Write();

    bias.push_back(gaus->GetParameter(1));
    pbeta.push_back(ipbeta);
    bias_err.push_back(gaus->GetParError(1));
    pbeta_err.push_back(0.);

  }

  TGraphErrors *ge = new TGraphErrors(bias.size(),
				      &pbeta[0], &bias[0],
				      &pbeta_err[0], &bias_err[0]);
  ge->SetName("ge");


  ofile->cd();
  ge->Write();
  ofile->Close();  

}
