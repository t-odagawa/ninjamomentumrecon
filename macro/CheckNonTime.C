void CheckNonTime() {

  gStyle->SetOptStat(0);

  //TString filename = "/home/t2k/odagawa/data/momch_file/ecc5/non-time/mip_non_time_ecc5_recon.root";
  TString filename = "/home/t2k/odagawa/data/momch_file/ecc5/non-time/hip_non_time_ecc5_recon.root";
  TFile *file = new TFile(filename, "read");
  TTree *tree = (TTree*)file->Get("tree");

  Int_t entry_in_daily_file;
  Int_t npl;
  std::vector<Double_t> *ax = 0;
  std::vector<Double_t> *ay = 0;
  std::vector<Double_t> *vph = 0;
  std::vector<Double_t> *pixel_count = 0;
  Double_t recon_pbeta;
  Double_t recon_pbeta_err;
  Double_t log_likelihood[3][2];
  Double_t recon_pbeta_candidate[3][2];
  Double_t recon_pbeta_err_candidate[3][2];

  tree->SetBranchAddress("entry_in_daily_file", &entry_in_daily_file);
  tree->SetBranchAddress("npl", &npl);
  tree->SetBranchAddress("ax", &ax);
  tree->SetBranchAddress("ay", &ay);
  tree->SetBranchAddress("vph", &vph);
  tree->SetBranchAddress("pixel_count", &pixel_count);
  tree->SetBranchAddress("recon_pbeta", &recon_pbeta);
  tree->SetBranchAddress("recon_pbeta_err", &recon_pbeta_err);
  tree->SetBranchAddress("log_likelihood", log_likelihood);
  tree->SetBranchAddress("recon_pbeta_candidate", recon_pbeta_candidate);
  tree->SetBranchAddress("recon_pbeta_err_candidate", recon_pbeta_err_candidate);

  TString pdfname = filename + ".pdf";
  TCanvas *c = new TCanvas("c", "c");
  //const Double_t vph_min = 0.;
  const Double_t vph_min = 100.;
  //const Double_t vph_max = 150;
  const Double_t vph_max = 800;
  TH1D *hist_recon_pbeta = new TH1D("hist_recon_pbeta","Reconstructed p#beta (all);p#beta [MeV/c]", 100, 0, 5000);
  TH2D *hist_recon_pbeta_vph = new TH2D("hist_recon_pbeta_vph", "All Events;p#beta [MeV/c];VPH (average)",
					100, 0, 2000, 50, vph_min, vph_max);
  TH2D *hist_recon_pbeta_pixel = new TH2D("hist_recon_pbeta_pixel","All Events;p#beta [MeV/c];Pixel count (average)",
					  100, 0, 2000, 100, 0, 400);
  TH2D *hist_recon_pbeta_pixel_correct = new TH2D("hist_recon_pbeta_pixel_correct","All Events;p#beta [MeV/c];Pixel count (average correction)",
						  100, 0, 2000, 100, 0, 200);
  TH1D *hist_log_likelihood_diff = new TH1D("hist_log_likelihood_diff","PID w/o VPH;NLL(MIP) - NLL(HIP)",100,-10,10);
  TH1D *hist_recon_pbeta_tangent[10];
  TH2D *hist_recon_pbeta_vph_tangent[10];
  TH2D *hist_recon_pbeta_pixel_tangent[10];
  TH2D *hist_recon_pbeta_pixel_correct_tangent[10];
  TH1D *hist_log_likelihood_diff_tangent[10];
  TH1D *hist_log_likelihood_diff_pbeta[10];
  for ( Int_t ihist = 0; ihist < 10; ihist++ ) {
    Double_t ang_min = ihist * 0.1;
    Double_t ang_max;
    if (ihist < 9) ang_max = (ihist + 1) * 0.1;
    else ang_max = 10;
    hist_recon_pbeta_tangent[ihist] = new TH1D(Form("hist_recon_pbeta%d", ihist),Form("Reconstructed p#beta (%.1f < tan#theta < %.1f);p#beta [MeV/c]",
										      ang_min, ang_max), 100, 0, 5000);
    hist_recon_pbeta_vph_tangent[ihist] = new TH2D(Form("hist_recon_pbeta_vph%d", ihist), Form("%.1f < tan#theta < %.1f;p#beta [MeV/c];VPH (average)",
											       ang_min, ang_max),
						   100, 0, 2000, 50, vph_min, vph_max);
    hist_recon_pbeta_pixel_tangent[ihist] = new TH2D(Form("hist_recon_pbeta_pixel%d", ihist), Form("%.1f < tan#theta < %.1f;p#beta [MeV/c];Pixel count (average)",
												   ang_min, ang_max),
						     100, 0, 2000, 100, 0, 400);
    hist_recon_pbeta_pixel_correct_tangent[ihist] = new TH2D(Form("hist_recon_pbeta_pixel_correct%d", ihist), Form("%.1f < tan#theta < %.1f;p#beta [MeV/c];Pixel count (average correction)",
														   ang_min, ang_max),
							     100, 0, 2000, 100, 0, 200);
    hist_log_likelihood_diff_tangent[ihist] = new TH1D(Form("hist_log_likelihood_diff%d", ihist), Form("%.1f < tan#theta < %.1f;NLL(MIP) - NLL(HIP)",
												       ang_min, ang_max),
						       100, -10, 10);
  }


  for ( Int_t ientry = 0; ientry < tree->GetEntries(); ientry++ ) {
    if (ientry % 1000 == 0)
      std::cout << ientry << "\r";
    tree->GetEntry(ientry);
    Double_t angle = std::sqrt(ax->at(0) * ax->at(0) + ay->at(0) * ay->at(0));
    Int_t binid = std::floor(angle * 10);
    if (binid > 9) binid = 9;
    std::vector<Double_t> pixel_count_correct;
    for ( Int_t j = 0; j < pixel_count->size(); j++ ) {
      if ( pixel_count->at(j) == -1 ) {
	pixel_count_correct.push_back(-1);
      }
      else
	pixel_count_correct.push_back(pixel_count->at(j) / std::sqrt(1 + ax->at(j) * ax->at(j) + ay->at(j) * ay->at(j)));
    }
    Double_t vph_average = std::accumulate(vph->begin(), vph->end(), 0);
    vph_average /= vph->size();
    pixel_count->erase(std::remove(pixel_count->begin(), pixel_count->end(), -1), pixel_count->end());
    Double_t pixel_count_average = std::accumulate(pixel_count->begin(), pixel_count->end(), 0);
    pixel_count_average /= pixel_count->size();
    pixel_count_correct.erase(std::remove(pixel_count_correct.begin(), pixel_count_correct.end(), -1),
			      pixel_count_correct.end());    
    Double_t pixel_count_correct_average = std::accumulate(pixel_count_correct.begin(),
							   pixel_count_correct.end(), 0);
    pixel_count_correct_average /= pixel_count_correct.size();
    hist_recon_pbeta->Fill(recon_pbeta);
    hist_recon_pbeta_vph->Fill(recon_pbeta, vph_average);
    hist_recon_pbeta_pixel->Fill(recon_pbeta, pixel_count_average);
    hist_recon_pbeta_pixel_correct->Fill(recon_pbeta, pixel_count_correct_average);
    hist_log_likelihood_diff->Fill(log_likelihood[0][0] - log_likelihood[2][0]);
    hist_recon_pbeta_tangent[binid]->Fill(recon_pbeta);
    hist_recon_pbeta_vph_tangent[binid]->Fill(recon_pbeta, vph_average);
    hist_recon_pbeta_pixel_tangent[binid]->Fill(recon_pbeta, pixel_count_average);
    hist_recon_pbeta_pixel_correct_tangent[binid]->Fill(recon_pbeta, pixel_count_correct_average);
    hist_log_likelihood_diff_tangent[binid]->Fill(log_likelihood[0][0] - log_likelihood[2][0]);
  }

  std::cout << std::endl;

  c->cd();
  c->Print(pdfname + "[", "pdf");
  hist_recon_pbeta->Draw("hist");
  c->Print(pdfname, "pdf");
  for ( Int_t ibin = 0; ibin < 10; ibin++ ) {
    hist_recon_pbeta_tangent[ibin]->Draw("hist");
    c->Print(pdfname, "pdf");
  }
  hist_recon_pbeta_vph->Draw("colz");
  c->Print(pdfname, "pdf");
  for ( Int_t ibin = 0; ibin < 10; ibin++ ) {
    hist_recon_pbeta_vph_tangent[ibin]->Draw("colz");
    c->Print(pdfname, "pdf");
  }
  hist_recon_pbeta_pixel->Draw("colz");
  c->Print(pdfname, "pdf");
  hist_recon_pbeta_pixel_correct->Draw("colz");
  c->Print(pdfname, "pdf");
  for ( Int_t ibin = 0; ibin < 10; ibin++ ) {
    hist_recon_pbeta_pixel_tangent[ibin]->Draw("colz");
    c->Print(pdfname, "pdf");
    hist_recon_pbeta_pixel_correct_tangent[ibin]->Draw("colz");
    c->Print(pdfname, "pdf");
  }
  hist_log_likelihood_diff->Draw("hist");
  c->Print(pdfname, "pdf");
  for ( Int_t ibin = 0; ibin < 10; ibin++ ) {
    hist_log_likelihood_diff_tangent[ibin]->Draw("hist");
    c->Print(pdfname, "pdf");
  }
  c->Print(pdfname + "]", "pdf");

}
