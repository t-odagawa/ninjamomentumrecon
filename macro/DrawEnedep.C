void DrawEnedep() {

  TString filename_wo = "/home/t2k/odagawa/muon_wo_enedep_graph.root";
  TString filename_w = "/home/t2k/odagawa/muon_w_angres_mom_graph.root";

  TFile *file_wo = new TFile(filename_wo, "read");
  TFile *file_w = new TFile(filename_w, "read");

  TLegend *l = new TLegend(0.25, 0.75, 0.7, 0.89);

  file_w->cd();
  TGraphErrors *ge_bias_w = (TGraphErrors*)file_w->Get("ge_bias");
  ge_bias_w->SetName("ge_bias_w");
  TGraphErrors *ge_width_w = (TGraphErrors*)file_w->Get("ge_width");
  ge_width_w->SetName("ge_width_w");

  file_wo->cd();
  TGraphErrors *ge_bias_wo = (TGraphErrors*)file_wo->Get("ge_bias");
  ge_bias_wo->SetName("ge_bias_wo");
  TGraphErrors *ge_width_wo = (TGraphErrors*)file_wo->Get("ge_width");
  ge_width_wo->SetName("ge_width_wo");

  TCanvas *c1 = new TCanvas("c1", "c1");
  
  ge_bias_w->GetYaxis()->SetRangeUser(-0.01, 0.25);
  ge_bias_w->Draw("AP");
  ge_bias_w->SetMarkerStyle(20);
  ge_bias_w->SetMarkerColor(kRed);

  ge_bias_w->SetTitle(";p#beta [MeV/c];Fractional bias");
  ge_bias_w->GetYaxis()->SetTitleOffset(1.3);
  ge_bias_w->GetYaxis()->CenterTitle();
  ge_bias_w->GetXaxis()->CenterTitle();

  ge_bias_wo->Draw("P SAME");
  ge_bias_wo->SetMarkerStyle(20);
  ge_bias_wo->SetMarkerColor(kBlue);

  l->AddEntry(ge_bias_wo, "#DeltaE = 0", "lp");
  l->AddEntry(ge_bias_w,  "#DeltaE follows the 4 order polynomial");
  l->SetName("l");
  l->Draw();

  gPad->SetGrid();

  TCanvas *c2 = new TCanvas("c2", "c2");

  ge_width_w->GetYaxis()->SetRangeUser(0., 0.3);
  ge_width_w->Draw("AP");
  ge_width_w->SetMarkerStyle(20);
  ge_width_w->SetMarkerColor(kRed);

  ge_width_w->SetTitle(";p#beta [MeV/c];Fractional width");
  ge_width_w->GetYaxis()->SetTitleOffset(1.3);
  ge_width_w->GetYaxis()->CenterTitle();
  ge_width_w->GetXaxis()->CenterTitle();

  ge_width_wo->Draw("P SAME");
  ge_width_wo->SetMarkerStyle(20);
  ge_width_wo->SetMarkerColor(kBlue);

  l->Draw();
  
  gPad->SetGrid();

  TFile *ofile = new TFile("~/enedep_graphs.root", "recreate");
  ofile->cd();
  ge_bias_w->Write();
  ge_bias_wo->Write();
  ge_width_w->Write();
  ge_width_wo->Write();
  l->Write();
  ofile->Close();

}
