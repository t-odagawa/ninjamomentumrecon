void DrawSTune() {

  TString filename_wo = "/home/t2k/odagawa/muon_wo_s_tune_graph.root";
  TString filename_w = "/home/t2k/odagawa/muon_w_s_tune_graph.root";
  
  TFile *file_wo = new TFile(filename_wo, "read");
  TFile *file_w = new TFile(filename_w, "read");

  TLegend *l = new TLegend(0.15, 0.75, 0.4, 0.89);

  file_w->cd();
  TGraphErrors *ge_w = (TGraphErrors*)file_w->Get("ge");
  
  ge_w->GetYaxis()->SetRangeUser(-0.01, 0.04);

  ge_w->Draw("AP");

  ge_w->SetMarkerStyle(20);
  ge_w->SetMarkerColor(kRed);

  ge_w->SetTitle(";p#beta [MeV/c];Fractional Bias");
  ge_w->GetYaxis()->SetTitleOffset(1.3);
  ge_w->GetYaxis()->CenterTitle();
  ge_w->GetXaxis()->CenterTitle();

  file_wo->cd();
  TGraphErrors *ge_wo = (TGraphErrors*)file_wo->Get("ge");
  ge_wo->SetMarkerStyle(20);
  ge_wo->SetMarkerColor(kBlue);

  ge_wo->Draw("P SAME");

  l->AddEntry(ge_wo, "w/o S tuning", "lp");
  l->AddEntry(ge_w,  "w/  S tuning", "lp");

  l->Draw();

  gPad->SetGrid();

}
