void DrawAngleResMom() {

  TString filename_wo = "/home/t2k/odagawa/muon_wo_angres_mom_graph.root";
  TString filename_w = "/home/t2k/odagawa/muon_w_angres_mom_graph.root";

  TFile *file_wo = new TFile(filename_wo, "read");
  TFile *file_w = new TFile(filename_w, "read");

  TLegend *l = new TLegend(0.15, 0.75, 0.4, 0.89);

  file_w->cd();
  TGraphErrors *ge_bias_w = (TGraphErrors*)file_w->Get("ge_bias");
  TGraphErrors *ge_width_w = (TGraphErrors*)file_w->Get("ge_width");

  file_wo->cd();
  TGraphErrors *ge_bias_wo = (TGraphErrors*)file_wo->Get("ge_bias");
  TGraphErrors *ge_width_wo = (TGraphErrors*)file_wo->Get("ge_width");

  TCanvas *c1 = new TCanvas("c1", "c1");

  ge_bias_w->GetYaxis()->SetRangeUser(-0.01, 0.02);
  ge_bias_w->Draw("AP");
  ge_bias_w->SetMarkerStyle(20);
  ge_bias_w->SetMarkerColor(kRed);
  
  ge_bias_w->SetTitle(";p#beta [MeV/c];Fractional Bias");
  ge_bias_w->GetYaxis()->SetTitleOffset(1.3);
  ge_bias_w->GetYaxis()->CenterTitle();
  ge_bias_w->GetXaxis()->CenterTitle();

  ge_bias_wo->Draw("P SAME");
  ge_bias_wo->SetMarkerStyle(20);
  ge_bias_wo->SetMarkerColor(kBlue);

  l->AddEntry(ge_bias_wo, "w/o angular resolution", "lp");
  l->AddEntry(ge_bias_w,  "w/  angular resolution", "lp");

  l->Draw();
  
  gPad->SetGrid();

  TCanvas *c2 = new TCanvas("c2", "c2");

  ge_width_w->GetYaxis()->SetRangeUser(0., 0.3);
  ge_width_w->Draw("AP");
  ge_width_w->SetMarkerStyle(20);
  ge_width_w->SetMarkerColor(kRed);
  
  ge_width_w->SetTitle(";p#beta [MeV/c];Fractional Width");
  ge_width_w->GetYaxis()->SetTitleOffset(1.3);
  ge_width_w->GetYaxis()->CenterTitle();
  ge_width_w->GetXaxis()->CenterTitle();

  ge_width_wo->Draw("P SAME");
  ge_width_wo->SetMarkerStyle(20);
  ge_width_wo->SetMarkerColor(kBlue);

  l->Draw();

  gPad->SetGrid();

}
