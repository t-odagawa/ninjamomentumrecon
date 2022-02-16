#include "RangeSpline.hpp"

#include <boost/log/core.hpp>
#include <boost/log/trivial.hpp>
#include <boost/log/expressions.hpp>
#include <boost/filesystem.hpp>

#include <string>

#include <TFile.h>
#include <TTree.h>
#include <TGraph.h>
#include <TSpline.h>

#include "McsConst.hpp"

const fs::path PROTON_DIRNAME("proton");
const fs::path PION_DIRNAME("pion");

const fs::path PROTON_IRON_FILENAME("iron_range.root");
const fs::path PROTON_WATER_FILENAME("water_range.root");
const fs::path PROTON_POLY_FILENAME("polystyrene_range.root");
const fs::path PROTON_EMULSION_FILENAME("emulsion_range.root");

const fs::path PION_IRON_FILENAME("iron_range.root");
const fs::path PION_WATER_FILENAME("water_range.root");
const fs::path PION_POLY_FILENAME("polystyrene_range.root");
const fs::path PION_EMULSION_FILENAME("emulsion_range.root");

RangeSpline::RangeSpline(const std::string &file_dir_path) {
  fs::path file_dir(file_dir_path);

  ReadProtonSplines(file_dir_path);
  ReadPionSplines(file_dir_path);

  BOOST_LOG_TRIVIAL(info) << "Range splines are initialized";
}

RangeSpline::RangeSpline(const fs::path &file_dir_path) : RangeSpline(file_dir_path.string()) {}

void RangeSpline::ReadProtonSplines(const fs::path &file_dir_path) {
  const std::string proton_dir_path = (file_dir_path/PROTON_DIRNAME).string();

  if ( !fs::exists(proton_dir_path) )
    throw std::runtime_error("Directory : " + proton_dir_path + "not found");

  ReadProtonIronSplines(proton_dir_path);
  ReadProtonWaterSplines(proton_dir_path);
  ReadProtonPolySplines(proton_dir_path);
  ReadProtonEmulsionSplines(proton_dir_path);
}

void RangeSpline::ReadPionSplines(const fs::path &file_dir_path) {
  const std::string pion_dir_path = (file_dir_path/PION_DIRNAME).string();

  if ( !fs::exists(pion_dir_path) )
    throw std::runtime_error("Directory : " + pion_dir_path + "not found");

  ReadPionIronSplines(pion_dir_path);
  ReadPionWaterSplines(pion_dir_path);
  ReadPionPolySplines(pion_dir_path);
  ReadPionEmulsionSplines(pion_dir_path);
}

void RangeSpline::ReadProtonIronSplines(const fs::path &file_dir_path) {
  const std::string iron_file_path = (file_dir_path/PROTON_IRON_FILENAME).string();
  
  if ( !fs::exists(iron_file_path) )
    throw std::runtime_error("File : " + iron_file_path + " not found");

  TFile iron_file(iron_file_path.c_str());
  auto *iron_tree = (TTree*)iron_file.Get("tree");
  TGraph *iron_energy_range_graph = new TGraph();
  TGraph *iron_range_energy_graph = new TGraph();
  Double_t range, energy;
  iron_tree->SetBranchAddress("range", &range);
  iron_tree->SetBranchAddress("energy", &energy);
  for ( Int_t ientry = 0; ientry < iron_tree->GetEntries(); ientry++ ) {
    iron_tree->GetEntry(ientry);
    iron_energy_range_graph->SetPoint(ientry, energy, range / IRON_DENSITY * 10.);
    iron_range_energy_graph->SetPoint(ientry, range / IRON_DENSITY * 10., energy);
  }

  proton_iron_energy_range_spline_ = new TSpline3("proton_iron_energy_range_spline_", iron_energy_range_graph);
  proton_iron_range_energy_spline_ = new TSpline3("proton_iron_range_energy_spline_", iron_range_energy_graph);

  iron_file.Close();
  return;

}

void RangeSpline::ReadProtonWaterSplines(const fs::path &file_dir_path) {
  const std::string water_file_path = (file_dir_path/PROTON_WATER_FILENAME).string();

  if ( !fs::exists(water_file_path))
    throw std::runtime_error("File : " + water_file_path + " not found");

  TFile water_file(water_file_path.c_str());
  auto *water_tree = (TTree*)water_file.Get("tree");
  TGraph *water_energy_range_graph = new TGraph();
  TGraph *water_range_energy_graph = new TGraph();
  Double_t range, energy;
  water_tree->SetBranchAddress("range", &range);
  water_tree->SetBranchAddress("energy", &energy);
  for ( Int_t ientry = 0; ientry < water_tree->GetEntries(); ientry++ ) {
    water_tree->GetEntry(ientry);
    water_energy_range_graph->SetPoint(ientry, energy, range / WATER_DENSITY * 10.);
    water_range_energy_graph->SetPoint(ientry, range / WATER_DENSITY * 10., energy);
  }
  proton_water_energy_range_spline_ = new TSpline3("proton_water_energy_range_spline_", water_energy_range_graph);
  proton_water_range_energy_spline_ = new TSpline3("proton_water_range_energy_spline_", water_range_energy_graph);

  water_file.Close();
  return;
  
}

void RangeSpline::ReadProtonPolySplines(const fs::path &file_dir_path) {
  const std::string poly_file_path = (file_dir_path/PROTON_POLY_FILENAME).string();

  if ( !fs::exists(poly_file_path))
    throw std::runtime_error("File : " + poly_file_path + " not found");

  TFile poly_file(poly_file_path.c_str());
  auto *poly_tree = (TTree*)poly_file.Get("tree");
  TGraph *poly_energy_range_graph = new TGraph();
  TGraph *poly_range_energy_graph = new TGraph();
  Double_t range, energy;
  poly_tree->SetBranchAddress("range", &range);
  poly_tree->SetBranchAddress("energy", &energy);
  for ( Int_t ientry = 0; ientry < poly_tree->GetEntries(); ientry++ ) {
    poly_tree->GetEntry(ientry);
    poly_energy_range_graph->SetPoint(ientry, energy, range / POLY_DENSITY * 10.);
    poly_range_energy_graph->SetPoint(ientry, range / POLY_DENSITY * 10., energy);
  }

  proton_poly_energy_range_spline_ = new TSpline3("proton_poly_energy_range_spline_", poly_energy_range_graph);
  proton_poly_range_energy_spline_ = new TSpline3("proton_poly_range_energy_spline_", poly_range_energy_graph);

  poly_file.Close();
  return;

}

void RangeSpline::ReadProtonEmulsionSplines(const fs::path &file_dir_path) {
  const std::string emulsion_file_path = (file_dir_path/PROTON_EMULSION_FILENAME).string();

  if ( !fs::exists(emulsion_file_path))
    throw std::runtime_error("File : " + emulsion_file_path + " not found");

  TFile emulsion_file(emulsion_file_path.c_str());
  auto *emulsion_tree = (TTree*)emulsion_file.Get("tree");
  TGraph *emulsion_energy_range_graph = new TGraph();
  TGraph *emulsion_range_energy_graph = new TGraph();
  Double_t range, energy;
  emulsion_tree->SetBranchAddress("range", &range);
  emulsion_tree->SetBranchAddress("energy", &energy);
  for ( Int_t ientry = 0; ientry < emulsion_tree->GetEntries(); ientry++ ) {
    emulsion_tree->GetEntry(ientry);
    emulsion_energy_range_graph->SetPoint(ientry, energy, range / GEL_DENSITY * 10.);
    emulsion_range_energy_graph->SetPoint(ientry, range / GEL_DENSITY * 10., energy);
  }

  proton_emulsion_energy_range_spline_ = new TSpline3("proton_emulsion_energy_range_spline_", emulsion_energy_range_graph);
  proton_emulsion_range_energy_spline_ = new TSpline3("proton_emulsion_range_energy_spline_", emulsion_range_energy_graph);

  emulsion_file.Close();
  return;

}

void RangeSpline::ReadPionIronSplines(const fs::path &file_dir_path) {
  const std::string iron_file_path = (file_dir_path/PION_IRON_FILENAME).string();
  
  if ( !fs::exists(iron_file_path) )
    throw std::runtime_error("File : " + iron_file_path + " not found");

  TFile iron_file(iron_file_path.c_str());
  auto *iron_tree = (TTree*)iron_file.Get("tree");
  TGraph *iron_energy_range_graph = new TGraph();
  TGraph *iron_range_energy_graph = new TGraph();
  Double_t range, energy;
  iron_tree->SetBranchAddress("range", &range);
  iron_tree->SetBranchAddress("energy", &energy);
  for ( Int_t ientry = 0; ientry < iron_tree->GetEntries(); ientry++ ) {
    iron_tree->GetEntry(ientry);
    iron_energy_range_graph->SetPoint(ientry, energy, range / IRON_DENSITY * 10.);
    iron_range_energy_graph->SetPoint(ientry, range / IRON_DENSITY * 10., energy);
  }

  pion_iron_energy_range_spline_ = new TSpline3("pion_iron_energy_range_spline_", iron_energy_range_graph);
  pion_iron_range_energy_spline_ = new TSpline3("pion_iron_range_energy_spline_", iron_range_energy_graph);

  iron_file.Close();
  return;

}

void RangeSpline::ReadPionWaterSplines(const fs::path &file_dir_path) {
  const std::string water_file_path = (file_dir_path/PION_WATER_FILENAME).string();

  if ( !fs::exists(water_file_path))
    throw std::runtime_error("File : " + water_file_path + " not found");

  TFile water_file(water_file_path.c_str());
  auto *water_tree = (TTree*)water_file.Get("tree");
  TGraph *water_energy_range_graph = new TGraph();
  TGraph *water_range_energy_graph = new TGraph();
  Double_t range, energy;
  water_tree->SetBranchAddress("range", &range);
  water_tree->SetBranchAddress("energy", &energy);
  for ( Int_t ientry = 0; ientry < water_tree->GetEntries(); ientry++ ) {
    water_tree->GetEntry(ientry);
    water_energy_range_graph->SetPoint(ientry, energy, range / WATER_DENSITY * 10.);
    water_range_energy_graph->SetPoint(ientry, range / WATER_DENSITY * 10., energy);
  }
  pion_water_energy_range_spline_ = new TSpline3("pion_water_energy_range_spline_", water_energy_range_graph);
  pion_water_range_energy_spline_ = new TSpline3("pion_water_range_energy_spline_", water_range_energy_graph);

  water_file.Close();
  return;
  
}

void RangeSpline::ReadPionPolySplines(const fs::path &file_dir_path) {
  const std::string poly_file_path = (file_dir_path/PION_POLY_FILENAME).string();

  if ( !fs::exists(poly_file_path))
    throw std::runtime_error("File : " + poly_file_path + " not found");

  TFile poly_file(poly_file_path.c_str());
  auto *poly_tree = (TTree*)poly_file.Get("tree");
  TGraph *poly_energy_range_graph = new TGraph();
  TGraph *poly_range_energy_graph = new TGraph();
  Double_t range, energy;
  poly_tree->SetBranchAddress("range", &range);
  poly_tree->SetBranchAddress("energy", &energy);
  for ( Int_t ientry = 0; ientry < poly_tree->GetEntries(); ientry++ ) {
    poly_tree->GetEntry(ientry);
    poly_energy_range_graph->SetPoint(ientry, energy, range / POLY_DENSITY * 10.);
    poly_range_energy_graph->SetPoint(ientry, range / POLY_DENSITY * 10., energy);
  }

  pion_poly_energy_range_spline_ = new TSpline3("pion_poly_energy_range_spline_", poly_energy_range_graph);
  pion_poly_range_energy_spline_ = new TSpline3("pion_poly_range_energy_spline_", poly_range_energy_graph);

  poly_file.Close();
  return;

}

void RangeSpline::ReadPionEmulsionSplines(const fs::path &file_dir_path) {
  const std::string emulsion_file_path = (file_dir_path/PION_EMULSION_FILENAME).string();

  if ( !fs::exists(emulsion_file_path))
    throw std::runtime_error("File : " + emulsion_file_path + " not found");

  TFile emulsion_file(emulsion_file_path.c_str());
  auto *emulsion_tree = (TTree*)emulsion_file.Get("tree");
  TGraph *emulsion_energy_range_graph = new TGraph();
  TGraph *emulsion_range_energy_graph = new TGraph();
  Double_t range, energy;
  emulsion_tree->SetBranchAddress("range", &range);
  emulsion_tree->SetBranchAddress("energy", &energy);
  for ( Int_t ientry = 0; ientry < emulsion_tree->GetEntries(); ientry++ ) {
    emulsion_tree->GetEntry(ientry);
    emulsion_energy_range_graph->SetPoint(ientry, energy, range / GEL_DENSITY * 10.);
    emulsion_range_energy_graph->SetPoint(ientry, range / GEL_DENSITY * 10., energy);
  }

  pion_emulsion_energy_range_spline_ = new TSpline3("pion_emulsion_energy_range_spline_", emulsion_energy_range_graph);
  pion_emulsion_range_energy_spline_ = new TSpline3("pion_emulsion_range_energy_spline_", emulsion_range_energy_graph);

  emulsion_file.Close();
  return;

}


void RangeSpline::GetProtonIronEnergyRangeSpline(TSpline3 &proton_iron_energy_range_spline) const {
  proton_iron_energy_range_spline = *proton_iron_energy_range_spline_;
  return;
}

void RangeSpline::GetProtonWaterEnergyRangeSpline(TSpline3 &proton_water_energy_range_spline) const {
  proton_water_energy_range_spline = *proton_water_energy_range_spline_;
  return;
}

void RangeSpline::GetProtonPolyEnergyRangeSpline(TSpline3 &proton_poly_energy_range_spline) const {
  proton_poly_energy_range_spline = *proton_poly_energy_range_spline_;
  return;
}

void RangeSpline::GetProtonEmulsionEnergyRangeSpline(TSpline3 &proton_emulsion_energy_range_spline) const {
  proton_emulsion_energy_range_spline = *proton_emulsion_energy_range_spline_;
  return;
}

void RangeSpline::GetProtonIronRangeEnergySpline(TSpline3 &proton_iron_range_energy_spline) const {
  proton_iron_range_energy_spline = *proton_iron_range_energy_spline_;
  return;
}

void RangeSpline::GetProtonWaterRangeEnergySpline(TSpline3 &proton_water_range_energy_spline) const {
  proton_water_range_energy_spline = *proton_water_range_energy_spline_;
  return;
}

void RangeSpline::GetProtonPolyRangeEnergySpline(TSpline3 &proton_poly_range_energy_spline) const {
  proton_poly_range_energy_spline = *proton_poly_range_energy_spline_;
  return;
}

void RangeSpline::GetProtonEmulsionRangeEnergySpline(TSpline3 &proton_emulsion_range_energy_spline) const {
  proton_emulsion_range_energy_spline = *proton_emulsion_range_energy_spline_;
  return;
}

void RangeSpline::GetPionIronEnergyRangeSpline(TSpline3 &pion_iron_energy_range_spline) const {
  pion_iron_energy_range_spline = *pion_iron_energy_range_spline_;
  return;
}

void RangeSpline::GetPionWaterEnergyRangeSpline(TSpline3 &pion_water_energy_range_spline) const {
  pion_water_energy_range_spline = *pion_water_energy_range_spline_;
  return;
}

void RangeSpline::GetPionPolyEnergyRangeSpline(TSpline3 &pion_poly_energy_range_spline) const {
  pion_poly_energy_range_spline = *pion_poly_energy_range_spline_;
  return;
}

void RangeSpline::GetPionEmulsionEnergyRangeSpline(TSpline3 &pion_emulsion_energy_range_spline) const {
  pion_emulsion_energy_range_spline = *pion_emulsion_energy_range_spline_;
  return;
}

void RangeSpline::GetPionIronRangeEnergySpline(TSpline3 &pion_iron_range_energy_spline) const {
  pion_iron_range_energy_spline = *pion_iron_range_energy_spline_;
  return;
}

void RangeSpline::GetPionWaterRangeEnergySpline(TSpline3 &pion_water_range_energy_spline) const {
  pion_water_range_energy_spline = *pion_water_range_energy_spline_;
  return;
}

void RangeSpline::GetPionPolyRangeEnergySpline(TSpline3 &pion_poly_range_energy_spline) const {
  pion_poly_range_energy_spline = *pion_poly_range_energy_spline_;
  return;
}

void RangeSpline::GetPionEmulsionRangeEnergySpline(TSpline3 &pion_emulsion_range_energy_spline) const {
  pion_emulsion_range_energy_spline = *pion_emulsion_range_energy_spline_;
  return;
}
