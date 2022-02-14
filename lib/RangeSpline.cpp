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

const fs::path IRON_FILENAME("iron_range.root");
const fs::path WATER_FILENAME("water_range.root");
const fs::path POLY_FILENAME("polystyrene_range.root");
const fs::path EMULSION_FILENAME("emulsion_range.root");

RangeSpline::RangeSpline(const std::string &file_dir_path) {
  fs::path file_dir(file_dir_path);

  ReadIronSplines(file_dir_path);
  ReadWaterSplines(file_dir_path);
  ReadPolySplines(file_dir_path);
  ReadEmulsionSplines(file_dir_path);

  BOOST_LOG_TRIVIAL(info) << "Range splines are initialized";
}

RangeSpline::RangeSpline(const fs::path &file_dir_path) : RangeSpline(file_dir_path.string()) {}

void RangeSpline::ReadIronSplines(const fs::path &file_dir_path) {
  const std::string iron_file_path = (file_dir_path/IRON_FILENAME).string();
  
  if ( !fs::exists(iron_file_path))
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

  iron_energy_range_spline_ = new TSpline3("iron_energy_range_spline_", iron_energy_range_graph);
  iron_range_energy_spline_ = new TSpline3("iron_range_energy_spline_", iron_range_energy_graph);

  iron_file.Close();
  return;

}

void RangeSpline::ReadWaterSplines(const fs::path &file_dir_path) {
  const std::string water_file_path = (file_dir_path/WATER_FILENAME).string();

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
  water_energy_range_spline_ = new TSpline3("water_energy_range_spline_", water_energy_range_graph);
  water_range_energy_spline_ = new TSpline3("water_range_energy_spline_", water_range_energy_graph);

  water_file.Close();
  return;
  
}

void RangeSpline::ReadPolySplines(const fs::path &file_dir_path) {
  const std::string poly_file_path = (file_dir_path/POLY_FILENAME).string();

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

  poly_energy_range_spline_ = new TSpline3("poly_energy_range_spline_", poly_energy_range_graph);
  poly_range_energy_spline_ = new TSpline3("poly_range_energy_spline_", poly_range_energy_graph);

  poly_file.Close();
  return;

}

void RangeSpline::ReadEmulsionSplines(const fs::path &file_dir_path) {
  const std::string emulsion_file_path = (file_dir_path/EMULSION_FILENAME).string();

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

  emulsion_energy_range_spline_ = new TSpline3("emulsion_energy_range_spline_", emulsion_energy_range_graph);
  emulsion_range_energy_spline_ = new TSpline3("emulsion_range_energy_spline_", emulsion_range_energy_graph);

  emulsion_file.Close();
  return;

}

void RangeSpline::GetIronEnergyRangeSpline(TSpline3 &iron_energy_range_spline) const {
  iron_energy_range_spline = *iron_energy_range_spline_;
  return;
}

void RangeSpline::GetWaterEnergyRangeSpline(TSpline3 &water_energy_range_spline) const {
  water_energy_range_spline = *water_energy_range_spline_;
  return;
}

void RangeSpline::GetPolyEnergyRangeSpline(TSpline3 &poly_energy_range_spline) const {
  poly_energy_range_spline = *poly_energy_range_spline_;
  return;
}

void RangeSpline::GetEmulsionEnergyRangeSpline(TSpline3 &emulsion_energy_range_spline) const {
  emulsion_energy_range_spline = *emulsion_energy_range_spline_;
  return;
}

void RangeSpline::GetIronRangeEnergySpline(TSpline3 &iron_range_energy_spline) const {
  iron_range_energy_spline = *iron_range_energy_spline_;
  return;
}

void RangeSpline::GetWaterRangeEnergySpline(TSpline3 &water_range_energy_spline) const {
  water_range_energy_spline = *water_range_energy_spline_;
  return;
}

void RangeSpline::GetPolyRangeEnergySpline(TSpline3 &poly_range_energy_spline) const {
  poly_range_energy_spline = *poly_range_energy_spline_;
  return;
}

void RangeSpline::GetEmulsionRangeEnergySpline(TSpline3 &emulsion_range_energy_spline) const {
  emulsion_range_energy_spline = *emulsion_range_energy_spline_;
  return;
}
