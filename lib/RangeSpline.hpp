#ifndef RANGE_SPLINE_HPP
#define RANGE_SPLINE_HPP

#include <boost/filesystem.hpp>

#include <string>

#include <TSpline.h>

namespace fs = boost::filesystem;

const extern fs::path IRON_FILENAME;
const extern fs::path WATER_FILENAME;
const extern fs::path POLY_FILENAME;
const extern fs::path EMULSION_FILENAME;

class RangeSpline {

private:
  TSpline3 *iron_energy_range_spline_;
  TSpline3 *iron_range_energy_spline_;
  TSpline3 *water_energy_range_spline_;
  TSpline3 *water_range_energy_spline_;
  TSpline3 *poly_energy_range_spline_;
  TSpline3 *poly_range_energy_spline_;
  TSpline3 *emulsion_energy_range_spline_;
  TSpline3 *emulsion_range_energy_spline_;

  void ReadIronSplines(const fs::path &file_dir_path);
  void ReadWaterSplines(const fs::path &file_dir_path);
  void ReadPolySplines(const fs::path &file_dir_path);
  void ReadEmulsionSplines(const fs::path &file_dir_path);

public:
  explicit RangeSpline(const std::string &file_dir_path);

  explicit RangeSpline(const fs::path &file_dir_path);

  void GetIronEnergyRangeSpline(TSpline3 &iron_energy_range_spline) const;
  void GetWaterEnergyRangeSpline(TSpline3 &water_energy_range_spline) const;
  void GetPolyEnergyRangeSpline(TSpline3 &poly_energy_range_spline) const;
  void GetEmulsionEnergyRangeSpline(TSpline3 &emulsion_energy_range_spline) const;
  void GetIronRangeEnergySpline(TSpline3 &iron_range_energy_spline) const;
  void GetWaterRangeEnergySpline(TSpline3 &water_range_energy_spline) const;
  void GetPolyRangeEnergySpline(TSpline3 &poly_range_energy_spline) const;
  void GetEmulsionRangeEnergySpline(TSpline3 &emulsion_range_energy_spline) const;

};

#endif
