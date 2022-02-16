#ifndef RANGE_SPLINE_HPP
#define RANGE_SPLINE_HPP

#include <boost/filesystem.hpp>

#include <string>

#include <TSpline.h>

namespace fs = boost::filesystem;

const extern fs::path PROTON_DIRNAME;
const extern fs::path PION_DIRNAME;
const extern fs::path PROTON_IRON_FILENAME;
const extern fs::path PROTON_WATER_FILENAME;
const extern fs::path PROTON_POLY_FILENAME;
const extern fs::path PROTON_EMULSION_FILENAME;
const extern fs::path PION_IRON_FILENAME;
const extern fs::path PION_WATER_FILENAME;
const extern fs::path PION_POLY_FILENAME;
const extern fs::path PION_EMULSION_FILENAME;

class RangeSpline {

private:
  TSpline3 *proton_iron_energy_range_spline_;
  TSpline3 *proton_iron_range_energy_spline_;
  TSpline3 *proton_water_energy_range_spline_;
  TSpline3 *proton_water_range_energy_spline_;
  TSpline3 *proton_poly_energy_range_spline_;
  TSpline3 *proton_poly_range_energy_spline_;
  TSpline3 *proton_emulsion_energy_range_spline_;
  TSpline3 *proton_emulsion_range_energy_spline_;

  TSpline3 *pion_iron_energy_range_spline_;
  TSpline3 *pion_iron_range_energy_spline_;
  TSpline3 *pion_water_energy_range_spline_;
  TSpline3 *pion_water_range_energy_spline_;
  TSpline3 *pion_poly_energy_range_spline_;
  TSpline3 *pion_poly_range_energy_spline_;
  TSpline3 *pion_emulsion_energy_range_spline_;
  TSpline3 *pion_emulsion_range_energy_spline_;

  void ReadSplines(const fs::path &file_dir_path);

  void ReadProtonSplines(const fs::path &file_dir_path);
  void ReadPionSplines(const fs::path &file_dir_path);

  void ReadProtonIronSplines(const fs::path &file_dir_path);
  void ReadProtonWaterSplines(const fs::path &file_dir_path);
  void ReadProtonPolySplines(const fs::path &file_dir_path);
  void ReadProtonEmulsionSplines(const fs::path &file_dir_path);

  void ReadPionIronSplines(const fs::path &file_dir_path);
  void ReadPionWaterSplines(const fs::path &file_dir_path);
  void ReadPionPolySplines(const fs::path &file_dir_path);
  void ReadPionEmulsionSplines(const fs::path &file_dir_path);

public:
  explicit RangeSpline(const std::string &file_dir_path);

  explicit RangeSpline(const fs::path &file_dir_path);

  void GetProtonIronEnergyRangeSpline(TSpline3 &proton_iron_energy_range_spline) const;
  void GetProtonWaterEnergyRangeSpline(TSpline3 &proton_water_energy_range_spline) const;
  void GetProtonPolyEnergyRangeSpline(TSpline3 &proton_poly_energy_range_spline) const;
  void GetProtonEmulsionEnergyRangeSpline(TSpline3 &proton_emulsion_energy_range_spline) const;
  void GetProtonIronRangeEnergySpline(TSpline3 &proton_iron_range_energy_spline) const;
  void GetProtonWaterRangeEnergySpline(TSpline3 &proton_water_range_energy_spline) const;
  void GetProtonPolyRangeEnergySpline(TSpline3 &proton_poly_range_energy_spline) const;
  void GetProtonEmulsionRangeEnergySpline(TSpline3 &proton_emulsion_range_energy_spline) const;

  void GetPionIronEnergyRangeSpline(TSpline3 &pion_iron_energy_range_spline) const;
  void GetPionWaterEnergyRangeSpline(TSpline3 &pion_water_energy_range_spline) const;
  void GetPionPolyEnergyRangeSpline(TSpline3 &pion_poly_energy_range_spline) const;
  void GetPionEmulsionEnergyRangeSpline(TSpline3 &pion_emulsion_energy_range_spline) const;
  void GetPionIronRangeEnergySpline(TSpline3 &pion_iron_range_energy_spline) const;
  void GetPionWaterRangeEnergySpline(TSpline3 &pion_water_range_energy_spline) const;
  void GetPionPolyRangeEnergySpline(TSpline3 &pion_poly_range_energy_spline) const;
  void GetPionEmulsionRangeEnergySpline(TSpline3 &pion_emulsion_range_energy_spline) const;

};

#endif
