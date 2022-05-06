#ifndef CONNECTION_DATA_HPP
#define CONNECTION_DATA_HPP

#include <boost/filesystem.hpp>

#include <string>
#include <map>

#include <B2EmulsionSummary.hh>

#include "ConnectionClass.hpp"

namespace fs = boost::filesystem;

const extern fs::path CONNECT_DIRNAME;
const extern fs::path FIDUCIAL_DIRNAME;
const extern fs::path EFFICIENCY_DIRNAME;

const extern fs::path ET_ECC_CONNECT_FILENAME;
const extern fs::path ET_ECC_FE_CONNECT_FILENAME;
const extern fs::path ET_ECC_FE_FE_CONNECT_FILENAME;
const extern fs::path BLACK_FE_CONNECT_FILENAME;
const extern fs::path BLACK_RE_FE_CONNECT_FILENAME;
const extern fs::path BLACK_RE_FE_FE_CONNECT_FILENAME;
const extern fs::path BLACK_RE_FE_WATER_CONNECT_FILENAME;
const extern fs::path BLACK_WATER_CONNECT_FILENAME;
const extern fs::path FE_CONNECT_FILENAME;
const extern fs::path FE_FE_CONNECT_FILENAME;
const extern fs::path FE_FE_FE_CONNECT_FILENAME;
const extern fs::path FE_WATER_CONNECT_FILENAME;
const extern fs::path FE_WATER_FE_CONNECT_FILENAME;
const extern fs::path RE_ET_ECC_FE_FE_FE_CONNECT_FILENAME;
const extern fs::path RE_ET_ECC_FE_FE_FE_FE_CONNECT_FILENAME;
const extern fs::path RE_ET_ECC_FE_FE_FE_WATER_CONNECT_FILENAME;
const extern fs::path RE_ET_ECC_FE_FE_WATER_CONNECT_FILENAME;
const extern fs::path RE_ET_ECC_FE_WATER_CONNECT_FILENAME;
const extern fs::path RE_FE_CONNECT_FILENAME;
const extern fs::path RE_FE_FE_CONNECT_FILENAME;
const extern fs::path RE_FE_FE_FE_CONNECT_FILENAME;
const extern fs::path RE_FE_FE_FE_FE_CONNECT_FILENAME;
const extern fs::path RE_FE_FE_FE_FE_FE_CONNECT_FILENAME;
const extern fs::path RE_FE_WATER_CONNECT_FILENAME;
const extern fs::path RE_FE_WATER_FE_CONNECT_FILENAME;
const extern fs::path RE_FE_WATER_FE_WATER_CONNECT_FILENAME;
const extern fs::path RE_FE_WATER_FE_WATER_FE_CONNECT_FILENAME;
const extern fs::path RE_WATER_FE_WATER_CONNECT_FILENAME;
const extern fs::path WATER_CONNECT_FILENAME;

const extern fs::path ECC_FIDUCIAL_FILENAMES[9];

const extern fs::path ECC_EFFICIENCY_FILENAMES[9];

const int EFFICIENCY_VECTOR_SIZE = 13;

class ConnectionData {

private:
  
  t2l_param et_ecc_param_;
  t2l_param et_ecc_fe_param_;
  t2l_param et_ecc_fe_fe_param_;
  t2l_param black_fe_param_;
  t2l_param black_re_fe_param_;
  t2l_param black_re_fe_fe_param_;
  t2l_param black_re_fe_water_param_;
  t2l_param black_water_param_;
  t2l_param fe_param_;
  t2l_param fe_fe_param_;
  t2l_param fe_fe_fe_param_;
  t2l_param fe_water_param_;
  t2l_param fe_water_fe_param_;
  t2l_param re_et_ecc_fe_fe_fe_param_;
  t2l_param re_et_ecc_fe_fe_fe_fe_param_;
  t2l_param re_et_ecc_fe_fe_fe_water_param_;
  t2l_param re_et_ecc_fe_fe_water_param_;
  t2l_param re_et_ecc_fe_water_param_;
  t2l_param re_fe_param_;
  t2l_param re_fe_fe_param_;
  t2l_param re_fe_fe_fe_param_;
  t2l_param re_fe_fe_fe_fe_param_;
  t2l_param re_fe_fe_fe_fe_fe_param_;
  t2l_param re_fe_water_param_;
  t2l_param re_fe_water_fe_param_;
  t2l_param re_fe_water_fe_water_param_;
  t2l_param re_fe_water_fe_water_fe_param_;
  t2l_param re_water_fe_water_param_;
  t2l_param water_param_;

  std::map<int, std::vector<FiducialArea> > ecc_fiducial_[9];
  std::map<int, std::vector<Efficiency > > ecc_efficiency_[9];

  void ReadConnectData(const fs::path &file_dir_path);
  void ReadFiducialData(const fs::path &file_dir_path);
  void ReadEfficiencyData(const fs::path &file_dir_path);
  
  void ReadETECCConnectData(const fs::path &file_dir_path);
  void ReadETECCFeConnectData(const fs::path &file_dir_path);
  void ReadETECCFeFeConnectData(const fs::path &file_dir_path);
  void ReadBlackFeConnectData(const fs::path &file_dir_path);
  void ReadBlackReFeConnectData(const fs::path &file_dir_path);
  void ReadBlackReFeFeConnectData(const fs::path &file_dir_path);
  void ReadBlackReFeWaterConnectData(const fs::path &file_dir_path);
  void ReadBlackWaterConnectData(const fs::path &file_dir_path);
  void ReadFeConnectData(const fs::path &file_dir_path);
  void ReadFeFeConnectData(const fs::path &file_dir_path);
  void ReadFeFeFeConnectData(const fs::path &file_dir_path);
  void ReadFeWaterConnectData(const fs::path &file_dir_path);
  void ReadFeWaterFeConnectData(const fs::path &file_dir_path);
  void ReadReETECCFeFeFeConnectData(const fs::path &file_dir_path);
  void ReadReETECCFeFeFeFeConnectData(const fs::path &file_dir_path);
  void ReadReETECCFeFeFeWaterConnectData(const fs::path &file_dir_path);
  void ReadReETECCFeFeWaterConnectData(const fs::path &file_dir_path);
  void ReadReETECCFeWaterConnectData(const fs::path &file_dir_path);
  void ReadReFeConnectData(const fs::path &file_dir_path);
  void ReadReFeFeConnectData(const fs::path &file_dir_path);
  void ReadReFeFeFeConnectData(const fs::path &file_dir_path);
  void ReadReFeFeFeFeConnectData(const fs::path &file_dir_path);
  void ReadReFeFeFeFeFeConnectData(const fs::path &file_dir_path);
  void ReadReFeWaterConnectData(const fs::path &file_dir_path);
  void ReadReFeWaterFeConnectData(const fs::path &file_dir_path);
  void ReadReFeWaterFeWaterConnectData(const fs::path &file_dir_path);
  void ReadReFeWaterFeWaterFeConnectData(const fs::path &file_dir_path);
  void ReadReWaterFeWaterConnectData(const fs::path &file_dir_path);
  void ReadWaterConnectData(const fs::path &file_dir_path);

  void ReadJsonData(t2l_param &param, const std::string json_file_path);

  void ReadFiducialAreaEcc(int ecc, const fs::path &file_dir_path);
  void ReadFiducialArea(std::map<int, std::vector<FiducialArea > > &fiducial_area_map,
			const std::string file_dir_path);

  void ReadEfficiencyEcc(int ecc, const fs::path &file_dir_path);
  void ReadEfficiency(std::map<int, std::vector<Efficiency > > &efficiency_map,
		      const std::string file_dir_path);

public:
  explicit ConnectionData(const std::string &file_dir_path);

  explicit ConnectionData(const fs::path &file_dir_path);

  void GetETECCConnectData(t2l_param &et_ecc_param) const;
  void GetETECCFeConnectData(t2l_param &et_ecc_fe_param) const;
  void GetETECCFeFeConnectData(t2l_param &et_ecc_fe_fe_param) const;
  void GetBlackFeConnectData(t2l_param &black_fe_param) const;
  void GetBlackReFeConnectData(t2l_param &black_re_fe_param) const;
  void GetBlackReFeFeConnectData(t2l_param &black_re_fe_fe_param) const;
  void GetBlackReFeWaterConnectData(t2l_param &black_re_fe_water_param) const;
  void GetBlackWaterConnectData(t2l_param &black_water_param) const;
  void GetFeConnectData(t2l_param &fe_param) const;
  void GetFeFeConnectData(t2l_param &fe_fe_param) const;
  void GetFeFeFeConnectData(t2l_param &fe_fe_fe_param) const;
  void GetFeWaterConnectData(t2l_param &fe_water_param) const;
  void GetFeWaterFeConnectData(t2l_param &fe_water_fe_param) const;
  void GetReETECCFeFeFeConnectData(t2l_param &re_et_ecc_fe_fe_fe_param) const;
  void GetReETECCFeFeFeFeConnectData(t2l_param &re_et_ecc_fe_fe_fe_fe_param) const;
  void GetReETECCFeFeFeWaterConnectData(t2l_param &re_et_ecc_fe_fe_fe_water_param) const;
  void GetReETECCFeFeWaterConnectData(t2l_param &re_et_ecc_fe_fe_water_param) const;
  void GetReETECCFeWaterConnectData(t2l_param &re_et_ecc_fe_water_param) const;
  void GetReFeConnectData(t2l_param &re_fe_param) const;
  void GetReFeFeConnectData(t2l_param &re_fe_fe_param) const;
  void GetReFeFeFeConnectData(t2l_param &re_fe_fe_fe_param) const;
  void GetReFeFeFeFeConnectData(t2l_param &re_fe_fe_fe_fe_param) const;
  void GetReFeFeFeFeFeConnectData(t2l_param &re_fe_fe_fe_fe_fe_param) const;
  void GetReFeWaterConnectData(t2l_param &re_fe_water_param) const;
  void GetReFeWaterFeConnectData(t2l_param &re_fe_water_fe_param) const;
  void GetReFeWaterFeWaterConnectData(t2l_param &re_fe_water_fe_water_param) const;
  void GetReFeWaterFeWaterFeConnectData(t2l_param &re_fe_water_fe_water_fe_param) const;
  void GetReWaterFeWaterConnectData(t2l_param &re_water_fe_water_param) const;
  void GetWaterConnectData(t2l_param &water_param) const;

  void GetFiducialAreaData(int ecc, std::map<int, std::vector<FiducialArea > > &ecc_fiducial) const;
  void GetEfficiencyData(int ecc, std::map<int, std::vector<Efficiency > > &ecc_efficiency) const;

  void DrawFiducialAreaData(int ecc, int plate,
			    std::vector<B2EmulsionSummary* > &emulsions,
			    int eventid) const;
  void DrawEfficiencyData(int ecc, int plate) const;
  
};

#endif
