#ifndef PID_DATA_HPP
#define PID_DATA_HPP

#include <string>
#include <vector>
#include <map>

#include <boost/filesystem.hpp>

#include <TH1D.h>

#include "PidClass.hpp"

namespace fs = boost::filesystem;

const extern fs::path PID_LIKELIHOOD_DIRNAME;
const extern fs::path PID_LIKELIHOOD_FILENAME;
const extern fs::path PID_VPH_FILENAME;
const extern fs::path DIST_PI_FILENAME;

class PidData {

private:

  std::map<int, std::map<double, Pid_data_ns::DataPoint > > likelihood_param_map_;
  std::map<int, Pid_data_ns::VphFuncParam > vph_func_param_map_;

  std::map<int, std::map<double, Pid_data_ns::VphPionMip > > vph_pion_mip_map_;
  std::map<int, double > vph_pion_mip_mean_map_;
  std::map<int, double > vph_pion_mip_thr_map_;
  std::map<int, double > vph_pion_mip_bin_width_map_;
  std::map<int, TH1D* > vph_pion_mip_hist_map_;

public:

  explicit PidData(const std::string &file_dir_path);
  explicit PidData(const fs::path &file_dir_path);
  
  void ReadPidData(std::map<int, std::map<double, Pid_data_ns::DataPoint > > &data,
		   const std::string file_dir_path);
  void ReadVphParamData(std::map<int, Pid_data_ns::VphFuncParam > &param,
			const std::string file_dir_path);

  void ReadPionPdfData(std::map<int, std::map<double, Pid_data_ns::VphPionMip > > &param_map,
		       const std::string file_dir_path);

  void GetLikelihoodParam(std::map<int, std::map<double, Pid_data_ns::DataPoint > > &likelihood_param_vec) const;
  void GetVphFuncParamMapData(std::map<int, Pid_data_ns::VphFuncParam > &vph_func_param_map) const;

  void GetPionMipParamMap(std::map<int, std::map<double, Pid_data_ns::VphPionMip > > &vph_pion_mip_map) const;
  void GetPionMipMeanMap(std::map<int, double > &vph_pion_mip_mean_map) const;
  void GetPionMipThrMap(std::map<int, double > &vph_pion_mip_thr_map) const;
  void GetPionMipBinWidthMap(std::map<int, double > &vph_pion_mip_bin_width_map) const;
  void GetPionPdfHistograms(std::map<int, TH1D* > &vph_pion_mip_hist_map) const;

private:

  void GeneratePionMipMeanMap();
  void GeneratePionMipThrMap();
  void GeneratePionMipProbMap();
  void GeneratePionMipBinWidthMap();
  void GeneratePionPdfHistograms();

};

#endif
