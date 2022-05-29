#ifndef PID_DATA_HPP
#define PID_DATA_HPP

#include <string>
#include <vector>
#include <map>

#include <boost/filesystem.hpp>

#include "PidClass.hpp"

namespace fs = boost::filesystem;

const extern fs::path PID_LIKELIHOOD_DIRNAME;
const extern fs::path PID_LIKELIHOOD_FILENAME;
const extern fs::path PID_VPH_FILENAME;

class PidData {

private:

  std::vector<Pid_data_ns::DataPoint > likelihood_param_vec_;
  std::vector<std::pair<double, double > > ang_bin_edge_vec_;
  std::map<std::pair<double, double >,
	   std::vector<std::pair<double, double> > > mom_bin_edge_map_;

  std::map<std::pair<double, double >, 
	   Pid_data_ns::VphFuncParam > vph_func_param_map_;
  

public:

  explicit PidData(const std::string &file_dir_path);
  explicit PidData(const fs::path &file_dir_path);
  
  void ReadPidData(std::vector<Pid_data_ns::DataPoint > &data,
		   std::vector<std::pair<double, double > > &ang_bin_edge_vec,
		   std::map<std::pair<double, double >, std::vector<std::pair<double, double > > > &mom_bin_edge_map,
		   const std::string file_dir_path);
  void ReadVphParamData(std::map<std::pair<double, double >, Pid_data_ns::VphFuncParam > &param,
			const std::string file_dir_path);

  void GetLikelihoodParam(std::vector<Pid_data_ns::DataPoint > &likelihood_param_vec) const;
  void GetAngBinVectorData(std::vector<std::pair<double, double > > &ang_bin_edge_vec) const;
  void GetMomBinMapData(std::map<std::pair<double, double >, std::vector<std::pair<double, double > >  > &mom_bin_edge_map) const;
  void GetVphFuncParamMapData(std::map<std::pair<double, double >, Pid_data_ns::VphFuncParam > &vph_func_param_map) const;

};

#endif
