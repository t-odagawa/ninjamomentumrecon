#include <boost/log/core.hpp>
#include <boost/log/trivial.hpp>
#include <boost/log/expressions.hpp>

#include <vector>
#include <map>

#include "PidData.hpp"
#include "PidClass.hpp"

namespace fs = boost::filesystem;

const fs::path PID_LIKELIHOOD_DIRNAME("pid");

const fs::path PID_LIKELIHOOD_FILENAME("likelihood_param.txt");
const fs::path PID_VPH_FILENAME("vph_func_param.txt");

PidData::PidData(const std::string &file_dir_path) {

  fs::path file_dir(file_dir_path);

  ReadPidData(likelihood_param_vec_, 
	      ang_bin_edge_vec_,
	      mom_bin_edge_map_,
	      file_dir.string());
  ReadVphParamData(vph_func_param_map_, file_dir.string());

  BOOST_LOG_TRIVIAL(info) << "Pid data are initialized";
  
}

PidData::PidData(const fs::path &file_dir_path) : PidData(file_dir_path.string()) {}

void PidData::ReadPidData(std::vector<Pid_data_ns::DataPoint > &data,
			  std::vector<std::pair<double, double > > &ang_bin_edge_vec,
			  std::map<std::pair<double, double >, std::vector<std::pair<double, double > > > &mom_bin_edge_map,
			  const std::string file_dir_path) {

  std::string file_path = (file_dir_path/PID_LIKELIHOOD_DIRNAME/PID_LIKELIHOOD_FILENAME).string();
  if ( !fs::exists(file_path) )
    throw std::runtime_error("File : " + file_path + " not found");

  std::ifstream ifs(file_path);
  BOOST_LOG_TRIVIAL(trace) << "Pid data file name : " << file_path;
  

  Pid_data_ns::DataPoint data_point_tmp;

  // 角度でソート，後運動量でソートされているファイルを仮定
  std::pair<double, double > ang_bin_edge_prev(-1., -1.);
  std::vector<std::pair<double, double > > mom_bin_edge_vec;

  while ( ifs >> data_point_tmp ) {
    auto ang_bin_edge_tmp = std::make_pair(data_point_tmp.input_ang_min,
					   data_point_tmp.input_ang_max);

    if ( data.size() != 0 && ang_bin_edge_tmp != ang_bin_edge_prev ) {
      ang_bin_edge_vec.push_back(ang_bin_edge_prev);
      mom_bin_edge_map.insert(std::make_pair(ang_bin_edge_prev, mom_bin_edge_vec));
      mom_bin_edge_vec.clear(); mom_bin_edge_vec.shrink_to_fit();
    }

    ang_bin_edge_prev = ang_bin_edge_tmp;
    
    mom_bin_edge_vec.push_back(std::make_pair(data_point_tmp.input_mom_min,
					      data_point_tmp.input_mom_max));

    data.push_back(data_point_tmp);
  }

  ang_bin_edge_vec.push_back(ang_bin_edge_prev);
  mom_bin_edge_map.insert(std::make_pair(ang_bin_edge_prev,
					 mom_bin_edge_vec));

  ifs.close();

  for ( auto itr = mom_bin_edge_map.begin(); itr != mom_bin_edge_map.end(); itr++ ) {
    BOOST_LOG_TRIVIAL(trace) << "Angle bin " << "(" << (*itr).first.first << ", " << (*itr).first.second << ")";
    for ( auto mom_bin : (*itr).second ) {
      BOOST_LOG_TRIVIAL(trace) << "(" << mom_bin.first << ", " << mom_bin.second << ")";
    }
  }

  for ( auto itr = data.begin(); itr != data.end(); itr++ ) {
    BOOST_LOG_TRIVIAL(trace) << (*itr);
  }

  return;

}

void PidData::ReadVphParamData(std::map<std::pair<double, double >, Pid_data_ns::VphFuncParam > &param,
			       const std::string file_dir_path) {

  std::string file_path = (file_dir_path/PID_LIKELIHOOD_DIRNAME/PID_VPH_FILENAME).string();
  if ( !fs::exists(file_path) )
    std::runtime_error("File : " + file_path + "not found");

  std::ifstream ifs(file_path);
  BOOST_LOG_TRIVIAL(trace) << "Pid VPH function parameter file name : " << file_path;

  Pid_data_ns::VphFuncParam vph_func_param_tmp_;
  while ( ifs >> vph_func_param_tmp_ ) {
    auto ang_bin = std::make_pair(vph_func_param_tmp_.input_ang_min,
				  vph_func_param_tmp_.input_ang_max);
    vph_func_param_map_.insert(std::make_pair(ang_bin, vph_func_param_tmp_));
  }

  for ( auto itr = vph_func_param_map_.begin(); itr != vph_func_param_map_.end(); itr++ ) {
    BOOST_LOG_TRIVIAL(trace) << (*itr).second;
  }

  return;

}

void PidData::GetLikelihoodParam(std::vector<Pid_data_ns::DataPoint > &likelihood_param_vec) const {
  likelihood_param_vec = likelihood_param_vec_;
  return;
}

void PidData::GetAngBinVectorData(std::vector<std::pair<double, double > > &ang_bin_edge_vec) const {
  ang_bin_edge_vec = ang_bin_edge_vec_;

  return;
}

void PidData::GetMomBinMapData(std::map<std::pair<double, double>, std::vector<std::pair<double, double > > > &mom_bin_edge_map) const {
  mom_bin_edge_map = mom_bin_edge_map_;
  return;
}

void PidData::GetVphFuncParamMapData(std::map<std::pair<double, double >, Pid_data_ns::VphFuncParam > &vph_func_param_map) const {
  vph_func_param_map = vph_func_param_map_;
  return;
}
