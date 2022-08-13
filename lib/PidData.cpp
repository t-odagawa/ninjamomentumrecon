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
const fs::path DIST_PI_FILENAME("vph_dist_pi.txt");

PidData::PidData(const std::string &file_dir_path) {

  fs::path file_dir(file_dir_path);

  ReadPidData(likelihood_param_map_, file_dir.string());
  ReadVphParamData(vph_func_param_map_, file_dir.string());
  ReadPionPdfData(vph_pion_mip_map_, file_dir.string());

  GeneratePionMipMeanMap();
  GeneratePionMipThrMap();
  GeneratePionMipBinWidthMap();
  GeneratePionPdfHistograms();
  
  BOOST_LOG_TRIVIAL(info) << "Pid data are initialized";
  
}

PidData::PidData(const fs::path &file_dir_path) : PidData(file_dir_path.string()) {}

void PidData::ReadPidData(std::map<int, std::map<double, Pid_data_ns::DataPoint> > &data_map,
			  const std::string file_dir_path) {
  std::string file_path = (file_dir_path/PID_LIKELIHOOD_DIRNAME/PID_LIKELIHOOD_FILENAME).string();
  if ( !fs::exists(file_path) )
    throw std::runtime_error("File : " + file_path + " not found");

  std::ifstream ifs(file_path);
  BOOST_LOG_TRIVIAL(trace) << "Pid data file name : " << file_path;

  Pid_data_ns::DataPoint data;

  while ( ifs >> data ) {

    BOOST_LOG_TRIVIAL(trace) << data;
    
    double pbeta_mean = (data.input_mom_min + data.input_mom_max) / 2.;

    int iang = data.input_ang_max * 10;
    auto res = data_map.find(iang);
    if ( res == data_map.end() ) {
      std::map<double, Pid_data_ns::DataPoint > map_tmp;
      map_tmp.insert(std::make_pair(pbeta_mean, data));
      data_map.insert(std::make_pair(iang, map_tmp));
    }
    else res->second.insert(std::make_pair(pbeta_mean, data));

  }

  return;

}

void PidData::ReadVphParamData(std::map<int, Pid_data_ns::VphFuncParam > &param, const std::string file_dir_path) {

  std::string file_path = (file_dir_path/PID_LIKELIHOOD_DIRNAME/PID_VPH_FILENAME).string();
  if ( !fs::exists(file_path) )
    std::runtime_error("File : " + file_path + " not found");

  std::ifstream ifs(file_path);
  BOOST_LOG_TRIVIAL(trace) << "Pid VPH function parameter file name : " << file_path;

  Pid_data_ns::VphFuncParam data;
  while ( ifs >> data ) {
    BOOST_LOG_TRIVIAL(trace) << data;
    int iang = data.input_ang_max * 10;
    param.insert(std::make_pair(iang, data));
  }

  return;

}

void PidData::ReadPionPdfData(std::map<int, std::map<double, Pid_data_ns::VphPionMip > > &param_map,
			       const std::string file_dir_path) {

  std::string file_path = (file_dir_path/PID_LIKELIHOOD_DIRNAME/DIST_PI_FILENAME).string();
  if ( !fs::exists(file_path) )
    std::runtime_error("File : " + file_path + "not found");
  
  std::ifstream ifs(file_path);
  Pid_data_ns::VphPionMip data;
  while ( ifs >> data ) {

    BOOST_LOG_TRIVIAL(trace) << data;

    int iang = data.ang_max * 10;
    auto res = param_map.find(iang);

    if ( res == param_map.end() ) {
      std::map<double, Pid_data_ns::VphPionMip > map_tmp;
      map_tmp.insert(std::make_pair(data.vph, data));
      param_map.insert(std::make_pair(iang, map_tmp));
    }
    else {
      res->second.insert(std::make_pair(data.vph, data));
    }
    
  }

  return;

}

void PidData::GeneratePionMipMeanMap() {

  for ( auto itr = vph_pion_mip_map_.begin(); itr != vph_pion_mip_map_.end(); itr++ ) {
    auto map = itr->second;
    auto param = map.begin()->second;
    
    vph_pion_mip_mean_map_.insert(std::make_pair(itr->first, param.expect));
  }

  return;

}

void PidData::GeneratePionMipThrMap() {

  for ( auto itr = vph_pion_mip_map_.begin(); itr != vph_pion_mip_map_.end(); itr++ ) {
    auto map = itr->second;
    auto param = map.begin()->second;

    vph_pion_mip_thr_map_.insert(std::make_pair(itr->first, param.pb_max));
  }

  return;

}

void PidData::GeneratePionMipBinWidthMap() {

  for ( auto itr = vph_pion_mip_map_.begin(); itr != vph_pion_mip_map_.end(); itr++ ) {
    auto map = itr->second;
    auto param1 = map.begin()->second;
    auto param2 = std::next(map.begin(), 1)->second;

    vph_pion_mip_bin_width_map_.insert(std::make_pair(itr->first, param2.vph - param1.vph));

  }

  return;

}

void PidData::GeneratePionPdfHistograms() {
  
  for ( auto itr = vph_pion_mip_map_.begin(); itr != vph_pion_mip_map_.end(); itr++ ) {
    auto map = itr->second;
    auto param1 = map.begin()->second;
    auto param2 = map.rbegin()->second;
    double start_vph = param1.vph;
    double end_vph = param2.vph;
    double bin_width = vph_pion_mip_bin_width_map_.at(itr->first);

    TH1D *tmp = new TH1D(Form("tmp_%d", itr->first),
			 Form("%.1f < tan#theta < %.1f, %3.0f < p#beta < %3.0f",
			      param1.ang_min, param1.ang_max, param1.pb_max, param1.pb_max + 200.),
			 map.size(), start_vph - 0.5 * bin_width, end_vph + 0.5 * bin_width);

    for ( auto itr1 = map.begin(); itr1 != map.end(); itr1++ ) {
      tmp->Fill(itr1->second.vph, itr1->second.entry);
    }
    
    vph_pion_mip_hist_map_.insert(std::make_pair(itr->first, tmp));
  }
  
  return;

}

void PidData::GetLikelihoodParam(std::map<int, std::map<double, Pid_data_ns::DataPoint > > &likelihood_param_map) const {
  likelihood_param_map = likelihood_param_map_;
  return;
}

void PidData::GetVphFuncParamMapData(std::map<int, Pid_data_ns::VphFuncParam > &vph_func_param_map) const {
  vph_func_param_map = vph_func_param_map_;
  return;
}

void PidData::GetPionMipParamMap(std::map<int, std::map<double, Pid_data_ns::VphPionMip > > &vph_pion_mip_map) const {
  vph_pion_mip_map = vph_pion_mip_map_;
  return;
}

void PidData::GetPionMipMeanMap(std::map<int, double > &vph_pion_mip_mean_map) const {
  vph_pion_mip_mean_map = vph_pion_mip_mean_map_;
  return;
}

void PidData::GetPionMipThrMap(std::map<int, double > &vph_pion_mip_thr_map) const {
  vph_pion_mip_thr_map = vph_pion_mip_thr_map_;
  return;
}

void PidData::GetPionMipBinWidthMap(std::map<int, double > &vph_pion_mip_bin_width_map) const {
  vph_pion_mip_bin_width_map = vph_pion_mip_bin_width_map_;
  return;
}

void PidData::GetPionPdfHistograms(std::map<int, TH1D* > &vph_pion_mip_hist_map) const {
  vph_pion_mip_hist_map = vph_pion_mip_hist_map_;
  return;
}
