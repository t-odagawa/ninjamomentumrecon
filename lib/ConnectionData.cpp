#include "ConnectionData.hpp"

#include <boost/log/core.hpp>
#include <boost/log/trivial.hpp>
#include <boost/log/expressions.hpp>
#include <boost/filesystem.hpp>

#include <string>

#include <TCanvas.h>
#include <TGraph.h>
#include <TGraphErrors.h>

#include <B2EmulsionSummary.hh>

#include "picojson.h"

#include "ConnectionClass.hpp"

const fs::path CONNECT_DIRNAME("connection");
const fs::path FIDUCIAL_DIRNAME("fiducial");
const fs::path EFFICIENCY_DIRNAME("efficiency");

const fs::path ET_ECC_CONNECT_FILENAME("t2l_ET_ECC.json");
const fs::path ET_ECC_FE_CONNECT_FILENAME("t2l_ET_ECC_fe.json");
const fs::path ET_ECC_FE_FE_CONNECT_FILENAME("t2l_ET_ECC_fe_fe.json");
const fs::path BLACK_FE_CONNECT_FILENAME("t2l_black_fe.json");
const fs::path BLACK_RE_FE_CONNECT_FILENAME("t2l_black_re_fe.json");
const fs::path BLACK_RE_FE_FE_CONNECT_FILENAME("t2l_black_re_fe_fe.json");
const fs::path BLACK_RE_FE_WATER_CONNECT_FILENAME("t2l_black_re_fe_water.json");
const fs::path BLACK_WATER_CONNECT_FILENAME("t2l_black_water.json");
const fs::path FE_CONNECT_FILENAME("t2l_fe.json");
const fs::path FE_FE_CONNECT_FILENAME("t2l_fe_fe.json");
const fs::path FE_FE_FE_CONNECT_FILENAME("t2l_fe_fe_fe.json");
const fs::path FE_WATER_CONNECT_FILENAME("t2l_fe_water.json");
const fs::path FE_WATER_FE_CONNECT_FILENAME("t2l_fe_water_fe.json");
const fs::path RE_ET_ECC_FE_FE_FE_CONNECT_FILENAME("t2l_re_ET_ECC_fe_fe_fe.json");
const fs::path RE_ET_ECC_FE_FE_FE_FE_CONNECT_FILENAME("t2l_re_ET_ECC_fe_fe_fe_fe.json");
const fs::path RE_ET_ECC_FE_FE_FE_WATER_CONNECT_FILENAME("t2l_re_ET_ECC_fe_fe_fe_water.json");
const fs::path RE_ET_ECC_FE_FE_WATER_CONNECT_FILENAME("t2l_re_ET_ECC_fe_fe_water.json");
const fs::path RE_ET_ECC_FE_WATER_CONNECT_FILENAME("t2l_re_ET_ECC_fe_water.json");
const fs::path RE_FE_CONNECT_FILENAME("t2l_re_fe.json");
const fs::path RE_FE_FE_CONNECT_FILENAME("t2l_re_fe_fe.json");
const fs::path RE_FE_FE_FE_CONNECT_FILENAME("t2l_re_fe_fe_fe.json");
const fs::path RE_FE_FE_FE_FE_CONNECT_FILENAME("t2l_re_fe_fe_fe_fe.json");
const fs::path RE_FE_FE_FE_FE_FE_CONNECT_FILENAME("t2l_re_fe_fe_fe_fe_fe.json");
const fs::path RE_FE_WATER_CONNECT_FILENAME("t2l_re_fe_water.json");
const fs::path RE_FE_WATER_FE_CONNECT_FILENAME("t2l_re_fe_water_fe.json");
const fs::path RE_FE_WATER_FE_WATER_CONNECT_FILENAME("t2l_re_fe_water_fe_water.json");
const fs::path RE_FE_WATER_FE_WATER_FE_CONNECT_FILENAME("t2l_re_fe_water_fe_water_fe.json");
const fs::path RE_WATER_FE_WATER_CONNECT_FILENAME("t2l_re_water_fe_water.json");
const fs::path WATER_CONNECT_FILENAME("t2l_water.json");

const fs::path ECC_FIDUCIAL_FILENAMES[9] = {"fa_new_mfile_local_ecc1.txt",
					    "fa_new_mfile_local_ecc2.txt",
					    "fa_new_mfile_local_ecc3.txt",
					    "fa_new_mfile_local_ecc4.txt",
					    "fa_new_mfile_local_ecc5.txt",
					    "fa_new_mfile_local_ecc6.txt",
					    "fa_new_mfile_local_ecc7.txt",
					    "fa_new_mfile_local_ecc8.txt",
					    "fa_new_mfile_local_ecc9.txt"};

const fs::path ECC_EFFICIENCY_FILENAMES[9] = {"efficiency_ecc1.txt",
					      "efficiency_ecc2.txt",
					      "efficiency_ecc3.txt",
					      "efficiency_ecc4.txt",
					      "efficiency_ecc5.txt",
					      "efficiency_ecc6.txt",
					      "efficiency_ecc7.txt",
					      "efficiency_ecc8.txt",
					      "efficiency_ecc9.txt"};

ConnectionData::ConnectionData(const std::string &file_dir_path) {
  fs::path file_dir(file_dir_path);

  ReadConnectData(file_dir);
  ReadFiducialData(file_dir);
  ReadEfficiencyData(file_dir);

  BOOST_LOG_TRIVIAL(info) << "Connection data are initialized";
}

ConnectionData::ConnectionData(const fs::path &file_dir_path) : ConnectionData(file_dir_path.string()) {}

void ConnectionData::ReadConnectData(const fs::path &file_dir_path) {
  const std::string connect_dir_path = (file_dir_path/CONNECT_DIRNAME).string();

  if ( !fs::exists(connect_dir_path) )
    throw std::runtime_error("Directory : " + connect_dir_path + " not found");

  ReadETECCConnectData(connect_dir_path);
  ReadETECCFeConnectData(connect_dir_path);
  ReadETECCFeFeConnectData(connect_dir_path);
  ReadBlackFeConnectData(connect_dir_path);
  ReadBlackReFeConnectData(connect_dir_path);
  ReadBlackReFeFeConnectData(connect_dir_path);
  ReadBlackReFeWaterConnectData(connect_dir_path);
  ReadBlackWaterConnectData(connect_dir_path);
  ReadFeConnectData(connect_dir_path);
  ReadFeFeConnectData(connect_dir_path);
  ReadFeFeFeConnectData(connect_dir_path);
  ReadFeWaterConnectData(connect_dir_path);
  ReadFeWaterFeConnectData(connect_dir_path);
  ReadReETECCFeFeFeConnectData(connect_dir_path);
  ReadReETECCFeFeFeFeConnectData(connect_dir_path);
  ReadReETECCFeFeFeWaterConnectData(connect_dir_path);
  ReadReETECCFeFeWaterConnectData(connect_dir_path);
  ReadReETECCFeWaterConnectData(connect_dir_path);
  ReadReFeConnectData(connect_dir_path);
  ReadReFeFeConnectData(connect_dir_path);
  ReadReFeFeFeConnectData(connect_dir_path);
  ReadReFeFeFeFeConnectData(connect_dir_path);
  ReadReFeFeFeFeFeConnectData(connect_dir_path);
  ReadReFeWaterConnectData(connect_dir_path);
  ReadReFeWaterFeConnectData(connect_dir_path);
  ReadReFeWaterFeWaterConnectData(connect_dir_path);
  ReadReFeWaterFeWaterFeConnectData(connect_dir_path);
  ReadReWaterFeWaterConnectData(connect_dir_path);
  ReadWaterConnectData(connect_dir_path);

  return;

}

void ConnectionData::ReadFiducialData(const fs::path &file_dir_path) {
  const fs::path fiducial_dir_path(file_dir_path/FIDUCIAL_DIRNAME);
  
  if ( !fs::exists(fiducial_dir_path) )
    throw std::runtime_error("Directory : " + fiducial_dir_path.string() + " not found");

  for ( int i = 0; i < 9; i++ ) {
    if ( i != 4 ) continue;
    ReadFiducialAreaEcc(i, fiducial_dir_path);
  }
  
}

void ConnectionData::ReadEfficiencyData(const fs::path &file_dir_path) {
  const fs::path efficiency_dir_path(file_dir_path/EFFICIENCY_DIRNAME);

  if ( !fs::exists(efficiency_dir_path) )
    throw std::runtime_error("Directory :" + efficiency_dir_path.string() + " not found");

  for ( int i = 0; i < 9; i++ ) {
    if ( i != 4 ) continue;
    ReadEfficiencyEcc(i, efficiency_dir_path);
  }
}

void ConnectionData::ReadETECCConnectData(const fs::path &file_dir_path) {
  const std::string json_file_path = (file_dir_path/ET_ECC_CONNECT_FILENAME).string();
  ReadJsonData(et_ecc_param_, json_file_path);
  return;
}

void ConnectionData::ReadETECCFeConnectData(const fs::path &file_dir_path) {
  const std::string json_file_path = (file_dir_path/ET_ECC_FE_CONNECT_FILENAME).string();
  ReadJsonData(et_ecc_fe_param_, json_file_path);
  return;
}

void ConnectionData::ReadETECCFeFeConnectData(const fs::path &file_dir_path) {
  const std::string json_file_path = (file_dir_path/ET_ECC_FE_FE_CONNECT_FILENAME).string();
  ReadJsonData(et_ecc_fe_fe_param_, json_file_path);
  return;
}

void ConnectionData::ReadBlackFeConnectData(const fs::path &file_dir_path) {
  const std::string json_file_path = (file_dir_path/BLACK_FE_CONNECT_FILENAME).string();
  ReadJsonData(black_fe_param_, json_file_path);
  return;
}

void ConnectionData::ReadBlackReFeConnectData(const fs::path &file_dir_path) {
  const std::string json_file_path = (file_dir_path/BLACK_RE_FE_CONNECT_FILENAME).string();
  ReadJsonData(black_re_fe_param_, json_file_path);
  return;
}

void ConnectionData::ReadBlackReFeFeConnectData(const fs::path &file_dir_path) {
  const std::string json_file_path = (file_dir_path/BLACK_RE_FE_FE_CONNECT_FILENAME).string();
  ReadJsonData(black_re_fe_fe_param_, json_file_path);
  return;
}

void ConnectionData::ReadBlackReFeWaterConnectData(const fs::path &file_dir_path) {
  const std::string json_file_path = (file_dir_path/BLACK_RE_FE_WATER_CONNECT_FILENAME).string();
  ReadJsonData(black_re_fe_water_param_, json_file_path);
  return;
}

void ConnectionData::ReadBlackWaterConnectData(const fs::path &file_dir_path) {
  const std::string json_file_path = (file_dir_path/BLACK_WATER_CONNECT_FILENAME).string();
  ReadJsonData(black_water_param_, json_file_path);
  return;
}

void ConnectionData::ReadFeConnectData(const fs::path &file_dir_path) {
  const std::string json_file_path = (file_dir_path/FE_CONNECT_FILENAME).string();
  ReadJsonData(fe_param_, json_file_path);
  return;
}

void ConnectionData::ReadFeFeConnectData(const fs::path &file_dir_path) {
  const std::string json_file_path = (file_dir_path/FE_FE_CONNECT_FILENAME).string();
  ReadJsonData(fe_fe_param_, json_file_path);
  return;
}

void ConnectionData::ReadFeFeFeConnectData(const fs::path &file_dir_path) {
  const std::string json_file_path = (file_dir_path/FE_FE_FE_CONNECT_FILENAME).string();
  ReadJsonData(fe_fe_fe_param_, json_file_path);
  return;
}

void ConnectionData::ReadFeWaterConnectData(const fs::path &file_dir_path) {
  const std::string json_file_path = (file_dir_path/FE_WATER_CONNECT_FILENAME).string();
  ReadJsonData(fe_water_param_, json_file_path);
  return;
}

void ConnectionData::ReadFeWaterFeConnectData(const fs::path &file_dir_path){
  const std::string json_file_path = (file_dir_path/FE_WATER_FE_CONNECT_FILENAME).string();
  ReadJsonData(fe_water_fe_param_, json_file_path);
  return;
}

void ConnectionData::ReadReETECCFeFeFeConnectData(const fs::path &file_dir_path){
  const std::string json_file_path = (file_dir_path/RE_ET_ECC_FE_FE_FE_CONNECT_FILENAME).string();
  ReadJsonData(re_et_ecc_fe_fe_fe_param_, json_file_path);
  return;
}

void ConnectionData::ReadReETECCFeFeFeFeConnectData(const fs::path &file_dir_path){
  const std::string json_file_path = (file_dir_path/RE_ET_ECC_FE_FE_FE_FE_CONNECT_FILENAME).string();
  ReadJsonData(re_et_ecc_fe_fe_fe_fe_param_, json_file_path);
  return;
}

void ConnectionData::ReadReETECCFeFeFeWaterConnectData(const fs::path &file_dir_path){
  const std::string json_file_path = (file_dir_path/RE_ET_ECC_FE_FE_FE_WATER_CONNECT_FILENAME).string();
  ReadJsonData(re_et_ecc_fe_fe_fe_water_param_, json_file_path);
  return;
}

void ConnectionData::ReadReETECCFeFeWaterConnectData(const fs::path &file_dir_path){
  const std::string json_file_path = (file_dir_path/RE_ET_ECC_FE_FE_WATER_CONNECT_FILENAME).string();
  ReadJsonData(re_et_ecc_fe_fe_water_param_, json_file_path);
  return;
}

void ConnectionData::ReadReETECCFeWaterConnectData(const fs::path &file_dir_path){
  const std::string json_file_path = (file_dir_path/RE_ET_ECC_FE_WATER_CONNECT_FILENAME).string();
  ReadJsonData(re_et_ecc_fe_water_param_, json_file_path);
  return;
}

void ConnectionData::ReadReFeConnectData(const fs::path &file_dir_path){
  const std::string json_file_path = (file_dir_path/RE_FE_CONNECT_FILENAME).string();
  ReadJsonData(re_fe_param_, json_file_path);
  return;
}

void ConnectionData::ReadReFeFeConnectData(const fs::path &file_dir_path){
  const std::string json_file_path = (file_dir_path/RE_FE_FE_CONNECT_FILENAME).string();
  ReadJsonData(re_fe_fe_param_, json_file_path);
  return;
}

void ConnectionData::ReadReFeFeFeConnectData(const fs::path &file_dir_path){
  const std::string json_file_path = (file_dir_path/RE_FE_FE_FE_CONNECT_FILENAME).string();
  ReadJsonData(re_fe_fe_fe_param_, json_file_path);
  return;
}

void ConnectionData::ReadReFeFeFeFeConnectData(const fs::path &file_dir_path){
  const std::string json_file_path = (file_dir_path/RE_FE_FE_FE_FE_CONNECT_FILENAME).string();
  ReadJsonData(re_fe_fe_fe_fe_param_, json_file_path);
  return;
}

void ConnectionData::ReadReFeFeFeFeFeConnectData(const fs::path &file_dir_path){
  const std::string json_file_path = (file_dir_path/RE_FE_FE_FE_FE_FE_CONNECT_FILENAME).string();
  ReadJsonData(re_fe_fe_fe_fe_fe_param_, json_file_path);
  return;
}

void ConnectionData::ReadReFeWaterConnectData(const fs::path &file_dir_path){
  const std::string json_file_path = (file_dir_path/RE_FE_WATER_CONNECT_FILENAME).string();
  ReadJsonData(re_fe_water_param_, json_file_path);
  return;
}

void ConnectionData::ReadReFeWaterFeConnectData(const fs::path &file_dir_path){
  const std::string json_file_path = (file_dir_path/RE_FE_WATER_FE_CONNECT_FILENAME).string();
  ReadJsonData(re_fe_water_fe_param_, json_file_path);
  return;
}

void ConnectionData::ReadReFeWaterFeWaterConnectData(const fs::path &file_dir_path){
  const std::string json_file_path = (file_dir_path/RE_FE_WATER_FE_WATER_CONNECT_FILENAME).string();
  ReadJsonData(re_fe_water_fe_water_param_, json_file_path);
  return;
}

void ConnectionData::ReadReFeWaterFeWaterFeConnectData(const fs::path &file_dir_path){
  const std::string json_file_path = (file_dir_path/RE_FE_WATER_FE_WATER_FE_CONNECT_FILENAME).string();
  ReadJsonData(re_fe_water_fe_water_fe_param_, json_file_path);
  return;
}

void ConnectionData::ReadReWaterFeWaterConnectData(const fs::path &file_dir_path){
  const std::string json_file_path = (file_dir_path/RE_WATER_FE_WATER_CONNECT_FILENAME).string();
  ReadJsonData(re_water_fe_water_param_, json_file_path);
  return;
}

void ConnectionData::ReadWaterConnectData(const fs::path &file_dir_path){
  const std::string json_file_path = (file_dir_path/WATER_CONNECT_FILENAME).string();
  ReadJsonData(water_param_, json_file_path);
  return;
}

void ConnectionData::ReadJsonData(t2l_param &param, const std::string json_file_path) {

  if ( !fs::exists(json_file_path) )
    throw std::runtime_error("File : " + json_file_path + " not found");

  std::ifstream ifs(json_file_path);
  const std::string json((std::istreambuf_iterator<char >(ifs)), std::istreambuf_iterator<char >());
  ifs.close();

  picojson::value v;
  const std::string err = picojson::parse(v, json);
  if ( !err.empty() ) 
    throw std::runtime_error("JSON parse error : " + err);

  picojson::object &all = v.get<picojson::object >();
  picojson::object &connect_param_angle = all["connect_param_angle"].get<picojson::object >();
  picojson::object &connect_param_position = all["connect_param_position"].get<picojson::object >();

  param.intercept_ax = connect_param_angle["intercept_x"].get<double>();
  param.intercept_ay = connect_param_angle["intercept_y"].get<double>();
  param.intercept_ar = connect_param_angle["intercept_r"].get<double>();
  param.intercept_al = connect_param_angle["intercept_l"].get<double>();
  param.slope_ax = connect_param_angle["slope_x"].get<double>();
  param.slope_ay = connect_param_angle["slope_y"].get<double>();
  param.slope_ar = connect_param_angle["slope_r"].get<double>();
  param.slope_al = connect_param_angle["slope_l"].get<double>();
  param.slope2_ax = connect_param_angle["slope2_x"].get<double>();
  param.slope2_ay = connect_param_angle["slope2_y"].get<double>();
  param.slope2_ar = connect_param_angle["slope2_r"].get<double>();
  param.slope2_al = connect_param_angle["slope2_l"].get<double>();

  param.intercept_px = connect_param_position["intercept_x"].get<double>();
  param.intercept_py = connect_param_position["intercept_y"].get<double>();
  param.intercept_pr = connect_param_position["intercept_r"].get<double>();
  param.intercept_pl = connect_param_position["intercept_l"].get<double>();
  param.slope_px = connect_param_position["slope_x"].get<double>();
  param.slope_py = connect_param_position["slope_y"].get<double>();
  param.slope_pr = connect_param_position["slope_r"].get<double>();
  param.slope_pl = connect_param_position["slope_l"].get<double>();
  param.slope2_px = connect_param_position["slope2_x"].get<double>();
  param.slope2_py = connect_param_position["slope2_y"].get<double>();
  param.slope2_pr = connect_param_position["slope2_r"].get<double>();
  param.slope2_pl = connect_param_position["slope2_l"].get<double>();

  BOOST_LOG_TRIVIAL(trace) << "Filename : "      << json_file_path << "\n"
			   << "intercept_ax : "  << param.intercept_ax << ", "
			   << "intercept_ay : "  << param.intercept_ay << ", "
			   << "intercept_ar : "  << param.intercept_ar << ", "
			   << "intercept_al : "  << param.intercept_al << "\n"
			   << "slope_ax : "      << param.slope_ax << ", "
			   << "slope_ay : "      << param.slope_ay << ", "
			   << "slope_ar : "      << param.slope_ar << ", "
			   << "slope_al : "      << param.slope_al << "\n"
			   << "slope2_ax : "     << param.slope2_ax << ", "
			   << "slope2_ay : "     << param.slope2_ay << ", "
			   << "slope2_ar : "     << param.slope2_ar << ", "
			   << "slope2_al : "     << param.slope2_al << "\n"
			   << "intercept_px : "  << param.intercept_px << ", "
			   << "intercept_py : "  << param.intercept_py << ", "
			   << "intercept_pr : "  << param.intercept_pr << ", "
			   << "intercept_pl : "  << param.intercept_pl << "\n"
			   << "slope_px : "      << param.slope_px << ", "
			   << "slope_py : "      << param.slope_py << ", "
			   << "slope_pr : "      << param.slope_pr << ", "
			   << "slope_pl : "      << param.slope_pl << "\n"
			   << "slope2_px : "     << param.slope2_px << ", "
			   << "slope2_py : "     << param.slope2_py << ", "
			   << "slope2_pr : "     << param.slope2_pr << ", "
			   << "slope2_pl : "     << param.slope2_pl;

  return;

}


void ConnectionData::ReadFiducialAreaEcc(int ecc, 
					 const fs::path &file_dir_path) {
  if ( ecc < 0 || ecc > 8 )
    throw std::invalid_argument("ECC id not valid : " + ecc);
  const std::string fiducial_file_path = (file_dir_path/ECC_FIDUCIAL_FILENAMES[ecc]).string();
  ReadFiducialArea(ecc_fiducial_[ecc], fiducial_file_path);
}

void ConnectionData::ReadFiducialArea(std::map<int, std::vector<FiducialArea > > &fiducial_area_map,
				      const std::string file_dir_path) {

  std::ifstream ifs(file_dir_path);
  BOOST_LOG_TRIVIAL(trace) << "Fiducial file name : " << file_dir_path;
  std::multimap<int, FiducialArea > fa_multi;
  FiducialArea fa;
  while ( ifs >> fa.pl
	  >> fa.p[0].x >> fa.p[0].y >> fa.p[0].z
	  >> fa.p[1].x >> fa.p[1].y >> fa.p[1].z ) {
    fa.p[0].y += 20.e3; fa.p[1].y += 20.e3; // vertical offset
    fa_multi.insert(std::make_pair(fa.pl, fa));
  }

  for ( auto itr = fa_multi.begin(); itr != fa_multi.end(); itr++ ) {
    auto range = fa_multi.equal_range(itr->first);
    std::vector<FiducialArea > vec;
    for ( auto itr2 = range.first; itr2 != range.second; itr2++ ) {
      vec.push_back(itr2->second);
    }
    fiducial_area_map.insert(std::make_pair(itr->first, vec));
  }

  ifs.close();


  for ( auto itr = fiducial_area_map.begin(); itr != fiducial_area_map.end(); itr++ ) {
    auto range = fiducial_area_map.equal_range(itr->first);
    for ( auto itr2 = range.first; itr2 != range.second; itr2++ ) {
      for ( auto fa2 : itr2->second )
	BOOST_LOG_TRIVIAL(trace) << "pl : " << fa2.pl << ", "
				 << "x0 : " << fa2.p[0].x << ", "
				 << "y0 : " << fa2.p[0].y << ", "
				 << "z0 : " << fa2.p[0].z << ", "
				 << "x1 : " << fa2.p[1].x << ", "
				 << "y1 : " << fa2.p[1].y << ", "
				 << "z1 : " << fa2.p[1].z;
    }
  }

}

void ConnectionData::ReadEfficiencyEcc(int ecc, 
				       const fs::path &file_dir_path) {
  if ( ecc < 0 || ecc > 8 )
    throw std::invalid_argument("ECC id not valid : " + ecc);
  const std::string efficiency_file_path = (file_dir_path/ECC_EFFICIENCY_FILENAMES[ecc]).string();
  ReadEfficiency(ecc_efficiency_[ecc], efficiency_file_path);
}

void ConnectionData::ReadEfficiency(std::map<int, std::vector<Efficiency > > &efficiency_map,
				    const std::string file_dir_path) {

  std::ifstream ifs(file_dir_path);
  BOOST_LOG_TRIVIAL(trace) << "Efficiency file name : " << file_dir_path;
  Efficiency eff;
  int tmp;
  std::vector<Efficiency > eff_vec;
  eff_vec.reserve(EFFICIENCY_VECTOR_SIZE);
  int i_eff_vec = 0;
  while ( ifs >> eff.pl >> eff.range[0] >> eff.range[1]
	  >> tmp >> tmp >> eff.efficiency >> eff.efficiency_err ) {
    eff_vec.push_back(eff);
    i_eff_vec++;
    if ( i_eff_vec == EFFICIENCY_VECTOR_SIZE ) {
      efficiency_map.insert(std::make_pair(eff.pl, eff_vec));
      eff_vec.clear(); eff_vec.shrink_to_fit();
      eff_vec.reserve(EFFICIENCY_VECTOR_SIZE);
      i_eff_vec = 0;
    }
  }

  ifs.close();

  for ( auto itr = efficiency_map.begin(); itr != efficiency_map.end(); itr++ ) {
    auto range = efficiency_map.equal_range(itr->first);
    for ( auto itr2 = range.first; itr2 != range.second; itr2++ ) {
      for ( auto eff2 : itr2->second)
	BOOST_LOG_TRIVIAL(trace) << "pl : " << eff2.pl << ", "
				 << "Angle : (" << eff2.range[0] << ", " << eff2.range[1] << "), "
				 << "Efficiency : " << eff2.efficiency << ", "
				 << " +/- " << eff2.efficiency_err;
    }
  }

}

void ConnectionData::GetETECCConnectData(t2l_param &et_ecc_param) const {
  et_ecc_param = et_ecc_param_;
  return;
}

void ConnectionData::GetETECCFeConnectData(t2l_param &et_ecc_fe_param) const {
  et_ecc_fe_param = et_ecc_fe_param_;
  return;
}

void ConnectionData::GetETECCFeFeConnectData(t2l_param &et_ecc_fe_fe_param) const {
  et_ecc_fe_fe_param = et_ecc_fe_fe_param_;
  return;
}

void ConnectionData::GetBlackFeConnectData(t2l_param &black_fe_param) const {
  black_fe_param = black_fe_param_;
  return;
}

void ConnectionData::GetBlackReFeConnectData(t2l_param &black_re_fe_param) const {
  black_re_fe_param = black_re_fe_param_;
  return;
}

void ConnectionData::GetBlackReFeFeConnectData(t2l_param &black_re_fe_fe_param) const {
  black_re_fe_fe_param = black_re_fe_fe_param_;
  return;
}

void ConnectionData::GetBlackReFeWaterConnectData(t2l_param &black_re_fe_water_param) const {
  black_re_fe_water_param = black_re_fe_water_param_;
  return;
}

void ConnectionData::GetBlackWaterConnectData(t2l_param &black_water_param) const {
  black_water_param = black_water_param_;
  return;
}

void ConnectionData::GetFeConnectData(t2l_param &fe_param) const {
  fe_param = fe_param_;
  return;
}

void ConnectionData::GetFeFeConnectData(t2l_param &fe_fe_param) const {
  fe_fe_param = fe_fe_param_;
  return;
}

void ConnectionData::GetFeFeFeConnectData(t2l_param &fe_fe_fe_param) const {
  fe_fe_fe_param = fe_fe_fe_param_;
  return;
}

void ConnectionData::GetFeWaterConnectData(t2l_param &fe_water_param) const {
  fe_water_param = fe_water_param_;
  return;
}

void ConnectionData::GetFeWaterFeConnectData(t2l_param &fe_water_fe_param) const {
  fe_water_fe_param = fe_water_fe_param_;
  return;
}

void ConnectionData::GetReETECCFeFeFeConnectData(t2l_param &re_et_ecc_fe_fe_fe_param) const {
  re_et_ecc_fe_fe_fe_param = re_et_ecc_fe_fe_fe_param_;
  return;
}

void ConnectionData::GetReETECCFeFeFeFeConnectData(t2l_param &re_et_ecc_fe_fe_fe_fe_param) const {
  re_et_ecc_fe_fe_fe_fe_param = re_et_ecc_fe_fe_fe_fe_param_;
  return;
}

void ConnectionData::GetReETECCFeFeFeWaterConnectData(t2l_param &re_et_ecc_fe_fe_fe_water_param) const {
  re_et_ecc_fe_fe_fe_water_param = re_et_ecc_fe_fe_fe_water_param_;
  return;
}

void ConnectionData::GetReETECCFeFeWaterConnectData(t2l_param &re_et_ecc_fe_fe_water_param) const {
  re_et_ecc_fe_fe_water_param = re_et_ecc_fe_fe_water_param_;
  return;
}

void ConnectionData::GetReETECCFeWaterConnectData(t2l_param &re_et_ecc_fe_water_param) const {
  re_et_ecc_fe_water_param = re_et_ecc_fe_water_param_;
  return;
}

void ConnectionData::GetReFeConnectData(t2l_param &re_fe_param) const {
  re_fe_param = re_fe_param_;
  return;
}

void ConnectionData::GetReFeFeConnectData(t2l_param &re_fe_fe_param) const {
  re_fe_fe_param = re_fe_fe_param_;
  return;
}

void ConnectionData::GetReFeFeFeConnectData(t2l_param &re_fe_fe_fe_param) const {
  re_fe_fe_fe_param = re_fe_fe_fe_param_;
  return;
}

void ConnectionData::GetReFeFeFeFeConnectData(t2l_param &re_fe_fe_fe_fe_param) const {
  re_fe_fe_fe_fe_param = re_fe_fe_fe_fe_param_;
  return;
}

void ConnectionData::GetReFeFeFeFeFeConnectData(t2l_param &re_fe_fe_fe_fe_fe_param) const {
  re_fe_fe_fe_fe_fe_param = re_fe_fe_fe_fe_fe_param_;
  return;
}

void ConnectionData::GetReFeWaterConnectData(t2l_param &re_fe_water_param) const {
  re_fe_water_param = re_fe_water_param_;
  return;
}

void ConnectionData::GetReFeWaterFeConnectData(t2l_param &re_fe_water_fe_param) const {
  re_fe_water_fe_param = re_fe_water_fe_param_;
  return;
}

void ConnectionData::GetReFeWaterFeWaterConnectData(t2l_param &re_fe_water_fe_water_param) const {
  re_fe_water_fe_water_param = re_fe_water_fe_water_param_;
  return;
}

void ConnectionData::GetReFeWaterFeWaterFeConnectData(t2l_param &re_fe_water_fe_water_fe_param) const {
  re_fe_water_fe_water_fe_param = re_fe_water_fe_water_fe_param_;
  return;
}

void ConnectionData::GetReWaterFeWaterConnectData(t2l_param &re_water_fe_water_param) const {
  re_water_fe_water_param = re_fe_water_fe_water_param_;
  return;
}

void ConnectionData::GetWaterConnectData(t2l_param &water_param) const {
  water_param = water_param_;
  return;
}

void ConnectionData::GetFiducialAreaData(int ecc, std::map<int, std::vector<FiducialArea > > &ecc_fiducial) const {
  if ( ecc < 0 || ecc > 8 )
    throw std::invalid_argument("Ecc id not valid : " + ecc);
  ecc_fiducial = ecc_fiducial_[ecc];
  return;
}

void ConnectionData::GetEfficiencyData(int ecc, std::map<int, std::vector<Efficiency > > &ecc_efficiency) const {
  if ( ecc < 0 || ecc > 8 )
    throw std::invalid_argument("Ecc id not valid : " + ecc);
  ecc_efficiency = ecc_efficiency_[ecc];
  return;
}

void ConnectionData::DrawFiducialAreaData(int ecc, int plate,
					  std::vector<B2EmulsionSummary* > &emulsions,
					  int eventid) const {

  std::vector<Double_t > edge_x, edge_y; // fiducial 外周の座標
  std::vector<Double_t > point_x, point_y; // emulsion の座標
  std::vector<FiducialArea > fa = ecc_fiducial_[ecc].at(plate);
  for ( auto fa_points : fa ) {
    edge_x.push_back(fa_points.p[0].x);
    edge_y.push_back(fa_points.p[0].y);
  }

  TGraph *edges = new TGraph(edge_x.size(),
			     &edge_x[0],
			     &edge_y[0]);

  for ( auto emulsion : emulsions ) {
    if ( emulsion->GetEcc() != ecc ) continue;
    if ( emulsion->GetPlate() + 1 != plate ) continue;
    point_x.push_back(emulsion->GetFilmPosition().GetValue().X() * 1000.);
    point_y.push_back(emulsion->GetFilmPosition().GetValue().Y() * 1000.);
  }

  TCanvas *c = new TCanvas("c", "c");
  edges->SetTitle(Form("Fiducial area (ECC%d, PL%d);x[um];y[um]", ecc+1, plate));
  edges->Draw("AP");

  if ( !point_x.empty() ) {
    TGraph *points = new TGraph(point_x.size(),
				&point_x[0],
				&point_y[0]);
    points->SetMarkerColor(kRed);
    points->SetMarkerStyle(kStar);
    points->Draw("SAME P");
  }

  c->SaveAs(Form("/home/t2k/odagawa/data/plots/fiducial_test/fiducial_ecc%d_plate%d_event%d.pdf", ecc+1, plate, eventid));
}

void ConnectionData::DrawEfficiencyData(int ecc, int plate) const {
  std::vector<double > ang_value, ang_error, eff_value, eff_error;
  std::vector<Efficiency > eff = ecc_efficiency_[ecc].at(plate);
  for ( auto ieff : eff ) {
    double iang_value = (ieff.range[0] + ieff.range[1]) / 2.;
    double iang_error = (ieff.range[1] - ieff.range[0]) / 2.;
    ang_value.push_back(iang_value);
    ang_error.push_back(iang_error);
    eff_value.push_back(ieff.efficiency);    
    eff_error.push_back(ieff.efficiency_err);
  }

  TGraphErrors *g_eff = new TGraphErrors(ang_value.size(),
					 &ang_value[0],
					 &eff_value[0],
					 &ang_error[0],
					 &eff_error[0]);
  TCanvas *c = new TCanvas("c_eff", "c_eff");
  c->DrawFrame(0., 0., 250.e3, 250.e3);
  g_eff->SetTitle(Form("Efficiency (ECC%d, PL%d);tan#theta;Efficiency", ecc+1, plate));
  g_eff->Draw("AP");

  c->SaveAs(Form("/home/t2k/odagawa/data/plots/efficiency_test/efficiency_ecc%d_plate%d.pdf", ecc+1, plate));

}
