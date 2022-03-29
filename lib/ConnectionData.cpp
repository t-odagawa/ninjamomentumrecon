#include "ConnectionData.hpp"

#include <boost/log/core.hpp>
#include <boost/log/trivial.hpp>
#include <boost/log/expressions.hpp>
#include <boost/filesystem.hpp>

#include <string>

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

ConnectionData::ConnectionData(const std::string &file_dir_path) {
  fs::path file_dir(file_dir_path);

  ReadConnectData(file_dir);

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

  return;

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
