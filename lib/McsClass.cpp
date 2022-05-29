#include <fstream>
#include <iostream>
#include <iomanip>
#include <string>
#include <vector>
#include <map>
#include <cmath>

#include "McsClass.hpp"

namespace Momentum_recon
{

  void WriteEventInformationHeader(std::ofstream &ofs, Event_information &ev) {
    int num_chain = ev.chains.size();
    int num_true_chain = ev.true_chains.size();

    ofs.write((char*)& ev.groupid, sizeof(int));
    ofs.write((char*)& ev.unixtime, sizeof(int));
    ofs.write((char*)& ev.tracker_track_id, sizeof(int));
    ofs.write((char*)& ev.entry_in_daily_file, sizeof(int));
    ofs.write((char*)& ev.vertex_pl, sizeof(int));
    ofs.write((char*)& ev.ecc_id, sizeof(int));
    ofs.write((char*)& ev.vertex_material, sizeof(int));
    ofs.write((char*)& ev.recon_vertex_position, sizeof(double) * 3);
    ofs.write((char*)& ev.true_vertex_position, sizeof(double) * 3);
    ofs.write((char*)& ev.weight, sizeof(double));
    ofs.write((char*)& ev.nu_energy, sizeof(double));
    ofs.write((char*)& ev.nu_ax, sizeof(double));
    ofs.write((char*)& ev.nu_ay, sizeof(double));
    ofs.write((char*)& num_chain, sizeof(int));
    ofs.write((char*)& num_true_chain, sizeof(int));

  }

  bool ReadEventInformationHeader(std::ifstream &ifs, Event_information &ev, int &num_chain , int &num_true_chain) {
    if (!ifs.read((char*)& ev.groupid, sizeof(int))) return false;
    if (!ifs.read((char*)& ev.unixtime, sizeof(int))) return false;
    if (!ifs.read((char*)& ev.tracker_track_id, sizeof(int))) return false;
    if (!ifs.read((char*)& ev.entry_in_daily_file, sizeof(int))) return false;
    if (!ifs.read((char*)& ev.vertex_pl, sizeof(int))) return false;
    if (!ifs.read((char*)& ev.ecc_id, sizeof(int))) return false;
    if (!ifs.read((char*)& ev.vertex_material, sizeof(int))) return false;
    if (!ifs.read((char*)& ev.recon_vertex_position, sizeof(double) * 3)) return false;
    if (!ifs.read((char*)& ev.true_vertex_position, sizeof(double) * 3)) return false;
    if (!ifs.read((char*)& ev.weight, sizeof(double))) return false;
    if (!ifs.read((char*)& ev.nu_energy, sizeof(double))) return false;
    if (!ifs.read((char*)& ev.nu_ax, sizeof(double))) return false;
    if (!ifs.read((char*)& ev.nu_ay, sizeof(double))) return false;
    if (!ifs.read((char*)& num_chain, sizeof(int))) return false;
    if (!ifs.read((char*)& num_true_chain, sizeof(int))) return false;
    return true;
  }

  void WriteMomChainHeader(std::ofstream &ofs, Mom_chain &mom_chain) {
    int num_base = mom_chain.base.size();
    int num_link = mom_chain.base_pair.size();
    
    ofs.write((char*)& mom_chain.chainid, sizeof(int));
    ofs.write((char*)& mom_chain.stop_flag, sizeof(int));
    ofs.write((char*)& mom_chain.particle_flag, sizeof(int));
    ofs.write((char*)& mom_chain.direction, sizeof(int));
    ofs.write((char*)& mom_chain.charge_sign, sizeof(int));
    ofs.write((char*)& mom_chain.ecc_range_mom, sizeof(double) * 2);
    ofs.write((char*)& mom_chain.ecc_mcs_mom, sizeof(double) * 2);
    ofs.write((char*)& mom_chain.bm_range_mom, sizeof(double));
    ofs.write((char*)& mom_chain.bm_curvature_mom, sizeof(double));
    ofs.write((char*)& mom_chain.ecc_range_mom_error, sizeof(double) * 2 * 2);
    ofs.write((char*)& mom_chain.ecc_mcs_mom_error, sizeof(double) * 2 * 2);
    ofs.write((char*)& mom_chain.bm_range_mom_error, sizeof(double) * 2);
    ofs.write((char*)& mom_chain.bm_curvature_mom_error, sizeof(double) * 2);
    ofs.write((char*)& mom_chain.muon_likelihood, sizeof(double));
    ofs.write((char*)& mom_chain.proton_likelihood, sizeof(double));
    ofs.write((char*)& num_base, sizeof(int));
    ofs.write((char*)& num_link, sizeof(int));
  }

  bool ReadMomChainHeader(std::ifstream &ifs, Mom_chain &mom_chain, int &num_base, int &num_link) {

    if (!ifs.read((char*)& mom_chain.chainid, sizeof(int))) return false;
    if (!ifs.read((char*)& mom_chain.stop_flag, sizeof(int))) return false;
    if (!ifs.read((char*)& mom_chain.particle_flag, sizeof(int))) return false;
    if (!ifs.read((char*)& mom_chain.direction, sizeof(int))) return false;
    if (!ifs.read((char*)& mom_chain.charge_sign, sizeof(int))) return false;
    if (!ifs.read((char*)& mom_chain.ecc_range_mom, sizeof(double) * 2)) return false;
    if (!ifs.read((char*)& mom_chain.ecc_mcs_mom, sizeof(double) * 2)) return false;
    if (!ifs.read((char*)& mom_chain.bm_range_mom, sizeof(double))) return false;
    if (!ifs.read((char*)& mom_chain.bm_curvature_mom, sizeof(double))) return false;
    if (!ifs.read((char*)& mom_chain.ecc_range_mom_error, sizeof(double) * 2 * 2)) return false;
    if (!ifs.read((char*)& mom_chain.ecc_mcs_mom_error, sizeof(double) * 2 * 2)) return false;
    if (!ifs.read((char*)& mom_chain.bm_range_mom_error, sizeof(double) * 2)) return false;
    if (!ifs.read((char*)& mom_chain.bm_curvature_mom_error, sizeof(double) * 2)) return false;
    if (!ifs.read((char*)& mom_chain.muon_likelihood, sizeof(double))) return false;
    if (!ifs.read((char*)& mom_chain.proton_likelihood, sizeof(double))) return false;
    if (!ifs.read((char*)& num_base, sizeof(int))) return false;
    if (!ifs.read((char*)& num_link, sizeof(int))) return false;
    return true;
    
  }

  void WriteEventInformationBin(std::string filename, std::vector<Event_information > &ev_vector) {
    std::ofstream ofs(filename, std::ios::binary);
    if ( !ofs ) {
      std::cerr << "File : " << filename << " not found!" << std::endl;
      std::exit(1);
    }
    else if ( ev_vector.empty() ) {
      std::cerr << "Target event data is null" << std::endl;
      std::exit(1);
    }
    else {
      for ( auto ev : ev_vector ) {
	WriteEventInformationHeader(ofs, ev);
	for ( auto mom_chain : ev.chains ) {
	  WriteMomChainHeader(ofs, mom_chain);
	  for ( auto base : mom_chain.base ) {
	    ofs.write((char*)& base, sizeof(Mom_basetrack));
	  } // base
	  for ( auto link : mom_chain.base_pair ) {
	    ofs.write((char*)& link.first, sizeof(Mom_basetrack));
	    ofs.write((char*)& link.second, sizeof(Mom_basetrack));
	  } // link
	} // mom_chain
	
	for ( auto mom_true_chain : ev.true_chains ) {
	  WriteMomChainHeader(ofs, mom_true_chain);
	  for ( auto true_base : mom_true_chain.base ) {
	    ofs.write((char*)& true_base, sizeof(Mom_basetrack));
	  } // true_base
	  for ( auto true_link : mom_true_chain.base_pair ) {
	    ofs.write((char*)& true_link.first, sizeof(Mom_basetrack));
	    ofs.write((char*)& true_link.second, sizeof(Mom_basetrack));
	  } // true_link		  
	} // true_mom_chain

      } // ev
    }
  }

  std::vector<Event_information > ReadEventInformationBin(std::string filename) {

    std::vector<Event_information > ret;
    std::ifstream ifs(filename, std::ios::binary);
    Event_information ev;
    Mom_chain mom_chain;
    Mom_basetrack basetrack;
    std::pair<Mom_basetrack, Mom_basetrack> link;
    int num_chain = 0;
    int num_true_chain = 0;
    int num_base = 0;
    int num_link = 0;

    while ( ReadEventInformationHeader(ifs, ev, num_chain, num_true_chain) ) {
      // chains
      for ( int ichain = 0; ichain < num_chain; ichain++ ) {
	if ( ReadMomChainHeader(ifs, mom_chain, num_base, num_link) ) {
	  mom_chain.base.clear();
	  mom_chain.base_pair.clear();
	  mom_chain.base.reserve(num_base);
	  mom_chain.base_pair.reserve(num_link);
	  for ( int j = 0; j < num_base; j++ ) {
	    ifs.read((char*)& basetrack, sizeof(Mom_basetrack));
	    mom_chain.base.push_back(basetrack);
	  }
	  for ( int j = 0; j < num_link; j++ ) {
	    ifs.read((char*)& link.first, sizeof(Mom_basetrack));
	    ifs.read((char*)& link.second, sizeof(Mom_basetrack));
	    mom_chain.base_pair.push_back(link);
	  }
	  ev.chains.push_back(mom_chain);
	}
      }

      // true chains

      for ( int itrue_chain = 0; itrue_chain < num_true_chain; itrue_chain++ ) {
	if ( ReadMomChainHeader(ifs, mom_chain, num_base, num_link) ) {
	  mom_chain.base.clear();
	  mom_chain.base_pair.clear();
	  mom_chain.base.reserve(num_base);
	  mom_chain.base_pair.reserve(num_link);
	  for ( int j = 0; j < num_base; j++ ) {
	    ifs.read((char*)& basetrack, sizeof(Mom_basetrack));
	    mom_chain.base.push_back(basetrack);
	  }
	  for ( int j = 0; j < num_link; j++ ) {
	    ifs.read((char*)& link.first, sizeof(Mom_basetrack));
	    ifs.read((char*)& link.second, sizeof(Mom_basetrack));
	    mom_chain.base_pair.push_back(link);
	  }
	  ev.true_chains.push_back(mom_chain);
	}
      }

      ret.push_back(ev);
      ev.chains.clear();
      ev.true_chains.clear();
    }

    return ret;

  }

  void WriteEventInformationTxt(std::string filename, std::vector<Event_information > &ev_vector) {
    std::ofstream ofs(filename);

    for ( auto &ev : ev_vector ) {
      //event inforamtion header 書き出し
      ofs << std::right << std::fixed
	  << std::setw(10) << std::setprecision(0) << ev.groupid << " "
	  << std::setw(10) << std::setprecision(0) << ev.unixtime << " "
	  << std::setw(3)  << std::setprecision(0) << ev.tracker_track_id << " "
	  << std::setw(5)  << std::setprecision(0) << ev.entry_in_daily_file << " "
	  << std::setw(3)  << std::setprecision(0) << ev.vertex_pl << " "
	  << std::setw(2)  << std::setprecision(0) << ev.ecc_id << " "
	  << std::setw(2)  << std::setprecision(0) << ev.vertex_material << " "
	  << std::setw(8)  << std::setprecision(1) << ev.recon_vertex_position[0] << " "
	  << std::setw(8)  << std::setprecision(1) << ev.recon_vertex_position[1] << " "
	  << std::setw(8)  << std::setprecision(1) << ev.recon_vertex_position[2] << " "
	  << std::setw(8)  << std::setprecision(1) << ev.true_vertex_position[0] << " "
	  << std::setw(8)  << std::setprecision(1) << ev.true_vertex_position[1] << " "
	  << std::setw(8)  << std::setprecision(1) << ev.true_vertex_position[2] << " "
	  << std::setw(7)  << std::setprecision(2) << ev.weight << " "
	  << std::setw(7)  << std::setprecision(1) << ev.nu_energy << " "
	  << std::setw(7)  << std::setprecision(4) << ev.nu_ax << " "
	  << std::setw(7)  << std::setprecision(4) << ev.nu_ay << " "
	  << std::setw(5)  << std::setprecision(0) << ev.chains.size() << " "
	  << std::setw(5)  << std::setprecision(0) << ev.true_chains.size() << std::endl;
      for ( int i = 0; i < ev.chains.size(); i++ ) {
	ofs << std::right << std::fixed
	    << std::setw(10) << std::setprecision(0) << ev.chains[i].chainid << " "
	    << std::setw(5)  << std::setprecision(0) << ev.chains[i].stop_flag << " "
	    << std::setw(5)  << std::setprecision(0) << ev.chains[i].particle_flag << " "
	    << std::setw(3)  << std::setprecision(0) << ev.chains[i].direction << " "
	    << std::setw(2)  << std::setprecision(0) << ev.chains[i].charge_sign << " "
	    << std::setw(7)  << std::setprecision(1) << ev.chains[i].ecc_range_mom[0] << " "
	    << std::setw(7)  << std::setprecision(1) << ev.chains[i].ecc_range_mom[1] << " "
	    << std::setw(7)  << std::setprecision(1) << ev.chains[i].ecc_mcs_mom[0] << " "
	    << std::setw(7)  << std::setprecision(1) << ev.chains[i].ecc_mcs_mom[1] << " "
	    << std::setw(7)  << std::setprecision(1) << ev.chains[i].bm_range_mom << " "
	    << std::setw(7)  << std::setprecision(1) << ev.chains[i].bm_curvature_mom << " "
	    << std::setw(7)  << std::setprecision(1) << ev.chains[i].ecc_range_mom_error[0][0] << " "
	    << std::setw(7)  << std::setprecision(1) << ev.chains[i].ecc_range_mom_error[0][1] << " "
	    << std::setw(7)  << std::setprecision(1) << ev.chains[i].ecc_range_mom_error[1][0] << " "
	    << std::setw(7)  << std::setprecision(1) << ev.chains[i].ecc_range_mom_error[1][1] << " "
	    << std::setw(7)  << std::setprecision(1) << ev.chains[i].ecc_mcs_mom_error[0][0] << " "
	    << std::setw(7)  << std::setprecision(1) << ev.chains[i].ecc_mcs_mom_error[0][1] << " "
	    << std::setw(7)  << std::setprecision(1) << ev.chains[i].ecc_mcs_mom_error[1][0] << " "
	    << std::setw(7)  << std::setprecision(1) << ev.chains[i].ecc_mcs_mom_error[1][1] << " "
	    << std::setw(7)  << std::setprecision(1) << ev.chains[i].bm_range_mom_error[0] << " "
	    << std::setw(7)  << std::setprecision(1) << ev.chains[i].bm_range_mom_error[1] << " "
	    << std::setw(7)  << std::setprecision(1) << ev.chains[i].bm_curvature_mom_error[0] << " "
	    << std::setw(7)  << std::setprecision(1) << ev.chains[i].bm_curvature_mom_error[1] << " "
	    << std::setw(7)  << std::setprecision(4) << ev.chains[i].muon_likelihood << " "
	    << std::setw(7)  << std::setprecision(4) << ev.chains[i].proton_likelihood << " "
	    << std::setw(3)  << std::setprecision(0) << ev.chains[i].base.size() << " "
	    << std::setw(3)  << std::setprecision(0) << ev.chains[i].base_pair.size() << std::endl;

	for ( auto &b : ev.chains[i].base ) {
	  ofs << std::right << std::fixed
	      << std::setw(4)  << std::setprecision(0) << b.pl << " "
	      << std::setw(10) << std::setprecision(0) << b.rawid << " "
	      << std::setw(7)  << std::setprecision(4) << b.ax << " "
	      << std::setw(7)  << std::setprecision(4) << b.ay << " "
	      << std::setw(8)  << std::setprecision(1) << b.x << " "
	      << std::setw(8)  << std::setprecision(1) << b.y << " "
	      << std::setw(8)  << std::setprecision(1) << b.z << " "
	      << std::setw(7)  << std::setprecision(0) << b.m[0].zone << " "
	      << std::setw(4)  << std::setprecision(0) << b.m[0].view << " "
	      << std::setw(3)  << std::setprecision(0) << b.m[0].imager << " "
	      << std::setw(5)  << std::setprecision(0) << b.m[0].ph << " "
	      << std::setw(5)  << std::setprecision(0) << b.m[0].pixelnum << " "
	      << std::setw(5)  << std::setprecision(0) << b.m[0].hitnum << " "
	      << std::setw(7)  << std::setprecision(0) << b.m[1].zone << " "
	      << std::setw(4)  << std::setprecision(0) << b.m[1].view << " "
	      << std::setw(3)  << std::setprecision(0) << b.m[1].imager << " "
	      << std::setw(5)  << std::setprecision(0) << b.m[1].ph << " "
	      << std::setw(5)  << std::setprecision(0) << b.m[1].pixelnum << " "
	      << std::setw(5)  << std::setprecision(0) << b.m[1].hitnum << std::endl;
	}
	for ( auto &p : ev.chains[i].base_pair ) {
	  ofs << std::right << std::fixed
	      << std::setw(4)  << std::setprecision(0) << p.first.pl << " "
	      << std::setw(10) << std::setprecision(0) << p.first.rawid << " "
	      << std::setw(4)  << std::setprecision(0) << p.second.pl << " "
	      << std::setw(10) << std::setprecision(0) << p.second.rawid << " "
	      << std::setw(7)  << std::setprecision(4) << p.first.ax << " "
	      << std::setw(7)  << std::setprecision(4) << p.first.ay << " "
	      << std::setw(8)  << std::setprecision(1) << p.first.x << " "
	      << std::setw(8)  << std::setprecision(1) << p.first.y << " "
	      << std::setw(8)  << std::setprecision(1) << p.first.z << " "
	      << std::setw(7)  << std::setprecision(4) << p.second.ax << " "
	      << std::setw(7)  << std::setprecision(4) << p.second.ay << " "
	      << std::setw(8)  << std::setprecision(1) << p.second.x << " "
	      << std::setw(8)  << std::setprecision(1) << p.second.y << " "
	      << std::setw(8)  << std::setprecision(1) << p.second.z << std::endl;
	}

      }

      for ( int i = 0; i < ev.true_chains.size(); i++ ) {
	ofs << std::right << std::fixed
	    << std::setw(10) << std::setprecision(0) << ev.true_chains[i].chainid << " "
	    << std::setw(5)  << std::setprecision(0) << ev.true_chains[i].stop_flag << " "
	    << std::setw(5)  << std::setprecision(0) << ev.true_chains[i].particle_flag << " "
	    << std::setw(3)  << std::setprecision(0) << ev.true_chains[i].direction << " "
	    << std::setw(2)  << std::setprecision(0) << ev.true_chains[i].charge_sign << " "
	    << std::setw(7)  << std::setprecision(1) << ev.true_chains[i].ecc_range_mom[0] << " "
	    << std::setw(7)  << std::setprecision(1) << ev.true_chains[i].ecc_range_mom[1] << " "
	    << std::setw(7)  << std::setprecision(1) << ev.true_chains[i].ecc_mcs_mom[0] << " "
	    << std::setw(7)  << std::setprecision(1) << ev.true_chains[i].ecc_mcs_mom[1] << " "
	    << std::setw(7)  << std::setprecision(1) << ev.true_chains[i].bm_range_mom << " "
	    << std::setw(7)  << std::setprecision(1) << ev.true_chains[i].bm_curvature_mom << " "
	    << std::setw(7)  << std::setprecision(1) << ev.true_chains[i].ecc_range_mom_error[0][0] << " "
	    << std::setw(7)  << std::setprecision(1) << ev.true_chains[i].ecc_range_mom_error[0][1] << " "
	    << std::setw(7)  << std::setprecision(1) << ev.true_chains[i].ecc_range_mom_error[1][0] << " "
	    << std::setw(7)  << std::setprecision(1) << ev.true_chains[i].ecc_range_mom_error[1][1] << " "
	    << std::setw(7)  << std::setprecision(1) << ev.true_chains[i].ecc_mcs_mom_error[0][0] << " "
	    << std::setw(7)  << std::setprecision(1) << ev.true_chains[i].ecc_mcs_mom_error[0][1] << " "
	    << std::setw(7)  << std::setprecision(1) << ev.true_chains[i].ecc_mcs_mom_error[1][0] << " "
	    << std::setw(7)  << std::setprecision(1) << ev.true_chains[i].ecc_mcs_mom_error[1][1] << " "
	    << std::setw(7)  << std::setprecision(1) << ev.true_chains[i].bm_range_mom_error[0] << " "
	    << std::setw(7)  << std::setprecision(1) << ev.true_chains[i].bm_range_mom_error[1] << " "
	    << std::setw(7)  << std::setprecision(1) << ev.true_chains[i].bm_curvature_mom_error[0] << " "
	    << std::setw(7)  << std::setprecision(1) << ev.true_chains[i].bm_curvature_mom_error[1] << " "
	    << std::setw(7)  << std::setprecision(4) << ev.true_chains[i].muon_likelihood << " "
	    << std::setw(7)  << std::setprecision(4) << ev.true_chains[i].proton_likelihood << " "
	    << std::setw(3)  << std::setprecision(0) << ev.true_chains[i].base.size() << " "
	    << std::setw(3)  << std::setprecision(0) << ev.true_chains[i].base_pair.size() << std::endl;

	for ( auto &b : ev.true_chains[i].base ) {
	  ofs << std::right << std::fixed
	      << std::setw(4)  << std::setprecision(0) << b.pl << " "
	      << std::setw(10) << std::setprecision(0) << b.rawid << " "
	      << std::setw(7)  << std::setprecision(4) << b.ax << " "
	      << std::setw(7)  << std::setprecision(4) << b.ay << " "
	      << std::setw(8)  << std::setprecision(1) << b.x << " "
	      << std::setw(8)  << std::setprecision(1) << b.y << " "
	      << std::setw(8)  << std::setprecision(1) << b.z << " "
	      << std::setw(7)  << std::setprecision(0) << b.m[0].zone << " "
	      << std::setw(4)  << std::setprecision(0) << b.m[0].view << " "
	      << std::setw(3)  << std::setprecision(0) << b.m[0].imager << " "
	      << std::setw(5)  << std::setprecision(0) << b.m[0].ph << " "
	      << std::setw(5)  << std::setprecision(0) << b.m[0].pixelnum << " "
	      << std::setw(5)  << std::setprecision(0) << b.m[0].hitnum << " "
	      << std::setw(7)  << std::setprecision(0) << b.m[1].zone << " "
	      << std::setw(4)  << std::setprecision(0) << b.m[1].view << " "
	      << std::setw(3)  << std::setprecision(0) << b.m[1].imager << " "
	      << std::setw(5)  << std::setprecision(0) << b.m[1].ph << " "
	      << std::setw(5)  << std::setprecision(0) << b.m[1].pixelnum << " "
	      << std::setw(5)  << std::setprecision(0) << b.m[1].hitnum << std::endl;
	}
	for ( auto &p : ev.true_chains[i].base_pair ) {
	  ofs << std::right << std::fixed
	      << std::setw(4)  << std::setprecision(0) << p.first.pl << " "
	      << std::setw(10) << std::setprecision(0) << p.first.rawid << " "
	      << std::setw(4)  << std::setprecision(0) << p.second.pl << " "
	      << std::setw(10) << std::setprecision(0) << p.second.rawid << " "
	      << std::setw(7)  << std::setprecision(4) << p.first.ax << " "
	      << std::setw(7)  << std::setprecision(4) << p.first.ay << " "
	      << std::setw(8)  << std::setprecision(1) << p.first.x << " "
	      << std::setw(8)  << std::setprecision(1) << p.first.y << " "
	      << std::setw(8)  << std::setprecision(1) << p.first.z << " "
	      << std::setw(7)  << std::setprecision(4) << p.second.ax << " "
	      << std::setw(7)  << std::setprecision(4) << p.second.ay << " "
	      << std::setw(8)  << std::setprecision(1) << p.second.x << " "
	      << std::setw(8)  << std::setprecision(1) << p.second.y << " "
	      << std::setw(8)  << std::setprecision(1) << p.second.z << std::endl;
	}

      }

    }

  }

  std::vector<Event_information > ReadEventInformationTxt(std::string filename){
    std::vector<Event_information > ret;

    std::ifstream ifs(filename);
    Event_information ev;
    Mom_chain t;
    Mom_basetrack m;
    std::pair<Mom_basetrack, Mom_basetrack> p;
    int chain_num, true_chain_num;
    int base_num, base_pair_num;

    while ( ifs >> ev.groupid >> ev.unixtime >> ev.tracker_track_id >> ev.entry_in_daily_file >> ev.vertex_pl >> ev.ecc_id >> ev.vertex_material
	   >> ev.recon_vertex_position[0] >> ev.recon_vertex_position[1] >> ev.recon_vertex_position[2]
	   >> ev.true_vertex_position[0] >> ev.true_vertex_position[1] >> ev.true_vertex_position[2]
	   >> ev.nu_energy >> ev.nu_ax >> ev.nu_ay >> chain_num >> true_chain_num ) {
      for ( int i = 0; i < chain_num; i++ ) {
	t.base.clear();
	t.base_pair.clear();
	t.base.reserve(base_num);
	t.base_pair.reserve(base_pair_num);
	ifs >> t.chainid >> t.stop_flag >> t.particle_flag >> t.direction >> t.charge_sign
	    >> t.ecc_range_mom[0] >> t.ecc_range_mom[1] >> t.ecc_mcs_mom[0] >> t.ecc_mcs_mom[1]
	    >> t.bm_range_mom >> t.bm_curvature_mom
	    >> t.ecc_range_mom_error[0][0] >> t.ecc_range_mom_error[0][1] >> t.ecc_range_mom_error[1][0] >> t.ecc_range_mom_error[1][1]
	    >> t.ecc_mcs_mom_error[0][0] >> t.ecc_mcs_mom_error[0][1] >> t.ecc_mcs_mom_error[1][0] >> t.ecc_mcs_mom_error[1][1]
	    >> t.bm_range_mom_error[0] >> t.bm_range_mom_error[1]
	    >> t.bm_curvature_mom_error[0] >> t.bm_curvature_mom_error[1]
	    >> t.muon_likelihood >> t.proton_likelihood
	    >> base_num >> base_pair_num;
	for ( int i = 0; i < base_num; i++ ) {
	  ifs >> m.pl >> m.rawid >> m.ax >> m.ay >> m.x >> m.y >> m.z
	      >> m.m[0].zone >> m.m[0].view >> m.m[0].imager >> m.m[0].ph >> m.m[0].pixelnum >> m.m[0].hitnum
	      >> m.m[1].zone >> m.m[1].view >> m.m[1].imager >> m.m[1].ph >> m.m[1].pixelnum >> m.m[1].hitnum;
	  t.base.push_back(m);
	}
	for ( int i = 0; i < base_pair_num; i++ ) {
	  ifs >> p.first.pl >> p.first.rawid >> p.second.pl >> p.second.rawid
	      >> p.first.ax >> p.first.ay >> p.first.x >> p.first.y >> p.first.z
	      >> p.second.ax >> p.second.ay >> p.second.x >> p.second.y >> p.second.z;
	  t.base_pair.push_back(p);
	}
	InputBasetrackInformation(t.base, t.base_pair);
	ev.chains.push_back(t);
      }
      for ( int i = 0; i < true_chain_num; i++ ) {
	t.base.clear();
	t.base_pair.clear();
	t.base.reserve(base_num);
	t.base_pair.reserve(base_pair_num);
	ifs >> t.chainid >> t.stop_flag >> t.particle_flag >> t.direction >> t.charge_sign
	    >> t.ecc_range_mom[0] >> t.ecc_range_mom[1] >> t.ecc_mcs_mom[0] >> t.ecc_mcs_mom[1]
	    >> t.bm_range_mom >> t.bm_curvature_mom
	    >> t.ecc_range_mom_error[0][0] >> t.ecc_range_mom_error[0][1] >> t.ecc_range_mom_error[1][0] >> t.ecc_range_mom_error[1][1]
	    >> t.ecc_mcs_mom_error[0][0] >> t.ecc_mcs_mom_error[0][1] >> t.ecc_mcs_mom_error[1][0] >> t.ecc_mcs_mom_error[1][1]
	    >> t.bm_range_mom_error[0] >> t.bm_range_mom_error[1]
	    >> t.bm_curvature_mom_error[0] >> t.bm_curvature_mom_error[1]
	    >> t.muon_likelihood >> t.proton_likelihood
	    >> base_num >> base_pair_num;
	for ( int i = 0; i < base_num; i++ ) {
	  ifs >> m.pl >> m.rawid >> m.ax >> m.ay >> m.x >> m.y >> m.z
	      >> m.m[0].zone >> m.m[0].view >> m.m[0].imager >> m.m[0].ph >> m.m[0].pixelnum >> m.m[0].hitnum
	      >> m.m[1].zone >> m.m[1].view >> m.m[1].imager >> m.m[1].ph >> m.m[1].pixelnum >> m.m[1].hitnum;
	  t.base.push_back(m);
	}
	for ( int i = 0; i < base_pair_num; i++ ) {
	  ifs >> p.first.pl >> p.first.rawid >> p.second.pl >> p.second.rawid
	      >> p.first.ax >> p.first.ay >> p.first.x >> p.first.y >> p.first.z
	      >> p.second.ax >> p.second.ay >> p.second.x >> p.second.y >> p.second.z;
	  t.base_pair.push_back(p);
	}
	InputBasetrackInformation(t.base, t.base_pair);
	ev.true_chains.push_back(t);
      }
      ret.push_back(ev);

      ev.chains.clear();
      ev.true_chains.clear();
    }

    return ret;

  }

  void InputBasetrackInformation(std::vector<Mom_basetrack > &base, std::vector<std::pair<Mom_basetrack, Mom_basetrack > > &base_pair) {
    std::multimap<std::pair<int, int>, Mom_basetrack*> base_pair_p;
    for ( auto itr = base_pair.begin(); itr != base_pair.end(); itr++ ) {
      base_pair_p.insert(std::make_pair(std::make_pair(itr->first.pl, itr->first.rawid), &(itr->first)));
      base_pair_p.insert(std::make_pair(std::make_pair(itr->second.pl, itr->second.rawid), &(itr->second)));
    }
    for ( auto itr = base.begin(); itr != base.end(); itr++ ) {
      std::pair<int, int> id = std::make_pair(itr->pl, itr->rawid);
      if ( base_pair_p.count(id) == 0 ) continue;
      auto range = base_pair_p.equal_range(id);
      for ( auto res = range.first; res != range.second; res++ ) {
	res->second->m[0] = itr->m[0];
	res->second->m[1] = itr->m[1];
      }
    }

  }

  microtrack_minimum::microtrack_minimum() {
    ph = -1;
    zone = -1;
    imager = -1;
    pixelnum = -1;
    hitnum = -1;
  }

  Mom_basetrack::Mom_basetrack() {
    pl = -1;
    rawid = -1;
    ax = NAN;
    ay = NAN;
    x = NAN;
    y = NAN;
    z = NAN;
  }

  Mom_chain::Mom_chain() {
    chainid = -1;
    stop_flag = -1;
    particle_flag = -1;
    direction = 0;
    charge_sign = 0;
    ecc_range_mom[0] = -1;
    ecc_range_mom[1] = -1;
    ecc_mcs_mom[0] = -1;
    ecc_mcs_mom[1] = -1;
    bm_range_mom = -1;
    charge_sign = -1;
    ecc_range_mom_error[0][0] = -1;
    ecc_range_mom_error[0][1] = -1;
    ecc_range_mom_error[1][0] = -1;
    ecc_range_mom_error[1][1] = -1;
    ecc_mcs_mom_error[0][0] = -1;
    ecc_mcs_mom_error[0][1] = -1;
    ecc_mcs_mom_error[1][0] = -1;
    ecc_mcs_mom_error[1][1] = -1;
    bm_range_mom_error[0] = -1;
    bm_range_mom_error[1] = -1;
    bm_curvature_mom_error[0] = -1;
    bm_curvature_mom_error[1] = -1;
    muon_likelihood = NAN;
    proton_likelihood = NAN;
  }

  Event_information::Event_information() {
    groupid = -1;
    unixtime = -1;
    tracker_track_id = -1;
    entry_in_daily_file = -1;
    vertex_pl = -1;
    ecc_id = -1;
    vertex_material = -1;
    recon_vertex_position[0] = NAN;
    recon_vertex_position[1] = NAN;
    recon_vertex_position[2] = NAN;
    true_vertex_position[0] = NAN;
    true_vertex_position[1] = NAN;
    true_vertex_position[2] = NAN;
    weight = -1;
    nu_energy = -1;
    nu_ax = NAN;
    nu_ay = NAN;
  }

}
