#include <fstream>
#include <iostream>
#include <iomanip>
#include <string>
#include <vector>

#include "McsClass.hpp"

namespace Momentum_recon
{

  std::vector<Mom_chain > ReadMomChain(std::string filename) {

    std::vector<Mom_chain > ret;
    std::ifstream ifs(filename, std::ios::binary);
    Mom_chain mom_chain;
    Mom_basetrack basetrack;
    std::pair<Mom_basetrack, Mom_basetrack> link;
    int num_base, num_link;

    while ( ReadMomChainHeader(ifs, mom_chain, num_base, num_link) ) {
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

      ret.emplace_back(mom_chain);
    }

    return ret;

  }

  void WriteMomChain(std::string filename, std::vector<Mom_chain > &mom_chain_vector) {

    std::ofstream ofs(filename, std::ios::binary);
    if ( !ofs ) {
      std::cerr << "File "<< filename << " cannot be open !" << std::endl;
      std::exit(1);
    }
    else if ( mom_chain_vector.empty() ) {
      std::cerr << "Target mom_chain data...null" << std::endl;
      std::exit(1);
    }
    else {
      int num_base = 0;
      int num_link = 0;
      for ( auto mom_chain : mom_chain_vector ) {
	num_base = mom_chain.base.size();
	num_link = mom_chain.base_pair.size();
	WriteMomChainHeader(ofs, mom_chain);
	for ( auto base : mom_chain.base )
	  ofs.write((char*)& base, sizeof(Mom_basetrack));
	for ( auto link : mom_chain.base_pair )	{
	  ofs.write((char*)& link.first, sizeof(Mom_basetrack));
	  ofs.write((char*)& link.second, sizeof(Mom_basetrack));
	}
      }
    }
  }

  bool ReadMomChainHeader(std::ifstream &ifs, Mom_chain &mom_chain, int &num_base, int &num_link) {
    if (!ifs.read((char*)& mom_chain.groupid, sizeof(int))) return false;
    if (!ifs.read((char*)& mom_chain.chainid, sizeof(int))) return false;
    if (!ifs.read((char*)& mom_chain.unixtime, sizeof(int))) return false;
    if (!ifs.read((char*)& mom_chain.tracker_track_id, sizeof(int))) return false;
    if (!ifs.read((char*)& mom_chain.entry_in_daily_file, sizeof(int))) return false;
    if (!ifs.read((char*)& mom_chain.stop_flag, sizeof(int))) return false;
    if (!ifs.read((char*)& mom_chain.particle_flag, sizeof(int))) return false;
    if (!ifs.read((char*)& mom_chain.ecc_range_mom, sizeof(double))) return false;
    if (!ifs.read((char*)& mom_chain.ecc_mcs_mom, sizeof(double))) return false;
    if (!ifs.read((char*)& mom_chain.bm_range_mom, sizeof(double))) return false;
    if (!ifs.read((char*)& mom_chain.bm_curvature_mom, sizeof(double))) return false;
    if (!ifs.read((char*)& num_base, sizeof(int))) return false;
    if (!ifs.read((char*)& num_link, sizeof(int))) return false;
    return true;
  }

  void WriteMomChainHeader(std::ofstream &ofs, Mom_chain &mom_chain) {
    int num_base = mom_chain.base.size();
    int num_link = mom_chain.base_pair.size();
    
    ofs.write((char*)& mom_chain.groupid, sizeof(int));
    ofs.write((char*)& mom_chain.chainid, sizeof(int));
    ofs.write((char*)& mom_chain.unixtime, sizeof(int));
    ofs.write((char*)& mom_chain.tracker_track_id, sizeof(int));
    ofs.write((char*)& mom_chain.entry_in_daily_file, sizeof(int));
    ofs.write((char*)& mom_chain.stop_flag, sizeof(int));
    ofs.write((char*)& mom_chain.particle_flag, sizeof(int));
    ofs.write((char*)& mom_chain.ecc_range_mom, sizeof(double));
    ofs.write((char*)& mom_chain.ecc_mcs_mom, sizeof(double));
    ofs.write((char*)& mom_chain.bm_range_mom, sizeof(double));
    ofs.write((char*)& mom_chain.bm_curvature_mom, sizeof(double));
    ofs.write((char*)& num_base, sizeof(int));
    ofs.write((char*)& num_link, sizeof(int));
  }

  void WriteMomChainText(std::string filename, std::vector<Mom_chain > &mom_chain_vector) {
    std::ofstream ofs(filename);
    for ( auto &mom_chain : mom_chain_vector ) {
      ofs << std::right << std::fixed
	  << std::setw(3) << std::setprecision(0) << mom_chain.groupid << " "
	  << std::setw(3) << std::setprecision(0) << mom_chain.chainid << " "
	  << std::setw(10) << std::setprecision(0) << mom_chain.unixtime << " "
	  << std::setw(3) << std::setprecision(0) << mom_chain.tracker_track_id << " "
	  << std::setw(5) << std::setprecision(0) << mom_chain.entry_in_daily_file << " "
	  << std::setw(3) << std::setprecision(0) << mom_chain.stop_flag << " "
	  << std::setw(3) << std::setprecision(0) << mom_chain.particle_flag << " "
	  << std::setw(7) << std::setprecision(1) << mom_chain.ecc_range_mom << " "
	  << std::setw(7) << std::setprecision(1) << mom_chain.ecc_mcs_mom << " "
	  << std::setw(7) << std::setprecision(1) << mom_chain.bm_range_mom << " "
	  << std::setw(7) << std::setprecision(1) << mom_chain.bm_curvature_mom << " "
	  << std::setw(4) << std::setprecision(0) << mom_chain.base.size() << " "
	  << std::setw(4) << std::setprecision(0) << mom_chain.base_pair.size() << std::endl;
      
      for ( auto &b : mom_chain.base) {
	ofs << std::right << std::fixed
	    << std::setw(4) << std::setprecision(0) << b.pl << " "
	    << std::setw(10) << std::setprecision(0) << b.rawid << " "
	    << std::setw(7) << std::setprecision(4) << b.ax << " "
	    << std::setw(7) << std::setprecision(4) << b.ay << " "
	    << std::setw(8) << std::setprecision(1) << b.x << " "
	    << std::setw(8) << std::setprecision(1) << b.y << " "
	    << std::setw(8) << std::setprecision(1) << b.z << " "
	    << std::setw(3) << std::setprecision(0) << b.m[0].zone << " "
	    << std::setw(4) << std::setprecision(0) << b.m[0].view << " "
	    << std::setw(3) << std::setprecision(0) << b.m[0].imager << " "
	    << std::setw(7) << std::setprecision(0) << b.m[0].ph << " "
	    << std::setw(5) << std::setprecision(0) << b.m[0].pixelnum << " "
	    << std::setw(5) << std::setprecision(0) << b.m[0].hitnum << " "
	    << std::setw(3) << std::setprecision(0) << b.m[1].zone << " "
	    << std::setw(4) << std::setprecision(0) << b.m[1].view << " "
	    << std::setw(3) << std::setprecision(0) << b.m[1].imager << " "
	    << std::setw(7) << std::setprecision(0) << b.m[1].ph << " "
	    << std::setw(5) << std::setprecision(0) << b.m[1].pixelnum << " "
	    << std::setw(5) << std::setprecision(0) << b.m[1].hitnum << std::endl;
      }
      
      for (auto &p : mom_chain.base_pair) {
	ofs << std::right << std::fixed
	    << std::setw(4) << std::setprecision(0) << p.first.pl << " "
	    << std::setw(10) << std::setprecision(0) << p.first.rawid << " "
	    << std::setw(4) << std::setprecision(0) << p.second.pl << " "
	    << std::setw(10) << std::setprecision(0) << p.second.rawid << " "
	    << std::setw(7) << std::setprecision(4) << p.first.ax << " "
	    << std::setw(7) << std::setprecision(4) << p.first.ay << " "
	    << std::setw(8) << std::setprecision(1) << p.first.x << " "
	    << std::setw(8) << std::setprecision(1) << p.first.y << " "
	    << std::setw(8) << std::setprecision(1) << p.first.z << " "
	    << std::setw(7) << std::setprecision(4) << p.second.ax << " "
	    << std::setw(7) << std::setprecision(4) << p.second.ay << " "
	    << std::setw(8) << std::setprecision(1) << p.second.x << " "
	    << std::setw(8) << std::setprecision(1) << p.second.y << " "
	    << std::setw(8) << std::setprecision(1) << p.second.z << std::endl;
      }
    }

  }

}

