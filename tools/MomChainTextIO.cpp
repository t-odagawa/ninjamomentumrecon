#include <fstream>
#include <iostream>
#include <string>
#include <iomanip>
#include <vector>

namespace Momentum_recon {

  class microtrack_minimum {
  public:
    int ph, zone, view, imager, pixelnum, hitnum;
  };
  class Mom_basetrack {
  public:
    int pl, rawid;
    float ax, ay, x, y, z;
    microtrack_minimum m[2];
  };
  class Mom_chain {
  public:
    //stop_flag 0:penetarte / 1:bm stop / 2:ecc stop
    int chainid, groupid, unixtime, tracker_track_id, entry_in_daily_file, stop_flag, particle_flag;
    double ecc_range_mom, ecc_mcs_mom, bm_range_mom, bm_curvature_mom;
    std::vector<Mom_basetrack> base;
    std::vector<std::pair<Mom_basetrack, Mom_basetrack>> base_pair;
  };

}

bool ReadMomChainHeader(std::ifstream &ifs, Momentum_recon::Mom_chain &mom_chain, int &num_base, int &num_link) {

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

std::vector<Momentum_recon::Mom_chain> ReadMomChainBin(std::string filename) {

  std::ifstream ifs(filename, std::ios::binary);

  std::vector<Momentum_recon::Mom_chain> mom_chain_vec;
  Momentum_recon::Mom_chain mom_chain;
  Momentum_recon::Mom_basetrack mom_basetrack;
  std::pair<Momentum_recon::Mom_basetrack, Momentum_recon::Mom_basetrack> mom_basetrack_pair;
  int num_base, num_link;

  while ( ReadMomChainHeader(ifs, mom_chain, num_base, num_link) ) {
    mom_chain.base.clear();
    mom_chain.base_pair.clear();
    mom_chain.base.reserve(num_base);
    mom_chain.base_pair.reserve(num_link);
    for ( int ibase = 0; ibase < num_base; ibase++ ) {
      ifs.read((char*)& mom_basetrack, sizeof(Momentum_recon::Mom_basetrack));
      mom_chain.base.push_back(mom_basetrack);
    }
    for ( int ilink = 0; ilink < num_link; ilink++ ) {
      ifs.read((char*)& mom_basetrack_pair.first, sizeof(Momentum_recon::Mom_basetrack));
      ifs.read((char*)& mom_basetrack_pair.second, sizeof(Momentum_recon::Mom_basetrack));
      mom_chain.base_pair.push_back(mom_basetrack_pair);
    }
    mom_chain_vec.emplace_back(mom_chain);
  }

  return mom_chain_vec;

}

void WriteMomChainText(std::string filename, std::vector<Momentum_recon::Mom_chain> &mom_chain_vector) {

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


int main ( int argc, char* argv[]) {

  if ( argc != 3 ) {
    std::cout << "Usage : " << argv[0]
	      << " <input momch binary file path> <output momch text file path>" << std::endl;
    std::exit(1);
  }

  std::string ifilename = argv[1];
  std::string ofilename = argv[2];

  std::vector<Momentum_recon::Mom_chain> mom_chain_vector = ReadMomChainBin(ifilename);

  WriteMomChainText(ofilename, mom_chain_vector);

  std::exit(0);

}
