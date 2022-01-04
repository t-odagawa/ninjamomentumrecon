#ifndef MCS_CLASS_HPP
#define MCS_CLASS_HPP

#include <fstream>
#include <string>
#include <vector>

namespace Momentum_recon
{

  class microtrack_minimum {
  public :
    int ph, zone, view, imager, pixelnum, hitnum;
  };

  class Mom_basetrack {
  public :
    int pl, rawid;
    float ax, ay, x, y, z;
    microtrack_minimum m[2];
  };

  class Mom_chain {
  public :
    int groupid, chainid, unixtime, tracker_track_id, entry_in_daily_file;
    int stop_flag, particle_flag;
    double ecc_range_mom, ecc_mcs_mom, bm_range_mom, bm_curvature_mom;
    // basetrack coordinate
    std::vector<Mom_basetrack > base;
    // local coordinate
    std::vector<std::pair<Mom_basetrack, Mom_basetrack > > base_pair;
  };

  std::vector<Mom_chain > ReadMomChain(std::string filename);
  void WriteMomChain(std::string filename, std::vector<Mom_chain > &mom_chain_vector);
  bool ReadMomChainHeader(std::ifstream &ifs, Mom_chain &mom_chain, int &num_base, int &num_link);
  void WriteMomChainHeader(std::ofstream &ofs, Mom_chain &mom_chain);
  void WriteMomChainText(std::string filename, std::vector<Mom_chain > &mom_chain_vector);
}

class SimpleMCBaseTrackCollection{
public:
  int particle_id, plate_id, ecc_id;
  double energy_deposit_1, energy_deposit_2, ax, ay, x, y, z;
  int track_id, event_id, file_id;
  double weight;
};


#endif
