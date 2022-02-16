#ifndef MCS_CLASS_HPP
#define MCS_CLASS_HPP

#include <fstream>
#include <string>
#include <vector>

namespace Momentum_recon
{

  class microtrack_minimum {
  public :
    int zone, view, imager, ph, pixelnum, hitnum;
    microtrack_minimum();
  };

  class Mom_basetrack {
  public :
    int pl, rawid;
    float ax, ay, x, y, z;
    microtrack_minimum m[2];
    Mom_basetrack();
  };

  class Mom_chain {
  public :
    int chainid;
    // stop_flag =  0 : penetrate, 1 : bm stop, 2 : ecc stop
    int stop_flag, particle_flag, direction, charge_sign;
    // [0] : muon (pion), [1] : proton
    double ecc_range_mom[2], ecc_mcs_mom[2], bm_range_mom, bm_curvature_mom;
    // best fit からの距離 (絶対値)
    double ecc_range_mom_error[2][2], ecc_mcs_mom_error[2][2], bm_range_mom_error[2], bm_curvature_mom_error[2];
    // likelihood = (vph - mean) / sigma (符号付き)
    double muon_likelihood, proton_likelihood;
    // basetrack coordinate
    std::vector<Mom_basetrack > base;
    // local coordinate
    std::vector<std::pair<Mom_basetrack, Mom_basetrack > > base_pair;
    Mom_chain();
  };

  class Event_information {
  public :
    int groupid, unixtime, tracker_track_id, entry_in_daily_file, vertex_pl, ecc_id;
    // 0 : water, 2 : iron
    int vertex_material;
    double recon_vertex_position[3], true_vertex_position[3];
    double weight, nu_energy, nu_ax, nu_ay; // energy : MeV
    std::vector<Mom_chain > chains;
    std::vector<Mom_chain > true_chains;
    Event_information();
  };

  void WriteEventInformationHeader(std::ofstream &ofs, Event_information &ev);
  bool ReadEventInformationHeader(std::ifstream &ifs, Event_information &ev, int &num_chain, int &num_true_chain);
  void WriteMomChainHeader(std::ofstream &ofs, Mom_chain &mom_chain);
  bool ReadMomChainHeader(std::ifstream &ifs, Mom_chain &mom_chain, int &num_base, int &num_link);
  void WriteEventInformationBin(std::string filename, std::vector<Event_information > &ev_vector);
  std::vector<Event_information > ReadEventInformationBin(std::string filename);
  void WriteEventInformationTxt(std::string filename, std::vector<Event_information > &ev_vector);
  std::vector<Event_information > ReadEventInformationTxt(std::string filename);
  void InputBasetrackInformation(std::vector<Mom_basetrack > &base, std::vector<std::pair<Mom_basetrack, Mom_basetrack > > &base_pair);
}

class SimpleMCBaseTrackCollection{
public:
  int particle_id, plate_id, ecc_id;
  double energy_deposit_1, energy_deposit_2, ax, ay, x, y, z;
  int track_id, event_id, file_id;
  double weight;
};


#endif
