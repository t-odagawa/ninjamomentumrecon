#include <fstream>

namespace Momentum_recon
{
  
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
    int chainid, groupid, unixtime, tracker_track_id, entry_in_daily_file, stop_flag, particle_flag;
    double ecc_range_mom, ecc_mcs_mom, bm_range_mom, bm_curvature_mom;
    //basetrack座標系
    std::vector<Mom_basetrack > base;
    //localな座標系
    std::vector<std::pair<Mom_basetrack, Mom_basetrack > > base_pair;
  };

}
