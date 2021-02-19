#include "McsCoordinateMethod.hpp"

// system includes
#include <vector>

// boost include
#include <boost/log/trivial.hpp>

#include <B2SpillSummary.hh>
#include <B2TrackSummary.hh>
#include <B2EmulsionSummary.hh>

void mcs_coordinate_method(const B2SpillSummary &spill_summary, int particle_id) {

  // Get emulsion tracks
  std::vector<const B2EmulsionSummary*> emulsions;
  auto it_emulsion = spill_summary.BeginEmulsion();
  while (const auto *emulsion = it_emulsion.Next()) {
    if (emulsion->GetParentTrackId() == 0) continue;
    emulsions.push_back(emulsion);
    int particle_id = emulsion->GetParentTrack().GetParticlePdg();
    BOOST_LOG_TRIVIAL(debug) << "Particle ID : " << particle_id;
  }

  if (emulsions.size() <= 0) return;
  BOOST_LOG_TRIVIAL(debug) << "# of emulsion tracks : " << emulsions.size();
}
