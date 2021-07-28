#include "McsAngleMethod.hpp"
#include "McsCommon.hpp"

// system includes
#include <vector>
#include <algorithm>
#include <numeric>

// boost include
#include <boost/log/trivial.hpp>

#include <B2SpillSummary.hh>
#include <B2TrackSummary.hh>
#include <B2EmulsionSummary.hh>

bool mcs_angle_method(const B2SpillSummary &spill_summary,
				     std::vector<double> &recon_pb,
				     std::vector<double> &angle_difference, 
				     std::vector<double> &path_length,
				     int particle_id) {

  // Get emulsion tracks
  std::vector<const B2EmulsionSummary*> emulsions;
  auto it_emulsion = spill_summary.BeginEmulsion();
  while (const auto *emulsion = it_emulsion.Next()) {
    if (emulsion->GetParentTrackId() == 0) continue;
    if (emulsion->GetParentTrack().GetParticlePdg() != particle_id) continue;
    if (emulsion->GetFilmType() != B2EmulsionType::kECC) continue;
    emulsions.push_back(emulsion);
  }

  if (emulsions.size() <= 0) return false;

  BOOST_LOG_TRIVIAL(debug) << "# of emulsion tracks : " << emulsions.size();

  std::sort(emulsions.begin(), emulsions.end(), emulsion_compare);

  int number_of_tracks = 0;
  std::vector<const B2EmulsionSummary*> emulsions_one_track;

  int track_id_tmp_ = emulsions.at(0)->GetParentTrackId();
  int ecc_tmp_ = emulsions.at(0)->GetEcc();

  for (const auto &emulsion : emulsions) {

    int track_id_ = emulsion->GetParentTrackId();
    int ecc_ = emulsion->GetEcc();

    if (track_id_ != track_id_tmp_ || ecc_ != ecc_tmp_ || emulsion == emulsions.back()) {
      // All base tracks consisting of the same particle's track
      if (emulsion == emulsions.back()) emulsions_one_track.push_back(emulsion);
      number_of_tracks++;

      int vertex_plate = emulsions_one_track.at(0)->GetPlate();
      TVector3 vertex_tangent = emulsions_one_track.at(0)->GetTangent().GetValue();
      double dz = std::sqrt(vertex_tangent.X() * vertex_tangent.X()
			    + vertex_tangent.Y() * vertex_tangent.Y() + 1.);
      path_length.push_back(dz);
      BOOST_LOG_TRIVIAL(debug) << "Vertex Plate ID : " << vertex_plate
			       << " Vertex Tangent X = " << vertex_tangent.X()
			       << ", Vertex Tangent Y = " << vertex_tangent.Y();
      std::array<std::pair<double, int>, MAX_NUM_SKIP> theta_rms;

      for (int i = 0; i < MAX_NUM_SKIP; i++) theta_rms.at(i) = std::make_pair(0., 0);

      // Calculate RMS of angle differences
      for (int iemulsion_up = 0; iemulsion_up < emulsions_one_track.size(); iemulsion_up++) {
	const auto emulsion_up = emulsions_one_track.at(iemulsion_up);
	if (emulsion_up->GetPlate() % 2 == 1
	    || emulsion_up->GetPlate() < 16) continue;
	TVector3 tangent_up = emulsion_up->GetTangent().GetValue();
	for (int iemulsion_down = iemulsion_up + 1; iemulsion_down < emulsions_one_track.size(); iemulsion_down++) {

	  const auto emulsion_down = emulsions_one_track.at(iemulsion_down);
	  if (emulsion_down->GetPlate() < 15) continue;
	  int plate_difference = emulsion_up->GetPlate() - emulsion_down->GetPlate();
	  if (plate_difference % 2 == 1 && plate_difference / 2 < MAX_NUM_SKIP - 1) {
	    BOOST_LOG_TRIVIAL(debug) << "up plate : " << emulsion_up->GetPlate()
				     << " down plate : " << emulsion_down->GetPlate();

	    TVector3 tangent_down = emulsion_down->GetTangent().GetValue();
	    theta_rms.at((plate_difference + 1) / 2).first += get_angle_difference_lateral(tangent_up, tangent_down, vertex_tangent)
	      * get_angle_difference_lateral(tangent_up, tangent_down, vertex_tangent);
	    theta_rms.at((plate_difference + 1 ) / 2).second++;
	    angle_difference.push_back(get_angle_difference_lateral(tangent_up, tangent_down, vertex_tangent));

	  }
	}
      }

      for (int iskip = 1; iskip < MAX_NUM_SKIP; iskip++) {
	theta_rms.at(iskip).first /= (double) theta_rms.at(iskip).second;
	theta_rms.at(iskip).first = std::sqrt(theta_rms.at(iskip).first);
      }

      // Reconstruct pbeta from RMS of angle difference distribution
      recon_pb.push_back(calculate_pb_average(dz, theta_rms));
      emulsions_one_track.clear();
    }

    emulsions_one_track.push_back(emulsion);
    track_id_tmp_ = track_id_;
    ecc_tmp_ = ecc_;
  }

  BOOST_LOG_TRIVIAL(debug) << "# of particles which make tracks in emulsions : " << number_of_tracks;
  BOOST_LOG_TRIVIAL(debug) << "Reconstructed pbeta = " << recon_pb.front();

  return true;

}

double calculate_pb_average(double dz, std::array<std::pair<double, int>, MAX_NUM_SKIP> theta_rms){ 

  std::array<double, MAX_NUM_SKIP> recon_pb = {};

  for (int skip = MIN_NUM_SKIP; skip < MAX_NUM_SKIP; skip++) {
    if (theta_rms.at(skip).second < 5) continue;
    BOOST_LOG_TRIVIAL(debug) << "theta_rms.at(" << skip << ").first = "
			     << theta_rms.at(skip).first
    			     << " theta_rms.at(" << skip << ").second = "
			     << theta_rms.at(skip).second;
    recon_pb.at(skip) = calculate_pb(theta_rms.at(skip).first, calculate_radiation_length(skip, dz));
  }

  return std::accumulate(recon_pb.begin(), recon_pb.end(), 0.) / recon_pb.size();

}

double calculate_pb(double theta_rms, double rad_length) {
  return 13.6 * std::sqrt(rad_length) / theta_rms * (1 + 0.038 * std::log(rad_length));
  // assume beta = 1
}

double calculate_radiation_length(int skip, double dz) {
  if (skip >= MAX_NUM_SKIP)
    throw std::out_of_range("skip should be less than MAX_NUM_SKIP");

  double rad_length = 0.;

  for (int material = 0; material < kNumberOfNinjaMaterials; material++) {
    int num_layers = 0;
    switch (material) {
    case kNinjaIron : 
      num_layers = skip;
      break;
    case kNinjaWater :
      num_layers = skip - 1;
      break;
    case kNinjaGel :
      num_layers = 4 * skip - 2;
      break;
    case kNinjaBase :
      num_layers = 2 * skip - 1;
      break;
    case kNinjaPacking :
      num_layers = 2 * skip - 2;
      break;
    }
    rad_length += num_layers * dz * MATERIAL_THICK[material] / RAD_LENGTH[material];
  }

  // BOOST_LOG_TRIVIAL(debug) << "Radiation length = " << rad_length << "X0";
  return rad_length;

}

double get_tangent_difference_lateral(TVector3 tangent_up, TVector3 tangent_down, TVector3 vertex_tangent) {

  /*  
  return (1. / std::sqrt(vertex_tangent.X() * vertex_tangent.X() + vertex_tangent.Y() * vertex_tangent.Y()))
    * ((tangent_down.X() - tangent_up.X()) * vertex_tangent.Y()
       - (tangent_down.Y() - tangent_up.Y()) * vertex_tangent.X());
  */
  
  return (1. / std::sqrt(tangent_up.X() * tangent_up.X() + tangent_down.Y() * tangent_down.Y()))
    * ((tangent_down.X() - tangent_up.X()) * tangent_up.Y()
       - (tangent_down.Y() - tangent_up.Y()) * tangent_up.X());
  
}

double get_angle_difference_lateral(TVector3 tangent_up, TVector3 tangent_down, TVector3 vertex_tangent) {
  return std::atan(tangent_down.Y() / tangent_down.X()) - std::atan(tangent_up.Y() / tangent_up.X());
}


double get_angle_difference_radial(TVector3 tangent_up, TVector3 tangent_down) {
  return std::atan((tangent_up.X() * tangent_down.X() + tangent_up.Y() * tangent_down.Y())
		   / std::sqrt(tangent_up.X() * tangent_up.X() + tangent_up.Y() * tangent_up.Y()))
    - std::atan(std::sqrt(tangent_up.X() * tangent_up.X() + tangent_up.Y() * tangent_up.Y()));
}

double get_tangent_difference_radial(TVector3 tangent_up, TVector3 tangent_down) {
  return (1. / std::sqrt(tangent_up.X() * tangent_up.X() + tangent_up.Y() * tangent_up.Y()))
    * ((tangent_down.X() - tangent_up.X()) * tangent_up.X()
       + (tangent_down.Y() - tangent_up.Y()) * tangent_up.Y());
}


double get_angle_difference_radial_new(TVector3 tangent_up, TVector3 tangent_down) {
  if (std::fabs(tangent_up.X()) < 1e-2 &&
      std::fabs(tangent_up.Y()) < 1e-2)
    return tangent_down.Y();

  double a = ( - tangent_up.X() * tangent_down.X() - tangent_up.Y() * tangent_down.Y()
	       + tangent_up.X() * tangent_up.X()   + tangent_up.Y() * tangent_up.Y())
    / std::sqrt( ( tangent_up.X() * tangent_up.X() + tangent_up.Y() * tangent_up.Y() ) 
		 * tangent_up.Mag2() );
  double b = (tangent_up.X() * tangent_down.X() + tangent_up.Y() * tangent_down.Y() + 1.)
    / tangent_up.Mag();
  return std::atan(a / b);
}

double get_angle_difference_lateral_new(TVector3 tangent_up, TVector3 tangent_down) {
  if (std::fabs(tangent_up.X()) < 1e-2 &&
      std::fabs(tangent_up.Y()) < 1e-2)
    return tangent_down.X();

  double a = ( - tangent_up.Y() * tangent_down.X() + tangent_up.X() * tangent_down.Y() )
    / std::sqrt(tangent_up.X() * tangent_up.X() + tangent_up.Y() * tangent_up.Y());
  double b = ( tangent_up.X() * tangent_down.X() + tangent_up.Y() * tangent_down.Y() + 1.)
    / tangent_up.Mag();
  return std::atan(a / b);
}
