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

std::vector<double> mcs_angle_method(const B2SpillSummary &spill_summary, 
				     std::vector<double> &angle_difference, 
				     std::vector<double> &path_length,
				     int particle_id) {

  std::vector<double> recon_pb = {};

  // Get emulsion tracks
  std::vector<const B2EmulsionSummary*> emulsions;
  auto it_emulsion = spill_summary.BeginEmulsion();
  while (const auto *emulsion = it_emulsion.Next()) {
    if (emulsion->GetParentTrackId() == 0) continue;
    if (emulsion->GetParentTrack().GetParticlePdg() != particle_id) continue;
    if (emulsion->GetFilmType() != B2EmulsionType::kECC) continue;
    emulsions.push_back(emulsion);
  }

  if (emulsions.size() <= 0) return recon_pb;

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
	  if (emulsion_down->GetPlate() < 16) continue;
	  int plate_difference = emulsion_up->GetPlate() - emulsion_down->GetPlate();
	  if (plate_difference % 2 == 1 && plate_difference / 2 < MAX_NUM_SKIP - 1) {
	    BOOST_LOG_TRIVIAL(debug) << "up plate : " << emulsion_up->GetPlate()
				     << " down plate : " << emulsion_down->GetPlate();

	    TVector3 tangent_down = emulsion_down->GetTangent().GetValue();
	    theta_rms.at((plate_difference + 1) / 2).first += get_angle_difference_lateral(tangent_up, tangent_down)
	      * get_angle_difference_lateral(tangent_up, tangent_down);
	    theta_rms.at((plate_difference + 1 ) / 2).second++;
	    angle_difference.push_back(get_angle_difference_lateral(tangent_up, tangent_down));

	  }
	}
      }

      for (int iskip = 1; iskip < MAX_NUM_SKIP; iskip++) {
	theta_rms.at(iskip).first /= (double) theta_rms.at(iskip).second;
	theta_rms.at(iskip).first = std::sqrt(theta_rms.at(iskip).first);
      }

      // Reconstruct pbeta from RMS of angle difference distribution
      recon_pb.push_back(reconstruct_pbeta(dz, theta_rms));
      emulsions_one_track.clear();
    }

    emulsions_one_track.push_back(emulsion);
    track_id_tmp_ = track_id_;
    ecc_tmp_ = ecc_;
  }

  BOOST_LOG_TRIVIAL(debug) << "# of particles which make tracks in emulsions : " << number_of_tracks;
  BOOST_LOG_TRIVIAL(debug) << "Reconstructed pbeta = " << recon_pb.front();
  return recon_pb;

}

double get_pb_oneskip(int skip, double dz, double theta_rms) {

  if (skip >= MAX_NUM_SKIP)
    throw std::out_of_range("skip should be less than MAX_NUM_SKIP");

  double pb2 = 0;

  for (int material = 0; material < kNumberOfNinjaMaterials; material++) {
    switch (material) {
    case kNinjaIron :
      pb2 += calculate_one_material(skip, dz, material) * calculate_one_material(skip, dz, material);
      break;
    case kNinjaWater :
      pb2 += calculate_one_material(skip - 1, dz, material) * calculate_one_material(skip - 1, dz, material);
      break;
    case kNinjaGel :
      pb2 += calculate_one_material(4 * skip - 2, dz, material) * calculate_one_material(4 * skip - 2, dz, material);
      break;
    case kNinjaBase :
      pb2 += calculate_one_material(2 * skip - 2, dz, material) * calculate_one_material(2 * skip - 2, dz, material);
      break;
    case kNinjaPacking :
      pb2 += calculate_one_material(2 * skip - 2, dz, material) * calculate_one_material(2 * skip - 2, dz, material);
      break;
    }
  }

  return std::sqrt(pb2) / theta_rms;
}

double calculate_one_material(int num_layer, double dz, int material) {
  if (num_layer == 0) return 0.;
  return 13.6 * std::sqrt(num_layer * dz * MATERIAL_THICK[material] / RAD_LENGTH[material])
    * (1 + 0.038 * std::log(num_layer * dz *  MATERIAL_THICK[material] / RAD_LENGTH[material]));
}

double reconstruct_pbeta(double dz, std::array<std::pair<double, int>, MAX_NUM_SKIP> theta_rms) {

  std::vector<double> pbeta = {};

  for (int iskip = MIN_NUM_SKIP; iskip < MAX_NUM_SKIP; iskip++) {
    if (theta_rms.at(iskip).second < 5) continue;
    BOOST_LOG_TRIVIAL(debug) << "theta_rms.at(" << iskip << ").first = " << theta_rms.at(iskip).first
    			     << " theta_rms.at(" << iskip << ").second = " << theta_rms.at(iskip).second;
    pbeta.push_back(get_pb_oneskip(iskip, dz, theta_rms.at(iskip).first));
  }

  return std::accumulate(pbeta.begin(), pbeta.end(), 0.)
    / pbeta.size();
}

double get_angle_difference_lateral(TVector3 tangent_up, TVector3 tangent_down) {
  return (1. / std::sqrt(tangent_up.X() * tangent_up.X() + tangent_up.Y() * tangent_up.Y()))
    * ((tangent_down.X() - tangent_up.X()) * tangent_up.Y()
       - (tangent_down.Y() - tangent_up.Y()) * tangent_up.X());
}

double get_angle_difference_radial(TVector3 tangent_up, TVector3 tangent_down) {
  return (1. / std::sqrt(tangent_up.X() * tangent_up.X() + tangent_up.Y() * tangent_up.Y()))
    * ((tangent_down.X() - tangent_up.X()) * tangent_up.X()
       + (tangent_down.Y() - tangent_up.Y()) * tangent_up.Y());
}
