#include <TRandom3.h>

#include <vector>

#include <B2SpillSummary.hh>
#include <B2EmulsionSummary.hh>

#include "McsCommon.hpp"

std::vector<Double_t> get_true_pbeta(const B2SpillSummary &spill_summary, Int_t particle_id) {

  std::vector<Double_t> true_pbeta = {};

}

bool emulsion_compare(const B2EmulsionSummary *lhs, const B2EmulsionSummary *rhs) {
  if (lhs->GetEcc() != rhs->GetEcc())
    return lhs->GetEcc() < rhs->GetEcc();
  else if (lhs->GetParentTrackId() != rhs->GetParentTrackId())
    return lhs->GetParentTrackId() < rhs->GetParentTrackId();
  else 
    return lhs->GetPlate() > rhs->GetPlate();
}

TVector3 smear_distance_vector(TVector3 &distance, Int_t material) {

  Double_t delta_x = gRandom->Gaus(0., xy_align_accuracy[material]);
  Double_t delta_y = gRandom->Gaus(0., xy_align_accuracy[material]);
  Double_t delta_z = gRandom->Gaus(0., z_align_accuracy[material]);

  distance.SetX(distance.X() + delta_x);
  distance.SetY(distance.Y() + delta_y);
  distance.SetZ(distance.Z() + delta_z);

  return distance;

}

TVector3 smear_tangent_vector(TVector3 &tangent, Int_t material) {

  Double_t delta_ay, delta_ax;

  if (std::fabs(tangent.X()) < lateral_tangent_accuracy(material)) {
    delta_ax = - gRandom->Gaus(0., lateral_tangent_accuracy(material));
    delta_ay =   gRandom->Gaus(0., radial_tangent_accuracy(tangent.Y(), material));
  } else {
    
    Double_t tan_theta = std::hypot(tangent.X(), tangent.Y());
    
    Double_t tan_phi = tangent.Y() / tangent.X();
    Double_t cos_phi = tangent.X() / std::fabs(tangent.X())
      * std::sqrt(1. / (1 + tan_phi * tan_phi));
    Double_t sin_phi = tan_phi * cos_phi;
    
    Double_t delta_radial = gRandom->Gaus(0., radial_tangent_accuracy(tan_theta, material));
    Double_t delta_lateral = gRandom->Gaus(0., lateral_tangent_accuracy(material));
    
    delta_ax = delta_radial * cos_phi - delta_lateral * sin_phi;
    delta_ay = delta_radial * sin_phi + delta_lateral * cos_phi;    
  }
  
  
  tangent.SetX(tangent.X() + delta_ax);
  tangent.SetY(tangent.Y() + delta_ay);
  
  return tangent;

}

Double_t radial_tangent_accuracy(Double_t tangent_theta, Int_t material) {
  return std::sqrt(2) / 210.e-3 * std::sqrt(xy_position_accuracy[material] * xy_position_accuracy[material]
					    + tangent_theta * tangent_theta * z_position_accuracy[material] * z_position_accuracy[material]);

}

Double_t lateral_tangent_accuracy(Int_t material) {
  return std::sqrt(2) / 210.e-3 * xy_position_accuracy[material];
}

Double_t new_radial_tangent_accuracy(Double_t tangent_theta, Int_t material) {
  return 5.e-3;
  // return std::sqrt(2) / 210.e-3 / std::sqrt(1 + tangent_theta * tangent_theta) * xy_position_accuracy[material];
}

Double_t new_lateral_tangent_accuracy(Double_t tangent_theta, Int_t material) {
  return 2.e-3;
  /*
  return std::sqrt(2) / 210.e-3 / std::sqrt(1 + xy_position_accuracy[material] * xy_position_accuracy[material]
					    / z_position_accuracy[material] / z_position_accuracy[material]
					    * tangent_theta * tangent_theta)
    * xy_position_accuracy[material];
    */
}
