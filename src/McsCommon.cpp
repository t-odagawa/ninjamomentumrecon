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

TVector3 smear_position_vector(TVector3 &position, Int_t material) {

  Double_t delta_x = gRandom->Gaus(0., xy_position_accuracy[material]);
  Double_t delta_y = gRandom->Gaus(0., xy_position_accuracy[material]);
  Double_t delta_z = gRandom->Gaus(0., z_position_accuracy[material]);

  position.SetX(position.X() + delta_x);
  position.SetY(position.Y() + delta_y);
  position.SetZ(position.Z() + delta_z);

  return position;

}

TVector3 smear_tangent_vector(TVector3 &tangent, Int_t material) {

  Double_t tan_theta = std::hypot(tangent.X(), tangent.Y());
  Double_t tan_phi = tangent.Y() / tangent.X();
  Double_t cos_phi = std::sqrt(1. / (1 + tan_phi * tan_phi));
  Double_t sin_phi = tan_phi * cos_phi;

  Double_t delta_radial = gRandom->Gaus(0., radial_angle_accuracy(tan_theta, material));
  Double_t delta_lateral = gRandom->Gaus(0., lateral_angle_accuracy(material));

  Double_t delta_ay = delta_radial * cos_phi + delta_lateral * sin_phi;
  Double_t delta_ax = delta_radial * sin_phi - delta_lateral * cos_phi;

  tangent.SetX(tangent.X() + delta_ax);
  tangent.SetY(tangent.Y() + delta_ay);

  return tangent;

}

Double_t radial_angle_accuracy(Double_t tangent_theta, Int_t material) {
  return std::sqrt(2) / 210.e-3 * std::sqrt(xy_position_accuracy[material] * xy_position_accuracy[material]
					    + tangent_theta * tangent_theta * z_position_accuracy[material] * z_position_accuracy[material]);

}

Double_t lateral_angle_accuracy(Int_t material) {
  return std::sqrt(2) / 210.e-3 * xy_position_accuracy[material];
}
