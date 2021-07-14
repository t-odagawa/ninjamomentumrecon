#include <vector>

#include <B2SpillSummary.hh>
#include <B2EmulsionSummary.hh>

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

TVector3 smear_position_vector(TVector3 &position) {

}
