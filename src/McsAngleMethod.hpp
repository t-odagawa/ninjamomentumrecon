#ifndef NINJA_ANGLE_METHOD_HPP
#define NINJA_ANGLE_METHOD_HPP

#include <vector>

#include <B2SpillSummary.hh>
#include "McsCommon.hpp"

std::vector<double> mcs_angle_method(const B2SpillSummary &spill_summary, int particle_id);

double get_pb_oneskip(int skip, double dz, double theta_rms);

double calculate_one_material(int num_layer, double dz, int material);

double reconstruct_pbeta(double dz, std::array<std::pair<double, int>, MAX_NUM_SKIP> theta_rms);

/**
 * Get angle difference in lateral coordinate of upstream base track
 * @param tangent_up tangent vector of upstream film
 * @param tangent_down tangent vector of downstream film
 * @return angle difference in lateral coordinate of upstream base track
 */
double get_angle_difference(TVector3 tangent_up, TVector3 tangent_down);
#endif
