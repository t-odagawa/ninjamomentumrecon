#ifndef NINJA_ANGLE_METHOD_HPP
#define NINJA_ANGLE_METHOD_HPP

#include <vector>

#include <B2SpillSummary.hh>
#include "McsCommon.hpp"

bool mcs_angle_method(const B2SpillSummary &spill_summary, 
		      std::vector<double> &recon_pb,
		      std::vector<double> &angle_difference,
		      std::vector<double> &path_length,
		      int particle_id);

double calculate_pb_average(double dz, std::array<std::pair<double, int>, MAX_NUM_SKIP> theta_rms);

double calculate_pb(double theta_rms, double rad_length);

double calculate_radiation_length(int skip, double dz);


/**
 * Get angle difference in lateral coordinate of upstream base track
 * @param tangent_up tangent vector of upstream film
 * @param tangent_down tangent vector of downstream film
 * @return angle difference in lateral coordinate of upstream base track
 */
double get_tangent_difference_lateral(TVector3 tangent_up, TVector3 tangent_down, TVector3 vertex_tangent);

double get_angle_difference_lateral(TVector3 tangent_up, TVector3 tangent_down, TVector3 vertex_tangent);

double get_tangent_difference_radial(TVector3 tangent_up, TVector3 tangent_down);
#endif
