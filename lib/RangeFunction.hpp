#ifndef RANGE_FUNCTION_HPP
#define RANGE_FUNCTION_HPP

#include <vector>

double IronRangeFromMomentum(double momentum);

double IronMomentumFromRange(double range);

double WaterRangeFromMomentum(double momentum);

double WaterMomentumFromRange(double range);

double PolystyreneRangeFromMomentum(double momentum);

double PolystyreneMomentumFromRange(double range);

double GelRangeFromMomentum(double momentum);

double GelMomentumFromRange(double range);

void ModifyVectors(std::vector<double> &ax, std::vector<double> &ay, std::vector<int> &pl);

double CalculateMomentumFromRange(std::vector<double> ax, std::vector<double> ay, std::vector<int> pl);

#endif
