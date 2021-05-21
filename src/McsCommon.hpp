#ifndef MCS_COMMON_HPP
#define MCS_COMMON_HPP

enum {kNinjaIron = 0,
      kNinjaWater = 1,
      kNinjaGel = 2,
      kNinjaBase = 3,
      kNinjaPacking = 4,
      kNumberOfNinjaMaterials = 5
} NinjaMaterial;

///> Radiation lengths of each material
//const static double IRON_RAD_LENGTH = 17.58; // mm
const static double IRON_RAD_LENGTH = 17.18; // mm
const static double WATER_RAD_LENGTH = 360.8; // mm
const static double GEL_RAD_LENGTH = 30.3; // mm
const static double POLY_RAD_LENGTH = 413.1; // mm
const static double RAD_LENGTH[kNumberOfNinjaMaterials] = { IRON_RAD_LENGTH,
							    WATER_RAD_LENGTH,
							    GEL_RAD_LENGTH,
							    POLY_RAD_LENGTH,
							    POLY_RAD_LENGTH
};

///> Thickness of each material
const static double IRON_THICK = 0.5; // mm
const static double WATER_THICK = 2.3; // mm
const static double GEL_THICK = 0.07; // mm
const static double BASE_THICK = 0.21; // mm
const static double PACK_THICK = 0.109; // mm
const static double MATERIAL_THICK[kNumberOfNinjaMaterials] = { IRON_THICK,
								WATER_THICK,
								GEL_THICK,
								BASE_THICK,
								PACK_THICK
};

///> Error of angle difference
const static double D_ANG_ERROR = 0.0022;
///> Error of position difference
const static double D_POS_ERROR = 0.006; // mm

///> Minimum number of intervals for plate pairs
const static int MIN_NUM_SKIP = 1;
///> Maximum number of intervals for plate pairs
const static int MAX_NUM_SKIP = 2;
///> Number of plates in one NINJA ECC
const static int MAX_NUM_ECC_PLATE = 133;

std::vector<Double_t> get_true_pbeta(const B2SpillSummary &spill_summary, int particle_id);

bool emulsion_compare(const B2EmulsionSummary *lhs, const B2EmulsionSummary *rhs);

#endif
