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

enum { kNinjaMcsMuon = 0,
       kNinjaMcsPion = 1,
       kNinjaMcsProton = 2,
       kNumberOfNinjaMcsParticles = 3
} NinjaMcsParticle;

///> Mass of particles in interest
const static double MCS_MUON_MASS = 105.658; // MeV/c2
const static double MCS_PION_MASS = 139.570; // MeV/c2
const static double MCS_PROTON_MASS = 938.262; // MeV/c2

const static double PARTICLE_MASS[kNumberOfNinjaMcsParticles] = { MCS_MUON_MASS,
								  MCS_PION_MASS,
								  MCS_PROTON_MASS
};

enum { kNinjaMcsForward = 1,
       kNinjaMcsBackward = -1,
       kNumberOfNinjaMcsDirections = 2
} NinjaMcsDirection;

const static int MCS_DIRECTION[kNumberOfNinjaMcsDirections] = { kNinjaMcsForward,
								kNinjaMcsBackward
};

///> Error of angle difference
const static double D_ANG_ERROR = 3.e-3;
///> Error of x/y position difference
const static double D_POS_XY_ERROR = 6.e-3; // mm
///> Error of z position difference
const static double D_POS_Z_ERROR = 4.e-3; // mm

///> Minimum number of intervals for plate pairs
const static int MIN_NUM_SKIP = 1;
///> Maximum number of intervals for plate pairs
const static int MAX_NUM_SKIP = 10;
///> Number of plates in one NINJA ECC
const static int MAX_NUM_ECC_PLATE = 133;

std::vector<Double_t> get_true_pbeta(const B2SpillSummary &spill_summary, int particle_id);

bool emulsion_compare(const B2EmulsionSummary *lhs, const B2EmulsionSummary *rhs);

TVector3 smear_position_vector(TVector3 &position);

TVector3 smear_tangent_vector(TVector3 &tangent);
/*
tangent (x, y, z) -> tangent (rad, lat, z)
rad, lat をsmear -> rad', lat' (新しい座標系ではない)
tangent (rad', lat') -> tangent(x', y', z)

これをそのまま新しい座標系で delta theta_rad, delta theta_lat に変換すれば
smear された角度差になる

radial/lateral angle accuracy はどうやって求める？
Data driven? 前の座標系の radial/lateral の幅を信じて，計算してみる？

iron  xy position accuracy : ~ .3 um (alignment?)
iron  z  position accuracy : ~ 4 um
water xy position accuracy : 
water z  position accuracy : 
radial angle accuracy : sqrt(2) / 210 um * sqrt (xy pos acc^2 + (tan theta)^2 * z pos acc^2)
lateral angle accuracy : sqrt(2) / 210 um * xy pos acc
*/


#endif
