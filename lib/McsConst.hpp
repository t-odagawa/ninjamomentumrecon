#ifndef MCS_CONST_HPP
#define MCS_CONST_HPP

#include <array>

///> Material enum in NINJA MCS
enum { kNinjaIron = 0,
       kNinjaWater = 1,
       kNinjaGel = 2,
       kNinjaBase = 3,
       kNinjaPacking = 4,
       kNumberOfNinjaMaterials = 5
} NinjaMaterial;

///> Particle enum in NINJA MCS
enum { kNinjaMcsMuon = 0,
       kNinjaMcsPion = 1,
       kNinjaMcsProton = 2,
       kNumberOfNinjaMcsParticles = 3
} NinjaMcsParticle;

///> Direction enum in NINJA MCS
enum { kNinjaMcsForward = 1,
       kNinjaMcsBackward = -1,
       kNumberOfNinjaMcsDirections = 2
} NinjaMcsDirection;

///> Radiation lengths of each material (mm)
const static double NJ_IRON_RAD_LENGTH = 17.18;
const static double NJ_WATER_RAD_LENGTH = 360.8;
const static double NJ_GEL_RAD_LENGTH = 30.3;
const static double NJ_POLY_RAD_LENGTH = 413.1;
const static double RAD_LENGTH[kNumberOfNinjaMaterials] = { NJ_IRON_RAD_LENGTH,
							    NJ_WATER_RAD_LENGTH,
							    NJ_GEL_RAD_LENGTH,
							    NJ_POLY_RAD_LENGTH,
							    NJ_POLY_RAD_LENGTH
};

///> Thickness of each material (mm)
const static double NJ_IRON_THICK = 0.5;
const static double NJ_WATER_THICK = 2.3;
const static double NJ_GEL_THICK = 0.07;
const static double NJ_BASE_THICK = 0.21;
const static double NJ_PACK_THICK = 0.109;
const static double MATERIAL_THICK[kNumberOfNinjaMaterials] = { NJ_IRON_THICK,
								NJ_WATER_THICK,
								NJ_GEL_THICK,
								NJ_BASE_THICK,
								NJ_PACK_THICK
};

///> Density of each material (g/cm3)
const static double NJ_IRON_DENSITY = 8.03;
const static double NJ_WATER_DENSITY = 1.00;
const static double NJ_GEL_DENSITY = 3.816;
const static double NJ_POLY_DENSITY = 1.032;
const static double MATERIAL_DENSITY[kNumberOfNinjaMaterials] = { NJ_IRON_DENSITY,
								  NJ_WATER_DENSITY,
								  NJ_GEL_DENSITY,
								  NJ_POLY_DENSITY,
								  NJ_POLY_DENSITY
};

///> Mass of particles in interest (MeV/c2)
const static double MCS_MUON_MASS = 105.658;
const static double MCS_PION_MASS = 139.570;
const static double MCS_PROTON_MASS = 938.272;
const static double PARTICLE_MASS[kNumberOfNinjaMcsParticles] = { MCS_MUON_MASS,
								  MCS_PION_MASS,
								  MCS_PROTON_MASS
};

///> Directions
const static int MCS_DIRECTION[kNumberOfNinjaMcsDirections] = { kNinjaMcsForward,
								kNinjaMcsBackward
};

///>  Minimum number of nseg
const static int MIN_NUM_SKIP = 1;
///> Maximum number of nseg
const static int MAX_NUM_SKIP = 10;
///> Number of plates in one NINJA ECC
const static int MAX_NUM_ECC_PLATE = 133;

///> Positional alignment accuracy in MC smearing (mm)
///> 0 : iron 1 : water
const static double XY_ALIGN_ACCURACY[2] = {.3e-3, .3e-3};
const static double Z_ALIGN_ACCURACY[2] = {6.e-3, 10.e-3};

///> Position determination accuracy in MC smearing (mm)
// const static double XY_POSITION_ACCURACY = .3e-3;
const static double XY_POSITION_ACCURACY = .245e-3;
//const static double Z_POSITION_ACCURACY = 4.e-3;
const static double Z_POSITION_ACCURACY = 1.64e-3;

///> Scale factor for MCS correction
const static double MCS_SCALE_FACTOR = 1.0265;

///> Energy deposit function parameters for iron
const std::array<double, 5> IRON_ENEDEP_FUNC_PAR = {  11.3616,
						       -41.0421,
						        63.9158,
						       -47.3527,
						        13.7393
};

const std::array<double, 2> IRON_BB_FUNC_PAR = {6.619e-2,
						7.994
};

///> Energy deposit function parameters for water
const std::array<double, 5> WATER_ENEDEP_FUNC_PAR = {  16.677,
						        -68.1774,
						        114.435,
						        -88.6891,
						         26.3072
};

const std::array<double, 2> WATER_BB_FUNC_PAR = {5.083e-2,
						 9.343
};

#endif
