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
const static double IRON_RAD_LENGTH = 17.18;
const static double WATER_RAD_LENGTH = 360.8;
const static double GEL_RAD_LENGTH = 30.3;
const static double POLY_RAD_LENGTH = 413.1;
const static double RAD_LENGTH[kNumberOfNinjaMaterials] = { IRON_RAD_LENGTH,
							    WATER_RAD_LENGTH,
							    GEL_RAD_LENGTH,
							    POLY_RAD_LENGTH,
							    POLY_RAD_LENGTH
};

///> Thickness of each material (mm)
const static double IRON_THICK = 0.5;
const static double WATER_THICK = 2.3;
const static double GEL_THICK = 0.07;
const static double BASE_THICK = 0.21;
const static double PACK_THICK = 0.109;
const static double MATERIAL_THICK[kNumberOfNinjaMaterials] = { IRON_THICK,
								WATER_THICK,
								GEL_THICK,
								BASE_THICK,
								PACK_THICK
};

///> Mass of particles in interest (MeV/c2)
const static double MCS_MUON_MASS = 105.658;
const static double MCS_PION_MASS = 139.570;
const static double MCS_PROTON_MASS = 938.262;
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
const static double XY_POSITION_ACCURACY = .3e-3;
const static double Z_POSITION_ACCURACY = 4.e-3;

///> Scale factor for MCS correction
const static double MCS_SCALE_FACTOR = 1.0265;

///> Energy deposit function parameters for iron
const std::array<double, 5> IRON_ENEDEP_FUNC_PAR = {  11.3616,
						       -41.0421,
						        63.9158,
						       -47.3527,
						        13.7393
};

///> Energy deposit function parameters for water
const std::array<double, 5> WATER_ENEDEP_FUNC_PAR = {  16.677,
						        -68.1774,
						        114.435,
						        -88.6891,
						         26.3072
};

#endif
