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
const static Double_t IRON_RAD_LENGTH = 17.18;
const static Double_t WATER_RAD_LENGTH = 360.8;
const static Double_t GEL_RAD_LENGTH = 30.3;
const static Double_t POLY_RAD_LENGTH = 413.1;
const static Double_t RAD_LENGTH[kNumberOfNinjaMaterials] = { IRON_RAD_LENGTH,
							    WATER_RAD_LENGTH,
							    GEL_RAD_LENGTH,
							    POLY_RAD_LENGTH,
							    POLY_RAD_LENGTH
};

///> Thickness of each material (mm)
const static Double_t IRON_THICK = 0.5;
const static Double_t WATER_THICK = 2.3;
const static Double_t GEL_THICK = 0.07;
const static Double_t BASE_THICK = 0.21;
const static Double_t PACK_THICK = 0.109;
const static Double_t MATERIAL_THICK[kNumberOfNinjaMaterials] = { IRON_THICK,
								WATER_THICK,
								GEL_THICK,
								BASE_THICK,
								PACK_THICK
};

///> Mass of particles in interest (MeV/c2)
const static Double_t MCS_MUON_MASS = 105.658;
const static Double_t MCS_PION_MASS = 139.570;
const static Double_t MCS_PROTON_MASS = 938.262;
const static Double_t PARTICLE_MASS[kNumberOfNinjaMcsParticles] = { MCS_MUON_MASS,
								  MCS_PION_MASS,
								  MCS_PROTON_MASS
};

///> Directions
const static Int_t MCS_DIRECTION[kNumberOfNinjaMcsDirections] = { kNinjaMcsForward,
								kNinjaMcsBackward
};

///>  Minimum number of nseg
const static Int_t MIN_NUM_SKIP = 1;
///> Maximum number of nseg
const static Int_t MAX_NUM_SKIP = 10;
///> Number of plates in one NINJA ECC
const static Int_t MAX_NUM_ECC_PLATE = 133;

///> Positional alignment accuracy in MC smearing (mm)
///> 0 : iron 1 : water
const static Double_t XY_ALIGN_ACCURACY[2] = {.3e-3, .3e-3};
const static Double_t Z_ALIGN_ACCURACY[2] = {6.e-3, 10.e-3};

///> Position determination accuracy in MC smearing (mm)
const static Double_t XY_POSITION_ACCURACY = .3e-3;
const static Double_t Z_POSITION_ACCURACY = 4.e-3;

///> Scale factor for MCS correction
const static Double_t MCS_SCALE_FACTOR = 1.0265;

///> Energy deposit function parameters for iron
const std::array<Double_t, 5> IRON_ENEDEP_FUNC_PAR = {  11.3616,
						       -41.0421,
						        63.9158,
						       -47.3527,
						        13.7393
};

///> Energy deposit function parameters for water
const std::array<Double_t, 5> WATER_ENEDEP_FUNC_PAR = {  16.677,
						        -68.1774,
						        114.435,
						        -88.6891,
						         26.3072
};

#endif
