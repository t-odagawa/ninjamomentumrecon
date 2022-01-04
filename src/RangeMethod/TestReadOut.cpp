#include <iostream>
#include <iomanip>
#include <fstream>

#include "McsClass.hpp"

int main (int argc, char* argv[]) {

  std::ifstream ifs(argv[0], std::ios::binary);
  SimpleMCBaseTrackCollection smbtc;
  while (ifs.read((char*)& smbtc, sizeof(SimpleMCBaseTrackCollection))) {

    std::cout << " Particle ID : " << smbtc.particle_id
	      << " Plate ID : " << smbtc.plate_id
	      << " ECC ID : " << smbtc.ecc_id
	      << " dE/dx 1 : " << smbtc.energy_deposit_1
	      << " dE/dx 2 : " << smbtc.energy_deposit_2
	      << " ax : " << smbtc.ax
	      << " ay : " << smbtc.ay
	      << " x : " << smbtc.x
	      << " y : " << smbtc.y
	      << " z : " << smbtc.z
	      << " Track ID : " << smbtc.track_id
	      << " Event ID :" << smbtc.event_id
	      << " File ID: " << smbtc.file_id
	      << " Weight : " << smbtc.weight << std::endl;
  }

  std::exit(0);
}
