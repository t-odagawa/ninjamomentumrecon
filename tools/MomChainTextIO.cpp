#include <fstream>
#include <iostream>
#include <string>
#include <vector>

#include "McsClass.hpp"
#include "McsClass.cpp"

int main ( int argc, char* argv[]) {

  if ( argc != 3 ) {
    std::cout << "Usage : " << argv[0]
	      << " <input momch binary file path> <output momch text file path>" << std::endl;
    std::exit(1);
  }

  std::string ifilename = argv[1];
  std::string ofilename = argv[2];

  std::vector<Momentum_recon::Event_information > ev_vector = Momentum_recon::ReadEventInformationBin(ifilename);
  std::cout << ev_vector.size() << std::endl;
  Momentum_recon::WriteEventInformationTxt(ofilename, ev_vector);

  std::exit(0);

}
