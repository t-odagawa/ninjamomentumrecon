#include <fstream>
#include <iostream>
#include <string>
#include <vector>

#include "McsClass.cpp"

int main ( int argc, char* argv[]) {

  if ( argc != 3 ) {
    std::cout << "Usage : " << argv[0]
	      << " <input momch binary file path> <output momch text file path>" << std::endl;
    std::exit(1);
  }

  std::string ifilename = argv[1];
  std::string ofilename = argv[2];

  std::vector<Momentum_recon::Mom_chain > mom_chain_vector = Momentum_recon::ReadMomChain(ifilename);

  Momentum_recon::WriteMomChainText(ofilename, mom_chain_vector);

  std::exit(0);

}
