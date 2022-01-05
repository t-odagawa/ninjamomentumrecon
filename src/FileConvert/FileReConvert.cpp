// boost includes
#include <boost/log/core.hpp>
#include <boost/log/trivial.hpp>
#include <boost/log/expressions.hpp>

// B2 includes
#include <B2Reader.hh>
#include <B2Writer.hh>
#include <B2SpillSummary.hh>
#include <B2EmulsionSummary.hh>

// ROOT includes

// system includes
#include <fstream>
#include <iostream>
#include <vector>
#include <algorithm>

// my include
#include "McsConst.hpp"
#include "McsClass.hpp"
#include "McsFunction.hpp"

namespace logging = boost::log;

int main (int argc, char* argv[]) {

  logging::core::get()->set_filter
    (
     logging::trivial::severity >= logging::trivial::info
     //logging::trivial::severity >= logging::trivial::debug
     );
  
  BOOST_LOG_TRIVIAL(info) << "==========MCS ReConvert Start==========";
  
  if (argc != 5) {
    BOOST_LOG_TRIVIAL(error) << "Usage : " << argv[0]
			     << " <Input momch file> <Input momrecon file> <Output momch file> <ECC id (0-8)> ";
    std::exit(1);
  }

  try {

    // input momch file read
    std::ifstream ifs(argv[1], std::ios::binary);
    Momentum_recon::Mom_chain mom_chain;
    Momentum_recon::Mom_basetrack mom_basetrack;
    std::pair<Momentum_recon::Mom_basetrack, Momentum_recon::Mom_basetrack > mom_basetrack_pair;

    // pbeta recon file read
    TFile *pbeta_recon_file = new TFile(argv[2], "read");
    TTree *pbeta_recon_tree = (TTree*)pbeta_recon_file->Get("tree");
    Double_t recon_pbeta;
    Int_t entry_in_daily_file;
    pbeta_recon_tree->SetBranchAddress("recon_pbeta", &recon_pbeta);
    pbeta_recon_tree->SetBranchAddress("entry_in_daily_file", &entry_in_daily_file);

    std::ofstream ofs(argv[3], std::ios::binary);

    int num_base, num_link;
    int num_entry = 0;
    while ( Momentum_recon::ReadMomChainHeader(ifs, mom_chain, num_base, num_link) ) {     
      mom_chain.base.clear();
      mom_chain.base_pair.clear();
      mom_chain.base.reserve(num_base);
      mom_chain.base_pair.reserve(num_link);
      for ( int ibase = 0; ibase < num_base; ibase++ ) {
	ifs.read((char*)& mom_basetrack, sizeof(Momentum_recon::Mom_basetrack));
	mom_chain.base.push_back(mom_basetrack);
      }
      for ( int ilink = 0; ilink < num_link; ilink++ ) {
	ifs.read((char*)& mom_basetrack_pair.first, sizeof(Momentum_recon::Mom_basetrack));
	ifs.read((char*)& mom_basetrack_pair.second, sizeof(Momentum_recon::Mom_basetrack));
	mom_chain.base_pair.push_back(mom_basetrack_pair);
      }

      pbeta_recon_tree->GetEntry(num_entry);
      BOOST_LOG_TRIVIAL(debug) << "Entry in daily file in momch : " << mom_chain.entry_in_daily_file
			       << "\nEntry in daily file in root : " << entry_in_daily_file;
      mom_chain.ecc_mcs_mom = CalculateMomentumFromPBeta(recon_pbeta, MCS_MUON_MASS);

      Momentum_recon::WriteMomChainHeader(ofs, mom_chain);
      for ( int ibase = 0; ibase < num_base; ibase++ ) {
	ofs.write((char*)& mom_chain.base.at(ibase), sizeof(Momentum_recon::Mom_basetrack));
      }
      for (int ipair = 0; ipair < num_link; ipair++ ) {
	ofs.write((char*)& mom_chain.base_pair.at(ipair).first, sizeof(Momentum_recon::Mom_basetrack));
	ofs.write((char*)& mom_chain.base_pair.at(ipair).second, sizeof(Momentum_recon::Mom_basetrack));
      }

      num_entry++;
    }

    ofs.close();

  } catch (const std::runtime_error &error) {
    BOOST_LOG_TRIVIAL(fatal) << "Runtime error : " << error.what();
    std::exit(1);
  } catch (const std::invalid_argument &error) {
    BOOST_LOG_TRIVIAL(fatal) << "Invalid argument error : " << error.what();
    std::exit(1);
  }
    
  BOOST_LOG_TRIVIAL(info) << "==========MCS ReConvert Finish==========";
  std::exit(0);

}
