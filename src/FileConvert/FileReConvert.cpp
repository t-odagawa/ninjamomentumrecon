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
#include "FileConvert.hpp"

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
			     << " <Input B2 file> <Input momrecon file> <Output momch file> <ECC id (0-8)> ";
    std::exit(1);
  }

  try {

    // B2 file を読み込む
    B2Reader reader(argv[1]);

    // 運動量recon file を読み込む
    TFile *pbeta_recon_file = new TFile(argv[2], "read");
    TTree *pbeta_recon_tree = (Ttree*)pbeta_recon_Tree->Get("tree");

    // mom_pbeta -> recon_mom? どうするか
    // とりあえず MIP はmu, HIP はp だと思ってつけた値を入れる

    std::ofstream ofs(argv[3], std::ios::binary);

    while ( reader.ReadNextSpill() > 0 ) {

      auto &spill_summary = reader.GetSpillSummary();

      if ( spill_summary.GetNumEmulsions() <= 0 ) continue;



      // momch binary に書き込む
      int num_base = mom_chain.base.size();
      int num_pair = mom_chain.base_pair.size();
      ofs.write((char*)& mom_chain.groupid, sizeof(int));
      ofs.write((char*)& mom_chain.chainid, sizeof(int));
      ofs.write((char*)& mom_chain.unixtime, sizeof(int));
      ofs.write((char*)& mom_chain.entry_in_daily_file, sizeof(int));
      ofs.write((char*)& mom_chain.mom_recon, sizeof(double));
      ofs.write((char*)& num_base, sizeof(int));
      ofs.write((char*)& num_pair, sizeof(int));
      for ( int ibase = 0; ibase < num_base; ibase++ ) {
	ofs.write((char*)& mom_chain.base.at(ibase), sizeof(Momentum_recon::Mom_basetrack));
      }
      for (int ipair = 0; ipair < num_pair; ipair++ ) {
	ofs.write((char*)& mom_chain.base_pair.at(ibase).first, sizeof(Momentum_recon::Mom_basetrack));
	ofs.write((char*)& mom_chain.base_pair.at(ibase).second, sizeof(Momentum_recon::Mom_basetrack));
      }
    }

    ofs.close();

  } catch (cosnt std::runtime_error &error) {
    BOOST_LOG_TRIVIAL(fatal) << "Runtime error : " << error.what();
    std::exit(1);
  } catch (const std::invalid_argument &error) {
    BOOST_LOG_TRIVIAL(fatal) << "Invalid argument error : " << error.what();
    std::exit(1);
  }
    
  BOOST_LOG_TRIVIAL(info) << "==========MCS ReConvert Finish==========";
  std::exit(0);

}
