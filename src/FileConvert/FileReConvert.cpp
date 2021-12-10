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
			     << " <Input momch file> <Input momrecon file> <Output momch file> <ECC id (0-8)> ";
    std::exit(1);
  }

  try {

    // input momch file を読み込む
    std::ifstream ifs(argv[1], std::ios::binary);
    Momentum_recon::Mom_chain mom_chain;
    Momentum_recon::Mom_basetrack mom_basetrack;
    std::pair<Momentum_recon::Mom_basetrack, Momentum_recon::Mom_basetrack > mom_basetrack_pair;
    int header[5];
    double mom_recon;
    int num_base, num_link;

    // pbeta recon file を読み込む
    TFile *pbeta_recon_file = new TFile(argv[2], "read");
    TTree *pbeta_recon_tree = (Ttree*)pbeta_recon_Tree->Get("tree");
    Double_t recon_pbeta;
    pbeta_recon_tree->SetBranchAddress("recon_pbeta", &recon_pbeta);

    // mom_pbeta -> recon_mom? どうするか
    // とりあえず MIP はmu, HIP はp だと思ってつけた値を入れる

    std::ofstream ofs(argv[3], std::ios::binary);

    while ( ifs.read((char*)& header, sizeof(int)*4) ) {

      mom_chain.groupid = header[0];
      mom_chain.chainid = header[1];
      mom_chain.unixtime = header[2];
      mom_chain.tracker_track_id = header[3];
      mom_chain.entry_in_daily_file = header[4];
      ifs.read((char*)& mom_recon, sizeof(double));
      mom_chain.mom_recon = mom_recon;
      ifs.read((char*)& header, sizeof(int)*2);
      num_base = header[0];
      num_link = header[1];
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

      pbeta_recon_tree->GetEntry();
      mom_chain.mom_recon = ConvertPbetaToMomentum(recon_pbeta);

      // momch binary に書き込む
      ofs.write((char*)& mom_chain.groupid, sizeof(int));
      ofs.write((char*)& mom_chain.chainid, sizeof(int));
      ofs.write((char*)& mom_chain.unixtime, sizeof(int));
      ofs.write((char*)& mom_chain.tracker_track_id, sizeof(int));
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
