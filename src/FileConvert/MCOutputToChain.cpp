// boost includes
#include <boost/log/core.hpp>
#include <boost/log/trivial.hpp>
#include <boost/log/expressions.hpp>

// B2 includes
#include <B2Reader.hh>
#include <B2Writer.hh>
#include <B2SpillSummary.hh>
#include <B2EventSummary.hh>
#include <B2VertexSummary.hh>
#include <B2EmulsionSummary.hh>

// ROOT includes

// system includes
#include <string>
#include <fstream>
#include <iostream> 
#include <vector>
#include <algorithm>
#include <map>
#include <utility>

// my includes
#include "McsClass.hpp"
#include "McsFunction.hpp"

namespace logging = boost::log;
namespace fs = boost::filesystem;

int main ( int argc, char* argv[] ) {
  
  logging::core::get()->set_filter
    (
     logging::trivial::severity >= logging::trivial::debug
     );

  BOOST_LOG_TRIVIAL(info) << "==========MC output to chain like data conversion start==========";

  if ( argc != 4 ) {
    BOOST_LOG_TRIVIAL(error) << "Usage : " << argv[0]
			     << " <Input B2 file> <Output event info file> <ECC id>";
    std::exit(1);
  }

  try {

    B2Reader reader((std::string)argv[1]);
    std::string ofilename = argv[2];

    std::vector<Momentum_recon::Event_information > ev_vec = {};

    const Int_t ecc_id = std::atoi(argv[3]);

    while ( reader.ReadNextSpill() > 0 ) {

      Momentum_recon::Event_information ev;

      auto &spill_summary = reader.GetSpillSummary();

      // MC weight
      auto it_event = spill_summary.BeginTrueEvent();
      const auto *event = it_event.Next();
      auto &primary_vertex_summary = event->GetPrimaryVertex();
      ev.weight = primary_vertex_summary.GetMCWeight();

      // Get emulsion tracks
      std::vector<const B2EmulsionSummary* > emulsions = GetTrueEmulsionTracks(spill_summary, ecc_id);

      if ( emulsions.empty() ) continue;
      std::sort(emulsions.begin(), emulsions.end(), EmulsionCompare);

      // Get true chains
      std::vector<std::vector<const B2EmulsionSummary* > > chains = GetTrueEmulsionChains(emulsions);

      // Add true information
      AddTrueChainsToEventInfo(ev, chains);

      // Smear
      SmearChains(chains);
      // Apply detection efficiency
      ApplyDetectionEfficiency(chains);
      // Fiducial volume cut
      ApplyFVCut(chains, ecc_id);
      // Linklet
      std::set<std::pair<const B2EmulsionSummary*, const B2EmusionSummary > > linklets = 
      // Multi-link erase is not necessary
      // Group reconstruction
      // 連結成分を取り出す -> それ以上上流に行けないtrack, 下流に行けないtrack を引っ張ってくる
      // 引っ張ってきたtrack をつかって再接続
      // Black だけでallowance を変えてやり直す (VPH が3倍以上で切る -> PID が違ったら切る？)
      // どちらも使ってgroup を作り直す
      // group をほどくアルゴリズムに突っ込む

      ev_vec.push_back(ev);

    }

    Momentum_recon::WriteEventInformationBin(ofilename, ev_vec);

  } catch (const std::runtime_error &error) {
    BOOST_LOG_TRIVIAL(fatal) << "Runtime error : " << error.what();
    std::exit(1);
  } catch (const std::invalid_argument &error) {
    BOOST_LOG_TRIVIAL(fatal) << "Invalid argument error : " << error.what();
    std::exit(1);
  }

  BOOST_LOG_TRIVIAL(info) << "==========MC output to chain like data conversion finish==========";
  std::exit(0);

}
