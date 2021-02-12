//boost includes
#include <boost/log/core.hpp>
#include <boost/log/trivial.hpp>
#include <boost/log/expressions.hpp>

// B2 includes
#include "B2Reader.hh"

#ifdef COORDINATE_METHOD
#include "McsCoordinateMethod.hpp"
#else
#include "McsAngleMethod.hpp"
#endif

namespace loggin = boost::log;

int main (int argc, char *argv[]) {

  logging::core::get()->set_filter
    (
     logging::trivial::severity >= logging::trivial::debug
     );

  BOOST_LOG_TRIVIAL(info) << "==========NINJA Momentum Reconstruction Start==========";

  if (argc != 3) {
    BOOST_LOG_TRIVIAL(error) << "Usage : " << argv[0]
			     << " <input B2 file path> <output file path>";
    std::exit(1);
  }

  try {

    B2Reader reader(argv[1]);

    TFile *file = new TFile(argv[2], "recreate");
    TTree *tree = new TTree("tree", "ECC reconstructed momentum");

    while(reader.ReadNextSpill() > 0) {
      auto &spill_summary = reader.GetSpillSummary();

#ifdef COORDINATE_METHOD
      mcs_coordinate_method(spill_summary);
#else
      mcs_angle_method(spill_summary);
#endif

    }

  } catch (const std::runtime_error &error) {
    BOOST_LOG_TRIVIAL(fatal) << "Runtime error : " << error.what();
    std::exit(1);
  } catch (const std::invalid_argument &error) {
    BOOST_LOG_TRIVIAL(fatal) << "Invalid argument error : " << error.what();
    std::exit(1);
  }

  BOOST_LOG_TRIVIAL(info) << "==========NINJA Momentum Reconstruction Finish==========";
  std::exit(0);

}
