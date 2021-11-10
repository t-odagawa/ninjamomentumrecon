// boost includes
#include <boost/log/core.hpp>
#include <boost/log/trivial.hpp>
#include <boost/log/expressions.hpp>

// B2 includes
#include <B2Writer.hh>
#include <B2SpillSummary.hh>
#include <B2EmulsionSummary.hh>

// ROOT includes

// system includes
#include <fstream>
#include <iostream>
#include <vector>
#include <algorithm>

namespace logging = boost::log;

int main (int argc, char *argv[]) {

  logging::core::get()->set_filter
    (
     logging::trivial::severity >= logging::trivial::info
     //logging::trivial::severity >= logging::trivial::debug
     );

  BOOST_LOG_TRIVIAL(info) << "==========MCS Convert Start==========";

  if (argc != 3) {
    BOOST_LOG_TRIVIAL(error) << "Usage : " << argv[0]
			     << " <Input text file> <Output B2 file>";
    std::exit(1);
  }

  try {

    std::ifstream ifs(argv[1]);
    int num_track = 0;

    B2Writer writer(argv[2]);

    int group, chain, num_base, num_link;

    while ( ifs >> group >> chain >> num_base >> num_link ) {
      num_track++;
      BOOST_LOG_TRIVIAL(debug) << "Number of tracks : " << num_track;
      auto &spill_summary = writer.GetSpillSummary();

      std::vector<int> base_id_vec(num_base);
      std::vector<int> plate_vec(num_base);
      std::vector<TVector3> film_position_vec(num_base);
      std::vector<TVector3> tangent_vec(num_base);
      std::vector<Double_t> vph_up_vec(num_base), vph_down_vec(num_base);
      std::vector<Double_t> pixel_count_up_vec(num_base), pixel_count_down_vec(num_base);

      std::vector<TVector3> film_position_in_down_coordinate_vec(num_base);
      std::vector<TVector3> tangent_in_down_coordinate_vec(num_base);

      int plate, base_id;
      double tangent_x, tangent_y;
      double position_x, position_y, position_z;
      int tmp;
      int ph[2], pixel_area[2], pixel_hit[2];

      TVector3 film_position, tangent;
      TVector3 film_position_in_down_coordinate, tangent_in_down_coordinate;
      Double_t vph[2];

      BOOST_LOG_TRIVIAL(debug) << "Basetrack information get";
      for ( int ibase = 0; ibase < num_base; ibase++ ) {
	ifs >> plate >> base_id
	    >> tangent_x >> tangent_y >> position_x >> position_y >> position_z
	    >> tmp >> tmp >> tmp
	    >> ph[0] >> pixel_area[0] >> pixel_hit[0]
	    >> tmp >> tmp >> tmp
	    >> ph[1] >> pixel_area[1] >> pixel_hit[1];

	film_position.SetX(position_x / 1e3);
	film_position.SetY(position_y / 1e3);
	film_position.SetZ(position_z / 1e3);
	tangent.SetX(tangent_x);
	tangent.SetY(tangent_y);
	tangent.SetZ(1.);
	vph[0] = ph[0] % 10000; vph[1] = ph[1] % 10000;

	base_id_vec.at(ibase) = base_id;
	plate_vec.at(ibase) = plate - 1;
	film_position_vec.at(ibase) = film_position;
	tangent_vec.at(ibase) = tangent;
	vph_up_vec.at(ibase) = vph[0];
	vph_down_vec.at(ibase) = vph[1];
	pixel_count_up_vec.at(ibase) = pixel_hit[0];
	pixel_count_down_vec.at(ibase) = pixel_hit[1];

      }

      BOOST_LOG_TRIVIAL(debug) << "Linklet information get";
      int link_plate[2], link_base_id[2];
      double link_tangent_x[2], link_tangent_y[2];
      double link_position_x[2], link_position_y[2], link_position_z[2];
      for ( int ilink = 0; ilink < num_link; ilink++ ) {
	ifs >> link_plate[0] >> link_base_id[0]
	    >> link_plate[1] >> link_base_id[1]
	    >> link_tangent_x[0] >> link_tangent_y[0]
	    >> link_position_x[0] >> link_position_y[0] >> link_position_z[0]
	    >> link_tangent_x[1] >> link_tangent_y[1]
	    >> link_position_x[1] >> link_position_y[1] >> link_position_z[1];

	film_position_in_down_coordinate.SetX(link_position_x[1] / 1e3);
	film_position_in_down_coordinate.SetY(link_position_y[1] / 1e3);
	film_position_in_down_coordinate.SetZ(link_position_z[1] / 1e3);
	tangent_in_down_coordinate.SetX(link_tangent_x[1]);
	tangent_in_down_coordinate.SetY(link_tangent_y[1]);
	tangent_in_down_coordinate.SetZ(1.);
	auto find_itr = std::find(base_id_vec.begin(), base_id_vec.end(), link_base_id[1]);
	int index = find_itr - base_id_vec.begin();
	film_position_in_down_coordinate_vec.at(index) = film_position_in_down_coordinate;
	tangent_in_down_coordinate_vec.at(index) = tangent_in_down_coordinate;
      }

      std::reverse(plate_vec.begin(), plate_vec.end());
      std::reverse(film_position_vec.begin(), film_position_vec.end());
      std::reverse(tangent_vec.begin(), tangent_vec.end());
      std::reverse(vph_up_vec.begin(), vph_up_vec.end());
      std::reverse(vph_down_vec.begin(), vph_down_vec.end());
      std::reverse(pixel_count_up_vec.begin(), pixel_count_up_vec.end());
      std::reverse(pixel_count_down_vec.begin(), pixel_count_down_vec.end());
      std::reverse(film_position_in_down_coordinate_vec.begin(),
		   film_position_in_down_coordinate_vec.end());
      std::reverse(tangent_in_down_coordinate_vec.begin(),
		   tangent_in_down_coordinate_vec.end());

      for ( int ibase = 0; ibase < num_base; ibase++ ) {
	auto &emulsion_summary = spill_summary.AddEmulsion();
	emulsion_summary.SetEmulsionTrackId((UInt_t)base_id_vec.at(ibase));
	emulsion_summary.SetParentTrackId(chain);
	emulsion_summary.SetFilmPosition(film_position_vec.at(ibase));
	emulsion_summary.SetTangent(tangent_vec.at(ibase));
	emulsion_summary.SetFilmPositionInDownCoordinate(film_position_in_down_coordinate_vec.at(ibase));
	emulsion_summary.SetTangentInDownCoordinate(tangent_in_down_coordinate_vec.at(ibase));
	emulsion_summary.SetVphUp(vph_up_vec.at(ibase));
	emulsion_summary.SetVphDown(vph_down_vec.at(ibase));
	emulsion_summary.SetPixelCountUp(pixel_count_up_vec.at(ibase));
	emulsion_summary.SetPixelCountDown(pixel_count_down_vec.at(ibase));
	emulsion_summary.SetFilmType(B2EmulsionType::kECC);
	emulsion_summary.SetEcc(4);
	emulsion_summary.SetPlate(plate_vec.at(ibase));
      }

      BOOST_LOG_TRIVIAL(debug) << "Total number of tracks : " << num_track;
      writer.Fill();

    }


  } catch (const std::runtime_error &error) {
    BOOST_LOG_TRIVIAL(fatal) << "Runtime error : " << error.what();
    std::exit(1);
  } catch (const std::invalid_argument &error) {
    BOOST_LOG_TRIVIAL(fatal) << "Invalid argument error : " << error.what();
    std::exit(1);
  }

  BOOST_LOG_TRIVIAL(info) << "==========MCS Convert Finish==========";
  std::exit(0);

}
