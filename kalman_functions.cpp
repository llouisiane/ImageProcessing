/**
 * @file   kalman_functions.cpp
 * @author Louisiane Lemaire
 * @date   2016-08
 * @brief  Functions for kalman.cpp and visualization.cpp
 * @see    some functions inspired from https://github.com/SwarmingSoft/ImageProcessing
 *
 * functions documentation in kalman_functions.hpp
 */

#include "kalman_functions.hpp"

#include <boost/numeric/ublas/lu.hpp>
#include <core/mat.hpp>
#include <core/mat.inl.hpp>
#include <core/types.hpp>
#include <imgcodecs/imgcodecs_c.h>
#include <imgcodecs.hpp>
#include <imgproc.hpp>
#include <video/tracking.hpp>
#include <algorithm>
#include <cassert>
#include <cmath>
#include <cstddef>
#include <cstdio>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <utility>

using namespace std;
using namespace cv;

configuration_parameters_t read_parameters(std::string conf_file) {
	configuration_parameters_t result;
	std::ifstream infile(conf_file);
	std::string line;
	int count_param_positions = 0;
	while (std::getline(infile, line)) {
		line = line.substr(0, line.find("//"));
		if (!line.empty()) {
			switch (count_param_positions) {
			case 0:
				result.working_dir = line.c_str();
				break;
			case 1:
				result.data_dir = line.c_str();
				break;
			case 2:
				result.positions_angles_data = line.c_str();
				break;
			case 3:
				result.length_width_data = line.c_str();
				break;
			case 4:
				result.pictures_filename_tpl = line.c_str();
				break;

			case 5:
				result.subregion_x = atoi(line.c_str());
				break;
			case 6:
				result.subregion_y = atoi(line.c_str());
				break;
			case 7:
				result.dx = atoi(line.c_str());
				break;
			case 8:
				result.dy = atoi(line.c_str());
				break;

			case 9:
				result.time_start = atoi(line.c_str());
				break;
			case 10:
				result.time_stop = atoi(line.c_str());
				break;
			case 11:
				result.time_step = atoi(line.c_str());
				break;

			case 12:
				result.width_min = atof(line.c_str());
				break;

			case 13:
				result.v_init = atof(line.c_str());
				break;
			case 14:
				result.initial_noise_x = atof(line.c_str());
				break;
			case 15:
				result.initial_noise_y = atof(line.c_str());
				break;
			case 16:
				result.initial_noise_vx = atof(line.c_str());
				break;
			case 17:
				result.initial_noise_vy = atof(line.c_str());
				break;

			case 18:
				result.observe_velocity = string2bool(line.c_str());
				break;

			case 19:
				result.weight_q_pos = atof(line.c_str());
				break;
			case 20:
				result.weight_q_vel = atof(line.c_str());
				break;
			case 21:
				result.weight_r_pos = atof(line.c_str());
				break;
			case 22:
				result.weight_r_vel = atof(line.c_str());
				break;

			case 23:
				result.edge_region = atof(line.c_str());
				break;

			case 24:
				result.angfactor = atof(line.c_str());
				break;
			case 25:
				result.rating_threshold = atof(line.c_str());
				break;

			case 26:
				result.multiple_pred_threshold = atoi(line.c_str());
				break;

			case 27:
				result.print_index = string2bool(line.c_str());
				break;
			case 28:
				result.plot_cov = string2bool(line.c_str());
				break;

			case 29:
				result.input_kind = atoi(line.c_str());
				break;

			default:
				std::cout << "too many parameters" << std::endl;
			}
			count_param_positions++;

		}
	}
	return result;
}

void read_data_len(std::string filename, size_dict_t &sizes, int time_start, int time_stop, int time_step) {
	std::ifstream in;
	in.open(filename);
	assert(!in.fail() && "Could not open file");

	double length, width;
	int linenum = 0;
	int time = time_start;

	std::string firstline, line;
	std::vector<Vector> sizes_v; // vector containing the length and the width of all rectangles for one frame

	std::getline(in, firstline); // skip the first line
	while (in.good() && linenum <= time_stop) {
		std::getline(in, line);
		if (linenum == time) {
			std::stringstream linestream(line);
			sizes_v.clear();
			//if line ends with " ", we read the last value 2 times
			while (!linestream.eof()) {
				linestream >> length;
				linestream >> width;
				sizes_v.push_back(Vector(length, width));
			}
			sizes[linenum] = sizes_v;
			time += time_step;
		}
		++linenum;
	}
}

void read_data_single(std::string filename, pos_dict_t &positions, ang_dict_t &angles, int time_start, int time_stop, int time_step, size_dict_t &sizes,
		double width_min) {
	std::ifstream in;
	in.open(filename);
	assert(!in.fail() && "Could not open file");

	double x, y, angle;
	int linenum = 0;
	int time = time_start;
	int counter = 0;

	std::string firstline, line;
	std::vector<Vector> positions_v;
	std::vector<double> angles_v;

	std::getline(in, firstline); // skip first line
	while (in.good() && linenum <= time_stop) {
		std::getline(in, line);
		if (linenum == time) {
			std::stringstream linestream(line);
			positions_v.clear();
			//if line ends with " ", we read the last value 2 times
			while (!linestream.eof()) {
				linestream >> x;
				linestream >> y;
				if (sizes[linenum][counter].y > width_min) {
					positions_v.push_back(Vector(x, y));
				}
				counter++;
			}
			positions[linenum] = positions_v;
			counter = 0;
		}
		std::getline(in, line);
		if (linenum == time) {
			std::stringstream linestream(line);
			angles_v.clear();
			//if line ends with " ", we read the last value 2 times
			while (!linestream.eof()) {
				linestream >> angle;
				if (sizes[linenum][counter].y > width_min) {
					angles_v.push_back(angle);
				}
				counter++;
			}
			angles[linenum] = angles_v;
			time += time_step;

		}
		++linenum;
		counter = 0;
	}
}

void add_new_bacteria(int t, list<int> pos_indices, pos_dict_t &positions, ang_dict_t &angles, matrix_dict_t &state_updated_dict, matrix_dict_t &updated_p_dict,
		int &highest_index, vector<string> filenames, configuration_parameters_t params, set<int> &new_bac) {

	// matrix to store one updated state (x, y, vx, vy)
	boost::numeric::ublas::matrix<double> updated_state(4, 1);

	// read the two pictures
	Mat orig_1 = imread(params.working_dir + params.data_dir + filenames[t - params.time_start], CV_LOAD_IMAGE_GRAYSCALE);
	Mat orig_2 = imread(params.working_dir + params.data_dir + filenames[t + params.time_step - params.time_start], CV_LOAD_IMAGE_GRAYSCALE);

	// use only the subregion we are working on
	Rect subrect = Rect(params.subregion_x, params.subregion_y, params.dx, params.dy);
	orig_1 = Mat(orig_1, subrect);
	orig_2 = Mat(orig_2, subrect);

	// histogram equalizer
	Mat norm_1, norm_2, flow;
	equalizeHist(orig_1, norm_1);
	equalizeHist(orig_2, norm_2);

	calcOpticalFlowFarneback(norm_1, norm_2, flow, 0.5, 4, FLOW_WINSIZE, 3, 5, 1.1, 0);

	double angle_flow;  // angle of the vector field at this position (x y)
	double angle_diff;   // difference with the orientation contained in the angles dictionary at time t=t_start

	for (list<int>::iterator it = pos_indices.begin(); it != pos_indices.end(); ++it) {

		// allow new bacteria to appear only in an edge region
		/*if(t<10 || !(positions[t+ params.time_step][*it].x>params.edge_region*params.dx
		 && positions[t+ params.time_step][*it].x<params.dx-params.edge_region*dx
		 && positions[t+ params.time_step][*it].y>params.edge_region*params.dy
		 && positions[t+ params.time_step][*it].y<params.dy-params.edge_region*params.dy)){
		 */

		// because angles are observed only modulo pi, use the picture at this time step and the one at the next time step to correct them
		// by comparing with the vector field (use code from opticalflow.cpp, https://github.com/SwarmingSoft/ImageProcessing)
		// put observed values from positions and angles dictionaries in updated_state (to initialize we consider the observed values as updated state estimates)
		updated_state(0, 0) = positions[t][*it].x;
		updated_state(1, 0) = positions[t][*it].y;
		updated_state(2, 0) = cos(angles[t][*it]) * params.v_init;
		updated_state(3, 0) = sin(angles[t][*it]) * params.v_init;

		// flow vector at the position of the bacteria
		// why is the y axis upside down ?
		// TODO: check where is the problem, maybe don't need to fix here
		const Point2f& fxy = flow.at<Point2f>(updated_state(0, 0), updated_state(1, 0));
		// angle of this vector
		angle_flow = atan2(fxy.x, fxy.y);
		angle_diff = angles[t][*it] - angle_flow;
		// when absolute difference is more than 90°, correct the vx and vy in the dictionary for updated states at time t=t_start
		if ((angle_diff > PI / 2 && angle_diff < 3 * PI / 2) || angle_diff < -PI / 2) {
			updated_state(2, 0) *= -1;
			updated_state(3, 0) *= -1;
		}
		state_updated_dict[t][highest_index + 1] = updated_state;

		// initialize error covariance matrices
		updated_p_dict[t][highest_index + 1] = boost::numeric::ublas::zero_matrix<double>(4, 4);
		updated_p_dict[t][highest_index + 1](0, 0) = params.initial_noise_x;
		updated_p_dict[t][highest_index + 1](1, 1) = params.initial_noise_y;
		updated_p_dict[t][highest_index + 1](2, 2) = params.initial_noise_vx;
		updated_p_dict[t][highest_index + 1](3, 3) = params.initial_noise_vy;

		new_bac.insert(highest_index + 1);
		highest_index++;   // increment the counter
		//}
	}

}

double angle_difference90(double ang1, double ang2) {
	return std::min(mod(ang1 - ang2, PI), mod(ang2 - ang1, PI));
}

double angle_difference180(double ang1, double ang2) {
	return std::min(mod(ang1 - ang2, 2 * PI), mod(ang2 - ang1, 2 * PI));
}

std::tuple<double, double> para_orth_distance(Vector pos1, Vector pos2, double ang1, double ang2) {

	double com_distance_para, com_distance_orth;

	Vector heading1(std::cos(ang1), std::sin(ang1));
	Vector distance1 = pos2 - pos1;

	Vector parallel1 = heading1 * (distance1 * heading1);
	Vector orthogonal1 = distance1 - parallel1;

	Vector heading2(std::cos(ang2), std::sin(ang2));
	Vector distance2 = pos1 - pos2;

	Vector parallel2 = heading2 * (distance2 * heading2);
	Vector orthogonal2 = distance2 - parallel2;

	com_distance_para = (parallel1.Norm() + parallel2.Norm()) / 2.;
	com_distance_orth = (orthogonal1.Norm() + orthogonal2.Norm()) / 2.;

	return std::make_tuple(com_distance_para, com_distance_orth);
}

double rate_match(Vector pos1, Vector pos2, double ang1, double ang2, double p_x, double p_y, double pos_para_factor, double pos_orth_factor, double ang_factor,
		double dx) {

	// vertical and parallel distance
	std::tuple<double, double> com_distance;
	com_distance = para_orth_distance(pos1, pos2, ang1, ang2);

	double pospara_rating = pos_para_factor * std::get < 0 > (com_distance) / dx;
	double posorth_rating = pos_orth_factor * std::get < 1 > (com_distance) / dx;

	// when only distance
	// double com_distance;
	// com_distance = pos1.GetDistance(pos2);
	// double posrating = ((pospara_factor+posorth_factor)/2) * com_distance / (std::sqrt(double(2))*dx);

	// angle difference cost
	double ang_rating = ang_factor * angle_difference90(ang1, ang2) / (PI / 2);

	// covariance cost
	double cov_rating = 0.005 * std::sqrt(p_x * p_x + p_y * p_x);

	// return posrating + angrating + cov_rating; // when only distance
	return pospara_rating + posorth_rating + ang_rating + cov_rating;
}

bool invert_matrix(const boost::numeric::ublas::matrix<double>& input, boost::numeric::ublas::matrix<double>& inverse) {
	typedef boost::numeric::ublas::permutation_matrix<std::size_t> pmatrix;

	// create a working copy of the input
	boost::numeric::ublas::matrix<double> A(input);

	// create a permutation matrix for the LU-factorization
	pmatrix pm(A.size1());

	// perform LU-factorization
	int res = lu_factorize(A, pm);
	if (res != 0)
		return false;

	// create identity matrix "inverse" of good size
	inverse.assign(boost::numeric::ublas::identity_matrix<double>(A.size1()));

	// backsubstitute to get the inverse
	lu_substitute(A, pm, inverse);

	return true;
}

std::vector<std::string> create_filename_list(const char *tpl, int start, int stop) {
	char buffer[64];
	std::vector<std::string> ret;

	for (int i = start; i <= stop; ++i) {
		snprintf(buffer, 64, tpl, i);
		ret.push_back(std::string(buffer));
	}

	return ret;
}

bool string2bool(std::string var) {
	if (var == "true" || var == "TRUE" || var == "1") {
		return true;
	} else if (var == "false" || var == "FALSE" || var == "0") {
		return false;
	}
	return false; // TODO: how to deal with that ?
}

