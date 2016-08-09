// 2016-06
// Louisiane Lemaire
// inspired from tracking.cpp (https://github.com/SwarmingSoft/ImageProcessing)

/* HOW TO USE
 * This code performs the tracking of bacteria trajectories over multiple time steps, using a Kalman filter.
 * It uses as input the output from the segmentation code that fit rectangles around bacteria on pictures at each time step:
 *
 * data_len.txt contains for each time step the length and the width of each fitted rectangle
 * header
 * length_time1_bact1 width_time1_bact1 length_time1_bact2 width_time1_bact2 ...
 * length_time2_bact1' width_time2_bact1' length_time2_bact2' width_time2_bact2' ...
 * ...
 *
 * data_exp.txt contains for each time step the angle (modulo pi) and the coordinates (x and y) of each rectangle
 * header
 * x_time1_bact1 y_time1_bact1 x_time1_bact2 y_time2_bact2 ...
 * angle_time1_bact1 angle_time1_bact2 ...
 * x_time2_bact1' y_time2_bact1' x_time2_bact2' y_time2_bact2' ...
 * angle_time2_bact1' angle_time2_bact2' ...
 * ...
 *
 * This code returns:
 * updated_states.txt, contains for each time step the index and the updated state (position and velocity) of each tracked bacteria
 * time id x y vx vy
 * ...
 *
 * updated_error_covariance.txt, contains for each time step the index and the diagonal of the updated covariance matrix of each tracked bacteria
 * time id cov_x cov_y cov_vx cov_vy
 *...
 *
 */

// because we need a lot of matrix operations in this code, we use boost (www.boost.org/doc/libs/1_61_0/libs/numeric/ublas/doc/)
#include <boost/numeric/ublas/expression_types.hpp>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/matrix_expression.hpp>
#include <core/mat.hpp>
#include <core/mat.inl.hpp>
#include <core/types.hpp>
#include <imgcodecs/imgcodecs_c.h>
#include <imgcodecs.hpp>
#include <imgproc.hpp>
#include <video/tracking.hpp>
#include <cmath>
#include <fstream>
#include <iostream>
#include <map>
#include <set>
#include <string>
#include <utility>
#include <vector>

#include "kalman_functions.hpp"
#include "Vector2D.hpp"

using namespace std;
using namespace cv;

int main(int argc, char **argv) {

	// check that the correct number of parameters (just one, the directory where the configuration file is) is given
	if (argc != 2) {
		cerr << "Wrong number of parameters. Usage: " << endl << "\t" << argv[0] << " CONFIG_FILE_DIRECTORY" << endl;
		return -1;
	}

	//////////////// 1) Initialization ////////////////

	// read parameters from the configuration file config.txt
	configuration_parameters_t params; // params contains all the configuration parameters
	params = read_parameters(string(argv[1]) + "/config.txt");

	/* input files (result of the segmentation) */

	// in_file_positions_angles contains for each frame the position (x and y) of each rectangle, and their angle (in radian modulo pi)
	string in_file_positions_angles = params.working_dir + params.data_dir + params.positions_angles_data;
	// in_file_sizes contains for each frame the length and the width of each rectangle
	string in_file_length_width = params.working_dir + params.data_dir + params.length_width_data;

	// put those observed data in dictionaries (time as a key)
	pos_dict_t positions; // positions (x and y)
	ang_dict_t angles;    // angles
	size_dict_t sizes; // lengths and widths
	read_data_len(in_file_length_width, sizes, params.time_start, params.time_stop, params.time_step);
	read_data_single(in_file_positions_angles, positions, angles, params.time_start, params.time_stop, params.time_step, sizes, params.width_min);

	// names of the original pictures
	//(those pictures are needed to compute the motion flow vector field to correct the angle values from the segmentation that are only given modulo pi)
	vector<string> filenames = create_filename_list(params.pictures_filename_tpl.c_str(), params.time_start, params.time_stop);


	/* Initialize predicted states and error covariance matrices, declare predicted states and error covariance matrices */

	// 4 by 4 zeros matrix
	boost::numeric::ublas::matrix<double> zeros = boost::numeric::ublas::zero_matrix<double>(4, 4);
	// 4 by 4 identity
	boost::numeric::ublas::identity_matrix<double> identity(4);

	// predicted states dictionary (map of map: time step, index of bacteria and predicted state(position and velocity))
	matrix_dict_t pred_state_dict;
	// matrix to store one predicted state (x, y, vx, vy)
	boost::numeric::ublas::matrix<double> pred_state(4, 1);

	// updated states dictionary (time step, index of bacteria and updated state (position and velocity))
	matrix_dict_t updated_state_dict;

	// predicted error covariance matrix dictionary
	matrix_dict_t pred_p_dict;
	// matrix to store one predicted error covariance matrix
	boost::numeric::ublas::matrix<double> pred_p(4, 4);

	matrix_dict_t updated_p_dict; // updated error covariance matrices dictionary (map of map, with times step and indices)

	int size_init = positions[params.time_start].size(); // number of detected rectangles at the first considered time step (time_start)
	int highest_index = -1;  // counter for highest used bacteria index

	// initialize updated states at time t=t_start
	// indices of the observations we want to put in updated_state_dict: all available observations at t=time_start (range from 0 to the highest one)
	list<int> pos_indices_init;
	for (int i = 0; i < size_init; i++) {
		pos_indices_init.push_back(i);
	}
	set<int> new_bac; // indices of the new bacteria at one time step (the ones in updated_state_dict[t])
	// to remove the one time step trajectory if we can't find a match at the next time step
	// put observed values in dictionary for updated steps at time t_start, correct the angles
	add_new_bacteria(params.time_start, pos_indices_init, positions, angles, updated_state_dict, updated_p_dict, highest_index, filenames, params, new_bac);

	/* compute some weight coefficients needed to compute the values of the energy function */
	//(they are used to do the assignment between predictions and observations)

	// weight for parallel distance (pos_para_factor)= 1/average length
	// weight for orthogonal distance (pos_orth_factor)= 1/average width
	int count = 0;
	double pos_para_factor = 0.;
	double pos_orth_factor = 0.;
	for (int t = params.time_start; t < params.time_stop; t += params.time_step) {
		for (Vector bact : sizes[t]) {
			pos_para_factor += bact.x;
			pos_orth_factor += bact.y;
			count++;
		}
	}
	pos_para_factor /= count;
	pos_para_factor = (1 / pos_para_factor) * params.dx;
	pos_orth_factor /= count;
	pos_orth_factor = (1 / pos_orth_factor) * params.dx;


	/* Kalman filter matrices */

	// matrix f for updated_state transition model
	boost::numeric::ublas::matrix<double> f = identity;
	f(0, 2) = DELTA_T;
	f(1, 3) = DELTA_T;

	// matrix h for observation model
	boost::numeric::ublas::matrix<double> h = identity;
	if (!params.observe_velocity) {
		h(2, 2) = 0;
		h(3, 3) = 0;
	}

	// matrix q for covariance of process noise
	boost::numeric::ublas::matrix<double> q = zeros;
	q(0, 0) = q(1, 1) = params.weight_q_pos;
	q(2, 2) = q(3, 3) = params.weight_q_vel;

	// matrix r for covariance of observation noise
	boost::numeric::ublas::matrix<double> r = zeros;
	r(0, 0) = r(1, 1) = params.weight_r_pos;
	if (!params.observe_velocity) {
		r(2, 2) = r(3, 3) = 0;
	} else {
		r(2, 2) = r(3, 3) = params.weight_r_vel;
	}

	// observed updated_state
	boost::numeric::ublas::matrix<double> z(4, 1);
	// updated_state prediction
	boost::numeric::ublas::matrix<double> x(4, 1);

	// innovation residual
	boost::numeric::ublas::matrix<double> y(4, 1);
	//innovation covariance
	boost::numeric::ublas::matrix<double> s(4, 4);
	// its inverse
	boost::numeric::ublas::matrix<double> s_inv(4, 4);

	// Kalman gain
	boost::numeric::ublas::matrix<double> k(4, 4);

	// updated updated_state estimate
	boost::numeric::ublas::matrix<double> updated_x;
	//updated estimate covariance
	boost::numeric::ublas::matrix<double> updated_p;

	/* keep track of how many multiple prediction steps we performed to stop continuing predicting
	    when bigger than the threshold (multiple_pred_threshold) */
	map<int, int> number_pred;


	/* output files */

	// updated states
	string out_file_s = params.working_dir + params.data_dir + "updated_states.txt";
	ofstream out_updated_states;
	out_updated_states.open(out_file_s);

	// updated error covariance matrices p (just the diagonal)
	string out_file_p = params.working_dir + params.data_dir + "updated_error_covariance.txt";
	ofstream out_updated_covariance;
	out_updated_covariance.open(out_file_p);


	cout << "Initialization complete." << endl;

	// loop over all frames
	for (int t = params.time_start; t < params.time_stop - params.time_step; t += params.time_step) {

		//////////////// 2) prediction of updated_state estimates at time step t+params.time_step using updated predicted state at time t ////////////////

		int prevnum = updated_state_dict[t].size(); // number of observed bacteria at time t (=number of predictions for time t+params.time_step)
		int nextnum = positions[t + params.time_step].size(); // number of observed bacteria at time t+params.time_step

		// loop over all the bacteria at a certain frame
		for (std::map<int, boost::numeric::ublas::matrix<double>>::iterator it = updated_state_dict[t].begin(); it != updated_state_dict[t].end(); ++it) {

			pred_state_dict[t + params.time_step][it->first] = prod(f, it->second);
		}

		//////////////// 3) prediction of error covariance matrices ////////////////
		for (std::map<int, boost::numeric::ublas::matrix<double>>::iterator it = updated_p_dict[t].begin(); it != updated_p_dict[t].end(); ++it) {
			pred_p = prod(f, it->second);
			pred_p_dict[t + params.time_step][it->first] = prod(pred_p, boost::numeric::ublas::trans(f)) + q;
		}

		//////////////// 4) assignment between predictions and observations, at time t+params.time_step ////////////////
		// using minimization of energy function
		// to deal with conflict use the same idea as in tracking.cpp (https://github.com/SwarmingSoft/ImageProcessing)

		cout << "Number of compared particles: " << prevnum << " " << nextnum << endl;

		multimap<double, pred_obs_mapping_t> matches_queue;

		for (map<int, boost::numeric::ublas::matrix<double>>::iterator it_pred = pred_state_dict[t + params.time_step].begin();
				it_pred != pred_state_dict[t + params.time_step].end(); ++it_pred) {
			for (int obs_particle = 0; obs_particle < nextnum; obs_particle++) {
				double rate = rate_match(Vector(it_pred->second(0, 0), it_pred->second(1, 0)), positions[t + params.time_step][obs_particle],
						atan2(it_pred->second(2, 0), it_pred->second(3, 0)), angles[t + params.time_step][obs_particle],
						pred_p_dict[t + params.time_step][it_pred->first](0, 0), pred_p_dict[t + params.time_step][it_pred->first](1, 1), pos_para_factor,
						pos_orth_factor, params.angfactor, params.dx);
				if (rate <= params.rating_threshold) {
					pred_obs_mapping_t pred_obs_mapping;
					pred_obs_mapping.prediction = it_pred->first;
					pred_obs_mapping.observation = obs_particle;
					matches_queue.insert(make_pair(rate, pred_obs_mapping));
				}
			}
		}

		// do the assignment (to deal with conflicts, assign starting by the lowest rate values, and check if observation has already been used)
		set<int> matched_predictions; // what prediction have a corresponding observation
		map<int, int> matches_result;  // map between observation and prediction
		for (multimap<double, pred_obs_mapping_t>::iterator it = matches_queue.begin(); it != matches_queue.end(); ++it) {
			if ((matches_result.find(it->second.observation) == matches_result.end())
					&& (matched_predictions.find(it->second.prediction) == matched_predictions.end())) {
				matches_result[it->second.observation] = it->second.prediction;
				matched_predictions.insert(it->second.prediction);
			}
		}

		// when we were doing multiple prediction and now have an observation available, remove the index of the trajectory from the dictionary counting the number of prediction
		for (map<int, int>::iterator it = number_pred.begin(); it != number_pred.end(); ++it) {
			if (matched_predictions.find(it->first) != matched_predictions.end()) {
				number_pred.erase(it->first);
			}
		}

		// if we assigned a new index to a bacteria and now no observation is available, remove the one time step long trajectory
		for (set<int>::iterator it = new_bac.begin(); it != new_bac.end(); ++it) {
			if (matched_predictions.find(*it) == matched_predictions.end()) {
				updated_state_dict[t].erase(*it);
				updated_p_dict[t].erase(*it);
			}
		}

		// for observations that were not close enough to any prediction, consider them as a new bacteria, with a new index
		list<int> pos_indices;  // contains indices of such observations
		for (int i = 0; i < nextnum; i++) {
			// look if this observation has not been assigned to a prediction
			if (matches_result.find(i) == matches_result.end()) {
				pos_indices.push_back(i);
			}
		}
		// indices of the new bacteria at one time step, to remove the one time step trajectory if we can't find a match at the next time step
		set<int> new_bac;
		// put observed values in dictionary for updated steps at time t+time_step, correct the angles
		add_new_bacteria(t + params.time_step, pos_indices, positions, angles, updated_state_dict, updated_p_dict, highest_index, filenames, params, new_bac);

		////////////////// 5) correction of predicted updated_state and of predicted error covariance matrices////////////////

		cout << "Number of matches: " << matches_result.size() << endl;

		for (map<int, int>::iterator m = matches_result.begin(); m != matches_result.end(); ++m) {

			int index_pred = m->second;
			int index_obs = m->first;

			// observed updated_state at indice_pred+1
			z(0, 0) = positions[t + params.time_step][index_obs].x;
			z(1, 0) = positions[t + params.time_step][index_obs].y;
			// estimate "observed" velocity
			// here use updated last position
			z(2, 0) = (positions[t + params.time_step][index_obs].x - updated_state_dict[t][index_pred](0, 0)) / DELTA_T; // vx
			z(3, 0) = (positions[t + params.time_step][index_obs].y - updated_state_dict[t][index_pred](1, 0)) / DELTA_T; // vy

			// predicted state for time t+params.time_step
			x = pred_state_dict[t + params.time_step][index_pred];

			// innovation residual
			y = z - prod(h, x);

			// innovation covariance
			pred_p = pred_p_dict[t + params.time_step][index_pred];
			s = prod(h, pred_p);
			s = prod(s, boost::numeric::ublas::trans(h)) + r;

			// Kalman gain
			k = prod(pred_p, boost::numeric::ublas::trans(h));
			invert_matrix(s, s_inv);
			k = prod(k, s_inv);

			// updated updated_state estimate
			updated_x = x + prod(k, y);

			// updated estimate covariance
			updated_p = identity - prod(k, h);
			updated_p = prod(updated_p, pred_p);

			// put in dictionary
			updated_state_dict[t + params.time_step][index_pred] = updated_x;
			updated_p_dict[t + params.time_step][index_pred] = updated_p;

		}

		// for the predictions that don't have any matching observation, continue predicting (only if number of steps below threshold)
		// when we don't find any matching observation after that number of steps, remove the trajectory
		if (params.multiple_pred_threshold > 0) {
			for (std::map<int, boost::numeric::ublas::matrix<double>>::iterator it = pred_state_dict[t + params.time_step].begin();
					it != pred_state_dict[t + params.time_step].end(); it++) {
				if (matched_predictions.find(it->first) == matched_predictions.end()) {
					// only if not on an edge region
					if ((it->second(0, 0) > params.edge_region * params.dx) && (it->second(0, 0) < params.dx - params.edge_region * params.dx)
							&& (it->second(1, 0) > params.edge_region * params.dy) && (it->second(1, 0) < params.dy - params.edge_region * params.dy)) {
						// if the first step of multiple prediction, initialize the number of steps to 1
						if (number_pred.find(it->first) == number_pred.end()) {
							number_pred[it->first] = 1;
						}
						if (number_pred[it->first] == params.multiple_pred_threshold + 1) {
							//remove trajectory
							for (int j = t + params.time_step - params.multiple_pred_threshold; j <= t + params.time_step; j++) {
								updated_state_dict[j].erase(it->first);
								// erase also covariance
								updated_p_dict[j].erase(it->first);
							}

						}
						if (number_pred[it->first] <= params.multiple_pred_threshold) {
							// put directly predicted updated_state and predicted error covariance without updating in the dictionaries for updated states and error covariance
							updated_state_dict[t + params.time_step][it->first] = it->second;
							updated_p_dict[t + params.time_step][it->first] = pred_p_dict[t + params.time_step][it->first];
							number_pred[it->first]++;

						}
					}
				}
			}
		}

		////////////////// 6) output ////////////////

		// there is a delay in the output, because when we do multiple prediction and we don't find a matching observation
		// after a certain number of steps (multiple_pred_threshol) we want to remove all what was predicted without available observation
		int output_delay = params.multiple_pred_threshold * params.time_step;
		if (t >= output_delay + params.time_start) {
			for (std::map<int, boost::numeric::ublas::matrix<double>>::iterator it = updated_state_dict[t - output_delay].begin();
					it != updated_state_dict[t - output_delay].end(); ++it) {
				// output the updated states
				out_updated_states << t - output_delay << " " << it->first << " " << it->second(0, 0) << " " << it->second(1, 0) << " " << it->second(2, 0)
						<< " " << it->second(3, 0) << endl;
				// output diagonal of updated error covariance matrices
				out_updated_covariance << t - output_delay << " " << it->first << " " << updated_p_dict[t - output_delay][it->first](0, 0) << " "
						<< updated_p_dict[t - output_delay][it->first](1, 1) << " " << updated_p_dict[t - output_delay][it->first](2, 2) << " "
						<< updated_p_dict[t - output_delay][it->first](3, 3) << endl;
			}
			// remove what we just output because we don't need anymore
			updated_state_dict.erase(t - output_delay);
			updated_p_dict.erase(t - output_delay);
		}
	}

	out_updated_states.close(); // output for updated predicted states
	out_updated_covariance.close(); // output for error covariance matrices (just output the diagonal)

	cout << "Evaluation complete." << endl;

}
