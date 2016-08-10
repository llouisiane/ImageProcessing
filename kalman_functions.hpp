/**
 * @file   kalman_functions.hpp
 * @author Louisiane Lemaire
 * @date   2016-08
 * @see    some functions inspired from https://github.com/SwarmingSoft/ImageProcessing
 */

#ifndef __KALMAN_FUNCTIONS
#define __KALMAN_FUNCTIONS

#include <boost/numeric/ublas/matrix.hpp>
#include <list>
#include <map>
#include <set>
#include <string>
#include <vector>

#include "Vector2D.hpp"

using namespace std;

#define FLOW_WINSIZE 16 // needed to compute the optical flow
#define DELTA_T 1       // express the time in unit of the difference between 2 time steps (for velocity values)

const double PI = 3.141592653589793238463;

typedef Vector2DT<double> Vector;
typedef std::map<int, std::vector<Vector>> pos_dict_t;  // store observed positions
typedef std::map<int, std::vector<double>> ang_dict_t;  // store observed angles
typedef std::map<int, std::vector<Vector>> size_dict_t; // store observed length and width
typedef std::map<int, std::map<int, boost::numeric::ublas::matrix<double>>>matrix_dict_t;  // store states (x, y, vx, vy) or error covariance matrices

/**
 * @brief used to store the index of a prediction and the index of an observation when do assignment between predictions and observations
 */
typedef struct {
	int prediction;
	int observation;
} pred_obs_mapping_t;

/**
 * @brief store configuration parameters
 */
typedef struct {
	// working directory
	std::string working_dir;
	// subdirectory where we have the data from the segmentation ("data_exp.txt" and "data_len.txt" ) and where will be the output
	std::string data_dir;
	// input files
	std::string positions_angles_data;
	std::string length_width_data;
	// names of original pictures (we need them to compute vector fields to correct angles (we have them only modulo pi), and for visualization of the results)
	std::string pictures_filename_tpl;

	// defines the part of the images that we are using (has to be the same as in imageprocessing.cpp) // TODO: make so that imageproc; take it from here
	int subregion_x;
	int subregion_y;
	int dx;
	int dy;

	// to define which frames we want to use
	int time_start;
	int time_stop;
	int time_step; // TODO: try when we use not 1

	// minimum width for a rectangle from segmentation to be considered as an observation (enables to detect some false positive)
	double width_min;

	// initial velocity
	double v_init;
	// initial error covariance matrices P
	double initial_noise_x;
	double initial_noise_y;
	double initial_noise_vx;
	double initial_noise_vy;

	// is the velocity part of the observation model?
	bool observe_velocity;

	// covariance of position process noise, values from kalman_parameter.cpp (matrix Q)
	double weight_q_pos;
	// covariance of velocity process noise (matrix q)
	double weight_q_vel;
	// covariance of position observation noise (matrix r)
	double weight_r_pos;
	// covariance of position observation noise (matrix r)
	double weight_r_vel;

	// proportion of the edge considered as edge region (e.g.: if edge=0.1 and dx=10, we take a band of width 1 all around which is the edge region)
	double edge_region;

	// weight coefficients for energy function (used to do the assignment between predictions and observations)
	double angfactor; // angular displacement
	// rating above means that cannot be the same particle (1.)
	double rating_threshold;

	// maximum number of step of multiple prediction
	int multiple_pred_threshold; // TODO, make so that if 0 then don't do multiple prediction

	// do we print indices next to bacteria (visualization.cpp)
	bool print_index;
	// do we plot covariance ellipses around bacteria (visualization.cpp)
	bool plot_cov;

	// input kind for visualization of the results, (1: use as input output of tracking.cpp, 2: use as input output of kalman.cpp)
	int input_kind;

} configuration_parameters_t;

/**
 * @brief read parameters from configuration file
 * @return a struct containing the parameters
 * @param conf_file name of the configuration file
 */
configuration_parameters_t read_parameters(std::string conf_file);

/**
 * @brief read input file containing the observed length and width
 * @see inspired from https://github.com/SwarmingSoft/ImageProcessing
 * @param filename name of the input file
 * @param sizes store the observed length and width in a  map between time step and vector with the length and width of each rectangle at this time step
 * @param time_start first considered time step
 * @param time_stop last considered time step
 * @param time_step step
 */
void read_data_len(std::string filename, size_dict_t &sizes, int time_start, int time_stop, int time_step);

/**
 * @brief read input file containing observed positions and angles
 * @see inspired from https://github.com/SwarmingSoft/ImageProcessing
 * @param filename name of the input file
 * @param positions store the observed positions (x and y), map between time step and vector with position of each rectangle at this time step
 * @param angles store the observed angles (in radian, modulo pi) in a map between time step and vector with the angle of each rectangle at this time step
 * @param time_start first considered time step
 * @param time_stop last considered time step
 * @param time_step step
 * @param sizes contains observed length and width in a map between time step and vector with the length and width of each rectangle at this time step
 * @param width_min minimum width for a rectangle to be considered as a bacteria (and not false positive)
 */
void read_data_single(std::string filename, pos_dict_t &positions, ang_dict_t &angles, int time_start, int time_stop, int time_step, size_dict_t &sizes,
		double width_min);

/**
 * @brief when start a trajectory, initialize updated predicted state estimate and updated error covariance matrix, and increment the maximum trajectory index
 * @param t time step of the beginning of this trajectory
 * @param pos_indices index of the observed positions (in positions at this time step) we want to consider as new bacteria
 * @param positions contains observed positions (x and y) in a map between time step and vector with position of each rectangle at this time step
 * @param angles observed angles (in radian, modulo pi) in a map between time step and vector with the angle of each rectangle at this time step
 * @param updated_state_dict contains the updated state estimates, in a map of map between the time step, the index of the trajectory and the state(x, y, vx, vy)
 * @param updated_p_dict contains the updated error covariance matrices, in a map of map between the time step, the index of the trajectory and the covariance matrix
 * @param highest_index highest trajectory index so far
 * @param filenames names of the original pictures, we need them to correct the angles (observed only modulo pi) using a motion flow computation between two consecutive pictures
 * @param params struct containing the configuration parameters
 * @param new_bac the indices (in updated_state_dict and updated_p_dict) of the newly considered trajectories
 */
void add_new_bacteria(int t, list<int> pos_indices, pos_dict_t &positions, ang_dict_t &angles, matrix_dict_t &updated_state_dict, matrix_dict_t &updated_p_dict,
		int &highest_index, vector<string> filenames, configuration_parameters_t params, set<int> &new_bac);

/**
 * @brief compute angle difference modulo pi
 * @see function from https://github.com/SwarmingSoft/ImageProcessing
 * @param ang1 angle in radian
 * @param ang2 other angle in radian
 * @return the angle difference in radian
 */
double angle_difference90(double ang1, double ang2);

/**
 * @brief compute angle difference modulo 2 pi
 * @see function from https://github.com/SwarmingSoft/ImageProcessing
 * @param ang1 angle in radian
 * @param ang2 other angle in radian
 * @return the angle difference in radian
 */
double angle_difference180(double ang1, double ang2);

/**
 * @brief compute orthogonal and parallel distances knowing coordinates and angles of two rectangles
 * @see function from https://github.com/SwarmingSoft/ImageProcessing
 * @param pos1 coordinates (x and y) of the first rectangle
 * @param pos2 coordinates (x and y) of the second rectangle
 * @param ang1 coordinate of the first rectangle
 * @param ang2 coordinate of the second rectangle
 * @return the parallel and the orthogonal distance
 */
std::tuple<double, double> para_orth_distance(Vector pos1, Vector pos2, double ang1, double ang2);

/**
 * @brief compute energy cost function value between two states (e.g.:one prediction and one observation)
 * @see strongly inspired from https://github.com/SwarmingSoft/ImageProcessing
 * @param pos1 coordinates(x and y) of the first state
 * @param pos2 coordinates(x and y) of the second state
 * @param ang1 angle of the first state (in radian)
 * @param ang1 angle of the second state (in radian)
 * @param p_x error covariance in x of the prediction (if not a prediction put 0)
 * @param p_y error covariance in y of the prediction (if not a prediction put 0)
 * @param pospara_factor weight of the parallel distance
 * @param posorth_factor weight of the orthogonal distance
 * @param ang_factor weight of the angle difference
 * @param dx edge in pixel of the considered subregion
 * @return the total cost value
 */
double rate_match(Vector pos1, Vector pos2, double ang1, double ang2, double p_x, double p_y, double pos_para_factor, double pos_orth_factor, double ang_factor,
		double dx);

/**
 * @brief matrix inversion using LU factorization
 * @see inspired from https://savingyoutime.wordpress.com/2009/09/21/c-matrix-inversion-boostublas/
 * @param input matrix to invert
 * @param inverse inverse
 * @return true if it worked, false if the LU factorization failed
 */
bool invert_matrix(const boost::numeric::ublas::matrix<double>& input, boost::numeric::ublas::matrix<double>& inverse);

/**
 * @brief create the names of the input picture files described by a template, and the range to use
 * @see function from https://github.com/SwarmingSoft/ImageProcessing
 * @param tpl template for the names
 * @param start first considered time step
 * @param stop last considered time step
 * @return a vector containing all the names
 */
std::vector<std::string> create_filename_list(const char *tpl, int start, int stop);

/**
 * @brief convert a string into a boolean (the string can be: "1", "0", "true", "false")
 * @param var the string
 * @return the corresponding boolean
 */
bool string2bool(std::string var);

#endif
