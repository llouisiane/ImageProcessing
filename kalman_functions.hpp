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

typedef Vector2DT<double> Vector;

typedef std::map<int, std::vector<Vector>> pos_dict_t;
typedef std::map<int, std::vector<double>> ang_dict_t;
typedef std::map<int, std::vector<Vector>> size_dict_t;
typedef std::map<int, std::map<int, boost::numeric::ublas::matrix<double>>> matrix_dict_t;

// struct used when doing the assignment between predictions and observations
typedef struct{
	int prediction;
	int observation;
} pred_obs_mapping_t;

// struct for the configuration parameters
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

	// input kind for visualisation of the results, (1: use as input output of tracking.cpp, 2: use as input output of kalman.cpp)
	int input_kind;

} configuration_parameters_t;

const double PI = 3.141592653589793238463;

configuration_parameters_t read_parameters(std::string conf_file);

void read_data_len(std::string filename, size_dict_t &sizes, int time_start, int time_stop, int time_step);
void read_data_single(std::string filename, pos_dict_t &positions, ang_dict_t &angles, int time_start, int time_stop, int time_step, size_dict_t &sizes, double width_min);

void add_new_bacteria(int t, list<int> pos_indices, pos_dict_t &positions, ang_dict_t &angles, matrix_dict_t &state_updated_dict, matrix_dict_t &updated_p_dict, int &highest_index, vector<string> filenames, configuration_parameters_t params, set<int> &new_bac);


double AngleDifference90(double a, double b);
double AngleDifference180(double a, double b);

std::tuple<double, double> ParaOrthDistance(Vector pos1, Vector pos2, double ang1, double ang2);

double rate_match(Vector pos1, Vector pos2, double ang1, double ang2, double p_x, double p_y, double pospara_factor, double
posorth_factor, double ang_factor, double dx);


bool invert_matrix(const boost::numeric::ublas::matrix<double>& input, boost::numeric::ublas::matrix<double>& inverse);

std::vector<std::string> create_filename_list(const char *tpl, int start, int stop);

bool string2bool(std::string var);

#endif
