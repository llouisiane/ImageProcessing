// 2016-06
// inspired from tracking.cpp of Markus and Janis

/* HOWTO USE
 * TODO
 */

#include <iostream>
#include <fstream>
#include <map>
#include <tuple>
#include <set>
#include <stdlib.h>
#include <string>

// for matrices
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/io.hpp>
#include <boost/numeric/ublas/operation.hpp>
#include <boost/numeric/ublas/operation_blocked.hpp>
#include <boost/numeric/ublas/lu.hpp>

#include "core.hpp"
#include "highgui.hpp"
#include "imgproc.hpp"
#include "imgcodecs.hpp"
#include "imgcodecs/imgcodecs_c.h"

#include "video/tracking.hpp"


#include "Vector2D.hpp"
#include "constantsImgPro.hpp"

typedef unsigned int uindex;


const real sqrt2 = std::sqrt(real(2));

std::string read_data_len(std::string filename, Sizedict &sizes, std::vector<unsigned int> &times)
{
    std::ifstream in;
    in.open(filename);
    assert(!in.fail() && "Could not open file");

    real size1, size2;
    unsigned int linenum = 1;
    unsigned int time = 0;

    std::string firstline, line;
    std::vector<Vector> sizesv; //()

    std::getline(in, firstline); //skip first line
    while (in.good())
    {
        //std::cout << time << " " << linenum << std::endl;
        std::getline(in, line);
        if (linenum == times[time])
        {
            std::stringstream linestream(line);
            sizesv.clear();
            //if line ends with " ", we read the last value 2 times
            while (!linestream.eof())
            {
                linestream >> size1;
                linestream >> size2;
                sizesv.push_back(Vector(size1, size2));
            }
            sizes[linenum] = sizesv;
            ++time;
        }
        ++linenum;
    }
    assert(time == times.size());
    return firstline;
}


real AngleDifference90(real a, real b)
{
    return std::min(mod(a - b, PI), mod(b - a, PI));
}

real AngleDifference180(real a, real b)
{
    return std::min(mod(a - b, 2*PI), mod(b - a, 2*PI));
}

std::tuple<real, real> ParaOrthDistance(Vector pos1, Vector pos2, real ang1, real ang2)
{
    real com_distance_para, com_distance_orth;

    Vector heading1(std::cos(ang1), std::sin(ang1));
    Vector distance1 = pos2 - pos1;

    Vector parallel1 = heading1*(distance1 * heading1);
    Vector orthogonal1 = distance1 - parallel1;

    Vector heading2(std::cos(ang2), std::sin(ang2));
    Vector distance2 = pos1 - pos2;

    Vector parallel2 = heading2*(distance2 * heading2);
    Vector orthogonal2 = distance2 - parallel2;

    com_distance_para = (parallel1.Norm() + parallel2.Norm()) / 2.;
    com_distance_orth = (orthogonal1.Norm() + orthogonal2.Norm()) / 2.;

    return std::make_tuple(com_distance_para, com_distance_orth);
}

void ParaOrthDistance_test(void)
{
    auto dist = ParaOrthDistance(Vector(1, 1), Vector(1, 2), 0, 0);
    assert_almost_equal(std::get<0>(dist), real(0));
    assert_almost_equal(std::get<1>(dist), real(1));

    dist = ParaOrthDistance(Vector(1, 1), Vector(1, 2), PI/2, PI/2);
    assert_almost_equal(std::get<0>(dist), real(1));
    assert_almost_equal(std::get<1>(dist), real(0));

    dist = ParaOrthDistance(Vector(1, 1), Vector(1, 2), 0, PI/2);
    assert_almost_equal(std::get<0>(dist), real(0.5));
    assert_almost_equal(std::get<1>(dist), real(0.5));

    dist = ParaOrthDistance(Vector(1, 1), Vector(1, 2), PI/4, PI/4);
    assert_almost_equal(std::get<0>(dist), real(1./std::sqrt(2)));
    assert_almost_equal(std::get<1>(dist), real(1./std::sqrt(2)));

    assert_almost_equal(0.1, 0.); //fails ofcourse
}

//imagesize in px
real rate_match(Vector pos1, Vector pos2, real ang1, real ang2, Vector size1, Vector size2, real p_x, real p_y, real pospara_factor, real posorth_factor, real ang_factor, real area_factor, real dx, bool periodic=false)
{
    //real com_distance;
    std::tuple<real, real> com_distance;
    if (!periodic)
    {
        //distance only
        //com_distance = pos1.GetDistance(pos2);

        //vertical and parallel distance
        com_distance = ParaOrthDistance(pos1, pos2, ang1, ang2);
    }
    else
    {
        //com_distance = distance_periodic(pos1, pos2, L);
        //not yet implemented!!!
        assert(false);
    }

    //real posrating = posfactor * com_distance / (sqrt2*L); //sqrt2 is not important
    real pospara_rating = pospara_factor * std::get<0>(com_distance) / dx;
    real posorth_rating = posorth_factor * std::get<1>(com_distance) / dx;

    //real angabsdiff = std::abs(ang1-ang2);
    //real angrating = angfactor * std::min(angabsdiff, PI - angabsdiff) / PI; //CAN BE NEGATIVE; FIX
    real ang_rating = ang_factor * AngleDifference90(ang1, ang2) / (PI/2); //not very very well tested...


    real area1 = size1.x * size1.y;
    real area2 = size2.x * size2.y;

    real area_rating;

    if (area1 >= area2)
    {
        area_rating = area_factor * (area1 - area2) / area1;
    }
    else
    {
        area_rating = area_factor * (area2 - area1)  /area2;
    }

    // covariance
    real cov_rating=0.0025*(p_x+p_y);

    //return +posrating + angrating + arearating;
    return pospara_rating + posorth_rating + ang_rating + area_rating + cov_rating;
}



real calcKal(unsigned int t, Posdict &positions, Matrixdict &pred, Matrixdict &p_pred_dict, Angdict &angles, Sizedict &sizes, auto particle, uindex other_particle, real L, bool periodic, real pospara_factor, real posorth_factor, real angfactor, real areafactor)
{
    return rate_match(Vector(particle.second(0,0), particle.second(1,0)), positions[t][other_particle], atan2(particle.second(2,0), particle.second(3,0)), angles[t][other_particle], sizes[t-1][particle.first], sizes[t][other_particle], p_pred_dict[t][particle.first](0,0), p_pred_dict[t][particle.first](1,1), pospara_factor, posorth_factor, angfactor, areafactor, L, periodic);
}
template <typename T> std::vector<T> range(T n)
{
    std::vector<T> ret(n);
    for (T i = 0; i < n; ++i)
    {
        ret[i] = i;
    }
    return ret;
}


// function for matrix inversion using LU factorisation (inspired from https://savingyoutime.wordpress.com/2009/09/21/c-matrix-inversion-boostublas/)
bool InvertMatrix(const boost::numeric::ublas::matrix<real>& input, boost::numeric::ublas::matrix<real>& inverse)
{
	typedef boost::numeric::ublas::permutation_matrix<std::size_t> pmatrix;

	// create a working copy of the input
	boost::numeric::ublas::matrix<real> A(input);

	// create a permutation matrix for the LU-factorization
	pmatrix pm(A.size1());
	//std::cout<<"pm "<<pm<<std::endl;

	// perform LU-factorization
	int res = lu_factorize(A, pm);
	if (res != 0)
		return false;

	//std::cout<<"A "<<A<<std::endl;

	// create identity matrix "inverse" of good size
	inverse.assign(boost::numeric::ublas::identity_matrix<real> (A.size1()));
	//std::cout<<"inverse"<<inverse<<std::endl;


	// backsubstitute to get the inverse
	lu_substitute(A, pm, inverse);

	return true;
}

//function to create a vector containing all the names of the files described by ...
std::vector<std::string> CreateFilenameList(const char *tpl, unsigned int start, unsigned int stop)
{
    char buffer[64];
    std::vector<std::string> ret;

    for (unsigned int i = start; i < stop; ++i)
    {
        snprintf(buffer, 64, tpl, i);
        ret.push_back(std::string(buffer));
    }

    return ret;
}


bool string2bool(std::string var){
    if(var == "true" || var == "TRUE" || var=="1" ){
  	    return true;
    }
    else if(var == "false" || var == "FALSE" || var=="0" ){
        return false;
    }
    return false; // TODO: how to deal with that ?
}


void readParameters(std::string conf_file, int& X, int& Y, int& DX, int& DY, bool& plot_cov, unsigned int& time_start, unsigned int& time_stop, unsigned int& time_step, std::string& data_directory, std::string& filenametpl, std::string& data_positions_angles, std::string& data_sizes, real& v_mean, real& delta_t, bool& observe_velocity, real& initialNoise_x, real& initialNoise_y, real& initialNoise_vx, real& initialNoise_vy, real& weight_q_pos, real& weight_q_vel, real& weight_r, real& edge_region, real& angfactor, real& areafactor, real& rating_threshold, bool& print_index, int& multiple_pred_threshold, bool& multiple_pred, real& width_min){
    std::ifstream infile(conf_file);
    std::string line;
    int count_param_positions=0;
    while(std::getline(infile,line)){
    	line=line.substr(0, line.find("//"));
    	if(!line.empty()){
    		switch(count_param_positions){
    		    case 0: X=atoi(line.c_str());
    		            break;
    		    case 1: Y=atoi (line.c_str());
    		            break;
    		    case 2: DX=atoi (line.c_str());
    		        	break;
    		    case 3: DY=atoi (line.c_str());
    		        	break;
    		    case 4: plot_cov=string2bool(line.c_str());
    		        	break;
    		    case 5: time_start=atoi(line.c_str());
    		        	break;
    		    case 6: time_stop=atoi(line.c_str());
    		        	break;
    		    case 7: time_step=atoi(line.c_str());
    		        	break;
                //TODO: deal with fact that there can be a some spaces at the end
    		    case 8: data_directory=line.c_str();
    		        	break;
    		    case 9: filenametpl=line.c_str();
    		        	break;
    		    case 10: data_positions_angles=line.c_str();
    		        	 break;
    		    case 11: data_sizes=line.c_str();
    		        	 break;
    		    case 12: v_mean=atof(line.c_str());
    		        	 break;
    		    case 13: delta_t=atof(line.c_str());
    		        	 break;
    		    case 14: observe_velocity=string2bool(line.c_str());
    		        	 break;
    		    case 15: initialNoise_x=atof(line.c_str());
    		        	 break;
    		    case 16: initialNoise_y=atof(line.c_str());
    		        	 break;
    		    case 17: initialNoise_vx=atof(line.c_str());
    		        	 break;
    		    case 18: initialNoise_vy=atof(line.c_str());
    		        	 break;
    		    case 19: weight_q_pos=atof(line.c_str());
    		        	 break;
    		    case 20: weight_q_vel=atof(line.c_str());
    		        	 break;
    		    case 21: weight_r=atof(line.c_str());
    		        	 break;
    		    case 22: edge_region=atof(line.c_str());
    		        	 break;
    		    case 23: angfactor=atof(line.c_str());
    		        	 break;
    		    case 24: areafactor=atof(line.c_str());
    		        	break;
    		    case 25: rating_threshold=atof(line.c_str());
    		        	 break;
    		    case 26: print_index=string2bool(line.c_str());
    		        	break;
    		    case 27: multiple_pred_threshold=atoi(line.c_str());
    		             break;
    		    case 28: multiple_pred=string2bool(line.c_str());
    		             break;
    		    case 29: width_min=atof(line.c_str());
    		             break;
    		    case 30: break;
    		    default: std::cout<<"too many parameters"<<std::endl;
    		}
            count_param_positions++;
    	}
    }
}

int main(int argc, char **argv)
{

	if(argc!=2){
		std::cerr<<"Wrong number of parameter. Usage: "<<std::endl<<"\t"<<argv[0]<<" DATA_DIRECTORY"<<std::endl;
		return -1;
	}

	//////////////// 1) Initialization ////////////////

	// read arguments from configuration file config.txt
    int subregion_x, subregion_y, dx, dy; // defines the part of the images that we are using (has to be the same as in imageprocessing.cpp)
	bool plot_cov;    // do we plot covariance ellipses
    unsigned int time_start, time_stop, time_step; // to define which frames we want to use
    // input files
    std::string dir=""; // working directory // TODO: how to deal with an empty parameter
    std::string data_directory; // subdirectory where we put the data and where will be the output
    std::string filenametpl;    // names of original pictures (we need them to compute vector fields to correct angles, and for visualization of the results)
    std::string data_positions_angles;
    std::string data_sizes;
    real v_mean;
    // time difference between 2 consecutive frames, what is it? say 1 (time in unit of delta t)
    real delta_t;
    bool observe_velocity;
	// initialize error covariance matrices P
    real initialNoise_x;
    real initialNoise_y;
    real initialNoise_vx;
    real initialNoise_vy;
    real weight_q_pos; // covariance of position process noise // values from kalman_parameter.cpp
    real weight_q_vel; // covariance of velocity process noise
    real weight_r;         // covariance of observation noise
    real edge_region;      // proportion of the edge considered as edge (e.g.: if edge=0.1 and dx=10, we take a band of width 1 all around which is the edge)
    // weight coefficients for energy function (used to do the assignment between predictions and observations)
    real angfactor;		// for angular displacement
    real areafactor;		// for area change
    real rating_threshold; // (1.)  rating above means that cannot be the same particle
    bool print_index;
    int multiple_pred_threshold;
    bool multiple_pred;
    real width_min; // minimum width for a rectangle from the segmentation to be considered as a bacteria
    readParameters(std::string(argv[1]) + "/config.txt", subregion_x, subregion_y, dx, dy, plot_cov, time_start, time_stop, time_step, data_directory, filenametpl, data_positions_angles, data_sizes, v_mean, delta_t, observe_velocity, initialNoise_x, initialNoise_y, initialNoise_vx, initialNoise_vy, weight_q_pos, weight_q_vel, weight_r, edge_region, angfactor, areafactor, rating_threshold, print_index, multiple_pred_threshold, multiple_pred, width_min);

    // times is a vector containing all the indices of the frames we want to use (time_stop not included)
    std::vector<unsigned int> times = MakeTimes(time_start, time_stop, time_step);

    std::vector<std::string> filenames = CreateFilenameList(filenametpl.c_str(), time_start-1, time_stop);
    // TODO: when time_step is not 1 ?

    // in_file_positions_angles contains x and y positions for each box, and angle of orientation (in radians, +/- pi) and that for each frame
    // in_file_sizes contains length and width of each box, and that for each frame
    std::string in_file_positions_angles = dir + data_directory + data_positions_angles;
    std::string in_file_sizes = dir + data_directory + data_sizes;

    // put observed data in dictionaries: one for the positions (x and y), one for the angles and one for the sizes (length and width) of the fitted rectangles
    Posdict positions;
    Angdict angles;
    Sizedict sizes;
    std::string header = read_data_len(in_file_sizes, sizes, times);
    read_data_single(in_file_positions_angles, positions, angles, times, sizes, width_min);

    // because in this program we need a lot of matrix operations use a library: boost
    // (www.boost.org/doc/libs/1_61_0/libs/numeric/ublas/doc/)

    // dictionary for predictions of states (dictionary of dictionary: time step, identity of bacteria and state)
    Matrixdict state_pred_dict;
    // vector of one predicted state (x, y, vx, vy)
    boost::numeric::ublas::matrix<real> state_pred(4,1);

    // dictionary for updated states (time, identity of bacteria and state)
    Matrixdict state_updated_dict;

    // 4 by 4 matrix with only zeros
	boost::numeric::ublas::matrix<real> zeros = boost::numeric::ublas::zero_matrix<real>(4,4);

	// identity 4
    boost::numeric::ublas::identity_matrix<real> identity(4);

	Matrixdict p_dict; // dictionary of dictionaries, with times step and identities, for updated error covariance matrices
    unsigned int sizeInit=positions[1].size();
	for (unsigned int i=0; i<sizeInit; i++){
		p_dict[1][i]=zeros;
		p_dict[1][i](0,0)=initialNoise_x;
		p_dict[1][i](1,1)=initialNoise_y;
		p_dict[1][i](2,2)=initialNoise_vx;
		p_dict[1][i](3,3)=initialNoise_vy;
	}

	// predicted error covariance matrix dictionary
	Matrixdict p_pred_dict;
	// one error covariance matrix
    boost::numeric::ublas::matrix<real> p_pred(4,4);

	// matrix f for state transition model
	boost::numeric::ublas::matrix<real> f=identity;
    f(0,2)=delta_t;
    f(1,3)=delta_t;

    // matrix h for observation model
    boost::numeric::ublas::matrix<real> h=identity;
    if(!observe_velocity){
        h(2,2)=0;
    	h(3,3)=0;
    }

    // matrix q for covariance of process noise
    boost::numeric::ublas::matrix<real> q=zeros;
    q(0,0)=q(1,1)=weight_q_pos;
  	q(2,2)=q(3,3)=weight_q_vel;

    // matrix r for covariance of observation noise
    boost::numeric::ublas::matrix<real> r=zeros;
    r(0,0)=r(1,1)=weight_r;
    if(!observe_velocity){
    	r(2,2)=r(3,3)=0;
    }
    else{
    	r(2,2)=r(3,3)=weight_r;
    }

    // counter for highest used bacteria index;
    int highest_index=-1;
    // initialize updated states at time t=1
    boost::numeric::ublas::matrix<real> state(4,1);
    for(unsigned int i=0; i<sizeInit; i++){
    	state(0,0)=positions[1][i].x;
    	state(1,0)=positions[1][i].y;
    	state(2,0)=std::cos(angles[1][i])*v_mean;
    	state(3,0)=std::sin(angles[1][i])*v_mean;
    	state_updated_dict[1][i]=state;
    	highest_index ++;   // increment the counter
    }

    // using the two first frames, correct the orientation (because +/- pi) by comparing with the vector field
    // use code from opticalflow.cpp

    // two first pictures
    cv::Mat orig_1 = cv::imread(data_directory +filenames[0], CV_LOAD_IMAGE_GRAYSCALE);
    cv::Mat orig_2 = cv::imread(data_directory +filenames[1], CV_LOAD_IMAGE_GRAYSCALE);

    // use only the subregion we are working on
    cv::Rect subrect = cv::Rect(subregion_x, subregion_y, dx, dy);
	orig_1=cv::Mat(orig_1, subrect);
	orig_2=cv::Mat(orig_2, subrect);;

    // histogram equalizer
    cv::Mat norm_1, norm_2, flow;
    unsigned int winsize=16;   // ?
    cv::equalizeHist(orig_1, norm_1);
    cv::equalizeHist(orig_2, norm_2);

    // compute vector field
    cv::calcOpticalFlowFarneback(norm_1, norm_2, flow, 0.5, 4, winsize, 3, 5, 1.1, 0);


    // do correction
    real angle_flow;  // angle of the vector field at this position (x y)
    real angle_dif;   // difference with the orientation contained in the angles dictionary at time t=1
    // loop for all the bacteria found at time t=1
    for (auto state_entry:state_updated_dict[1]){
    	// flow vector at the position of the bacteria
    	// why is the y axis upside down ?
    	// TODO: check where is the problem, maybe don't need to fix here
    	const cv::Point2f& fxy=flow.at<cv::Point2f>(state_entry.second(0,0), state_entry.second(1,0));
    	// angle of this vector
        angle_flow=atan2(fxy.x,fxy.y);
        angle_dif=angles[1][state_entry.first]-angle_flow;
        // when absolute difference is more than 90°, correct the vx and vy in the dictionary for updated states at time t=1
        if ((angle_dif> PI/2 && angle_dif< 3*PI/2) || angle_dif<-PI/2){
            state_entry.second(2,0)*= -1;
            state_entry.second(3,0)*= -1;
        }
    }

    // weight coefficients for energy function (used to do the assignment between predictions and observations)
    // for parallel distance: 1/average length  (10.)
    int count=0;
    real pospara_factor = 0.;
    for (unsigned int t: times){
    	for (Vector bact: sizes[t]){
    	    pospara_factor+=bact.x;
    	    count++;
    	}
    }
    pospara_factor/=count;
    pospara_factor=(1/pospara_factor)*dx;

    //TODO collapse in one loop
    // for orthogonal distance: 1/average width (50.)
    count=0;
    real posorth_factor = 0.;
    for (unsigned int t: times){
    	for (Vector bact: sizes[t]){
    	    posorth_factor+=bact.y;
    	    count++;;
    	}
    }
    posorth_factor/=count;
    posorth_factor=(1/posorth_factor)*dx;


    // periodic boundary condition or not
    bool periodic_boundary_cond=false; // true not implemented yet

    // observed state
    boost::numeric::ublas::matrix<real> z(4,1);
    // state prediction
    boost::numeric::ublas::matrix<real> x(4,1);

    // innovation residual
    boost::numeric::ublas::matrix<real> y(4,1);
    //innovation covariance
    boost::numeric::ublas::matrix<real> s(4,4);
    // its inverse
    boost::numeric::ublas::matrix<real> s_inv(4,4);

    // Kalman gain
    boost::numeric::ublas::matrix<real> k(4,4);

    // updated state estimate
    boost::numeric::ublas::matrix<real> updated_x;
    //updated estimate covariance
    boost::numeric::ublas::matrix<real> updated_p;

    // updated states
    std::string out_file_s = data_directory + "updated_states.txt";
    std::ofstream out1;
    out1.open(out_file_s);

    // updated error covariance matrices p (just the diagonal)
    std::string out_file_p = data_directory + "updated_error_covariance.txt";
    std::ofstream out2;
    out2.open(out_file_p);

    // to stop continuing predicting when no obervation available for too much time
    std::map<int, int> number_pred;

    std::set<int> new_bac;

    std::cout << "Initialisation complete." << std::endl;


    // loop for all frames
    for (unsigned int t:times){

    	//////////////// 2) Prediction of state estimate at time step t+1 using updated predicted state at time t ////////////////


        unsigned int prevnum = state_updated_dict[t].size();     // number of observed bacteria at time t (=number of predictions for time t+1)
        unsigned int nextnum = positions[t+1].size();	// number of observed bacteria at time t+1

// #pragma unroll
        // loop for all the bacteria at a certain frame
        for (auto i: state_updated_dict[t]){
        // can do that in parallel ??
    	    state_pred_dict[t+1][i.first]=prod(f,i.second);
        }

        //////////////// 3) prediction of error covariance matrices ////////////////
        for(auto i:p_dict[t]){
            p_pred=prod(f, i.second);
            p_pred_dict[t+1][i.first]=prod(p_pred,boost::numeric::ublas::trans(f))+q;
        }


        //////////////// 4) assignment between predictions and observations, at time t+1 ////////////////
        // using minimization of energy function
        // to deal with conflict use the same algorithm as in tracking.cpp (copy code from there)

        std::cout << "Number of compared particles: " << prevnum << " " << nextnum << std::endl;

        std::vector<uindex> nextbox = range(nextnum);

        std::vector<std::vector<uindex>> neighbouring_nextboxes; // including self, needs to be a function
        neighbouring_nextboxes.push_back(nextbox);

        std::map<uindex, std::map<real,uindex >> rates;


// TODO: make less ugly

        for (auto particle : state_pred_dict[t+1])
        {
            for (auto nbox : neighbouring_nextboxes)
            {
                for (uindex other_particle : nbox)
                {
                    //map should contain only tracking partners with rating below threshold
                	// TODO: declare the function somewhere, change the name
                    real rate = calcKal(t+1, positions, state_pred_dict, p_pred_dict, angles, sizes, particle, other_particle, dx, periodic_boundary_cond, pospara_factor, posorth_factor, angfactor, areafactor);
                    if (rate <= rating_threshold)
                    {
                        rates[particle.first][rate] = other_particle;
                    }
                }
            }
        }


        // use the computed rates to do assignment between predicted states and observed states, at t+1
        std::multimap<real, std::tuple<uindex, uindex, uindex>> matches;
        std::multimap<real, std::tuple<uindex, uindex, uindex>> newmatches;

        for (uindex i = 0; i < rates.size(); ++i)
        {
            auto p = rates[i].begin();
            if (p != rates[i].end()) //threshold not to small
            {
                matches.insert(std::make_pair(p->first, std::make_tuple(i, p->second, 0)));
            }
        }

        std::set<uindex> used_matches;
        std::set<uindex> matched_observations;  // to know what observation were used, to give new identity to the ones that were not used
        std::set<uindex> matched_predictions;  // same for predictions, to continue predicting without updated for those
        bool tracking_changes = true;
        while (tracking_changes) //keep iterating
        {
            tracking_changes = false;

            for (auto m : matches)
            {
                auto nextmatch = std::get<1>(m.second);
                auto nextrate = m.first;
                auto nextindex = std::get<2>(m.second);

                auto ratelist = rates[std::get<0>(m.second)];
                auto rateiter = ratelist.begin();
                std::advance(rateiter, nextindex);

                while (true)
                {
                    if (used_matches.insert(nextmatch).second == false) // already in set
                    {
                        //finished all particles in rating map rates[p] -> assume
                        //the first detection was false positive
                        //(next match would be below rating_threshold)
                        if (rateiter == ratelist.end())
                        {
                            break;
                        }

                        tracking_changes = true;

                        nextmatch = rateiter->second;
                        nextrate = rateiter->first;
                        ++nextindex;
                        ++rateiter;
                    }
                    else
                    {
                        newmatches.insert(std::make_pair(nextrate, std::make_tuple(std::get<0>(m.second), nextmatch, nextindex)));
                        break;
                    }
                }
            }
            matches = newmatches;
            matched_observations=used_matches;
            used_matches.clear();

            newmatches.clear(); //hm really? ka, apparently this fixes everything (at least something)
        }

        for(auto m: matches){
            auto indice_pred=std::get<0>(m.second);
            matched_predictions.insert(indice_pred);
        }

        // when we were doing multiple prediction and now have an observation available, remove the index of the trajectory from the dictionary counting the number of prediction
        for(auto i: number_pred){
        	//
        	if(matched_predictions.find(i.first)!=matched_predictions.end()){
        	    number_pred.erase(i.first);
        	}

        }

        // if we assigned a new index to a bacteria and now no observation is available, remove the one time step long trajectory

        for(std::set<int>::iterator it = new_bac.begin(); it != new_bac.end(); ++it){
        	if(matched_predictions.find(*it)==matched_predictions.end()){
        		state_updated_dict[t].erase(*it);
        		p_dict[t].erase(*it);

        	}
        }


        // what to do with observation that were not close enough to any prediction ? Say it's a new bacteria ??

        // two last pictures
        cv::Mat img_1 = cv::imread(data_directory +filenames[t-1], CV_LOAD_IMAGE_GRAYSCALE);
        cv::Mat img_2 = cv::imread(data_directory +filenames[t], CV_LOAD_IMAGE_GRAYSCALE);

        // use only the subregion we are working on
        cv::Rect subrect = cv::Rect(subregion_x, subregion_y, dx, dy);
    	cv::Mat subImage_img_1(img_1, subrect);
    	img_1=subImage_img_1;
    	cv::Mat subImage_img_2(img_2, subrect);
    	img_2=subImage_img_1;

        // histogram equalizer
        unsigned int winsize=16;   // ?
        cv::equalizeHist(img_1, norm_1);
        cv::equalizeHist(img_2, norm_2);

        // compute vector field
        cv::calcOpticalFlowFarneback(norm_1, norm_2, flow, 0.5, 4, winsize, 3, 5, 1.1, 0);

        //
        new_bac.clear();

        for (unsigned int i=0; i<nextnum; i++){
        	// look if this observation has been assigned to a prediction
            if (matched_observations.find(i)==matched_observations.end()){
	    	    //if(t<10 || !(positions[t+1][i].x>edge_region*DX && positions[t+1][i].x<DX-edge_region*DX && positions[t+1][i].y>edge_region*DY && positions[t+1][i].y<DY-edge_region*DY)){
                	new_bac.insert(highest_index+1);
                	// consider it as a new bacteria, put an updated state for it at time t+1
                    state(0,0)=positions[t+1][i].x;
                    state(1,0)=positions[t+1][i].y;
                    state(2,0)=std::cos(angles[t+1][i])*v_mean;
                    state(3,0)=std::sin(angles[t+1][i])*v_mean;
                    // correct the orientation (because +/- pi) by comparing with the vector field
                    // use code from opticalflow.cpp

                    // do correction
            	    const cv::Point2f& fxy=flow.at<cv::Point2f>(state(0,0), state(1,0));
            	    // angle of this vector
                    angle_flow=atan2(fxy.x,fxy.y);
                    angle_dif=angles[t+1][i]-angle_flow;
                    // when absolute difference is more than 90°, correct the vx and vy in the dictionary for updated states at time t=1
                    if ((angle_dif> PI/2 && angle_dif< 3*PI/2) || angle_dif<-PI/2){
                    	state(2,0)*= -1;
                	    state(3,0)*= -1;
                    }

            	    state_updated_dict[t+1][highest_index+1]=state;
            	    // TODO: correct angle using vector field ? or not ? (if do it change the energy function to have different computation of ange costs)
        		    p_dict[t+1][highest_index+1]=zeros;
        		    p_dict[t+1][highest_index+1](0,0)=initialNoise_x;
        		    p_dict[t+1][highest_index+1](1,1)=initialNoise_y;
        		    p_dict[t+1][highest_index+1](2,2)=initialNoise_vx;
        		    p_dict[t+1][highest_index+1](3,3)=initialNoise_vy;
            	    highest_index++;
	    	    //}

            }
        }




        ////////////////// 5) correction of predicted state and of predicted error covariance matrices////////////////


        std::cout<<"Number of match: "<<matches.size()<<std::endl;

        for(auto m: matches){

            auto indice_pred=std::get<0>(m.second);
            auto indice_obs=std::get<1>(m.second);



            // observed state at indice_pred+1
            z(0,0)=positions[t+1][indice_obs].x;
            z(1,0)=positions[t+1][indice_obs].y;
            // estimate "observed" velocity
            // to compute pseudo observed velocity, use updated or observed last position ?
            // here use updated last position
            z(2,0)=(positions[t+1][indice_obs].x-state_updated_dict[t][indice_pred](0,0))/delta_t;   // vx
            z(3,0)=(positions[t+1][indice_obs].y-state_updated_dict[t][indice_pred](1,0))/delta_t;   // vy

            // predicted state for time t+1
            x=state_pred_dict[t+1][indice_pred];

            // innovation residual
            y=z-prod(h,x);

            // innovation covariance
            p_pred=p_pred_dict[t+1][indice_pred];
            s=prod(h,p_pred);
            s=prod(s,boost::numeric::ublas::trans(h))+r;

            // Kalman gain
            k=prod(p_pred, boost::numeric::ublas::trans(h));
            InvertMatrix(s, s_inv);
            k=prod(k,s_inv);

            // updated state estimate
            updated_x=x+prod(k,y);

            // updated estimate covariance
            updated_p=identity-prod(k,h);
            updated_p=prod(updated_p, p_pred);

            // put in dictionary
            state_updated_dict[t+1][indice_pred]=updated_x;
            p_dict[t+1][indice_pred]=updated_p;

        }


        // for the predictions that don't have any matching observation, continue predicting (only if number of steps below threshold)
        // when don't find any matching observation, remove the trajectory
        if(multiple_pred){
            for(auto i: state_pred_dict[t+1]){
    	        if(matched_predictions.find(i.first)==matched_predictions.end()){
    	    	    // only if not on an edge region
    	    	    if(i.second(0,0)>edge_region*dx && i.second(0,0)<dx-edge_region*dx && i.second(1,0)>edge_region*dy && i.second(1,0)<dy-edge_region*dy){
    	    		    if(number_pred.find(i.first)==number_pred.end()){
        	    		    number_pred[i.first]=1;
    	    		    }
    	    		    if(number_pred[i.first]==multiple_pred_threshold+1){
    	    		    	//remove trajectory
    	    		    	for (unsigned int j=t+1-multiple_pred_threshold; j<=t+1; j++){
    	    		    	    state_updated_dict[j].erase(i.first);
    	    		    	    // erase also covariance
    	    		    	    p_dict[j].erase(i.first);
    	    		    	}

    	    		    }
    	    		    if(number_pred[i.first]<=multiple_pred_threshold){
            	    	    // put directly predicted state and predicted error covariance without updating in the dictionaries for updated stuff
            		        state_updated_dict[t+1][i.first]=i.second;
    	        	        p_dict[t+1][i.first]=p_pred_dict[t+1][i.first];
    	        	        number_pred[i.first]++;

    	    		    }
    	    	    }
        	    }
            }
        }

        // output
        if(t>multiple_pred_threshold){
        	for(auto i:state_updated_dict[t-multiple_pred_threshold]){
        		out1<< t-multiple_pred_threshold<< " "<<i.first<< " "<< i.second(0,0)<< " "<< i.second(1,0)<< " "<< i.second(2,0)<< " "<< i.second(3,0)<< std::endl;
        		out2<< t-multiple_pred_threshold<< " "<<i.first<< " "<< p_dict[t-multiple_pred_threshold][i.first](0,0)<< " "<< p_dict[t-multiple_pred_threshold][i.first](1,1)<< " "<< p_dict[t-multiple_pred_threshold][i.first](2,2)<< " "<< p_dict[t-multiple_pred_threshold][i.first](3,3)<<std::endl;
        	}
        	state_updated_dict.erase(t-multiple_pred_threshold);
        	p_dict.erase(t-multiple_pred_threshold);
        }

    }

    ////////////////// 6) output ////////////////

    out1.close(); // output for updated predicted states
    out2.close(); // output for error covariance matrices (just output the diagonal)

    std::cout << "Evaluation complete." << std::endl;

}


// TODO
// class for matrices/ not everything in the main
// not store everything for every timestep
// no need for the time vector
// store reference to map elements and dont access the map 1000000000000000000 times



