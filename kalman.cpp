// 2015-02-19

/* TODO
-optional periodic boundary conditions: TEST!
-generalize to more than two frames? test python script to 1:1 find trajectories
*/

/* HOWTO USE
needs the data_exp.txt and data_len.txt from the ImageProcessing.cbp/exe
where there's a 1:1 conversion between the indices in those two files

returns the tracked particles in a Entropy.cbp/exe compatible
format data_exp_tracked.txt: (x,y) position, (a) angle
"""
header
x01 y01 ... x0n x0n
a01 ... a0n
x11 y11 ... x1n x1n
a11 ... a1n
x21 y21 ... x2n x2n
a21 ... a2n
...
""" -> 4 * times lines
where the "0" are from time=1, the "1" are from time=2 and there's a 1:1
conversion between those two, while the "2" are from time=2 and "3" from time=3 and so on
notice that time=2 is double, because in the {0,1} set we track the 0 in 1
and in the {2,3} set we track the 3 in 2;
returns also data_len_tracked.txt: (h,w) height and width
"""
h01 w01 ... h0n w0n
h21 w21 ... h2n w2n
""" -> 1 * times lines
notice that since there is a 1:1 conversion between each two lines in the _exp_tracked
here only every second two line set has a line in _exp_tracked
returns also data_vel_tracked.txt: (vx, vy) velocity in x and y direction
"""
vx01 vy01 ... vx0n vy0n
vx21 vy21 ... vx2n vy2n
...
""" -> 1 * times lines
see the _len_tracked.txt notes, only that here now the difference between {0,1}
from time 1 and 2 is a "0" and the difference between {2,3} from time 2 and 3
is a "2"
*/

#include <iostream>
#include <fstream>
#include <map>
#include <tuple>
#include <set>
#include <stdlib.h>

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

//using namespace boost::numeric::ublas;

typedef std::vector<cv::RotatedRect> RectangleList;
typedef std::map<unsigned int, std::vector<Vector>> Sizedict;
typedef unsigned int uindex;

const cv::Scalar WHITE = cv::Scalar(255);
const cv::Scalar GRAY = cv::Scalar(127);
const cv::Scalar BLACK = cv::Scalar(0);
const cv::Scalar RED = cv::Scalar(0, 0, 255);
const cv::Scalar BLUE = cv::Scalar(255, 0, 0);
const cv::Scalar GREEN = cv::Scalar(0, 255, 0);
const cv::Scalar YELLOW = cv::Scalar(0, 255, 255);
const cv::Scalar ORANGE = cv::Scalar(0, 128, 255);

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

//obsolete
real AngleDifference(real a, real b) //might be wrong if angles are are not from the same interval (eg. [-pi,pi] and [0, 2pi])
{
    //copied from clustering.cpp
    real diff = std::abs(a - b);
    return std::min(diff, 2*PI - diff);
}

/*
ang = Flatten[Table[{i, j}, {i, -4*Pi, 4*Pi, Pi/8}, {j, -4*Pi, 4*Pi, Pi/8}], 1];
angdif[a_, b_] := Min[{Mod[a - b, Pi], Mod[b - a, Pi]}]
Map[angdif[#[[1]], #[[2]]] &, ang];
ListPlot[%, PlotRange -> All]
*/
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

//reference, pointer?????
//imagesize in px
real rate_match(Vector pos1, Vector pos2, real ang1, real ang2, Vector size1, Vector size2, real pospara_factor, real posorth_factor, real ang_factor, real area_factor, real L, bool periodic=false)
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
    real pospara_rating = pospara_factor * std::get<0>(com_distance) / L;
    real posorth_rating = posorth_factor * std::get<1>(com_distance) / L;

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

    //return +posrating + angrating + arearating;
    return pospara_rating + posorth_rating + ang_rating + area_rating;
}

real calc(unsigned int t, unsigned int t2, Posdict &positions, Angdict &angles, Sizedict &sizes, uindex particle, uindex other_particle, real L, bool periodic, real pospara_factor, real posorth_factor, real angfactor, real areafactor)
{
    return rate_match(positions[t][particle], positions[t2][other_particle], angles[t][particle], angles[t2][other_particle], sizes[t][particle], sizes[t2][other_particle], pospara_factor, posorth_factor, angfactor, areafactor, L, periodic);
}


real calcKal(unsigned int t, Posdict &positions, Matrixdict &pred, Angdict &angles, Sizedict &sizes, auto particle, uindex other_particle, real L, bool periodic, real pospara_factor, real posorth_factor, real angfactor, real areafactor)
{
    return rate_match(Vector(particle.second(0,0),particle.second(1,0)), positions[t][other_particle], atan2(particle.second(2,0),particle.second(3,0)), angles[t][other_particle], sizes[t-1][particle.first], sizes[t][other_particle], pospara_factor, posorth_factor, angfactor, areafactor, L, periodic);
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


real FixAngle(Posdict &positions, Angdict &angles, std::tuple<uindex, uindex, uindex> &msecond, unsigned int t, unsigned int t2, unsigned int whichparticle, bool last_step = false)
{
    unsigned int vt1, vt2;
    if (last_step)
    {
        vt1 = t-1; //is this correct?
        vt2 = t2-1;
    }
    else
    {
        vt1 = t;
        vt2 = t2;
    }
    //at1 = t
    //at2 = t2

    // does atan2 give us "fixed" angles from fixed coords?
    real velocityangle = std::atan2(positions[vt2][std::get<1>(msecond)].y - positions[vt1][std::get<0>(msecond)].y, positions[vt2][std::get<1>(msecond)].x - positions[vt1][std::get<0>(msecond)].x);
    //not a nice hack : ( for std::get<"MUST HAVE CONSTANT">()
    real angleangle;
    if (whichparticle == 0)
    {
        angleangle = angles[t][std::get<0>(msecond)];
    }
    else
    {
        angleangle = angles[t][std::get<1>(msecond)];
    }

    real ang_distance = AngleDifference180(velocityangle, angleangle); //not very very well tested....

    if (ang_distance > PI/2)
    {
        if (angleangle > 0)
        {
            angleangle -= PI;
        }
        else
        {
            angleangle += PI;
        }
    }
    return angleangle;
}
// TODO: declare it somewhere ??
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


int main(void)
{

	//////////////// 1) Initialization ////////////////

	// TODO: make that the arguments can be entered from command line or from configuration file

	// defines the part of the images that we are using (has to be the same as in imageprocessing.cpp)
    int X = 512;
    int Y = 0;
    int DX = 512;
    int DY = 512;
    // TODO: somehow write this information in the output of imageprocessing.cpp so that we don't need to enter it again

    // times is a vector containing all the indices of the frames we want to use (time_stop not included)
    unsigned int time_start = 1, time_stop = 10, time_step = 1;
    std::vector<unsigned int> times = MakeTimes(time_start, time_stop, time_step);

    // input files:
    // in_file contains x and y positions for each box, and angle of orientation (in radians, +/- pi) and that for each frame
    // in_file2 contains length and width of each box, and that for each frame
    std::string dir = "";
    std::string in_file = dir + "data_exp_noduplicates.txt"; //data_exp.txt
    std::string in_file2 = dir + "data_len_noduplicates.txt"; //data_len.txt
    // TODO: check if maybe we could not use directly data_exp.txt, ... because we do not seem to have the duplication bug here

    // put observed data in dictionaries: one for the positions (x and y), one for the angles and one for the sizes (length and width) of the fitted rectangles
    Posdict positions;
    Angdict angles;
    Sizedict sizes;
    read_data_single(in_file, positions, angles, times);
    std::string header = read_data_len(in_file2, sizes, times);

    // dictionary for predictions of states (dictionary of dictionary: time step, identity of bacteria and state)
    Matrixdict state_pred_dict;
    // vector of one predicted state (x, y, vx, vy)
    boost::numeric::ublas::matrix<real> state_pred(4,1);

    // dictionary for updated states (time, identity of bacteria and state)
    Matrixdict state_updated_dict;

    // TODO
    // somehow find some velocity values between the 2 first frames
    // for now will just use the mean velocity

    // TODO: find the mean velocity with the good unit
    //mean velocity: 34.5 micro meters per second
    real v_mean=0.1;

    // time difference between 2 consecutive frames
    // what is it? say 1 (time in unit of delta t)
    int delta_t=1;

    // because in this program we need a lot of matrix operations use a library: boost
    // (www.boost.org/doc/libs/1_61_0/libs/numeric/ublas/doc/)

    // 4 by 4 matrix with only zeros
	boost::numeric::ublas::matrix<real> zeros = boost::numeric::ublas::zero_matrix<real>(4,4);

	// identity 4
    boost::numeric::ublas::identity_matrix<real> identity(4);


	// initialize error covariance matrices P
    real initialNoise_x=0.01;
    real initialNoise_y=0.01;
    real initialNoise_vx=0.01;
    real initialNoise_vy=0.01;

	// TODO: change and not put only zero: put something in the diagonal
	Matrixdict p_dict; // dictionary of dictionaries, with times step and identities
    unsigned int sizeInit=positions[1].size();
	for (unsigned int i=0; i<sizeInit; i++){
		p_dict[1][i]=zeros;
		p_dict[1][i](0,0)=initialNoise_x;
		p_dict[1][i](1,0)=initialNoise_y;
		p_dict[1][i](2,0)=initialNoise_vx;
		p_dict[1][i](3,0)=initialNoise_vy;

	}

	std::cout<<p_dict[1][0]<<std::endl;

	// predicted error covariance matrix dictionary
	Matrixdict p_pred_dict;
	// one error covariance matrix
    boost::numeric::ublas::matrix<real> p_pred(4,4);

	// matrix f for state transition model
	boost::numeric::ublas::matrix<real> f=identity;
    f(0,2)=delta_t;
    f(1,3)=delta_t;

    // matrix h for observation model
    boost::numeric::ublas::identity_matrix<real> h(4);

    // matrix q for covariance of process noise
    // identity with a weight
     // TODO: what weight to put for q and r?
    int weight_q=10;
    boost::numeric::ublas::matrix<real> q=weight_q*identity;

    // matrix r for covariance of observation noise
    int weight_r=10;
    boost::numeric::ublas::matrix<real> r=weight_r*identity;

    // counter for highest used bacteria index;
    int highest_index=-1;

    // initialize updated states at time t=1
    boost::numeric::ublas::matrix<real> state(4,1);
    for(unsigned int i=0; i<sizeInit; i++){
    	state(0,0)=positions[1][i].x;
    	state(1,0)=positions[1][i].y;
    	state(2,0)=std::cos(angles[1][i])*v_mean;
    	state(3,0)=std::sin(angles[1][i])*v_mean;
    	//state(2,0)=0;    	//or initialize the velocity to zero and hope that will be corrected by the filter
    	//state(3,0)=0;
    	state_updated_dict[1][i]=state;
    	highest_index ++;   // increment the counter
    }

    std::cout<<highest_index<<std::endl;
    // using the two first frames, correct the orientation (because +/- pi) by comparing with the vector field
    // use code from opticalflow.cpp

    // two first pictures
    // TODO: make that an argument or sth
    cv::Mat orig_1 = cv::imread("img_000000000_00-BF_EGFP_000.tif", CV_LOAD_IMAGE_GRAYSCALE);
    cv::Mat orig_2 = cv::imread("img_000000001_00-BF_EGFP_000.tif", CV_LOAD_IMAGE_GRAYSCALE);

    // use only the subregion we are working on
    cv::Rect subrect = cv::Rect(X, Y, DX, DY);
	cv::Mat subImage_orig_1(orig_1, subrect);
	orig_1=subImage_orig_1;
	cv::Mat subImage_orig_2(orig_2, subrect);
	orig_2=subImage_orig_1;

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
    for (auto i:state_updated_dict[1]){
    	// flow vector at the position of the bacteria
    	// why is the y axis upside down ?
    	// TODO: check where is the problem, maybe don't need to fix here
    	const cv::Point2f& fxy=flow.at<cv::Point2f>(i.second(0,0), DY-i.second(1,0));
    	// angle of this vector
        angle_flow=atan2(fxy.x,fxy.y);
        angle_dif=angles[1][i.first]-angle_flow;
        // when absolute difference is more than 90°, correct the vx and vy in the dictionary for updated states at time t=1
        if ((angle_dif> PI/2 && angle_dif< 3*PI/2) || angle_dif<-PI/2){
            i.second(2,0)*= -1;
            i.second(3,0)*= -1;
        }
    }

    // edge of whole considered area ?
    // TODO: look at what happens when not square
    real L=DX;

    // weight coefficients for energy function (used to do the assignment between predictions and observations)

    // for parallel distance: 1/average length  (10.)
    int count=0;
    real pospara_factor = 0.;
    for (unsigned int t: times){
    	for (Vector bact: sizes[t]){
    	    pospara_factor+=bact.x;
    	    count++;;
    	}
    }
    pospara_factor/=count;
    pospara_factor=(1/pospara_factor)*L;

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
    posorth_factor=(1/posorth_factor)*L;

    real angfactor = 1.;		// for angular displacement
    real areafactor = 1.;		// for area change
    real rating_threshold = 2.; // (1.)  rating above means that cannot be the same particle


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


    // TODO
    // open output files (also put some headers ?)
    // updated states
    std::string out_file_s = "updated_states.txt";
    std::ofstream out1;
    out1.open(out_file_s);

    // updated error covariance matrices p
    std::string out_file_p = "updated_error_covariance.txt";
    std::ofstream out2;
    out2.open(out_file_p);

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


        //////////////// 3) assignment between predictions and observations, at time t+1 ////////////////
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
                    real rate = calcKal(t+1, positions, state_pred_dict, angles, sizes, particle, other_particle, L, periodic_boundary_cond, pospara_factor, posorth_factor, angfactor, areafactor);
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

        // TODO:

        // what to do with observation that were not close enough to any prediction ? Say it's a new bacteria ??
        for (unsigned int i=0; i<nextnum; i++){
        	// look if this observation has been assigned to a prediction
            if (matched_observations.find(i)==matched_observations.end()){
            	// consider it as a new bacteria, put an updated state for it at time t+1
                state(0,0)=positions[t+1][i].x;
                state(1,0)=positions[t+1][i].y;
                state(2,0)=std::cos(angles[t+1][i])*v_mean;
                state(3,0)=std::sin(angles[t+1][i])*v_mean;
            	state_updated_dict[t+1][highest_index+1]=state;
            	// TODO: correct angle using vector field
        		p_dict[t+1][highest_index+1]=zeros;
        		p_dict[t+1][highest_index+1](0,0)=initialNoise_x;
        		p_dict[t+1][highest_index+1](1,0)=initialNoise_y;
        		p_dict[t+1][highest_index+1](2,0)=initialNoise_vx;
        		p_dict[t+1][highest_index+1](3,0)=initialNoise_vy;
            	highest_index++;

            }
        }





        //////////////// 4) prediction of error covariance matrices ////////////////
        for(auto i:p_dict[t]){
            p_pred=prod(f, i.second);
            p_pred_dict[t+1][i.first]=prod(p_pred,boost::numeric::ublas::trans(f))+q;
        }

        ////////////////// 5) correction of predicted state and of predicted error covariance matrices////////////////

        // (idea: could we use the values of the energy function to help estimating confidence ?)

        //output time in ouput files
        out1<<t+1<<" ";
        out2<<t+1<<" ";

        for(auto m: matches){
            auto indice_pred=std::get<0>(m.second);
            auto indice_obs=std::get<1>(m.second);

            // observed state at tindice_pred+1
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

            //output
            out1<<indice_pred<<" ";
            out1<<updated_x(0,0)<<" "<<updated_x(1,0)<<" "<<updated_x(2,0)<<" "<<updated_x(3,0)<<" ";

            out2<<updated_p(0,0)<<" "<<updated_p(1,1)<<" "<<updated_p(2,2)<<" "<<updated_p(3,3)<<" ";
        }
        out1<<std::endl;
        out2<<std::endl;

    }

    // TODO:
    ////////////////// 6) output ////////////////

    out1.close(); // output for updated predicted states
    out2.close(); // output for error covariance matrices

    ////////////////// 7) visualization of results ////////////////


    // for each time step, plot the updated position of the bacteria on the original picture
    // keep the same color for the same bacteria over multiple frames

    cv::Mat original;
    cv::Point bac;    // contains coordinate x and y of the point
    cv::Scalar col;   //color of the point (blue, yellow, red)

    // TODO: take that as argument or sth
    std::string FILENAMETPL = "img_%09d_00-BF_EGFP_000.tif";
	int START = 0;
	int STOP = 10;
    std::vector<std::string> filenames = CreateFilenameList(FILENAMETPL.c_str(), START, STOP);

    std::string output_name; // name of an output file

    // assign a different random blue, yellow and red amount to each bacteria identity, keep it in vectors
	srand(1); //seed for random choice of the color
	std::vector<int> bleu, jaune, rouge;
	for (int i=0; i<highest_index; i++){
		bleu.push_back(rand()%200);
	}
	for (int i=0; i<highest_index; i++){
		jaune.push_back(rand()%200);
	}
	for (int i=0; i<highest_index; i++){
		rouge.push_back(rand()%200);
	}

    real x_cov, y_cov;
    cv::Size axes;
    double angle=0, startAngle=0, endAngle=359;
	for(unsigned int t:times){
		// import the original picture
		original= cv::imread(filenames[t-1]);
		// keep only the part we worked on
		cv::Mat subImage(original, subrect);
		original=subImage;
		for(auto i: state_updated_dict[t]){
			// color
	    	col[0]=bleu[i.first];
	    	col[1]=jaune[i.first];
	    	col[2]=rouge[i.first];
	    	// coordinate x and y of the point
	    	bac=cv::Point(i.second(0,0), DY-i.second(1,0));
	    	// make the point
	        cv::circle(original, bac, 2, col, -1);
	        // plot covariance for x and y -> ellipse
	        x_cov=p_dict[t][i.first](0,0);
	        y_cov=p_dict[t][i.first](1,1);
            axes=cv::Size(x_cov, y_cov);
	        cv::ellipse(original, bac, axes, angle, startAngle, endAngle, col);

	    }
	    output_name="pred_"+filenames[t-1];
	    cv::imwrite(output_name, original);
	}


	// try following one bacteria and plot the error covariance each time



    std::cout << "Evaluation complete." << std::endl;

}


// TODO
// class for matrices/ not everything in the main
// not store everything for every timestep
// no need for the time vector
// put arguments in configuration file
// store reference to map elements and dont access the map 1000000000000000000 times



