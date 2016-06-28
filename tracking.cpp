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

#include "core.hpp"
#include "highgui.hpp"
#include "imgproc.hpp"
#include "imgcodecs.hpp"
#include "imgcodecs/imgcodecs_c.h"

#include "Vector2D.hpp"
#include "constantsImgPro.hpp"

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

int main(void)
{
    //ParaOrthDistance_test();

    std::string dir = "";

    std::string in_file = dir + "data_exp_noduplicates.txt"; //data_exp.txt
    std::string in_file2 = dir + "data_len_noduplicates.txt"; //data_len.txt
    std::string out_file = dir + "data_exp_tracked.txt";
    std::string out_file2 = dir + "data_len_tracked.txt";
    std::string out_file3 = dir + "data_vel_tracked.txt";
    std::string out_file4 = dir + "data_extras_tracked.txt";
    unsigned int time_start = 1, time_stop = 11, time_step = 1;

    real L = 512.; //assert quadratic

    //0.25 4 4 0.4 0.2
    real rating_threshold = 1.0; // rating above means: cannot be same particle
    //real posfactor = 4.;
    real pospara_factor = 10.;
    real posorth_factor = 50.;
    real angfactor = 1.;
    real areafactor = 1.;

    bool periodic_boundary_cond = false; //true not yet implemented

    Posdict positions;
    Angdict angles;
    Sizedict sizes;
    std::vector<unsigned int> times = MakeTimes(time_start, time_stop, time_step);


    std::ofstream out;
    out.open(out_file);
    std::ofstream out2;
    out2.open(out_file2);
    std::ofstream out3;
    out3.open(out_file3);
    std::ofstream out4;
    out4.open(out_file4);

    //read in detection file from ImageProcessing.cbp
    read_data_single(in_file, positions, angles, times);
    std::string header = read_data_len(in_file2, sizes, times); //header in_file1 == in_file2??

    std::cout << "Initialisation complete." << std::endl;

    //output header
    out << header << std::endl;
    out2 << header << std::endl;
    out3 << header << std::endl;

    //DO TRACKING

    //CARE: angles ambigious +-180°: have some way to tell AFTER tracking!! (-> displacement)

    for (auto t : times)
    {
        //otherwise out of bounds
        if (t == times.back())
        {
            break;
        }
        auto t2 = t+1;
        unsigned int prevnum = sizes[t].size();
        unsigned int nextnum = sizes[t2].size();

        std::cout << "Number of compared particles: " << prevnum << " " << nextnum << std::endl;

        std::vector<uindex> prevbox = range(prevnum); // all indices
        std::vector<uindex> nextbox = range(nextnum);

        std::vector<std::vector<uindex>> boxes;
        boxes.push_back(prevbox);
        std::vector<std::vector<uindex>> neighbouring_nextboxes; //including self, needs to be a function
        neighbouring_nextboxes.push_back(nextbox);

        /*for (uindex i = 0; i < prevnum; ++i)
        {
            assert(positions[t][i] == positions[t2][i]);
            assert(angles[t][i] == angles[t2][i]);
            if (sizes[t][i] != sizes[t2][i])
            {
                std::cout << sizes[t][i].x << " " << sizes[t][i].y << ", " << sizes[t2][i].x << " " << sizes[t2][i].y << std::endl;
            }
        }*/

        //std::vector< std::map<real,uindex,std::greater<real>> > rates;
        //vector with dictionary (map) for each particle p that maps all ratings
        //of all particles to their index p2 (including self) and is ordered
        //from low rating (good match) to high rating
        //length is varying since all ratings > rating_threshold are cut
        std::vector< std::map<real,uindex> > rates(prevnum);

        /*
        for all combinations of particles of two frames:
            calculate rating based parallel and orthogonal distance, relative angle and relative size
            low rating is better
            combinations with rating higher than threshold will be ignored
        */
        for (auto box : boxes)
        {
            for (uindex particle : box)
            {
                for (auto nbox : neighbouring_nextboxes)
                {
                    for (uindex other_particle : nbox)
                    {
                        //map should contain only tracking partners with rating below threshold
                        real rate = calc(t, t2, positions, angles, sizes, particle, other_particle, L, periodic_boundary_cond, pospara_factor, posorth_factor, angfactor, areafactor);
                        if (rate <= rating_threshold)
                        {
                            rates[particle][rate] = other_particle;
                        }
                    }
                }
            }
        }

        //console output for debugging
        /*
        for (auto particle : rates)
        {
            //for (auto p : particle)
            auto p = particle.begin();

            for (uindex j = 0; j < 10 && j < particle.size(); ++j)
            {
                std::cout << p->first << " " << p->second << " ";
                ++p;
            }
            std::cout << std::endl;
        }
        std::cout << std::endl;
        */

        //populate structure that maps rating on tuple (p, p2, ranking)
        //of best matching tracking partner p->p2 for each p
        //where ranking is the index in the rating map rates[p]
        //e.g. 0 -> 0, 0; 0 -> 5, 0.1; 0 -> 8, 0.2 it would mean {0.2 : (0, 8, 2)}
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

        //iterative scheme to find tracking partners

        /*

        /*
        make list of mappings with highest ratings for each particle sorted by rating, called prev mappings (particles may have same successor)

        create new mapping:
        for each particle:
            find successor with highest rating
            if successor was already used by another particle, find the one with the next highest rating
            if no soccessor can be found for particle: ignore
        set prev mappings to new mappings
        if another successor had to be found for particle,
        repeat this process to ensure the particle with the highest rating is found
        in case this successor is (again) used by multiple particles

        */

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
            used_matches.clear();
            newmatches.clear(); //hm really? ka, apparently this fixes everything (at least something)
        }

        // calculate sum
        real ratesum = 0;
        for (auto m : matches)
        {
            ratesum += m.first;
        }
        std::cout << "Matches: " << matches.size() << ", Quality: " << ratesum / matches.size() << std::endl;
        out4 << "matches: "<< matches.size() << " ratesum: " << ratesum;

        //output in console
        /*
        for (auto m : matches)
        {
            std::cout << std::get<0>(m.second) << " " << std::get<1>(m.second)
            << " (" << m.first << ")" << std::endl;
        }
        std::cout << std::endl;
        */

        //CARE: the output has " " in end of line -> so reading in this data
        //again with our c++ read in functions will result in an error (but no msg)
        //output in file
        for (auto m: matches)
        {
            //data_exp_tracked.txt
            out << positions[t][std::get<0>(m.second)].x << " " << positions[t][std::get<0>(m.second)].y << " ";
            //data_len_tracked.txt
            out2 << sizes[t][std::get<0>(m.second)].x << " " << sizes[t][std::get<0>(m.second)].y << " ";
            //data_vel_tracked.txt
            out3 << positions[t2][std::get<1>(m.second)].x - positions[t][std::get<0>(m.second)].x << " " << positions[t2][std::get<1>(m.second)].y - positions[t][std::get<0>(m.second)].y << " ";
        }


        out << std::endl;
        out2 << std::endl;
        out3 << std::endl;
        out4 << std::endl;

        // correct the angles, because the detection leaves them arbitrary up to 180°
        // if angle between velocity and vector from angles is > 90° => rotate angle by 180°
        for (auto m: matches)
        {

            //out << angles[t][std::get<0>(m.second)] << " ";
            out << FixAngle(positions, angles, m.second, t, t2, 0) << " ";
        }
        out << std::endl;
        for (auto m: matches)
        {
            out << positions[t2][std::get<1>(m.second)].x << " " << positions[t2][std::get<1>(m.second)].y << " ";
        }
        out << std::endl;

        auto t3 = t2+1;
        //for last time velocity: take the velocity from before (in all other cases it is between t2 and t3
        for (auto m: matches)
        {
            //out << angles[t2][std::get<1>(m.second)] << " ";
            out << FixAngle(positions, angles, m.second, t2, t3, 1, t2 == times.back()) << " ";
        }
        out << std::endl;
    }
    out.close();
    out2.close();
    out3.close();
    out4.close();
    std::cout << "Evaluation complete." << std::endl;
}





