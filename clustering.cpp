// 2015-12-15

/* TODO
-optional periodic boundary conditions: TEST!
*/

/*
HOW TO USE
-needs 1:1 particle conversion, meaning same number of positions and angles per line
and particle 0 in line 1 is particle 0 in line 2 and so on
-so it should work on Vicsek.exe and SPR.exe (but not on Tracking.exe yet)
    -> do a python script that does that 1:1 conversion
*/

#include <iostream>
#include <fstream>
#include <map>
#include <tuple>
#include <algorithm>
#include <set>
#include <numeric>

#include "core.hpp"
#include "highgui.hpp"
#include "imgproc.hpp"
#include "imgcodecs.hpp"
#include "imgcodecs/imgcodecs_c.h"

#include "Vector2D.hpp"
#include "constantsImgPro.hpp"

typedef double real; //float, double, ... //not needed
typedef Vector2DT<real> Vector; //not needed

typedef std::vector<cv::RotatedRect> RectangleList;

const cv::Scalar WHITE = cv::Scalar(255);
const cv::Scalar GRAY = cv::Scalar(127);
const cv::Scalar BLACK = cv::Scalar(0);
const cv::Scalar RED = cv::Scalar(0, 0, 255);
const cv::Scalar BLUE = cv::Scalar(255, 0, 0);
const cv::Scalar GREEN = cv::Scalar(0, 255, 0);
const cv::Scalar YELLOW = cv::Scalar(0, 255, 255);
const cv::Scalar ORANGE = cv::Scalar(0, 128, 255);

//COPIED from tracking.cpp
real AngleDifference180(real a, real b)
{
    return std::min(mod(a - b, 2*PI), mod(b - a, 2*PI));
}

bool cluster_condition(Vector pos1, Vector pos2, real ang1, real ang2, real max_com_distance, real max_ang_distance, real L, bool periodic)
{
    //center of mass distance
    real com_distance;
    //PERIODIC!????????
    //com_distance = distance_periodic(positions[i], positions[j], L);
    //NON-PERIODIC??????? -> more yes than no (periodic for simulation, nonperiodic for exp?)
    if (!periodic)
    {   
        //real xdiff = pos1.x-pos2.x;
        //real ydiff = pos1.y-pos2.y;
        //com_distance = std::sqrt(xdiff*xdiff + ydiff*ydiff);
        com_distance = pos1.GetDistance(pos2);
    }
    else
    {
        com_distance = distance_periodic(pos1, pos2, L);
    }

    //angle difference CORRECT????
    //real ang_distance = std::abs(ang1 - ang2);
    //ang_distance = std::min(ang_distance, 2*PI - ang_distance);
    real ang_distance = AngleDifference180(ang1, ang2); //not very very well tested....

    return (com_distance <= max_com_distance) && (ang_distance <= max_ang_distance);
}

struct Counter
{
    size_t count = 0;
    struct value_type { template<typename T> value_type(const T&) { } };
    void push_back(const value_type&) { ++count; }
};

template <typename T> size_t set_intersect_len(const std::set<T> &a, const std::set<T> &b)
{
    Counter c;
    std::set_intersection(a.begin(), a.end(), b.begin(), b.end(), std::back_inserter(c));
    return c.count;
}

void test_set_intersect_len(void)
{
    std::set<int> a = {1, 2, 3};
    std::set<int> b = {2, 3, 4};
    std::set<int> c = {4, 5, 6};

    assert(set_intersect_len(a, a) == 3);
    assert(set_intersect_len(a, b) == 2);
    assert(set_intersect_len(b, c) == 1);
    assert(set_intersect_len(a, c) == 0);
}

int main(void)
{
    //pass these over cmd args
    std::string in_file = "data.txt";
    std::string out_file = "data_clustered.txt";
    unsigned int time_start = 1, time_stop = 5001, time_step = 1;
    real max_com_distance = 5.*1.2; //maximum center of mass distance
    real max_ang_distance = PI/12.; //maximum angle difference
    real L = 1024.; //square for periodic boundary conditions
    bool periodic_boundary_cond = false;

    std::ofstream out;
    out.open(out_file);

    Posdict positions;
    Angdict angles;
    std::vector<unsigned int> times = MakeTimes(time_start, time_stop, time_step);
    //read in file
    std::string header = read_data_single(in_file, positions, angles, times);

    std::cout << "Initialisation complete." << std::endl;

    out << header << std::endl;

    /*pseudo code for 1 time cluster*/
    //for each time
        //make boxes
        //make cluster vector
        //for each particle
            //if not in cluster: new cluster (id?)
            //for each particle in 9 boxes:
                //if not in cluster and cluster conditions match: 2nd particle in cluster (with 1st's id)
                //if in cluster and cluster conditions match: join clusters (keep 1st's id)
                //else: nothing
        //output info
            //size, avg heading, std of heading of each cluster; which other infos?????

    std::map<unsigned int, std::vector<int>> clusters; //maps times to vectors with index posid and value clusterid
    int delete_index, cluster_index = 0;
    bool cluster_cond;
    for (unsigned int t : times)
    {
        std::cout << "single frame (part 1/2), t: " << t << "/" << times.size() << std::endl;
        //make boxes and distribute particles into them

        //initialize with -1s
        std::vector<int> cluster(positions[t].size(), -1);

        //later: iterate over all particles in 9 boxes
        for (unsigned int p = 0; p < positions[t].size(); ++p)
        {
            //particle not in cluster yet
            if (cluster[p] == -1)
            {
                cluster[p] = cluster_index;
                cluster_index++;
            }

            for (unsigned int p2 = 0; p2 < positions[t].size(); ++p2)
            {
                if (p == p2)
                {
                    continue;
                }
                cluster_cond = cluster_condition(positions[t][p], positions[t][p2], angles[t][p], angles[t][p2], max_com_distance, max_ang_distance, L, periodic_boundary_cond);
                //p and p2 are cluster connected
                if (cluster_cond)
                {

                    //p2 already in cluster
                    if (cluster[p2] != -1)
                    {
                        //remember id of p2
                        delete_index = cluster[p2];
                        //put p2 in cluster of p
                        cluster[p2] = cluster[p];

                        //join clusters (now the id of p2 is free?????!!!!!!)
                        //naiv: faster method??? like for example bookkeeping of all clusters
                        for (unsigned int c = 0; c < positions[t].size(); ++c)
                        {
                            if (cluster[c] == delete_index)
                            {
                                cluster[c] = cluster[p];
                            }
                        }
                    }
                    //p2 not in cluster yet
                    else
                    {
                        cluster[p2] = cluster[p];
                    }
                }
            }
        }
        clusters[t] = cluster;
    }

    /*pseudo code for multiple time cluster identification*/
    //track consecutive time clusters here
    //for each cluster1 in time1
        //for each cluster2 in time2
            //count intersection both clusters -> 1-2 for loops
            //if count(intersection) > 0.5 * count(groesserer cluster)
                //cluster1 === cluster2, ie cluster in time2 bekommt id von time1

    //identify clusters in consecutive times:
    //cluster from t2 = clsuter from t1 if either contains
    //more than 50% of the others particles in an intersection
    for (auto t : times)
    {
        std::cout << "cluster over time (part 2/2), t: " << t << "/" << times.size() << std::endl;
        if (t == times.back())
        {
            break;
        }
        auto t2 = t+1;
        //otherwise out of bounds

        //build maps: mapping cluster id to set of particle ids (their
        //index in positions and angles)
        std::map<int, std::set<unsigned int>> clustersets1;
        std::map<int, std::set<unsigned int>> clustersets2;
        for (unsigned int i = 0; i < clusters[t].size(); ++i)
        {
            clustersets1[clusters[t][i]].insert(i);
        }
        for (unsigned int i = 0; i < clusters[t2].size(); ++i)
        {
            clustersets2[clusters[t2][i]].insert(i);
        }

        auto mapcolsum = [](unsigned int a, const auto &b) -> unsigned int {return a+b.second.size();};
        assert(std::accumulate(clustersets1.begin(), clustersets1.end(), 0, mapcolsum) == positions[t].size());
        assert(std::accumulate(clustersets2.begin(), clustersets2.end(), 0, mapcolsum) == positions[t2].size());

        for (auto item1 : clustersets1)
        {
            for (auto item2 : clustersets2)
            {
                int len = set_intersect_len(item1.second, item2.second);
                unsigned int bigger_cluster_len = std::max(item1.second.size(), item2.second.size());

                //ids from t are copied to t2
                if (len > 0.5 * bigger_cluster_len)
                {
                    for (auto &clusterid : clusters[t2])
                    {
                        if (clusterid == item2.first)
                        {
                            clusterid = item1.first;
                        }
                    }
                }
            }
        }
    }

    //do all the output here
    for (auto t : times)
    {
        for (unsigned int cl = 0; cl < clusters[t].size() - 1; ++cl)
        {
            out << clusters[t][cl] << " ";
        }
        out << clusters[t].back() << std::endl;
    }

    out.close();
    std::cout << "Evaluation complete." << std::endl;
}
