#ifndef CORRELATION_HEADER
#define CORRELATION_HEADER

#include <vector>
#include <cstdlib>

#include "core.hpp"

#include "Vector2D.hpp"

typedef double real; //float, double, ...
//const real PI = 3.141592653589793238463;

typedef Vector2DT<real> Vector;

real centralDifference(real a, real b, real delta);
real localCurl(const cv::Mat &flow, unsigned int x, unsigned int y, unsigned int flowstep);

real correlationCisnerosXV2M(const cv::Mat &flow, unsigned int flowstep);
real correlationCisnerosXVM2(const cv::Mat &flow, unsigned int flowstep);
real correlationCisnerosXVRM(const cv::Mat &flow, real r, real rstep, unsigned int flowstep);
std::vector<real> correlationCisnerosXVel(const cv::Mat &flow, real flowstep, real rstep);
std::vector<real> correlationPoneXVel(const cv::Mat &flow, real flowstep, real rstep);

real correlationCisnerosXO2M(const cv::Mat &flow, unsigned int flowstep);
real correlationCisnerosXOM2(const cv::Mat &flow, unsigned int flowstep);
real correlationCisnerosOXRM(const cv::Mat &flow, real r, real rstep, unsigned int flowstep);
std::vector<real> correlationCisnerosXVor(const cv::Mat &flow, real flowstep, real rstep);
std::vector<real> correlationPoneXVor(const cv::Mat &flow, real flowstep, real rstep);

real correlationCisnerosTV2M(const std::vector<cv::Mat> &flows, unsigned int x, unsigned int y);
real correlationCisnerosTVM2(const std::vector<cv::Mat> &flows, unsigned int x, unsigned int y);
real correlationCisnerosVTSM(const std::vector<cv::Mat> &flows, unsigned int dt, unsigned int x, unsigned int y);
std::vector<real> correlationCisnerosTVel(const std::vector<cv::Mat> &flows, unsigned int x, unsigned int y);
std::vector<real> correlationPoneTVel(const std::vector<cv::Mat> &flows, unsigned int x, unsigned int y);

real correlationCisnerosTO2M(const std::vector<cv::Mat> &flows, unsigned int x, unsigned int y, unsigned int flowstep);
real correlationCisnerosTOM2(const std::vector<cv::Mat> &flows, unsigned int x, unsigned int y, unsigned int flowstep);
real correlationCisnerosOTSM(const std::vector<cv::Mat> &flows, unsigned int dt, unsigned int x, unsigned int y, unsigned int flowstep);
std::vector<real> correlationCisnerosTVor(const std::vector<cv::Mat> &flows, unsigned int x, unsigned int y, unsigned int flowstep);
std::vector<real> correlationPoneTVor(const std::vector<cv::Mat> &flows, unsigned int x, unsigned int y, unsigned int flowstep);

real correlationCisnerosXV2Mvic(const std::vector<Vector> &velocity);
real correlationCisnerosXVM2vic(const std::vector<Vector> &velocity);
real correlationCisnerosXVRMvic(const std::vector<Vector> &position, const std::vector<Vector> &velocity, real r, real rstep);
std::vector<real> correlationCisnerosXVelvic(const std::vector<Vector> &position, const std::vector<Vector> &velocity, int size, real rstep);
std::vector<real> correlationPoneXVelvic(const std::vector<Vector> &position, const std::vector<Vector> &velocity, int size, real rstep);

#endif //CORRELATION_HEADER
