// 2015-10-09

#include <iostream>
#include <fstream>
#include <map>
#include <tuple>
#include <iomanip>

#include "core.hpp"
#include "highgui.hpp"
#include "imgproc.hpp"
#include "imgcodecs.hpp"
#include "imgcodecs/imgcodecs_c.h"
#include "video/tracking.hpp"
#include "videoio.hpp"

#include "correlation.hpp"

const cv::Scalar WHITE = cv::Scalar(255);
const cv::Scalar GRAY = cv::Scalar(127);
const cv::Scalar BLACK = cv::Scalar(0);
const cv::Scalar RED = cv::Scalar(0, 0, 255);
const cv::Scalar BLUE = cv::Scalar(255, 0, 0);
const cv::Scalar GREEN = cv::Scalar(0, 255, 0);
const cv::Scalar YELLOW = cv::Scalar(0, 255, 255);
const cv::Scalar ORANGE = cv::Scalar(0, 128, 255);

static void drawOptFlowMap(const cv::Mat& flow, cv::Mat& cflowmap, int step, double factor, const cv::Scalar& color)
{
    for(int y = 0; y < cflowmap.rows; y += step)
    {
        for(int x = 0; x < cflowmap.cols; x += step)
        {
            const cv::Point2f& fxy = flow.at<cv::Point2f>(y, x);
            cv::line(cflowmap, cv::Point(x,y), cv::Point(cvRound(x+factor*fxy.x), cvRound(y+factor*fxy.y)), color);
            //cv::circle(cflowmap, cv::Point(x,y), 2, color, -1);
        }
    }
}

void saveFlowMap(const std::string &filename, const cv::Mat& flow, int step)
{
    std::ofstream out;
    out.open(filename);

    for (int y = 0; y < flow.rows; y += step)
    {
        for (int x = 0; x < flow.cols; x += step)
        {
            const cv::Point2f& fxy = flow.at<cv::Point2f>(y, x);
            out << fxy.x << " " << -fxy.y << " ";
        }
        out << std::endl;
    }
    out.close();
}

void convertFlowMatrix(const cv::Mat &flow, cv::Mat &newflow, unsigned int reducesize)
{
    unsigned int x, y;
    newflow.create(flow.rows / reducesize, flow.cols / reducesize, flow.type());

    for (x = 0; x < flow.cols; x += reducesize)
    {
        for (y = 0; y < flow.rows; y += reducesize)
        {
            newflow.at<cv::Point2f>(y / reducesize, x / reducesize) = flow.at<cv::Point2f>(y, x);
        }
    }
}

int main(void)
{
    cv::Mat orig_1, orig_2, norm_1, norm_2, flow, smallflow, orig_wflow;

    unsigned int i, winsize = 16, smallwinsize = 1;
    unsigned int no_of_flows = 999; // 10

    std::string indir = "./"; // "rawpics"
    std::string outdir = "./"; // "flow"

    std::ofstream corxvelout, corxvorout, cortvelout, cortvorout;
    corxvelout.open(outdir+"correlations-x-vel.txt");
    corxvorout.open(outdir+"correlations-x-vor.txt");
    cortvelout.open(outdir+"correlations-t-vel.txt");
    cortvorout.open(outdir+"correlations-t-vor.txt");

    std::cout << "Initialization complete." << std::endl;

    std::vector<cv::Mat> flows(no_of_flows);

    for (i = 0; i < no_of_flows; ++i)
    {
        std::cout << "Optical flow: " << i << std::endl;
        std::stringstream ss;
        //load original image
        ss << std::setw(9) << std::setfill('0') << i;
        std::string s = ss.str();
        orig_1 = cv::imread(indir+"img_" + s + "_02-BF1_000.tif", CV_LOAD_IMAGE_GRAYSCALE); //img_000000000_02-BF1_000.tif
        ss.str("");
        ss << std::setw(9) << std::setfill('0') << i+1;
        s = ss.str();
        orig_2 = cv::imread(indir+"img_" + s + "_02-BF1_000.tif", CV_LOAD_IMAGE_GRAYSCALE);

        //histogram equalizer
        cv::equalizeHist(orig_1, norm_1);
        cv::equalizeHist(orig_2, norm_2);
        //cv::imwrite("optFlow_Normalized.png", normalized);

        //bilateral filter
        //cv::bilateralFilter(normalized, blurred, 20, 100, 100); //20, 100, 100
        //cv::imwrite("test_Blurred.png", blurred);


        cv::calcOpticalFlowFarneback(norm_1, norm_2, flow, 0.5, 4, winsize, 3, 5, 1.1, 0);
        convertFlowMatrix(flow, smallflow, winsize);

        flows[i] = smallflow.clone();
        //cv::calcOpticalFlowSF(norm_1, norm_2, flow, 3, 2, 4);
        //cvCalcOpticalFlowHS(const CvArr* prev, const CvArr* curr, int use_previous, CvArr* velx, CvArr* vely, double lambda, CvTermCriteria criteria);
        cv::cvtColor(orig_1, orig_wflow, cv::COLOR_GRAY2BGR);
        //uflow.copyTo(flow);
        drawOptFlowMap(flow, orig_wflow, winsize, 1., cv::Scalar(0, 255, 0));

        ss.str("");
        ss << std::setw(9) << std::setfill('0') << i;
        s = ss.str();
        cv::imwrite(outdir+"flow_farneback" + s + ".png", orig_wflow);
        saveFlowMap(outdir+"flow_" + s + ".txt", flow, winsize);
    }

    for (i = 0; i < flows.size(); ++i)
    {
        std::cout << "spatial correlation (velocity): " << i << std::endl;
        //spatial velocity
        for (auto c : correlationCisnerosXVel(flows[i], smallwinsize, smallwinsize))
        {
            corxvelout << c << " ";
        }
        corxvelout << std::endl;
        std::cout << "spatial correlation (vorticity): " << i << std::endl;
        //spatial vorticity
        for (auto c : correlationCisnerosXVor(flows[i], smallwinsize, smallwinsize))
        {
            corxvorout << c << " ";
        }
        corxvorout << std::endl;
    }

    for (unsigned int x = 0; x < flows[0].cols; x += smallwinsize)
    {
        for (unsigned int y = 0; y < flows[0].rows; y += smallwinsize)
        {
            std::cout << "temporal correlation (velocity): " << x << " " << y << std::endl;
            //temporal velocity
            for (auto c : correlationCisnerosTVel(flows, x, y))
            {
                cortvelout << c << " ";
            }
            cortvelout << std::endl;
        }
    }

    for (unsigned int x = smallwinsize; x < flows[0].cols - smallwinsize; x += smallwinsize)
    {
        for (unsigned int y = smallwinsize; y < flows[0].rows - smallwinsize; y += smallwinsize)
        {
            std::cout << "temporal correlation (vorticity): " << x << " " << y << std::endl;
            //temporal vorticity
            for (auto c : correlationCisnerosTVor(flows, x, y, smallwinsize))
            {
                cortvorout << c << " ";
            }
            cortvorout << std::endl;
        }
    }


    std::cout << "Evaluation complete." << std::endl;

    corxvelout.close();
    corxvorout.close();
    cortvelout.close();
    cortvorout.close();
}
