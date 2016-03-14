// 2015-09-04

#include <iostream>
#include <map>
#include <tuple>

#include "core.hpp"
#include "highgui.hpp"
#include "imgproc.hpp"
#include "imgcodecs.hpp"
#include "imgcodecs/imgcodecs_c.h"

int main(int argc, char** argv)
{
    std::string filename = argv[1];

    cv::Mat original, normalized;

    original = cv::imread(filename, CV_LOAD_IMAGE_GRAYSCALE); //cv::imread("img_000000000_02-BF1_000.tif", CV_LOAD_IMAGE_GRAYSCALE);

    cv::equalizeHist(original, normalized);

    //cv::Ptr<cv::CLAHE> clahe = cv::createCLAHE();
    //clahe->setClipLimit(4);
    //cv::Mat dst;
    //clahe->apply(original, normalized);

    cv::imwrite(filename+".norm.png", normalized);

}
