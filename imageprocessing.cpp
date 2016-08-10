#include <iostream>
#include <cstdio>
#include <fstream>
#include <map>
#include <tuple>
#include <iomanip>
#include <sstream>
#include <stdexcept>

#include "core.hpp"
#include "highgui.hpp"
#include "imgproc.hpp"
#include "imgcodecs.hpp"
#include "imgcodecs/imgcodecs_c.h"

#include "Vector2D.hpp"

//A: test each rectangle only with the first in each group
//B: test each rectangle with each other: if any two are found to not belong to the same particle, it is considered as such //B is senseless.... use A or C....
//C: test each rectangle with each other: if all pairs are found to not belong to the same particle, it is considered as such
const int GROUP_MODE_A = 1;
const int GROUP_MODE_B = 2;
const int GROUP_MODE_C = 3;
const int GROUP_MODE = GROUP_MODE_C;

typedef double real; //float, double, ...

typedef Vector2DT<real> Vector;
typedef std::vector<cv::RotatedRect> RectangleList;

const cv::Scalar WHITE = cv::Scalar(255);
const cv::Scalar GRAY = cv::Scalar(127);
const cv::Scalar BLACK = cv::Scalar(0);
const cv::Scalar RED = cv::Scalar(0, 0, 255);
const cv::Scalar BLUE = cv::Scalar(255, 0, 0);
const cv::Scalar GREEN = cv::Scalar(0, 255, 0);
const cv::Scalar YELLOW = cv::Scalar(0, 255, 255);
const cv::Scalar ORANGE = cv::Scalar(0, 128, 255);

void drawRotatedRectangle(cv::Mat &image, const cv::RotatedRect& rect, const cv::Scalar& color, bool filled = false, unsigned int thickness = 1)
{
    cv::Point2f vertices2f[4];
    cv::Point vertices[4];
    rect.points(vertices2f);

    if (filled)
    {
        for (int i = 0; i < 4; ++i)
        {
            vertices[i] = vertices2f[i];
        }

        cv::fillConvexPoly(image, vertices, 4, color); //was: v::floodFill(image, rect.center, BLACK);
    }
    else
    {
        for (int i = 0; i < 4; ++i)
        {
            cv::line(image, vertices2f[i], vertices2f[(i+1)%4], color, thickness);
        }
    }
}

cv::Mat CreateRectangle(int width, int height, real angle)
{
    cv::RotatedRect rrect(cv::Point2f(0,0), cv::Size2f(width,height), angle);
    cv::Rect brect = rrect.boundingRect();

    rrect.center = cv::Point2f(brect.width/2, brect.height/2);

    cv::Mat image(brect.height, brect.width, CV_8UC1, WHITE);
    drawRotatedRectangle(image, rrect, BLACK, true); //opencv cannot draw rotated rectangles :(

    return image;
}

int MatchRectangles(const cv::Mat &image, const cv::Mat &rect, int x, int y, real ignore_ratio, int width, int height)
{
    int subx, suby;
    int false_pixels = 0;

    //int area_pixels = rect.cols*rect.rows;
    //int ignore_pixels = int(ignore_ratio*width*height - width*height + area_pixels);
    int ignore_pixels = int(ignore_ratio*width*height);

    for (suby = 0; suby < rect.rows; ++suby)
    {
        for (subx = 0; subx < rect.cols; ++subx)
        {
            false_pixels += int(!(image.at<unsigned char>(y+suby, x+subx) || rect.at<unsigned char>(suby, subx))); //both black pixels

            if (false_pixels > ignore_pixels)
            {
                return -1;
            }
        }
    }
    //std::cout << "ignore: " << ignore_pixels << ", false: " << false_pixels << std::endl;
    return false_pixels;
}

RectangleList FindRectangles(const cv::Mat &image, int width, int height, real delta_angle, real ignore_ratio, int delta_x, int delta_y)
{
    int x, y;
    real angle;
    //int ignore_pixels = int(width*height*ignore_ratio); //was: int(image.cols*image.rows*ignore_ratio)
    RectangleList ret;
    bool found;
    //bool last;

    int uneq;

    cv::Mat show;
    //PositionedRectangle posrect;

    std::cout << "Initialization complete." << std::endl;

    for (angle = 0.; angle < 180.; angle += delta_angle)
    {
        std::cout << "angle " << angle << std::endl;
        cv::Mat rect = CreateRectangle(width, height, angle);


        for (y = 0; y < image.rows - rect.rows + 1; y += delta_y)
        {
            for (x = 0; x < image.cols - rect.cols + 1; x += delta_x)
            {

                //std::cout << "x: " << x << ", y: " << y << ", rect.cols: " << rect.cols << ", rect.rows: " << rect.rows << std::endl;
                /*cv::bitwise_and(image(cv::Range(y, y+rect.rows), cv::Range(x, x+rect.cols)), rect, show);
                cv::imshow("Image", show);
                cv::waitKey();*/


                /*cv::Mat img_both = cv::Mat(rect.rows, 2 * rect.cols + 50, CV_8UC1, cv::Scalar(127));
                cv::Mat left(img_both, cv::Rect(0, 0, rect.cols, rect.rows)); // Copy constructor
                rect.copyTo(left);
                cv::Mat right(img_both, cv::Rect(rect.cols + 50, 0, rect.cols, rect.rows)); // Copy constructor
                image(cv::Range(y, y+rect.rows), cv::Range(x, x+rect.cols)).copyTo(right);
                cv::imshow("Both", img_both);
                uneq = MatchRectangles(image, rect, x, y, ignore_ratio, width, height);
                //std::cout << "ignore: " << int(ignore_ratio*width*height - width*height + rect.cols*rect.rows) << " uneq: " << uneq << std::endl;
                cv::waitKey();*/


                uneq = MatchRectangles(image, rect, x, y, ignore_ratio, width, height);
                found = uneq >= 0;


                /*cv::Mat image_copy = image.clone();
                cvtColor(image_copy, image_copy, CV_GRAY2BGR);
                cv::Point2f vertices[4];
                cv::RotatedRect(cv::Point2f(x+rect.cols/2, y+rect.rows/2), cv::Size2f(width, height), angle).points(vertices);
                for (int i = 0; i < 4; i++)
                {
                    cv::line(image_copy, vertices[i], vertices[(i+1)%4], cv::Scalar(0,0,255));
                }
                uneq = MatchRectangles(image, rect, x, y, ignore_ratio, width, height);
                found = uneq >= 0;
                cv::imshow("Orig_with_scanning_rect", image_copy);
                cv::waitKey();*/


                //CARE: found = MatchRectangles(image, rect, x, y, ignore_ratio, width, height) >= 0;
                /*if (found)
                {
                    if (last)
                    {
                        posrect = MergePositionedRectangles(posrect, PositionedRectangle{x, y, width, height, angle});
                    }
                    else
                    {
                        posrect = PositionedRectangle{x, y, width, height, angle};
                    }
                }
                else // found == false
                {
                    if (last)
                    {
                        ret.push_back(posrect); //no touching rectangle found: append the last merged rectangle
                    }
                }
                last = found;*/

                /*if (found)
                {
                    ret.push_back(PositionedRectangle{x, y, width, height, angle});
                }*/
                if (found)
                {
                    ret.push_back(cv::RotatedRect(cv::Point2f(x+rect.cols/2, y+rect.rows/2), cv::Size2f(width, height), angle));
                }
            }
            //if rectangle in the last corner: append
            /*if (last)
            {
                ret.push_back(posrect);
                last = false;
            }*/
        }

    }
    return ret;
}

void NumpySubtract(const cv::Mat &src1, const cv::Mat &src2, cv::Mat &dest)
{
    assert(src1.cols == src2.cols && src1.rows == src2.rows && src1.type() == src2.type());
    dest.create(src1.rows, src1.cols, src1.type());

    int x, y;

    for (y = 0; y < src1.rows; ++y)
    {
        for (x = 0; x < src1.cols; ++x)
        {
            dest.at<unsigned char>(y, x) = src1.at<unsigned char>(y, x) - src2.at<unsigned char>(y, x);
        }
    }
}

typedef std::vector<cv::Point> Contour;

bool RectInContour(const Contour &cont, const cv::RotatedRect &rect)
{
    //works because rect is either completly inside or completely outside of contour
    return cv::pointPolygonTest(cont, rect.center, false) != -1;
}

void drawWhiteRectanglesInBlackPicture(const RectangleList &rects, cv::Mat &image, int rows, int cols)
{
    //image.create(src1.rows, src1.cols, CV_8UC1, BLACK);
    image = cv::Mat(rows, cols, CV_8UC1, BLACK);

    for (auto rect : rects)
    {
        drawRotatedRectangle(image, rect, WHITE, true);
    }
}

real angleDiff(real angle1, real angle2)
{
    //case 0° and 359°
    return std::min(std::abs(angle1 - angle2), std::abs(angle1 - angle2 - 360.));
}

bool AreRectanglesInSameParticle(const cv::RotatedRect &rect1, const cv::RotatedRect &rect2, real max_allowed_angle_diff, real max_allowed_center_diff_factor, real max_allowed_perp_center_diff_factor)
{
    assert(rect1.size == rect2.size);
    real x_diff = rect1.center.x - rect2.center.x;
    real y_diff = rect1.center.y - rect2.center.y;
    Vector v(x_diff, y_diff);
    real center_distance = v.Norm();
    real width = rect1.size.width;
    real height = rect1.size.height;
    //angledifference of rectangle is too high to belong to one particle
    bool angleDiffTooHigh = (angleDiff(rect1.angle, rect2.angle) > max_allowed_angle_diff);
    if (angleDiffTooHigh)
    {
        return false;
    }
    //center of rectangles is too far appart to belong to one particle
    if (center_distance > height*max_allowed_center_diff_factor)
    {
        return false;
    }
    //angle difference is low (already first if), but perpendicular distance is high; idea: U - shaped 3 particle contours have to be differenciated
    Vector n(std::sin(rect1.angle/360.*2.*M_PI), -std::cos(rect1.angle/360.*2.*M_PI)); //assume rect1.angle \approx rect2.angle
    real perp_center_distance = (v - n*(n*v)).Norm();

    /*cv::Mat img_test = cv::Mat(1024, 1024, CV_8UC1, WHITE);

    drawRotatedRectangle(img_test, rect1, GRAY);
    drawRotatedRectangle(img_test, rect2, GRAY);
    std::cout << "n: " << n.x << ", " << n.y << std::endl;
    cv::line(img_test, cv::Point2f(rect1.center.x, rect1.center.y), cv::Point2f(rect1.center.x + v.x, rect1.center.y + v.y), BLACK);
    cv::line(img_test, cv::Point2f(rect1.center.x, rect1.center.y), cv::Point2f(rect1.center.x + n.x, rect1.center.y + n.y), BLACK);

    cv::imshow("Test", img_test);
    cv::waitKey();*/

    //std::cout << perp_center_distance << " / " << width*max_allowed_perp_center_diff_factor << std::endl;
    if (perp_center_distance > width*max_allowed_perp_center_diff_factor) //&& !angleDiffTooHigh
    {
        return false;
    }
    return true;
}

bool GroupIsSingleParticle(const RectangleList &rects, real max_allowed_angle_diff, real max_allowed_center_diff_factor, real max_allowed_perp_center_diff_factor)
{
    unsigned int i, j;
    for (i = 0; i < rects.size(); ++i)
    {
        for (j = i+1; j < rects.size(); ++j)
        {
            if (!AreRectanglesInSameParticle(rects[i], rects[j], max_allowed_angle_diff, max_allowed_center_diff_factor, max_allowed_perp_center_diff_factor))
            {
                return false;
            }
        }
    }
    return true;
}

std::vector<RectangleList> GroupParticles(const RectangleList &rects, real max_allowed_angle_diff, real max_allowed_center_diff_factor, real max_allowed_perp_center_diff_factor)
{
    std::vector<RectangleList> ret;
    unsigned int i, index;
    bool found;

    for (auto rect : rects)
    {
        found = false;
        for (i = 0; i < ret.size(); ++i)
        {
            if (GROUP_MODE == GROUP_MODE_A) //MODE a) test each rectangle only with the first in each group
            {
                if (AreRectanglesInSameParticle(rect, ret[i][0], max_allowed_angle_diff, max_allowed_center_diff_factor, max_allowed_perp_center_diff_factor))
                {
                    index = i;
                    found = true;
                }
            }
            else if (GROUP_MODE == GROUP_MODE_B) //MODE b) test each rectangle with each other: if any two are found to not belong to the same particle, it is considered as such
            {
                for (unsigned int j = 0; j < ret[i].size(); ++j)
                {
                    if (AreRectanglesInSameParticle(rect, ret[i][j], max_allowed_angle_diff, max_allowed_center_diff_factor, max_allowed_perp_center_diff_factor))
                    {
                        index = i;
                        found = true;
                    }
                    else
                    {
                        found = false;
                        break;
                    }
                }
            }
            else if (GROUP_MODE == GROUP_MODE_C) //MODE c) test each rectangle with each other: if all pairs are found to not belong to the same particle, it is considered as such
            {
                for (unsigned int j = 0; j < ret[i].size(); ++j)
                {
                    if (AreRectanglesInSameParticle(rect, ret[i][j], max_allowed_angle_diff, max_allowed_center_diff_factor, max_allowed_perp_center_diff_factor))
                    {
                        index = i;
                        found = true;
                        break;
                    }
                }
            }
        }
        if (found)
        {
            ret[index].push_back(rect);
        }
        else
        {
            RectangleList newgroup;
            newgroup.push_back(rect);
            ret.push_back(newgroup);
        }
    }

    return ret;
}

cv::Rect_<float> BoundingRect(const cv::RotatedRect &rect1, const cv::RotatedRect &rect2) //RotatedRect is float only
{
    std::vector<Vector> points(8);
    real a1 = rect1.angle/360.*2.*M_PI;
    real a2 = rect2.angle/360.*2.*M_PI;
    Vector c1 = Vector(rect1.center.x, rect1.center.y);
    Vector c2 = Vector(rect2.center.x, rect2.center.y);
    Vector vh1 = Vector(std::sin(a1), -std::cos(a1)) * rect1.size.height/2.;
    Vector vw1 = Vector(std::cos(a1), +std::sin(a1)) * rect1.size.width/2.;
    Vector vh2 = Vector(std::sin(a2), -std::cos(a2)) * rect2.size.height/2.;
    Vector vw2 = Vector(std::cos(a2), +std::sin(a2)) * rect2.size.width/2.;
    points[0] = c1 + vw1 + vh1;
    points[1] = c1 - vw1 + vh1;
    points[2] = c1 + vw1 - vh1;
    points[3] = c1 - vw1 - vh1;
    points[4] = c2 + vw2 + vh2;
    points[5] = c2 - vw2 + vh2;
    points[6] = c2 + vw2 - vh2;
    points[7] = c2 - vw2 - vh2;

    real left = points[0].x, right = points[0].x, top = points[0].y, bottom = points[0].y;
    for (unsigned int i = 1; i < 8; ++i)
    {
        left = std::min(left, points[i].x);
        right = std::max(right, points[i].x);
        top = std::min(top, points[i].y);
        bottom = std::max(bottom, points[i].y);
    }

    return cv::Rect_<float>(cv::Point_<float>(left, top), cv::Point_<float>(right, bottom));
}

unsigned int CountPixels(const cv::Mat &image, unsigned char color)
{
    int x, y;
    unsigned int ret = 0;

    for (y = 0; y < image.rows; ++y)
    {
        for (x = 0; x < image.cols; ++x)
        {
            if (image.at<unsigned char>(y, x) == color)
            {
                ret += 1;
            }
        }
    }

    return ret;
}

std::tuple<real, real, real> IntersectionArea(cv::RotatedRect rect1, cv::RotatedRect rect2) //copies because we change .center fields
{
    cv::Rect_<float> bbox_ = BoundingRect(rect1, rect2);
    cv::Rect_<float> bbox(std::floor(bbox_.x), std::floor(bbox_.y), std::ceil(bbox_.width), std::ceil(bbox_.height));

    rect1.center -= bbox.tl();
    rect2.center -= bbox.tl();

    cv::Mat box1(bbox.height, bbox.width, CV_8UC1, BLACK);
    cv::Mat box2(bbox.height, bbox.width, CV_8UC1, BLACK);
    cv::Mat andbox(bbox.height, bbox.width, CV_8UC1, BLACK);
    drawRotatedRectangle(box1, rect1, WHITE, true);
    drawRotatedRectangle(box2, rect2, WHITE, true);
    cv::bitwise_and(box1, box2, andbox); //intersection, _or: union

    unsigned int area1 = CountPixels(box1, 255); //WHITE
    unsigned int area2 = CountPixels(box2, 255); //WHITE
    unsigned int area_and = CountPixels(andbox, 255); //WHITE

    //real areareal = rect1.size.height * rect1.size.width +rect2.size.height * rect2.size.width - area_and;
    //int areaint = area1 + area2 - area_and;

    /*std::cout << "intersection: " << area_and << " / " << rect1.size.height*rect1.size.width << " (" << area1 << "), "<< rect2.size.height*rect2.size.width << " (" << area2 << ")" << std::endl;
    cv::imshow("intersection1", box1);
    cv::waitKey();
    cv::imshow("intersection2", box2);
    cv::waitKey();
    cv::imshow("intersection3", andbox);
    cv::waitKey();*/

    return std::make_tuple(area1, area2, area_and);
}

std::vector<Contour> MakeContoursFromRectangles(const cv::Mat &normalized, const RectangleList &rects, int rows, int cols)
{
    std::vector<Contour> contours;

    cv::Mat customcontourimage;
    drawWhiteRectanglesInBlackPicture(rects, customcontourimage, rows, cols);
    ///cv::imwrite("test_FilledRects.png", customcontourimage);

    cv::findContours(customcontourimage, contours, cv::noArray(), CV_RETR_LIST, CV_CHAIN_APPROX_SIMPLE);
    cv::Mat img_with_contours = normalized.clone();
    cvtColor(img_with_contours, img_with_contours, CV_GRAY2BGR);
    cv::drawContours(img_with_contours, contours, -1, RED);
    ///cv::imwrite("test_Orig_wContours.png", img_with_contours);

    return contours;
}

std::vector<RectangleList> GroupRectangles(const RectangleList &rects, const std::vector<Contour> &contours)
{
    std::vector<RectangleList> groups;

    real rects_per_contour_avg = 0;
    for (auto cont : contours)
    {
        RectangleList rects_in_cont;
        for (auto rect : rects)
        {
            if (RectInContour(cont, rect))
            {
                rects_in_cont.push_back(rect);
            }
        }
        rects_per_contour_avg += rects_in_cont.size();
        groups.push_back(rects_in_cont);
    }
    rects_per_contour_avg /= groups.size();
    std::cout << "Average number of rectangles per contour: " << rects_per_contour_avg << std::endl;

    return groups;
}

unsigned int FindNonSingleGroups(const std::vector<RectangleList> &groups, real max_allowed_angle_diff, real max_allowed_center_diff_factor, real  max_allowed_perp_center_diff_factor)
{
    unsigned int num_not_single = 0;
    for (auto group : groups)
    {
        if (!GroupIsSingleParticle(group, max_allowed_angle_diff, max_allowed_center_diff_factor, max_allowed_perp_center_diff_factor))
        {
            num_not_single += 1;
        }
    }

    return num_not_single;
}

cv::Mat DrawRotatedRectList(cv::Mat &image, const RectangleList &rects, const cv::Scalar &color, unsigned int thickness = 1)
{
    cv::Mat output = image.clone();
    cvtColor(output, output, CV_GRAY2BGR);

    for (auto rect : rects)
    {
        drawRotatedRectangle(output, rect, color, false, thickness);
    }

    return output;
}

cv::Mat DrawRotatedRectListList(cv::Mat &image, const std::vector<RectangleList> &rectslist, const cv::Scalar &color)
{
    cv::Mat output = image.clone();
    cvtColor(output, output, CV_GRAY2BGR);

    for (auto rects : rectslist)
    {
        for (auto rect : rects)
        {
            drawRotatedRectangle(output, rect, color, false);
        }
    }

    return output;
}

template<class InputIt, class UnaryFunction> std::ostream &PrintIterable(std::ostream &stream, InputIt iter, UnaryFunction f)
{
    auto it = iter.begin();
    auto last = iter.end();
    for (;;)
    {
        f(stream, *it);
        stream << " ";
        ++it;
        if (std::next(it) == last)
        {
            break;
        }
    }
    f(stream, *it);
    return stream;
}

void test_PrintIterable(void)
{
    std::vector<int> asd = {1, 2, 3};
    std::cout << "\""; PrintIterable(std::cout, asd, [](std::ostream &stream, int i)
        {
            stream << i;
        }
    ) << "\"" << std::endl;
    std::cout << "\""; PrintIterable(std::cout, asd, [](std::ostream &stream, int i)
        {
            stream << i;
        }
    ) << "\"" << std::endl;
}

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

int main(int argc, char **argv)
{
//    if (argc != 10)
//    {
//        std::cerr << "inpath(string)\tfilenametpl(string)\tstart(int)\tstop(int)\tX\tY\tDX\tDY\toutpath(string)" << std::endl;
//        throw std::runtime_error("Parameters missing");
//    }

	std::string INPATH = std::string("/Users/nikita/Documents/spr/20160809/100x_10BF_4um_100/Pos0/");argv[1]; // path with trailing slash
	std::string FILENAMETPL = std::string("img_%09d_01-BF0_000.tif");argv[2]; // "img_%09d_00-BF_EGFP_000.tif"
	int START = 0;std::atoi(argv[3]); // 0
	int STOP = 15;std::atoi(argv[4]); // 1000
	int X = 0;512;std::atoi(argv[5]); //512
	int Y = 575;325;std::atoi(argv[6]); //325
	int DX = 1024;std::atoi(argv[7]); //512
	int DY = 1024;std::atoi(argv[8]); //512
	std::string OUTPATH = std::string("/Users/nikita/Documents/spr/20160809/100x_10BF_4um_100/Pos0_output/");argv[9]; // path with trailing slash
	
	std::vector<std::string> filenames = CreateFilenameList(FILENAMETPL.c_str(), START, STOP);

    cv::Rect subrect = cv::Rect(X, Y, DX, DY);

    std::ofstream out, out_len, out_params;
    out.open(OUTPATH + "data_exp.txt");
    out_len.open(OUTPATH + "data_len.txt");
    out_params.open(OUTPATH + "data_params.txt");

/*for (real anglee = 0.; anglee < 360.; anglee+=1.)
{
    cv::RotatedRect rect1(cv::Point2f(800, 200), cv::Size(100, 200), anglee), rect2(cv::Point2f(500, 200), cv::Size(100, 200), anglee);
    real x_diff = rect1.center.x - rect2.center.x;
    real y_diff = rect1.center.y - rect2.center.y;
    Vector v(x_diff, y_diff);
    Vector n(std::sin(rect1.angle/360.*2.*M_PI),-std::cos(rect1.angle/360.*2.*M_PI)); //assume rect1.angle \approx rect2.angle
    Vector perp = (v - n*(n*v));

    cv::Mat img_test = cv::Mat(1024, 1024, CV_8UC1, WHITE);

    drawRotatedRectangle(img_test, rect1, GRAY);
    drawRotatedRectangle(img_test, rect2, GRAY);
    std::cout << "n: " << n.x << ", " << n.y << std::endl;
    cv::line(img_test, cv::Point2f(rect2.center.x, rect2.center.y), cv::Point2f(rect2.center.x + v.x, rect2.center.y + v.y), BLACK);
    cv::line(img_test, cv::Point2f(rect2.center.x, rect2.center.y), cv::Point2f(rect2.center.x + 100*n.x, rect2.center.y + 100*n.y), BLACK);
    cv::line(img_test, cv::Point2f(rect2.center.x, rect2.center.y), cv::Point2f(rect2.center.x + perp.x, rect2.center.y + perp.y), BLACK);

    cv::imshow("Test", img_test);
    cv::waitKey();
}*/

	//parameters for finding the rectangles //1 20 1 0.08 2 2 //1 30 1 0.12 2 2
	int width = 10;2; //3; //5; //of test rectangle >=1
	int height = 50; //40; //80; //of test rectangle =>1
	real delta_angle = 1; //1.; //20.; //1.; //rotation of rectangle >0, <180 // in degree
	real ignore_ratio = 0.026; //0.00; //0.015; //maximum area error in binary picture <=1, >=0
	int delta_x = 1; //2; //1; //translation of rectangle in x >=1
	int delta_y = 1; //2; //1; //translation of rectangle in y >=1
	//parameters for rectangle grouping in contours
	real max_allowed_angle_diff = 2.; //18.; //in degree
	real max_allowed_center_diff_factor = 0.6; // 1.2; //factor of "height"
	real max_allowed_perp_center_diff_factor = 0.4; // 8.; //factor of "width"
	//parameters for filtering out additional overlapping rects in some contour
	real max_area_overlapp = 0.52; // 0.7; //factor of area of smaller of 2 compared rects
	//GroupParticles has 2-modes (a, b) where a is better but slower and needs more tweaking (mainly overlapp parameter)
	
	cv::Mat original, grey, normalized, blurred, subtracted, subtracted_normalized, edge_detected, binary, dilatederoded;

    //cv::Rect subrect(360,  90, 512, 512); // subrect(512, 325, 512, 512) for 60x_1.0_BF_2016_01_27_2, subrect(380, 150, 512, 512) for 60x_1.0_BF_2016_01_27_4, subrect(360,  90, 512, 512) for 60x_1.0_BF_2016_01_27_6
    bool invert_from_binary_onwards = false;

    out << "W: " << subrect.width << " H: " << subrect.height << std::endl;
    out_len << "W: " << subrect.width << " H: " << subrect.height << std::endl;

    out_params << "group_mode = " << GROUP_MODE << std::endl
    << "subrect = " << subrect.x << ", " << subrect.y << ", " << subrect.width << ", " << subrect.height << std::endl
    << "width = " << width << std::endl
    << "height = " << height << std::endl
    << "delta_angle = " << delta_angle << std::endl
    << "ignore_ratio = " << ignore_ratio << std::endl
    << "delta_x = " << delta_x << std::endl
    << "delta_y = " << delta_y << std::endl
    << "max_allowed_angle_diff = " << max_allowed_angle_diff << std::endl
    << "max_allowed_center_diff_factor = " << max_allowed_center_diff_factor << std::endl
    << "max_allowed_perp_center_diff_factor = " << max_allowed_perp_center_diff_factor << std::endl
    << "max_area_overlapp = " << max_area_overlapp << std::endl
    << "invert = " << invert_from_binary_onwards << std::endl;
    out_params.close();

    unsigned int imagenum = 0;
    for (std::string filename : filenames)
    {
        std::stringstream ss;
        ss << std::setfill('0') << std::setw(5);
        ss << imagenum;

        //load original image
        original = cv::imread(INPATH + filename, CV_LOAD_IMAGE_GRAYSCALE); //("1.0ROIwhite9.9.15-gut.tif", CV_LOAD_IMAGE_GRAYSCALE); //cv::imread("img_000000000_02-BF1_000.tif", CV_LOAD_IMAGE_GRAYSCALE);

        //subimage
        cv::Mat subImage(original, subrect);
        original = subImage;

        //convert to greyscale
        //cvtColor(original, grey, cv::CV_BGR2GRAY);
        //cv::imwrite(OUTPATH + ss.str() + "test_Grey.png", grey);

        //histogram equalizer
        cv::equalizeHist(original, normalized);
        //cv::Ptr<cv::CLAHE> clahe = cv::createCLAHE();
        //clahe->setClipLimit(4);
        //clahe->apply(original, normalized);
        cv::imwrite(OUTPATH + ss.str() + "_Normalized.png", normalized);

        //bilateral filter
        cv::bilateralFilter(normalized, blurred, 20, 100, 100); //20, 100, 100
        ///cv::imwrite(OUTPATH + ss.str() + "test_Blurred.png", blurred);

        //subtract
        //subtracted = blurred - original;
        //cv::subtract(normalized, blurred, subtracted);
        NumpySubtract(blurred, normalized, subtracted);
        //NumpySubtract(normalized, blurred, subtracted);
        ///cv::imwrite(OUTPATH + ss.str() + "test_Subtracted.png", subtracted);

        //histogram equalizer
        cv::equalizeHist(subtracted, subtracted_normalized);
        ///cv::imwrite(OUTPATH + ss.str() + "test_Subtracted_Normalized.png", subtracted_normalized);

        //edge detection
        //cv::Canny(subtracted_normalized, edge_detected, 0, 50, 7);
        //cv::imwrite(OUTPATH + ss.str() + "test_Edge_Detected.png", edge_detected);
        //CONTINUE /W EDGE_DETECTED?

        //binary
        //CARE: because numpy-subtract was used, we binarize subtracted instead of subtracted_normalized!
        cv::threshold(subtracted, binary, 127, 255, 0); //127, 255, 0
        //cv::adaptiveThreshold(normalized, binary, 255, cv::ADAPTIVE_THRESH_MEAN_C, cv::THRESH_BINARY, blocksize, C);

        //optional invert
        if (invert_from_binary_onwards)
        {
            cv::bitwise_not(binary, binary);
        }
        ///cv::imwrite(OUTPATH + ss.str() + "test_Binary.png", binary);

        //dilate/erode
        //cv::dilate(binary, dilatederoded, cv::Mat::ones(3,3,CV_8U), cv::Point(-1, -1), 1);
        //cv::imwrite(OUTPATH + ss.str() + "test_Binary_eroded_dilated.png", dilatederoded);
        dilatederoded = binary;

        RectangleList rects = FindRectangles(dilatederoded, width, height, delta_angle, ignore_ratio, delta_x, delta_y);

        cv::Mat img_with_rectangles = DrawRotatedRectList(normalized, rects, RED);
        ///cv::imwrite(OUTPATH + ss.str() + "test_Orig_wRectangles.png", img_with_rectangles);

        std::vector<Contour> contours = MakeContoursFromRectangles(normalized, rects, original.rows, original.cols);
        std::cout << "Found rectangles: " << rects.size() << ", contours: " << contours.size() << std::endl;

        std::vector<RectangleList> groups = GroupRectangles(rects, contours);

        unsigned int num_not_single = FindNonSingleGroups(groups, max_allowed_angle_diff, max_allowed_center_diff_factor, max_allowed_perp_center_diff_factor);
        std::cout << "Not single particle contours (maybe more than 1 new particle per new group): " << num_not_single << std::endl;

        int num_particles = 0;
        std::vector<Contour> subcontours;
        cv::Mat image_ellipse_fit;

        std::vector<RectangleList> fit_rects;
        std::vector<std::vector<real>> fit_rects_weights;
        //RectangleList fit_ellipses;

        //fit rectangles in grouped contours
        for (auto group : groups)
        {
            std::vector<RectangleList> subgroups = GroupParticles(group, max_allowed_angle_diff, max_allowed_center_diff_factor, max_allowed_perp_center_diff_factor);
            //num_particles += subgroups.size();

            //each subgroup is now a single particle
            RectangleList fit_subrects;
            std::vector<real> fit_subrects_weights;
            for (auto subgroup : subgroups)
            {
                //TODO: merge rectangles and calc bounding box: this is now a single particle

                //very unperformant alternativ

                drawWhiteRectanglesInBlackPicture(subgroup, image_ellipse_fit, original.rows, original.cols);
                cv::findContours(image_ellipse_fit, subcontours, cv::noArray(), CV_RETR_LIST, CV_CHAIN_APPROX_SIMPLE); //modifies image
                //cvtColor(image_ellipse_fit, image_ellipse_fit, CV_GRAY2BGR);
                //cv::drawContours(image_ellipse_fit, subcontours, -1, RED);
                /*if (subcontours.size() > 1.)
                {
                    drawWhiteRectanglesInBlackPicture(subgroup, image_ellipse_fit, original.rows, original.cols);
                    cv::imshow("image_ellipse_fit", image_ellipse_fit);
                    std::cout << subcontours.size() << std::endl;
                    cv::waitKey();
                }*/
                //assert(subcontours.size() == 1); //does fail sometimes!

                for (Contour subcont : subcontours)
                {
                    //filters out little points that come from superposition of 3 "bonfire" particles with hole in middle
                    if (cv::contourArea(subcont) >= 0.9 * width * height)
                    {
                        num_particles += 1;
                        //minRect = cv::minAreaRect( cv::Mat(subcont) );
                        fit_subrects.push_back(cv::minAreaRect(subcont));
                        fit_subrects_weights.push_back(real(subgroup.size())/subcontours.size()); //weight
                    }
                }

                //BUG: if RotatedRect 0°, 90°, 180°, 270° -> consists only of 4 points -> will not get fitted
                //if (subcontours[0].size() > 5)
                //{
                //    minEllipse = cv::fitEllipse( cv::Mat(subcontours[0]) );
                //}

                //fit_rects.push_back(minRect);
                //fit_ellipses.push_back(minEllipse);
            }
            fit_rects.push_back(fit_subrects);
            fit_rects_weights.push_back(fit_subrects_weights);
        }
        //assert will fail if there are any RotatedRect 0°, 90°, ... -> consists only of 4 points -> will not get fitEllipse-fitted
        //assert(fit_rects.size() == fit_ellipses.size());
        std::cout << "Particles (fitted rectangles): " << num_particles << std::endl;
        //std::cout << "Fitted Rectangles: " << fit_rects.size() << std::endl;

        // Draw rotated rects (non-filtered)
        cv::Mat img_with_fits_before_filter = DrawRotatedRectListList(normalized, fit_rects, RED);
        ///cv::imwrite(OUTPATH + ss.str() + "test_Orig_wfittedRects_beforeFilter.png", img_with_fits_before_filter);


        // filter subrects in group that overlapp
        unsigned int i, j, index, num_filtered = 0;
        real area_smaller_rect; //area, area_i, area_j
        RectangleList fit_rects_filtered;
        for (auto group : fit_rects)
        {
            std::vector<int> dismissed(group.size(), 0);
            for (i = 0; i < group.size(); ++i)
            {
                for (j = i + 1; j < group.size(); ++j)
                {
                    std::tuple<real, real, real> area_tuple = IntersectionArea(group[i], group[j]);
                    if (std::get<2>(area_tuple) > 0.)
                    {
                        //area_i = group[i].size.width*group[i].size.height;
                        //area_j = group[j].size.width*group[j].size.height;
                        if (std::get<0>(area_tuple) < std::get<1>(area_tuple))
                        {
                            area_smaller_rect = std::get<0>(area_tuple);
                            index = i;
                        }
                        /*//fix for bug: area equal, then both got dismissed; with this only the one with the smaller index gets dismissed; lower index rect gets dismissed increased by 2
                        else if (std::get<0>(area_tuple) == std::get<1>(area_tuple))
                        {
                            area_smaller_rect = std::get<0>(area_tuple); //or std::get<1>()
                            index = std::min(i, j); // or std::max()
                            std::cout << "Same rectangle area: here was a bug before" << std::endl;
                        }*/
                        else
                        {
                            area_smaller_rect = std::get<1>(area_tuple);
                            index = j;
                        }
                        assert(area_smaller_rect > 0);
                        if (std::get<2>(area_tuple) > max_area_overlapp * area_smaller_rect)
                        {
                            dismissed[index] += 1;
                        }
                    }
                }
            }

            //new list with only not-dismissed rects
            for (i = 0; i < group.size(); ++i)
            {
                if (dismissed[i] == 0) //was not even once kicked out
                {
                    fit_rects_filtered.push_back(group[i]);
                }
                else
                {
                    num_filtered += 1;
                }
            }
        }

        PrintIterable(out, fit_rects_filtered, [&](std::ostream &stream, cv::RotatedRect rect)
            {
                stream << rect.center.x << " " << original.rows - rect.center.y;
            }
        ) << std::endl;
        PrintIterable(out_len, fit_rects_filtered, [](std::ostream &stream, cv::RotatedRect rect)
            {
                if (rect.size.width < rect.size.height)
                {
                    stream << rect.size.height << " " << rect.size.width;
                }
                else
                {
                    stream << rect.size.width << " " << rect.size.height;
                }
            }
        ) << std::endl;
        PrintIterable(out, fit_rects_filtered, [](std::ostream &stream, cv::RotatedRect rect)
            {
                if (rect.size.width < rect.size.height)
                {
                    stream << ((90. - rect.angle)/360.*2*M_PI);
                }
                else
                {
                    stream << ((-rect.angle)/360.*2*M_PI);
                }
            }
        ) << std::endl;

        std::cout << "Filtered additional rectangles due to overlap: " << num_filtered << std::endl;
        std::cout << "Final count of particles: " << fit_rects_filtered.size() << std::endl;

        // Draw rotated rects (filtered)
        cv::Mat img_with_fits = DrawRotatedRectList(normalized, fit_rects_filtered, RED, 2);
        cv::imwrite(OUTPATH + ss.str() + "_Orig_wfittedRects.png", img_with_fits);
        ++imagenum;
    }

    out.close();
    out_len.close();
    std::cout << "Evaluation complete." << std::endl;
    return 0;
}
