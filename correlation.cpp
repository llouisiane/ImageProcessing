#include "correlation.hpp"

real centralDifference(real a, real b, real delta)
{
    return (a - b)/delta/2.;
}

//z component of curl (rotation)
real localCurl(const cv::Mat &flow, unsigned int x, unsigned int y, unsigned int flowstep)
{
    const cv::Point2f& vt = flow.at<cv::Point2f>(y-flowstep, x);
    const cv::Point2f& vb = flow.at<cv::Point2f>(y+flowstep, x);
    const cv::Point2f& vl = flow.at<cv::Point2f>(y, x-flowstep);
    const cv::Point2f& vr = flow.at<cv::Point2f>(y, x+flowstep);
    // dvy/dx - dvx/dy
    return centralDifference(vl.y, vr.y, flowstep) - centralDifference(vt.x, vb.x, flowstep);
}

////////// Fluid dynamics of self-propelled microorganisms, from individuals to concentrated populations (Cisneros)
////////// Collective Motion of Spherical Bacteria (Pone)
/// spatial velocity
// <v(x)**2>
real correlationCisnerosXV2M(const cv::Mat &flow, unsigned int flowstep)
{
    unsigned int num = 0, x, y;
    real sum = 0.;
    for (x = 0; x < flow.cols; x += flowstep)
    {
        for (y = 0; y < flow.rows; y += flowstep)
        {
            const cv::Point2f& v = flow.at<cv::Point2f>(y, x);
            sum += Vector(v.x, v.y)*Vector(v.x, v.y);
            num += 1;
        }
    }
    return num == 0 ? 0. : sum / num;
}

// <v(x)>**2
real correlationCisnerosXVM2(const cv::Mat &flow, unsigned int flowstep)
{
    unsigned int num = 0, x, y;
    Vector sum(0., 0.);
    for (x = 0; x < flow.cols; x += flowstep)
    {
        for (y = 0; y < flow.rows; y += flowstep)
        {
            const cv::Point2f& v = flow.at<cv::Point2f>(y, x);
            sum += Vector(v.x, v.y);
            num += 1;
        }
    }
    return num == 0 ? 0. : (sum / num)*(sum / num);
}

// <v(x)v(x+r)>
real correlationCisnerosXVRM(const cv::Mat &flow, real r, real rstep, unsigned int flowstep)
{
    unsigned int num = 0, x, y, xr, yr;
    real n, sum = 0.;
    for (x = 0; x < flow.cols; x += flowstep)
    {
        for (y = 0; y < flow.rows; y += flowstep)
        {
            for (xr = 0; xr < flow.cols; xr += flowstep)
            {
                for (yr = 0; yr < flow.rows; yr += flowstep)
                {
                    n = (Vector(x, y) - Vector(xr, yr)).Norm();
                    if (n >= r && n < r + rstep)
                    {
                        const cv::Point2f& v = flow.at<cv::Point2f>(y, x);
                        const cv::Point2f& vr = flow.at<cv::Point2f>(yr, xr);
                        sum += Vector(v.x, v.y)*Vector(vr.x, vr.y);
                        num += 1;
                    }
                }
            }
        }
    }
    return num == 0 ? 0. : sum / num;
}

std::vector<real> correlationCisnerosXVel(const cv::Mat &flow, real flowstep, real rstep)
{
    real r, c;

    std::vector<real> ret;

    real v2m = correlationCisnerosXV2M(flow, flowstep);
    real vm2 = correlationCisnerosXVM2(flow, flowstep);

    for (r = 0.; r < flow.cols; r += rstep) //r should be able to go higher (\sqrt(2)*flow.clos)
    {
        c = (correlationCisnerosXVRM(flow, r, rstep, flowstep) - vm2) / (v2m - vm2);
        ret.push_back(c);
    }

    return ret;
}

std::vector<real> correlationPoneXVel(const cv::Mat &flow, real flowstep, real rstep)
{
    real r, c;

    std::vector<real> ret;

    for (r = 0.; r < flow.cols; r += rstep) //r should be able to go higher (\sqrt(2)*flow.clos)
    {
        unsigned int num = 0;
        real x, y, xr, yr;
        real n, sum = 0.;
        Vector sum1(0.,0.), sum2(0.,0.);
        for (x = 0; x < flow.cols; x += flowstep)
        {
            for (y = 0; y < flow.rows; y += flowstep)
            {
                for (xr = 0; xr < flow.cols; xr += flowstep)
                {
                    for (yr = 0; yr < flow.rows; yr += flowstep)
                    {
                        n = (Vector(x, y) - Vector(xr, yr)).Norm();
                        if (n >= r && n < r + rstep)
                        {
                            const cv::Point2f& v = flow.at<cv::Point2f>(y, x);
                            const cv::Point2f& vr = flow.at<cv::Point2f>(yr, xr);
                            sum += Vector(v.x, v.y)*Vector(vr.x, vr.y);
                            sum1 += Vector(v.x, v.y);
                            sum2 += Vector(vr.x, vr.y);
                            num += 1;
                        }
                    }
                }
            }
        }
        c = (num == 0 ? 0. : sum / num - sum1 * sum2 / (num * num));
        ret.push_back(c);
    }

    for (unsigned int i = 0; i < ret.size(); ++i)
    {
        ret[i] = ret[i] / ret[0]; // doesn't work, overwrites
    }

    return ret;
}

/// spatial vorticity
// <o(x)**2>
real correlationCisnerosXO2M(const cv::Mat &flow, unsigned int flowstep)
{
    unsigned int num = 0, x, y;
    real sum = 0., curl;
    for (x = flowstep; x < flow.cols - flowstep; x += flowstep)
    {
        for (y = flowstep; y < flow.rows - flowstep; y += flowstep)
        {
            curl = localCurl(flow, x, y, flowstep);
            sum += curl*curl;
            num += 1;
        }
    }
    return num == 0 ? 0. : sum / num;
}

// <o(x)>**2
real correlationCisnerosXOM2(const cv::Mat &flow, unsigned int flowstep)
{
    unsigned int num = 0, x, y;
    real sum = 0.;
    for (x = flowstep; x < flow.cols - flowstep; x += flowstep)
    {
        for (y = flowstep; y < flow.rows - flowstep; y += flowstep)
        {
            sum += localCurl(flow, x, y, flowstep);
            num += 1;
        }
    }
    return num == 0 ? 0. : (sum / num)*(sum / num);
}

// <o(x)o(x+r)>
real correlationCisnerosOXRM(const cv::Mat &flow, real r, real rstep, unsigned int flowstep)
{
    unsigned int num = 0, x, y, xr, yr;
    real n, sum = 0.;
    for (x = flowstep; x < flow.cols - flowstep; x += flowstep)
    {
        for (y = flowstep; y < flow.rows - flowstep; y += flowstep)
        {
            for (xr = flowstep; xr < flow.cols - flowstep; xr += flowstep)
            {
                for (yr = flowstep; yr < flow.rows - flowstep; yr += flowstep)
                {
                    n = (Vector(x, y) - Vector(xr, yr)).Norm();
                    if (n >= r && n < r + rstep)
                    {
                        sum += localCurl(flow, x, y, flowstep)*localCurl(flow, xr, yr, flowstep);
                        num += 1;
                    }
                }
            }
        }
    }
    return num == 0 ? 0. : sum / num;
}

std::vector<real> correlationCisnerosXVor(const cv::Mat &flow, real flowstep, real rstep)
{
    real r, c;

    std::vector<real> ret;

    real o2m = correlationCisnerosXO2M(flow, flowstep);
    real om2 = correlationCisnerosXOM2(flow, flowstep);

    for (r = 0.; r < flow.cols; r += rstep)
    {
        c = (correlationCisnerosOXRM(flow, r, rstep, flowstep) - om2) / (o2m - om2);
        ret.push_back(c);
    }

    return ret;
}

std::vector<real> correlationPoneXVor(const cv::Mat &flow, real flowstep, real rstep)
{
    real r, c;

    std::vector<real> ret;

    for (r = 0.; r < flow.cols; r += rstep)
    {
        unsigned int num = 0, x, y, xr, yr;
        real n, sum = 0., sum1 = 0., sum2 = 0.;
        for (x = flowstep; x < flow.cols - flowstep; x += flowstep)
        {
            for (y = flowstep; y < flow.rows - flowstep; y += flowstep)
            {
                for (xr = flowstep; xr < flow.cols - flowstep; xr += flowstep)
                {
                    for (yr = flowstep; yr < flow.rows - flowstep; yr += flowstep)
                    {
                        n = (Vector(x, y) - Vector(xr, yr)).Norm();
                        if (n >= r && n < r + rstep)
                        {
                            sum += localCurl(flow, x, y, flowstep)*localCurl(flow, xr, yr, flowstep);
                            sum1 += localCurl(flow, x, y, flowstep);
                            sum2 += localCurl(flow, xr, yr, flowstep);
                            num += 1;
                        }
                    }
                }
            }
        }
        c = (num == 0 ? 0. : sum / num - sum1 * sum2 / (num * num));
        ret.push_back(c);
    }

    for (unsigned int i = 0; i < ret.size(); ++i)
    {
        ret[i] = ret[i] / ret[0]; // doesn't work, overwrites
    }

    return ret;
}

/// temporal velocity

// <v(t)**2>
real correlationCisnerosTV2M(const std::vector<cv::Mat> &flows, unsigned int x, unsigned int y)
{
    unsigned int t, num = 0;
    real sum = 0.;
    for (t = 0; t < flows.size(); ++t)
    {
        const cv::Point2f& v = flows[t].at<cv::Point2f>(y, x);
        sum += Vector(v.x, v.y)*Vector(v.x, v.y);
        num += 1;
    }
    return num == 0 ? 0. : sum / num;
}

// <v(t)>**2
real correlationCisnerosTVM2(const std::vector<cv::Mat> &flows, unsigned int x, unsigned int y)
{
    unsigned int num = 0, t;
    Vector sum(0., 0.);
    for (t = 0; t < flows.size(); ++t)
    {
        const cv::Point2f& v = flows[t].at<cv::Point2f>(y, x);
        sum += Vector(v.x, v.y);
        num += 1;
    }
    return num == 0 ? 0. : (sum / num)*(sum / num);
}

// <v(t)v(t+s)>
real correlationCisnerosVTSM(const std::vector<cv::Mat> &flows, unsigned int dt, unsigned int x, unsigned int y)
{
    unsigned int num = 0, t, ts;
    real sum = 0.;
    for (t = 0; t < flows.size(); ++t)
    {
        for (ts = 0; ts < flows.size(); ++ts)
        {
            if (std::abs(ts-t) == dt)
            {
                const cv::Point2f& v = flows[t].at<cv::Point2f>(y, x);
                const cv::Point2f& vs = flows[ts].at<cv::Point2f>(y, x);
                sum += Vector(v.x, v.y)*Vector(vs.x, vs.y);
                num += 1;
            }
        }
    }
    return num == 0 ? 0. : sum / num;
}

std::vector<real> correlationCisnerosTVel(const std::vector<cv::Mat> &flows, unsigned int x, unsigned int y)
{
    real c;
    unsigned int t;

    std::vector<real> ret;

    real t2m = correlationCisnerosTV2M(flows, x, y);
    real tm2 = correlationCisnerosTVM2(flows, x, y);

    for (t = 0; t < flows.size(); ++t)
    {
        c = (correlationCisnerosVTSM(flows, t, x, y) - tm2) / (t2m - tm2);
        ret.push_back(c);
    }

    return ret;
}

std::vector<real> correlationPoneTVel(const std::vector<cv::Mat> &flows, unsigned int x, unsigned int y)
{
    real c;
    unsigned int dt;

    std::vector<real> ret;


    for (dt = 0; dt < flows.size(); ++dt)
    {
        unsigned int num = 0, t, ts;
        real sum = 0.;
        Vector sum1(0.,0.), sum2(0.,0.);
        for (t = 0; t < flows.size(); ++t)
        {
            for (ts = 0; ts < flows.size(); ++ts)
            {
                if (std::abs(ts-t) == dt)
                {
                    const cv::Point2f& v = flows[t].at<cv::Point2f>(y, x);
                    const cv::Point2f& vs = flows[ts].at<cv::Point2f>(y, x);
                    sum += Vector(v.x, v.y)*Vector(vs.x, vs.y);
                    sum1 += Vector(v.x, v.y);
                    sum2 += Vector(vs.x, vs.y);
                    num += 1;
                }
            }
        }

        c = (num == 0 ? 0. : sum / num - sum1 * sum2 / (num * num));
        ret.push_back(c);
    }

    for (unsigned int i = 0; i < ret.size(); ++i)
    {
        ret[i] = ret[i] / ret[0]; // doesn't work, overwrites
    }

    return ret;
}

/// temporal vorticity

// <o(t)**2>
real correlationCisnerosTO2M(const std::vector<cv::Mat> &flows, unsigned int x, unsigned int y, unsigned int flowstep)
{
    unsigned int num = 0, t;
    real sum = 0., curl;
    for (t = 0; t < flows.size(); ++t)
    {
        curl = localCurl(flows[t], x, y, flowstep);
        sum += curl*curl;
        num += 1;
    }
    return num == 0 ? 0. : sum / num;
}

// <o(t)>**2
real correlationCisnerosTOM2(const std::vector<cv::Mat> &flows, unsigned int x, unsigned int y, unsigned int flowstep)
{
    unsigned int num = 0, t;
    real sum = 0.;
    for (t = 0; t < flows.size(); ++t)
    {
        sum += localCurl(flows[t], x, y, flowstep);
        num += 1;
    }
    return num == 0 ? 0. : (sum / num)*(sum / num);
}

// <o(t)o(t+s)>
real correlationCisnerosOTSM(const std::vector<cv::Mat> &flows, unsigned int dt, unsigned int x, unsigned int y, unsigned int flowstep)
{
    unsigned int t, ts, num = 0;
    real sum = 0.;
    for (t = 0; t < flows.size(); ++t)
    {
        for (ts = 0; ts < flows.size(); ++ts)
        {
            if (std::abs(ts-t) == dt)
            {
                sum += localCurl(flows[t], x, y, flowstep)*localCurl(flows[ts], x, y, flowstep);
                num += 1;
            }
        }
    }
    return num == 0 ? 0. : sum / num;
}

std::vector<real> correlationCisnerosTVor(const std::vector<cv::Mat> &flows, unsigned int x, unsigned int y, unsigned int flowstep)
{
    real c;
    unsigned int t;

    std::vector<real> ret;

    real t2m = correlationCisnerosTO2M(flows, x, y, flowstep);
    real tm2 = correlationCisnerosTOM2(flows, x, y, flowstep);

    for (t = 0; t < flows.size(); ++t)
    {
        c = (correlationCisnerosOTSM(flows, t, x, y, flowstep) - tm2) / (t2m - tm2);
        ret.push_back(c);
    }

    return ret;
}

std::vector<real> correlationPoneTVor(const std::vector<cv::Mat> &flows, unsigned int x, unsigned int y, unsigned int flowstep)
{
    real c;
    unsigned int dt;

    std::vector<real> ret;

    for (dt = 0; dt < flows.size(); ++dt)
    {
        unsigned int t, ts, num = 0;
        real sum = 0., sum1 = 0., sum2 = 0.;
        for (t = 0; t < flows.size(); ++t)
        {
            for (ts = 0; ts < flows.size(); ++ts)
            {
                if (std::abs(ts-t) == dt)
                {
                    sum += localCurl(flows[t], x, y, flowstep)*localCurl(flows[ts], x, y, flowstep);
                    sum1 += localCurl(flows[t], x, y, flowstep);
                    sum2 += localCurl(flows[ts], x, y, flowstep);
                    num += 1;
                }
            }
        }

        c = (num == 0 ? 0. : sum / num - sum1 * sum2 / (num * num));
        ret.push_back(c);
    }

    for (unsigned int i = 0; i < ret.size(); ++i)
    {
        ret[i] = ret[i] / ret[0]; // doesn't work, overwrites
    }

    return ret;
}


//obsolete: uses manhattan metric (wrongly implemented from journal pone 0083760)
std::vector<real> correlationPone0083760(const cv::Mat &flow, const int stepsize)
{
    assert(flow.cols == flow.rows);
    //what? assert(flow.cols % 0 == 0); //otherwise not restlessly dividable
    std::vector<real> rs(flow.cols/stepsize);
    real phi2;
    std::vector<real> phix(2);
    std::vector<real> phiy(2);
    int num, isum, jsum;
    cv::Point2f rp, rp2;

    for (int r = 0; r < flow.cols/stepsize; ++r)
    {
        num = 0;
        phi2 = 0;
        phix[0] = 0; phix[1] = 0;
        phiy[0] = 0; phiy[1] = 0;
        for (int i = 0; i < flow.cols; i += stepsize)
        {
            for (int j = 0; j < flow.rows; j += stepsize)
            {
                for (int i2 = -(r + 1) * stepsize; i2 < (r + 1) * stepsize; i2 += stepsize)
                {
                    for (int j2 = -(r + 1) * stepsize; j2 < (r + 1) * stepsize; j2 += stepsize)
                    {
                        isum = i2 + i;
                        jsum = j2 + j;
                        if (isum >= 0 && isum < flow.cols && jsum >=0 && jsum < flow.rows)
                        {
                            //std::cout << i, j, i2, j2;
                            //std::cout << std::endl;
                            rp = flow.at<cv::Point2f>(j, i);
                            rp2 = flow.at<cv::Point2f>(jsum, isum);
                            phi2 += rp.x * rp2.x + rp.y * rp2.y;
                            phix[0] += rp.x;
                            phix[1] += rp.y;
                            phiy[0] +=  rp2.x;
                            phiy[1] +=  rp2.y;
                            num += 1;
                        }
                    }
                }
            }
        }
        if (num == 0)
        {
            //redundant
            phi2 /= 1;
            phix[0] /= 1;
            phix[1] /= 1;
            phiy[0] /= 1;
            phiy[1] /= 1;
        }
        else
        {
            phi2 /= num;
            phix[0] /= num;
            phix[1] /= num;
            phiy[0] /= num;
            phiy[1] /= num;
        }
        rs[r] = phi2 - (phix[0] * phiy[0] + phix[1] * phiy[1]);
        if (r > 0)
        {
            rs[r] /= rs[0];
        }
    }
    rs[0] = 1;
    return rs;
}

// VICSEK

real correlationCisnerosXV2Mvic(const std::vector<Vector> &velocity)
{
    real total = 0.;
    int num = 0;
    for (auto v : velocity)
    {
        total += v*v;
        num++;
    }
    return num == 0 ? 0. : total / num;
}

real correlationCisnerosXVM2vic(const std::vector<Vector> &velocity)
{
    Vector total(0., 0.);
    int num = 0;
    for (auto v : velocity)
    {
        total += v;
        num++;
    }
    return num == 0 ? 0. : (total / num)*(total / num);
}

real correlationCisnerosXVRMvic(const std::vector<Vector> &position, const std::vector<Vector> &velocity, real r, real rstep)
{
    assert(position.size() == velocity.size());
    real total = 0.;
    unsigned int num = 0;
    for (unsigned int i = 0; i < position.size(); ++i)
    {
        for (unsigned int j = 0; j < position.size(); ++j)
        {
            real n = (position[i] - position[j]).Norm(); // evtl periodic
            if (n >= r && n < r + rstep)
            {
                total += velocity[i]*velocity[j];
                num++;
            }
        }
    }

    return num == 0 ? 0. : total / num;
}

std::vector<real> correlationCisnerosXVelvic(const std::vector<Vector> &position, const std::vector<Vector> &velocity, int size, real rstep)
{
    std::vector<real> ret;
    real v2m = correlationCisnerosXV2Mvic(velocity);
    real vm2 = correlationCisnerosXVM2vic(velocity);
    for (real r = 0.; r < size * std::sqrt(2.); r += rstep)
    {
        ret.push_back((correlationCisnerosXVRMvic(position, velocity, r, rstep) - vm2) / (v2m - vm2));
    }

    return ret;
}

std::vector<real> correlationPoneXVelvic(const std::vector<Vector> &position, const std::vector<Vector> &velocity, int size, real rstep)
{
    std::vector<real> ret;

    for (real r = 0.; r < size * std::sqrt(2.); r += rstep)
    {
        assert(position.size() == velocity.size());
        real total = 0.;
        Vector total1(0.,0.), total2(0.,0.);
        unsigned int num = 0;
        for (unsigned int i = 0; i < position.size(); ++i)
        {
            for (unsigned int j = 0; j < position.size(); ++j)
            {
                real n = (position[i] - position[j]).Norm(); // evtl periodic
                if (n >= r && n < r + rstep)
                {
                    total += velocity[i]*velocity[j];
                    total1 += velocity[i];
                    total2 += velocity[j];
                    num++;
                }
            }
        }

        ret.push_back(num == 0 ? 0. : total / num - total1 * total2 / (num * num));
    }

    return ret;
}
