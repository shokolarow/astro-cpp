#include "tools.h"

using namespace EARTH;

datetime::datetime(double yr, double mn, double d, double h, double m, double s) :
        year(yr), month(mn), day(d), hour(h), minute(m), second(s)
{}

std::ostream& operator <<(std::ostream &o, const datetime& date){
    o << date.day << "/" << date.month << "/" << date.year << " " << date.hour
    << ":" << date.minute << ":" << date.second << std::endl;

    return o;
}

Eigen::RowVectorXd koe_to_csv(double a, double e, double i, double om, double w, double nu){
    i *= d2r;
    w *= d2r;
    om *= d2r;
    nu *= d2r;

    double r = a * (1 - e*e) / (1 + e*std::cos(nu));
    double h = std::sqrt(mu * a * (1 - e*e));

    double X = r * (std::cos(om) * std::cos(w+nu) - std::sin(om) * std::sin(w+nu) * std::cos(i));
    double Y = r * (std::sin(om) * std::cos(w+nu) + std::cos(om) * std::sin(w+nu) * std::cos(i));
    double Z = r * std::sin(i) * std::sin(w+nu);

    double p = a * (1 - e*e);

    double Vx = ((X*h*e) / (r*p)) * std::sin(nu) - (h/r) * (std::cos(om) * std::sin(w+nu) \
                                                            + std::sin(om) * std::cos(w+nu) * std::cos(i));
    double Vy = ((Y*h*e) / (r*p)) * std::sin(nu) - (h/r) * (std::sin(om) * std::sin(w+nu) \
                                                            - std::cos(om) * std::cos(w+nu) * std::cos(i));
    double Vz = ((Z*h*e) / (r*p)) * std::sin(nu) + (h/r) * std::sin(i) * std::cos(w+nu);

    Eigen::RowVectorXd csv(6);
    csv << X,Y,Z,Vx,Vy,Vz;
    return csv;
}

double MA2v(double MA, double e){
    double tolerance = 1E-3;
    double EA = MA;
    while ((EA - e*std::sin(EA) - MA) > tolerance){
        EA = (EA + (MA - EA + e*std::sin(e))/(1 - e*std::cos(e)));
        EA = std::fmod(EA, 2*M_PI);
    }

    return 2 * atan(std::sqrt((1+e)/(1-e)) * std::tan(EA/2));
}

std::pair<Eigen::Matrix<double,1,6>, datetime> parse_tle(const std::string& s1, const std::string& s2){
    std::istringstream it1(s1);
    std::vector<std::string> tle1{
        std::istream_iterator<std::string>(it1), {}
    };

    std::string yd = tle1[3];
    int year = std::stoi("20"+yd.substr(0,2));
    int day = std::stoi(yd.substr(2, 5)) - 1;

    time_t loctime;
    struct tm timeinfo, *loctimeinfo;

    bzero(&timeinfo, sizeof(struct tm));
    timeinfo.tm_isdst = -1;
    timeinfo.tm_mon = 0;
    timeinfo.tm_mday = day;
    timeinfo.tm_year = year - 1900;
    loctime = mktime(&timeinfo);
    loctimeinfo = localtime(&loctime);

    int month = loctimeinfo->tm_mon + 1;
    day = loctimeinfo->tm_mday + 1;

    double fr = std::stod(yd.substr(5, yd.length()));
    double hour = fr * 24;
    double minute = (hour - static_cast<int>(hour)) * 60;
    double second = (minute - static_cast<int>(minute)) * 60;

    datetime date(year, month, day, static_cast<int>(hour), static_cast<int>(minute),
                  static_cast<int>(second));

    std::istringstream it2(s2);
    std::vector<std::string> tle2{
            std::istream_iterator<std::string>(it2), {}
    };

    double i = std::stod(tle2[2]);
    double om = std::stod(tle2[3]);
    double e = std::stod("0."+tle2[4]);
    double w = std::stod(tle2[5]);
    double MA = std::stod(tle2[6]) * d2r;
    double v = MA2v(MA, e) * r2d;
    double n = std::stod(tle2[7]) * (2 * M_PI / 24 / 60 / 60);
    double a = std::cbrt(mu/n/n);

    Eigen::Matrix<double, 1, 6> koe;
    koe << a, e, i, om, w, v;

    return {koe, date};
}