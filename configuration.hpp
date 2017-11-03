#include <cmath>
#include "common_globals.hpp"

using namespace std;

int no_of_pts=1000;
double low_lim = -6.0;
double up_lim = 6.0;
double dx = (up_lim - low_lim)/double(no_of_pts);
const int no_of_sps = 3;
double tolerance = 0.0001;

const double omega=1.0;

double V(double);

double V(double x)
{
  return 0.5*pow(omega*x,2);
}
