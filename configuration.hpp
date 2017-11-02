#include <cmath>
#include "common_globals.hpp"

using namespace std;

int no_of_pts=1000;
const int number_of_mesh=100;
double low_lim = -6.0;
double up_lim = 6.0;
double dx = (up_lim - low_lim)/double(no_of_pts);
const int no_of_sps = 3;

const double a = 1.0;
const double omega=1.0;
double alpha = 1/sqrt(omega);

double V(double);

double V(double x)
{
  return 0.5*pow(omega*x,2);
}
