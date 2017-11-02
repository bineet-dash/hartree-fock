#include <cmath>
#include "common_globals.hpp"

using namespace std;

int no_of_pts=1000;
const int number_of_mesh=100;
double low_lim = -6.0;
double up_lim = 6.0;
double dx = (up_lim - low_lim)/double(no_of_pts);
const int no_of_sps = 3;
double tolerance = 0.0001;


double separation;

double V(double);

double inverted_gaussian(double x, double beta)
{
  double A0, alpha, gamma;
  A0 = 600.0;
  alpha = 15.0;
  gamma = 1.0;
  double v = -A0*exp(-(alpha*pow((x-beta),2))/(2*pow(gamma,2)));
  return v;
}

double V(double x)
{
  double beta;
  double v=0.0;

  for(int i=-2; i <=2 ; i++)
   {
     beta = i*separation;
     v += inverted_gaussian(x,beta);
   }
  return v;
}
