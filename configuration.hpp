#include <cmath>
#include "common_globals.hpp"

using namespace std;

/*Declaration of parameters fixed during runtime. */

/*Essential parameters.*/
/*Don't add your own parameter in this block */

int no_of_pts=1000;       //Number of grid points
double low_lim = -6.0;    //Lower limit of x. Grid starts here
double up_lim = 6.0;      //Upper Limit of x. Grid ends here.
double dx = (up_lim - low_lim)/double(no_of_pts); //Step size
const int no_of_sps = 3;  //Number of electrons (i.e. no of trial basis statesa)
double tolerance = 1e-4;  //Tolerance for the convergence of hartree_fock iterations

double V(double);
/*Essential parameters end*/

/*Custom parameters*/
/*Define your own parameters here*/
const double omega=1.0;

/*Custom parameters end*/

/*Definition of Potential*/
double V(double x)
{
  return 0.5*pow(omega*x,2);
}
