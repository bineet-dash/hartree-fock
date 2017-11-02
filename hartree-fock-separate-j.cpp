/*
 * hartree_fock_separate_j.cpp
 *
 * Copyright 2017 Bineet Dash <bineet@bineet-ubuntu>
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston,
 * MA 02110-1301, USA.
 *
 *
 */
#include <cstdlib>
#include <iomanip>
#include <ctime>
#include <cstring>
#include <iostream>
#include <cmath>
#include <fstream>
#include <chrono>
#include <string>
#include <Eigen/Dense>
#include <lapacke.h>

using namespace std;
using namespace std::chrono;
using namespace Eigen;

typedef std::complex <double> cd;

#include "hermite_polynomial.hpp"

double Sqr(cd x){return (x*conj(x)).real();}

bool diagonalize(MatrixXcd Ac, VectorXcd& lambdac, MatrixXcd& vc)
{
  int N;
  if(Ac.cols()==Ac.rows())  N = Ac.cols(); else return false;

  MatrixXd A = Ac.real();
  lambdac.resize(N);
  vc.resize(N,N);
  VectorXd lambda = lambdac.real();

  int LDA = A.outerStride();
  int INFO = 0;
  char Uchar = 'U';
  char Vchar = 'V';

  int LWORK = 5*(2*LDA*LDA+6*LDA+1);
  int LIWORK = 5*(3+5*LDA);

  VectorXd WORK(LWORK);
  VectorXi IWORK(IWORK);

  dsyevd_(&Vchar, &Uchar, &N, A.data(), &LDA, lambda.data(),  WORK.data(), &LWORK, IWORK.data(), &LIWORK, &INFO);
  vc.real() = A;
  lambdac.real() = lambda;

  return INFO==0;
}

bool cEP(MatrixXcd A, VectorXcd& lambda, MatrixXcd& v)
{
  int N = A.cols();
  if (A.rows()!=N)  return false;
  v.resize(N,N);
  lambda.resize(N);

  int LDA = A.outerStride();
  int LDV = v.outerStride();
  int INFO = 0;
  cd* w = const_cast<cd*>(lambda.data());
  char Nchar = 'N';
  char Vchar = 'V';
  int LWORK = int(A.size())*4;
  VectorXcd WORK(LWORK);
  VectorXd RWORK(2*LDA);

  zgeev_(&Nchar, &Vchar, &N, reinterpret_cast <__complex__ double*> (A.data()), &LDA, reinterpret_cast <__complex__ double*> (w), 0, &LDV, reinterpret_cast <__complex__ double*> (v.data()), &LDV,  reinterpret_cast <__complex__ double*> (WORK.data()), &LWORK, RWORK.data(), &INFO);

  for(int i=0; i<N; i++)
 	 v.col(i)=v.col(i)/v.col(i).unaryExpr(&Sqr).sum();

  return INFO==0;
}


int no_of_pts=1000;
const int number_of_mesh=100;
const double a = 1.0;
double low_lim = -6.0;
double up_lim = 6.0;
double dx = (up_lim - low_lim)/double(no_of_pts);
const double omega=1.0;
double alpha = 1/sqrt(omega);
const int no_of_sps = 3;
VectorXd point(no_of_pts+1);
MatrixXcd states(point.size(),no_of_sps);

void show_time(milliseconds begin_ms, milliseconds end_ms);
double fac(int n) {double prod=1.0; for(int i=n; i>1;i--) prod*=n; return prod;}
double hermite(int n, double y)
{
    double x_vec[1];
    x_vec[0]=y;
    double* fx2_vec = h_polynomial_value ( 1, n, x_vec );
    double fx2 = fx2_vec[n];
    return fx2;
}

double V(double x) {return 0.5*pow(omega*x,2);}

double psi(int n, double x)
{
  n++;
  if(abs(x)>1) return 0;
  else return (n%2==0)? sqrt(2/a)*sin(n*M_PI*x/(2*a)):sqrt(2/a)*cos(n*M_PI*x/(2*a));
}

// double psi(int n, double x)
// {
//   long double y = alpha*x;
//   long double result = 1/sqrt(pow(2,n)*fac(n))*sqrt(alpha)/pow(M_PI,0.25)*exp(-y*y/2)*hermite(n,y);
//   return result;
// }

bool compare(const pair<double, VectorXcd>&i, const pair<double, VectorXcd>&j) {return i.first < j.first;}

double rho_H(double r_prime)
{
  int n_prime = int((r_prime - low_lim)/dx);
  double rho=0.0;
  for(int i=0; i< no_of_sps; i++)
  {
    rho += (conj(states(n_prime,i))*states(n_prime,i)).real();
  }          //phi_i*(r')phi_i(r')
  return rho;
}

double rho_HF(int n , int n_prime, int j)
{
  double num=0.0; double denom=0.0;

  for(int k=0; k< no_of_sps; k++)
  {
    num += (conj(states(n_prime,k))*states(n_prime,j)*conj(states(n,j))*states(n,k)).real();
           //phi_k*(r')phi_j(r')phi_j*(r)phi_k(r)
  }
  denom = (states(n,j)*conj(states(n,j))).real();
  if(denom != 0) return num/denom; else return 0.0;
}

double rho_HF(double r, double r_prime, int j)
{
  double num=0.0; double denom=0.0;
  int n = (r - low_lim)/dx;
  int n_prime = (r_prime - low_lim)/dx;

  for(int k=0; k< no_of_sps; k++)
  {
    num += (conj(states(n_prime,k))*states(n_prime,j)*conj(states(n,j))*states(n,k)).real();
           //phi_k*(r')phi_j(r')phi_j*(r)phi_k(r)
  }
  denom = (states(n,j)*conj(states(n,j))).real();
  if(denom != 0) return num/denom; else return 0.0;
}

double integrand(double r, double r_prime, int j)
{
  return (rho_H(r_prime)-rho_HF(r,r_prime,j))/(abs(r - r_prime)+1/(2.0*double(number_of_mesh)));
}

double integrate_rho(double r, int j, double (*func_x)(double, double, int)) //func_x=integrand
{
  double trapez_sum;
  double fa, fb,x, step;
  step=(up_lim - low_lim)/((double) number_of_mesh);
  fa=(*func_x)(r,low_lim,j);
  fb=(*func_x)(r,up_lim,j);
  trapez_sum=0.;
  for (int k=1; k < number_of_mesh; k++)
  {
    x=j*step+low_lim;
    trapez_sum+=(*func_x)(r,x,j);
  }
  trapez_sum=(trapez_sum+fb+fa)*step;
  return trapez_sum;
}

int main()
{

  double tolerance;
  cout << "Enter tolerance: "; cin >> tolerance;

  for(int i=0; i<= no_of_pts; i++) {point(i)=low_lim+i*dx;}
  for(int i=0; i< point.size(); i++)
  {
    for(int j=0; j<no_of_sps; j++)
      states(i,j) = psi(j,point(i));
  }

  for(int j=0; j<no_of_sps; j++)
      states.col(j) = states.col(j)/sqrt(states.col(j).unaryExpr(&Sqr).sum());

  MatrixXcd H = MatrixXcd::Zero(point.size(),point.size());
  for(int i=0; i<point.size(); i++)
  {
      int j = (i==point.size()-1)? 0 : i+1;
      H(i,j)= -1/(2*dx*dx);
      H(j,i)= -1/(2*dx*dx);
  }

  VectorXcd v; MatrixXcd eigenvectors; VectorXd eigenvalues;
  int master_loop = 1;
  VectorXd oldeival= VectorXd::Zero(no_of_sps);
  VectorXd neweival= VectorXd::Zero(no_of_sps);
  vector < pair<double,VectorXcd> > selected_spectrum;

   ofstream fout("data/initialstate.txt");
   for(int i=0; i<point.size(); i++)
   {
     fout << point(i) << " ";
     for(int j=0; j<no_of_sps; j++) fout << (states(i,j)).real() << " ";
     fout << endl;
   }
   fout.close();

   cout.precision(8);

  //  int n, n_prime;
  //  while(1==1)
  // { cin >> n >> n_prime;
  //  for(int j=0; j<no_of_sps; j++) cout << rho_HF(n,n_prime,j) << endl;}
  //  exit(1);

  for(; ; )
  {
    cout << "Loop-" << master_loop << "\n============================\n";

    for(int j=0; j<no_of_sps; j++)
      {
        for(int i=0; i<point.size(); i++)
          H(i,i) = 1/(dx*dx)+ V(point(i))+ integrate_rho(point(i),j,&integrand);

        milliseconds begin_ms = duration_cast< milliseconds >(system_clock::now().time_since_epoch());
        diagonalize(H,v,eigenvectors); eigenvalues = v.real();
        milliseconds end_ms = duration_cast< milliseconds >(system_clock::now().time_since_epoch());
        show_time(begin_ms,end_ms);

        vector < pair<double,VectorXcd> > eigenspectrum;
        for(int k=0; k<point.size(); k++) eigenspectrum.push_back(make_pair(eigenvalues(k),eigenvectors.col(k)));

        sort(eigenspectrum.begin(),eigenspectrum.end(),compare);
        selected_spectrum.push_back(eigenspectrum.at(j));
        eigenspectrum.clear();
      }


    for(int i=0; i<no_of_sps; i++) states.col(i)= selected_spectrum.at(i).second;
    oldeival = neweival;
    for(int i=0; i<neweival.size(); i++) neweival(i) = selected_spectrum.at(i).first;
    cout << "Eigenvalues are: " << neweival.transpose() << endl << endl;

    string filename = "data/loop"+to_string(master_loop)+".txt";
    fout.open(filename);
    for(int i=0; i<point.size(); i++)
    {
      fout << point(i) << " ";
      for(int j=0; j<no_of_sps; j++) fout << (states(i,j)).real() << " ";
      fout << endl;
    }
    fout.close();

    // cout << "Normalization check: " << endl;
    // for(int i=0; i< states.cols(); i++)
  	// cout << "Column-" << i << " Normalization=" << sqrt(states.col(i).unaryExpr(&Sqr).sum()) << endl;

    double max_deviation = (neweival - oldeival).cwiseAbs().maxCoeff();
    if(max_deviation < tolerance) break;

    selected_spectrum.clear();
    master_loop++; cout << endl;
  }

  	// cout << "Final Normalization: " << endl;
    // for(int i=0; i< states.cols(); i++)
  	// cout << "Column-" << i << " Normalization=" << sqrt(states.col(i).unaryExpr(&Sqr).sum()) << endl;
}

void show_time(milliseconds begin_ms, milliseconds end_ms)
{
   long int t = (end_ms.count()-begin_ms.count())/1000;
    if (t<=60)
    { cout << "Diagonalization took " << t << " seconds." << endl; }
    else if (t>60 && t<3600)
    {
      int minutes=int(t/60); int seconds= t%minutes;
      cout << "Diagonalization took " << minutes << " minute and " << seconds << " seconds." << endl;
    }
    else if(t>=3600)
    {
      int hrs= int(t/3600); int minutes=int((t-3600*hrs)/60); int seconds= int(t-3600*hrs-60*minutes);
      cout << "Diagonalization took " << hrs << " hour, " << minutes << " minutes and " << seconds << " seconds. ";
    }
    else
    {cout << "Diagonalization took " << t << "time. Wrong t received.\n"; }
}
