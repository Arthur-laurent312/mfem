//Compile with: make cas_troue
//
//test de remaillage autour d'untre arc de cercle

#include "mfem.hpp"
#include <fstream>
#include <iostream>

using namespace std;
using namespace mfem;

void Transform(DenseMatrix &, Vector &, Vector &);
double L1(const Vector &);
double L2(const Vector &);
double L3(const Vector &);
double N1(const Vector &);
double N2(const Vector &);
double N3(const Vector &);
double N4(const Vector &);
double N5(const Vector &);
double N6(const Vector &);

int main(int argc, char *argv[])
{
  double step = 0.001;	//pas de calcul pour le tracé
  DenseMatrix coord(2,3);
  coord(0,0) = -9, coord(1,0) = 0;
  coord(0,1) = 0, coord(1,1) = -9;
  //coord(0,1) = -6.36396, coord(1,1) = -6.36396;
  coord(0,2) = -6.36396, coord(1,2) = -6.36396;
  //coord(1,2) = -8.31492 , coord(0,2) = -3.44415;

  //coord1(1,2) = -1.78762 ,coord1(0,2) = -8.81569;

  string const parabole1("Trace_parabole.txt");
  ofstream parabole_flux1(parabole1.c_str());
  string const parabole2("Trace_parabole2.txt");
  ofstream parabole_flux2(parabole2.c_str());

  if (!parabole_flux2.is_open()) {
    cout << "Problem in openning file" << endl;
    exit(0);
  }
  else{
    cout <<"File openning"<<endl;
    Vector varx(2), vary(2);
	Vector Tx(2), Ty(2);
    for(double i=0; i<=1 ; i+=step)
      {
	varx(0) = i, varx(1) = 0.;
	vary(0) = i, vary(1) = -(i-1);
	Transform(coord, varx, Tx);
	Transform(coord, vary, Ty);
	parabole_flux1<<Tx(0)<<" "<<Tx(1)<<" "<<Ty(0)<<" "<<Ty(1)<<endl;
      }
	//for(double j=-6.36396; j<=0; j+=step)
	for(double j=-9; j<=0; j+=step)
	{
	 parabole_flux2<<j<<" "<<-pow(9*9-j*j,0.5)<<endl;
	} 
 }
  return 0;
}

void Transform(DenseMatrix &coord, Vector &variable, Vector &T)
{
  T.SetSize(2);
  T(0) = N1(variable)*coord(0,0) + N2(variable)*coord(0,1) +
    N3(variable)*coord(0,2);
  T(1) = N1(variable)*coord(1,0) + N2(variable)*coord(1,1) +
    N3(variable)*coord(1,2);
}

double L1(const Vector &var)
{
  return 1 - var(0) - var(1);
}
double L2(const Vector &var)
{
  return var(0);
}
double L3(const Vector &var)
{
  return var(1);
}

double N1(const Vector &var)
{
  double L=L1(var);
  return L*(2*L-1);
  //return (1-var(0))*(1-var(1));
}
double N2(const Vector &var)
{
  double L=L2(var);
  return L*(2*L-1);
  //return (1+var(0))*(1-var(1));
}
double N3(const Vector &var)
{
  double L1_ = L1(var);
  double L2_ = L2(var);
  return 4*L1_*L2_;
  //return (1+var(0))*(1+var(1));
}
double N4(const Vector &var)
{
  return (1-var(0))*(1+var(1));
}
double N5(const Vector &var)
{
  double L1_ = L1(var);
  double L3_ = L3(var);
  return 4*L1_*L3_;
}
double N6(const Vector &var)
{
  double L1_ = L1(var);
  double L2_ = L2(var);
  return 4*L1_*L2_;
}
