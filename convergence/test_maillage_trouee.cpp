//Compile with: make cas_troue
//
//test de remaillage autour d'untre arc de cercle

#include "mfem.hpp"
#include <fstream>
#include <iostream>

using namespace std;
using namespace mfem;
Vector param(4);
int main(int argc, char *argv[])
{
  // Parse command-line options.
  int order = 1;
  int rep = 1;
  const char *mesh_file_ = "test_remaillage.msh";
  OptionsParser args(argc, argv);
  args.AddOption(&order, "-o", "--order",
		 "Finite element order (polynomial degree).");
  args.AddOption(&rep, "-r", "--repet",
		 "Repetition de malliage.");
  args.AddOption(&mesh_file_, "-m", "--mesh",
 		 "mesh file to use.");
  args.Parse();
  args.PrintOptions(cout);

  Mesh *mesh = new Mesh(mesh_file_, 1, 1);
  const int dim = mesh->Dimension();

  for (int r = 0; r < rep; r++){
    cout<<endl;
    cout<<"Maillage "<<r<<endl<<"Coordonées sur le demi cercle:"<<endl;;
    //Save in Paraview format
    string  name_paraview = "Test_remaillage";
    string R_(to_string(r));
    name_paraview.append(R_);
    char* name_praview_ = &name_paraview[0];
    ParaViewDataCollection paraview_dc(name_praview_, mesh);
    paraview_dc.SetPrefixPath("ParaView");
    paraview_dc.SetLevelsOfDetail(order+1);
    paraview_dc.SetHighOrderOutput(true);
    paraview_dc.SaveMesh();
    paraview_dc.Save();

    //Coordonées sur le quart de cercle:
	const double* coord;
cout<<mesh->GetNE()<<endl;
for(int i=0;i<mesh->GetNE();i++){
      Element *bdr;
	bdr = mesh->GetElement(i);

	Array<int> v;
	bdr->GetVertices(v);
for (int j=0;j<v.Size();j++){
	coord = mesh->GetVertex(v[j]);
if(-6.3639<coord[0] && -6.3639>coord[1] && -8.1<coord[1]){
	cout<<"x: "<<coord[0]<<" y: "<<coord[1]<<endl;}
}
}
    mesh->UniformRefinement();
  }
  return 0;
}

void eqParabol(Vector &coord1, Vector &coord2, Vector &coord3,
	       Vector &coord4, Vector &param)
{
  param.SetSize(4);

  double b = coord1(0), y = coord1(1), a = coord1(0)*coord1(0);
  double s = coord2(0), f = coord2(1), q = coord2(0)*coord2(0);
  double x = coord3(0), v = coord3(1), w = coord3(0)*coord3(0);
  /*
    double r = coord1(0), t = coord1(1), z = coord1(0)*coord1(0), 
    a = coord1(0)*coord1(0)*coord1(0);
    double d = coord2(0), f = coord2(1), s = coord2(0)*coord2(0),
    q= coord2(0)*coord2(0)*coord2(0);
    double c = coord3(0), v = coord3(1), x= coord3(0)*coord3(0),
    w = coord3(0)*coord3(0)*coord3(0);
    double o = coord4(0), p = coord4(1), u = coord4(0)*coord4(0),
    y = coord4(0)*coord4(0)*coord4(0);
  */
  param(0) = (b*f-b*v+s*v-f*x-s*y+x*y)/(b*q-a*s-b*w+s*w+a*x-q*x);
  param(1) = (-a*f+a*v-q*v+f*w+q*y-w*y)/(b*q-a*s-b*w+s*w+a*x-q*x);
  param(2) = (b*q*v-a*s*v-b*f*w+a*f*x+s*w*y-q*x*y)/(b*q-a*s-b*w+s*w+a*x-q*x);
  /*
    param(0)=(-c*p*s+p*r*s+c*s*t-o*s*t+c*f*u-f*r*u-c*t*u+d*t*u+o*s*v-r*s*v-d*u*v+r*u*v-f*o*x+d*p*x+f*r*x-p*r*x-d*t*x+o*t*x-c*f*z+f*o*z+c*p*z-d*p*z+d*v*z-o*v*z)/(a*c*s-a*o*s-a*c*u+a*d*u+c*q*u-q*r*u+o*s*w-r*s*w-d*u*w+r*u*w-a*d*x+a*o*x-o*q*x+q*r*x-c*s*y+r*s*y+d*x*y-r*x*y-c*q*z+o*q*z+d*w*z-o*w*z+c*y*z-d*y*z);
    param(1)=(a*c*f-a*f*o-a*c*p+a*d*p+c*p*q-p*q*r-c*q*t+o*q*t-a*d*v+a*o*v-o*q*v+q*r*v+f*o*w-d*p*w-f*r*w+p*r*w+d*t*w-o*t*w-c*f*y+f*r*y+c*t*y-d*t*y+d*v*y-r*v*y)/(a*c*s-a*o*s-a*c*u+a*d*u+c*q*u-q*r*u+o*s*w-r*s*w-d*u*w+r*u*w-a*d*x+a*o*x-o*q*x+q*r*x-c*s*y+r*s*y+d*x*y-r*x*y-c*q*z+o*q*z+d*w*z-o*w*z+c*y*z-d*y*z);
    param(2)=(-a*p*s+a*f*u-q*t*u+a*s*v-a*u*v+q*u*v+p*s*w-s*t*w-f*u*w+t*u*w-a*f*x+a*p*x-p*q*x+q*t*x+s*t*y-s*v*y+f*x*y-t*x*y+p*q*z-q*v*z+f*w*z-p*w*z-f*y*z+v*y*z)/(a*c*s-a*o*s-a*c*u+a*d*u+c*q*u-q*r*u+o*s*w-r*s*w-d*u*w+r*u*w-a*d*x+a*o*x-o*q*x+q*r*x-c*s*y+r*s*y+d*x*y-r*x*y-c*q*z+o*q*z+d*w*z-o*w*z+c*y*z-d*y*z);
    param(3)=(a*c*p*s-a*c*f*u+c*q*t*u-a*o*s*v+a*d*u*v-q*r*u*v-p*r*s*w+o*s*t*w+f*r*u*w-d*t*u*w+a*f*o*x-a*d*p*x+p*q*r*x-o*q*t*x-c*s*t*y+r*s*v*y-f*r*x*y+d*t*x*y-c*p*q*z+o*q*v*z-f*o*w*z+d*p*w*z+c*f*y*z-d*v*y*z)/(a*c*s-a*o*s-a*c*u+a*d*u+c*q*u-q*r*u+o*s*w-r*s*w-d*u*w+r*u*w-a*d*x+a*o*x-o*q*x+q*r*x-c*s*y+r*s*y+d*x*y-r*x*y-c*q*z+o*q*z+d*w*z-o*w*z+c*y*z-d*y*z);
  */
}
