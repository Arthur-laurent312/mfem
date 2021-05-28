//Compile with: make cas_troue
//
//Cas test elsatique linéaire d'une plaque trouée. Un 
//bord est en traction avec la contrainte Sigmainf imposée
//Par symétrie on peut considérer seulement un quart de la plaque
//
//Ref: https://perso.ensta-paris.fr/~mbonnet/codes_donnees.pdf

#include "mfem.hpp"
#include <fstream>
#include <iostream>

using namespace std;
using namespace mfem;

static constexpr double R2 = 0.25*0.25;	//Rayon troue au carré
static constexpr double Sigmainf = 10.0;	//contrainte à l'infinie
static constexpr double E = 1000.;	//Module de Young
static constexpr double nu = 0.25;	//Coef de poisson
static constexpr double lambda = E*nu/((1.+nu)*(1.-2.*nu)); //Coef de lamé
static constexpr double mu = E/(2.*(1.+nu));

//Solution exacte
void conversion(const double, const double , double &, double &);
void Stress_exacte(const Vector &, Vector &);
void Strain_(DenseMatrix &, Vector &);

//Chargement
double F(const Vector &);

//Erreur en norme energy
double ComputeEnergyNorm(GridFunction &,
			 Coefficient &, Coefficient &,
			 VectorFunctionCoefficient &);

void Elasticy_mat(ElementTransformation &,const IntegrationPoint &, int,
		  Coefficient &, Coefficient &, DenseMatrix &);
void ComputeStrain(ElementTransformation &,const IntegrationPoint &,
		   GridFunction &, int,  Vector &);

void SaveStrain(GridFunction &,Coefficient &, Coefficient &,
		Mesh &, QuadratureFunction );
int main(int argc, char *argv[])
{
  DenseMatrix slope_ener;
  double err_tmp_ener = 0;
  double h_tmp = 0.;
  int iter = 0;
  // Parse command-line options.
  bool static_cond = false;
  int order = 1;
  bool iterative = true;
  int rep = 1;
  OptionsParser args(argc, argv);
  args.AddOption(&order, "-o", "--order",
		 "Finite element order (polynomial degree).");
  args.AddOption(&rep, "-r", "--repet",
		 "Repetition de malliage.");
  args.AddOption(&static_cond, "-sc", "--static-condensation", "-no-sc",
		 "--no-static-condensation", "Enable static condensation.");
  args.AddOption(&iterative, "-it", "--iterative", "-di",
		 "--direct", "Enable direct or iterative solver.");
  args.Parse();
  if (!args.Good())
    {
      args.PrintUsage(cout);
      return 1;
    }
  args.PrintOptions(cout);
  slope_ener.SetSize(rep-1,3);
  for (int r = 0; r < rep; r++){
    string  mesh_file = "hole_mesh/quarter_phole";
    string buf(mesh_file);
    string X_(to_string(r));
    buf.append(X_);
    string buff(buf);
    string mesh_file2 = ".msh";
    buff.append(mesh_file2);

    char* mesh_file_ = &buff[0];
    cout<<"Nom du maillage: "<<mesh_file_<<endl;
    Mesh *mesh = new Mesh(mesh_file_, 1, 1);
    // 2. Read the mesh from the given mesh file. We can handle triangular,
    //    quadrilateral, tetrahedral or hexahedral elements with the same code.
    int dim = mesh->Dimension();
  
    // 5. Define a finite element space on the mesh. Here we use vector finite
    //    elements, i.e. dim copies of a scalar finite element space. The vector
    //    dimension is specified by the last argument of the FiniteElementSpace
    //    constructor. For NURBS meshes, we use the (degree elevated) NURBS space
    //    associated with the mesh nodes.
    FiniteElementCollection *fec;
    FiniteElementSpace *fespace;
    fec = new H1_FECollection(order, dim);
    fespace = new FiniteElementSpace(mesh, fec, dim);
    cout << "Numbers of elements: " << mesh->GetNE() <<endl;
    cout << "Number of finite element unknowns: " << fespace->GetTrueVSize()
	 << endl << "Assembling: " << flush;
  
    // 6. Determine the list of true (i.e. conforming) essential boundary dofs.
    //    In this example, the boundary conditions are defined by marking only
    //    boundary attribute 1 from the mesh as essential and converting it to a
    //    list of true dofs.
  
    // List of True DoFs : Define (here) Dirichlet conditions
    Array<int> ess_tdof_list, tmp_tdof, ess_bdr(mesh->bdr_attributes.Max());
    cout << "bdr_attributes.Max()" << std::endl;
    ess_bdr = 0;
    // ess_bdr[0] refers to dof on the "xz0" plane
    ess_bdr[1] = 1;
    // grab all dof on "xz0" and put them into tmp_dof.
    // last parameter is 1 to set "y" direction for dirichlet condition.
    fespace->GetEssentialTrueDofs(ess_bdr, tmp_tdof, 1); 
    // ess_tof_list accumulates all needed dof
    ess_tdof_list.Append(tmp_tdof);

    ess_bdr = 0;
    // ess_bdr[1] refer to dof on the "xy0" plane
    ess_bdr[3] = 1;
    // last parameter is 2 to set "z" direction for dirichlet condition.
    fespace->GetEssentialTrueDofs(ess_bdr, tmp_tdof, 0);
    // ess_tof_list accumulates all needed dof
    ess_tdof_list.Append(tmp_tdof);
  
    // 7. Set up the linear form b(.) which corresponds to the right-hand side of
    //    the FEM linear system. In this case, b_i equals the boundary integral
    //    of f*phi_i where f represents a "pull down" force on the Neumann part
    //    of the boundary and phi_i are the basis functions in the finite element
    //    fespace. The force is defined by the VectorArrayCoefficient object f,
    //    which is a vector of Coefficient objects.
    //    The fact that f is non-zero
    //    on boundary attribute 2 is indicated by the use of piece-wise constants
    //    coefficient for its last component.
  
    VectorArrayCoefficient f(dim);
    for (int i = 0; i < dim-1; i++)
      {
	f.Set(i, new ConstantCoefficient(0.0));
      }
    f.Set(dim-1, new FunctionCoefficient(F));
  
    LinearForm *b = new LinearForm(fespace);
    b->AddBoundaryIntegrator(new VectorBoundaryLFIntegrator(f));
    cout << "r.h.s. ... " << flush;

    // 8. Define the solution vector x as a finite element grid function
    //    corresponding to fespace. Initialize x with initial guess of zero,
    //    which satisfies the boundary conditions.
    GridFunction x(fespace);
    x = 0.0;
  
    // 9. Set up the bilinear form a(.,.) on the finite element space
    //    corresponding to the linear elasticity integrator with piece-wise
    //    constants coefficient lambda and mu.
  
    ConstantCoefficient mu_func(mu);
    ConstantCoefficient lambda_func(lambda);
    BilinearForm *a = new BilinearForm(fespace);
    BilinearFormIntegrator *integ = new ElasticityIntegrator(lambda_func, mu_func);
    a->AddDomainIntegrator(integ);
  
    // 10. Assemble the bilinear form and the corresponding linear system,
    //     applying any necessary transformations such as: eliminating boundary
    //     conditions, applying conforming constraints for non-conforming AMR,
    //     static condensation, etc.
  
    cout << "matrix ... " << flush;
  
    a->Assemble();
    b->Assemble();
  
    SparseMatrix A;
    Vector B, X;
  
    a->FormLinearSystem(ess_tdof_list, x, *b, A, X, B);  
    cout << "done." << endl;
  
    cout << "Size of linear system: " << A.Height() << endl;
  
    if(iterative){
      GSSmoother M(A);
      PCG(A, M, B, X, 2, 10000, 1e-20, 0.0);
    } else {
#ifdef MFEM_USE_SUITESPARSE	  
      // 11. If MFEM was compiled with SuiteSparse, use UMFPACK to solve the system.
      UMFPackSolver umf_solver;
      umf_solver.Control[UMFPACK_ORDERING] = UMFPACK_ORDERING_METIS;
      umf_solver.SetOperator(A);
      umf_solver.Mult(B, X);
#else
      cout<<"Direct solver not implemented" << endl;
      exit(0);
#endif
    }
    // 12. Recover the solution as a finite element grid function.
    a->RecoverFEMSolution(X, *b, x);

    // Compute norms of error
    int tdim = dim*(dim+1)/2; // num. entries in a symmetric tensor
    VectorFunctionCoefficient Stress_exacte_coef(tdim,Stress_exacte);
    double ener_error = ComputeEnergyNorm(x, lambda_func, mu_func,
					  Stress_exacte_coef);
    cout << "Energy norm of error: " << ener_error << endl;
    cout << endl;
    double h = mesh->GetElementSize(0);
    //Compute the slope
    slope_ener(iter,0) = log(err_tmp_ener/ener_error) / log(h_tmp/h);
    slope_ener(iter,1) = h;
    slope_ener(iter,2) = ener_error;
    err_tmp_ener = ener_error;
    h_tmp = h;
    iter++;
    //Save in Praview format
    QuadratureFunction Strain;
    SaveStrain(x, lambda_func, mu_func, *mesh, Strain);
    if (!mesh->NURBSext)
      {
        mesh->SetNodalFESpace(fespace);
      }
    GridFunction *nodes = mesh->GetNodes();
    *nodes += x;
    x *= -1;
    ParaViewDataCollection paraview_dc("Troue", mesh);
    paraview_dc.SetPrefixPath("ParaView");
    paraview_dc.SetLevelsOfDetail(order+1);
    paraview_dc.SetCycle(0);
    paraview_dc.SetDataFormat(VTKFormat::BINARY);
    paraview_dc.SetHighOrderOutput(true);
    paraview_dc.SetTime(0.0); // set the time
    paraview_dc.RegisterField("numerical_solution",&x);  
    //paraview_dc.RegisterQField("Contrainte",&Strain);
    paraview_dc.Save();	
    /*
      delete a;
      delete b;
      if (fec) {
      delete fespace;
      delete fec;
      }*/
  }//Fin de boucle maillage
  //Affichage des normes et pentes.
  cout<<endl;
  cout<<"Erreur en norme:"<<endl;
  for (int i=1; i<iter; i++){
    cout << " Energie: " << slope_ener(i,2)<<" Taille de maille= "
	 <<slope_ener(i,1)<<endl;}

  cout<<endl;
  cout<<"Pente de convergence:"<<endl;
  for (int i=1; i<iter; i++){
    cout << " Energie: " << slope_ener(i,0)<<" Taille de maille= "<<
      slope_ener(i,1)<<endl;}
  return 0;
}
void conversion(const double x, const double y , double &r, double &theta)
{
  double x1=x-0.5;
  double y1=y-0.5;
  r = sqrt(x1*x1 + y1*y1);
  //theta = atan2(y1,x1);
  if ( x > 0 && y > 0)
    theta = atan(y/x);
  if ( x > 0 && y == 0)
    theta = 2.0*M_PI;
  if ( x > 0 && y < 0)
    theta = atan(y/x) + 2*M_PI;
  if ( x < 0)
    theta = atan(y/x) + M_PI;
  if ( x == 0 && y > 0)
    theta = M_PI/2.0;
  if ( x == 0 && y < 0)
    theta = 3*M_PI/2.0;
  if ( x == 0 && y == 0)
    theta = 0.0;

}
//===================== Stress exacte =====================
void Stress_exacte(const Vector &x, Vector &stress)
{
  double r, theta;
  conversion(x(0), x(1), r, theta);
  double r2 = r*r;
  stress(0) = Sigmainf*0.5*((1-R2/r2) + (1+3*R2*R2/(r2*r2)-4*R2/r2)*cos(2*theta));
  stress(1) = Sigmainf*0.5*((1+R2/r2) - (1+3*R2*R2/(r2*r2))*cos(2*theta));
  stress(2) = -Sigmainf*0.5*(1-3*R2*R2/(r2*r2)+2*R2/r2)*sin(2*theta);
}
//===================== Strain exacte =====================
void Strain_(DenseMatrix &grad, Vector &strain)
{
  int dim = grad.Size();
  strain.SetSize(dim*(dim+1)/2);
  strain(0)=grad(0,0);
  strain(1)=grad(1,1);
  if(dim==2){
    strain(2)=0.5*(grad(1,0)+grad(0,1));
  }
  else if(dim==3){
    strain(2)=grad(2,2);
    strain(3)=0.5*(grad(1,0)+grad(0,1));
    strain(4)=0.5*(grad(2,0)+grad(0,2));
    strain(5)=0.5*(grad(2,1)+grad(1,2));
  }
  else{
    cout<<"Dimention not suported"<<endl;}
}
//===================== Charge sur le bord droit =====================
double F(const Vector &x)
{
  double force;
  if(x(0) == 0){
    force = Sigmainf;
  } else {
    force = 0.;
  }
  return force;
}
//===================== Strain dans QuadratureFunction =====================
void SaveStrain(GridFunction &x,Coefficient &lambdah, Coefficient &muh,
		Mesh &mesh, QuadratureFunction Strain){
  FiniteElementSpace *fes = x.FESpace();
  int dim = fes->GetMesh()->SpaceDimension();
  int tdim = dim*(dim+1)/2; // num. entries in a symmetric tensor
  ElementTransformation *Trans;
  const FiniteElement *fe = fes->GetFE(0);
  const int order = 2*fe->GetOrder() + 3;   //<----------
  const IntegrationRule *ir = &(IntRules.Get(fe->GetGeomType(), order));
  int Ni = ir->GetNPoints();
  int Ne = fes->GetNE();
  Vector data(tdim*Ni*Ne);
  Vector stressh(tdim), strainh(tdim);
  DenseMatrix C;
  for (int i = 0; i < Ne ; i++)
    {
      const FiniteElement *fe = fes->GetFE(i);
      const int order = 2*fe->GetOrder() + 3;   //<----------
      const IntegrationRule *ir = &(IntRules.Get(fe->GetGeomType(), order));
      Trans = fes->GetElementTransformation(i);
      const int dof = fe->GetDof();
      for (int k = 0; k < Ni; k++)
	{
	  const IntegrationPoint &ip = ir->IntPoint(k);
	  Trans->SetIntPoint(&ip);
	  double w = Trans->Weight() * ip.weight;
	  ComputeStrain(*Trans, ip, x, i, strainh);
	  Elasticy_mat(*Trans,ip,dim,lambdah,muh,C);
	  C.Mult(strainh,stressh);	
	  for(int j = 0; j < tdim; j++)
	    data(k+i+j)=stressh(j);
	}
    }
  QuadratureSpace StrainSpace(&mesh, order);
  Strain.SetSpace(&StrainSpace, data, tdim);
}

//==============Erreur en Norme energie ===================
double ComputeEnergyNorm(GridFunction &x, Coefficient &lambdah, Coefficient &muh,
			 VectorFunctionCoefficient &Stress_exacte_coef)
{
  FiniteElementSpace *fes = x.FESpace();
  int dim = fes->GetMesh()->SpaceDimension();
  const int tdim = dim*(dim+1)/2; // num. entries in a symmetric tensor
  ElementTransformation *Trans;

  double energy = 0.0;
  for (int i = 0; i < fes->GetNE() ; i++)
    {
      const FiniteElement *fe = fes->GetFE(i);
      const int order = 2*fe->GetOrder()+3; // <----------
      const IntegrationRule *ir = &IntRules.Get(fe->GetGeomType(), order);
      const int dof = fe->GetDof();
      Trans = fes->GetElementTransformation(i);
      Vector stressh(tdim), strainh(tdim);	//approché
      Vector stress(tdim), strain(tdim);	//exacte
      DenseMatrix Ch;
      for (int j = 0; j < ir->GetNPoints(); j++)
	{
	  const IntegrationPoint &ip = ir->IntPoint(j);
	  Trans->SetIntPoint(&ip);
	  double w = Trans->Weight() * ip.weight;

	  ComputeStrain(*Trans, ip, x, i, strainh);//approx
	  Stress_exacte_coef.Eval(stress,*Trans,ip);//exact
	  Elasticy_mat(*Trans,ip,dim,lambdah,muh,Ch);
	  Ch.Mult(strainh,stressh);	
	  Ch.Invert();
	  Ch.Mult(stress,strain);

	  strainh -= strain;
	  stressh -= stress;

	  double pdc=0.0;
	  for (int k = 0; k< dim; k++)
	    pdc += strainh(k)*strainh(k);
	  for (int k = dim; k < dim*(dim+1)/2; k++)
	    pdc += 2*strainh(k)*strainh(k); 

	  energy += w * pdc;
	}
    }
  if(energy>0.0){
    return sqrt(energy);}
  else{
    cout<<"Negative Energy error"<<endl;
    exit(0);
  }
}

//===================== Matrice élasticité =====================
void Elasticy_mat(ElementTransformation &T,const IntegrationPoint &ip, 
		  int dim, Coefficient &lambda, Coefficient &mu_func, DenseMatrix &C){
  double M = mu_func.Eval(T, ip);
  double L = lambda.Eval(T, ip);

  C.SetSize(dim*(dim+1)/2,dim*(dim+1)/2);
  C = 0.;
  for (int k = 0; k< dim; k++)
    {
      // Extra-diagonal terms
      for (int l = 0; l< dim; l++)
	C(k,l) = L;
      // Diagonal terms
      C(k,k) = L+2*M;
    }
  // Diagonal terms
  for (int k = dim; k < dim*(dim+1)/2; k++)
    C(k,k) = 2*M;
}

//===================== Déformation =====================
void ComputeStrain(ElementTransformation &T,const IntegrationPoint &ip,
		   GridFunction &x, int elem,  Vector &strain){
  FiniteElementSpace *fes = x.FESpace();
  Array<int> udofs;
  const FiniteElement *fe = fes->GetFE(elem);
  const int dof = fe->GetDof();
  const int dim = fe->GetDim();
  const int tdim = dim*(dim+1)/2; // num. entries in a symmetric tensor
  DenseMatrix dshape(dof, dim);
  DenseMatrix gh(dim, dim),grad (dim, dim);
  fes->GetElementVDofs(elem, udofs);
  DenseMatrix loc_data(dof, dim);
  for (int s=0 ; s<dim ; s++)
    {
      Array<int> udofs_tmp(dof);
      udofs_tmp = 0.;
      for(int j=0 ; j<dof ; j++){
	udofs_tmp[j] = udofs[j+dof*s];}

      Vector loc_data_tmp;
      x.GetSubVector(udofs_tmp, loc_data_tmp);

      for (int j=0 ; j<dof ; j++)
	loc_data(j,s) = loc_data_tmp(j);
    }
  double w = T.Weight() * ip.weight;
  fe->CalcDShape(ip, dshape);
  MultAtB(loc_data, dshape, gh);
  Mult(gh, T.InverseJacobian(), grad);
  Strain_(grad, strain);
}
