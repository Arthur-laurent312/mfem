//Compile with: make cas_flexion
//
//Cas test elsatique linéaire d'une poutre isotrope en flexion.
//Résolution de l'équation: -div(\sigma(u))=f
//La charge P est appliqué suivant un profil parabolique F.
//Condition de dirichlet u_h=u_{exacte} à l'encastrement.
//Solution analytique disponible.
//
//Ref: https://www.sciencedirect.com/science/article
//		/abs/pii/S0168874X08000140?via%3Dihub#bib1

//Calculs avec plusieurs maillage possible
//GL: quels maillages ? Exemple ...

#include "mfem.hpp"
#include <fstream>
#include <iostream>

using namespace std;
using namespace mfem;

//GL: commenter ces constantes 
static constexpr double pull_force = -1.;
static constexpr double L = 8.0;
static constexpr double D =1.0;
static constexpr double E_old = 1000.;
static constexpr double nu_old = 0.25;
static constexpr double I = D*D*D*D/12.;
static constexpr double lambda_old  = E_old*nu_old/((1.+nu_old)*(1.-2.*nu_old));
static constexpr double mu_old = E_old/(2.*(1.+nu_old));
static constexpr double E = E_old/(1.-nu_old*nu_old);
static constexpr double nu = nu_old/(1.-nu_old);		//Passage en plane stress
static constexpr double lambda = E*nu/((1.+nu)*(1.-2.*nu));
static constexpr double mu = E/(2.*(1.+nu));

//Solution exacte
void sol_exact(const Vector &, Vector &);
void grad_exact(const Vector &, DenseMatrix &);
void Stress_exacte(const Vector &, Vector &);
void Strain_exacte(const Vector &, Vector &);

//Chargement
double F(const Vector &);

//Erreur en norme H1
double ComputeH1Norm(GridFunction &);
//Erreur en norme energy
double ComputeEnergyNorm(GridFunction &,
			 Coefficient &, Coefficient &);

void Elasticy_mat(ElementTransformation &,const IntegrationPoint &, int,
		  Coefficient &, Coefficient &, DenseMatrix &);
void ComputeStrain(ElementTransformation &,const IntegrationPoint &,
		   GridFunction &, int,  Vector &);

int main(int argc, char *argv[])
{
  DenseMatrix slope_l2, slope_ener, slope_grad;
  double err_tmp_ener = 0, err_tmp_l2=0, err_tmp_grad=0;
  double h_tmp = 0.;
  int iter = 0;

  // Parse command-line options.
  const char *mesh_file = "beam-bar.msh";
  bool static_cond = false;
  int order = 2;
  int rep = 7 ;
  bool iterative = true;

  OptionsParser args(argc, argv);
  args.AddOption(&mesh_file, "-m", "--mesh",
		 "Mesh file to use.");
  args.AddOption(&order, "-o", "--order",
		 "Finite element order (polynomial degree).");
  args.AddOption(&rep, "-r", "--Repetition",
		 "Nombre de maillage.");
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
  slope_l2.SetSize(rep-1,3);
  slope_grad.SetSize(rep-1,3);

  string const err_energy("err_flexion_ordre2.txt");
  ofstream err_energy_flux(err_energy.c_str());
  if (!err_energy_flux.is_open()) {
    cout << "Problem in openning file" << endl;
    exit(0);
  }

  {
    cout << "Openning file" << endl;
    Mesh *mesh = new Mesh(mesh_file, 1, 1);
    // 2. Read the mesh from the given mesh file. We can handle triangular,
    //    quadrilateral, tetrahedral or hexahedral elements with the same code.
    int dim = mesh->Dimension();
  
    for (int ref_levels=1; ref_levels<rep; ref_levels++){ 
  
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
      Array<int> ess_tdof_list;
      Array<int> ess_bdr(mesh->bdr_attributes.Max());
      ess_bdr = 0;
      ess_bdr[0] =1;
      fespace->GetEssentialTrueDofs(ess_bdr, ess_tdof_list);
  
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
      VectorFunctionCoefficient Boundary_Dirichlet_coef(dim, sol_exact);
      x = 0.0;
      // To use if there are different Dirichlet conditions.
      // Beware, the values of dirichlet boundary conditions are set here !
      x.ProjectBdrCoefficient(Boundary_Dirichlet_coef, ess_bdr);
  
      // 9. Set up the bilinear form a(.,.) on the finite element space
      //    corresponding to the linear elasticity integrator with piece-wise
      //    constants coefficient lambda and mu.
  
      ConstantCoefficient mu_func(mu_old);
      ConstantCoefficient lambda_func(lambda_old);
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
  	PCG(A, M, B, X, 2, 50000, 1e-20, 0.0);
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

      // Compute error
      double ener_error = ComputeEnergyNorm(x, lambda_func, mu_func);
      double err_H1 = ComputeH1Norm(x);
      VectorFunctionCoefficient sol_exact_coef(dim, sol_exact);
      double L2_error = x.ComputeL2Error(sol_exact_coef);
  	
      cout << "\nL2 norm of error: " << L2_error << endl;
      cout << "Energy norm of error: " << ener_error << endl;
      cout << "H1 norm of error: " << err_H1 << endl << endl;
  
      double h = mesh->GetElementSize(0);
      //Compute the slope
      slope_l2(iter,0) = log(err_tmp_l2/L2_error) / log(h_tmp/h);
      slope_l2(iter,1) = h;
      slope_l2(iter,2) = L2_error;
      err_tmp_l2 = L2_error;
      //Compute the slope
      slope_ener(iter,0) = log(err_tmp_ener/ener_error) / log(h_tmp/h);
      slope_ener(iter,1) = h;
      slope_ener(iter,2) = ener_error;
      err_tmp_ener = ener_error;
      //Compute the slope
      slope_grad(iter,0) = log(err_tmp_grad/err_H1) / log(h_tmp/h);
      slope_grad(iter,1) = h;
      slope_grad(iter,2) = err_H1;
      err_tmp_grad = err_H1;
      h_tmp = h;
      iter++;
  
      //Save in Praview format
      GridFunction ex(fespace);
      ex.ProjectCoefficient(sol_exact_coef);
      GridFunction diff(fespace);
      diff.ProjectCoefficient(sol_exact_coef);
      diff -= x;
  
      ParaViewDataCollection paraview_dc("Flexion", mesh);
      paraview_dc.SetPrefixPath("ParaView");
      paraview_dc.SetLevelsOfDetail(order+1);
      paraview_dc.SetCycle(0);
      paraview_dc.SetDataFormat(VTKFormat::BINARY);
      paraview_dc.SetHighOrderOutput(true);
      paraview_dc.SetTime(iter*1.0); // set the time
      paraview_dc.RegisterField("numerical_solution",&x);
      paraview_dc.RegisterField("diff-exact_solution",&diff);
      paraview_dc.RegisterField("exact_solution",&ex);
      paraview_dc.Save();	
  
      delete a;
      delete b;
      if (fec) {
  	delete fespace;
  	delete fec;
      }
      //    refine the mesh to increase the resolution. In this example we do
      //    'ref_levels' of uniform refinement. We choose 'ref_levels' to be the
      //    largest number that gives a final mesh with no more than 5,000
      //    elements.
      mesh->UniformRefinement();
    }  //end loop mesh
    //Affichage des normes et pentes.
    cout<<endl;
    cout<<"Erreur en norme:"<<endl;
    err_energy_flux << "Erreur en norme:"<<endl;
    for (int i=1; i<iter; i++){
      err_energy_flux << "L2: " << slope_l2(i,2)<<" H1: "<<slope_grad(i,2) 
		      << " Energie: " << slope_ener(i,2)<<" Taille de maille= "
		      <<slope_l2(i,1)<<endl;
      cout << "L2: " << slope_l2(i,2)<<" H1: "<<slope_grad(i,2) 
  	   << " Energie: " << slope_ener(i,2)<<" Taille de maille= "
  	   <<slope_l2(i,1)<<endl;}

    cout<<endl;
    err_energy_flux<<endl;
    cout<<"Pente de convergence:"<<endl;
    err_energy_flux << "Pente de convergence:"<<endl;
    for (int i=1; i<iter; i++){
      err_energy_flux<<"Pente L2: " << slope_l2(i,0)<<" H1: "<<slope_grad(i,0) 
		     << " Energie: " << slope_ener(i,0)<<" Taille de maille= "<<
  	slope_l2(i,1)<<endl;
      cout << "Pente L2: " << slope_l2(i,0)<<" H1: "<<slope_grad(i,0) 
  	   << " Energie: " << slope_ener(i,0)<<" Taille de maille= "<<
  	slope_l2(i,1)<<endl;}
    cout<<endl;
  }//end flux .txt
  return 0;
}

//===================== Solution exacte =====================
void sol_exact(const Vector &x, Vector &u)
{
  double y = x(1)-D*0.5;


  u(0) = -pull_force*y/(6.*E*I) * ((6.*L-3.*x(0))*x(0) + (2.+nu)*(y*y - D*D/4.));
  u(1) = pull_force/(6.*E*I) * (3.*nu*y*y*(L-x(0)) + (3.*L-x(0))*x(0)*x(0)+ 
				(4.+5.*nu)*D*D*x(0)/4.);
}
//===================== Gradient exacte =====================
void grad_exact(const Vector &x, DenseMatrix &grad)
{
  double y = x(1)-D*0.5;

  grad(0,0) = -pull_force*y/(E*I)*(L-x(0));
  grad(0,1) = -pull_force/(6.*E*I)*((6*L-3*x(0))*x(0) + (2+nu)*(y*y-D*D/4)
				    + 2*y*y*(2+nu));
  grad(1,1) = pull_force*y*nu/(E*I)*(L-x(0));
  grad(1,0) = pull_force/(6.*E*I)*(2*x(0)*(3*L-x(0))-x(0)*x(0) -
				   3*nu*y*y+(4+5*nu)*D*D/4);
}
//===================== Stress exacte =====================
void Stress_exacte(const Vector &x, Vector &stress)
{
  double y = x(1)-D*0.5;
  stress(0) = -pull_force*(L-x(0))*y/I;
  stress(1) = 0.0;
  stress(2) = pull_force/(2*I)*(D*D/4-y*y);
}
//===================== Strain exacte =====================
void Strain_exacte(const Vector &x, Vector &strain)
{
  int dim = x.Size();
  DenseMatrix grad(dim,dim);
  double y = x(1)-D*0.5;
  grad(0,0) = -pull_force*y/(E*I)*(L-x(0));
  grad(0,1) = -pull_force/(6.*E*I)*((6*L-3*x(0))*x(0) + (2+nu)*(y*y-D*D/4)
				    + 2*y*y*(2+nu));
  grad(1,1) = pull_force*y*nu/(E*I)*(L-x(0));
  grad(1,0) = pull_force/(6.*E*I)*(2*x(0)*(3*L-x(0))-x(0)*x(0) -
				   3*nu*y*y+(4+5*nu)*D*D/4);
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
  if(x(0) >= 8.-1.e-6){
    force = pull_force/(2.*I)*(pow(D*0.5,2) - pow((x(1)-D*0.5),2));
  } else {
    force = 0.;
  }

  return force;
}

//===================== Erreur en norme H1 =====================
double ComputeH1Norm(GridFunction &x){
  FiniteElementSpace *fes = x.FESpace();
  int dim = fes->GetMesh()->SpaceDimension();

  MatrixFunctionCoefficient grad_exact_coef (dim, grad_exact);
  ElementTransformation *Trans;
  DenseMatrix grad, gradh;
  double error = 0.0;
  Array<int> udofs;
  for (int i = 0; i < fes->GetNE() ; i++)
    {
      const FiniteElement *fe = fes->GetFE(i);
      const int order = 2*fe->GetOrder() + 3;   //<----------
      const IntegrationRule *ir = &(IntRules.Get(fe->GetGeomType(), order));
      Trans = fes->GetElementTransformation(i);
      const int dof = fe->GetDof();

      DenseMatrix dshape(dof, dim);
      DenseMatrix gh(dim, dim),gradh (dim, dim),grad(dim,dim);
      fes->GetElementVDofs(i, udofs);
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
      for (int j = 0; j < ir->GetNPoints(); j++)
	{
	  const IntegrationPoint &ip = ir->IntPoint(j);
	  Trans->SetIntPoint(&ip);
	  double w = Trans->Weight() * ip.weight;
	  fe->CalcDShape(ip, dshape);
	  MultAtB(loc_data, dshape, gh);
	  Mult(gh, Trans->InverseJacobian(), gradh);
	  grad_exact_coef.Eval(grad,*Trans,ip);
	  grad -= gradh;
	  error += w * grad.FNorm2();
	}			
    }
  if(error>0.0){
    return sqrt(error);}
  else{
    cout<<"Negative H1 error"<<endl;
    exit(0);}
}

//==============Erreur en Norme energie ===================
double ComputeEnergyNorm(GridFunction &x,
			 Coefficient &lambdah, Coefficient &muh)
{
  FiniteElementSpace *fes = x.FESpace();
  int dim = fes->GetMesh()->SpaceDimension();
  const int tdim = dim*(dim+1)/2; // num. entries in a symmetric tensor
  VectorFunctionCoefficient Strain_exacte_coef(tdim,Strain_exacte);
  VectorFunctionCoefficient Stress_exacte_coef(tdim,Stress_exacte);
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

	  //======= strain vectors ========
	  ComputeStrain(*Trans, ip, x, i, strainh);//approx
	  Strain_exacte_coef.Eval(strain,*Trans,ip);//exact

	  //======= Stress vectors ========
	  Elasticy_mat(*Trans,ip,dim,lambdah,muh,Ch);
	  Ch.Mult(strainh,stressh);			//approx
	  Stress_exacte_coef.Eval(stress,*Trans,ip);//exact

	  strain -= strainh;
	  stress -= stressh;

	  double pdc=0.0;
	  for (int k = 0; k< dim; k++)
	    pdc += strain(k)*stress(k);
	  for (int k = dim; k < dim*(dim+1)/2; k++)
	    pdc += 2*strain(k)*stress(k);

	  energy += w * pdc;
	}
    }
  if(energy>0.0){
    return sqrt(energy);}
  else{
    cout<<"Negative Energy error"<<endl;
    exit(0);}
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
    C(k,k) = M;
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
    cout<<"dimention not suported"<<endl;}
}

