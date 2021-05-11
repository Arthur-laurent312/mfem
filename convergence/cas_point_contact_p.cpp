// Compile with: make cas_point_contact
//
//Cas test avec solution analytique. Plaque 2D avec point de contact en (0,0).
//Calcul de la déformation en élasticité linéaire. 
//Version parallèle

//Ref: http://e6.ijs.si/medusa/wiki/index.php/Point_contact

#include "mfem.hpp"
#include <fstream>
#include <iostream>

using namespace std;
using namespace mfem;

static constexpr double Force = -1.;
static constexpr double E = 1000.;	//Module de young
static constexpr double nu = 0.25;	//Coef de poisson
static constexpr double lambda = E*nu/((1.+nu)*(1.-2.*nu));
static constexpr double mu = E/(2.*(1.+nu));	//coef de Lamé

void sol_exact(const Vector &, Vector &);
double ux_exact(const Vector &);
double uy_exact(const Vector &);
void Grad_Exact(const Vector &, DenseMatrix &);
void GradExact_x(const Vector &, Vector &);
void GradExact_y(const Vector &, Vector &);

double ComputeEnergyNorm(GridFunction &,
			 Coefficient &, Coefficient &);

void Elasticy_mat(ElementTransformation &,const IntegrationPoint &, 
		  int, Coefficient &, Coefficient &, DenseMatrix &);

double ComputeH1Norm(GridFunction &);

void ComputeStress(ElementTransformation &,const IntegrationPoint &,
		   GridFunction &, int,  Vector &);

int main(int argc, char *argv[])
{


  // Initialize MPI.
  int num_procs, myid;
  MPI_Init(&argc, &argv);
  MPI_Comm_size(MPI_COMM_WORLD, &num_procs);
  MPI_Comm_rank(MPI_COMM_WORLD, &myid);
 // Parse command-line options.
  int order=2;
  const char *mesh_file = "carre.msh";
  bool static_cond = false;
  bool amg_elast = 0;
  bool reorder_space = false;
  bool solver=true;
  OptionsParser args(argc, argv);
  args.AddOption(&mesh_file, "-m", "--mesh",
		 "Mesh file to use.");
  args.AddOption(&order, "-o", "--order",
		 "Finite element order (polynomial degree).");
  args.AddOption(&amg_elast, "-elast", "--amg-for-elasticity", "-sys",
		 "--amg-for-systems",
		 "Use the special AMG elasticity solver (GM/LN approaches), "
		 "or standard AMG for systems (unknown approach).");
  args.AddOption(&static_cond, "-sc", "--static-condensation", "-no-sc",
		 "--no-static-condensation", "Enable static condensation.");
  args.AddOption(&reorder_space, "-nodes", "--by-nodes", "-vdim", "--by-vdim",
		 "Use byNODES ordering of vector space instead of byVDIM");
  args.AddOption(&solver, "-sol", "--Itératif", "-Direct",
		 "--Solver Itératif", "Solver direct.");   

  int ref_levels = 6;

  args.Parse();
  if (!args.Good())
    {
      if (myid == 0)
	{
	  args.PrintUsage(cout);
	}
      MPI_Finalize();
      return 1;
    }
  if (myid == 0)
    {
      args.PrintOptions(cout);
    }

  //Lecture du malliage
  Mesh *mesh = new Mesh(mesh_file, 1, 1);
  int dim = mesh->Dimension();

  //  refine the mesh to increase the resolution. In this example we do
  //    'ref_levels' of uniform refinement. We choose 'ref_levels' to be the
  //    largest number that gives a final mesh with no more than 5,000
  //    elements.
  //   int ref_levels = (int)floor(log(5000./mesh->GetNE())/log(2.)/dim);

  for (int l = 0; l < ref_levels; l++)
    {
      mesh->UniformRefinement();
    }

  //Define parallel mesh by a partitioning of the serial mesh.
  ParMesh *pmesh = new ParMesh(MPI_COMM_WORLD, *mesh);
/*
    int par_ref_levels = 1;
    for (int l = 0; l < par_ref_levels; l++)
      {
	pmesh->UniformRefinement();
      }
*/
  //  Define a finite element space on the mesh. Here we use vector finite
  //    elements, i.e. dim copies of a scalar finite element space. The vector
  //    dimension is specified by the last argument of the FiniteElementSpace
  //    constructor. For NURBS meshes, we use the (degree elevated) NURBS space
  //    associated with the mesh nodes.
  FiniteElementCollection *fec;
  ParFiniteElementSpace *fespace;
  const bool use_nodal_fespace = pmesh->NURBSext && !amg_elast;
  if (use_nodal_fespace)
    {
      fec = NULL;
      fespace = (ParFiniteElementSpace *)pmesh->GetNodes()->FESpace();
    }
  else
    {
      fec = new H1_FECollection(order, dim);
      if (reorder_space)
	{
	  fespace = new ParFiniteElementSpace(pmesh, fec, dim, Ordering::byNODES);
	}
      else
	{
	  fespace = new ParFiniteElementSpace(pmesh, fec, dim, Ordering::byVDIM);
	}
    }
  HYPRE_Int size = fespace->GlobalTrueVSize();
  if (myid == 0)
    {
      cout << "Number of finite element unknowns: " << size << endl
           << "Assembling: " << flush;
    }

  // Determine the list of true (i.e. conforming) essential boundary dofs.
  //    In this example, the boundary conditions are defined by marking only
  //    boundary attribute 1 from the mesh as essential and converting it to a
  //    list of true dofs.


  // List of True DoFs : Define (here) Dirichlet conditions
  Array<int> ess_tdof_list;
  Array<int> ess_bdr(pmesh->bdr_attributes.Max());
  ess_bdr = 1;		//On sélectionne tous les bords
  fespace->GetEssentialTrueDofs(ess_bdr, ess_tdof_list);

  // Set up the linear form b(.) which corresponds to the right-hand side of
  //    the FEM linear system. In this case, b_i equals the boundary integral
  //    of f*phi_i where f represents a "pull down" force on the Neumann part
  //    of the boundary and phi_i are the basis functions in the finite element
  //    fespace. The force is defined by the VectorArrayCoefficient object f,
  //    which is a vector of Coefficient objects. The fact that f is non-zero
  //    on boundary attribute 2 is indicated by the use of piece-wise constants
  //    coefficient for its last component.

  VectorArrayCoefficient f(dim);
  for (int i = 0; i < dim-1; i++)
    {
      f.Set(i, new ConstantCoefficient(0.0));
    }

  ParLinearForm *b = new ParLinearForm(fespace);
  b->Assemble();
  if (myid == 0)
    {
      cout << "r.h.s. ... " << flush;
    }

  // Define the solution vector x as a finite element grid function
  //    corresponding to fespace. Initialize x with initial guess of zero,
  //    which satisfies the boundary conditions.
  ParGridFunction x(fespace);
  VectorFunctionCoefficient Boundary_Dirichlet_coef(dim, sol_exact);
  x = 0.0;
  // To use if there are different Dirichlet conditions.
  // Beware, the values of dirichlet boundary conditions are set here !
  //Projète solution exacte sur les bords
  x.ProjectBdrCoefficient(Boundary_Dirichlet_coef, ess_bdr);	

  // Set up the bilinear form a(.,.) on the finite element space
  //    corresponding to the linear elasticity integrator with piece-wise
  //    constants coefficient lambda and mu.

  ConstantCoefficient mu_func(mu);
  ConstantCoefficient lambda_func(lambda);

  ParBilinearForm *a = new ParBilinearForm(fespace);
  BilinearFormIntegrator *integ = new ElasticityIntegrator(lambda_func, mu_func);
  a->AddDomainIntegrator(integ);

  // Assemble the parallel bilinear form and the corresponding linear
  //     system, applying any necessary transformations such as: parallel
  //     assembly, eliminating boundary conditions, applying conforming
  //     constraints for non-conforming AMR, static condensation, etc.
  if (myid == 0) { cout << "matrix ... " << flush; }
  if (static_cond) { a->EnableStaticCondensation(); }
  a->Assemble();

  HypreParMatrix A;
  Vector B, X;

  a->FormLinearSystem(ess_tdof_list, x, *b, A, X, B);
  if (myid == 0)
    {
      cout << "done." << endl;
      cout << "Size of linear system: " << A.GetGlobalNumRows() << endl;
    }

  // 13. Define and apply a parallel PCG solver for A X = B with the BoomerAMG
  //     preconditioner from hypre.

  if(solver){
    HypreBoomerAMG *amg = new HypreBoomerAMG(A);
    if (amg_elast && !a->StaticCondensationIsEnabled())
      {
	amg->SetElasticityOptions(fespace);
      }
    else
      {
	amg->SetSystemsOptions(dim, reorder_space);
      }
    HyprePCG *pcg = new HyprePCG(A);
    pcg->SetTol(1e-15);
    pcg->SetMaxIter(5000);
    pcg->SetPrintLevel(2);
    pcg->SetPreconditioner(*amg);
    pcg->Mult(B, X);
    delete pcg;
    delete amg;
  }
  else{
    cout<<"Solver direct non implémenté"<<endl;
  }   


  // 12. Recover the solution as a finite element grid function.
  a->RecoverFEMSolution(X, *b, x);


  // 13. For non-NURBS meshes, make the mesh curved based on the finite element
  //     space. This means that we define the mesh elements through a fespace
  //     based transformation of the reference element. This allows us to save
  //     the displaced mesh as a curved mesh when using high-order finite
  //     element displacement field. We assume that the initial mesh (read from
  //     the file) is not higher order curved mesh c<ompared to the chosen FE
  //     space.
  if (!use_nodal_fespace)
    {
      pmesh->SetNodalFESpace(fespace);
    }

      // Compute errors
      double ener_error = ComputeEnergyNorm(x, lambda_func, mu_func);
      VectorFunctionCoefficient sol_exact_coef(dim, sol_exact);
      double L2_error = x.ComputeL2Error(sol_exact_coef);
      double H1_error = ComputeH1Norm(x);
  if (myid == 0)
    {
      cout<<"Erreur en Norme L2: "<<L2_error<<endl;
      cout<<"Erreur en Norme Énergie: "<<ener_error<<endl;
      cout<<"Erreur en Norme H1: "<<H1_error<<endl;  
      cout<<"Taille de maille: "<<h<<endl;    
      cout << "numbers of elements: " << pmesh->GetNE() <<endl;
    }
  //  Free the used memory.
  delete a;
  delete b;
  if (fec)
    {
      delete fec;
    }
  delete mesh;
  delete pmesh;

  MPI_Finalize();
  /*
  //Save in Praview format
  if (ref1==0){
  ParGridFunction diff(fespace);
  ParGridFunction ex1(fespace);
  diff.ProjectCoefficient(sol_exact_coef);
  ex1.ProjectCoefficient(sol_exact_coef);
  diff-= x;
 
  ParaViewDataCollection paraview_dc("Example2", mesh);
  paraview_dc.SetPrefixPath("ParaView");
  paraview_dc.SetLevelsOfDetail(order+1);
  paraview_dc.SetCycle(0);
  paraview_dc.SetDataFormat(VTKFormat::BINARY);
  paraview_dc.SetHighOrderOutput(true);
  paraview_dc.SetTime(0.0); // set the time
  paraview_dc.RegisterField("numerical_solution",&x);
  paraview_dc.RegisterField("diff-exact_solution",&diff);
  paraview_dc.RegisterField("exact_solution",&ex1);
  paraview_dc.Save();	
  }
  */

  return 0;
}

//===================== Solution exacte =====================
void sol_exact(const Vector &x, Vector &u)
{
  double r2 = x(1)*x(1)+x(0)*x(0);
  double pi = M_PI;
  u(0) = -Force/(4*pi*mu)*(2*x(0)*x(1)/r2 + 2*mu/(lambda+mu)*atan2(x(1),x(0)));
  u(1) = -Force/(4*pi*mu)*((x(1)*x(1)-x(0)*x(0))/r2 - 
			   (lambda+2*mu)*log(r2)/(lambda+mu));
}
double ux_exact(const Vector &x)
{
  double r2 = x(1)*x(1)+x(0)*x(0);
  double pi = M_PI;
  return -Force/(4*pi*mu)*(2*x(0)*x(1)/r2 + 2*mu/(lambda+mu)*atan2(x(1),x(0)));
}

double uy_exact(const Vector &x)
{
  double r2 = x(1)*x(1)+x(0)*x(0);
  double pi = M_PI;
  return -Force/(4*pi*mu)*((x(1)*x(1)-x(0)*x(0))/r2 - 
			   (lambda+2*mu)*log(r2)/(lambda+mu));
}
//===================== Grad exacte =====================
void Grad_Exact(const Vector &x, DenseMatrix &grad)
{
  double r2 = x(1)*x(1)+x(0)*x(0);
  double pi = M_PI;

  grad(0,0) = -Force/(4*pi*mu)*(2*x(1)*(x(1)*x(1)-x(0)*x(0))/(r2*r2)-2*mu/(lambda+mu)*x(1)/r2);
  grad(0,1) = -Force/(4*pi*mu)*(2*x(0)*(x(0)*x(0)-x(1)*x(1))/(r2*r2) + 2*mu/(lambda+mu)*x(0)/r2);
  grad(1,1) = -Force/(4*pi*mu)*(4*x(1)*x(0)*x(0)/(r2*r2) - 
			(lambda+2*mu)/(lambda+mu)*2*x(1)/r2);
  grad(1,0) = -Force/(4*pi*mu)*(-4*x(1)*x(1)*x(0)/(r2*r2) - 
			(lambda+2*mu)/(lambda+mu)*2*x(0)/r2);
}
void GradExact_x(const Vector &x, Vector &grad)
{
  double r2 = x(1)*x(1)+x(0)*x(0);
  double pi = M_PI;
  grad(0) = -Force/(4*pi*mu)*(2*x(1)*(x(1)*x(1)-x(0)*x(0))/(r2*r2)-2*mu/(lambda+mu)*x(1)/r2);
  grad(1) = -Force/(4*pi*mu)*(2*x(0)*(x(0)*x(0)-x(1)*x(1))/(r2*r2) + 2*mu/(lambda+mu)*x(0)/r2);
}
void GradExact_y(const Vector &x, Vector &grad)
{
  double r2 = x(1)*x(1)+x(0)*x(0);
  double pi = M_PI;
  grad(1) = -Force/(4*pi*mu)*(4*x(1)*x(0)*x(0)/(r2*r2) - 
		(lambda+2*mu)/(lambda+mu)*2*x(1)/r2);
  grad(0) = -Force/(4*pi*mu)*(-4*x(1)*x(1)*x(0)/(r2*r2) - 
		(lambda+2*mu)/(lambda+mu)*2*x(0)/r2);
}

//===================== Erreur en norme H1 =====================
double ComputeH1Norm(GridFunction &x){
  FiniteElementSpace *fes = x.FESpace();
  int dim = fes->GetMesh()->SpaceDimension();

  MatrixFunctionCoefficient grad_exact_coef (dim, Grad_Exact);
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
      const int tdim = dim*(dim+1)/2; // num. entries in a symmetric tensor

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
  return (error < 0.0) ? -sqrt(-error) : sqrt(error);
}

//==============Erreur en Norme energie ===================
double ComputeEnergyNorm(GridFunction &x,
			 Coefficient &lambdah, Coefficient &muh)
{
  FiniteElementSpace *fes = x.FESpace();
  int dim = fes->GetMesh()->SpaceDimension();
  VectorFunctionCoefficient sol_exact_coef (dim, sol_exact);
  GridFunction ex(fes);
  ex.ProjectCoefficient(sol_exact_coef);
  
  ConstantCoefficient lambda_func(lambda);
  ConstantCoefficient mu_func(mu);
  ElementTransformation *Trans;

  double energy = 0.0;
  for (int i = 0; i < fes->GetNE() ; i++)
    {
      const FiniteElement *fe = fes->GetFE(i);
      const int order = 2*fe->GetOrder()+3; // <----------
      const IntegrationRule *ir = &IntRules.Get(fe->GetGeomType(), order);
      const int dof = fe->GetDof();
      const int tdim = dim*(dim+1)/2; // num. entries in a symmetric tensor
      Trans = fes->GetElementTransformation(i);
      Vector stressh(tdim), strainh(tdim);	//approché
      Vector stress(tdim), strain(tdim);	//exacte
      DenseMatrix C,Ch;
      for (int j = 0; j < ir->GetNPoints(); j++)
	{
	  const IntegrationPoint &ip = ir->IntPoint(j);
	  Trans->SetIntPoint(&ip);
	  double w = Trans->Weight() * ip.weight;

	  ComputeStress(*Trans, ip, x, i, stressh);
	  ComputeStress(*Trans, ip, ex, i, stress);

	  //======= Stress vectors ========
	  Elasticy_mat(*Trans,ip,dim,lambda_func,mu_func,C);
	  Elasticy_mat(*Trans,ip,dim,lambdah,muh,Ch);

	  Ch.Mult(stressh,strainh);	//approx
	  C.Mult(stress,strain);	//exacte

	  strainh -= strain;
	  stressh -= stress;

	  double pdc=0.0;
	  for (int k = 0; k< dim; k++)
	    pdc += strainh(k)*stressh(k);

	  for (int k = dim; k < dim*(dim+1)/2; k++)
	    pdc += 2*strainh(k)*stressh(k);

	  energy += w * pdc;
	}
    }
  return (energy < 0.0) ? -sqrt(-energy) : sqrt(energy);
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
void ComputeStress(ElementTransformation &T,const IntegrationPoint &ip,
		   GridFunction &x, int elem,  Vector &stress){
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
  stress(0)=grad(0,0);
  stress(1)=grad(1,1);
  if(dim==2){
    stress(2)=0.5*(grad(1,0)+grad(0,1));
  }
  else if(dim==3){
    stress(2)=grad(2,2);
    stress(3)=0.5*(grad(1,0)+grad(0,1));
    stress(4)=0.5*(grad(2,0)+grad(0,2));
    stress(5)=0.5*(grad(2,1)+grad(1,2));
  }
  else{
    cout<<"dimention not suported"<<endl;}
}

