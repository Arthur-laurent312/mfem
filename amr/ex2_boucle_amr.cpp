//
// Compile with: make ex2_boucle_amr
//
// Sample runs:  ex2 -m ../data/beam-tri.mesh
//               ex2 -m ../data/beam-quad.mesh
//               ex2 -m ../data/beam-tet.mesh
//               ex2 -m ../data/beam-hex.mesh
//               ex2 -m ../data/beam-wedge.mesh
//               ex2 -m ../data/beam-quad.mesh -o 3 -sc
//               ex2 -m ../data/beam-quad-nurbs.mesh
//               ex2 -m ../data/beam-hex-nurbs.mesh
//
//			Problème élastique a plusieurs chargements et rafinement de maillage 
//				entre les chargements

#include "mfem.hpp"
#include <fstream>
#include <iostream>

using namespace std;
using namespace mfem;

int main(int argc, char *argv[])
{
   // 1. Parse command-line options.
   const char *mesh_file = "../data/beam-tri.mesh";
   int order = 1;
   int max_dofs = 50000;
   bool static_cond = false;
   bool visualization = 1;
   int boucle = 0;
   OptionsParser args(argc, argv);
   args.AddOption(&mesh_file, "-m", "--mesh",
                  "Mesh file to use.");
   args.AddOption(&order, "-o", "--order",
                  "Finite element order (polynomial degree).");
   args.AddOption(&static_cond, "-sc", "--static-condensation", "-no-sc",
                  "--no-static-condensation", "Enable static condensation.");
   args.AddOption(&max_dofs, "-md", "--max-dofs",
                  "Stop after reaching this many degrees of freedom.");
   args.AddOption(&visualization, "-vis", "--visualization", "-no-vis",
                  "--no-visualization",
                  "Enable or disable GLVis visualization.");
   args.AddOption(&boucle, "-r", "--rep", "-Nb-répét",
                  "Boucle de répétition du calcul.");
   args.Parse();
   if (!args.Good())
   {
      args.PrintUsage(cout);
      return 1;
   }
   args.PrintOptions(cout);

   // 2. Read the mesh from the given mesh file. We can handle triangular,
   //    quadrilateral, tetrahedral or hexahedral elements with the same code.
   Mesh *mesh = new Mesh(mesh_file, 1, 1);
   int dim = mesh->Dimension();
   int sdim = mesh->SpaceDimension();

   // 3. Select the order of the finite element discretization space. For NURBS
   //    meshes, we increase the order by degree elevation.
   if (mesh->NURBSext)
   {
      mesh->DegreeElevate(order, order);
   }

   // 4. Refine the mesh to increase the resolution. In this example we do
   //    'ref_levels' of uniform refinement. We choose 'ref_levels' to be the
   //    largest number that gives a final mesh with no more than 5,000
   //    elements.
   {
      int ref_levels =
         (int)floor(log(5000./mesh->GetNE())/log(2.)/dim);
      for (int l = 0; l < ref_levels; l++)
      {
         mesh->UniformRefinement();
      }
   }

   // 5. Define a finite element space on the mesh. Here we use vector finite
   //    elements, i.e. dim copies of a scalar finite element space. The vector
   //    dimension is specified by the last argument of the FiniteElementSpace
   //    constructor. For NURBS meshes, we use the (degree elevated) NURBS space
   //    associated with the mesh nodes.
   FiniteElementCollection *fec;
   FiniteElementSpace *fespace;
   if (mesh->NURBSext)
   {
      fec = NULL;
      fespace = mesh->GetNodes()->FESpace();
   }
   else
   {
      fec = new H1_FECollection(order, dim);
      fespace = new FiniteElementSpace(mesh, fec, dim);
   }
   cout << "Number of finite element unknowns: " << fespace->GetTrueVSize()
        << endl << "Assembling: " << flush;

   // 6. Determine the list of true (i.e. conforming) essential boundary dofs.
   //    In this example, the boundary conditions are defined by marking only
   //    boundary attribute 1 from the mesh as essential and converting it to a
   //    list of true dofs.


   Array<int> ess_tdof_list, ess_bdr(mesh->bdr_attributes.Max());
   ess_bdr = 0;
   ess_bdr[0] = 1;
   // 7. Set up the linear form b(.) which corresponds to the right-hand side of
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
   {
      double pull_force;
      pull_force = -1.0e-2;
			pull_force = pull_force/boucle;
      f.Set(dim-1, new ConstantCoefficient(pull_force));
   }

   LinearForm *b = new LinearForm(fespace);
   b->AddBoundaryIntegrator(new VectorBoundaryLFIntegrator(f));
   cout << "r.h.s. ... " << flush;


   // 8. Define the solution vector x as a finite element grid function
   //    corresponding to fespace. Initialize x with initial guess of zero,
   //    which satisfies the boundary conditions.
   GridFunction x(fespace);
   x = 0.0;
   GridFunction x1(fespace);		// somme des résultats
   x1 = 0.0;

   // 9. Set up the bilinear form a(.,.) on the finite element space
   //    corresponding to the linear elasticity integrator with piece-wise
   //    constants coefficient lambda and mu.
   double lambda;
   lambda = 1.0;
   ConstantCoefficient lambda_func(lambda);
   double mu;
   mu = 1.0;
   ConstantCoefficient mu_func(mu);

   BilinearForm *a = new BilinearForm(fespace);
	 BilinearFormIntegrator *integ = new ElasticityIntegrator(lambda_func,mu_func);
   a->AddDomainIntegrator(new ElasticityIntegrator(lambda_func,mu_func));

   // 10. Assemble the bilinear form and the corresponding linear system,
   //     applying any necessary transformations such as: eliminating boundary
   //     conditions, applying conforming constraints for non-conforming AMR,
   //     static condensation, etc.
   const int tdim = dim*(dim+1)/2; // num. entries in a symmetric tensor
   ConstantCoefficient zero(0.0);
for(int j=0;j<boucle;j++){		//répétition du calcul

//===============Définition estimateur==========================
   FiniteElementSpace flux_fespace(mesh, fec, tdim);
   ZienkiewiczZhuEstimator estimator(*integ, x, flux_fespace);
   estimator.SetAnisotropic();

   ThresholdRefiner refiner(estimator);
   refiner.SetTotalErrorFraction(1.);


   for (int it = 0; ; it++)		// Boucle raffinement
   {
	 int cdofs = fespace->GetTrueVSize();
   cout << "matrix ... " << flush;
   if (static_cond) { a->EnableStaticCondensation(); }
   a->Assemble();
   b->Assemble();

   fespace->GetEssentialTrueDofs(ess_bdr, ess_tdof_list);

   SparseMatrix A;
   Vector B, X;


   a->FormLinearSystem(ess_tdof_list, x, *b, A, X, B);
   cout << "done." << endl;

   cout << "Size of linear system: " << A.Height() << endl;

#ifndef MFEM_USE_SUITESPARSE
   // 11. Define a simple symmetric Gauss-Seidel preconditioner and use it to
   //     solve the system Ax=b with PCG.
   GSSmoother M(A);
   PCG(A, M, B, X, 2, 5000, 1e-8, 0.0);
#else
   // 11. If MFEM was compiled with SuiteSparse, use UMFPACK to solve the system.
   UMFPackSolver umf_solver;
   umf_solver.Control[UMFPACK_ORDERING] = UMFPACK_ORDERING_METIS;
   umf_solver.SetOperator(A);
   umf_solver.Mult(B, X);
#endif

   // 12. Recover the solution as a finite element grid function after calcul.
   a->RecoverFEMSolution(X, *b, x);
      if (cdofs > max_dofs)
      {
         cout << "Reached the maximum number of dofs. Stop." << endl;
         break;
      }

      // 13. Call the refiner to modify the mesh. The refiner calls the error
      //     estimator to obtain element errors, then it selects elements to be
      //     refined and finally it modifies the mesh. The Stop() method can be
      //     used to determine if a stopping criterion was met.
     refiner.Apply(*mesh);
      if (refiner.Stop())
      {
         cout << "Stopping criterion satisfied. Stop." << endl;
         break;
      }

      // 14. Update the space to reflect the new state of the mesh. Also,
      //     interpolate the solution x so that it lies in the new space but
      //     represents the same function. This saves solver iterations later
      //     since we'll have a good initial guess of x in the next step.
	      //     Internally, FiniteElementSpace::Update() calculates an
      //     interpolation matrix which is then used by GridFunction::Update().
      fespace->Update();
      x.Update();

      // 15. Inform also the bilinear and linear forms that the space has
      //     changed.
      a->Update();
      b->Update();
}
  x1+=x;  //somme résultat
  } 
 //fin de boucle de répétition du calcul


   // 16. For non-NURBS meshes, make the mesh curved based on the finite element
   //     space. This means that we define the mesh elements through a fespace
   //     based transformation of the reference element. This allows us to save
   //     the displaced mesh as a curved mesh when using high-order finite
   //     element displacement field. We assume that the initial mesh (read from
   //     the file) is not higher order curved mesh compared to the chosen FE
   //     space.
   if (!mesh->NURBSext)
   {
      mesh->SetNodalFESpace(fespace);
   }

   // 17. Save the displaced mesh and the inverted solution (which gives the
   //     backward displacements to the original grid). This output can be
   //     viewed later using GLVis: "glvis -m displaced.mesh -g sol.gf".
   {
      GridFunction *nodes = mesh->GetNodes();
      *nodes += x;
      x *= -1;
      ofstream mesh_ofs("displaced.mesh");
      mesh_ofs.precision(8);
      mesh->Print(mesh_ofs);
      ofstream sol_ofs("sol.gf");
      sol_ofs.precision(8);
      x.Save(sol_ofs);
   }


		// save dans un fichier .txt
   {
      GridFunction *nodes = mesh->GetNodes();
      *nodes += x1;
      ofstream sol_ofs("sol_boucle_amr.txt");
      sol_ofs.precision(8);
      x1.Save(sol_ofs);
   }

   // 18. Send the above data by socket to a GLVis server. Use the "n" and "b"
   //     keys in GLVis to visualize the displacements.
   if (visualization)
   {
      char vishost[] = "localhost";
      int  visport   = 19916;
      socketstream sol_sock(vishost, visport);
      sol_sock.precision(8);
      sol_sock << "solution\n" << *mesh << x1 << flush;
   }

   // 19. Free the used memory.
   delete a;
   delete b;
   if (fec)
   {
      delete fespace;
      delete fec;
   }
   delete mesh;

   return 0;
}
