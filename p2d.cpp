//                       MFEM Example 16 - Parallel Version
//
// Compile with: make ex16p
//
// Sample runs:  mpirun -np 4 ex16p
//               mpirun -np 4 ex16p -m ../data/inline-tri.mesh
//               mpirun -np 4 ex16p -m ../data/disc-nurbs.mesh -tf 2
//               mpirun -np 4 ex16p -s 1 -a 0.0 -k 1.0
//               mpirun -np 4 ex16p -s 2 -a 1.0 -k 0.0
//               mpirun -np 8 ex16p -s 3 -a 0.5 -k 0.5 -o 4
//               mpirun -np 4 ex16p -s 14 -dt 1.0e-4 -tf 4.0e-2 -vs 40
//               mpirun -np 16 ex16p -m ../data/fichera-q2.mesh
//               mpirun -np 16 ex16p -m ../data/fichera-mixed.mesh
//               mpirun -np 16 ex16p -m ../data/escher-p2.mesh
//               mpirun -np 8 ex16p -m ../data/beam-tet.mesh -tf 10 -dt 0.1
//               mpirun -np 4 ex16p -m ../data/amr-quad.mesh -o 4 -rs 0 -rp 0
//               mpirun -np 4 ex16p -m ../data/amr-hex.mesh -o 2 -rs 0 -rp 0
//
// Description:  This example solves a time dependent nonlinear heat equation
//               problem of the form du/dt = C(u), with a non-linear diffusion
//               operator C(u) = \nabla \cdot (\kappa + \alpha u) \nabla u.
//
//               The example demonstrates the use of nonlinear operators (the
//               class EquationOperator defining C(u)), as well as their
//               implicit time integration. Note that implementing the method
//               EquationOperator::ImplicitSolve is the only requirement for
//               high-order implicit (SDIRK) time integration. In this example,
//               the diffusion operator is linearized by evaluating with the
//               lagged solution from the previous timestep, so there is only
//               a linear solve.

#include "mfem.hpp"
#include <fstream>
#include <iostream>

using namespace std;
using namespace mfem;

const real_t C0 = 0.0;
const real_t D = 1.0;

double  function1(const Vector & x){ return D * x(0) * x(0); }
double  function2(const Vector & x){ return x(0) * x(0); }

/** After spatial discretization, the conduction model can be written as:
 *
 *     du/dt = M^{-1}(-Ku)
 *
 *  where u is the vector representing the temperature, M is the mass matrix,
 *  and K is the diffusion operator with diffusivity depending on u:
 *  (\kappa + \alpha u).
 *
 *  Class EquationOperator represents the right-hand side of the above ODE.
 */
class EquationOperator : public TimeDependentOperator
{
protected:
   ParFiniteElementSpace &fespace;
   Array<int> ess_tdof_list; // this list remains empty for pure Neumann b.c.
   Array<int> nbc_bdr; // this list remains empty for pure Neumann b.c.

   ParBilinearForm *M;
   ParBilinearForm *K;
   ParLinearForm *Q;

   HypreParMatrix Mmat;
   HypreParMatrix Kmat;
   HypreParVector Qvec;

   HypreParMatrix *C; // C = M + dt K
   real_t current_dt;

   CGSolver M_solver;    // Krylov solver for inverting the mass matrix M
   HypreSmoother M_prec; // Preconditioner for the mass matrix M

   CGSolver C_solver;    // Implicit solver for T = M + dt K
   HypreSmoother C_prec; // Preconditioner for the implicit solver

   mutable Vector z; // auxiliary vector

public:
   EquationOperator(ParFiniteElementSpace &f, const Vector &u, const Array<int> &etl, const Array<int> &nb);

   virtual void Mult(const Vector &u, Vector &du_dt) const;
   /** Solve the Backward-Euler equation: k = f(u + dt*k, t), for the unknown k.
       This is the only requirement for high-order SDIRK implicit integration.*/
   virtual void ImplicitSolve(const real_t dt, const Vector &u, Vector &k);

   /// Update the diffusion BilinearForm K using the given true-dof vector `u`.
   virtual void SetParameters(const Vector &u) = 0;

   virtual ~EquationOperator() {}
};

class ParticleConcentrationOperator : public EquationOperator
{
   public:
      using EquationOperator::EquationOperator;
      virtual void SetParameters(const Vector &u);
      virtual ~ParticleConcentrationOperator();
};

class ElectrolyteConcentrationOperator : public EquationOperator
{
   public:
      using EquationOperator::EquationOperator;
      virtual void SetParameters(const Vector &u);
      virtual ~ElectrolyteConcentrationOperator();
};

int main(int argc, char *argv[])
{
   // 1. Initialize MPI and HYPRE.
   Mpi::Init(argc, argv);
   int num_procs = Mpi::WorldSize();
   int myid = Mpi::WorldRank();
   Hypre::Init();

   // 2. Parse command-line options.
   int ser_ref_levels = 0;
   int par_ref_levels = 0;
   int order = 1;
   int ode_solver_type = 3;
   real_t t_final = 1.0;
   real_t dt = 1.0e-2;
   bool visualization = true;
   int vis_steps = 5;

   int precision = 8;
   cout.precision(precision);

   OptionsParser args(argc, argv);
   args.AddOption(&ser_ref_levels, "-rs", "--refine-serial",
                  "Number of times to refine the mesh uniformly in serial.");
   args.AddOption(&par_ref_levels, "-rp", "--refine-parallel",
                  "Number of times to refine the mesh uniformly in parallel.");
   args.AddOption(&order, "-o", "--order",
                  "Order (degree) of the finite elements.");
   args.AddOption(&ode_solver_type, "-s", "--ode-solver",
                  "ODE solver: 1 - Backward Euler, 2 - SDIRK2, 3 - SDIRK3,\n\t"
                  "\t   11 - Forward Euler, 12 - RK2, 13 - RK3 SSP, 14 - RK4.");
   args.AddOption(&t_final, "-tf", "--t-final",
                  "Final time; start time is 0.");
   args.AddOption(&dt, "-dt", "--time-step",
                  "Time step.");
   args.AddOption(&visualization, "-vis", "--visualization", "-no-vis",
                  "--no-visualization",
                  "Enable or disable ParaView visualization.");
   args.AddOption(&vis_steps, "-vs", "--visualization-steps",
                  "Visualize every n-th timestep.");
   args.Parse();
   if (!args.Good())
   {
      args.PrintUsage(cout);
      return 1;
   }

   if (myid == 0)
      args.PrintOptions(cout);

   // 3. Read the serial mesh from the given mesh file on all processors. We can
   //    handle triangular, quadrilateral, tetrahedral and hexahedral meshes
   //    with the same code.
   Mesh serial_mesh = Mesh::MakeCartesian1D(10);
   int dim = serial_mesh.Dimension();

   // 4. Define the ODE solver used for time integration. Several implicit
   //    singly diagonal implicit Runge-Kutta (SDIRK) methods, as well as
   //    explicit Runge-Kutta methods are available.
   ODESolver *ode_solver;
   switch (ode_solver_type)
   {
      // Implicit L-stable methods
      case 1:  ode_solver = new BackwardEulerSolver; break;
      case 2:  ode_solver = new SDIRK23Solver(2); break;
      case 3:  ode_solver = new SDIRK33Solver; break;
      // Explicit methods
      case 11: ode_solver = new ForwardEulerSolver; break;
      case 12: ode_solver = new RK2Solver(0.5); break; // midpoint method
      case 13: ode_solver = new RK3SSPSolver; break;
      case 14: ode_solver = new RK4Solver; break;
      case 15: ode_solver = new GeneralizedAlphaSolver(0.5); break;
      // Implicit A-stable methods (not L-stable)
      case 22: ode_solver = new ImplicitMidpointSolver; break;
      case 23: ode_solver = new SDIRK23Solver; break;
      case 24: ode_solver = new SDIRK34Solver; break;
      default:
         cout << "Unknown ODE solver type: " << ode_solver_type << '\n';
         return 1;
   }

   // 5. Refine the mesh in serial to increase the resolution. In this example
   //    we do 'ser_ref_levels' of uniform refinement, where 'ser_ref_levels' is
   //    a command-line parameter.
   for (int lev = 0; lev < ser_ref_levels; lev++)
      serial_mesh.UniformRefinement();

   // 6. Define a parallel mesh by a partitioning of the serial mesh. Refine
   //    this mesh further in parallel to increase the resolution. Once the
   //    parallel mesh is defined, the serial mesh can be deleted.
   ParMesh *pmesh = new ParMesh(MPI_COMM_WORLD, serial_mesh);
   serial_mesh.Clear(); // the serial mesh is no longer needed
   for (int lev = 0; lev < par_ref_levels; lev++)
      pmesh->UniformRefinement();

   // 7. Define the vector finite element space representing the current and the
   //    initial temperature, u_ref.
   H1_FECollection fe_coll(order, dim);
   ParFiniteElementSpace fespace(pmesh, &fe_coll);

   HYPRE_BigInt fe_size_global = fespace.GlobalTrueVSize();
   if (myid == 0)
      cout << "Unknowns (total): " << fe_size_global << endl;

   HYPRE_BigInt fe_size_owned = fespace.TrueVSize();
   cout << "Unknowns (rank " << myid << "): " << fe_size_owned << endl;

   ParGridFunction u_gf(&fespace);

   // 8. Set the initial/boundary conditions for u.
   Array<int> ess_tdof_list;
   Array<int> ess_bdr(pmesh->bdr_attributes.Max());
   ess_bdr = 0;
   fespace.GetEssentialTrueDofs(ess_bdr, ess_tdof_list);

   Array<int> nbc_bdr(pmesh->bdr_attributes.Max());
   nbc_bdr = 0; nbc_bdr[1] = 1;

   ConstantCoefficient u_0(C0);
   u_gf.ProjectCoefficient(u_0);
   ConstantCoefficient u_b(C0);
   u_gf.ProjectBdrCoefficient(u_b, ess_bdr);
   Vector u;
   u_gf.GetTrueDofs(u);

   // 9. Initialize the conduction operator and the VisIt visualization.
   ParticleConcentrationOperator oper(fespace, u, ess_tdof_list, nbc_bdr);

   //u_gf.SetFromTrueDofs(u);
   ParaViewDataCollection pd("particle", pmesh);
   pd.SetPrefixPath("ParaView");
   pd.SetLevelsOfDetail(order);
   pd.SetHighOrderOutput(true);
   pd.SetDataFormat(VTKFormat::BINARY);
   pd.RegisterField("u", &u_gf);
   pd.SetCycle(0);
   pd.SetTime(0.0);
   pd.Save();

   // 10. Perform time-integration (looping over the time iterations, ti, with a
   //     time-step dt).
   ode_solver->Init(oper);
   real_t t = 0.0;

   oper.SetParameters(u);

   bool last_step = false;
   for (int ti = 1; !last_step; ti++)
   {
      last_step = t + dt >= t_final - dt/2;

      ode_solver->Step(u, t, dt);

      if (last_step || (ti % vis_steps) == 0)
      {
         if (myid == 0)
            cout << "step " << ti << ", t = " << t << endl;

         u_gf.SetFromTrueDofs(u);

         LinearForm sum(&fespace);
         GridFunctionCoefficient u_gfc(&u_gf);
         FunctionCoefficient r2(function1);
         ProductCoefficient ur2(u_gfc,r2);
         sum.AddDomainIntegrator(new DomainLFIntegrator(ur2));
         sum.Assemble();
   
         std::cout << "To total flux accumulated = " << sum.Sum() << std::endl;

         if (last_step || visualization)
         {
            pd.SetCycle(ti);
            pd.SetTime(t);
            pd.Save();
         }
      }
      oper.SetParameters(u);
   }

   // 11. Free the used memory.
   delete ode_solver;
   delete pmesh;

   return 0;
}

EquationOperator::EquationOperator(ParFiniteElementSpace &f, const Vector &u, const Array<int> &etl, const Array<int> &nb)
   : TimeDependentOperator(f.GetTrueVSize(), (real_t) 0.0), fespace(f),
     ess_tdof_list(etl), nbc_bdr(nb), M(NULL), K(NULL), Q(NULL), C(NULL), current_dt(0.0),
     M_solver(f.GetComm()), C_solver(f.GetComm()), z(height)
{
   const real_t rel_tol = 1e-8;

   M_solver.iterative_mode = false;
   M_solver.SetRelTol(rel_tol);
   M_solver.SetAbsTol(0.0);
   M_solver.SetMaxIter(100);
   M_solver.SetPrintLevel(0);
   M_prec.SetType(HypreSmoother::Jacobi);
   M_solver.SetPreconditioner(M_prec);
   M_solver.SetOperator(Mmat);

   C_solver.iterative_mode = false;
   C_solver.SetRelTol(rel_tol);
   C_solver.SetAbsTol(0.0);
   C_solver.SetMaxIter(100);
   C_solver.SetPrintLevel(0);
   C_solver.SetPreconditioner(C_prec);
}

void EquationOperator::Mult(const Vector &u, Vector &du_dt) const
{
   // Compute:
   //    du_dt = M^{-1}*-Ku
   // for du_dt, where K is linearized by using u from the previous timestep
   Kmat.Mult(u, z);
   z.Neg(); // z = -z
   z += Qvec;
   HypreParMatrix A; Vector X, Z;
   ParGridFunction uu(&fespace), zz(&fespace);
   uu.SetFromTrueDofs(u);
   const SparseMatrix &R = *(fespace.GetRestrictionMatrix());
   R.MultTranspose(z, zz);
   K->FormLinearSystem(ess_tdof_list, uu, zz, A, X, Z);
   M_solver.Mult(Z, du_dt);
   du_dt.SetSubVector(ess_tdof_list, 0.0);
}

void EquationOperator::ImplicitSolve(const real_t dt,
                                       const Vector &u, Vector &du_dt)
{
   // Solve the equation:
   //    du_dt = M^{-1}*[-K(u + dt*du_dt)]
   // for du_dt, where K is linearized by using u from the previous timestep
   if (!C)
   {
      C = Add(1.0, Mmat, dt, Kmat);
      current_dt = dt;
      C_solver.SetOperator(*C);
   }
   MFEM_VERIFY(dt == current_dt, ""); // SDIRK methods use the same dt
   Kmat.Mult(u, z);
   z.Neg();
   z += Qvec;
   HypreParMatrix A; Vector X, Z;
   ParGridFunction uu(&fespace), zz(&fespace);
   uu.SetFromTrueDofs(u);
   const SparseMatrix &R = *(fespace.GetRestrictionMatrix());
   R.MultTranspose(z, zz);
   K->FormLinearSystem(ess_tdof_list, uu, zz, A, X, Z);
   C_solver.Mult(Z, du_dt);
   du_dt.SetSubVector(ess_tdof_list, 0.0);
}

void ParticleConcentrationOperator::SetParameters(const Vector &u)
{
   ParGridFunction u_gf(&fespace);
   u_gf.SetFromTrueDofs(u);

   IntegrationRule ir = IntRules.Get(Geometry::SEGMENT, 6);

   //std::cout << "Number of points " << ir.GetNPoints() << std::endl;

   FunctionCoefficient diff(function1);
   FunctionCoefficient coeff(function2);
   ConstantCoefficient nbcCoef(1.0);

   delete M;
   M = new ParBilinearForm(&fespace);
   M->AddDomainIntegrator(new MassIntegrator(coeff, &ir));
   M->Assemble(0); // keep sparsity pattern of M and K the same
   M->FormSystemMatrix(ess_tdof_list, Mmat);

   delete K;
   K = new ParBilinearForm(&fespace);
   K->AddDomainIntegrator(new DiffusionIntegrator(diff, &ir));
   K->Assemble(0); // keep sparsity pattern of M and K the same
   K->FormSystemMatrix(ess_tdof_list, Kmat);

   delete Q;
   Q = new ParLinearForm(&fespace);
   Q->AddBoundaryIntegrator(new BoundaryLFIntegrator(nbcCoef), nbc_bdr);
   Q->Assemble();
   Qvec = std::move(*(Q->ParallelAssemble()));
   Qvec.SetSubVector(ess_tdof_list, 0.0); // do we need this?
   
   delete C;
   C = NULL; // re-compute C on the next ImplicitSolve
}

ParticleConcentrationOperator::~ParticleConcentrationOperator()
{
   delete C;
   delete M;
   delete K;
   delete Q;
}


void ElectrolyteConcentrationOperator::SetParameters(const Vector &u)
{
   ParGridFunction u_gf(&fespace);
   u_gf.SetFromTrueDofs(u);

   real_t a, tplus;
   real_t temp = (1 - tplus) * a; 

   IntegrationRule ir = IntRules.Get(Geometry::SEGMENT, 6);

   //std::cout << "Number of points " << ir.GetNPoints() << std::endl;

   ConstantCoefficient one(1.0);
   ConstantCoefficient nbcCoef(1.0);
   ConstantCoefficient Coeff(temp);

   delete M;
   M = new ParBilinearForm(&fespace);
   M->AddDomainIntegrator(new MassIntegrator(one, &ir));
   M->Assemble(0); // keep sparsity pattern of M and K the same
   M->FormSystemMatrix(ess_tdof_list, Mmat);

   delete K;
   K = new ParBilinearForm(&fespace);
   K->AddDomainIntegrator(new DiffusionIntegrator(one, &ir));
   K->Assemble(0); // keep sparsity pattern of M and K the same
   K->FormSystemMatrix(ess_tdof_list, Kmat);

   delete Q;
   Q = new ParLinearForm(&fespace);
   Q->AddDomainIntegrator(new DomainLFIntegrator(Coeff));
   Q->Assemble();
   Qvec = std::move(*(Q->ParallelAssemble()));
   Qvec.SetSubVector(ess_tdof_list, 0.0); // do we need this?
   
   delete C;
   C = NULL; // re-compute C on the next ImplicitSolve
}

ElectrolyteConcentrationOperator::~ElectrolyteConcentrationOperator()
{
   delete C;
   delete M;
   delete K;
   delete Q;
}