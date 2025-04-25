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
#include "EquationOperator.hpp"
#include "ParticleConcentrationOperator.hpp"

using namespace std;
using namespace mfem;

const real_t T0 = 298.15;
const real_t C0 = 0.0;

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

