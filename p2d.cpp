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
//               problem of the form du/dt = C(x), with a non-linear diffusion
//               operator C(x) = \nabla \cdot (\kappa + \alpha x) \nabla x.
//
//               The example demonstrates the use of nonlinear operators (the
//               class EquationOperator defining C(x)), as well as their
//               implicit time integration. Note that implementing the method
//               EquationOperator::ImplicitSolve is the only requirement for
//               high-order implicit (SDIRK) time integration. In this example,
//               the diffusion operator is linearized by evaluating with the
//               lagged solution from the previous timestep, so there is only
//               a linear solve.

#include "mfem.hpp"
#include <fstream>
#include <iostream>
#include "P2DOperator.hpp"
#include <cmath>

using namespace std;
using namespace mfem;

int main(int argc, char *argv[])
{
   // 1. Initialize MPI and HYPRE.
   Mpi::Init(argc, argv);
   int num_procs = Mpi::WorldSize();
   int myid = Mpi::WorldRank();
   Hypre::Init();

   // 2. Parse command-line options.
   Method method = SPM;
   int ser_ref_levels = 0;
   int par_ref_levels = 0;
   int order = 1;
   int ode_solver_type = 1;
   real_t t_final = 1.0;
   real_t dt = 1.0/3600.0;
   bool visualization = true;
   int vis_steps = 5;

   int precision = 8;
   cout.precision(precision);

   OptionsParser args(argc, argv);
   args.AddOption(reinterpret_cast<int*>(&method), "-m", "--method",
                  "Electrochemical method: 0) SPM, 1) SPMe, 2) P2D.");
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

   if (Mpi::Root())
      args.PrintOptions(cout);

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

   // Initialise grid and layout properties dependent on the EC method and FE order
   init_params(method, order);

   // 3. Read the serial mesh from the given mesh file on all processors. We can
   //    handle triangular, quadrilateral, tetrahedral and hexahedral meshes
   //    with the same code.
   Mesh x_smesh = Mesh::MakeCartesian1D(NX);
 
   x_smesh.EnsureNodes();
   GridFunction *vert_coord = x_smesh.GetNodes();

   int numnodes = (*vert_coord).Size();

   MFEM_ASSERT(numnodes == order * NX + 1, "Total number of nodes in x_mesh is wrong");

   //Vector jubat_mesh({0.0, 0.049305555555555554, 0.09861111111111111, 0.14791666666666664, 0.19722222222222222, 0.24652777777777776, 0.2958333333333333, 0.34513888888888883, 0.39444444444444443, 0.44375, 0.4930555555555555, 0.5, 0.5069444444444444, 0.5138888888888888, 0.5208333333333333, 0.5277777777777778, 0.5347222222222222, 0.5416666666666666, 0.5486111111111112, 0.5555555555555556, 0.5625, 0.60625, 0.65, 0.69375, 0.7375, 0.78125, 0.825, 0.8687499999999999, 0.9125000000000001, 0.95625, 1.0});
   /*Vector jubat_mesh({0.0,    0.0437,    0.0875,    0.1313,    0.1750,    0.2188,    0.2625,    0.3063,    0.3500,   0.3938,    0.4375,    0.4444,    0.4514,    0.4583,    0.4653,    0.4722,   0.4792,    0.4861,    0.4931,    0.5000,    0.5069,    0.5562,    0.6056,    0.6549,    0.7042,    0.7535,    0.8028,    0.8521,    0.9014,    0.9507,    1.0000});
   for (int i = 0; i < numnodes; i++)
   {
      std::cout << "jubat_mesh(" << i << ") = " << jubat_mesh(i) << std::endl;
      (*vert_coord)(i) = jubat_mesh(i); // RHS has to be changed with right coordinatess
   }*/

   for (unsigned i = 0; i < NX; i++)
      x_smesh.SetAttribute(i, i < NPE ? PE : i < NPE + NSEP ? SEP : NE);
   x_smesh.SetAttributes();

   Mesh r_smesh[NPAR];
   for (unsigned p = 0; p < NPAR; p++)
      r_smesh[p] = Mesh::MakeCartesian1D(NR);

   // 6. Define a parallel mesh by a partitioning of the serial mesh. Refine
   //    this mesh further in parallel to increase the resolution. Once the
   //    parallel mesh is defined, the serial mesh can be deleted.
   ParMesh *x_pmesh = new ParMesh(MPI_COMM_WORLD, x_smesh);
   x_smesh.Clear(); // the serial mesh is no longer needed
   ParMesh *r_pmesh[NPAR];
   for (unsigned p = 0; p < NPAR; p++)
   {
      r_pmesh[p] = new ParMesh(MPI_COMM_WORLD, r_smesh[p]);
      r_smesh[p].Clear(); // the serial mesh is no longer needed
   }

   // 7. Define the vector finite element space representing the current and the
   //    initial temperature, u_ref.
   H1_FECollection fe_coll(order, /*dim*/ 1);
   ParFiniteElementSpace * x_fespace = new ParFiniteElementSpace(x_pmesh, &fe_coll);
   Array<ParFiniteElementSpace *> r_fespace(NPAR);
   for (unsigned p = 0; p < NPAR; p++)
      r_fespace[p] = new ParFiniteElementSpace(r_pmesh[p], &fe_coll);

   // 8. Get the total number of dofs in the system (including boundaries)
   {
      HYPRE_BigInt fe_size_global = SC * x_fespace->GlobalTrueVSize();
      for (unsigned p = 0; p < NPAR; p++)
         fe_size_global += r_fespace[p]->GlobalTrueVSize();

      if (Mpi::Root())
         cout << "Unknowns (total): " << fe_size_global << endl;
   }

   // 8.5 Get the total number of dofs _owned_ by this processor
   HYPRE_BigInt fe_size_owned = SC * x_fespace->GetTrueVSize();
   for (unsigned p = 0; p < NPAR; p++)
      fe_size_owned += r_fespace[p]->GetTrueVSize();

   // 9. Initialize the conduction operator and the VisIt visualization.
   BlockVector x;
   P2DOperator oper(x_fespace, r_fespace, fe_size_owned, x);

   // X. Visualization for the 0th particle
   ParGridFunction u_gf(r_fespace[0]);
   ParaViewDataCollection pd("particle", r_pmesh[0]);
   pd.SetPrefixPath("ParaView");
   pd.SetLevelsOfDetail(order);
   pd.SetHighOrderOutput(true);
   pd.SetDataFormat(VTKFormat::BINARY);
   pd.RegisterField("x", &u_gf);
   pd.SetCycle(0);
   pd.SetTime(0.0);
   pd.Save();

   // Visualisation for the electrolyte
   ParGridFunction ec_gf(x_fespace);
   ParaViewDataCollection pd_e("electrolyte", x_pmesh);
   pd_e.SetPrefixPath("ParaView");
   pd_e.SetLevelsOfDetail(order);
   pd_e.SetHighOrderOutput(true);
   pd_e.SetDataFormat(VTKFormat::BINARY);
   pd_e.RegisterField("ec", &ec_gf);
   pd_e.SetCycle(0);
   pd_e.SetTime(0.0);
   pd_e.Save();

   // 10. Perform time-integration (looping over the time iterations, ti, with a
   //     time-step dt).
   ode_solver->Init(oper);
   real_t t = 0.0;

   oper.Update(x, dt);

   bool last_step = false;
   for (int ti = 1; !last_step; ti++)
   {
      std::cout << "---------------- Time step " << ti << ", t = " << t
                << " ----------------" << std::endl;
      last_step = t + dt >= t_final - dt/2;

      ode_solver->Step(x, t, dt);
      oper.ComputeVoltage(x, t, dt);

      if (last_step || (ti % vis_steps) == 0)
      {
         if (Mpi::Root())
            cout << "step " << ti << ", t = " << t << endl;

         // TODO: Stop sim at cutoff voltage
         if (last_step || visualization)
         {
            u_gf.SetFromTrueDofs(x.GetBlock(SC + 0));
            ParLinearForm sum(r_fespace[0]);
            GridFunctionCoefficient u_gfc(&u_gf);
            FunctionCoefficient r2([](const Vector & x){ return x(0) * x(0); });
            ProductCoefficient ur2(u_gfc,r2);
            sum.AddDomainIntegrator(new DomainLFIntegrator(ur2));
            sum.Assemble();
            std::cout << "[Rank " << Mpi::WorldRank() << "]"
                      << " Total flux accumulated (" << 0 << ") = " << sum.Sum() << std::endl;

            ec_gf.SetFromTrueDofs(x.GetBlock(EC));

            //std::cout << "ec(0) = " << ec_gf(0) << std::endl;

            pd.SetCycle(ti);
            pd.SetTime(t);
            pd.Save();

            pd_e.SetCycle(ti);
            pd_e.SetTime(t);
            pd_e.Save();
         }
      }
      oper.Update(x, dt);
   }

   // 11. Free the used memory.
   delete ode_solver;
   delete x_pmesh;
   for (unsigned p = 0; p < NPAR; p++)
      delete r_pmesh[p];

   return 0;
}
