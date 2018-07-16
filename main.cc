//LIC//======================================================================
///LIC// This file forms part of oomph-lib,the object-oriented,
//LIC// multi-physics finite-element library,available
//LIC// at http://www.oomph-lib.org.
//LIC//
//LIC//           Version 0.90. August 3,2009.
//LIC//
//LIC// Copyright (C) 2006-2009 Matthias Heil and Andrew Hazel
//LIC//
//LIC// This library is free software; you can redistribute it and/or
//LIC// modify it under the terms of the GNU Lesser General Public
//LIC// License as published by the Free Software Foundation; either
//LIC// version 2.1 of the License,or (at your option) any later version.
//LIC//
//LIC// This library is distributed in the hope that it will be useful,
//LIC// but WITHOUT ANY WARRANTY; without even the implied warranty of
//LIC// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
//LIC// Lesser General Public License for more details.
//LIC//
//LIC// You should have received a copy of the GNU Lesser General Public
//LIC// License along with this library; if not,write to the Free Software
//LIC// Foundation,Inc.,51 Franklin Street,Fifth Floor,Boston,MA
//LIC// 02110-1301  USA.
//LIC//
//LIC// The authors may be contacted at oomph-lib@maths.man.ac.uk.
//LIC//
//LIC//
//LIC//======================================================================
// Helmoholtz equation with one PML and forcing on the opposite side and
// periodic BC on all other sides
#include <fenv.h>
#include "math.h"
#include <complex>

// Generic routines
#include "generic.h"

#include "oomph_crbond_bessel.h"

// The Helmholtz equations and complex-valued multigrid machinery
#include "SourceFiles/pml_helmholtz.h"

// The mesh
#include "meshes.h"


using namespace std;
using namespace oomph;

//////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////

//================================================start_of_namespace======
/// Namespace for the Helmholtz problem parameters
//========================================================================
namespace GlobalParameters
{
 /// \short Solver specific parameters:
 ///----------------------------------- 
 /// The number of nodes in one direction (default=2)
 unsigned Nnode_1d=2;
 
 /// The minimum level of uniform refinement
 unsigned Min_refinement_level=1;

 /// The additional levels of uniform refinement 
 unsigned Add_refinement_level=0;
 
 /// The number of adaptations allowed by the Newton solver
 unsigned N_adaptations=1;

 /// \short The choice of whether or not to use adaptation
 ///    0 = Uniform refinement
 ///    1 = Adaptive refinement
 unsigned Use_adaptation_flag=0;

 /// \short The choice of pre-smoother:
 ///    0 = Automatic (GMRES as a smoother on levels where kh>0.5)
 ///    1 = Damped Jacobi on all levels with a constant omega value
 unsigned Pre_smoother_flag=0;
 
 /// \short The choice of post-smoother:
 ///    0 = Automatic (GMRES as a smoother on levels where kh>0.5)
 ///    1 = Damped Jacobi on all levels with a constant omega value
 unsigned Post_smoother_flag=0;

 /// \short The choice of linear solver
 ///    0 = SuperLU
 ///    1 = Multigrid
 unsigned Linear_solver_flag=0;

  /// \short The MG solver allows for five different levels of output:
 ///    0 = Outputs everything
 ///    1 = Outputs everything except the smoother timings 
 ///    2 = Outputs setup information but no V-cycle timings
 ///    3 = Suppresses all output
 unsigned Output_management_flag=0;
  
 /// \short Variable used to decide whether or not convergence information
 /// is displayed:
 ///    0 = Don't display convergence information
 ///    1 = Display convergence information
 unsigned Doc_convergence_flag=0;
 
 /// DocInfo object used for documentation of the solution
 DocInfo Doc_info;
 
 // Pointer to the output stream -- defaults to oomph_info
 std::ostream* Stream_pt;
  
 /// \short Problem specific parameters:
 ///------------------------------------
 /// Dimensions of the mesh
 
 // Mesh lengths in each direction
 double Lx=1.0;
 double Ly=1.0;

 // x coordinates
 double Xmin = 1.0;
 double Xmax = Xmin + Lx;

 // y coordinates
 double Ymin = 1.0;
 double Ymax = Ymin + Ly;

 /// Number of elements in each direction (used by SimpleCubicMesh)
 unsigned Nx=27;
 unsigned Ny=27;
 
 /// Store the value of Pi
 double Pi=MathematicalConstants::Pi;

 /// Choose the value of the shift to create the complex-shifted
 /// Laplacian preconditioner (CSLP)
 double Alpha_shift=0.0;
 
 /// Square of the wavenumber (also known as k^2)
 double K_squared=1.0;

 /// Wavenumber (also known as k),k=omega/c
 double Wavenumber=sqrt(K_squared);

 /// Update the parameters passed in at the command line
 void update_parameters()
 {
  /// Wavenumber (also known as k), k=omega/c
  Wavenumber=sqrt(K_squared);
 }

 /// Fourier wavenumber
 // If N=0, then we are solving the standard Helmholtz problem
 int N_fourier_wavenumber=4;

 double L_squared = 0.0;
 double L = sqrt(L_squared);

 /// Exact solution as a Vector of size 2, containing real and imag parts
 void get_exact_u(const Vector<double>& x, Vector<double>& u)
 {
   // Convert to cylindrical coordinates
   double r = x[0];
   double z = x[1];

   // Argument for Bessel function, scale radial distance by wavenumber
   double r_times_ksq_lsq = r*sqrt(K_squared+L_squared);

   // Calculate one more term than the Fourier wavenumber
   int n_terms = N_fourier_wavenumber + 1;

   Vector<double> jv(n_terms);
   Vector<double> djv(n_terms);
   Vector<double> yv(n_terms);
   Vector<double> dyv(n_terms);
   
   double n_actual = 0;
   CRBond_Bessel::bessjyv(n_terms, 
                          r_times_ksq_lsq,
                          n_actual,
                          &jv[0],&yv[0],
                          &djv[0],&dyv[0]);

  // double j0,j1,y0,y1,j0p,j1p,y0p,y1p;
  // CRBond_Bessel::bessjy01a(ksq_lsq*r,j0,j1,y0,y1,j0p,j1p,y0p,y1p);
  
  complex<double> u_ex(0.0, 0.0);
  u_ex+= jv[N_fourier_wavenumber]*exp(L*z);
  
  u[0]=u_ex.real();
  u[1]=u_ex.imag();
  std::cout << "Do we even get here?" << std::endl;
 
  // u[0] = 0.0;
  // u[1] = 0.0;
 }
 
 FiniteElement::SteadyExactSolutionFctPt exact_u_pt=&get_exact_u;

} // End of namespace

/////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////

//======start_of_namespace================================================
/// Returns a pointer to a smoother of the appropriate type
//========================================================================
namespace Smoother_Factory_Function_Helper
{
 /// The value of the damping factor for the damped Jacobi smoother
 double Omega=0.4;
 
 /// \short Returns a pointer to a Smoother object which is to be used as
 /// the pre-smoother
 HelmholtzSmoother* set_pre_smoother()
 {
  // Create a new DampedJacobi object
  return new ComplexDampedJacobi<CRDoubleMatrix>(Omega);
 } 
 
 /// \short Returns a pointer to a Smoother object which is to be used as
 /// the post-smoother
 HelmholtzSmoother* set_post_smoother()
 {
  // Create a new DampedJacobi object
  return new ComplexDampedJacobi<CRDoubleMatrix>(Omega);
 }
} // End of Smoother_Factory_Function_Helper

//////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////

//============================================start_of_problem_class======
/// Problem class
//========================================================================
template<class ELEMENT>
class PMLFourierDecomposedHelmholtzProblem : public HelmholtzMGProblem
{

public:

 /// Constructor
 PMLFourierDecomposedHelmholtzProblem();

 /// Destructor (empty)
 ~PMLFourierDecomposedHelmholtzProblem();

 /// Doc the solution
 void doc_solution();

 /// Update the problem specs before solve (empty)
 void actions_before_newton_solve(){}

 /// Update the problem specs after solve (empty)
 void actions_after_newton_solve()
 {
  // Document the solution
  doc_solution();
 }

 /// Actions before adapt: (empty)
 void actions_before_adapt(){}

 /// Actions after adapt:(empty)
 void actions_after_adapt();

 /// Set GMRES preconditioner by multigrid as the linear solver
 void set_gmres_multigrid_solver();
 
 /// Enable the PML mapping function for all nodes in the PML region
 void enable_pmls();

 // Apply boundary conditions
 void apply_boundary_conditions();

private:

 /// Pointer to the "bulk" mesh
 // RefineableRectangularQuadMesh<ELEMENT>* Bulk_mesh_pt;
 RectangularQuadMesh<ELEMENT>* Bulk_mesh_pt;

 /// Overload the make_new_problem function to return an object of this class
 HelmholtzMGProblem* make_new_problem()
 {
  // Return a new problem pointer
  return new PMLFourierDecomposedHelmholtzProblem<ELEMENT>;
 }

 /// \short Overload the mg_bulk_mesh_pt function to return a pointer to the
 /// "refineable" portion of the mesh
 TreeBasedRefineableMeshBase* mg_bulk_mesh_pt()
 {
  // Return the pointer to the bulk mesh
  // return Bulk_mesh_pt;
  return nullptr;
 }
 
 /// Trace file
 ofstream Trace_file;
}; // End of PMLFourierDecomposedHelmholtzProblem class

//==============================================start_of_constructor======
/// Constructor for Helmholtz problem
//========================================================================
template<class ELEMENT>
PMLFourierDecomposedHelmholtzProblem<ELEMENT>::PMLFourierDecomposedHelmholtzProblem()
{
 // Indicate that the problem is nonlinear to ensure the residual is
 // calculated at the end of the iteration
 problem_is_nonlinear(true);

 // Set the number of Newton iterations to one
 max_newton_iterations()=10;

 // Set up solver specific information:
 //------------------------------------
 // If we're choosing to use GMRES & MG as our linear solver
 if (GlobalParameters::Linear_solver_flag==1)
 {
  // Set the solver
  set_gmres_multigrid_solver();
 }
 
 // Open trace file
 Trace_file.open("RESLT/trace.dat");

 // Build the mesh using the specified parameters:
 //-----------------------------------------------
 // Build the "bulk" mesh
 Bulk_mesh_pt=new RectangularQuadMesh<ELEMENT>(
  GlobalParameters::Nx,GlobalParameters::Ny,
  GlobalParameters::Xmin,GlobalParameters::Xmax,
  GlobalParameters::Ymin,GlobalParameters::Ymax);

 // Create the main mesh
 add_sub_mesh(Bulk_mesh_pt);

 // Build the entire mesh from its submeshes
 build_global_mesh();

 // Complete the build of all elements so they are fully functional: 
 //-----------------------------------------------------------------
 // How many elements in the mesh?
 unsigned n_element=mesh_pt()->nelement();

 // Loop over the elements and pass a pointer to the value of k^2
 for (unsigned e=0;e<n_element;e++)
 {
  // Upcast from GeneralisedElement to Helmholtz bulk element
  ELEMENT* el_pt=dynamic_cast<ELEMENT*>(Bulk_mesh_pt->element_pt(e));

  // If the upcast was successful
  if (el_pt!=0)
  {
   // Set the wavenumber function pointer
   el_pt->k_squared_pt()=&GlobalParameters::K_squared;

   // Set the Fourier wavenumber
   el_pt->n_fourier_wavenumber_pt()=&GlobalParameters::N_fourier_wavenumber;
   
  } // if (el_pt!=0)
 } // for (unsigned e=0;e<n_element;e++)

 // Apply the boundary conditions, both in the central region and on the
 // outer boundary (since these nodes are PML nodes)
 apply_boundary_conditions();
 
 // Setup equation numbering scheme
 assign_eqn_numbers(); 
} // End of constructor

//===============================================start_of_destructor======
/// Destructor for Helmholtz problem
//========================================================================
template<class ELEMENT>
PMLFourierDecomposedHelmholtzProblem<ELEMENT>::~PMLFourierDecomposedHelmholtzProblem()
{   

 // If we're using GMRES & MG as the linear solver
 if (GlobalParameters::Linear_solver_flag==1)
 {
  // Delete the MG solver pointers
  delete dynamic_cast<HelmholtzFGMRESMG<CRDoubleMatrix>* >
   (linear_solver_pt())->preconditioner_pt();

  // Set the pointer to null
  dynamic_cast<HelmholtzFGMRESMG<CRDoubleMatrix>* >
   (linear_solver_pt())->preconditioner_pt()=0;
    
  // Delete the MG solver pointers
  delete linear_solver_pt();

  // Set the pointer to null
  linear_solver_pt()=0;    
 }
   
 // Delete the error estimator
 // delete Bulk_mesh_pt->spatial_error_estimator_pt();
 
 // Set the pointer to null
 // Bulk_mesh_pt->spatial_error_estimator_pt()=0;
 // Delete the "bulk" mesh
 delete Bulk_mesh_pt;

 // Set the pointer to null
 Bulk_mesh_pt=0;
} // End of ~PMLFourierDecomposedHelmholtzProblem

//================================start_of_apply_boundary_conditions======
/// Apply boundary conditions
//========================================================================
template<class ELEMENT>
void PMLFourierDecomposedHelmholtzProblem<ELEMENT>::apply_boundary_conditions()
{
 // Find the number of elements in the mesh
 unsigned n_element=Bulk_mesh_pt->nelement();

 // Vector to hold the local coordinates of a point in any given element 
 Vector<double> s(2,0.0);
 
 // Vector to hold the (Eulerian) coordinates of a point
 // Vector<double> x(3,0.0);
 Vector<double> x(2,0.0);

 // Vector to hold the real and imaginary part of the solution
 Vector<double> u(2,0.0);
    
 // Find the number of boundaries in the mesh
 unsigned n_bound=Bulk_mesh_pt->nboundary();
 
 // Loop over all boundaries
 for (unsigned b=0;b<n_bound;b++)
 {
  // Find the number of nodes on the b-th boundary
  unsigned n_node=Bulk_mesh_pt->nboundary_node(b);

  // Loop over the nodes on the b-th boundary
  for(unsigned n=0;n<n_node;n++)
  {
   // All of these nodes sides are PMLs so pin to 0
   Node* boundary_node_pt=Bulk_mesh_pt->boundary_node_pt(b,n);

   // Pin the (real) dof at this node
   boundary_node_pt->pin(0);
   
   // Pin the (imaginary) dof at this node
   boundary_node_pt->pin(1);

   GlobalParameters::get_exact_u(x, u);

   // Set the solution value at this point (real part)
   boundary_node_pt->set_value(0,u[0]);
   
   // Set the solution value at this point (imaginary part)
   boundary_node_pt->set_value(1,u[1]);
  }
 } // for(unsigned b=0;b<n_bound;b++)
} // End of apply_boundary_conditions

//=======set_gmres_multigrid_solver=======================================
/// Build and set GMRES preconditioner by multigrid as the linear solver
//========================================================================
template<class ELEMENT>
void PMLFourierDecomposedHelmholtzProblem<ELEMENT>::set_gmres_multigrid_solver()
{
 // Create linear solver
 HelmholtzFGMRESMG<CRDoubleMatrix>* solver_pt=
  new HelmholtzFGMRESMG<CRDoubleMatrix>;

 // Set the number of iterations
 solver_pt->max_iter()=200;

 // Set the tolerance (to ensure the Newton solver converges in one step)
 solver_pt->tolerance()=1.0e-10;
   
 // If the user wishes to document the convergence information
 if (GlobalParameters::Doc_convergence_flag)
 {
  // Create a file to record the convergence history
  solver_pt->open_convergence_history_file_stream("RESLT/conv.dat");
 }
 
 // Create linear solver 
 linear_solver_pt()=solver_pt;
 
 // This preconditioner uses multigrid on the block version of the full
 // matrix. 2 V-cycles will be used here per preconditioning step
 HelmholtzMGPreconditioner<2>* prec_pt=new HelmholtzMGPreconditioner<2>(this);

 // Set preconditioner
 solver_pt->preconditioner_pt()=prec_pt;
  
 // Set the shift
 prec_pt->alpha_shift()=GlobalParameters::Alpha_shift;
   
 // If the user wants to use damped Jacobi on every level as a smoother
 if (GlobalParameters::Pre_smoother_flag==1)
 {
  // Set the pre-smoother factory function
  prec_pt->set_pre_smoother_factory_function
   (Smoother_Factory_Function_Helper::set_pre_smoother);
 }

 // If the user wants to use damped Jacobi on every level as a smoother
 if (GlobalParameters::Post_smoother_flag==1)
 {
  // Set the post-smoother factory function
  prec_pt->set_post_smoother_factory_function
   (Smoother_Factory_Function_Helper::set_post_smoother);
 }
 
 // Suppress certain timings
 if (GlobalParameters::Output_management_flag==1)
 {
  prec_pt->disable_doc_time();
 }
 else if (GlobalParameters::Output_management_flag==2)
 {
  prec_pt->disable_v_cycle_output();
 }
 else if (GlobalParameters::Output_management_flag==3)
 {
  prec_pt->disable_output();
 }
} // End of set_gmres_multigrid_solver

//======================================start_of_actions_after_adapt======
/// Actions after adapt: Re-apply the boundary conditions
//========================================================================
template<class ELEMENT>
void PMLFourierDecomposedHelmholtzProblem<ELEMENT>::actions_after_adapt()
{
 // Complete the build of all elements so they are fully functional: 
 //-----------------------------------------------------------------
 // How many elements in the mesh?
 unsigned n_element=mesh_pt()->nelement();

 // Loop over the elements and pass a pointer to the value of k^2
 for (unsigned e=0;e<n_element;e++)
 {
  // Upcast from GeneralisedElement to Helmholtz bulk element
  ELEMENT* el_pt=dynamic_cast<ELEMENT*>(Bulk_mesh_pt->element_pt(e));

  // If the upcast was successful
  if (el_pt!=0)
  {
   // Set the wavenumber function pointer
   el_pt->k_squared_pt()=&GlobalParameters::K_squared;

   // Set the Fourier wavenumber
   el_pt->n_fourier_wavenumber_pt()=&GlobalParameters::N_fourier_wavenumber;
  }
 } // for (unsigned e=0;e<n_element;e++)
 
 // Re-apply boundary conditions
 apply_boundary_conditions();

 // Rebuild the mesh
 rebuild_global_mesh();
} // End of actions_after_adapt

//======================================================start_of_doc======
/// Doc the solution: doc_info contains labels/output directory etc.
//========================================================================
template<class ELEMENT>
void PMLFourierDecomposedHelmholtzProblem<ELEMENT>::doc_solution()
{
 // Tell the user
 oomph_info << "\nDocumentation step: "
	    << GlobalParameters::Doc_info.number() << std::endl;
 
 // Create an output stream
 ofstream some_file;

 // Create space for the file name
 char filename[100];

 // Number of plot points
 unsigned npts=2;
 
 // Number of plot points in the coarse solution
 unsigned npts_coarse=2;

 // Output solution
 //-----------------
 sprintf(filename,"%s/soln%i.dat",
	 GlobalParameters::Doc_info.directory().c_str(),
	 GlobalParameters::Doc_info.number());
 some_file.open(filename);
 Bulk_mesh_pt->output(some_file,npts);
 some_file.close();
 
 // Ouput exact solution
 //---------------------
 sprintf(filename,"%s/exact_soln%i.dat",
	 GlobalParameters::Doc_info.directory().c_str(),
	 GlobalParameters::Doc_info.number());
 some_file.open(filename);
 Bulk_mesh_pt->output_fct(some_file,npts,GlobalParameters::exact_u_pt);
 some_file.close();

 // Output coarse solution
 //-----------------------
 sprintf(filename,"%s/coarse_soln%i.dat",
	 GlobalParameters::Doc_info.directory().c_str(),
	 GlobalParameters::Doc_info.number());
 some_file.open(filename);
 Bulk_mesh_pt->output(some_file,npts_coarse);
 some_file.close();

 // Compute error
 //--------------
 sprintf(filename,"%s/error%i.dat",
	 GlobalParameters::Doc_info.directory().c_str(),
	 GlobalParameters::Doc_info.number());
 some_file.open(filename);
 
 //---------------------------------------------------------------------
 // To compute the norm of the error norm we first need to loop over all
 // of the elements in the mesh. Before we compute the norm of the error
 // in any element we need to make sure it doesn't lie in the PML region
 // or in the pinned region
 //--------------------------------------------------------------------- 
 // Variables to hold the L2 norm of the error in the solution
 double error=0.0;
 
 // Variable to hold the L2 norm of the solution
 double norm=0.0;
 
 Bulk_mesh_pt->compute_error(some_file,
   GlobalParameters::get_exact_u,
   error, norm);

 // Now close the file
 some_file.close();

 // Output the L2 norm of the error and the solution and then output
 // the relative error of the solution
 oomph_info << "\nSolution norm : " << norm
	    << "\nAbsolute error: " << error	   
	    << "\nRelative error: " << error/norm
	    << std::endl;
 
 // Write the L2 norm of the solution to the trace file
 Trace_file << norm << std::endl;

 // Increment the documentation number
 GlobalParameters::Doc_info.number()++;
} // End of doc_solution


//=====================================================start_of_main======
/// Solve 3D Helmholtz problem for a point source in a unit cube
//========================================================================
int main(int argc,char **argv)
{
 //------------------------
 // Command line arguments
 //------------------------
 // Store command line arguments
 CommandLineArgs::setup(argc,argv);

 // Choose the number of nodes in one direction of an element;
 // Values of nnode_1d:
 //        2: Bilinear interpolation 
 //        3: Biquadratic interpolation 
 //        4: Bicubic interpolation 
 CommandLineArgs::specify_command_line_flag(
  "--nnode_1d",&GlobalParameters::Nnode_1d);

 // Choose the minimum level of uniform refinement
 CommandLineArgs::specify_command_line_flag(
  "--min_ref",&GlobalParameters::Min_refinement_level);
 
 // Choose the value of k^2
 CommandLineArgs::specify_command_line_flag(
  "--k_sq",&GlobalParameters::K_squared);
 
 // Parse command line
 CommandLineArgs::parse_and_assign();
  
 // Document what has been specified on the command line
 CommandLineArgs::doc_specified_flags();

 // Update any parameters that need to be updated
 GlobalParameters::update_parameters();
 
 //--------------------------------
 // Set the documentation directory
 //--------------------------------
 // Set output directory
 GlobalParameters::Doc_info.set_directory("RESLT");

 //-------------------
 // Set up the problem
 //-------------------
 // Initialise a null pointer to the class Problem 
 Problem* problem_pt=0;

 // Typedef element name
//  typedef QPMLHelmholtzElement<2,2> ELEMENT;
 typedef QPMLHelmholtzElement<2,2> ELEMENT;
 
 // Set the problem pointer
 problem_pt=new PMLFourierDecomposedHelmholtzProblem<ELEMENT>;

 problem_pt->newton_solve();

return 0;
//------------------ 
 // Solve the problem
 //------------------ 
 // If the user wishes to use adaptive refinement then we use the Newton
 // solver with a given argument to indicate how many adaptations to use
 if (GlobalParameters::Use_adaptation_flag)
 {
  // If the user wishes to silence everything
  if (GlobalParameters::Output_management_flag==3)
  {
   // Store the output stream pointer
   GlobalParameters::Stream_pt=oomph_info.stream_pt();

   // Now set the oomph_info stream pointer to the null stream to
   // disable all possible output
   oomph_info.stream_pt()=&oomph_nullstream;
  }
    
  // Keep refining until the minimum refinement level is reached
  for (unsigned i=0;i<GlobalParameters::Min_refinement_level;i++)
  { 
   oomph_info << "\n===================="
	      << "Initial Refinement"
	      << "====================\n"
	      << std::endl;

   // Add additional refinement
   problem_pt->refine_uniformly();
  }

  // If we silenced the adaptation, allow output again
  if (GlobalParameters::Output_management_flag==3)
  {
   // Now set the oomph_info stream pointer to the null stream to
   // disable all possible output
   oomph_info.stream_pt()=GlobalParameters::Stream_pt;
  }
  
  // Solve the problem
  problem_pt->newton_solve();
   
  // Keep refining until the minimum refinement level is reached
  for (unsigned i=0;i<GlobalParameters::N_adaptations;i++)
  { 
   // If the user wishes to silence everything
   if (GlobalParameters::Output_management_flag==3)
   {
    // Store the output stream pointer
    GlobalParameters::Stream_pt=oomph_info.stream_pt();

    // Now set the oomph_info stream pointer to the null stream to
    // disable all possible output
    oomph_info.stream_pt()=&oomph_nullstream;
   }
   
   // Adapt the problem
   problem_pt->adapt();
   
   // If we silenced the adaptation, allow output again
   if (GlobalParameters::Output_management_flag==3)
   {
    // Now set the oomph_info stream pointer to the null stream to
    // disable all possible output
    oomph_info.stream_pt()=GlobalParameters::Stream_pt;
   }
  
   // Solve the problem
   problem_pt->newton_solve();
  }
 }
 // If the user instead wishes to use uniform refinement
 else
 {
  // If the user wishes to silence everything
  if (GlobalParameters::Output_management_flag==3)
  {
   // Store the output stream pointer
   GlobalParameters::Stream_pt=oomph_info.stream_pt();

   // Now set the oomph_info stream pointer to the null stream to
   // disable all possible output
   oomph_info.stream_pt()=&oomph_nullstream;
  }
  
  // Keep refining until the minimum refinement level is reached
  for (unsigned i=0;i<GlobalParameters::Min_refinement_level;i++)
  { 
   oomph_info << "\n===================="
	      << "Initial Refinement"
	      << "====================\n"
	      << std::endl;

   // Add additional refinement
   problem_pt->refine_uniformly();
  }
 
  // If we silenced the adaptation, allow output again
  if (GlobalParameters::Output_management_flag==3)
  {
   // Now set the oomph_info stream pointer to the null stream to
   // disable all possible output
   oomph_info.stream_pt()=GlobalParameters::Stream_pt;
  }
  
  // Solve the problem
  problem_pt->newton_solve();
 
  // Refine and solve until the additional refinements have been completed
  for (unsigned i=0;i<GlobalParameters::Add_refinement_level;i++)
  {
   // If the user wishes to silence everything
   if (GlobalParameters::Output_management_flag==3)
   {
    // Store the output stream pointer
    GlobalParameters::Stream_pt=oomph_info.stream_pt();

    // Now set the oomph_info stream pointer to the null stream to
    // disable all possible output
    oomph_info.stream_pt()=&oomph_nullstream;
   }
  
   oomph_info << "==================="
	      << "Additional Refinement"
	      << "==================\n"
	      << std::endl;
 
   // Add additional refinement
   problem_pt->refine_uniformly();
  
   // If we silenced the adaptation, allow output again
   if (GlobalParameters::Output_management_flag==3)
   {
    // Now set the oomph_info stream pointer to the null stream to
    // disable all possible output
    oomph_info.stream_pt()=GlobalParameters::Stream_pt;
   }
  
   // Solve the problem
   problem_pt->newton_solve();
  }
 } // if (GlobalParameters::Use_adaptation_flag)
 
 // Delete the problem pointer
 delete problem_pt;

 // Make it a null pointer
 problem_pt=0;
} // End of main
