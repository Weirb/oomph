/*

2d solving of the same problem
Need to change the coordinates and stuff
Also need to change the equations
- This means change the element

*/



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

// The Helmholtz equations and complex-valued multigrid machinery
#include "SourceFiles/pml_helmholtz.h"

// The mesh
#include "meshes/simple_cubic_mesh.h"
#include "meshes/simple_rectangular_quadmesh.h"

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
 
 /// DocInfo object used for documentation of the solution
 DocInfo Doc_info;
 
 // Pointer to the output stream -- defaults to oomph_info
 std::ostream* Stream_pt;
  
 /// \short Problem specific parameters:
 ///------------------------------------
 /// Length of the cube in each direction
 double Lx=1.0;
 double Ly=1.0;
 // double Lz=1.0;

 /// Number of elements in each direction (used by SimpleCubicMesh)
 unsigned Nx=27;
 unsigned Ny=27;
 // unsigned Nz=7;
 
 /// The element width
 double Element_width=Lx/double(Nx);
 
 /// Length of cube in each direction
 double Pml_thickness=Element_width;
 
 /// Store the value of Pi
 double Pi=MathematicalConstants::Pi;
 
 /// Square of the wavenumber (also known as k^2)
 double K_squared=20.0;

 /// Wavenumber (also known as k),k=omega/c
 double Wavenumber=sqrt(K_squared);

 /// Update the parameters passed in at the command line
 void update_parameters()
 {
  /// Wavenumber (also known as k), k=omega/c
  Wavenumber=sqrt(K_squared);
 }

 /// Fourier wavenumber
 int N_fourier_wavenumber=0;
 
 /// The x and y coordinate of the centre of the cube 
 double Centre=Lx/2.0;

 /// Get the exact solution, u, at the spatial position, x
 void get_simple_exact_u(const Vector<double>& x,Vector<double>& u)
 {
  // Initialise a variable to store the radial distance
  double r=std::sqrt((x[0]-Centre)*(x[0]-Centre)
		     +(x[1]-Centre)*(x[1]-Centre));

  // Scale the radial distance by the wavenumber
  double kr=Wavenumber*r;

  // The solution is singular at the centre so set it to zero 
  if (r==0.0)
  {
   // Set the real part of the solution value
   u[0]=0.0;
   
   // Set the imaginary part of the solution value
   u[1]=0.0;
  }
  // Otherwise set the correct solution value
  else
  {
   // Set the real part of the solution value
   u[0]=cos(kr)/kr;
   
   // Set the imaginary part of the solution value
   u[1]=sin(kr)/kr;
  }
 } // End of get_simple_exact_u

 // Set the exact solution pointer to the get_simple_exact_u function above
 FiniteElement::SteadyExactSolutionFctPt simple_exact_u_pt=&get_simple_exact_u;

 /// \short New mapping function that makes the mapping independent of the
 /// PML thickness
 class TestPMLMapping : public virtual PMLMapping
 {
 public:

  /// Default constructor (empty)
  TestPMLMapping(){};

  /// \short Overwrite the pure PML mapping coefficient function to return the
  /// coeffcients proposed by Bermudez et al
  std::complex<double> gamma(const double& nu_i,
			     const double& pml_width_i,
			     const double& k_squared_local,
			     const double& alpha_shift)
   {
    // The "effective k^2" is shifted, so we shift the k used in the
    // transformation too
    std::complex<double> k_shifted=
     sqrt(k_squared_local*std::complex<double>(1.0,alpha_shift));

    // Return the gamma in J++, with the shifted k
    return (1.0/k_shifted)*std::complex<double>
     (0.0,1.0/(std::fabs(pml_width_i-nu_i)));    
   } // End of gamma
 }; // End of TestPMLMapping

 /// Set the new PML mapping
 TestPMLMapping* Test_pml_mapping_pt=new TestPMLMapping;
 
 /// \short The choice of whether or not to enable the new test mapping
 ///    1 = Enable test mapping
 ///    0 = Disable test mapping
 unsigned Enable_test_pml_mapping_flag=1;
 
 /// The tolerance for a point relative to the bounding inner square
 double Eps=1.0e-12;
   
 /// \short Function to determine whether or not a point lies in the centre
 /// of the mesh (in the pinned region)
 bool is_in_pinned_region(const Vector<double>& x)
 {
  // Check if the element lies in the central cube region
  // return (abs(x[0]-GlobalParameters::Centre)<
	 //  (0.5*GlobalParameters::Element_width+Eps)&&
	 //  abs(x[1]-GlobalParameters::Centre)<
	 //  (0.5*GlobalParameters::Element_width+Eps)&&
	 //  abs(x[2]-GlobalParameters::Centre)<
	 //  (0.5*GlobalParameters::Element_width+Eps));
  return (abs(x[0]-GlobalParameters::Centre)<
	  (0.5*GlobalParameters::Element_width+Eps)&&
	  abs(x[1]-GlobalParameters::Centre)<
	  (0.5*GlobalParameters::Element_width+Eps));
 } // End of is_in_pinned_region 
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
class PMLFourierDecomposedHelmholtzProblem : public Problem
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
 RectangularQuadMesh<ELEMENT>* Bulk_mesh_pt;
 
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
 
 // Open trace file
 Trace_file.open("RESLT/trace.dat");

 // Build the mesh using the specified parameters:
 //-----------------------------------------------
 // Build the "bulk" mesh
 Bulk_mesh_pt=new RectangularQuadMesh<ELEMENT>(
  GlobalParameters::Nx,GlobalParameters::Ny,
  GlobalParameters::Lx,GlobalParameters::Ly);

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
   
   // If we're using Jonathon's new test mapping
   if (GlobalParameters::Enable_test_pml_mapping_flag)
   {
    // Set the PML mapping function
    el_pt->pml_mapping_pt()=GlobalParameters::Test_pml_mapping_pt;
   }
  } // if (el_pt!=0)
 } // for (unsigned e=0;e<n_element;e++)

 // Apply the boundary conditions, both in the central region and on the
 // outer boundary (since these nodes are PML nodes)
 apply_boundary_conditions();
 
 // Enable the PML mapping in elements in the PML region
 enable_pmls();
 
 // Setup equation numbering scheme
 assign_eqn_numbers(); 
} // End of constructor

//===============================================start_of_destructor======
/// Destructor for Helmholtz problem
//========================================================================
template<class ELEMENT>
PMLFourierDecomposedHelmholtzProblem<ELEMENT>::~PMLFourierDecomposedHelmholtzProblem()
{   
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
    
 // Loop over the elements in the mesh
 for (unsigned e=0;e<n_element;e++)
 {
  // Upcast from GeneralisedElement to Helmholtz bulk element
  ELEMENT* el_pt=dynamic_cast<ELEMENT*>(Bulk_mesh_pt->element_pt(e));

  // If the upcast was successful
  if (el_pt!=0)
  {
   // Get the (Eulerian) coordinates of the centre of the element
   el_pt->get_x(s,x);

   // Check if the element lies in the central cube region
   if (GlobalParameters::is_in_pinned_region(x))
   {
    // Calculate the number of nodes in the element
    unsigned nnode=el_pt->nnode();

    // Loop over all of the nodes in the element
    for (unsigned i=0;i<nnode;i++)
    {
     // Create a node pointer to store the i-th node in the element
     Node* node_pt=el_pt->node_pt(i);

     // Get the spatial position of this node
     for (unsigned k=0;k<2;k++)
     {
      // Store the k-th coordinate value in the vector, x
      x[k]=node_pt->x(k);
     }

     // Get the exact solution at this (Eulerian) position
     GlobalParameters::get_simple_exact_u(x,u);

     // Make sure each dof at this point is pinned (real and imaginary)
     node_pt->pin(0);
     node_pt->pin(1);

     // Set the solution value at this point
     node_pt->set_value(0,u[0]);
     node_pt->set_value(1,u[1]);
    }
   } // if(abs(x[0]-GlobalParameters::Centre) < 0.51 ... 
  } // if (el_pt!=0)
 } // for (unsigned e=0;e<n_element;e++)
          
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

   // Set the solution value at this point (real part)
   boundary_node_pt->set_value(0,0.0);
   
   // Set the solution value at this point (imaginary part)
   boundary_node_pt->set_value(1,0.0);
  }
 } // for(unsigned b=0;b<n_bound;b++)
} // End of apply_boundary_conditions

//==============================================start_of_enable_pmls======
/// Enable the PML mapping function for each node in the PML region
//========================================================================
template<class ELEMENT>
void PMLFourierDecomposedHelmholtzProblem<ELEMENT>::enable_pmls()
{
 // Find the number of elements in the mesh
 unsigned n_element=Bulk_mesh_pt->nelement();

 // Vector to hold the local coordinates of a point in any given element 
 Vector<double> s(2,0.0);
 
 // Vector to hold the (Eulerian) coordinates of a point
 Vector<double> x(2,0.0);

 // Vector to hold the real and imaginary part of the solution
 Vector<double> u(2,0.0);

 // Store the required coordinate of the inner boundary of the left PML; 
 // in any given direction this will be the value of Pml_thickness
 double left_boundary=GlobalParameters::Pml_thickness;
   
 // Store the required coordinate of the inner boundary of the right
 // PML; in the x-direction this will be the value of Lx-Pml_thickness
 // (or Ly-Pml_thickness in the y-direction and Lz-Pml_thickness in
 // the z-direction) but we assume the PML has the same thickness in
 // all directions
 double right_boundary=GlobalParameters::Lx-GlobalParameters::Pml_thickness;

 // Loop over the elements in the mesh
 for (unsigned e=0;e<n_element;e++)
 {
  // Upcast from GeneralisedElement to Helmholtz bulk element
  ELEMENT* el_pt=dynamic_cast<ELEMENT*>(Bulk_mesh_pt->element_pt(e));

  // If the upcast was successful
  if (el_pt!=0)
  {
   // If we're using Jonathon's new test mapping
   if (GlobalParameters::Enable_test_pml_mapping_flag)
   {
    // Set the PML mapping function
    el_pt->pml_mapping_pt()=GlobalParameters::Test_pml_mapping_pt;
   }
  
   // Get the (Eulerian) coordinates of the centre of the element
   el_pt->get_x(s,x);
   
   // If it's in the left (x-direction) PML region
   if (x[0]<=left_boundary)
    el_pt->enable_pml(0,left_boundary,0.0);

   // If it's in the right (x-direction) PML region
   if (x[0]>=right_boundary)
    el_pt->enable_pml(0,right_boundary,GlobalParameters::Lx);

   // If it's in the left (y-direction) PML region
   if (x[1]<=left_boundary)
    el_pt->enable_pml(1,left_boundary,0.0);

   // If it's in the right (y-direction) PML region
   if (x[1]>=right_boundary)
    el_pt->enable_pml(1,right_boundary,GlobalParameters::Ly);
  }
 } // for (unsigned e=0;e<n_element;e++)
} // End of enable_pmls

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

 // Re-enable the PML mapping in elements in the PML region
 enable_pmls();

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
 Bulk_mesh_pt->output_fct(some_file,npts,GlobalParameters::simple_exact_u_pt);
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
 
 // Vector to hold the local coordinates of a point in an element
 Vector<double> s(2,0.0);

 // Vector to hold the spatial position of a point in an element
 Vector<double> x(2,0.0);

 // Store the required coordinate of the inner boundary of the left PML; 
 // in any given direction this will be the value of Pml_thickness
 double left_boundary=GlobalParameters::Pml_thickness;
   
 // Store the required coordinate of the inner boundary of the right
 // PML; in the x-direction this will be the value of Lx-Pml_thickness
 // (or Ly-Pml_thickness in the y-direction and Lz-Pml_thickness in
 // the z-direction) but we assume the PML has the same thickness in
 // all directions
 double right_boundary=GlobalParameters::Lx-GlobalParameters::Pml_thickness;

 // Find out how many elements there are in the mesh
 unsigned n_element=Bulk_mesh_pt->nelement();

 // Loop over all of the elements in the mesh
 for (unsigned e=0;e<n_element;e++)
 { 
  // Variables to hold the L2 norm of the error in the elemental solution
  double el_error=0.0;
 
  // Variable to hold the L2 norm of the elemental solution
  double el_norm=0.0;
 
  // Upcast from GeneralisedElement to Helmholtz bulk element
  ELEMENT* el_pt=dynamic_cast<ELEMENT*>(Bulk_mesh_pt->element_pt(e));

  // If the upcast was successful
  if (el_pt!=0)
  {
   // Get the (Eulerian) coordinates of the centre of the element
   el_pt->get_x(s,x);

   // We only take the contribution from this element if it does
   // not lie in the PML region 
   if(x[0]<=left_boundary) continue;
   if(x[0]>=right_boundary) continue;
   if(x[1]<=left_boundary) continue;
   if(x[1]>=right_boundary) continue;
   // if(x[2]<=left_boundary) continue;
   // if(x[2]>=right_boundary) continue;

   // If it's in the (pinned) central region, ignore it 

   // Check if the element lies in the central cube region
   if (GlobalParameters::is_in_pinned_region(x))
   {
    // Skip to the next element
    continue;
   }

   // Otherwise, compute the L2 norm of the error over this element
   el_pt->compute_error(some_file,
			GlobalParameters::get_simple_exact_u,
			el_error,
			el_norm);

   // Update the global error norm value
   error+=el_error;
   
   // Update the global norm value
   norm+=el_norm;
  }
 } // for(unsigned e=0;e<n_element;e++)

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
 
 // Decide whether or not to display convergence information
 CommandLineArgs::specify_command_line_flag(
  "--test_pml_mapping",&GlobalParameters::Enable_test_pml_mapping_flag);
 
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
 typedef QPMLHelmholtzElement<2,2> ELEMENT;
 
 // Set the problem pointer
 problem_pt=new PMLFourierDecomposedHelmholtzProblem<ELEMENT>;

 //------------------ 
 // Solve the problem
 //------------------ 
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
  
 // Solve the problem
 problem_pt->newton_solve();
 
 // Delete the problem pointer
 delete problem_pt;

 // Make it a null pointer
 problem_pt=0;
} // End of main
