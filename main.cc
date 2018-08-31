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

#include <fenv.h>
#include "math.h"
#include <complex>

// Generic routines
#include "generic.h"

// For the Bessel functions
#include "oomph_crbond_bessel.h"

#include "SourceFiles/pml_helmholtz.h"

// The mesh
#include "meshes.h"

#include "fenv.h"
using namespace std;

using namespace oomph;


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
 unsigned Linear_solver_flag=1;
 
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

 // Specify the number of V-cycles to perform in the preconditioner
 unsigned Nvcycle=1;
 
 /// DocInfo object used for documentation of the solution
 DocInfo Doc_info;
 
 // Pointer to the output stream -- defaults to oomph_info
 std::ostream* Stream_pt;

 // Use the 1D test problem?
 // 	0 = standard problem, Helmholtz/Poisson
 //		1 = solve the test problem
 unsigned One_dimensional_test_problem = 0;
  
 /// Choose the value of the shift to create the complex-shifted
 /// Laplacian preconditioner (CSLP)
 double Alpha_shift=0.0;

 // Fourier wavenumber
 int N_fourier = 0;

 // Helmholtz frequency
 double K_squared = 0.0;
 double K = sqrt(K_squared);

 // Mesh length in the R direction
 double R_min = 1.0;
 double R_max = R_min + 1.0;
 
 // Mesh length in the Z direction
 double Z_min = 0.0;
 double Z_max = Z_min + 1.0;

 // Mesh length in the R direction (test problem)
 double R_min_test = 1.0;
 double R_max_test = R_min_test + 1.0;
 
 // Mesh length in the Z direction (test problem)
 double Z_min_test = -0.1;
 double Z_max_test = 0.1;

 // Number of terms in the solution
 unsigned N_terms = 6;

 /// Coefficients in the exact solution
 Vector<double> Coeff(N_terms,1.0);

 /// Imaginary unit 
 std::complex<double> I(0.0,1.0); 

 /// Exact solution as a Vector of size 2, containing real and imag parts
 void get_exact_u(const Vector<double>& x, Vector<double>& u)
 {

	 // Solve the one dimensional test problem
	 if (One_dimensional_test_problem == 1){
		
		// Ensure parameters are set for this problem
		N_fourier = 0;
		K_squared = 0.0;
		
		// radial coordinate
		double R=x[0];

		// exact solution
		complex<double> u_ex(0.0,0.0);

		// Solution to (r*u'(r))'/r=0 is u(r)=A*log(r)+B
		u_ex += 6.0*log(R) + 2.0;
		
		// Get the real & imaginary part of the result
		u[0]=u_ex.real();
		u[1]=u_ex.imag();
	 }

	 // We are solving the Poisson problem if k^2 = 0
	 // Avoid divide by 0 error with separate exact solution
	 // Also avoid == comparison with double
	 else {
		 if (K_squared != 0.0) {
			
      // K^2 != 0
      // This is the Helmholtz equation in Fourier decomposed coordinates
      // We are soling for the N_fourier mode

			double R = x[0];
			double Z = x[1];
			
      complex<double> u_ex(0.0, 0.0);

			// Argument for Bessel/Hankel functions
			double kr = sqrt(K_squared)*R;

      N_terms = N_fourier + 1;

			// Evaluate Bessel/Hankel functions
			Vector<double> jv(N_terms);
			Vector<double> yv(N_terms);
			Vector<double> djv(N_terms);
			Vector<double> dyv(N_terms);
			double order_max_in=double(N_terms-1);
			double order_max_out=0;
			
			// This function returns vectors containing 
			// J_k(x), Y_k(x) and their derivatives
			// up to k=order_max, with k increasing in
			// integer increments starting with smallest
			// positive value. So, e.g. for order_max=3.5
			// jv[0] contains J_{1/2}(x),
			// jv[1] contains J_{3/2}(x),
			// jv[2] contains J_{5/2}(x),
			// jv[3] contains J_{7/2}(x).
			CRBond_Bessel::bessjyv(order_max_in,
									kr,
									order_max_out,
									&jv[0],&yv[0],
									&djv[0],&dyv[0]);

      u_ex += (jv[N_fourier]+I*yv[N_fourier])*Z;

      // Get the real & imaginary part of the result
			u[0]=u_ex.real();
			u[1]=u_ex.imag();

		} else if (K_squared == 0.0 && N_fourier != 0) {

      // K^2 = 0, N_fourier != 0
      // This is the Poisson problem in Fourier decomposed coordinates

      // Variables for coordinates
			double R = x[0];
			double Z = x[1];
			
      // Exact solution
      complex<double> u_ex(0.0, 0.0);

      // Add the exact solution to the 
      u_ex += (cosh(N_fourier*log(R)) + I*sinh(N_fourier*log(R)))*Z;

			
			// Get the real & imaginary part of the result
			u[0]=u_ex.real();
			u[1]=u_ex.imag();

		} else if (K_squared == 0.0 && N_fourier == 0) {
			
      // This is the axisymmetric Poisson problem.

			double R = x[0];
			double Z = x[1];

			double m = 3.0;
			
      complex<double> u_ex(0.0, 0.0);

			// Argument for Bessel/Hankel functions
			double kr = m*R;

      N_terms = N_fourier + 1;

			// Evaluate Bessel/Hankel functions
			Vector<double> jv(N_terms);
			Vector<double> yv(N_terms);
			Vector<double> djv(N_terms);
			Vector<double> dyv(N_terms);
			double order_max_in=double(N_terms-1);
			double order_max_out=0;
			
			// This function returns vectors containing 
			// J_k(x), Y_k(x) and their derivatives
			// up to k=order_max, with k increasing in
			// integer increments starting with smallest
			// positive value. So, e.g. for order_max=3.5
			// jv[0] contains J_{1/2}(x),
			// jv[1] contains J_{3/2}(x),
			// jv[2] contains J_{5/2}(x),
			// jv[3] contains J_{7/2}(x).
			CRBond_Bessel::bessjyv(order_max_in,
									kr,
									order_max_out,
									&jv[0],&yv[0],
									&djv[0],&dyv[0]);

      u_ex += (jv[N_fourier]+I*yv[N_fourier])*exp(m*Z);

      // Get the real & imaginary part of the result
			u[0]=u_ex.real();
			u[1]=u_ex.imag();

		}

      // Need one more condition for when K^2=0 and N_fourier = 0.
      // But we may actually ignore this case since we are not interested.

	 }
	 
  
 }//end of get_exact_u

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



//========================================================================
// AnnularQuadMesh, derived from SimpleRectangularQuadMesh.
//========================================================================
template<class ELEMENT> 
class AnnularQuadMesh : public RefineableRectangularQuadMesh<ELEMENT>
{
 
  public:

 // \short Constructor for angular mesh with n_r x n_phi 
 // 2D quad elements. Calls constructor for the underlying 
 // SimpleRectangularQuadMesh; then deforms the mesh so that it fits 
 // into the annular region bounded by the radii r_min and r_max
 // and angles (in degree) of phi_min and phi_max.
 AnnularQuadMesh(const unsigned& n_r, const unsigned& n_phi,
                 const double& r_min, const double& r_max,
                 const double& phi_min, const double& phi_max) :
                 RectangularQuadMesh<ELEMENT>(n_r,n_phi,1.0,1.0,&Mesh::Default_TimeStepper),
                 RefineableRectangularQuadMesh<ELEMENT>(n_r,n_phi,1.0,1.0,&Mesh::Default_TimeStepper)
  {

   // The constructor for the  SimpleRectangularQuadMesh has
   // built the mesh with n_x x n_y = n_r x n_phi elements in the unit
   // square. Let's reposition the nodal points so that the mesh
   // gets mapped into the required annular region:

   // Find out how many nodes there are
   unsigned n_node=this->nnode();
   
   // Loop over all nodes
   for (unsigned n=0;n<n_node;n++)
    {
     // Pointer to node:
     Node* nod_pt=this->node_pt(n);
     
     // Get the x/y coordinates
     double x_old=nod_pt->x(0);
     double y_old=nod_pt->x(1);

     // Map from the old x/y to the new r/phi:
     double r=r_min+(r_max-r_min)*x_old;
     double phi=(phi_min+(phi_max-phi_min)*y_old)*
      MathematicalConstants::Pi/180.0;

     // Set new nodal coordinates
     nod_pt->x(0)=r*cos(phi);
     nod_pt->x(1)=r*sin(phi);
    }
  }
};



//////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////

//============================================start_of_problem_class======
/// Problem class
//========================================================================
template<class ELEMENT>
class PMLStructuredCubicHelmholtz : public HelmholtzMGProblem
{

public:

 /// Constructor
 PMLStructuredCubicHelmholtz();

 /// Destructor (empty)
 ~PMLStructuredCubicHelmholtz();

 /// Doc the solution
 void doc_solution();

 /// Update the problem specs before solve (empty)
 void actions_before_newton_solve();

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

 /// Overload the make_new_problem function to return an object of this class
 HelmholtzMGProblem* make_new_problem()
 {
  // Return a new problem pointer
  return new PMLStructuredCubicHelmholtz<ELEMENT>;
 }

 /// \short Overload the mg_bulk_mesh_pt function to return a pointer to the
 /// "refineable" portion of the mesh
 TreeBasedRefineableMeshBase* mg_bulk_mesh_pt()
 {
  // Return the pointer to the bulk mesh
  return dynamic_cast<TreeBasedRefineableMeshBase*>(mesh_pt());
 }
 
 /// Trace file
 ofstream Trace_file;
}; // End of PMLStructuredCubicHelmholtz class

//==============================================start_of_constructor======
/// Constructor for Helmholtz problem
//========================================================================
template<class ELEMENT>
PMLStructuredCubicHelmholtz<ELEMENT>::PMLStructuredCubicHelmholtz()
{
 // Indicate that the problem is nonlinear to ensure the residual is
 // calculated at the end of the iteration
 problem_is_nonlinear(true);

 // Set the number of Newton iterations to one
 max_newton_iterations()=1;

 // Increase the maximum residuals
 max_residuals()=15.0;

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

 if (GlobalParameters::One_dimensional_test_problem == 1) 
 {
	 // Mesh for the one dimensional test problem
	 mesh_pt() = new RefineableRectangularQuadMesh<ELEMENT>(10,10,
		GlobalParameters::R_min_test,GlobalParameters::R_max_test,
		GlobalParameters::Z_min_test,GlobalParameters::Z_max_test);
 }
 else
 {
	// Create mesh for the Laplacian problem
	mesh_pt() = new RefineableRectangularQuadMesh<ELEMENT>(10,10,
		GlobalParameters::R_min,GlobalParameters::R_max,
		GlobalParameters::Z_min,GlobalParameters::Z_max);
 }

 // Set the boundary conditions for this problem: All nodes are
 // free by default -- only need to pin the ones that have Dirichlet conditions
 // here. 
 unsigned n_bound = mesh_pt()->nboundary();
 for(unsigned i=0;i<n_bound;i++)
  {
   unsigned n_node = mesh_pt()->nboundary_node(i);
   for (unsigned n=0;n<n_node;n++)
    {
     mesh_pt()->boundary_node_pt(i,n)->pin(0); 
		   mesh_pt()->boundary_node_pt(i,n)->pin(1); 
    }
  }

 // Complete the build of all elements so they are fully functional

 // Loop over the elements to set up element-specific 
 // things that cannot be handled by the (argument-free!) ELEMENT 
 // constructor: Pass pointer to source function
 unsigned n_element = mesh_pt()->nelement();
 for(unsigned i=0;i<n_element;i++)
  {
   // Upcast from GeneralsedElement to the present element
   ELEMENT *el_pt = dynamic_cast<ELEMENT*>(mesh_pt()->element_pt(i));

   //Set the source function pointer
   el_pt->k_squared_pt()=&GlobalParameters::K_squared;
	  el_pt->n_fourier_wavenumber_pt()=&GlobalParameters::N_fourier;
  }


 // Setup equation numbering scheme
 cout <<"Number of equations: " << assign_eqn_numbers() << std::endl; 

} // End of constructor

//===============================================start_of_destructor======
/// Destructor for Helmholtz problem
//========================================================================
template<class ELEMENT>
PMLStructuredCubicHelmholtz<ELEMENT>::~PMLStructuredCubicHelmholtz()
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
   
 // // Delete the error estimator
 // delete mesh_pt()->spatial_error_estimator_pt();
 
 // // Set the pointer to null
 // mesh_pt()->spatial_error_estimator_pt()=0;

 // Delete the "bulk" mesh
 delete mesh_pt();

 // Set the pointer to null
 mesh_pt()=0;
 
} // End of ~PMLStructuredCubicHelmholtz

//=======set_gmres_multigrid_solver=======================================
/// Build and set GMRES preconditioner by multigrid as the linear solver
//========================================================================
template<class ELEMENT>
void PMLStructuredCubicHelmholtz<ELEMENT>::set_gmres_multigrid_solver()
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
 //SuperLUPreconditioner* prec_pt = new SuperLUPreconditioner();

 // Set preconditioner
 solver_pt->preconditioner_pt()=prec_pt;
  
 // Set the shift
 prec_pt->alpha_shift()=GlobalParameters::Alpha_shift;

 // Set the number of V-cycles
 prec_pt->n_v_cycle()=GlobalParameters::Nvcycle;
   
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


//========================================start_of_actions_before_newton_solve===
/// Update the problem specs before solve: (Re-)set boundary conditions
/// to the values from the exact solution.
//========================================================================
template<class ELEMENT>
void PMLStructuredCubicHelmholtz<ELEMENT>::actions_before_newton_solve()
{
 // How many boundaries are there?
 unsigned n_bound = mesh_pt()->nboundary();
 
 //Loop over the boundaries
 for(unsigned i=0;i<n_bound;i++)
  {
   // How many nodes are there on this boundary?
   unsigned n_node = mesh_pt()->nboundary_node(i);

   // Loop over the nodes on boundary
   for (unsigned n=0;n<n_node;n++)
    {
     // Get pointer to node
     Node* nod_pt=mesh_pt()->boundary_node_pt(i,n);

     // Extract nodal coordinates from node:
     Vector<double> x(2);
     x[0]=nod_pt->x(0);
     x[1]=nod_pt->x(1);

     // Compute the value of the exact solution at the nodal point
     Vector<double> u(2);
     GlobalParameters::get_exact_u(x,u);

     // Assign the value to the one (and only) nodal value at this node
     nod_pt->set_value(0,u[0]);
		 nod_pt->set_value(1,u[1]);
    }
  } 
}  // end of actions before solve


//================================start_of_apply_boundary_conditions======
/// Apply boundary conditions
//========================================================================
template<class ELEMENT>
void PMLStructuredCubicHelmholtz<ELEMENT>::apply_boundary_conditions()
{
 // How many boundaries are there?
 unsigned n_bound = mesh_pt()->nboundary();
 
 //Loop over the boundaries
 for(unsigned i=0;i<n_bound;i++)
  {
   // How many nodes are there on this boundary?
   unsigned n_node = mesh_pt()->nboundary_node(i);

   // Loop over the nodes on boundary
   for (unsigned n=0;n<n_node;n++)
    {
     // Get pointer to node
     Node* nod_pt=mesh_pt()->boundary_node_pt(i,n);

     // Extract nodal coordinates from node:
     Vector<double> x(2);
     x[0]=nod_pt->x(0);
     x[1]=nod_pt->x(1);

     // Compute the value of the exact solution at the nodal point
     Vector<double> u(2);
     GlobalParameters::get_exact_u(x,u);

     // Assign the value to the one (and only) nodal value at this node
     nod_pt->pin(0); 
		   nod_pt->pin(1); 

     nod_pt->set_value(0,u[0]);
		   nod_pt->set_value(1,u[1]);
    }
  } 
} // End of apply_boundary_conditions



//======================================start_of_actions_after_adapt======
/// Actions after adapt: Re-apply the boundary conditions
//========================================================================
template<class ELEMENT>
void PMLStructuredCubicHelmholtz<ELEMENT>::actions_after_adapt()
{
 // Complete the build of all elements so they are fully functional: 
 //-----------------------------------------------------------------
 // How many elements in the mesh?
 unsigned n_element=mesh_pt()->nelement();

 // Loop over the elements and pass a pointer to the value of k^2
 for (unsigned e=0;e<n_element;e++)
 {
  // Upcast from GeneralisedElement to Helmholtz bulk element
  ELEMENT* el_pt=dynamic_cast<ELEMENT*>(mesh_pt()->element_pt(e));

  // If the upcast was successful
  if (el_pt!=0)
  {
   // Set the wavenumber function pointer
   el_pt->k_squared_pt()=&GlobalParameters::K_squared;
   el_pt->n_fourier_wavenumber_pt()=&GlobalParameters::N_fourier;
  }
 } // for (unsigned e=0;e<n_element;e++)
 
 // Re-apply boundary conditions
 apply_boundary_conditions();

} // End of actions_after_adapt

//======================================================start_of_doc======
/// Doc the solution: doc_info contains labels/output directory etc.
//========================================================================
template<class ELEMENT>
void PMLStructuredCubicHelmholtz<ELEMENT>::doc_solution()
{
 // Tell the user
 oomph_info << "\nDocumentation step: "
	    << GlobalParameters::Doc_info.number() << std::endl;
 
 // Create an output stream
 ofstream some_file;

 // Create space for the file name
 char filename[100];

 // Number of plot points
 unsigned npts=5;
 
 // Number of plot points in the coarse solution
 unsigned npts_coarse=2;

 // Output solution
 //-----------------
 sprintf(filename,"%s/soln%i.dat",
	 GlobalParameters::Doc_info.directory().c_str(),
	 GlobalParameters::Doc_info.number());
 some_file.open(filename);
 mesh_pt()->output(some_file,npts);
 some_file.close();
 
 // Ouput exact solution
 //---------------------
 sprintf(filename,"%s/exact_soln%i.dat",
	 GlobalParameters::Doc_info.directory().c_str(),
	 GlobalParameters::Doc_info.number());
 some_file.open(filename);
 mesh_pt()->output_fct(some_file,npts,GlobalParameters::exact_u_pt);
 some_file.close();

 // Output coarse solution
 //-----------------------
 sprintf(filename,"%s/coarse_soln%i.dat",
	 GlobalParameters::Doc_info.directory().c_str(),
	 GlobalParameters::Doc_info.number());
 some_file.open(filename);
 mesh_pt()->output(some_file,npts_coarse);
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
 
 // Doc error and return of the square of the L2 error
 //---------------------------------------------------
 sprintf(filename,"%s/error%i.dat",
	 GlobalParameters::Doc_info.directory().c_str(),
	 GlobalParameters::Doc_info.number());
 some_file.open(filename);
//  std::cout << "The error is here" << std::endl;
 mesh_pt()->compute_error(some_file,
 													GlobalParameters::get_exact_u,
                          error,norm); 
 some_file.close();

 // Doc L2 error and norm of solution
 cout << "Norm of solution: " << sqrt(norm) << std::endl
 	 	  << "Norm of error   : " << sqrt(error) << std::endl
      << "Relative error  : " << sqrt(error/norm) << std::endl;

 // Increment the documentation number
 GlobalParameters::Doc_info.number()++;
} // End of doc_solution


//=====================================================start_of_main======
/// Solve 3D Helmholtz problem for a point source in a unit cube
//========================================================================
int main(int argc,char **argv)
{

	feenableexcept(FE_INVALID | FE_DIVBYZERO | FE_OVERFLOW | FE_UNDERFLOW);
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
 
 // Choose the additional levels of uniform refinement
 CommandLineArgs::specify_command_line_flag(
  "--add_ref",&GlobalParameters::Add_refinement_level);
 
 // Choose the maximum number of adaptive refinements
 CommandLineArgs::specify_command_line_flag(
  "--n_adapt",&GlobalParameters::N_adaptations);
 
 // Choose how many additional levels of uniform refinement to use
 CommandLineArgs::specify_command_line_flag(
  "--use_adapt",&GlobalParameters::Use_adaptation_flag);
  
 // Choose the value of k^2
 CommandLineArgs::specify_command_line_flag(
  "--k_sq",&GlobalParameters::K_squared);

 // Choose the value of Fourier wavenumber, N_fourier
 CommandLineArgs::specify_command_line_flag(
  "--n_fourier",&GlobalParameters::N_fourier);

 // Choose the number of V-cycles in the preconditioner 
 CommandLineArgs::specify_command_line_flag(
  "--vcycles",&GlobalParameters::Nvcycle);

 // Choose the value of the shift in the CSLP
 CommandLineArgs::specify_command_line_flag(
  "--alpha",&GlobalParameters::Alpha_shift);
 
 // Choose the value of the damping factor in the damped Jacobi solver
 CommandLineArgs::specify_command_line_flag(
  "--omega",&Smoother_Factory_Function_Helper::Omega);
  
 // Choose the pre-smoother
 CommandLineArgs::specify_command_line_flag(
  "--presmoother",&GlobalParameters::Pre_smoother_flag);
  
 // Choose the post-smoother
 CommandLineArgs::specify_command_line_flag(
  "--postsmoother",&GlobalParameters::Post_smoother_flag);
  
 // Choose the linear solver
 CommandLineArgs::specify_command_line_flag(
  "--linear_solver",&GlobalParameters::Linear_solver_flag);
 
 // Decide whether or not to suppress all or some of the MG solver output
 CommandLineArgs::specify_command_line_flag(
  "--output_flag",&GlobalParameters::Output_management_flag);
     
 // Decide whether or not to display convergence information
 CommandLineArgs::specify_command_line_flag(
  "--conv_flag",&GlobalParameters::Doc_convergence_flag);
 
 // Solve the 'thin' test problem, with no z dependence
  CommandLineArgs::specify_command_line_flag(
  "--1d_test",&GlobalParameters::One_dimensional_test_problem);
 
 // Parse command line
 CommandLineArgs::parse_and_assign();
  
 // Document what has been specified on the command line
 CommandLineArgs::doc_specified_flags();

 
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
 
 typedef RefineableQPMLHelmholtzElement<2,3> ELEMENT;
  
 // Set the problem pointer
 problem_pt=new PMLStructuredCubicHelmholtz<ELEMENT>;

 // problem_pt->refine_uniformly();

 // problem_pt->newton_solve();

 // delete problem_pt;
 // problem_pt = 0;

 // return 0;


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
