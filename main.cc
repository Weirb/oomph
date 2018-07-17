//LIC// ====================================================================
//LIC// This file forms part of oomph-lib, the object-oriented, 
//LIC// multi-physics finite-element library, available 
//LIC// at http://www.oomph-lib.org.
//LIC// 
//LIC//    Version 1.0; svn revision $LastChangedRevision: 1097 $
//LIC//
//LIC// $LastChangedDate: 2015-12-17 11:53:17 +0000 (Thu, 17 Dec 2015) $
//LIC// 
//LIC// Copyright (C) 2006-2016 Matthias Heil and Andrew Hazel
//LIC// 
//LIC// This library is free software; you can redistribute it and/or
//LIC// modify it under the terms of the GNU Lesser General Public
//LIC// License as published by the Free Software Foundation; either
//LIC// version 2.1 of the License, or (at your option) any later version.
//LIC// 
//LIC// This library is distributed in the hope that it will be useful,
//LIC// but WITHOUT ANY WARRANTY; without even the implied warranty of
//LIC// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
//LIC// Lesser General Public License for more details.
//LIC// 
//LIC// You should have received a copy of the GNU Lesser General Public
//LIC// License along with this library; if not, write to the Free Software
//LIC// Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA
//LIC// 02110-1301  USA.
//LIC// 
//LIC// The authors may be contacted at oomph-lib@maths.man.ac.uk.
//LIC// 
//LIC//====================================================================
//Driver for a simple 2D poisson problem

//Generic routines
#include "generic.h"

// The Poisson equations
// #include "helmholtz.h"
// #include "poisson.h"
#include "fourier_decomposed_helmholtz.h"

#include "SourceFiles/pml_helmholtz.h"

// The mesh
#include "meshes.h"

using namespace std;

using namespace oomph;


namespace GlobalParameters {

	int N_fourier = 3;

	double K_squared = 10.0;
	double K = sqrt(K_squared);

 unsigned N_terms = 6;

  /// Coefficients in the exact solution
 Vector<double> Coeff(N_terms,1.0);

 /// Imaginary unit 
 std::complex<double> I(0.0,1.0); 

 /// Exact solution as a Vector of size 2, containing real and imag parts
 void get_exact_u(const Vector<double>& x, Vector<double>& u)
 {
  // Switch to spherical coordinates
  double R=sqrt(x[0]*x[0]+x[1]*x[1]);
  
  double theta;
  theta=atan2(x[0],x[1]);
  
  // Argument for Bessel/Hankel functions
  double kr = sqrt(K_squared)*R;  
  
  // Need half-order Bessel functions
  double bessel_offset=0.5;

  // Evaluate Bessel/Hankel functions
  Vector<double> jv(N_terms);
  Vector<double> yv(N_terms);
  Vector<double> djv(N_terms);
  Vector<double> dyv(N_terms);
  double order_max_in=double(N_terms-1)+bessel_offset;
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
  
  // Assemble  exact solution (actually no need to add terms
  // below i=N_fourier as Legendre polynomial would be zero anyway)
  complex<double> u_ex(0.0,0.0);
  for(unsigned i=N_fourier;i<N_terms;i++)
   {
    //Associated_legendre_functions
    double p=Legendre_functions_helper::plgndr2(i,N_fourier,
                                                cos(theta));
    // Set exact solution
    u_ex+=Coeff[i]*sqrt(MathematicalConstants::Pi/(2.0*kr))*(jv[i]+I*yv[i])*p;
   }
  
  // Get the real & imaginary part of the result
  u[0]=u_ex.real();
  u[1]=u_ex.imag();
  
 }//end of get_exact_u

	FiniteElement::SteadyExactSolutionFctPt exact_u_pt=&get_exact_u;
}


//========================================================================
// AnnularQuadMesh, derived from SimpleRectangularQuadMesh.
//========================================================================
template<class ELEMENT> 
class AnnularQuadMesh : public RectangularQuadMesh<ELEMENT>
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
   RectangularQuadMesh<ELEMENT>(n_r,n_phi,1.0,1.0)
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


//====== start_of_problem_class=======================================
/// 2D Poisson problem on rectangular domain, discretised with
/// 2D QPoisson elements. The specific type of element is
/// specified via the template parameter.
//====================================================================
template<class ELEMENT> 
class PoissonProblem : public Problem
{

public:

 /// Constructor: Pass pointer to source function
 PoissonProblem();

 /// Destructor (empty)
 ~PoissonProblem(){}

 /// \short Update the problem specs before solve: Reset boundary conditions
 /// to the values from the exact solution.
 void actions_before_newton_solve();

 /// Update the problem after solve (empty)
 void actions_after_newton_solve(){}

 /// \short Doc the solution. DocInfo object stores flags/labels for where the
 /// output gets written to
 void doc_solution(DocInfo& doc_info);

}; // end of problem class




//=====start_of_constructor===============================================
/// Constructor for Poisson problem: Pass pointer to source function.
//========================================================================
template<class ELEMENT>
PoissonProblem<ELEMENT>::PoissonProblem()
{ 
 // Build annular mesh
 // # of elements in r-direction 
 unsigned n_r=10;
 
 // # of elements in theta-direction 
 unsigned n_theta=10;
 
 // Domain boundaries in theta-direction
 double theta_min=-90.0;
 double theta_max=90.0;
 
 // Domain boundaries in r-direction
 double r_min=1.0;
 double r_max=3.0;
 
 // Build and assign mesh
 mesh_pt() = 
  new AnnularQuadMesh<ELEMENT>(n_r,n_theta,r_min,r_max,theta_min,theta_max);

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

} // end of constructor




//========================================start_of_actions_before_newton_solve===
/// Update the problem specs before solve: (Re-)set boundary conditions
/// to the values from the exact solution.
//========================================================================
template<class ELEMENT>
void PoissonProblem<ELEMENT>::actions_before_newton_solve()
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



//===============start_of_doc=============================================
/// Doc the solution: doc_info contains labels/output directory etc.
//========================================================================
template<class ELEMENT>
void PoissonProblem<ELEMENT>::doc_solution(DocInfo& doc_info)
{ 

 ofstream some_file;
 char filename[100];

 // Number of plot points: npts x npts
 unsigned npts=5;

 // Output solution 
 //-----------------
//  sprintf(filename,"%s/soln%i.dat",doc_info.directory().c_str(),
//          doc_info.number());
//  some_file.open(filename);
//  mesh_pt()->output(some_file,npts);
//  some_file.close();


 // Doc error and return of the square of the L2 error
 //---------------------------------------------------
 double error,norm;
 sprintf(filename,"%s/error%i.dat",doc_info.directory().c_str(),
         doc_info.number());
 some_file.open(filename);
//  std::cout << "The error is here" << std::endl;
 mesh_pt()->compute_error(some_file,
 													GlobalParameters::get_exact_u,
                          error,norm); 
 some_file.close();

 // Doc L2 error and norm of solution
 cout << "Norm of error   : " << sqrt(error) << std::endl
 		  << "Norm of solution: " << sqrt(norm) << std::endl
 			<< "Relative error  : " << sqrt(error/norm) << std::endl;

} // end of doc

 




//===== start_of_main=====================================================
/// Driver code for 2D Poisson problem
//========================================================================
int main()
{

 //Set up the problem
 //------------------

 // Create the problem with 2D nine-node elements from the
 // QPoissonElement family. Pass pointer to source function. 
//  typedef QPoissonElement<2,2> ELEMENT;
 // typedef QFourierDecomposedHelmholtzElement<3> ELEMENT;
 typedef QPMLHelmholtzElement<2,3> ELEMENT;
//  typedef QHelmholtzElement<2,3> ELEMENT;

 PoissonProblem<ELEMENT> problem;

 // Create label for output
 //------------------------
 DocInfo doc_info;

 // Set output directory
 doc_info.set_directory("RESLT");

 // Step number
 doc_info.number()=0;

 // Check that we're ready to go:
 //----------------------------
 cout << "\n\n\nProblem self-test ";
 if (problem.self_test()==0) 
  {
   cout << "passed: Problem can be solved." << std::endl;
  }
 else 
  {
   throw OomphLibError("Self test failed",
                       OOMPH_CURRENT_FUNCTION,
                       OOMPH_EXCEPTION_LOCATION);
  }

	// Solve the problem
	problem.newton_solve();

	//Output the solution
	problem.doc_solution(doc_info);

} //end of main