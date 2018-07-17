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

// The mesh
#include "meshes.h"

using namespace std;

using namespace oomph;


namespace GlobalParameters {

	int N_fourier = 0;

	double K_squared = 3.0;
	double K = sqrt(K_squared);

	double M_squared = 1.0;
	double M = sqrt(M_squared);

  // Cartesian solution
	// void get_exact_u(const Vector<double> &x, Vector<double> &u){
	// 	u[0] = sin(sqrt(M_squared+K_squared)*x[1])*exp(-M*x[0]);
	// 	u[1] = 0.0;
	// }
	
	void get_exact_u(const Vector<double> &x, Vector<double> &u){

   int n_terms = N_fourier + 1;

   Vector<double> jv(n_terms);
   Vector<double> djv(n_terms);
   Vector<double> yv(n_terms);
   Vector<double> dyv(n_terms);
   
   double n_actual = 0;
   CRBond_Bessel::bessjyv(n_terms, 
                          x[0]*sqrt(K_squared-M_squared),
                          n_actual,
                          &jv[0],&yv[0],
                          &djv[0],&dyv[0]);

  // double j0,j1,y0,y1,j0p,j1p,y0p,y1p;
  // CRBond_Bessel::bessjy01a(ksq_lsq*r,j0,j1,y0,y1,j0p,j1p,y0p,y1p);
  
  complex<double> u_ex(0.0, 0.0);

	for (unsigned j=0; j<5; ++j){
		// u_ex+= jv[j]*exp(-M*x[1]);
		u_ex += sin(sqrt(M_squared+K_squared)*x[1])*exp(-M*x[0]);

	}
  
  u[0]=u_ex.real();
  u[1]=u_ex.imag();
 }

	FiniteElement::SteadyExactSolutionFctPt exact_u_pt=&get_exact_u;
}


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
 // Setup mesh

 // # of elements in x-direction
 unsigned n_x=4;

 // # of elements in y-direction
 unsigned n_y=4;

 // Domain length in x-direction
 double l_x=1.0;

 // Domain length in y-direction
 double l_y=1.0;

 // Build and assign mesh
 mesh_pt() = new RectangularQuadMesh<ELEMENT>(n_x,n_y,l_x,l_y);

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
	//  el_pt->fourier_wavenumber_pt()=&GlobalParameters::N_fourier;
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
 typedef QFourierDecomposedHelmholtzElement<3> ELEMENT;
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