//LIC// ====================================================================
//LIC// This file forms part of oomph-lib, the object-oriented, 
//LIC// multi-physics finite-element library, available 
//LIC// at http://www.oomph-lib.org.
//LIC// 
//LIC//    Version 1.0; svn revision $LastChangedRevision: 1259 $
//LIC//
//LIC// $LastChangedDate: 2016-11-08 23:52:03 +0000 (Tue, 08 Nov 2016) $
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
#include "refineable_pml_helmholtz_elements.h"


namespace oomph
{
 //======================================================================
 /// Compute element residual Vector and/or element Jacobian matrix
 ///
 /// flag=1: compute both
 /// flag=0: compute only residual Vector
 ///
 /// Pure version without hanging nodes
 //======================================================================
 template <unsigned DIM>
 void  RefineablePMLHelmholtzEquations<DIM>::
 fill_in_generic_residual_contribution_helmholtz(Vector<double> &residuals,
						 DenseMatrix<double> &jacobian,
						 const unsigned& flag)
 {
  //Find out how many nodes there are
  const unsigned n_node = nnode();

  //Set up memory for the shape and test functions
  Shape psi(n_node), test(n_node);
  DShape dpsidx(n_node,DIM), dtestdx(n_node,DIM);

  // Local storage for pointers to hang_info objects
  HangInfo *hang_info_pt=0, *hang_info2_pt=0;

  //Set the value of n_intpt
  const unsigned n_intpt = integral_pt()->nweight();

  //Integers to store the local equation and unknown numbers
  int local_eqn_real=0, local_unknown_real=0;
  int local_eqn_imag=0, local_unknown_imag=0;

  //Loop over the integration points
  for(unsigned ipt=0;ipt<n_intpt;ipt++)
  {
   //Get the integral weight
   double w = integral_pt()->weight(ipt);

   //Call the derivatives of the shape and test functions
   double J = this->dshape_and_dtest_eulerian_at_knot_helmholtz(ipt,psi,dpsidx,
                                                                test,dtestdx);

   //Premultiply the weights and the Jacobian
   double W = w*J;

   //Calculate local values of unknown
   //Allocate and initialise to zero
   std::complex<double> interpolated_u(0.0,0.0);
   Vector<double> interpolated_x(DIM,0.0);
   Vector< std::complex<double> > interpolated_dudx(DIM);

   //Calculate function value and derivatives:
   //-----------------------------------------
   // Loop over nodes
   for(unsigned l=0;l<n_node;l++)
   {
    // Loop over directions
    for(unsigned j=0;j<DIM;j++)
    {
     interpolated_x[j] += nodal_position(l,j)*psi(l);
    }

    //Get the nodal value of the helmholtz unknown
    const std::complex<double>
     u_value(this->nodal_value(l,this->u_index_helmholtz().real()),
	     this->nodal_value(l,this->u_index_helmholtz().imag()));

    //Add to the interpolated value
    interpolated_u += u_value*psi(l);

    // Loop over directions
    for(unsigned j=0;j<DIM;j++)
    {
     interpolated_dudx[j] += u_value*dpsidx(l,j);
    }
   }

    // Create complex variable for the real radial distance
    std::complex<double> r(interpolated_x[0],0.0);
    std::complex<double> rr = r*r;

    // Fourier wavenumber
    double n = (double)this->n_fourier_wavenumber();
    double nn = n*n;

    // Complex alpha contribution
    std::complex<double> alpha_factor(1.0, this->alpha());

    // Compute the Jacobians
    //----------------------

    // Jacobian for the Helmholtz term
    std::complex<double> k_squared_jacobian =
      -alpha_factor*(this->k_squared()-nn/rr)*r*W;

    // Jacobian for the Laplacian
    Vector<std::complex<double> > laplace_jacobian(2);
    for (unsigned k=0;k<2;k++)
    {
      laplace_jacobian[k] = r*W;
    }

    // Compute the residuals
    //----------------------

    // Residual for the Helmholtz term
    std::complex<double> k_squared_residual =
    k_squared_jacobian*interpolated_u;

    // Residual for the Laplacian
    Vector<std::complex<double> > laplace_residual(2);
    for (unsigned k=0;k<2;k++)
    {
      laplace_residual[k] = laplace_jacobian[k]*interpolated_dudx[k];
    }
   
   // Assemble residuals and Jacobian
   //--------------------------------
   // Loop over the test functions
   for(unsigned l=0;l<n_node;l++)
   {
     
    //Local variables used to store the number of master nodes and the
    //weight associated with the shape function if the node is hanging
    unsigned n_master=1; double hang_weight=1.0;
     
    //Local bool (is the node hanging)
    bool is_node_hanging = this->node_pt(l)->is_hanging();
     
    //If the node is hanging, get the number of master nodes
    if(is_node_hanging)
    {
     hang_info_pt = this->node_pt(l)->hanging_pt();
     n_master = hang_info_pt->nmaster();
    }
    //Otherwise there is just one master node, the node itself
    else
    {
     n_master = 1;
    }
     
    //Loop over the master nodes
    for(unsigned m=0;m<n_master;m++)
    {
     //Get the local equation number and hang_weight
     //If the node is hanging
     if(is_node_hanging)
     {
      //Read out the local equation number from the m-th master node
      local_eqn_real = 
       this->local_hang_eqn(hang_info_pt->master_node_pt(m),
			    this->u_index_helmholtz().real());
         
      local_eqn_imag = 
       this->local_hang_eqn(hang_info_pt->master_node_pt(m),
			    this->u_index_helmholtz().imag());
         
      //Read out the weight from the master node
      hang_weight = hang_info_pt->master_weight(m);
     }
     //If the node is not hanging
     else
     {
      //The local equation number comes from the node itself
      local_eqn_real = 
       this->nodal_local_eqn(l,this->u_index_helmholtz().real());         
      local_eqn_imag = 
       this->nodal_local_eqn(l,this->u_index_helmholtz().imag());
         
      //The hang weight is one
      hang_weight = 1.0;
     }
       
     // first, compute the real part contribution
     //-------------------------------------------
       
     /*IF it's not a boundary condition*/
     if(local_eqn_real >= 0)
     {
        // Add the Helmholtz term to the residuals
        residuals[local_eqn_real]+=k_squared_residual.real()*test(l)*hang_weight;

        // Add the Laplace term to the residuals
        for (unsigned k=0;k<2;k++)
        {
          residuals[local_eqn_real]+=laplace_residual[k].real()*dtestdx(l,k)*hang_weight;
        }
         
      // Calculate the jacobian
      //-----------------------
      if(flag)
      {
       //Local variables to store the number of master nodes
       //and the weights associated with each hanging node
       unsigned n_master2=1; double hang_weight2=1.0;
           
       //Loop over the nodes for the variables
       for(unsigned l2=0;l2<n_node;l2++)
            { 
        //Local bool (is the node hanging)
        bool is_node2_hanging = this->node_pt(l2)->is_hanging();
                  
        //If the node is hanging, get the number of master nodes
        if(is_node2_hanging)
        {
        hang_info2_pt = this->node_pt(l2)->hanging_pt();
        n_master2 = hang_info2_pt->nmaster();
        }
        //Otherwise there is one master node, the node itself
        else
        {
        n_master2 = 1;
        }
                  
        //Loop over the master nodes
        for(unsigned m2=0;m2<n_master2;m2++)
        {
        //Get the local unknown and weight
        //If the node is hanging
        if(is_node2_hanging)
        {
          //Read out the local unknown from the master node
          local_unknown_real = 
          this->local_hang_eqn(hang_info2_pt->master_node_pt(m2),
              this->u_index_helmholtz().real());
          local_unknown_imag = 
          this->local_hang_eqn(hang_info2_pt->master_node_pt(m2),
              this->u_index_helmholtz().imag());
                      
          //Read out the hanging weight from the master node
          hang_weight2 = hang_info2_pt->master_weight(m2);
        }
        //If the node is not hanging
        else
        {
          //The local unknown number comes from the node
          local_unknown_real = 
          this->nodal_local_eqn(l2,this->u_index_helmholtz().real());
                      
          local_unknown_imag = 
          this->nodal_local_eqn(l2,this->u_index_helmholtz().imag());
                      
          //The hang weight is one
          hang_weight2 = 1.0;
        }
                    
                    
        //If at a non-zero degree of freedom add in the entry
        if(local_unknown_real >= 0)
        {
          // Add the Helmholtz term to the Jacobian
          jacobian(local_eqn_real, local_unknown_real) +=
            k_squared_jacobian.real()*psi(l2)*test(l)*hang_weight*hang_weight2;
          
          // Add the Laplace term to the Jacobian
          for (unsigned k=0;k<2;k++)
          {
            jacobian(local_eqn_real,local_unknown_real) +=
              laplace_jacobian[k].real()*dpsidx(l2,k)*dtestdx(l,k)
              *hang_weight*hang_weight2;
          }
        }

        //If at a non-zero degree of freedom add in the entry
        if(local_unknown_imag >= 0)
        {
              // Add the Helmholtz term to the Jacobian
              jacobian(local_eqn_real, local_unknown_imag) +=
                -k_squared_jacobian.imag()*psi(l2)*test(l)
                *hang_weight*hang_weight2;
              
              // Add the Laplace term to the Jacobian
              for (unsigned k=0;k<2;k++)
              {
                jacobian(local_eqn_real,local_unknown_imag) +=
                  -laplace_jacobian[k].imag()*dpsidx(l2,k)*dtestdx(l,k)
                  *hang_weight*hang_weight2;
              }
        }
        }
       }
      }
     }
       
     // Second, compute the imaginary part contribution
     //------------------------------------------------
       
     /*IF it's not a boundary condition*/
     if(local_eqn_imag >= 0)
     {
        // Add the Helmholtz term to the residuals
        residuals[local_eqn_imag]+=k_squared_residual.imag()*test(l)*hang_weight;

        // Add the Laplace term to the residuals
        for (unsigned k=0;k<2;k++)
        {
          residuals[local_eqn_imag]+=laplace_residual[k].imag()*dtestdx(l,k)*hang_weight;
        }
         
      // Calculate the jacobian
      //-----------------------
      if(flag)
      {
       //Local variables to store the number of master nodes
       //and the weights associated with each hanging node
       unsigned n_master2=1; double hang_weight2=1.0;
           
       //Loop over the nodes for the variables
       for(unsigned l2=0;l2<n_node;l2++)
       {
        //Local bool (is the node hanging)
        bool is_node2_hanging = this->node_pt(l2)->is_hanging();
                    
        //If the node is hanging, get the number of master nodes
        if(is_node2_hanging)
        {
         hang_info2_pt = this->node_pt(l2)->hanging_pt();
         n_master2 = hang_info2_pt->nmaster();
        }
        //Otherwise there is one master node, the node itself
        else
        {
         n_master2 = 1;
        }
                    
        //Loop over the master nodes
        for(unsigned m2=0;m2<n_master2;m2++)
        {
         //Get the local unknown and weight
         //If the node is hanging
         if(is_node2_hanging)
         {
          //Read out the local unknown from the master node
          local_unknown_real = 
           this->local_hang_eqn(hang_info2_pt->master_node_pt(m2),
           this->u_index_helmholtz().real());
          local_unknown_imag = 
           this->local_hang_eqn(hang_info2_pt->master_node_pt(m2),
           this->u_index_helmholtz().imag());
                        
          //Read out the hanging weight from the master node
          hang_weight2 = hang_info2_pt->master_weight(m2);
         }
         //If the node is not hanging
         else
         {
          //The local unknown number comes from the node
          local_unknown_real = 
           this->nodal_local_eqn(l2,this->u_index_helmholtz().real());
                        
          local_unknown_imag = 
           this->nodal_local_eqn(l2,this->u_index_helmholtz().imag());
                        
          //The hang weight is one
          hang_weight2 = 1.0;
         }
                      
         //If at a non-zero degree of freedom add in the entry
         if(local_unknown_imag >= 0)
         {
              // Add the Helmholtz term to the Jacobian
              jacobian(local_eqn_imag, local_unknown_imag) +=
                k_squared_jacobian.real()*psi(l2)*test(l)*hang_weight*hang_weight2;
              
              // Add the Laplace term to the Jacobian
              for (unsigned k=0;k<2;k++)
              {
                jacobian(local_eqn_imag,local_unknown_imag) +=
                  laplace_jacobian[k].real()*dpsidx(l2,k)*dtestdx(l,k)*hang_weight*hang_weight2;
              }
         }
         if(local_unknown_real >= 0)
         {
              // Add the Helmholtz term to the Jacobian
              jacobian(local_eqn_imag, local_unknown_real) +=
                k_squared_jacobian.imag()*psi(l2)*test(l)*hang_weight*hang_weight2;
              
              // Add the Laplace term to the Jacobian
              for (unsigned k=0;k<2;k++)
              {
                jacobian(local_eqn_imag,local_unknown_real) +=
                  laplace_jacobian[k].imag()*dpsidx(l2,k)*dtestdx(l,k)*hang_weight*hang_weight2;
              }
         }
        }
       }
      }
     }
    }
   }

  } // End of loop over integration points
 }






//====================================================================
// Force build of templates
//====================================================================
 template class RefineableQPMLHelmholtzElement<1,2>;
 template class RefineableQPMLHelmholtzElement<1,3>;
 template class RefineableQPMLHelmholtzElement<1,4>;

 template class RefineableQPMLHelmholtzElement<2,2>;
 template class RefineableQPMLHelmholtzElement<2,3>;
 template class RefineableQPMLHelmholtzElement<2,4>;

 template class RefineableQPMLHelmholtzElement<3,2>;
 template class RefineableQPMLHelmholtzElement<3,3>;
 template class RefineableQPMLHelmholtzElement<3,4>;

}
