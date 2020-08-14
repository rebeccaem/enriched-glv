 /*------------------------------------------------------------------
 * Brief description of this file: 
 * 
 * This file contains the code for the forward model.
 *-----------------------------------------------------------------*/

#include "model.h"
/* #include "dynamics_info.h" */
#include <cmath>
#include <vector>
#include <iostream>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_odeiv2.h>
#include <gsl/gsl_integration.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_roots.h>
#include <assert.h>
#include "eigen3/Eigen/Dense"
/* //antioch */
/* #include <antioch/kinetics_evaluator.h> */
/* #include <antioch/cea_evaluator.h> */

#ifndef __EPS_ABS
#define __EPS_ABS 1e-8
#endif
#ifndef __EPS_REL
#define __EPS_REL 1e-8
#endif

//first define the function for the ODE solve
int glvFunction( double t,
                 const double Y[], // TESTING: REMOVED THE CONST
                 double dYdt[],
                 void* params)
{
  //here, params is sending the function all the reaction info
  dynamics_info dyn = *(dynamics_info *) params;

  const unsigned int n_s = dyn.N_s;
  const unsigned int inad_type = dyn.Inad_type;
  const unsigned int pf = dyn.Params_factor;
  const unsigned int n_deltas = pf * n_s;
  const double       delta_t = dyn.Delta_t; //time interval, different from deltas above!
  std::vector< std::vector<double> > matrixA = dyn.MatrixA;
  std::vector<double> vectorb = dyn.Vectorb;
  std::vector<double> delta = dyn.Deltas;
  std::vector<double> data = dyn.Data;

  std::vector<double> pops(n_s,0);
  for (unsigned int i = 0; i < n_s; i++){
    pops[i] = Y[i];
    if(pops[i] <= 0){
      pops[i] = 0;
    }
  }

  for (unsigned int i = 0; i < n_s; i++){
      double sum = 0;
      for (unsigned int j = 0; j < n_s; j++){
          sum += matrixA[i][j]*pops[j];
      }
      dYdt[i] = vectorb[i]*pops[i] + sum*pops[i];
//    std::cout << "dYdt["<<i<<"] = " << dYdt[i] <<"\n";
  }

//inadequacy formulation
  if ( inad_type == 0) {
    for (unsigned int i = 0; i < n_s; i++){
      dYdt[i] += 0.;
    }
  }
  else if ( inad_type == 1) {
    for (unsigned int i = 0; i < n_s; i++){
      dYdt[i] += delta[pf*i+0]*pops[i] + delta[pf*i+1]*std::abs(dYdt[i]);
      /* dYdt[i] += delta[pf*i+0] + delta[pf*i+1]*pops[i] + delta[pf*i+2]*std::abs(dYdt[i]); */
    }
  }
  else if ( inad_type == 2) {
    for (unsigned int i = 0; i < n_s; i++){
      dYdt[i] += delta[pf*i+0]*pops[i] + delta[pf*i+1]*std::pow(dYdt[i],2);
      /* dYdt[i] += delta[pf*i+0] + delta[pf*i+1]*pops[i] + delta[pf*i+2]*std::abs(dYdt[i]); */
    }
  }
  else if ( inad_type == 3) {
    for (unsigned int i = 0; i < n_s; i++){
      dYdt[i] += delta[pf*i+0]*pops[i] + delta[pf*i+1]*std::abs(dYdt[i]) + 
          delta[pf*i+2]*std::pow(pops[i],2) + delta[pf*i+3]*std::pow(dYdt[i],2);
    }
  }
  else if ( inad_type == 4) {
    for (unsigned int i = 0; i < n_s; i++){
        for (unsigned int j = 0; j < n_s; j++){
            dYdt[i] += delta[pf*i + 2*j +0]*pops[j] + delta[pf*i + 2*j + 1]*std::abs(dYdt[j]);
            dYdt[i] += delta[pf*i + 2*j +2]*std::pow(pops[j],2) + delta[pf*i + 2*j + 3]*std::pow(dYdt[j],2);
         }
    }
  }
  //memory method
  else if ( inad_type == 10) {
    double trap = 0.0;
    double t_int, t_frac; //integer time and fractional time
    //std::cout<<"time = "<<t<<"\n";
    //std::cout<<"delta_t = "<<delta_t<<"\n";
    t_int = floor(t/delta_t);
    t_frac = fmod(t, delta_t);
    //std::cout<<"fraction of time = "<<t_frac<<"\n";
    //std::cout<<"integer of time = "<<t_int<<"\n";
    double height;
    for (unsigned int i = 0; i < n_s; i++){
//        //TODO: FILL IN HERE!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        for (unsigned int j = 0; j < t_int; j++){
     //       std::cout<<"data[j*n_s+ i] = "<<data[j*n_s+ i]<<"\n";
      //      std::cout<<"data[(j+1)*n_s + i]"<<data[(j+1)*n_s+ i]<<"\n";
           trap  += 0.5*delta_t*(data[j*n_s+ i] + data[(j+1)*n_s + i]);}
      height = data[t_int*n_s + i] + (data[(t_int+1)*n_s + i] -data[t_int*n_s + i])/delta_t *( t_frac);
      //trap += 0.5*delta_t*(height + data[(t_int+1)*n_s + i]);
      //std::cout<<"i = " << i << " and trap = "<<trap<<"\n";
      //std::cout<<"height = "<<height<<"\n";
      dYdt[i] += delta[pf*i+0]*pops[i] + delta[pf*i+1]*trap;
      //dYdt[i] += delta[pf*i+0]*std::pow(pops[i],2) + delta[pf*i+1]*trap*pops[i];
//      //dYdt[i] += delta[pf*i+0] + delta[pf*i+1]*pops[i] + delta[pf*i+2]*std::abs(dYdt[i]);
    }
  }
  //algebraic
  else if ( inad_type == 11) {
    for (unsigned int i = 0; i < n_s; i++){
      dYdt[i] += delta[pf*i+0]*std::pow(pops[i],2) + delta[pf*i+1]*pops[i]*dYdt[i];
      //dYdt[i] += delta[pf*i+0] + delta[pf*i+1]*pops[i] + delta[pf*i+2]*std::abs(dYdt[i]);
    }
  }

  else {
     std::cout << "Choose valid inadequacy type.\n";
  }

  // set derivatives inside dyn
  for (unsigned int i = 0; i < n_s; i++){
    //std::cout << "Y["<<i<<"] = " << Y[i] <<"\n";
    //std::cout << "dYdt["<<i<<"] = " << dYdt[i] <<"\n";
      dyn.Derivatives[i] = dYdt[i];
  }
  
//  for (unsigned int i = n_s; i < 2*n_s; i++){
//      dYdt[i] = 0;
//  }

  return GSL_SUCCESS;
}

//jacobian for ode solve---------------------------------------------
int glvJacobian( double t, 
				const double Y[],
				double *dfdY,
				double dfdt[],
				void* params )
{
  return GSL_SUCCESS;
}

void glvComputeModel(
  std::vector<double>&  initialValues,
  std::vector<double>&  timePoints,
  dynamics_info*        dyn,
  std::vector<double>&  returnValues)
{  
  // Compute model
  // GSL prep
  unsigned int dim = initialValues.size();
  //unsigned int n_s = dim;
  gsl_odeiv2_system sys = { glvFunction, 
			   glvJacobian, 
			   dim, dyn };
  
  double h = 1e-10;    //initial step-size
  gsl_odeiv2_driver * d = gsl_odeiv2_driver_alloc_y_new( &sys, gsl_odeiv2_step_rkf45,h,1e-8,1e-4);   
  // initialize values
  double Y[dim];
  for (unsigned int i = 0; i < dim; ++i){
    Y[i] = initialValues[i];
  };
  double t = 0.0;
  double prevt;
  double sumY;
//  std::cout << "Starting Integration..." << std::endl;
  // GSL integration
  double finalTime;
  for (unsigned int i = 0; i < timePoints.size(); i++){
    finalTime = timePoints[i];
    while (t < finalTime)
      {
     //   prevt = t;
     //   sumY = 0;
        // t and y are updated and placed back into those variables
        int status = gsl_odeiv2_driver_apply( d, // necessary gsl vars
               &t,    // current time
               finalTime,    // maY time
               Y );   // current solution values (at time t)
//        std::cout<<"Time = "<<t<<"\nAfter integration Y values : \n";
//        for (unsigned int i = 0; i <= n_species; i++) std::cout<<Y[i]<<"\n";
        /* std::cout<<"t = "<<t<<"\n"; */
        /* std::cout<<"Y[0] = "<<Y[0]<<"\n"; */
        /* std::cout<<"Y[1] = "<<Y[1]<<"\n"; */
        //std::cout<<"T = "<<Y[n_species]<<"\n";
//        for( i = 0; i<7;i++) sumY+=Y[i];
//        std::cout<<"N = "<<sumY<<"\n";
        // check that the evolution was successful
        #ifdef UQ_FATAL_TEST_MACRO
          UQ_FATAL_TEST_MACRO( status != GSL_SUCCESS,
             0,
             "reduced Lotka-Volterra",
             "The status of GSL integration != GSL_SUCCESS" );
        #else 
          if ( status != GSL_SUCCESS )
          {
     //       std::cout<< "h approx = " << t-prevt<<"\n\n";
            std::cout << "ERROR: status of GSL integration != GSL_SUCCESS" <<
              std::endl;
            assert( status == GSL_SUCCESS );
          }
        #endif
      }
    //  std::cout << " h is " << h << std::endl;
    // save results, right now return values are all the species of the reduced
    // model
    for (unsigned int j = 0; j < dim; j++){
      returnValues[2 * dim * i + j] = Y[j];
      //std::cout << "Y["<<j<<"] = " << Y[j] <<"\n";
    }
    for (unsigned int j = 0; j < dim; j++){
      returnValues[2 * dim * i + (dim + j)] = dyn->Derivatives[j];
      //std::cout << "dYdt["<<j<<"] = " << dyn->Derivatives[j] <<"\n";
    }
  }
  /* std::cout << "x1 = " << Y[0] << std::endl; */
  //std::cout << "O2 = " << Y[1] << std::endl;
  // deallocate memory
  gsl_odeiv2_driver_free( d );
}
