/*-------------------------------------------------------------------
 *
 * Copyright (C) 2012 The PECOS Development Team
 *
 * Please see http://pecos.ices.utexas.edu for more information.
 *
 * This file is part of the QUESO Library (Quantification of Uncertainty
 * for Estimation, Simulation and Optimization).
 *
 * QUESO is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * QUESO is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with QUESO. If not, see <http://www.gnu.org/licenses/>.
 *
 *-------------------------------------------------------------------
 *
 */
 /*------------------------------------------------------------------
 * This file contains the code for the user defined likelihood data 
 * class and the user defined likelihood routine.
 *-----------------------------------------------------------------*/

#include "likelihood.h"
#include "dynamics_info.h"
#include "model.h"
#include <cmath>
#include <stdio.h>
#include <fstream>

// Constructor
likelihoodRoutine_Data::likelihoodRoutine_Data(
    const QUESO::BaseEnvironment& env,
    const std::vector<double> & phis, 
    const std::vector<double> & times, 
    const std::vector<double> & pops,
    double & var,
    dynamics_info * dynInfo)
: m_env(&env),
  m_phis(phis),
  m_times(times),
  m_pops(pops),
  m_var(var),
  m_dynMain(dynInfo)
{
}

// Destructor
likelihoodRoutine_Data::~likelihoodRoutine_Data()
{
}

//------------------------------------------------------
// The user defined likelihood routine
//------------------------------------------------------

double likelihoodRoutine(
  const QUESO::GslVector& paramValues,
  const QUESO::GslVector* paramDirection,
  const void*             functionDataPtr,
  QUESO::GslVector*       gradVector,
  QUESO::GslMatrix*       hessianMatrix,
  QUESO::GslVector*       hessianEffect)
{
  const QUESO::BaseEnvironment& env = *(((likelihoodRoutine_Data*) functionDataPtr)->m_env);
    
  if (paramDirection && functionDataPtr && gradVector && hessianMatrix && hessianEffect) 
  {
    // Just to eliminate INTEL compiler warnings
  }
  
  env.subComm().Barrier(); 
  
  // Compute likelihood 
  // get data from likelihood data structure
  const std::vector<double>&  phis
    = ((likelihoodRoutine_Data *) functionDataPtr)->m_phis;
  const std::vector<double>&  times
    = ((likelihoodRoutine_Data *) functionDataPtr)->m_times;
  const std::vector<double>&  pops
    = ((likelihoodRoutine_Data *) functionDataPtr)->m_pops;
  double & var
    = ((likelihoodRoutine_Data *) functionDataPtr)->m_var;
  dynamics_info *        dyn
    = ((likelihoodRoutine_Data *) functionDataPtr)->m_dynMain;

  const unsigned int n_S = dyn->N_S;          //the number of species included in the model
  const unsigned int n_s = dyn->N_s;          //the number of species included in the model
  const unsigned int n_times = dyn->N_times;  //the number of time points in every time series of data
  const unsigned int n_phis_cal = dyn->N_phis_cal;    //the number of initial conditions
  const unsigned int n_params = dyn->Params_factor * n_s;      //the number of parameters to be calibrated

  /* double variance = .001; */

  //set up lambda vector for loop, right now just one
  /* std::vector<double> phiPoints(n_phis,0.); */

  /* for (unsigned int j = 0; j < n_phis; j++){ */
  /*   phiPoints[j] = phis[n_times * j]; */
  /* } */

  std::vector<double> timePoints(n_times,0.);
  for (unsigned int i = 0; i < n_times; i++){
    timePoints[i] = times[i*n_S];
    /* std::cout << "timePoints[i] = " << timePoints[i]; */
  }

  //init initial conditions
  std::vector<double> initial_conditions(n_s, 0.0);

  //return all times of n_s species
  std::vector<double> returnValues(n_times * 2 * n_s, 0.);

  double misfitValue = 0.;
  double diff = 0.;

  /* for (unsigned int i = 0; i < n_params; i++){  dyn->Deltas[i] = -std::exp(paramValues[i]); } */
  // folded normal
  for (unsigned int i = 0; i < n_params; i++){  
      /* dyn->Deltas[i] = -std::abs(paramValues[i]); */ 
      dyn->Deltas[i] = paramValues[i]; 
      /* std::cout << "dyn->Deltas[i] = " << dyn->Deltas[i]<< "\n"; */
  }
  /* for (unsigned int i = 0; i < n_params; i++){  dyn->Deltas[i] = 0.; } */

  /* int count = 0; */
  for (unsigned int i = 0; i < n_phis_cal; i++) {
    //set initial conditions
    for (unsigned int ic =0; ic < n_s; ic++) {
      initial_conditions[ic] = pops[(i*n_times*n_S) + ic];}
      //set data for memory method
    for (unsigned int t = 0; t < n_times; t++) {
        for (unsigned int sp = 0; sp < n_s; sp++){
          dyn->Data[t*n_s + sp] = pops[i*n_times*n_S + t*n_S + sp];}
          //std::cout << "pops = " << pops[i*n_times*n_S + t*n_S + sp] << "\n";
          //std::cout << "index = " << (t*n_s + sp) << "\n";
          //std::cout << "data = " << dyn->Data[t*n_s + sp] << "\n";}
    }
      /* std::cout << "initial_conditions[ic] = " << pops[(i*n_times*n_S) + ic] << "\n"; */
    
    try
     {
      glvComputeModel(initial_conditions,timePoints,dyn,returnValues);
//      std::cout << "Finished compute model" << std::endl;
      for (unsigned int j = 0; j < n_times; j++){
        for (unsigned int k = 0; k < n_s; k++){
          diff = (returnValues[2*n_s * j + k] - pops[(n_times*n_S) * i + n_S * j + k]);
          //std::cout<<"likelihood: " << (n_tICs*n_times*n_species) * i + (n_times*n_species) *ii + n_species*j + k << "\n";
          misfitValue += diff * diff / var;
          /* std::cout << "misfit = " << misfitValue << "\n"; */
          /* count += 1; */
        }
      }
     } catch( int exception )
     {
      misfitValue = 1000000;
     }
   }

  /* std::cout << "likelihood count = " << count << "\n"; */
  /* std::cout << " the misfit is " << misfitValue << std::endl; */
  return (-0.5 * misfitValue);
}
