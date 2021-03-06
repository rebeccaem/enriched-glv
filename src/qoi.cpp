/*------------------------------------------------------------------------
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
 *------------------------------------------------------------------------
 *
 */
 /*------------------------------------------------------------------
 * Brief description of this file: 
 * 
 * This file contains the code for the user defined qoi routine.
 *-----------------------------------------------------------------*/

#include "qoi.h"
#include "model.h"
#include "dynamics_info.h"
#include <cmath>
//------------------------------------------------------
/// The actual (user-defined) qoi routine
//------------------------------------------------------
// Constructor
qoiRoutine_Data::qoiRoutine_Data(
    const QUESO::BaseEnvironment& env,
    const std::vector<double> & phis,
    const std::vector<double> & times,
    std::vector<double> & pops,
    dynamics_info * dynInfo)
: m_env(&env),
  m_phis(phis),
  m_times(times),
  m_pops(pops),
  m_dynMain(dynInfo)
{
}

// Destructor
qoiRoutine_Data::~qoiRoutine_Data()
{
}

void
qoiRoutine(
  const QUESO::GslVector&                    paramValues,
  const QUESO::GslVector*                    paramDirection,
  const void*                                functionDataPtr,
        QUESO::GslVector&                    qoiValues,
        QUESO::DistArray<QUESO::GslVector*>* gradVectors,
        QUESO::DistArray<QUESO::GslMatrix*>* hessianMatrices,
        QUESO::DistArray<QUESO::GslVector*>* hessianEffects)
{
  const QUESO::BaseEnvironment& env = paramValues.env();

  if (paramDirection && 
      gradVectors    &&
      hessianEffects &&
      hessianMatrices) {
    // Logic just to avoid warnings from INTEL compiler
  }

  // get data from qoi data structure
  const std::vector<double>&  phis
    = ((qoiRoutine_Data *) functionDataPtr)->m_phis;
  const std::vector<double>&  times
    = ((qoiRoutine_Data *) functionDataPtr)->m_times;
  const std::vector<double>&  pops
    = ((qoiRoutine_Data *) functionDataPtr)->m_pops;
  dynamics_info *       dyn
    = ((qoiRoutine_Data *) functionDataPtr)->m_dynMain;

  const unsigned int n_S = dyn->N_S;          //the number of species included in the model
  const unsigned int n_s = dyn->N_s;          //the number of species included in the model
  const unsigned int n_times = dyn->N_times;  //the number of time points in every time series of data
  const unsigned int n_phis_cal = dyn->N_phis_cal;    //the number of initial conditions
  const unsigned int n_phis_val = dyn->N_phis_val;    //the number of initial conditions
  const unsigned int pf = dyn->Params_factor;    //the number of initial conditions
  const unsigned int n_params = pf * n_s;      //the number of parameters to be calibrated

  unsigned int n_phis_tot = n_phis_cal + n_phis_val;
  //define initial conditions = state + temp
  std::vector<double> initial_conditions(n_s, 0.0);

  //set up lambda vector for loop, right now just one
  /* std::vector<double> phiPoints(n_phis,0.); */
  /* for (unsigned int j = 0; j < n_phis; j++){ */
  /*   phiPoints[j] = phis[n_times * j]; */
  /* } */

  std::vector<double> timePoints(n_times,0.);
  for (unsigned int i = 0; i < n_times; i++){
    timePoints[i] = times[i*n_S];
  }

  //return 10 time points of x1
  std::vector<double> returnValues(n_times * 2 * n_s, 0.);

  /* for (unsigned int i = 0; i < n_params; i++){  dyn->Deltas[i] = -std::exp(paramValues[i]); } */
  //folded normal
  for (unsigned int i = 0; i < n_params; i++){  
      /* dyn->Deltas[i] = -std::abs(paramValues[i]); */
      dyn->Deltas[i] = paramValues[i];
  }
  /* for (unsigned int i = 0; i < n_params; i++){  dyn->Deltas[i] = 0.; } */

  for (unsigned int i = 0; i < n_phis_tot; ++i) {
    //set initial conditions
    for (unsigned int ic =0; ic < n_s; ic++) {
      initial_conditions[ic] = pops[(i*n_times*n_S) + ic];
    }
    glvComputeModel(initial_conditions,timePoints,dyn,returnValues);
    //std::cout<< "qoi: ret val = " <<  returnValues[7 * 9 + 6] << std::endl; 
    for (unsigned int j = 0; j < returnValues.size(); j++){
      /* std::cout << "i = " << i << " and j = " << j << "\n"; */
      qoiValues[(n_times*2*n_s) * i + j] = returnValues[j];
      //std::cout << "returnValues[j] = " << returnValues[j] << "\n\n";
    }
  }

// for (unsigned int i = 0; i < 90; i++) qoiValues[i] = i;
  return;
}
