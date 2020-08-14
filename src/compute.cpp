/*-------------------------------------------------------------------
 * This file is divided in two parts:
 * - the first one handles the statistical inverse problem (SIP) for estimating
 *   the reaction rates 'k_i', i=1,...,5, where k = A * T^beta * exp(-E/RT)
 *   and the stochastic operator hyperparameters
 * - the second part handles the statistical forward problem (SFP) for
 *   predicting the QoI which has not yet been defined
 *-----------------------------------------------------------------*/
/* #include <fstream> */
#include "compute.h"
#include "likelihood.h"
#include "qoi.h"
#include "dynamics_info.h"
//queso
#include <queso/GslVector.h>
#include <queso/GslMatrix.h>
#include <queso/GenericScalarFunction.h>
#include <queso/GenericVectorFunction.h>
#include <queso/GenericVectorRV.h>
#include <queso/GaussianVectorRV.h>
#include <queso/LogNormalVectorRV.h>
#include <queso/UniformVectorRV.h>
#include <queso/StatisticalInverseProblem.h>
#include <queso/StatisticalForwardProblem.h>
#include <sys/time.h>
#include <cmath>
#include <vector>

void computeParams(const QUESO::FullEnvironment& env) {
  struct timeval timevalNow;
  
  gettimeofday(&timevalNow, NULL);
  if (env.fullRank() == 0) {
    std::cout << "\nBeginning run of 'GLV plus' example at "
              << ctime(&timevalNow.tv_sec)
              << "\n my fullRank = "         << env.fullRank()
              << "\n my subEnvironmentId = " << env.subId()
              << "\n my subRank = "          << env.subRank()
              << "\n my interRank = "        << env.inter0Rank()
               << std::endl << std::endl;
  }

  // Just examples of possible calls
  if ((env.subDisplayFile()       ) && 
      (env.displayVerbosity() >= 2)) {
    *env.subDisplayFile() << "Beginning run of 'GLV plus' example at "
                          << ctime(&timevalNow.tv_sec)
                          << std::endl;
  }
  env.fullComm().Barrier();
  env.subComm().Barrier();  // Just an example of a possible call
  
  //================================================================
  // Statistical inverse problem (SIP): find posterior PDF for k's
  //================================================================
  gettimeofday(&timevalNow, NULL);
  if (env.fullRank() == 0) {
    std::cout << "Beginning 'SIP -> all parameters estimation' at "
              << ctime(&timevalNow.tv_sec)
              << std::endl;
  }

  //------------------------------------------------------
  // SIP Step 0 of 6: Read in the data
  //------------------------------------------------------
  unsigned int n_S;  //number of species in detailed model
  unsigned int n_s;  //number of species in reduced model
  unsigned int n_phis_cal;  //number of initial conditions
  unsigned int n_phis_val;  //number of initial conditions
  unsigned int n_times;  //number of data points in every time series
  double var;            //variance in the data
  unsigned int inad_type;  //type of inadequacy model

  //read in info file
  FILE *infoFile;
  infoFile = fopen("./inputs/info.txt","r");
  if(fscanf(infoFile,"%u %u %u %u %u %lf %u", &n_S, &n_s, &n_phis_cal, &n_phis_val, &n_times, &var, &inad_type)){};
  fclose(infoFile);

  std::cout << "n_S = " << n_S << "\n";
  std::cout << "n_s = " << n_s << "\n";
  std::cout << "n_phis_cal = " << n_phis_cal << "\n";
  std::cout << "n_phis_val = " << n_phis_val << "\n";
  std::cout << "n_times = " << n_times << "\n";
  std::cout << "variance = " << var << "\n";
  std::cout << "inadequacy type = " << inad_type << "\n\n";

  unsigned int n_phis_tot = n_phis_cal + n_phis_val;
  unsigned int params_factor;
  if( inad_type == 0 ) { params_factor = 1;}
  if( inad_type == 1 ) { params_factor = 2;}
  if( inad_type == 2 ) { params_factor = 2;}
  if( inad_type == 3 ) { params_factor = 4;}
  if( inad_type == 4 ) { params_factor = 2 * n_s;}
  if( inad_type == 10 ) { params_factor = 2;}
  if( inad_type == 11 ) { params_factor = 2;}
  unsigned int n_delta = params_factor*n_s;         //the model inadequacy terms
  unsigned int n_params = n_delta;      //no hyperparameters, for now

  //read matrix data into temporary vector
  FILE *matrixFile;
  matrixFile = fopen("./inputs/matrix.txt","r");
  //this file is in format of one entry per line
  std::vector<double> tmpMatrix(n_S * (n_S + 1), 0.);
  double tmpMat;
  int numLines = 0;
  while (fscanf(matrixFile,"%lf", &tmpMat) != EOF) {
    tmpMatrix[numLines]  = tmpMat;
    numLines++;
  }
  fclose(matrixFile);
  /* std::cout << "numlines = " << numLines << "\n"; */
  //now place in b vector and A matrix
  std::vector<double> vectorb(n_S, 0.);
  for (unsigned int i=0; i<n_S; i++){vectorb[i] = tmpMatrix[i];}

  std::vector< std::vector<double> > matrixA(n_S,std::vector<double>(n_S, 0.));
  for (unsigned int i=0; i<n_S; i++){
      for (unsigned int j=0; j<n_S; j++){
          matrixA[i][j] = tmpMatrix[(i+1)*n_S + j];
      }
  }

  //read in data points
  double tmpPhis;
  double tmpTimes;
  double tmpx;
  double tmpdx;
  numLines = 0;

  FILE *dataFile;
  dataFile = fopen("./inputs/datafile.txt","r");

  std::vector<double> phis(n_phis_tot * n_times * n_S, 0.);
  std::vector<double> times(n_phis_tot * n_times * n_S, 0.);
  std::vector<double> pops(n_phis_tot * n_times * n_S, 0.);
  std::vector<double> derivs(n_phis_tot * n_times * n_S, 0.);
  while (fscanf(dataFile,"%lf %lf %lf %lf", &tmpPhis, &tmpTimes, &tmpx, &tmpdx) != EOF) {
    phis[numLines]     = tmpPhis;
    times[numLines]    = tmpTimes;
    pops[numLines]     = tmpx;
    derivs[numLines]   = tmpdx;
    numLines++;
  }
  fclose(dataFile);
  std::cout << "The number of data points is " << n_s*n_phis_cal*n_times << "\n\n";

  const double delta_t = times[n_S] - times[0]; //set time interval this way TODO: improve this?

  //condense phi and times vectors down so no repeats
  /* std::vector<double> phis(n_phis, 0.); */
  /* for (unsigned int i=0; i<n_phis; i++){ */
  /*     phis[i] = phisAll[n_times*i]; */
  /* } */
  /* std::vector<double> times(n_times, 0.); */
  /* for (unsigned int i=0; i<n_times; i++){ */
  /*     times[i] = timesAll[i]; */
  /* } */

  //create dummy vector to be filled with params inside likelihood
  std::vector<double> queso_params(n_params, 0.0); //alpha_{0,1,2}

  //------------------------------------------------------
  // SIP Step 1 of 6: Instantiate the parameter space
  //------------------------------------------------------
  QUESO::VectorSpace<QUESO::GslVector,QUESO::GslMatrix> paramSpace(env, "param_", n_params, NULL);

  //------------------------------------------------------
  // SIP Step 2 of 6: Instantiate the parameter domain
  //------------------------------------------------------
  QUESO::GslVector paramMinValues(paramSpace.zeroVector());
  QUESO::GslVector paramMaxValues(paramSpace.zeroVector());
  //mean of xi
  for (unsigned int i=0; i<n_params; ++i){
    paramMinValues[i] = -10;
    paramMaxValues[i] = 0;
}
  /* //variance of xi */
  /* for (unsigned int i=n_xi; i<2*n_xi; ++i){ */
  /*   paramMinValues[i] = -INFINITY; */
  /*   paramMaxValues[i] = INFINITY;} */
  /* //xi */
  /* for (unsigned int i=2*n_xi; i<3*n_xi; ++i){ */
  /*   paramMinValues[i] = -INFINITY; */
  /*   paramMaxValues[i] = INFINITY;} */
  //mean of a_0 and b_0

  QUESO::BoxSubset<QUESO::GslVector,QUESO::GslMatrix>
    paramDomain("param_", paramSpace, paramMinValues, paramMaxValues);

  //placeholder vector to return derivative values during ODE solve
  std::vector<double> derivatives(n_s, 0.);
  //placeholder for data to pass to ODE solver
  std::vector<double> data(n_s*n_times, 0.);
  dynamics_info dynMain(n_S, n_s, n_phis_cal, n_phis_val, n_times, inad_type, params_factor, delta_t, vectorb, matrixA, queso_params, derivatives, data);

  //------------------------------------------------------
  // SIP Step 3 of 6: Instantiate the likelihood function 
  // object to be used by QUESO.
  //------------------------------------------------------
  likelihoodRoutine_Data likelihoodRoutine_Data1(env, phis, times, pops, var, &dynMain);

  QUESO::GenericScalarFunction<>
    likelihoodFunctionObj(
        "like_",
			  paramDomain,
			  likelihoodRoutine,
        static_cast<void *> (&likelihoodRoutine_Data1),
			  true); // the routine computes [ln(function)]
    
  //------------------------------------------------------
  // SIP Step 4 of 6: Define the prior RV
  //------------------------------------------------------
  // UNIFORM PRIOR
  QUESO::UniformVectorRV<QUESO::GslVector,QUESO::GslMatrix>
    priorRv("prior_",paramDomain);

  // NORMAL PRIOR
  //QUESO::GslVector meanVec(paramSpace.zeroVector());
  //QUESO::GslVector diagVec(paramSpace.zeroVector());
  //for (unsigned int i = 0; i < n_params; i++) meanVec[i] = 1.;
  //for (unsigned int i = 0; i < n_params; i++) diagVec[i] = 1.;
  //QUESO::GslMatrix covMatrix(diagVec);

  //QUESO::GaussianVectorRV<QUESO::GslVector,QUESO::GslMatrix>
    //priorRv("prior_",paramDomain,meanVec,covMatrix);

  //------------------------------------------------------
  // SIP Step 5 of 6: Instantiate the inverse problem
  //------------------------------------------------------
  QUESO::GenericVectorRV<>
    postTotal("post_",  // Extra prefix before the default "rv_" prefix
           paramSpace);
        
  QUESO::StatisticalInverseProblem<>
    ip("",          // No extra prefix before the default "ip_" prefix
       NULL,
       priorRv,
       likelihoodFunctionObj,
       postTotal);

  //------------------------------------------------------
  // SIP Step 6 of 6: Solve the inverse problem, that is,
  // set the 'pdf' and the 'realizer' of the posterior RV //------------------------------------------------------
  std::cout << "Solving the SIP with Multi-Level Metropolis Hastings" 
	    << std::endl << std::endl;  

  //The following is set if use ip.solveWithBayesMetropolisHastings
  QUESO::GslVector paramInitials(paramSpace.zeroVector());
  for (unsigned int i = 0; i < n_params; i++) {paramInitials[i] = -1;}
   //
//priorRv.realizer().realization(paramInitials);

  QUESO::GslVector diagVec(paramSpace.zeroVector());
  QUESO::GslMatrix proposalCovMatrix(diagVec);
  for (unsigned int i = 0; i < n_params; i++) proposalCovMatrix(i,i) = 1;

  ip.solveWithBayesMetropolisHastings(NULL, paramInitials, &proposalCovMatrix);

  /* ip.seedWithMAPEstimator(); */
  /* ip.solveWithBayesMetropolisHastings(); */

  /* ip.solveWithBayesMLSampling(); */

  //================================================================
  // Statistical forward problem (SFP)
  //================================================================
  gettimeofday(&timevalNow, NULL);
  std::cout << "Beginning 'SFP -> Undecided QoI' at " 
            << ctime(&timevalNow.tv_sec)
            << std::endl;

  //------------------------------------------------------
  // SFP Step 1 of 6: Instantiate the parameter *and* qoi spaces. 
  // SFP input RV = FIP posterior RV, so SFP parameter space
  // has been already defined.
  //------------------------------------------------------
  // Return also model discrepancy piece so *2
  QUESO::VectorSpace<QUESO::GslVector,QUESO::GslMatrix> qoiSpace(env, "qoi_", n_phis_tot * n_times * n_s * 2, NULL);

  //------------------------------------------------------
  // SFP Step 2 of 6: Instantiate the parameter domain 
  //------------------------------------------------------
  
  // Not necessary because input RV of the SFP = output RV of SIP. 
  // Thus, the parameter domain has been already defined.
  
  //------------------------------------------------------ 
  // SFP Step 3 of 6: Instantiate the qoi function object 
  // to be used by QUESO.
  //------------------------------------------------------
  qoiRoutine_Data qoiRoutine_Data(env, phis, times, pops, &dynMain);
  
  QUESO::GenericVectorFunction<QUESO::GslVector,QUESO::GslMatrix,QUESO::GslVector,QUESO::GslMatrix>
    qoiFunctionObj("qoi_",
                   paramDomain,
                   qoiSpace,
                   qoiRoutine,
                   (void *) &qoiRoutine_Data);
      
  //------------------------------------------------------
  // SFP Step 4 of 6: Define the input RV
  //------------------------------------------------------
  
  // Not necessary because input RV of SFP = output RV of SIP 
  // (postRv).
      
  //------------------------------------------------------
  // SFP Step 5 of 6: Instantiate the forward problem
  //------------------------------------------------------
  QUESO::GenericVectorRV<QUESO::GslVector,QUESO::GslMatrix> qoiRv("qoi_", qoiSpace);
  
  QUESO::StatisticalForwardProblem<QUESO::GslVector,QUESO::GslMatrix,QUESO::GslVector,QUESO::GslMatrix>
    fp("",
       NULL,
       postTotal,
       qoiFunctionObj,
       qoiRv);

  //------------------------------------------------------
  // SFP Step 6 of 6: Solve the forward problem
  //------------------------------------------------------
  std::cout << "Solving the SFP with Monte Carlo" 
            << std::endl << std::endl;  
        std::cout << "hellooooooooooooooooooooooooooooooooooo\n";
  fp.solveWithMonteCarlo(NULL);

  //------------------------------------------------------
  gettimeofday(&timevalNow, NULL);
  if ((env.subDisplayFile()       ) && 
      (env.displayVerbosity() >= 2)) {
    *env.subDisplayFile() << "Ending run of 'GLV plus' example at "
                          << ctime(&timevalNow.tv_sec)
                          << std::endl;
  }
  if (env.fullRank() == 0) {
    std::cout << "Ending run of 'GLV plus' example at "
              << ctime(&timevalNow.tv_sec)
              << std::endl;
  }

  return;
}
