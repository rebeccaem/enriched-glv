/*-------------------------------------------------------------------
 * Brief description of this file: 
 * testing the model output of the reduced 5 rxn mechanism 
 *-----------------------------------------------------------------*/

#include <fstream>
#include <iostream>
#include <random>
#include "model.h"
#include "dynamics_info.h"

int main()
{

  unsigned int n_S;  //number of species in detailed model
  unsigned int n_s;  //number of species in reduced model
  unsigned int n_phis_cal; //= 1;  //number of initial conditions
  unsigned int n_phis_val; //= 1;  //number of initial conditions
  unsigned int n_times; // = 100;  //number of data points in every time series
  double var;           // the variance on each data point
  unsigned int inad_type;

  double finalTime = 1.;

  //read in info file
  FILE *infoFile;
  infoFile = fopen("./inputs/info.txt","r");
  if(fscanf(infoFile,"%u %u %u %u %u %lf %u", &n_S, &n_s, &n_phis_cal, &n_phis_val, &n_times, &var, &inad_type)){};
  fclose(infoFile);

  std::cout << "var = " << var << "\n";
  //use detailed model, i.e. n_s = n_S, what
  /* n_s = n_S; */
  //read matrix data into temporary vector
  FILE *matrixFile;
  matrixFile = fopen("./inputs/matrix.txt","r");
  std::vector<double> tmpMatrix(n_S * (n_S + 1), 0.);
  double tmpMat;
  int numLines = 0;
  while (fscanf(matrixFile,"%lf", &tmpMat) != EOF) {
    tmpMatrix[numLines]  = tmpMat;
    numLines++;
  }

  //now place in b vector and A matrix
  std::vector<double> vectorb(n_S, 0.);
  for (unsigned int i=0; i<n_S; i++){vectorb[i] = tmpMatrix[i];}
  //reduced
  std::vector<double> vectorbRed(n_s, 0.);
  for (unsigned int i=0; i<n_s; i++){
      vectorbRed[i] = tmpMatrix[i];
          /* std::cout << "i = " << i << "and vector = " << vectorbRed[i] << "\n"; */
  }

  std::vector< std::vector<double> > matrixA(n_S,std::vector<double>(n_S, 0.));
  for (unsigned int i=0; i<n_S; i++){
      for (unsigned int j=0; j<n_S; j++){
          matrixA[i][j] = tmpMatrix[(i+1)*n_S + j];
      }
  }
  //reduced
  std::vector< std::vector<double> > matrixARed(n_s,std::vector<double>(n_s, 0.));
  for (unsigned int i=0; i<n_s; i++){
      for (unsigned int j=0; j<n_s; j++){
          matrixARed[i][j] = tmpMatrix[(i+1)*n_S + j];
          /* std::cout << "i,j = " << i << "," << j << "and matrix = " << matrixARed[i][j] << "\n"; */
      }
  }

  unsigned int params_factor;
  if( inad_type == 1 ) { params_factor = 2;}
  if( inad_type == 2 ) { params_factor = 2;}
  if( inad_type == 10 ) { params_factor = 2;}
  if( inad_type == 11 ) { params_factor = 2;}
  unsigned int n_delta = params_factor*n_S;         //the model inadequacy terms, should be set to zero
  unsigned int n_params = n_delta;      //no hyperparameters, for now
  //reduced
  unsigned int n_deltaRed = params_factor*n_s;         //the model inadequacy terms, should be set to zero
  unsigned int n_paramsRed = n_deltaRed;      //no hyperparameters, for now

  //dummy vector for now
  std::vector<double> delta(n_delta,0.);
  std::vector<double> deltaRed(n_deltaRed,0.);
  std::vector<double> derivatives(n_S,0.);
  std::vector<double> derivativesRed(n_s,0.);
  /* dynamics_info dynMain(n_S, n_s, n_phis, n_times, vectorb, matrixA, delta); */
  unsigned int n_phis_total = n_phis_cal + n_phis_val;

  std::vector<double> timePoints(n_times,0.);
  for (unsigned int i = 0; i < timePoints.size(); i++){
    timePoints[i] = (i) * finalTime/n_times;
    /* timePoints[i] = (i) * 2.5e-5 + 3.0e-5; */
    /* timePoints[i] = (i+1) * 5.0e-6; */
  }

  const double delta_t = timePoints[1] - timePoints[0];
  std::vector<double> data(n_phis_total*n_s*n_times, 0.); //don't have this yet but need placeholder

  dynamics_info dynMain(n_S, n_S, n_phis_cal, n_phis_val, n_times, inad_type, params_factor, delta_t, vectorb, matrixA, delta, derivatives, data);
  dynamics_info dynRed(n_s, n_s, n_phis_cal, n_phis_val, n_times, inad_type, params_factor, delta_t, vectorbRed, matrixARed, deltaRed, derivativesRed, data);
  
  std::vector<double> phiPoints(n_phis_total,0.);
  for (unsigned int i = 0; i < n_phis_total; i++){
  phiPoints[i] = i; //phiPoints[1] = 1.0; phiPoints[2] = 1.1;
  }
  std::vector<double> initialValues(n_S,0.); 
  //double timePoint = 2.0e-5;
  std::vector<double> initialValuesRed(n_s,0.); 


  //open file now to write data
  std::ofstream datafile[2];
  /* datafile.open ("datafile-det-all-dense.txt"); */
  datafile[0].open ("./inputs/datafile.txt");
  datafile[1].open ("./inputs/datafile-reduced.txt");

  double IC;
  std::default_random_engine generator;
  std::vector<double> set_ICs(n_phis_total*n_S,0.);
  set_ICs[0] = 1.;
  set_ICs[1] = 1.;
  set_ICs[2] = 1.;
  set_ICs[3] = 2.;
  set_ICs[4] = 2.;
  set_ICs[5] = 1.;
  for (unsigned int l = 0; l < n_phis_total; l++){
    //use random values to set initial conditions ???????????????????????????????????
    //better way to set these????????????????????????????????????????????????????????
    std::normal_distribution<double> distributionIC(0.0,std::sqrt(1));
    for (unsigned int i = 0; i < n_S; i++){
      IC = distributionIC(generator);
      /* IC = 0; //if want to set ICs to 1*/
      initialValues[i] = std::exp(IC);
//      initialValues[i] = set_ICs[l*n_S+i]; //hard code ICs
      /* initialValues[i] = 1; //start everything at 1 for now */
    }
    for (unsigned int i = 0; i < n_s; i++){
      initialValuesRed[i] = initialValues[i];
    }

    //return 10 time points of all species and temp
    std::vector<double> returnValues(2*n_S*timePoints.size(),0.);
    std::vector<double> returnValuesRed(2*n_s*timePoints.size(),0.);

    glvComputeModel(initialValues,timePoints,&dynMain,returnValues);
    std::cout << "helloooooooooooooooooooooo?\n";
    glvComputeModel(initialValuesRed,timePoints,&dynRed,returnValuesRed);
   
    //create measurement error
    std::normal_distribution<double> distribution(0.0,std::sqrt(var));
    //TODO increase this appropriately for temperature measurements

    //right now not adding any error
    for (unsigned int i = 0; i < timePoints.size(); i++){
        /* std::cout << "i = " << i << "\n"; */
      
        for (unsigned int j = 0; j < n_S; j++){
      //    double error = distribution(generator);
          double measurement = returnValues[2*n_S*i + j];//take out error for now: + error;
          double m_deriv = returnValues[2*n_S*i + j + n_S];//take out error for now: + error;
          datafile[0] << phiPoints[l] << " " << timePoints[i] << " " << measurement << " " << m_deriv << "\n";
          //if (measurement>0){datafile[0] << phiPoints[l] << " " << timePoints[i] << " " << measurement << "\n";}
          //else {datafile[0] << phiPoints[l] << " " << timePoints[i] << " " <<  0.0 << "\n";}
      } //datafile << "\n";
        for (unsigned int j = 0; j < n_s; j++){
          double measurementRed = returnValuesRed[2*n_s*i + j];
          double mRed_deriv = returnValuesRed[2*n_s*i + j + n_s];
          datafile[1] << phiPoints[l] << " " << timePoints[i] << " " << measurementRed << " " << mRed_deriv << "\n";
        }
    }
  }
  return 0;
}
