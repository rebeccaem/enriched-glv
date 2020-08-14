#ifndef DYNAMICS_INFO_H
#define DYNAMICS_INFO_H

//c++
#include <vector>

// define struct that holds all reaction info, except params
struct dynamics_info { dynamics_info(
  const unsigned int & n_S,
  const unsigned int & n_s,
  const unsigned int & n_phis_cal,
  const unsigned int & n_phis_val,
  const unsigned int & n_times,
  const unsigned int & inad_type,
  const unsigned int & params_factor,
  const double       & delta_t,
  std::vector<double> & vectorb,
  std::vector< std::vector<double> > & matrixA,
  std::vector<double> & deltas,
  std::vector<double> & derivatives,
  std::vector<double> & data);
 ~dynamics_info();

  const unsigned int & N_S;
  const unsigned int & N_s;
  const unsigned int & N_phis_cal;
  const unsigned int & N_phis_val;
  const unsigned int & N_times;
  const unsigned int & Inad_type;
  const unsigned int & Params_factor;
  const double       & Delta_t;
  std::vector<double> & Vectorb;
  std::vector< std::vector<double> > & MatrixA;
  std::vector<double> & Deltas;
  std::vector<double> & Derivatives;
  std::vector<double> & Data;
};
#endif
