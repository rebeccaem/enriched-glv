#include "dynamics_info.h"

//Constructor
dynamics_info::dynamics_info(
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
    std::vector<double> & data)
:
  N_S(n_S),
  N_s(n_s),
  N_phis_cal(n_phis_cal),
  N_phis_val(n_phis_val),
  N_times(n_times),
  Inad_type(inad_type),
  Params_factor(params_factor),
  Delta_t(delta_t),
  Vectorb(vectorb),
  MatrixA(matrixA),
  Deltas(deltas),
  Derivatives(derivatives),
  Data(data)
{
}

//Destructor
dynamics_info::~dynamics_info()
{
}
