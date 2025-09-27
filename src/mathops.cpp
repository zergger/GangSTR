#include "mathops.h"

#include <gsl/gsl_math.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>

#include <nlopt.hpp>
#include <iostream>
using namespace std;
//////////
// int count = 0;
double myfunc_test(unsigned n, const double *x, double *grad, void *my_func_data)
{
  // ++count;
    if (grad) {
        grad[0] = 0.0;
        grad[1] = 0.5 / sqrt(x[1]);
    }
    return sqrt(x[1]);
}
typedef struct {
    double a, b;
} my_constraint_data_test;
double myconstraint_test(unsigned n, const double *x, double *grad, void *data)
{
    my_constraint_data_test *d = (my_constraint_data_test *) data;
    double a = d->a, b = d->b;
    if (grad) {
        grad[0] = 3 * a * (a*x[0] + b) * (a*x[0] + b);
        grad[1] = -1.0;
    }
    return ((a*x[0] + b) * (a*x[0] + b) * (a*x[0] + b) - x[1]);
 }
////////////

double TestNLOPT(){
  nlopt::opt opt(nlopt::LD_MMA, 2);
  std::vector<double> lb(2);
  lb[0] = -HUGE_VAL; lb[1] = 0;
  opt.set_lower_bounds(lb);

  opt.set_min_objective(myfunc_test, NULL);

  my_constraint_data_test data[2] = { {2,0}, {-1,1} };
  opt.add_inequality_constraint(myconstraint_test, &data[0], 1e-8);
  opt.add_inequality_constraint(myconstraint_test, &data[1], 1e-8);

  opt.set_xtol_rel(1e-4);

  std::vector<double> xx(2);
  xx[0] = 1.234; xx[1] = 5.678;
  double minf;
  nlopt::result result = opt.optimize(xx, minf);

  cout<<result;
  return 0.0;
}
double dummy_func (const gsl_vector *v, void *params)
{
  double A;
  double *p = (double *)params;

  A = gsl_vector_get(v, 0);

  return pow(A - p[0], p[1]);
}
double TestGSL() {
  return gsl_ran_gaussian_pdf(0, 1);
}

double fast_log_sum_exp(double log_v1, double log_v2){
  if (log_v1 > log_v2){
    double diff = log_v2-log_v1;
    return diff < LOG_THRESH ? log_v1 : log_v1 + fastlog(1 + fastexp(diff));
  }
  else {
    double diff = log_v1-log_v2;
    return diff < LOG_THRESH ? log_v2 : log_v2 + fastlog(1 + fastexp(diff));
  }
}

double log_sum_exp(const std::vector<double>& log_vals) {
    if (log_vals.empty()) {
        return -1.0/0.0; // Return negative infinity
    }
    double max_val = log_vals[0];
    for (size_t i = 1; i < log_vals.size(); ++i) {
        if (log_vals[i] > max_val) {
            max_val = log_vals[i];
        }
    }
    double sum = 0.0;
    for (double val : log_vals) {
        sum += exp(val - max_val);
    }
    return max_val + log(sum);
}

double log_sum_exp(const double* start, const double* end) {
    if (start == end) {
        return -1.0/0.0; // Return negative infinity
    }
    double max_val = *start;
    for (const double* p = start + 1; p != end; ++p) {
        if (*p > max_val) {
            max_val = *p;
        }
    }
    double sum = 0.0;
    for (const double* p = start; p != end; ++p) {
        sum += exp(*p - max_val);
    }
    return max_val + log(sum);
}


