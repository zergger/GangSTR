#include <cmath>
#include <iostream>
#include <vector>

#include "src/hipstr_models/hip_em_learner.h"
#include "src/hipstr_models/hip_stutter_model.h"

namespace {

bool IsStrictProbability(double value) {
  return std::isfinite(value) && value > 0.0 && value < 1.0;
}

}  // namespace

int main() {
  std::vector< std::vector<int> > observed_alleles_by_sample;
  observed_alleles_by_sample.push_back(std::vector<int>{10, 10, 10, 10, 11, 10});
  observed_alleles_by_sample.push_back(std::vector<int>{12, 12, 12, 11, 12, 12});
  observed_alleles_by_sample.push_back(std::vector<int>{10, 10, 12, 12, 11, 11});

  HipEMLearner learner(observed_alleles_by_sample, 2);
  if (!learner.train(20, 0.01)) {
    std::cerr << "HipEMLearner training failed on a stable synthetic sample-aware dataset" << std::endl;
    return 1;
  }

  HipStutterModel* model = learner.get_stutter_model();
  if (model == NULL) {
    std::cerr << "HipEMLearner did not return a learned stutter model" << std::endl;
    return 1;
  }

  const double in_geom = model->get_parameter(true, 'P');
  const double in_up = model->get_parameter(true, 'U');
  const double in_down = model->get_parameter(true, 'D');
  const double out_geom = model->get_parameter(false, 'P');
  const double out_up = model->get_parameter(false, 'U');
  const double out_down = model->get_parameter(false, 'D');

  if (!IsStrictProbability(in_geom) ||
      !IsStrictProbability(in_up) ||
      !IsStrictProbability(in_down) ||
      !IsStrictProbability(out_geom) ||
      !IsStrictProbability(out_up) ||
      !IsStrictProbability(out_down)) {
    std::cerr << "HipEMLearner produced invalid stutter model parameters" << std::endl;
    return 1;
  }

  if (in_up + in_down + out_up + out_down >= 1.0) {
    std::cerr << "HipEMLearner produced incompatible transition probabilities" << std::endl;
    return 1;
  }

  if (model->period() != 2) {
    std::cerr << "HipEMLearner returned a stutter model with the wrong motif length" << std::endl;
    return 1;
  }

  const double ll_equal = model->log_stutter_pmf(20, 20);
  const double ll_up = model->log_stutter_pmf(20, 22);
  const double ll_down = model->log_stutter_pmf(20, 18);
  if (!(std::isfinite(ll_equal) && std::isfinite(ll_up) && std::isfinite(ll_down))) {
    std::cerr << "HipEMLearner produced non-finite stutter likelihoods" << std::endl;
    return 1;
  }
  if (!(ll_equal > ll_up && ll_equal > ll_down)) {
    std::cerr << "HipEMLearner learned a model that does not prefer the exact allele over one-step stutter" << std::endl;
    return 1;
  }

  return 0;
}
