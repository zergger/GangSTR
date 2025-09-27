#ifndef HIP_EM_LEARNER_H_
#define HIP_EM_LEARNER_H_

#include <algorithm>
#include <iostream>
#include <map>
#include <math.h>
#include <set>
#include <string>
#include <vector>

#include "../common.h"
#include "hip_stutter_model.h"

class HipEMLearner {
 private:
  int motif_len_;
  HipStutterModel* stutter_model_;
  std::vector<int> alleles_;
  std::vector<int> allele_counts_;
  int num_alleles_;
  int num_reads_;

  double* log_gt_priors_;
  double* log_allele_likelihoods_;
  double* log_sample_posteriors_; // P(g | D) for a single sample

  void init_allele_likelihoods();
  void init_log_gt_priors();
  void init_stutter_model();
  
  void recalc_log_gt_priors();
  void recalc_stutter_model();
  
  double calc_log_likelihood();

  // Private unimplemented copy constructor and assignment operator to prevent operations
  HipEMLearner(const HipEMLearner& other);
  HipEMLearner& operator=(const HipEMLearner& other);

 public:
  HipEMLearner(const std::vector<int>& observed_alleles, int motif_len) {
    motif_len_ = motif_len;
    num_reads_ = observed_alleles.size();

    std::map<int, int> allele_counts_map;
    for (int allele : observed_alleles) {
        allele_counts_map[allele]++;
    }
    for (auto const& [allele, count] : allele_counts_map) {
        alleles_.push_back(allele);
        allele_counts_.push_back(count);
    }
    num_alleles_ = alleles_.size();

    log_gt_priors_ = new double[num_alleles_];
    log_allele_likelihoods_ = new double[num_alleles_ * num_alleles_];
    log_sample_posteriors_ = new double[num_alleles_ * num_alleles_];
    stutter_model_ = NULL;
  }

  ~HipEMLearner(){
    delete [] log_gt_priors_;
    delete [] log_allele_likelihoods_;
    delete [] log_sample_posteriors_;
    delete stutter_model_;
  }  
  
  bool train(int max_iter, double min_LL_abs_change);

  HipStutterModel* get_stutter_model() const {
    if (stutter_model_ == NULL) {
        PrintMessageDieOnError("No stutter model has been specified or learned", M_ERROR, false);
    }
    return stutter_model_;
  }
};

#endif
