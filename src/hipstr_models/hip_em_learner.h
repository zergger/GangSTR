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
  int num_alleles_;
  int num_reads_;
  int num_samples_;

  int* allele_index_;
  int* sample_label_;
  std::vector<int> reads_per_sample_;

  double* log_gt_priors_;
  double* log_allele_likelihoods_;
  double* log_sample_posteriors_;

  void init_log_gt_priors();
  void init_stutter_model();
  
  void recalc_log_gt_priors();
  void recalc_stutter_model();
  
  double calc_log_likelihood();

  // Private unimplemented copy constructor and assignment operator to prevent operations
  HipEMLearner(const HipEMLearner& other);
  HipEMLearner& operator=(const HipEMLearner& other);

 public:
  HipEMLearner(const std::vector< std::vector<int> >& observed_alleles_by_sample, int motif_len) {
    motif_len_ = motif_len;
    num_samples_ = observed_alleles_by_sample.size();
    num_reads_ = 0;

    std::map<int, int> allele_counts_map;
    for (const std::vector<int>& sample_alleles : observed_alleles_by_sample) {
      num_reads_ += sample_alleles.size();
      for (int allele : sample_alleles) {
        allele_counts_map[allele]++;
      }
    }
    for (auto const& [allele, count] : allele_counts_map) {
      (void)count;
      alleles_.push_back(allele);
    }
    num_alleles_ = alleles_.size();

    std::map<int, int> allele_to_index;
    for (int i = 0; i < num_alleles_; ++i) {
      allele_to_index[alleles_[i]] = i;
    }

    allele_index_ = new int[num_reads_];
    sample_label_ = new int[num_reads_];
    log_gt_priors_ = new double[num_alleles_];
    log_allele_likelihoods_ = new double[num_reads_ * num_alleles_];
    log_sample_posteriors_ = new double[num_samples_ * num_alleles_ * num_alleles_];
    stutter_model_ = NULL;

    int read_index = 0;
    for (int sample_index = 0; sample_index < num_samples_; ++sample_index) {
      reads_per_sample_.push_back(observed_alleles_by_sample[sample_index].size());
      for (int allele : observed_alleles_by_sample[sample_index]) {
        allele_index_[read_index] = allele_to_index[allele];
        sample_label_[read_index] = sample_index;
        ++read_index;
      }
    }
  }

  ~HipEMLearner(){
    delete [] allele_index_;
    delete [] sample_label_;
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
