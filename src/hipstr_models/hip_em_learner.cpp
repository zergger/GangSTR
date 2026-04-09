#include <algorithm>
#include <cfloat>
#include <cstring>
#include <sstream>

#include "hip_em_learner.h"
#include "../mathops.h"

const double TOLERANCE = 1e-6;
const double LOG_ONE_HALF = log(0.5);
const double INFRAME_UP_PSEUDOCOUNT = 0.5;
const double INFRAME_DOWN_PSEUDOCOUNT = 2.0;
const double OUTFRAME_UP_PSEUDOCOUNT = 1.0;
const double OUTFRAME_DOWN_PSEUDOCOUNT = 1.0;
const double EQ_PSEUDOCOUNT = 1.0;
const double DIFF_PSEUDOCOUNT = 1.0;
const double DIFF_TAIL_PSEUDOCOUNT = 1.1;

double HipEMLearner::calc_log_likelihood() {
    for (int read_index = 0; read_index < num_reads_; ++read_index) {
        for (int allele_index = 0; allele_index < num_alleles_; ++allele_index) {
            log_allele_likelihoods_[read_index * num_alleles_ + allele_index] =
                stutter_model_->log_stutter_pmf(alleles_[allele_index] * motif_len_,
                                                alleles_[allele_index_[read_index]] * motif_len_);
        }
    }

    double total_log_likelihood = 0.0;
    const int num_diplotypes = num_alleles_ * num_alleles_;
    std::fill(log_sample_posteriors_, log_sample_posteriors_ + num_samples_ * num_diplotypes, 0.0);

    for (int sample_index = 0; sample_index < num_samples_; ++sample_index) {
        double* sample_posteriors = log_sample_posteriors_ + sample_index * num_diplotypes;
        for (int index_1 = 0; index_1 < num_alleles_; ++index_1) {
            for (int index_2 = 0; index_2 < num_alleles_; ++index_2) {
                sample_posteriors[index_1 * num_alleles_ + index_2] =
                    log_gt_priors_[index_1] + log_gt_priors_[index_2];
            }
        }
    }

    for (int read_index = 0; read_index < num_reads_; ++read_index) {
        double* sample_posteriors = log_sample_posteriors_ + sample_label_[read_index] * num_diplotypes;
        double* read_likelihoods = log_allele_likelihoods_ + read_index * num_alleles_;
        for (int index_1 = 0; index_1 < num_alleles_; ++index_1) {
            for (int index_2 = 0; index_2 < num_alleles_; ++index_2) {
                sample_posteriors[index_1 * num_alleles_ + index_2] +=
                    fast_log_sum_exp(LOG_ONE_HALF + read_likelihoods[index_1],
                                     LOG_ONE_HALF + read_likelihoods[index_2]);
            }
        }
    }

    for (int sample_index = 0; sample_index < num_samples_; ++sample_index) {
        double* sample_posteriors = log_sample_posteriors_ + sample_index * num_diplotypes;
        double norm = log_sum_exp(sample_posteriors, sample_posteriors + num_diplotypes);
        total_log_likelihood += norm;
        for (int i = 0; i < num_diplotypes; ++i) {
            sample_posteriors[i] -= norm;
        }
    }

    return total_log_likelihood;
}

void HipEMLearner::recalc_stutter_model() {
    std::vector<double> in_log_up, in_log_down, in_log_eq, in_log_diffs;
    std::vector<double> out_log_up, out_log_down, out_log_diffs;

    // Use a weak directional prior that favors downward stutter over upward stutter.
    // This reflects typical PCR slippage behavior while still allowing the data to dominate
    // once enough enclosing reads are available.
    in_log_up.push_back(log(INFRAME_UP_PSEUDOCOUNT));
    in_log_down.push_back(log(INFRAME_DOWN_PSEUDOCOUNT));
    in_log_diffs.push_back(log(DIFF_PSEUDOCOUNT));
    in_log_diffs.push_back(log(DIFF_TAIL_PSEUDOCOUNT));

    out_log_up.push_back(log(OUTFRAME_UP_PSEUDOCOUNT));
    out_log_down.push_back(log(OUTFRAME_DOWN_PSEUDOCOUNT));
    out_log_diffs.push_back(log(DIFF_PSEUDOCOUNT));
    out_log_diffs.push_back(log(DIFF_TAIL_PSEUDOCOUNT));

    in_log_eq.push_back(log(EQ_PSEUDOCOUNT));

    const int num_diplotypes = num_alleles_ * num_alleles_;
    for (int read_index = 0; read_index < num_reads_; ++read_index) {
        const int observed_allele = allele_index_[read_index];
        double* sample_posteriors = log_sample_posteriors_ + sample_label_[read_index] * num_diplotypes;
        double* read_likelihoods = log_allele_likelihoods_ + read_index * num_alleles_;
        for (int index_1 = 0; index_1 < num_alleles_; ++index_1) {
            for (int index_2 = 0; index_2 < num_alleles_; ++index_2) {
                const double posterior = exp(sample_posteriors[index_1 * num_alleles_ + index_2]);
                const double p_obs_given_1 = exp(read_likelihoods[index_1]);
                const double p_obs_given_2 = exp(read_likelihoods[index_2]);
                const double denom = p_obs_given_1 + p_obs_given_2;

                double expected_from_1 = 0.0;
                double expected_from_2 = 0.0;
                if (index_1 == index_2) {
                    expected_from_1 = posterior;
                } else if (denom > 0.0) {
                    expected_from_1 = posterior * p_obs_given_1 / denom;
                    expected_from_2 = posterior * p_obs_given_2 / denom;
                }

                int bp_diff_1 = (alleles_[observed_allele] - alleles_[index_1]) * motif_len_;
                int bp_diff_2 = (alleles_[observed_allele] - alleles_[index_2]) * motif_len_;

                if (expected_from_1 > 0.0) {
                    if (bp_diff_1 == 0) {
                        in_log_eq.push_back(log(expected_from_1));
                    } else if (bp_diff_1 % motif_len_ != 0) {
                        out_log_diffs.push_back(log(expected_from_1 * abs(bp_diff_1)));
                        if (bp_diff_1 > 0) out_log_up.push_back(log(expected_from_1));
                        else out_log_down.push_back(log(expected_from_1));
                    } else {
                        in_log_diffs.push_back(log(expected_from_1 * abs(bp_diff_1 / motif_len_)));
                        if (bp_diff_1 > 0) in_log_up.push_back(log(expected_from_1));
                        else in_log_down.push_back(log(expected_from_1));
                    }
                }

                if (expected_from_2 > 0.0) {
                    if (bp_diff_2 == 0) {
                        in_log_eq.push_back(log(expected_from_2));
                    } else if (bp_diff_2 % motif_len_ != 0) {
                        out_log_diffs.push_back(log(expected_from_2 * abs(bp_diff_2)));
                        if (bp_diff_2 > 0) out_log_up.push_back(log(expected_from_2));
                        else out_log_down.push_back(log(expected_from_2));
                    } else {
                        in_log_diffs.push_back(log(expected_from_2 * abs(bp_diff_2 / motif_len_)));
                        if (bp_diff_2 > 0) in_log_up.push_back(log(expected_from_2));
                        else in_log_down.push_back(log(expected_from_2));
                    }
                }
            }
        }
    }

    double in_log_total_up = log_sum_exp(in_log_up);
    double in_log_total_down = log_sum_exp(in_log_down);
    double in_log_total_eq = log_sum_exp(in_log_eq);
    double in_log_total_diffs = log_sum_exp(in_log_diffs);
    double out_log_total_up = log_sum_exp(out_log_up);
    double out_log_total_down = log_sum_exp(out_log_down);
    double out_log_total_diffs = log_sum_exp(out_log_diffs);
    double out_log_total = fast_log_sum_exp(out_log_total_up, out_log_total_down);
    double in_pgeom_hat = std::min(0.999, exp(fast_log_sum_exp(in_log_total_up, in_log_total_down) - in_log_total_diffs));
    double out_pgeom_hat = std::min(0.999, exp(out_log_total - out_log_total_diffs));
    double log_total = fast_log_sum_exp(fast_log_sum_exp(in_log_total_up, in_log_total_down),
                                        fast_log_sum_exp(in_log_total_eq, out_log_total));
    double in_pup_hat = exp(in_log_total_up - log_total);
    double in_pdown_hat = exp(in_log_total_down - log_total);
    double out_pup_hat = exp(out_log_total_up - log_total);
    double out_pdown_hat = exp(out_log_total_down - log_total);

    delete stutter_model_;
    stutter_model_ = new HipStutterModel(in_pgeom_hat, in_pup_hat, in_pdown_hat,
                                         out_pgeom_hat, out_pup_hat, out_pdown_hat, motif_len_);
}

void HipEMLearner::recalc_log_gt_priors() {
    std::vector<double> allele_i_posteriors;
    for (int i = 0; i < num_alleles_; ++i) {
        allele_i_posteriors.clear();
        for (int sample_index = 0; sample_index < num_samples_; ++sample_index) {
            double* sample_posteriors = log_sample_posteriors_ + sample_index * num_alleles_ * num_alleles_;
            for (int j = 0; j < num_alleles_; ++j) {
                allele_i_posteriors.push_back(sample_posteriors[i * num_alleles_ + j]);
                allele_i_posteriors.push_back(sample_posteriors[j * num_alleles_ + i]);
            }
        }
        log_gt_priors_[i] = log_sum_exp(allele_i_posteriors);
    }
    double log_total = log_sum_exp(log_gt_priors_, log_gt_priors_ + num_alleles_);
    for (int i = 0; i < num_alleles_; ++i) {
        log_gt_priors_[i] -= log_total;
    }
}

void HipEMLearner::init_log_gt_priors() {
    std::vector<double> prior_mass(num_alleles_, 1.0);
    for (int read_index = 0; read_index < num_reads_; ++read_index) {
        const int sample = sample_label_[read_index];
        if (reads_per_sample_[sample] > 0) {
            prior_mass[allele_index_[read_index]] += 1.0 / reads_per_sample_[sample];
        }
    }
    double total = 0.0;
    for (int i = 0; i < num_alleles_; ++i) {
        total += prior_mass[i];
    }
    for (int i = 0; i < num_alleles_; ++i) {
        log_gt_priors_[i] = log(prior_mass[i]) - log(total);
    }
}

void HipEMLearner::init_stutter_model() {
    delete stutter_model_;
    stutter_model_ = new HipStutterModel(0.9, 0.05, 0.05, 0.8, 0.01, 0.01, motif_len_);
}

bool HipEMLearner::train(int max_iter, double min_LL_abs_change) {
    init_log_gt_priors();
    init_stutter_model();

    double LL = -DBL_MAX;
    for (int i = 0; i < max_iter; ++i) {
        double new_LL = calc_log_likelihood();

        if (i > 0 && (new_LL - LL) < min_LL_abs_change) {
            return true;
        }
        if (new_LL < LL - TOLERANCE) {
            return false;
        }
        LL = new_LL;

        recalc_log_gt_priors();
        recalc_stutter_model();
    }
    return true;
}
