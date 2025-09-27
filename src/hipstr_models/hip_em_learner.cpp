#include <algorithm>
#include <cfloat>
#include <cstring>
#include <sstream>

#include "hip_em_learner.h"
#include "../mathops.h"

const double TOLERANCE = 1e-6;

// E-Step: Calculate expected genotype counts given current stutter model
double HipEMLearner::calc_log_likelihood() {
    // P(A_i | g_j) - Likelihood of observed allele i given true genotype j
    for (int i = 0; i < num_alleles_; ++i) {
        for (int j = 0; j < num_alleles_; ++j) {
            log_allele_likelihoods_[i * num_alleles_ + j] = stutter_model_->log_stutter_pmf(alleles_[j] * motif_len_, alleles_[i] * motif_len_);
        }
    }

    double total_log_likelihood = 0;
    // P(g_j, g_k | D) ~ P(D | g_j, g_k) * P(g_j) * P(g_k)
    for (int j = 0; j < num_alleles_; ++j) {
        for (int k = j; k < num_alleles_; ++k) {
            double log_likelihood_diplotype = 0;
            for (int i = 0; i < num_alleles_; ++i) {
                double log_likelihood_allele_i = fast_log_sum_exp(
                    log_allele_likelihoods_[i * num_alleles_ + j],
                    log_allele_likelihoods_[i * num_alleles_ + k]
                );
                log_likelihood_diplotype += allele_counts_[i] * log_likelihood_allele_i;
            }
            double prior = (j == k) ? log_gt_priors_[j] * 2 : log_gt_priors_[j] + log_gt_priors_[k];
            log_sample_posteriors_[j * num_alleles_ + k] = log_likelihood_diplotype + prior;
            log_sample_posteriors_[k * num_alleles_ + j] = log_sample_posteriors_[j * num_alleles_ + k];
        }
    }

    // Normalize posteriors
    double norm = log_sum_exp(log_sample_posteriors_, log_sample_posteriors_ + num_alleles_ * num_alleles_);
    for (int i = 0; i < num_alleles_ * num_alleles_; ++i) {
        log_sample_posteriors_[i] -= norm;
    }
    return norm; // Return total log likelihood
}

// M-step: Recalculate stutter model parameters based on expected counts
void HipEMLearner::recalc_stutter_model() {
    std::vector<double> in_log_up, in_log_down, in_log_eq, in_log_diffs;
    std::vector<double> out_log_up, out_log_down, out_log_diffs;

    // Add pseudocounts
    in_log_up.push_back(0.0); in_log_down.push_back(0.0); in_log_diffs.push_back(0.0); in_log_diffs.push_back(log(1.1));
    out_log_up.push_back(0.0); out_log_down.push_back(0.0); out_log_diffs.push_back(0.0); out_log_diffs.push_back(log(1.1));
    in_log_eq.push_back(0.0);

    for (int i = 0; i < num_alleles_; ++i) { // Observed allele
        for (int j = 0; j < num_alleles_; ++j) { // True allele 1
            for (int k = j; k < num_alleles_; ++k) { // True allele 2
                double posterior = exp(log_sample_posteriors_[j * num_alleles_ + k]);
                double p_i_given_j = exp(log_allele_likelihoods_[i * num_alleles_ + j]);
                double p_i_given_k = exp(log_allele_likelihoods_[i * num_alleles_ + k]);
                
                double expected_from_j = posterior * p_i_given_j / (p_i_given_j + p_i_given_k);
                double expected_from_k = posterior * p_i_given_k / (p_i_given_j + p_i_given_k);

                if (j==k) expected_from_j = posterior;

                int bp_diff_j = (alleles_[i] - alleles_[j]) * motif_len_;
                int bp_diff_k = (alleles_[i] - alleles_[k]) * motif_len_;

                // Contribution from allele j
                if (bp_diff_j == 0) {
                    in_log_eq.push_back(log(expected_from_j * allele_counts_[i]));
                } else if (bp_diff_j % motif_len_ != 0) {
                    out_log_diffs.push_back(log(expected_from_j * allele_counts_[i] * abs(bp_diff_j)));
                    if (bp_diff_j > 0) out_log_up.push_back(log(expected_from_j * allele_counts_[i]));
                    else out_log_down.push_back(log(expected_from_j * allele_counts_[i]));
                } else {
                    in_log_diffs.push_back(log(expected_from_j * allele_counts_[i] * abs(bp_diff_j / motif_len_)));
                    if (bp_diff_j > 0) in_log_up.push_back(log(expected_from_j * allele_counts_[i]));
                    else in_log_down.push_back(log(expected_from_j * allele_counts_[i]));
                }

                if (j==k) continue;

                // Contribution from allele k
                 if (bp_diff_k == 0) {
                    in_log_eq.push_back(log(expected_from_k * allele_counts_[i]));
                } else if (bp_diff_k % motif_len_ != 0) {
                    out_log_diffs.push_back(log(expected_from_k * allele_counts_[i] * abs(bp_diff_k)));
                    if (bp_diff_k > 0) out_log_up.push_back(log(expected_from_k * allele_counts_[i]));
                    else out_log_down.push_back(log(expected_from_k * allele_counts_[i]));
                } else {
                    in_log_diffs.push_back(log(expected_from_k * allele_counts_[i] * abs(bp_diff_k / motif_len_)));
                    if (bp_diff_k > 0) in_log_up.push_back(log(expected_from_k * allele_counts_[i]));
                    else in_log_down.push_back(log(expected_from_k * allele_counts_[i]));
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
    double log_total = fast_log_sum_exp(fast_log_sum_exp(in_log_total_up, in_log_total_down), fast_log_sum_exp(in_log_total_eq, out_log_total));
    double in_pup_hat = exp(in_log_total_up - log_total);
    double in_pdown_hat = exp(in_log_total_down - log_total);
    double out_pup_hat = exp(out_log_total_up - log_total);
    double out_pdown_hat = exp(out_log_total_down - log_total);

    delete stutter_model_;
    stutter_model_ = new HipStutterModel(in_pgeom_hat, in_pup_hat, in_pdown_hat, out_pgeom_hat, out_pup_hat, out_pdown_hat, motif_len_);
}

// M-step: Recalculate priors based on expected genotype counts
void HipEMLearner::recalc_log_gt_priors() {
    for (int i = 0; i < num_alleles_; ++i) {
        std::vector<double> allele_i_posteriors;
        for (int j = 0; j < num_alleles_; ++j) {
            allele_i_posteriors.push_back(log_sample_posteriors_[i * num_alleles_ + j]);
            allele_i_posteriors.push_back(log_sample_posteriors_[j * num_alleles_ + i]);
        }
        log_gt_priors_[i] = log_sum_exp(allele_i_posteriors);
    }
    double log_total = log_sum_exp(log_gt_priors_, log_gt_priors_ + num_alleles_);
    for (int i = 0; i < num_alleles_; i++) {
        log_gt_priors_[i] -= log_total;
    }
}

void HipEMLearner::init_log_gt_priors() {
    for (int i = 0; i < num_alleles_; i++) {
        log_gt_priors_[i] = log(allele_counts_[i] + 1.0);
    }
    double log_total = log_sum_exp(log_gt_priors_, log_gt_priors_ + num_alleles_);
    for (int i = 0; i < num_alleles_; i++) {
        log_gt_priors_[i] -= log_total;
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
        // E-step
        double new_LL = calc_log_likelihood();

        // Check for convergence
        if (i > 0 && (new_LL - LL) < min_LL_abs_change) {
            return true;
        }
        if (new_LL < LL - TOLERANCE) {
            // This can happen due to pseudocounts, but a large drop is an issue
            return false;
        }
        LL = new_LL;

        // M-step
        recalc_log_gt_priors();
        recalc_stutter_model();
    }
    return true; // Reached max iterations
}
