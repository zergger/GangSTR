/*
Copyright (C) 2017 Melissa Gymrek <mgymrek@ucsd.edu>
and Nima Mousavi (mousavi@ucsd.edu)

This file is part of GangSTR.

GangSTR is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

GangSTR is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with GangSTR.  If not, see <http://www.gnu.org/licenses/>.
*/

#include "src/tests/LikelihoodMaximizer_test.h"

#include <gsl/gsl_cdf.h>
#include <gsl/gsl_randist.h>

#include <cmath>

using namespace std;

CPPUNIT_TEST_SUITE_REGISTRATION(LikelihoodMaximizerTest);

namespace {

SampleProfile BuildSampleProfile(const double mean, const double sdev, const int dist_size) {
  SampleProfile sp;
  sp.rg_sample = "test";
  sp.rg_id = "test_rg";
  sp.dist_mean = mean;
  sp.dist_sdev = sdev;
  sp.coverage = 30.0;
  sp.dist_pdf.resize(dist_size);
  sp.dist_cdf.resize(dist_size);
  sp.dist_integral.resize(dist_size);

  double running = 0.0;
  for (int i = 0; i < dist_size; ++i) {
    sp.dist_pdf[i] = gsl_ran_gaussian_pdf(i - mean, sdev);
    sp.dist_cdf[i] = gsl_cdf_gaussian_P(i - mean, sdev);
    running += i * sp.dist_pdf[i];
    sp.dist_integral[i] = running;
  }
  if (!sp.dist_cdf.empty()) {
    sp.dist_cdf.back() = 1.0;
  }
  return sp;
}

}  // namespace

void LikelihoodMaximizerTest::setUp() {
  options = Options();
  options.stutter_up = 0.01;
  options.stutter_down = 0.02;
  options.stutter_p = 0.95;
  options.flanklen = 2000;
  options.realignment_flanklen = 100;
  options.frr_weight = 1.0;
  options.enclosing_weight = 1.0;
  options.spanning_weight = 1.0;
  options.flanking_weight = 1.0;
  options.verbose = false;
  options.min_match = 0;
  options.read_len = 100;
  options.use_cov = false;
  options.hist_mode = false;

  read_len = 100;
  motif_len = 3;
  ref_count = 10;
  resampled = false;

  locus.chrom = "19";
  locus.start = 5000;
  locus.end = 5039;
  locus.motif = "CTG";
  locus.period = 3;

  sample_profile_ = BuildSampleProfile(400.0, 50.0, 2000);
  str_info_.exp_thresh = 0;
  str_info_.stutter_up = 0.01;
  str_info_.stutter_down = 0.02;
  str_info_.stutter_p = 0.95;

  likelihood_maximizer_ = new LikelihoodMaximizer(options, sample_profile_, read_len, "F");
  likelihood_maximizer_->Reset();
  likelihood_maximizer_->SetLocusParams(str_info_, sample_profile_.coverage, read_len, motif_len, ref_count, locus.chrom);
}

void LikelihoodMaximizerTest::tearDown() {
  delete likelihood_maximizer_;
  likelihood_maximizer_ = nullptr;
}

void LikelihoodMaximizerTest::test_Reset() {
  likelihood_maximizer_->AddEnclosingData(10);
  likelihood_maximizer_->AddSpanningData(20);
  likelihood_maximizer_->AddFRRData(30);
  likelihood_maximizer_->Reset();
  CPPUNIT_ASSERT_EQUAL(static_cast<int>(likelihood_maximizer_->GetEnclosingDataSize()), 0);
  CPPUNIT_ASSERT_EQUAL(static_cast<int>(likelihood_maximizer_->GetSpanningDataSize()), 0);
  CPPUNIT_ASSERT_EQUAL(static_cast<int>(likelihood_maximizer_->GetFRRDataSize()), 0);
  CPPUNIT_ASSERT_EQUAL(static_cast<int>(likelihood_maximizer_->GetReadPoolSize()), 0);
}

void LikelihoodMaximizerTest::test_AddEnclosingData() {
  likelihood_maximizer_->Reset();
  likelihood_maximizer_->SetLocusParams(str_info_, sample_profile_.coverage, read_len, motif_len, ref_count, locus.chrom);
  likelihood_maximizer_->AddEnclosingData(10);
  likelihood_maximizer_->AddEnclosingData(20);
  likelihood_maximizer_->AddEnclosingData(30);
  CPPUNIT_ASSERT_EQUAL(static_cast<int>(likelihood_maximizer_->GetEnclosingDataSize()), 3);
  CPPUNIT_ASSERT_EQUAL(static_cast<int>(likelihood_maximizer_->GetSpanningDataSize()), 0);
  CPPUNIT_ASSERT_EQUAL(static_cast<int>(likelihood_maximizer_->GetFRRDataSize()), 0);
}

void LikelihoodMaximizerTest::test_AddSpanningData() {
  likelihood_maximizer_->Reset();
  likelihood_maximizer_->SetLocusParams(str_info_, sample_profile_.coverage, read_len, motif_len, ref_count, locus.chrom);
  likelihood_maximizer_->AddSpanningData(10);
  CPPUNIT_ASSERT_EQUAL(static_cast<int>(likelihood_maximizer_->GetSpanningDataSize()), 1);
  CPPUNIT_ASSERT_EQUAL(static_cast<int>(likelihood_maximizer_->GetEnclosingDataSize()), 0);
  CPPUNIT_ASSERT_EQUAL(static_cast<int>(likelihood_maximizer_->GetFRRDataSize()), 0);
}

void LikelihoodMaximizerTest::test_AddFRRData() {
  likelihood_maximizer_->Reset();
  likelihood_maximizer_->SetLocusParams(str_info_, sample_profile_.coverage, read_len, motif_len, ref_count, locus.chrom);
  likelihood_maximizer_->AddFRRData(10);
  likelihood_maximizer_->AddFRRData(20);
  CPPUNIT_ASSERT_EQUAL(static_cast<int>(likelihood_maximizer_->GetFRRDataSize()), 2);
  CPPUNIT_ASSERT_EQUAL(static_cast<int>(likelihood_maximizer_->GetEnclosingDataSize()), 0);
  CPPUNIT_ASSERT_EQUAL(static_cast<int>(likelihood_maximizer_->GetSpanningDataSize()), 0);
}

void LikelihoodMaximizerTest::test_GetGenotypeNegLogLikelihood() {
  likelihood_maximizer_->Reset();
  likelihood_maximizer_->SetLocusParams(str_info_, sample_profile_.coverage, read_len, motif_len, ref_count, locus.chrom);

  likelihood_maximizer_->AddEnclosingData(25);
  likelihood_maximizer_->AddEnclosingData(25);
  likelihood_maximizer_->AddSpanningData(445);

  double ll_match = 0.0;
  double ll_mismatch = 0.0;
  CPPUNIT_ASSERT(likelihood_maximizer_->GetGenotypeNegLogLikelihood(25, 25, resampled, &ll_match));
  CPPUNIT_ASSERT(likelihood_maximizer_->GetGenotypeNegLogLikelihood(18, 18, resampled, &ll_mismatch));
  CPPUNIT_ASSERT(std::isfinite(ll_match));
  CPPUNIT_ASSERT(std::isfinite(ll_mismatch));
  CPPUNIT_ASSERT(ll_match < ll_mismatch);
}

void LikelihoodMaximizerTest::test_OptimizeLikelihood() {
  likelihood_maximizer_->Reset();
  likelihood_maximizer_->SetLocusParams(str_info_, sample_profile_.coverage, read_len, motif_len, ref_count, locus.chrom);

  likelihood_maximizer_->AddEnclosingData(25);
  likelihood_maximizer_->AddEnclosingData(25);
  likelihood_maximizer_->AddEnclosingData(25);
  likelihood_maximizer_->SetGridSize(20, 30);

  int32_t allele1 = -1;
  int32_t allele2 = -1;
  double min_neg_like = 0.0;
  CPPUNIT_ASSERT(likelihood_maximizer_->OptimizeLikelihood(false, 2, 0, 0.0, &allele1, &allele2, &min_neg_like));
  CPPUNIT_ASSERT_EQUAL(25, allele1);
  CPPUNIT_ASSERT_EQUAL(25, allele2);
  CPPUNIT_ASSERT(std::isfinite(min_neg_like));
}
