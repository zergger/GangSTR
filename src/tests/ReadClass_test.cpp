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

#include "src/tests/ReadClass_test.h"

#include <gsl/gsl_cdf.h>
#include <gsl/gsl_randist.h>

#include <cmath>
#include <string>

using namespace std;

CPPUNIT_TEST_SUITE_REGISTRATION(ReadClassTest);

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

void AssertAlmostEqual(const double observed, const double expected, const double tol = 1e-9) {
  CPPUNIT_ASSERT_DOUBLES_EQUAL(expected, observed, tol);
}

}  // namespace

void ReadClassTest::setUp() {
  sample_profile_ = BuildSampleProfile(400.0, 50.0, 2000);
  str_info_.exp_thresh = 0;
  str_info_.stutter_up = 0.01;
  str_info_.stutter_down = 0.02;
  str_info_.stutter_p = 0.95;

  encl_class_.SetGlobalParams(sample_profile_, 2000, false, false);
  span_class_.SetGlobalParams(sample_profile_, 2000, false, false);
  frr_class_.SetGlobalParams(sample_profile_, 2000, false, false);

  encl_class_.SetLocusParams(str_info_);
  span_class_.SetLocusParams(str_info_);
  frr_class_.SetLocusParams(str_info_);

  encl_class_.SetCoverage(30);
  span_class_.SetCoverage(30);
  frr_class_.SetCoverage(30);

  encl_class_.Reset();
  span_class_.Reset();
  frr_class_.Reset();

  read_len = 100;
  motif_len = 3;
  ref_count = 10;
  ploidy = 2;
}

void ReadClassTest::tearDown() {}

void ReadClassTest::test_AddData() {
  encl_class_.AddData(10);
  encl_class_.AddData(20);
  span_class_.AddData(10);
  span_class_.AddData(20);
  span_class_.AddData(30);
  frr_class_.AddData(10);
  CPPUNIT_ASSERT_EQUAL(static_cast<int>(encl_class_.GetDataSize()), 2);
  CPPUNIT_ASSERT_EQUAL(static_cast<int>(span_class_.GetDataSize()), 3);
  CPPUNIT_ASSERT_EQUAL(static_cast<int>(frr_class_.GetDataSize()), 1);
}

void ReadClassTest::test_Reset() {
  encl_class_.AddData(10);
  encl_class_.AddData(20);
  span_class_.AddData(10);
  span_class_.AddData(20);
  span_class_.AddData(30);
  frr_class_.AddData(10);
  encl_class_.Reset();
  span_class_.Reset();
  frr_class_.Reset();
  CPPUNIT_ASSERT_EQUAL(static_cast<int>(encl_class_.GetDataSize()), 0);
  CPPUNIT_ASSERT_EQUAL(static_cast<int>(span_class_.GetDataSize()), 0);
  CPPUNIT_ASSERT_EQUAL(static_cast<int>(frr_class_.GetDataSize()), 0);
}

void ReadClassTest::test_GetReadDictStr() {
  encl_class_.AddData(20);
  encl_class_.AddData(10);
  encl_class_.AddData(20);
  encl_class_.AddData(30);
  encl_class_.AddData(20);
  encl_class_.AddData(10);
  CPPUNIT_ASSERT_EQUAL(string("10,2|20,3|30,1"), encl_class_.GetReadDictStr());
}

void ReadClassTest::test_SpanClassProb() {
  double log_class_prob = 0.0;
  CPPUNIT_ASSERT(span_class_.GetLogClassProb(25, read_len, motif_len, &log_class_prob));
  AssertAlmostEqual(log_class_prob, log(0.0838726946383), 1e-12);
}

void ReadClassTest::test_SpanReadProb() {
  double log_allele_prob = 0.0;
  CPPUNIT_ASSERT(span_class_.GetLogReadProb(25, 450, read_len, motif_len, ref_count, &log_allele_prob));
  AssertAlmostEqual(log_allele_prob, log(0.00131231629549), 1e-11);
}

void ReadClassTest::test_FRRClassProb() {
  double log_class_prob = 0.0;
  CPPUNIT_ASSERT(frr_class_.GetLogClassProb(45, read_len, motif_len, &log_class_prob));
  AssertAlmostEqual(log_class_prob, log(0.0177896348168 / 2.0), 1e-11);
}

void ReadClassTest::test_FRRReadProb() {
  double log_allele_prob = 0.0;
  CPPUNIT_ASSERT(frr_class_.GetLogReadProb(45, 80, read_len, motif_len, ref_count, &log_allele_prob));
  AssertAlmostEqual(log_allele_prob, log(0.0363690786878), 1e-12);
}

void ReadClassTest::test_EnclosingClassProb() {
  double log_class_prob = 0.0;
  CPPUNIT_ASSERT(encl_class_.GetLogClassProb(25, read_len, motif_len, &log_class_prob));
  AssertAlmostEqual(log_class_prob, log(0.0129032258065 / 2.0), 1e-11);
}

void ReadClassTest::test_EnclosingReadProb() {
  double log_allele_prob = 0.0;
  CPPUNIT_ASSERT(encl_class_.GetLogReadProb(25, 25, motif_len, nullptr, &log_allele_prob));
  AssertAlmostEqual(log_allele_prob, log(0.97), 1e-12);
}

void ReadClassTest::test_GetClassLogLikelihood() {
  double class_ll = 0.0;

  encl_class_.AddData(10);
  encl_class_.AddData(20);
  encl_class_.AddData(30);
  encl_class_.AddData(40);
  CPPUNIT_ASSERT(encl_class_.GetClassLogLikelihood(20, 50, read_len, motif_len, ref_count, ploidy, nullptr, &class_ll));
  // Regression anchor for the current enclosing-likelihood implementation.
  AssertAlmostEqual(class_ll, -134.979508559486, 1e-9);

  span_class_.Reset();
  span_class_.AddData(20);
  span_class_.AddData(30);
  span_class_.AddData(40);
  CPPUNIT_ASSERT(span_class_.GetClassLogLikelihood(20, 50, read_len, motif_len, ref_count, ploidy, &class_ll));
  AssertAlmostEqual(class_ll, -62.392179166719, 1e-8);

  frr_class_.Reset();
  frr_class_.AddData(10);
  frr_class_.AddData(20);
  frr_class_.AddData(30);
  frr_class_.AddData(40);
  CPPUNIT_ASSERT(frr_class_.GetClassLogLikelihood(20, 50, read_len, motif_len, ref_count, ploidy, &class_ll));
  AssertAlmostEqual(class_ll, -40.8239143859, 1e-9);
}

void ReadClassTest::test_GetAlleleLogLikelihood() {
  double allele_ll = 0.0;
  CPPUNIT_ASSERT(encl_class_.GetAlleleLogLikelihood(25, 24, read_len, motif_len, ref_count, &allele_ll));
  AssertAlmostEqual(allele_ll, -9.00674141673, 1e-9);

  CPPUNIT_ASSERT(span_class_.GetAlleleLogLikelihood(55, 430, read_len, motif_len, ref_count, &allele_ll));
  AssertAlmostEqual(allele_ll, -13.1016086835, 1e-9);

  CPPUNIT_ASSERT(frr_class_.GetAlleleLogLikelihood(55, 80, read_len, motif_len, ref_count, &allele_ll));
  AssertAlmostEqual(allele_ll, -6.17069762893, 1e-9);
}
