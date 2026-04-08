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

#include "src/tests/ReadExtractor_test.h"

#include <gsl/gsl_cdf.h>
#include <gsl/gsl_randist.h>

#include <cmath>
#include <cstdlib>
#include <vector>

using namespace std;

CPPUNIT_TEST_SUITE_REGISTRATION(ReadExtractorTest);

namespace {

SampleProfile BuildGaussianProfile(const double mean, const double sdev, const int dist_size) {
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

string ReadExtractorTest::ResolveTestDir() const {
  const char* from_env = getenv("GANGSTR_TEST_DIR");
  if (from_env != nullptr && from_env[0] != '\0') {
    return string(from_env);
  }
#ifdef GANGSTR_TEST_DATA_DIR
  return string(GANGSTR_TEST_DATA_DIR);
#else
  return "unused/tests";
#endif
}

SampleProfile ReadExtractorTest::BuildSampleProfile() const {
  return BuildGaussianProfile(500.0, 50.0, 2000);
}

void ReadExtractorTest::RegisterFixtureProfile(const string& bam_path) {
  sample_info_.RegisterSampleProfile(bam_path, "test", BuildSampleProfile(), "F", 100, true);
}

void ReadExtractorTest::setUp() {
  test_dir = ResolveTestDir();

  locus.chrom = "3";
  locus.start = 63898362;
  locus.end = 63898391;
  locus.pre_flank = "taggagcggaaagaatgtcggagcgggccgcggatgacgtcaggggggagccgcgccgcgcggcggcggcggcgggcggagcagcggccgcggccgcccgg";
  locus.post_flank = "ccgccgcctccgcagccccagcggcagcagcacccgccaccgccgccacggcgcacacggccggaggacggcgggcccggcgccgcctccacctcggccgc";
  locus.motif = "cag";
  locus.period = 3;

  regionsize = 5000;
  min_match = 0;

  options = Options();
  options.flanklen = 3000;
  options.read_len = 100;
  options.realignment_flanklen = 100;
  options.min_match = min_match;
  options.min_score = 0;
  options.max_processed_reads_per_sample = 10000;
  options.realign_match_perc = 0.9;
  options.output_readinfo = false;

  RegisterFixtureProfile(test_dir + "/test.sorted.bam");
  RegisterFixtureProfile(test_dir + "/test.enclosing.single.bam");
  RegisterFixtureProfile(test_dir + "/test.spanning.single.bam");
  RegisterFixtureProfile(test_dir + "/test.frr.bam");

  read_extractor_ = new ReadExtractor(options, sample_info_);
}

void ReadExtractorTest::tearDown() {
  delete read_extractor_;
  read_extractor_ = nullptr;
}

void ReadExtractorTest::test_ExtractReads() {
  const string bam_file = test_dir + "/test.sorted.bam";
  const string fastafile = test_dir + "/test.fa";
  vector<string> files(1, bam_file);
  BamCramMultiReader bamreader(files, fastafile);

  SampleProfile sp = BuildSampleProfile();
  STRLocusInfo sli;
  sli.exp_thresh = 0;
  sli.stutter_up = options.stutter_up;
  sli.stutter_down = options.stutter_down;
  sli.stutter_p = options.stutter_p;

  LikelihoodMaximizer likmax(options, sp, options.read_len, "F");
  likmax.SetLocusParams(sli, sp.coverage, options.read_len, locus.period, (locus.end - locus.start + 1) / locus.period, locus.chrom);

  map<string, LikelihoodMaximizer*> lms;
  lms["test"] = &likmax;

  CPPUNIT_ASSERT(read_extractor_->ExtractReads(&bamreader, locus, regionsize, min_match, lms));
  CPPUNIT_ASSERT(likmax.GetEnclosingDataSize() > 0);
  CPPUNIT_ASSERT(likmax.GetSpanningDataSize() > 0);
  CPPUNIT_ASSERT(likmax.GetFlankingDataSize() > 0);
}

void ReadExtractorTest::test_ProcessSingleRead() {
  {
    const string bam_file = test_dir + "/test.enclosing.single.bam";
    const string fastafile = test_dir + "/test.fa";
    vector<string> files(1, bam_file);
    BamCramMultiReader bamreader(files, fastafile);
    const int32_t buffer = 50000;
    const int32_t chrom_ref_id = bamreader.bam_header()->ref_id(locus.chrom);
    bamreader.SetRegion(locus.chrom, locus.start - buffer, locus.end + buffer);

    BamAlignment aln;
    int seen = 0;
    while (bamreader.GetNextAlignment(aln)) {
      ReadType read_type = RC_UNKNOWN;
      SingleReadType srt = SR_UNKNOWN;
      int32_t data_value = 0;
      int32_t ncopy_value = 0;
      int32_t score_value = 0;
      if (!read_extractor_->ProcessSingleRead(aln, chrom_ref_id, locus, min_match, false,
                                              &data_value, &ncopy_value, &score_value, &read_type, &srt)) {
        continue;
      }
      int64_t true_value = -1;
      CPPUNIT_ASSERT(aln.GetIntTag("nc", true_value));
      CPPUNIT_ASSERT_EQUAL(RC_ENCL, read_type);
      CPPUNIT_ASSERT_EQUAL(static_cast<int32_t>(true_value), data_value);
      ++seen;
    }
    CPPUNIT_ASSERT(seen > 0);
  }

  {
    const string bam_file = test_dir + "/test.spanning.single.bam";
    const string fastafile = test_dir + "/test.fa";
    vector<string> files(1, bam_file);
    BamCramMultiReader bamreader(files, fastafile);
    const int32_t buffer = 50000;
    const int32_t chrom_ref_id = bamreader.bam_header()->ref_id(locus.chrom);
    bamreader.SetRegion(locus.chrom, locus.start - buffer, locus.end + buffer);

    BamAlignment aln;
    bool found_spanning = false;
    while (bamreader.GetNextAlignment(aln)) {
      ReadType read_type = RC_UNKNOWN;
      SingleReadType srt = SR_UNKNOWN;
      int32_t data_value = 0;
      int32_t ncopy_value = 0;
      int32_t score_value = 0;
      if (!read_extractor_->ProcessSingleRead(aln, chrom_ref_id, locus, min_match, false,
                                              &data_value, &ncopy_value, &score_value, &read_type, &srt)) {
        continue;
      }
      if (read_type == RC_SPAN) {
        CPPUNIT_ASSERT(data_value > 0);
        found_spanning = true;
      }
    }
    CPPUNIT_ASSERT(found_spanning);
  }
}

void ReadExtractorTest::test_RescueMate() {
  const string bam_file = test_dir + "/test.frr.bam";
  const string fastafile = test_dir + "/test.fa";
  vector<string> files(1, bam_file);
  BamCramMultiReader bamreader(files, fastafile);
  const int32_t buffer = 50000;
  bamreader.SetRegion(locus.chrom, locus.start - buffer, locus.end + buffer);

  BamAlignment aln;
  bool checked = false;
  while (bamreader.GetNextAlignment(aln)) {
    if (aln.Name() != "ATXN7_52_cov60_dist500_DIP_const70_70_constAllele_3065_3556_0:0:0_0:0:0_136" &&
        aln.Name() != "ATXN7_52_cov60_dist500_DIP_const70_70_constAllele_2667_3154_0:0:0_0:0:0_a8") {
      continue;
    }
    BamAlignment matealn;
    CPPUNIT_ASSERT(read_extractor_->RescueMate(&bamreader, aln, &matealn));
    CPPUNIT_ASSERT_EQUAL(aln.MateRefID(), matealn.RefID());
    if (aln.Name() == "ATXN7_52_cov60_dist500_DIP_const70_70_constAllele_3065_3556_0:0:0_0:0:0_136") {
      CPPUNIT_ASSERT(std::abs(matealn.Position() - 53253385) <= 1);
    } else {
      CPPUNIT_ASSERT(std::abs(matealn.Position() - 53253384) <= 1);
    }
    checked = true;
    break;
  }
  CPPUNIT_ASSERT(checked);
}
