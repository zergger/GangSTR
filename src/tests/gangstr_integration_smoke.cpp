#include <cmath>
#include <iostream>

#include "src/flanking_class.h"
#include "src/bam_info_extract.h"
#include "src/hipstr_models/hip_stutter_model.h"
#include "src/str_info.h"

namespace {

SampleProfile MakeSampleProfile() {
  SampleProfile sp;
  sp.dist_mean = 400;
  sp.dist_sdev = 50;
  sp.coverage = 30;
  sp.dist_pdf.assign(1000, 0.0);
  sp.dist_cdf.assign(1000, 1.0);
  return sp;
}

STRLocusInfo MakeLocusInfo() {
  STRLocusInfo info;
  info.exp_thresh = 0;
  info.stutter_up = 0.05;
  info.stutter_down = 0.05;
  info.stutter_p = 0.9;
  return info;
}

bool NearlyEqual(double lhs, double rhs, double tol) {
  return std::fabs(lhs - rhs) <= tol;
}

}  // namespace

int main() {
  SampleProfile sp = MakeSampleProfile();
  STRLocusInfo info = MakeLocusInfo();
  FlankingClass flanking_class;
  flanking_class.SetGlobalParams(sp, 2000, false, false);
  flanking_class.SetLocusParams(info);
  flanking_class.SetCoverage(30);
  flanking_class.AddData(17);
  flanking_class.AddData(25);

  const int read_len = 100;
  const int motif_len = 2;
  const int ref_count = 12;
  const int ploidy = 2;
  double ll_without_model = 0.0;
  double ll_with_model = 0.0;

  HipStutterModel stutter_model(0.9, 0.01, 0.01, 0.9, 0.01, 0.01, motif_len);

  if (!flanking_class.GetClassLogLikelihood(18, 26, read_len, motif_len, ref_count, ploidy, nullptr, &ll_without_model)) {
    std::cerr << "flanking likelihood without stutter model failed" << std::endl;
    return 1;
  }
  if (!flanking_class.GetClassLogLikelihood(18, 26, read_len, motif_len, ref_count, ploidy, &stutter_model, &ll_with_model)) {
    std::cerr << "flanking likelihood with stutter model failed" << std::endl;
    return 1;
  }
  if (NearlyEqual(ll_without_model, ll_with_model, 1e-6)) {
    std::cerr << "flanking likelihood did not change when stutter model was provided" << std::endl;
    std::cerr << ll_without_model << " vs " << ll_with_model << std::endl;
    return 1;
  }

  return 0;
}
