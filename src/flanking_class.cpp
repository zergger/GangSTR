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

#include "src/mathops.h"
#include "src/flanking_class.h"
#include <algorithm>
#include <iostream>
#include <iomanip>
#include <math.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_cdf.h>
using namespace std;


bool FlankingClass::GetAlleleLogLikelihood(const int32_t& allele,
				   const int32_t& data,
				   const int32_t& motif_len,
                   const HipStutterModel* stutter_model,
				   double* allele_ll){
    if (stutter_model == nullptr) {
        *allele_ll = NEG_INF; // Or some other default behavior
        return true;
    }
    *allele_ll = stutter_model->log_stutter_pmf(allele * motif_len, data * motif_len);
	return true;
}


bool FlankingClass::GetClassLogLikelihood(const int32_t& allele1,
				      const int32_t& allele2,
				      const int32_t& read_len, const int32_t& motif_len,
				      const int32_t& ref_count, const int32_t& ploidy,
                      const HipStutterModel* stutter_model,
				      double* class_ll) {
  *class_ll = 0;
  double samp_log_likelihood, a1_ll, a2_ll;
  for (std::vector<int32_t>::iterator data_it = read_class_data_.begin();
       data_it != read_class_data_.end();
       data_it++) {
    if (!FlankingClass::GetAlleleLogLikelihood(allele1, *data_it, motif_len, stutter_model, &a1_ll)) {
      return false;
    }
    if (!FlankingClass::GetAlleleLogLikelihood(allele2, *data_it, motif_len, stutter_model, &a2_ll)) {
      return false;
    }
    if (ploidy == 2){
      *class_ll += fast_log_sum_exp(log(allele1_weight_)+a1_ll, log(allele2_weight_)+a2_ll);
  	}
    else if (ploidy == 1){
      *class_ll += log(allele1_weight_) + a1_ll;
    }
  }
  return true;
}

bool FlankingClass::GetGridBoundaries(int32_t* min_allele, int32_t* max_allele) {
  if (read_class_data_.empty()) return false;
  // For flanking, only care about pushing the max allele size
  std::vector<int32_t>::iterator itmax = std::max_element(read_class_data_.begin(), read_class_data_.end());
  if (*itmax > *max_allele) {
    *max_allele = *itmax;
  }
  return true;
}

