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
#include "src/enclosing_class.h"

#include <math.h>
#include <iostream>
#include <algorithm>
#include <map>
using namespace std;

bool EnclosingClass::GetLogClassProb(const int32_t& allele,
				     const int32_t& read_len, const int32_t& motif_len,
				     double* log_class_prob) {
	int str_len = allele * motif_len;					// (L)
	double class_prob;
	if (double(2 * flank_len + str_len - 2 * read_len) == 0)
	{
	  cerr << "Divide by Zero prevented!" << endl;
	  *log_class_prob = NEG_INF;
	  return true;
	}
	if (read_len <= str_len)
		class_prob = 0;
	else
		class_prob = double(read_len - str_len) / double(2 * flank_len + str_len - 2 * read_len);

	if (class_prob > 0){
		*log_class_prob = log(class_prob);
		return true;
	}
	else if (class_prob == 0){
		*log_class_prob = NEG_INF;
		return true;
	}
	else
		return false;
}

// Original GetLogReadProb for backward compatibility if needed (e.g. base class calls)
bool EnclosingClass::GetLogReadProb(const int32_t& allele,
				    const int32_t& data,
				    const int32_t& read_len,
				    const int32_t& motif_len,
				    const int32_t& ref_count,
				    double* log_allele_prob) {
	double delta = data - allele;
	double allele_prob;
	if (delta == 0)
		allele_prob = 1 - stutter_up - stutter_down;
	else if (delta > 0)
		allele_prob = stutter_up * stutter_p * pow(1.0 - stutter_p, delta - 1.0);
	else
		allele_prob = stutter_down * stutter_p * pow(1.0 - stutter_p, -delta - 1.0);
	if (allele_prob > 0){
		*log_allele_prob = log(allele_prob);
		return true;
	}
	else if (allele_prob == 0){
		*log_allele_prob = NEG_INF;
		return true;
	}
	else
		return false;
}

// Overloaded GetLogReadProb that uses the new stutter model
bool EnclosingClass::GetLogReadProb(const int32_t& allele,
				    const int32_t& data,
				    const int32_t& motif_len,
				    const HipStutterModel* stutter_model,
				    double* log_allele_prob) {
	if (stutter_model == nullptr) {
		// Fallback to old model if no new model is provided
		return GetLogReadProb(allele, data, 0, motif_len, 0, log_allele_prob);
	}
	*log_allele_prob = stutter_model->log_stutter_pmf(allele * motif_len, data * motif_len);
	return true;
}

// Overloaded GetClassLogLikelihood that uses the new stutter model
bool EnclosingClass::GetClassLogLikelihood(const int32_t& allele1,
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
    double log_class_prob1, log_read_prob1;
    if (!GetLogClassProb(allele1, read_len, motif_len, &log_class_prob1)) return false;
    if (!GetLogReadProb(allele1, *data_it, motif_len, stutter_model, &log_read_prob1)) return false;
    a1_ll = log_class_prob1 + log_read_prob1;

    double log_class_prob2, log_read_prob2;
    if (!GetLogClassProb(allele2, read_len, motif_len, &log_class_prob2)) return false;
    if (!GetLogReadProb(allele2, *data_it, motif_len, stutter_model, &log_read_prob2)) return false;
    a2_ll = log_class_prob2 + log_read_prob2;

    if (ploidy == 2){
      *class_ll += fast_log_sum_exp(log(allele1_weight_)+a1_ll, log(allele2_weight_)+a2_ll);
    }
    else if (ploidy == 1){
      *class_ll += log(allele1_weight_) + a1_ll;
    }
  }
  return true;
}

bool EnclosingClass::GetGridBoundaries(int32_t* min_allele, int32_t* max_allele) {
  if (read_class_data_.empty()) return false;
  std::vector<int32_t>::iterator itmin = std::min_element(read_class_data_.begin(), read_class_data_.end());
  if (*itmin < *min_allele) {
    *min_allele = *itmin;
  }
  std::vector<int32_t>::iterator itmax = std::max_element(read_class_data_.begin(), read_class_data_.end());
  if (*itmax > *max_allele) {
    *max_allele = *itmax;
  }
  return true;
}

bool EnclosingClass::ExtractAllEnclosingAlleles(std::vector<int> *alleles) const {
  if (alleles == NULL) {
    return false;
  }
  size_t before_size = alleles->size();
  for (std::vector<int32_t>::const_iterator data_it = read_class_data_.begin();
       data_it != read_class_data_.end();
       ++data_it) {
    alleles->push_back(*data_it);
  }
  return alleles->size() > before_size;
}

int32_t EnclosingClass::GetAlleleCount(const int32_t& allele) const {
  int32_t count = 0;
  for (std::vector<int32_t>::const_iterator data_it = read_class_data_.begin();
       data_it != read_class_data_.end();
       ++data_it) {
    if (*data_it == allele) {
      count++;
    }
  }
  return count;
}

bool EnclosingClass::ExtractEnclosingAlleles(std::vector<int> *alleles){
    if (alleles == NULL) {
      return false;
    }
    size_t before_size = alleles->size();
    std::map<int32_t, int32_t> allele_repeats;

	for (std::vector<int32_t>::iterator data_it = this->read_class_data_.begin();
       data_it != this->read_class_data_.end();
       data_it++) { 
       	if (allele_repeats.find(*data_it) == allele_repeats.end()){
       		allele_repeats[*data_it] = 1;
       	}
       	else{
       		allele_repeats[*data_it]++;
       	}
  	}
	// Now refill read_class_data_ only with repeated enclosing reads
	read_class_data_.clear();
  	for (map<int32_t, int32_t>::iterator it = allele_repeats.begin(); it != allele_repeats.end(); it++){
  		if (it->second >= 2){
		    (*alleles).push_back(it->first);
  		    //cerr << it->first << "\t" << it -> second << endl;
		    for (int i = 0; i < it->second; i++){
		      read_class_data_.push_back(it->first);
		    }
  		}
  	}
    return alleles->size() > before_size;
}
