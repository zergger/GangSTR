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

#ifndef SRC_ENCLOSING_CLASS_H__
#define SRC_ENCLOSING_CLASS_H__

#include "src/read_class.h"
#include "src/hipstr_models/hip_stutter_model.h"

/*
  Type of ReadClass

  Enclosing reads have a single read completely spanning the STR

 */
class EnclosingClass: public ReadClass {
 public:
  bool GetLogClassProb(const int32_t& allele,
		       const int32_t& read_len, const int32_t& motif_len,
		       double* log_class_prob);
  bool GetClassLogLikelihood(const int32_t& allele1, const int32_t& allele2,
			     const int32_t& read_len, const int32_t& motif_len,
			     const int32_t& ref_count, const int32_t& ploidy,
			     const HipStutterModel* stutter_model,
			     double* class_ll);
  bool GetLogReadProb(const int32_t& allele, const int32_t& data,
		      const int32_t& motif_len,
		      const HipStutterModel* stutter_model,
		      double* log_allele_prob);
  bool GetGridBoundaries(int32_t* min_allele, int32_t* max_allele);
  // Function to extract all enclosing alleles present
  bool ExtractEnclosingAlleles(std::vector<int> *alleles);
private:
    bool GetLogReadProb(const int32_t& allele, const int32_t& data,
		      const int32_t& read_len,
		      const int32_t& motif_len,
		      const int32_t& ref_count,
		      double* log_allele_prob) override;
};

#endif  // SRC_ENCLOSING_CLASS_H__
