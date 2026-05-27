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

#include <cerrno>
#include <climits>
#include <cmath>
#include <cctype>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <sstream>

#include "src/genotyper.h"
#include "src/mathops.h"
#include "src/stringops.h"

#include <set>
using namespace std;

namespace {

std::string TrimAscii(const std::string& value) {
  size_t begin = 0;
  while (begin < value.size() &&
         std::isspace(static_cast<unsigned char>(value[begin]))) {
    begin++;
  }
  size_t end = value.size();
  while (end > begin &&
         std::isspace(static_cast<unsigned char>(value[end - 1]))) {
    end--;
  }
  return value.substr(begin, end - begin);
}

bool ParseStrictInt32(const std::string& raw_value,
                      const std::string& field_name,
                      int32_t* parsed_value,
                      std::string* error_message) {
  const std::string value = TrimAscii(raw_value);
  if (value.empty()) {
    *error_message = field_name + " is empty";
    return false;
  }
  char* parse_end = NULL;
  errno = 0;
  const long parsed = std::strtol(value.c_str(), &parse_end, 10);
  if (parse_end == value.c_str() || *parse_end != '\0' ||
      errno == ERANGE || parsed < INT_MIN || parsed > INT_MAX) {
    *error_message = field_name + " is not a valid integer: " + raw_value;
    return false;
  }
  *parsed_value = static_cast<int32_t>(parsed);
  return true;
}

bool ParseStrictDouble(const std::string& raw_value,
                       const std::string& field_name,
                       double* parsed_value,
                       std::string* error_message) {
  const std::string value = TrimAscii(raw_value);
  if (value.empty()) {
    *error_message = field_name + " is empty";
    return false;
  }
  char* parse_end = NULL;
  errno = 0;
  const double parsed = std::strtod(value.c_str(), &parse_end);
  if (parse_end == value.c_str() || *parse_end != '\0' ||
      errno == ERANGE || !std::isfinite(parsed)) {
    *error_message = field_name + " is not a valid finite number: " + raw_value;
    return false;
  }
  *parsed_value = parsed;
  return true;
}

bool ValidateExternalStutterModel(const std::string& chrom,
                                  int32_t start,
                                  int32_t end,
                                  int32_t period,
                                  const std::string& motif,
                                  double in_geom,
                                  double in_up,
                                  double in_down,
                                  double out_geom,
                                  double out_up,
                                  double out_down,
                                  std::string* error_message) {
  if (TrimAscii(chrom).empty()) {
    *error_message = "chrom cannot be empty";
    return false;
  }
  if (start < 0 || end < start) {
    *error_message = "start/end are invalid";
    return false;
  }
  if (period < 1 || period > 9) {
    *error_message = "period must be in [1, 9]";
    return false;
  }
  if (motif.empty()) {
    *error_message = "motif cannot be empty";
    return false;
  }
  if (static_cast<int32_t>(motif.size()) != period) {
    *error_message = "motif length must equal period";
    return false;
  }
  if (!(in_geom > 0.0 && in_geom < 1.0) ||
      !(out_geom > 0.0 && out_geom < 1.0)) {
    *error_message = "in_geom and out_geom must be in (0, 1)";
    return false;
  }
  if (!(in_up > 0.0 && in_down > 0.0 &&
        out_up > 0.0 && out_down > 0.0)) {
    *error_message = "stutter direction probabilities must be > 0";
    return false;
  }
  if (in_up + in_down + out_up + out_down >= 1.0) {
    *error_message = "stutter direction probabilities must sum to < 1";
    return false;
  }
  if (in_up > in_down) {
    *error_message = "in-frame P_UP must not exceed P_DOWN";
    return false;
  }
  return true;
}

}  // namespace

Genotyper::Genotyper(RefGenome& _refgenome,
		     Options& _options,
		     SampleInfo& _sample_info,
		     STRInfo& _str_info) {
  refgenome = &_refgenome;
  options = &_options;
  sample_info = &_sample_info;
  str_info = &_str_info;
  read_extractor = new ReadExtractor(_options, *sample_info);
  std::set<std::string> rg_samples = sample_info->GetSamples();
  for (std::set<std::string>::iterator it=rg_samples.begin();
       it != rg_samples.end(); it++) {
    SampleProfile sp;
    if (sample_info->GetSampleProfile(*it, &sp)) {
      sample_likelihood_maximizers[*it] = new LikelihoodMaximizer(_options, sp, sample_info->GetReadLength(), sample_info->GetSampleSex(*it));
    } else {
      PrintMessageDieOnError("Could not find sample profile for " + *it, M_ERROR, false);
    }
  }
}

bool Genotyper::SetFlanks(Locus* locus) {
  if (!refgenome->GetSequence(locus->chrom,
			      locus->start-options->realignment_flanklen-1,
			      locus->start-2,
			      &locus->pre_flank)) {
    return false;
  }
  if (!refgenome->GetSequence(locus->chrom,
			      locus->end,
			      locus->end+options->realignment_flanklen-1,
			      &locus->post_flank)) {
    return false;
  }
  return true;
}

std::string Genotyper::BuildLocusStutterKey(const std::string& chrom,
                                            int32_t start,
                                            int32_t end,
                                            int32_t period,
                                            const std::string& motif) const {
  std::stringstream key;
  key << chrom
      << ":" << start
      << ":" << end
      << ":" << period
      << ":" << lowercase(motif);
  return key.str();
}

std::string Genotyper::GetLocusStutterKey(const Locus& locus) const {
  return BuildLocusStutterKey(locus.chrom, locus.start, locus.end, locus.period, locus.motif);
}

bool Genotyper::SetGGL(Locus& locus, const std::string& samp) {
  std::map<std::pair<int32_t, int32_t>, double> gridlik;
  for (int32_t a1=locus.grid_min_allele; a1<=locus.grid_max_allele; a1++) {
    for (int32_t a2=a1; a2<=locus.grid_max_allele; a2++) {
      double gt_ll;
      if (!sample_likelihood_maximizers[samp]->GetGenotypeNegLogLikelihood(a1, a2, false, &gt_ll)) {
	return false;
      }
      std::pair<int32_t, int32_t> gt(a1, a2);
      gridlik[gt] = gt_ll/log(10)*-1;
    }
  }
  locus.grid_likelihoods[samp] = gridlik;
  return true;
}

bool Genotyper::LoadExternalStutterModels(const std::string& model_path) {
  if (model_path.empty()) return true;

  std::ifstream model_stream(model_path.c_str());
  if (!model_stream.good()) {
    PrintMessageDieOnError("Could not open external stutter model file: " + model_path, M_ERROR, false);
    return false;
  }

  std::string line;
  std::map<std::string, size_t> header_index;
  bool saw_header = false;
  int loaded_models = 0;
  int line_number = 0;
  const std::vector<std::string> required_columns = {
    "chrom", "start", "end", "period", "motif",
    "in_geom", "in_up", "in_down", "out_geom", "out_up", "out_down"
  };

  while (std::getline(model_stream, line)) {
    line_number++;
    if (line.empty()) continue;
    if (!line.empty() && line[0] == '#') continue;

    std::vector<std::string> fields;
    split_by_delim(line, '\t', fields);
    if (fields.empty()) continue;

    if (!saw_header && lowercase(TrimAscii(fields[0])) == "chrom") {
      for (size_t i = 0; i < fields.size(); i++) {
        header_index[lowercase(TrimAscii(fields[i]))] = i;
      }
      for (size_t i = 0; i < required_columns.size(); i++) {
        if (header_index.count(required_columns[i]) == 0) {
          PrintMessageDieOnError("External stutter model header is missing required column: " + required_columns[i], M_ERROR, false);
          return false;
        }
      }
      saw_header = true;
      continue;
    }

    if (fields.size() < 11 && header_index.empty()) {
      PrintMessageDieOnError("External stutter model row has fewer than 11 columns: " + line, M_ERROR, false);
      return false;
    }

    bool row_error = false;
    auto get_field = [&](const std::string& name, size_t fallback_index) -> std::string {
      std::map<std::string, size_t>::const_iterator it = header_index.find(name);
      size_t idx = fallback_index;
      if (it != header_index.end()) idx = it->second;
      if (idx >= fields.size()) {
        row_error = true;
        return "";
      }
      return TrimAscii(fields[idx]);
    };

    const std::string chrom = get_field("chrom", 0);
    int32_t start = 0;
    int32_t end = 0;
    int32_t period = 0;
    double in_geom = 0.0;
    double in_up = 0.0;
    double in_down = 0.0;
    double out_geom = 0.0;
    double out_up = 0.0;
    double out_down = 0.0;
    const std::string motif = lowercase(get_field("motif", 4));
    if (row_error) {
      PrintMessageDieOnError("Malformed external stutter model row in " + model_path + ": " + line, M_ERROR, false);
      return false;
    }
    std::string parse_error;
    if (!ParseStrictInt32(get_field("start", 1), "start", &start, &parse_error) ||
        !ParseStrictInt32(get_field("end", 2), "end", &end, &parse_error) ||
        !ParseStrictInt32(get_field("period", 3), "period", &period, &parse_error) ||
        !ParseStrictDouble(get_field("in_geom", 5), "in_geom", &in_geom, &parse_error) ||
        !ParseStrictDouble(get_field("in_up", 6), "in_up", &in_up, &parse_error) ||
        !ParseStrictDouble(get_field("in_down", 7), "in_down", &in_down, &parse_error) ||
        !ParseStrictDouble(get_field("out_geom", 8), "out_geom", &out_geom, &parse_error) ||
        !ParseStrictDouble(get_field("out_up", 9), "out_up", &out_up, &parse_error) ||
        !ParseStrictDouble(get_field("out_down", 10), "out_down", &out_down, &parse_error)) {
      std::stringstream err;
      err << "Invalid external stutter model row in " << model_path
          << " at line " << line_number << ": " << parse_error;
      PrintMessageDieOnError(err.str(), M_ERROR, false);
      return false;
    }
    if (!ValidateExternalStutterModel(chrom, start, end, period, motif,
                                      in_geom, in_up, in_down,
                                      out_geom, out_up, out_down,
                                      &parse_error)) {
      std::stringstream err;
      err << "Invalid external stutter model parameters in " << model_path
          << " at line " << line_number << ": " << parse_error;
      PrintMessageDieOnError(err.str(), M_ERROR, false);
      return false;
    }

    std::string locus_id = BuildLocusStutterKey(chrom, start, end, period, motif);
    std::map<std::string, HipStutterModel*>::iterator existing_model = locus_stutter_models.find(locus_id);
    if (existing_model != locus_stutter_models.end()) {
      PrintMessageDieOnError("Duplicate external stutter model for locus " + locus_id + "; keeping the last row", M_WARNING, options->quiet);
      delete existing_model->second;
      existing_model->second = nullptr;
    }
    locus_stutter_models[locus_id] = new HipStutterModel(in_geom, in_up, in_down,
                                                         out_geom, out_up, out_down,
                                                         period);
    loaded_models++;
  }

  if (options->verbose) {
    std::stringstream ss;
    ss << "Loaded " << loaded_models << " external stutter models from " << model_path;
    PrintMessageDieOnError(ss.str(), M_PROGRESS, options->quiet);
  }
  return true;
}

bool Genotyper::LearnStutterModels(BamCramMultiReader* bamreader, std::vector<Locus*>& loci) {
    if (options->verbose) {
        PrintMessageDieOnError("Learning locus-specific stutter models...", M_PROGRESS, options->quiet);
    }

    for (auto const& locus : loci) {
        const std::string locus_id = GetLocusStutterKey(*locus);
        std::map<std::string, HipStutterModel*>::const_iterator existing_model = locus_stutter_models.find(locus_id);
        if (existing_model != locus_stutter_models.end() && existing_model->second != nullptr) {
            if (options->verbose) {
                PrintMessageDieOnError("\tSkipping stutter model learning for locus " + locus_id + " (external model already loaded)", M_PROGRESS, options->quiet);
            }
            continue;
        }

        if (!SetFlanks(locus)) {
            PrintMessageDieOnError("Failed to set flanking sequence for stutter learning at locus " + locus->chrom + ":" + std::to_string(locus->start), M_WARNING, options->quiet);
            continue;
        }

        // We need a temporary map of likelihood maximizers just for extracting reads for this locus.
        std::map<std::string, LikelihoodMaximizer*> temp_lms;
        std::set<std::string> rg_samples = sample_info->GetSamples();
        for (const auto& samp : rg_samples) {
            SampleProfile sp;
            if (sample_info->GetSampleProfile(samp, &sp)) {
                temp_lms[samp] = new LikelihoodMaximizer(*options, sp, sample_info->GetReadLength(), sample_info->GetSampleSex(samp));
            } else {
                for (auto const& item : temp_lms) {
                    delete item.second;
                }
                PrintMessageDieOnError("Could not find sample profile for " + samp, M_ERROR, false);
                return false;
            }
        }

        // Extract reads for the current locus
        if (!read_extractor->ExtractReads(bamreader, *locus, options->regionsize, options->min_match, temp_lms)) {
            PrintMessageDieOnError("Failed to extract reads for stutter learning at locus " + locus->chrom + ":" + std::to_string(locus->start), M_WARNING, options->quiet);
            // Cleanup
            for(auto const& [key, val] : temp_lms) { delete val; }
            continue; // Skip to next locus
        }

        // Aggregate enclosing read alleles per sample for sample-aware EM training.
        std::vector< std::vector<int> > per_sample_enclosing_alleles;
        size_t total_enclosing_reads = 0;
        for (const auto& samp : rg_samples) {
            std::vector<int> sample_alleles;
            temp_lms[samp]->enclosing_class_.ExtractAllEnclosingAlleles(&sample_alleles);
            total_enclosing_reads += sample_alleles.size();
            per_sample_enclosing_alleles.push_back(sample_alleles);
        }

        // Cleanup the temporary LMs
        for(auto const& [key, val] : temp_lms) { delete val; }

        // Check if we have enough data to learn a model
        if (total_enclosing_reads < 20) { // Heuristic threshold
            if (options->verbose) {
                PrintMessageDieOnError("\tSkipping stutter model learning for locus " + locus->chrom + ":" + std::to_string(locus->start) + " (not enough enclosing reads)", M_PROGRESS, options->quiet);
            }
            continue;
        }

        // Learn the model
        HipEMLearner em_learner(per_sample_enclosing_alleles, locus->motif.length());
        bool success = em_learner.train(20, 0.01);

        std::map<std::string, HipStutterModel*>::iterator stored_model = locus_stutter_models.find(locus_id);
        if (stored_model != locus_stutter_models.end()) {
            delete stored_model->second;
            stored_model->second = nullptr;
        }
        if (success) {
            HipStutterModel* learned_model = em_learner.get_stutter_model()->copy();
            const double in_down = learned_model->get_parameter(true, 'D');
            const double in_up = learned_model->get_parameter(true, 'U');
            const bool suspicious_direction = in_up > in_down;

            if (options->verbose) {
                PrintMessageDieOnError("\tSuccessfully learned stutter model for locus " + locus_id, M_PROGRESS, options->quiet);
                std::ostringstream model_ss;
                model_ss << *learned_model;
                PrintMessageDieOnError(model_ss.str(), M_PROGRESS, options->quiet);
            }

            if (suspicious_direction) {
                if (options->verbose) {
                    std::ostringstream warn_ss;
                    warn_ss << "\tRejecting learned stutter model for locus " << locus_id
                            << " because in-frame P_UP (" << in_up
                            << ") exceeded P_DOWN (" << in_down
                            << "); falling back to default GangSTR behavior";
                    PrintMessageDieOnError(warn_ss.str(), M_PROGRESS, options->quiet);
                }
                delete learned_model;
                locus_stutter_models[locus_id] = nullptr;
            } else {
                locus_stutter_models[locus_id] = learned_model;
            }
        } else {
            PrintMessageDieOnError("\tFailed to learn stutter model for locus " + locus_id, M_WARNING, options->quiet);
            locus_stutter_models[locus_id] = nullptr;
        }
    }

    return true;
}

bool Genotyper::ProcessLocus(BamCramMultiReader* bamreader, Locus* locus) {
  int32_t read_len = sample_info->GetReadLength();
  // Load preflank and postflank to locus
  if (options->verbose) {
    PrintMessageDieOnError("\tSetting flanking regions", M_PROGRESS, options->quiet);
  }
  if (!SetFlanks(locus)) {
    return false;
  }

  for (std::map<std::string, LikelihoodMaximizer*>::iterator it = sample_likelihood_maximizers.begin();
       it != sample_likelihood_maximizers.end(); it++) {
    it->second->Reset();
  }

  // Infer GC bin
  int gcbin = -1;
  if (options->model_gc_cov) {
    std::string seq;
    if (refgenome->GetSequence(locus->chrom,
			       int((locus->start+locus->end)/2-options->gc_region_len/2),
			       int((locus->start+locus->end)/2+options->gc_region_len/2),
			       &seq)) {
      float gc = GetGC(seq);
      gcbin = int(floor(gc/options->gc_bin_size));
    }
  }
  // Load all read data
  if (options->verbose) {
    PrintMessageDieOnError("\tLoading read data", M_PROGRESS, options->quiet);
  }
  if (!read_extractor->ExtractReads(bamreader, *locus, options->regionsize,
				    options->min_match, sample_likelihood_maximizers)) {
    return false;
  }

  // Get the learned stutter model for this locus
  std::string locus_id = GetLocusStutterKey(*locus);
  const HipStutterModel* stutter_model = nullptr;
  if (locus_stutter_models.count(locus_id)) {
      stutter_model = locus_stutter_models.at(locus_id);
  }
  if (options->verbose && options->stutter_mode == "external" && stutter_model == nullptr) {
      PrintMessageDieOnError("\tNo external stutter model matched " + locus_id + "; falling back to default GangSTR stutter behavior", M_WARNING, options->quiet);
  }

  std::set<std::string> rg_samples = sample_info->GetSamples();
  // First set grid size
  int32_t sample_min_allele, sample_max_allele;
  int32_t min_allele = 100000;
  int32_t max_allele = 0;
  int32_t ref_count = (int32_t)((locus->end-locus->start+1)/locus->motif.size());
  for (std::set<std::string>::iterator it = rg_samples.begin();
       it != rg_samples.end(); it++) {
    const std::string samp = *it;
    // Default to original GangSTR behavior for the provisional per-sample call.
    sample_likelihood_maximizers[samp]->SetStutterModel(nullptr);

    if (gcbin != -1) {
      sample_likelihood_maximizers[samp]->SetLocusParams(str_info->GetSTRInfo(locus->chrom, locus->start),
							 sample_info->GetGCCoverage(samp, gcbin),
							 sample_info->GetReadLength(), (int32_t)(locus->motif.size()),
							 ref_count, locus->chrom);
    } else {
      sample_likelihood_maximizers[samp]->SetLocusParams(str_info->GetSTRInfo(locus->chrom, locus->start),
							 sample_info->GetCoverage(samp),
							 sample_info->GetReadLength(), (int32_t)(locus->motif.size()),
							 ref_count, locus->chrom);
    }
    if (!sample_likelihood_maximizers[samp]->InferGridSize() ) {
      PrintMessageDieOnError("Error inferring grid size", M_PROGRESS, options->quiet);
      return false;
    }
    sample_likelihood_maximizers[samp]->GetGridSize(&sample_min_allele, &sample_max_allele);
    if (sample_max_allele > max_allele) max_allele = sample_max_allele;
    if (sample_min_allele < min_allele) min_allele = sample_min_allele;
  }
  // Set locus info to output to VCF
  locus->grid_min_allele = min_allele;
  locus->grid_max_allele = max_allele;
  locus->expansion_threshold = str_info->GetExpansionThreshold(locus->chrom, locus->start);
  locus->stutter_up = str_info->GetStutterUp(locus->chrom, locus->start);
  locus->stutter_down = str_info->GetStutterDown(locus->chrom, locus->start);
  locus->stutter_p = str_info->GetStutterP(locus->chrom, locus->start);
  // Maximize the likelihood
  if (options->verbose) {
    PrintMessageDieOnError("\tMaximizing likelihood", M_PROGRESS, options->quiet);
  }
  int32_t allele1, allele2;
  double min_negLike, lob1, lob2, hib1, hib2, q_score;
  double a1_se, a2_se;
  bool resampled = false;
  for (std::set<std::string>::iterator it = rg_samples.begin();
       it != rg_samples.end(); it++) {
    const std::string samp = *it;
    locus->called[samp] = false;
    sample_likelihood_maximizers[samp]->SetGridSize(min_allele, max_allele);
    sample_likelihood_maximizers[samp]->SetStutterEnclosingWeightScale(1.0);
    sample_likelihood_maximizers[samp]->SetStutterFlankingWeightScale(1.0);
    try {
      sample_likelihood_maximizers[samp]->SetStutterModel(nullptr);
      int32_t base_allele1 = 0, base_allele2 = 0;
      double base_min_negLike = 0.0;
      if (!sample_likelihood_maximizers[samp]->OptimizeLikelihood(resampled, options->ploidy,
								  0,
								  locus->offtarget_share,
								  &base_allele1,
								  &base_allele2,
								  &base_min_negLike)) {
	continue;
      }
      allele1 = base_allele1;
      allele2 = base_allele2;
      min_negLike = base_min_negLike;

      bool apply_sample_stutter = false;
      bool use_large_probe_boost = false;
      bool use_upper_large_probe_boost = false;
      bool use_spanning_surface_boost = false;
      bool use_frr_surface_boost = false;
      if (options->verbose && stutter_model != nullptr) {
        stringstream gate_msg;
        gate_msg << "\tStutter Gate [" << samp << "]: "
                 << sample_likelihood_maximizers[samp]->GetStutterGateDebugSummary(
                        base_allele1,
                        base_allele2,
                        options->stutter_gate_max_adjacent_diff,
                        options->stutter_gate_min_adjacent_count);
        PrintMessageDieOnError(gate_msg.str(), M_PROGRESS, options->quiet);
      }
      if (stutter_model != nullptr) {
        use_upper_large_probe_boost =
            sample_likelihood_maximizers[samp]->ShouldApplyStutterModelForUpperLargeSideProbeCall(
                base_allele1,
                base_allele2,
                options->stutter_gate_min_adjacent_count);
        use_large_probe_boost =
            sample_likelihood_maximizers[samp]->ShouldApplyStutterModelForLargeSideProbeCall(
                base_allele1,
                base_allele2,
                options->stutter_gate_min_adjacent_count);
        use_spanning_surface_boost =
            sample_likelihood_maximizers[samp]->ShouldApplyStutterModelForSpanningCall(
                base_allele1,
                base_allele2);
        use_frr_surface_boost =
            sample_likelihood_maximizers[samp]->ShouldApplyStutterModelForFRRCall(
                base_allele1,
                base_allele2);
      }
      if (stutter_model != nullptr &&
          sample_likelihood_maximizers[samp]->ShouldApplyStutterModelForCall(
              base_allele1,
              base_allele2,
              options->stutter_gate_max_adjacent_diff,
              options->stutter_gate_min_adjacent_count)) {
        if (options->verbose) {
          stringstream apply_msg;
          apply_msg << "\tApplying learned stutter model for sample " << samp;
          if (use_upper_large_probe_boost) {
            apply_msg << " with upper-large flanking boost";
          } else if (use_large_probe_boost) {
            apply_msg << " with large-probe flanking boost";
          }
          if (use_spanning_surface_boost) {
            apply_msg << " with spanning-surface rerank";
          }
          if (use_frr_surface_boost) {
            apply_msg << " with FRR-surface rerank";
          }
          PrintMessageDieOnError(apply_msg.str(), M_PROGRESS, options->quiet);
        }
        int32_t model_allele1 = 0, model_allele2 = 0;
        double model_min_negLike = 0.0;
        double enclosing_scale = 1.0;
        double flanking_scale = 1.0;
        if (use_large_probe_boost) {
          flanking_scale = std::max(flanking_scale, 1.10);
        }
        if (use_upper_large_probe_boost) {
          flanking_scale = std::max(flanking_scale, 1.25);
          enclosing_scale = std::max(enclosing_scale, 1.05);
        }
        if (use_spanning_surface_boost) {
          enclosing_scale = std::max(enclosing_scale, 1.08);
          flanking_scale = std::max(flanking_scale, 1.15);
        }
        if (use_frr_surface_boost) {
          enclosing_scale = std::max(enclosing_scale, 1.10);
          flanking_scale = std::max(flanking_scale, 1.20);
        }
        sample_likelihood_maximizers[samp]->SetStutterEnclosingWeightScale(enclosing_scale);
        sample_likelihood_maximizers[samp]->SetStutterFlankingWeightScale(flanking_scale);
        sample_likelihood_maximizers[samp]->SetStutterModel(stutter_model);
        if (sample_likelihood_maximizers[samp]->OptimizeLikelihood(resampled, options->ploidy,
								   0,
								   locus->offtarget_share,
								   &model_allele1,
								   &model_allele2,
								   &model_min_negLike)) {
          allele1 = model_allele1;
          allele2 = model_allele2;
          min_negLike = model_min_negLike;
          apply_sample_stutter = true;
        } else {
          sample_likelihood_maximizers[samp]->SetStutterEnclosingWeightScale(1.0);
          sample_likelihood_maximizers[samp]->SetStutterFlankingWeightScale(1.0);
          sample_likelihood_maximizers[samp]->SetStutterModel(nullptr);
        }
      }
      if (!apply_sample_stutter) {
        sample_likelihood_maximizers[samp]->SetStutterEnclosingWeightScale(1.0);
        sample_likelihood_maximizers[samp]->SetStutterFlankingWeightScale(1.0);
        sample_likelihood_maximizers[samp]->SetStutterModel(nullptr);
      }

      std::vector<double> sample_prob_vec;
      if (!sample_likelihood_maximizers[samp]->GetExpansionProb(&sample_prob_vec, locus->expansion_threshold)) {
	sample_prob_vec.clear();
	sample_prob_vec.push_back(-1.0);
	sample_prob_vec.push_back(-1.0);
	sample_prob_vec.push_back(-1.0);
      }
      if (!options->skip_qscore){
	if (!sample_likelihood_maximizers[samp]->GetQScore(allele1, allele2, &q_score)){ 
	  q_score = -1;
	  PrintMessageDieOnError("\tProblem setting quality scores", M_WARNING, options->quiet);
	}
      }
      else{
	q_score = -1;
      }
      locus->q_scores[samp] = q_score;
      locus->expansion_probs[samp] = sample_prob_vec;
      locus->allele1[samp] = allele1;
      locus->allele2[samp] = allele2;
      locus->min_neg_lik[samp] = min_negLike;
      locus->enclosing_reads[samp] = sample_likelihood_maximizers[samp]->GetEnclosingDataSize();
      locus->spanning_reads[samp] = sample_likelihood_maximizers[samp]->GetSpanningDataSize();
      locus->frr_reads[samp] = sample_likelihood_maximizers[samp]->GetFRRDataSize();
      locus->flanking_reads[samp] = sample_likelihood_maximizers[samp]->GetFlankingDataSize();
      locus->enclosing_reads_dict[samp] = sample_likelihood_maximizers[samp]->GetEnclosingReadDictStr();
      locus->flanking_reads_dict[samp] = sample_likelihood_maximizers[samp]->GetFlankingReadDictStr();


      locus->depth[samp] = sample_likelihood_maximizers[samp]->GetReadPoolSize();
      locus->called[samp] = true;
      if (allele1 <= 0 and allele2 <= 0){
	PrintMessageDieOnError("\tProblem maximizing likelihood. Skipping locus", M_WARNING, options->quiet);
	locus->called[samp] = false;
	continue;
      }
      if (options->include_ggl && !SetGGL(*locus, samp)) {
	PrintMessageDieOnError("\tProblem setting genotype likelihoods", M_WARNING, options->quiet);
      }
      if (options->num_boot_samp > 0){
	if (options->verbose) {
	  PrintMessageDieOnError("\tGetting confidence intervals", M_PROGRESS, options->quiet);
	}
	try{
	  if (!sample_likelihood_maximizers[samp]->GetConfidenceInterval(allele1, allele2, *locus,
									 &lob1, &hib1, &lob2, &hib2, &a1_se, &a2_se)) {
	    locus->called[samp] = false;
	    continue;
	  }
	  locus->lob1[samp] = lob1;
	  locus->lob2[samp] = lob2;
	  locus->hib1[samp] = hib1;
	  locus->hib2[samp] = hib2;
	  locus->a1_se[samp] = a1_se;
	  locus->a2_se[samp] = a2_se;
	  
	  stringstream msg;
	  msg<<"\tGenotyper Results:  "<<allele1<<", "<<allele2<<"\tlikelihood = "<<min_negLike;
	  PrintMessageDieOnError(msg.str(), M_PROGRESS, options->quiet);
	  if (options->verbose) {
	    msg.clear();
	    msg.str(std::string());
	    msg<<"\tSmall Allele Bound: ["<<lob1<<", "<<hib1<<"]";
	    PrintMessageDieOnError(msg.str(), M_PROGRESS, options->quiet);
	    msg.clear();
	    msg.str(std::string());
	    msg<<"\tLarge Allele Bound: ["<<lob2<<", "<<hib2<<"]";
	    PrintMessageDieOnError(msg.str(), M_PROGRESS, options->quiet);
	  }
	}
	catch (std::exception &exc){
	  if (options->verbose) {
	    stringstream msg;
	    msg<<"\tEncountered error("<< exc.what() <<") in likelihood maximization for confidence interval. Skipping locus";
	    PrintMessageDieOnError(msg.str(), M_PROGRESS, options->quiet);
	  }
	  locus->called[samp] = false;
	}
      }
    }
    catch (std::exception &exc){
      if (options->verbose) {
	stringstream msg;
	msg<<"\tEncountered error("<< exc.what() <<") in likelihood maximization. Skipping locus";
	PrintMessageDieOnError(msg.str(), M_PROGRESS, options->quiet);
      }
      locus->called[samp] = false;
    }
  }
  return true;
}

void Genotyper::Debug(BamCramMultiReader* bamreader) {
  cerr << "testing refgenome" << endl;
  std::string seq;
  refgenome->GetSequence("3", 63898261, 63898360, &seq);
  cerr << seq << endl;
  cerr << "testing bam" << endl;
  bamreader->SetRegion("1", 0, 10000);
  BamAlignment aln;
  if (bamreader->GetNextAlignment(aln)) { // Requires SetRegion was called
    std::string testread = aln.QueryBases();
    cerr << testread << endl;
  } else {
    cerr << "testing bam failed" << endl;
  }
  cerr << "testing GSL" << endl;
  double x = TestGSL();
  cerr << "gsl_ran_gaussian_pdf(0, 1) " << x << endl;
  //  double y = TestNLOPT();
}

Genotyper::~Genotyper() {
  delete read_extractor;
  for (std::map<std::string, LikelihoodMaximizer*>::iterator it = sample_likelihood_maximizers.begin();
       it != sample_likelihood_maximizers.end(); it++) {
    delete it->second;
  }
  for (std::map<std::string, HipStutterModel*>::iterator it = locus_stutter_models.begin();
       it != locus_stutter_models.end(); it++) {
    delete it->second;
  }
}
