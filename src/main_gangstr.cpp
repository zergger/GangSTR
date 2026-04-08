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

#include <getopt.h>
#include <stdlib.h>
#include <iostream>
#include <set>
#include <sstream>

//#include "src/bam_reader.h"
#include "src/bam_info_extract.h"
#include "src/bam_io.h"
#include "src/common.h"
#include "src/genotyper.h"
#include "src/options.h"
#include "src/ref_genome.h"
#include "src/region_reader.h"
#include "src/sample_info.h"
#include "src/str_info.h"
#include "src/stringops.h"
#include "src/vcf_writer.h"
#include "GangSTRConfig.h"

using namespace std;

void show_help() {
  Options options;
  std::stringstream help_msg;
  help_msg << "\nUsage: GangSTR [OPTIONS] "
	   << "--bam <file1[,file2,...]> "
	   << "--ref <reference.fa> "
	   << "--regions <regions.bed> "
	   << "--out <outprefix> "
	   << "\n\n Required options:\n"
	   << "\t" << "--bam         <file.bam,[file2.bam]>" << "\t" << "Comma separated list of input BAM files" << "\n"
	   << "\t" << "--ref         <genome.fa>     " << "\t" << "FASTA file for the reference genome" << "\n"
	   << "\t" << "--regions     <regions.bed>   " << "\t" << "BED file containing TR coordinates" << "\n"
	   << "\t" << "--out         <outprefix>     " << "\t" << "Prefix to name output files" << "\n"
	   << "\n Additional general options:\n"
	   << "\t" << "--targeted                    " << "\t" << "Targeted mode" << "\n"
	   << "\t" << "--chrom                       " << "\t" << "Only genotype regions on this chromosome" << "\n"
           << "\t" << "--bam-samps   <string>        " << "\t" << "Comma separated list of sample IDs for --bam" << "\n"
	   << "\t" << "--samp-sex    <string>        " << "\t" << "Comma separated list of sample sex for each sample ID (--bam-samps must be provided)" << "\n"
	   << "\t" << "--str-info    <string>        " << "\t" << "Tab file with additional per-STR info (see docs)" << "\n"
	   << "\t" << "--period      <string>        " << "\t" << "Only genotype loci with periods (motif lengths) in this comma-separated list." << "\n"
	   << "\t" << "--skip-qscore                 " << "\t" << "Skip calculation of Q-score" << "\n"
	   << "\n Options for different sequencing settings\n"
	   << "\t" << "--readlength  <int>           " << "\t" << "Read length. Default: " << options.read_len << "\n"
	   << "\t" << "--coverage    <float>         " << "\t" << "Average coverage. must be set for exome/targeted data. Comma separated list to specify for each BAM" << "\n"
	   << "\t" << "--model-gc-coverage           " << "\t" << "Model coverage as a function of GC content. Requires genome-wide data" << "\n"
	   << "\t" << "--insertmean  <float>         " << "\t" << "Fragment length mean. Comma separated list to specify for each BAM separately." << "\n" 
	   << "\t" << "--insertsdev  <float>         " << "\t" << "Fragment length standard deviation. Comma separated list to specify for each BAM separately. " << "\n"
	   << "\t" << "--nonuniform                  " << "\t" << "Indicate whether data has non-uniform coverage (i.e., exome)" << "\n"
	   << "\t" << "--min-sample-reads <int>      " << "\t" << "Minimum number of reads per sample." << "\n"
	   << "\n Advanced paramters for likelihood model:\n"
	   << "\t" << "--frrweight   <float>         " << "\t" << "Weight for FRR reads. Default: " << options.frr_weight << "\n"
	   << "\t" << "--enclweight  <float>         " << "\t" << "Weight for enclosing reads. Default: " << options.enclosing_weight << "\n"
	   << "\t" << "--spanweight  <float>         " << "\t" << "Weight for spanning reads. Default: " << options.spanning_weight << "\n"
	   << "\t" << "--flankweight <float>         " << "\t" << "Weight for flanking reads. Default: " << options.flanking_weight << "\n"
	   << "\t" << "--ploidy      <int>           " << "\t" << "Indicate whether data is haploid (1) or diploid (2). Default: " << options.ploidy << "\n"
	   << "\t" << "--skipofftarget               " << "\t" << "Skip off target regions included in the BED file." << "\n"
	   << "\t" << "--read-prob-mode              " << "\t" << "Use only read probability (ignore class probability)" << "\n"
	   << "\t" << "--numbstrap   <int>           " << "\t" << "Number of bootstrap samples. Default: " << options.num_boot_samp << "\n"
	   << "\t" << "--bootstrap-alpha <float>     " << "\t" << "Bootstrap CI two-sided tail probability. Default: " << options.bootstrap_alpha << "\n"
	   << "\t" << "--grid-threshold <int>        " << "\t" << "Use optimization rather than grid search to find MLE if more than this many possible alleles. Default: " << options.grid_threshold << "\n"
	   << "\t" << "--rescue-count <int>          " << "\t" << "Number of regions that GangSTR attempts to rescue mates from (excluding off-target regions) Default: " << options.rescue_count << "\n"
	   << "\t" << "--rescue-match-perc <float>   " << "\t" << "Minimum match fraction to accept rescued reads. Default: " << options.rescue_match_perc << "\n"
	   << "\t" << "--max-proc-read <int>         " << "\t" << "Maximum number of processed reads per sample before a region is skipped. Default: " << options.max_processed_reads_per_sample << "\n"
	   << "\t" << "--opt-xtol-rel-2d <float>     " << "\t" << "Relative tolerance for 2D optimizer. Default: " << options.nlopt_xtol_rel_2d << "\n"
	   << "\t" << "--opt-xtol-rel-1d <float>     " << "\t" << "Relative tolerance for 1D optimizer. Default: " << options.nlopt_xtol_rel_1d << "\n"
	   << "\n Parameters for local realignment:\n"
	   << "\t" << "--minscore    <int>           " << "\t" << "Minimum alignment score (out of 100). Default: " << options.min_score << "\n"
	   << "\t" << "--minmatch    <int>           " << "\t" << "Minimum number of matching basepairs on each end of enclosing reads. Default: " << options.min_match<< "\n"
	   << "\t" << "--realign-match-perc <float>  " << "\t" << "Minimum match fraction for realignment classification. Default: " << options.realign_match_perc << "\n"
	   << "\t" << "--max-spanning-data <int>     " << "\t" << "Maximum spanning data value (0 disables). Default: " << options.max_spanning_data_value << "\n"
	   << "\n Default stutter model parameters:\n"
	   << "\t" << "--stutterup   <float>         " << "\t" << "Stutter insertion probability. Default: " << options.stutter_up << "\n"
	   << "\t" << "--stutterdown <float>         " << "\t" << "Stutter deletion probability. Default: " << options.stutter_down << "\n"
	   << "\t" << "--stutterprob <float>         " << "\t" << "Stutter step size parameter. Default: " << options.stutter_p << "\n"
	   << "\n Parameters for more detailed info about each locus:\n"
	   << "\t" << "--output-bootstraps           " << "\t" << "Output file with bootstrap samples" << "\n"
	   << "\t" << "--output-readinfo             " << "\t" << "Output read class info (for debugging)" << "\n"
	   << "\t" << "--include-ggl                 " << "\t" << "Output GGL (special GL field) in VCF" << "\n"
	   << "\n Additional optional paramters:\n"
	   << "\t" << "-h,--help                     " << "\t" << "display this help screen" << "\n"
	   << "\t" << "--seed                        " << "\t" << "Random number generator initial seed" << "\n"
	   << "\t" << "-v,--verbose                  " << "\t" << "Print out useful progress messages" << "\n"
	   << "\t" << "--very                        " << "\t" << "Print out more detailed progress messages for debugging" << "\n"
	   << "\t" << "--quiet                       " << "\t" << "Don't print anything" << "\n"
	   << "\t" << "--version                     " << "\t" << "Print out the version of this software.\n"
	   << "\n\nThis program takes in aligned reads in BAM format\n"
	   << "and outputs estimated genotypes at each TR in VCF format.\n\n";
  cerr << help_msg.str();
  exit(1);
}

void parse_commandline_options(int argc, char* argv[], Options* options) {
  enum LONG_OPTIONS {
    OPT_MAXPROCREAD,
    OPT_SKIPQ,
    OPT_MINREAD,
    OPT_PERIOD,
    OPT_GGL,
    OPT_GRIDTHRESH,
    OPT_RESCUE,
    OPT_BAMFILES,
    OPT_BAMSAMP,
    OPT_SAMPSEX,
    OPT_STRINFO,
    OPT_CHROM,
    OPT_REFFA,
    OPT_REGIONS,
    OPT_OUT,
    OPT_HELP,
    OPT_WFRR,
    OPT_WENCLOSE,
    OPT_WSPAN,
    OPT_WFLANK,
    OPT_TARGETED,
    OPT_PLOIDY,
    OPT_READLEN,
    OPT_COVERAGE,
    OPT_GCCOV,
    OPT_SKIPOFF,
    OPT_NONUNIF,
    OPT_INSMEAN,
    OPT_INSSDEV,
    OPT_MINSCORE,
    OPT_MINMATCH,
    OPT_STUTUP,
    OPT_STUTDW,
    OPT_STUTPR,
    OPT_NBSTRAP,
    OPT_BOOTSTRAP_ALPHA,
    OPT_RDPROB,
    OPT_OUTBS,
    OPT_OUTREADINFO,
    OPT_OPT_XTOL2D,
    OPT_OPT_XTOL1D,
    OPT_RESCUE_MATCH_PERC,
    OPT_REALIGN_MATCH_PERC,
    OPT_MAX_SPANNING_DATA,
    OPT_SEED,
    OPT_VERBOSE,
    OPT_VERYVERBOSE,
    OPT_QUIET,
    OPT_VERSION,
  };
  static struct option long_options[] = {
    {"max-proc-read", required_argument, NULL, OPT_MAXPROCREAD},
    {"skip-qscore", no_argument, NULL, OPT_SKIPQ},
    {"min-sample-reads", required_argument, NULL, OPT_MINREAD},
    {"period",      required_argument,  NULL, OPT_PERIOD},
    {"include-ggl", no_argument, NULL, OPT_GGL},
    {"grid-threshold", required_argument, NULL, OPT_GRIDTHRESH},
    {"rescue-count", required_argument, NULL, OPT_RESCUE},
    {"rescue-match-perc", required_argument, NULL, OPT_RESCUE_MATCH_PERC},
    {"bam",         required_argument,  NULL, OPT_BAMFILES},
    {"bam-samps",   required_argument,  NULL, OPT_BAMSAMP},
    {"samp-sex",    required_argument,  NULL, OPT_SAMPSEX},
    {"str-info",    required_argument,  NULL, OPT_STRINFO},
    {"chrom",       required_argument,  NULL, OPT_CHROM},
    {"ref",         required_argument,  NULL, OPT_REFFA},
    {"regions",     required_argument,  NULL, OPT_REGIONS},
    {"out",         required_argument,  NULL, OPT_OUT},
    {"help",        no_argument,        NULL, OPT_HELP},
    {"frrweight",   required_argument,  NULL, OPT_WFRR},
    {"enclweight",  required_argument,  NULL, OPT_WENCLOSE},
    {"spanweight",  required_argument,  NULL, OPT_WSPAN},
    {"flankweight", required_argument,  NULL, OPT_WFLANK},
    {"targeted",  no_argument, NULL, OPT_TARGETED},
    {"ploidy",      required_argument,  NULL, OPT_PLOIDY},
    {"readlength",  required_argument,  NULL, OPT_READLEN},
    {"coverage",    required_argument,  NULL, OPT_COVERAGE},
    {"model-gc-coverage", no_argument, NULL, OPT_GCCOV},
    {"nonuniform",  no_argument,  NULL, OPT_NONUNIF},
    {"skipofftarget",no_argument,  NULL, OPT_SKIPOFF},
    {"insertmean",  required_argument,  NULL, OPT_INSMEAN},
    {"insertsdev",  required_argument,  NULL, OPT_INSSDEV},
    {"minscore",    required_argument,  NULL, OPT_MINSCORE},
    {"minmatch",    required_argument,  NULL, OPT_MINMATCH},
    {"stutterup",   required_argument,  NULL, OPT_STUTUP},
    {"stutterdown", required_argument,  NULL, OPT_STUTDW},
    {"stutterprob", required_argument,  NULL, OPT_STUTPR},
    {"numbstrap",   required_argument,  NULL, OPT_NBSTRAP},
    {"bootstrap-alpha", required_argument, NULL, OPT_BOOTSTRAP_ALPHA},
    {"read-prob-mode",   no_argument,  NULL, OPT_RDPROB},
    {"opt-xtol-rel-2d", required_argument, NULL, OPT_OPT_XTOL2D},
    {"opt-xtol-rel-1d", required_argument, NULL, OPT_OPT_XTOL1D},
    {"output-bootstraps", no_argument,      NULL, OPT_OUTBS},
    {"output-readinfo", no_argument,        NULL, OPT_OUTREADINFO},
    {"realign-match-perc", required_argument, NULL, OPT_REALIGN_MATCH_PERC},
    {"max-spanning-data", required_argument, NULL, OPT_MAX_SPANNING_DATA},
    {"seed",        required_argument,  NULL, OPT_SEED},
    {"verbose",     no_argument,        NULL, OPT_VERBOSE},
    {"very",  no_argument, NULL, OPT_VERYVERBOSE},
    {"quiet", no_argument, NULL, OPT_QUIET},
    {"version",     no_argument,        NULL, OPT_VERSION},
    {NULL,          no_argument,        NULL, 0},
  };
  std::vector<std::string> dist_means_str, dist_sdev_str, coverage_str, pers;
  int ch;
  int option_index = 0;
  ch = getopt_long(argc, argv, "hv?",
                   long_options, &option_index);
  while (ch != -1) {
    switch (ch) {
    case OPT_MAXPROCREAD:
      options->max_processed_reads_per_sample = atoi(optarg);
      break;
    case OPT_SKIPQ:
      options->skip_qscore = true;
      break;
    case OPT_MINREAD:
      options->min_reads_per_sample = atoi(optarg);
      break;
    case OPT_PERIOD:
      split_by_delim(optarg, ',', pers);
      options->period.clear();
      for (size_t i=0; i<pers.size(); i++) {
	options->period.push_back(atoi(pers[i].c_str()));
      }
      break;
    case OPT_GGL:
      options->include_ggl = true;
      break;
    case OPT_GRIDTHRESH:
      options->grid_threshold = atoi(optarg);
      break;
    case OPT_RESCUE:
      options->rescue_count = atoi(optarg);
      break;
    case OPT_RESCUE_MATCH_PERC:
      options->rescue_match_perc = atof(optarg);
      break;
    case OPT_BAMFILES:
      options->bamfiles.clear();
      split_by_delim(optarg, ',', options->bamfiles);
      break;
    case OPT_BAMSAMP:
      options->rg_sample_string = optarg;
      break;
    case OPT_SAMPSEX:
      options->sample_sex_string = optarg;
      break;
    case OPT_STRINFO:
      options->str_info_file = optarg;
      break;
    case OPT_CHROM:
      options->chrom = optarg;
      break;
    case OPT_REFFA:
      options->reffa = optarg;
      break;
    case OPT_REGIONS:
      options->regionsfile = optarg;
      break;
    case OPT_OUT:
      options->outprefix = optarg;
      break;
    case OPT_HELP:
    case 'h':
      show_help();
    case OPT_WFRR:
      options->frr_weight = atof(optarg);
      break;
    case OPT_WENCLOSE:
      options->enclosing_weight = atof(optarg);
      break;
    case OPT_WSPAN:
      options->spanning_weight = atof(optarg);
      break;
    case OPT_WFLANK:
      options->flanking_weight = atof(optarg);
      break;
    case OPT_TARGETED:
      options->genome_wide = false;
      break;
    case OPT_PLOIDY:
      options->ploidy = atoi(optarg);
      break;
    case OPT_READLEN:
      options->read_len = atoi(optarg);
      break;
    case OPT_COVERAGE:
      split_by_delim(optarg, ',', coverage_str);
      options->coverage.clear();
      for (size_t i=0; i<coverage_str.size(); i++) {
	options->coverage.push_back(strtof(coverage_str[i].c_str(), NULL));
      }
      break;
    case OPT_GCCOV:
      options->model_gc_cov = true;
      break;
    case OPT_NONUNIF:
      options->use_cov = false;
      options->use_off = false;
      break;
    case OPT_SKIPOFF:
      options->use_off = false;
      break;
    case OPT_INSMEAN:
      split_by_delim(optarg, ',', dist_means_str);
      options->dist_mean.clear();
      for (size_t i=0; i<dist_means_str.size(); i++) {
	options->dist_mean.push_back(strtof(dist_means_str[i].c_str(), NULL));
      }
      options->dist_man_set = true;
      break;
    case OPT_INSSDEV:
      split_by_delim(optarg, ',', dist_sdev_str);
      options->dist_sdev.clear();
      for (size_t i=0; i<dist_sdev_str.size(); i++) {
	options->dist_sdev.push_back(strtof(dist_sdev_str[i].c_str(), NULL));
      }
      options->dist_man_set = true;
      break;
    case OPT_MINSCORE:
      options->min_score = atoi(optarg);
      break;
    case OPT_MINMATCH:
      options->min_match = atoi(optarg);
      break;
    case OPT_STUTUP:
      options->stutter_up = atof(optarg);
      break;
    case OPT_STUTDW:
      options->stutter_down = atof(optarg);
      break;
    case OPT_STUTPR:
      options->stutter_p = atof(optarg);
      break;
    case OPT_NBSTRAP:
      options->num_boot_samp = atoi(optarg);
      break;
    case OPT_BOOTSTRAP_ALPHA:
      options->bootstrap_alpha = atof(optarg);
      break;
    case OPT_OUTBS:
      options->output_bootstrap = true;
      break;
    case OPT_RDPROB:
      options->read_prob_mode = true;
      break;
    case OPT_OPT_XTOL2D:
      options->nlopt_xtol_rel_2d = atof(optarg);
      break;
    case OPT_OPT_XTOL1D:
      options->nlopt_xtol_rel_1d = atof(optarg);
      break;
    case OPT_OUTREADINFO:
      options->output_readinfo= true;
      break;
    case OPT_REALIGN_MATCH_PERC:
      options->realign_match_perc = atof(optarg);
      break;
    case OPT_MAX_SPANNING_DATA:
      options->max_spanning_data_value = atoi(optarg);
      break;
    case OPT_SEED:
      options->seed = atoi(optarg);
      break;
    case OPT_VERBOSE:
    case 'v':
      options->verbose = true;
      break;
    case OPT_VERYVERBOSE:
      options->very_verbose = true;
      break;
    case OPT_QUIET:
      options->quiet = true;
      break;
    case OPT_VERSION:
      cerr << GangSTR_VER << endl;
      exit(0);
    case '?':
      show_help();
    default:
      show_help();
    };
    ch = getopt_long(argc, argv, "hv?",
		     long_options, &option_index);
  };
  // Leftover arguments are errors
  if (optind < argc) {
    PrintMessageDieOnError("Unnecessary leftover arguments", M_ERROR, false);
  }
  // Perform other error checking
  if (options->bamfiles.empty()) {
    PrintMessageDieOnError("No --bam files specified", M_ERROR, false);
  }
  if (options->regionsfile.empty()) {
    PrintMessageDieOnError("No --regions option specified", M_ERROR, false);
  }
  if (options->reffa.empty()) {
    PrintMessageDieOnError("No --ref option specified", M_ERROR, false);
  }
  if (options->outprefix.empty()) {
    PrintMessageDieOnError("No --out option specified", M_ERROR, false);
  }
  if (options->min_match < 0 or (options->read_len != -1 and options->min_match > options->read_len)){
    PrintMessageDieOnError("--minmatch parameter must be >= 0 and <= read_len (when --readlength is set)", M_ERROR, false);
  }
  if (options->min_score < 0 || options->min_score > 100){
    PrintMessageDieOnError("--minscore parameter must be in [0, 100] range", M_ERROR, false);
  }
  if (options->bootstrap_alpha <= 0.0 || options->bootstrap_alpha >= 1.0){
    PrintMessageDieOnError("--bootstrap-alpha parameter must be in (0, 1) range", M_ERROR, false);
  }
  if (options->nlopt_xtol_rel_2d <= 0.0){
    PrintMessageDieOnError("--opt-xtol-rel-2d parameter must be > 0", M_ERROR, false);
  }
  if (options->nlopt_xtol_rel_1d <= 0.0){
    PrintMessageDieOnError("--opt-xtol-rel-1d parameter must be > 0", M_ERROR, false);
  }
  if (options->rescue_match_perc <= 0.0 || options->rescue_match_perc > 1.0){
    PrintMessageDieOnError("--rescue-match-perc parameter must be in (0, 1] range", M_ERROR, false);
  }
  if (options->realign_match_perc <= 0.0 || options->realign_match_perc > 1.0){
    PrintMessageDieOnError("--realign-match-perc parameter must be in (0, 1] range", M_ERROR, false);
  }
  if (options->max_spanning_data_value < 0){
    PrintMessageDieOnError("--max-spanning-data parameter must be >= 0", M_ERROR, false);
  }
  
  // Check if bam-samps provided when user uses sam-sex
  if (options->rg_sample_string.empty() and !options->sample_sex_string.empty()){
    PrintMessageDieOnError("--sample-sex requires input argument --bam-samps to be set", M_ERROR, false);
  }
}

int main(int argc, char* argv[]) {
  // Set up
  Options options;
  parse_commandline_options(argc, argv, &options);
  stringstream full_command_ss;
  full_command_ss << "GangSTR-" << GangSTR_VER;
  for (int i = 1; i < argc; i++) {
    full_command_ss << " " << argv[i];
  }
  std::string full_command = full_command_ss.str();
  RegionReader region_reader(options.regionsfile);
  Locus locus;
  int merge_type = BamCramMultiReader::ORDER_ALNS_BY_FILE;
  BamCramMultiReader bamreader(options.bamfiles, options.reffa, merge_type);

  // Extract sample info
  SampleInfo sample_info;
  if (!options.rg_sample_string.empty()) {
    if (!sample_info.SetCustomReadGroups(options)) {
      PrintMessageDieOnError("Error setting custom read groups", M_ERROR, false);
    }
  } else {
    if (!sample_info.LoadReadGroups(options, bamreader)) {
      PrintMessageDieOnError("Error loading read groups", M_ERROR, false);
    }
  }
  // Extract information from bam file (read length, insert size distribution, ..)
  RefGenome refgenome(options.reffa);
  
  if (!sample_info.ExtractBamInfo(options, bamreader, region_reader, refgenome)) {
    PrintMessageDieOnError("Error extracting info from BAM file", M_ERROR, false);
  }
  options.realignment_flanklen = sample_info.GetReadLength();
  // Write out values found for each sample
  sample_info.PrintSampleInfo(options.outprefix + ".samplestats.tab");
  
  // Process each region
  region_reader.Reset();
  STRInfo str_info(options);
  VCFWriter vcfwriter(options.outprefix + ".vcf", full_command,
		      refgenome, 
		      sample_info, options.include_ggl);
  Genotyper genotyper(refgenome, options, sample_info, str_info);

  // First pass: Read all loci from the file
  std::vector<Locus*> all_loci;
  while (region_reader.GetNextRegion(&locus)) {
      if (!options.chrom.empty() && locus.chrom != options.chrom) {continue;}
      if (!options.period.empty() &&
    std::find(options.period.begin(), options.period.end(), locus.period) == options.period.end()) {
        continue;
      }
      Locus* new_locus = new Locus();
      *new_locus = locus;
      all_loci.push_back(new_locus);
  }

  // Second pass: Learn stutter models
  if (!genotyper.LearnStutterModels(&bamreader, all_loci)) {
    PrintMessageDieOnError("Error learning stutter models", M_ERROR, false);
  }

  // Third pass: Genotype each locus
  stringstream ss;
  for (auto const& current_locus : all_loci) {
    ss.str("");
    ss.clear();
    ss << "Processing " << current_locus->chrom << ":" << current_locus->start;
    PrintMessageDieOnError(ss.str(), M_PROGRESS, options.quiet);

    if (options.use_off){
      current_locus->offtarget_share = 1.0;
    }
    else{
      current_locus->offtarget_share = 0.0;
    }

    if (genotyper.ProcessLocus(&bamreader, current_locus)) {
      vcfwriter.WriteRecord(*current_locus);
    }
    delete current_locus; // free memory
  };
}
