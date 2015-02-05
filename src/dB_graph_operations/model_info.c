/*
  model_info.c
*/

#include "model_info.h"

void initialise_model_info(GraphAndModelInfo* model_info, GraphInfo* ginfo, 
			   long long genome_len, double mu, //double seq_err_rate_per_base,
			   int ref_colour, int num_chroms, ExperimentType type, AssumptionsOnGraphCleaning assump)
{
  model_info->ginfo = ginfo;
  model_info->genome_len = genome_len;
  model_info->mu = mu;
  // model_info->seq_error_rate_per_base = seq_err_rate_per_base;
  model_info->ref_colour=ref_colour; //if -1, that means no ref colour
  model_info->num_haploid_chromosomes = num_chroms;
  model_info->expt_type=type;
  model_info->assump = assump;
}
