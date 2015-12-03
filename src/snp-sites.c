/*
 *  Wellcome Trust Sanger Institute
 *  Copyright (C) 2013  Wellcome Trust Sanger Institute
 *  
 *  This program is free software; you can redistribute it and/or
 *  modify it under the terms of the GNU General Public License
 *  as published by the Free Software Foundation; either version 3
 *  of the License, or (at your option) any later version.
 *  
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *  
 *  You should have received a copy of the GNU General Public License
 *  along with this program; if not, write to the Free Software
 *  Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301, USA.
 */



#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <regex.h>
#include "vcf.h"
#include "alignment-file.h"
#include "snp-sites.h"
#include "phylib-of-snp-sites.h"
#include "fasta-of-snp-sites.h"

int generate_snp_sites(char filename[],int output_multi_fasta_file, int output_vcf_file, int output_phylip_file, char output_filename[])
{
	int i;
	
	detect_snps(filename);

	char* bases_for_snps[get_number_of_snps()];
	
	for(i = 0; i < get_number_of_snps(); i++)
	{
		bases_for_snps[i] = calloc(get_number_of_samples()+1 ,sizeof(char));
	}
	
	get_bases_for_each_snp(filename, bases_for_snps);
	
	char output_filename_base[FILENAME_MAX];
	char filename_without_directory[FILENAME_MAX];
	strip_directory_from_filename(filename, filename_without_directory);
	strncpy(output_filename_base, filename_without_directory, FILENAME_MAX);
	
	if(output_filename != NULL && *output_filename != '\0')
	{
		strncpy(output_filename_base, output_filename, FILENAME_MAX);
	}

	if(output_vcf_file)
	{
		char vcf_output_filename[FILENAME_MAX];
		strncpy(vcf_output_filename, output_filename_base, FILENAME_MAX);
		if((output_vcf_file + output_phylip_file + output_multi_fasta_file) > 1 || (output_filename == NULL || *output_filename == '\0') )
		{
			strcat(vcf_output_filename, ".vcf");
		}
		
	  create_vcf_file(vcf_output_filename, get_snp_locations(), get_number_of_snps(), bases_for_snps, get_sequence_names(), get_number_of_samples(), get_length_of_genome());
  }

  if(output_phylip_file)
  {
		char phylip_output_filename[FILENAME_MAX];
		strncpy(phylip_output_filename, output_filename_base, FILENAME_MAX);
		if((output_vcf_file + output_phylip_file + output_multi_fasta_file) > 1 || (output_filename == NULL || *output_filename == '\0') )
		{
			strcat(phylip_output_filename, ".phylip");
		}
	  create_phylib_of_snp_sites(phylip_output_filename, get_number_of_snps(), bases_for_snps, get_sequence_names(), get_number_of_samples());
  }

  if((output_multi_fasta_file) || (output_vcf_file ==0 && output_phylip_file == 0 && output_multi_fasta_file == 0))
  {
		char multi_fasta_output_filename[FILENAME_MAX];
		strncpy(multi_fasta_output_filename, output_filename_base, FILENAME_MAX);
		if((output_vcf_file + output_phylip_file + output_multi_fasta_file) > 1 || (output_filename == NULL || *output_filename == '\0') )
		{
			strcat(multi_fasta_output_filename, ".snp_sites.aln");
		}
	  create_fasta_of_snp_sites(multi_fasta_output_filename, get_number_of_snps(), bases_for_snps, get_sequence_names(), get_number_of_samples());
  }

  // free memory
	free(get_snp_locations());
	for(i = 0; i < get_number_of_samples(); i++)
	{
	//	free(get_sequence_names().[i]);
	}
	for(i = 0; i < get_number_of_snps(); i++)
	{
		free(bases_for_snps[i]);
	}
	

	return 1;
}

// Inefficient
void strip_directory_from_filename(char * input_filename, char * output_filename)
{
  int i;
  int end_index = 0;
  int last_forward_slash_index = -1;
  for(i = 0; i< FILENAME_MAX; i++)
  {
    if(input_filename[i] == '/')
    {
      last_forward_slash_index = i;
    }
    
    if(input_filename[i] == '\0' || input_filename[i] == '\n')
    {
      end_index = i;
      break;
    }
  }
  
  int current_index = 0;
  for(i = last_forward_slash_index+1; i< end_index; i++)
  {
    output_filename[current_index] = input_filename[i];
    current_index++;
  }
  output_filename[current_index] = '\0';
}


