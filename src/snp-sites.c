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
#include "parse-phylip.h"
#include "string-cat.h"


void build_snp_locations(int snp_locations[], char reference_sequence[])
{
	int i;
	int snp_counter = 0;
	
	for(i = 0; reference_sequence[i]; i++)
    {
		if(reference_sequence[i] == '*')
		{
			snp_locations[snp_counter] = i;
			snp_counter++;
		}
	}
}


int generate_snp_sites(char filename[],int output_multi_fasta_file, int output_vcf_file, int output_phylip_file, char output_filename[])
{
	size_t length_of_genome;
	char * reference_sequence;
	int number_of_snps;
	int * snp_locations;
	int number_of_samples;
	int i;
	
	length_of_genome = genome_length(filename);
	reference_sequence = (char *) calloc((length_of_genome +1),sizeof(char));
	
	build_reference_sequence(reference_sequence,filename);
	number_of_snps = detect_snps(reference_sequence, filename, length_of_genome);
	
	snp_locations = (int *) calloc((number_of_snps+1),sizeof(int));
	build_snp_locations(snp_locations, reference_sequence);
	free(reference_sequence);
	
	number_of_samples = number_of_sequences_in_file(filename);
	
	// Find out the names of the sequences
	char* sequence_names[number_of_samples];
	sequence_names[number_of_samples-1] = '\0';
	for(i = 0; i < number_of_samples; i++)
	{
		sequence_names[i] = calloc(MAX_SAMPLE_NAME_SIZE,sizeof(char));
	}
	
	get_sample_names_for_header(filename, sequence_names, number_of_samples);
	
	char* bases_for_snps[number_of_snps];
	
	for(i = 0; i < number_of_snps; i++)
	{
		bases_for_snps[i] = calloc(number_of_samples+1 ,sizeof(char));
	}
	
	get_bases_for_each_snp(filename, snp_locations, bases_for_snps, length_of_genome, number_of_snps);
	
	char output_filename_base[MAX_FILENAME_SIZE];
	char filename_without_directory[MAX_FILENAME_SIZE];
	strip_directory_from_filename(filename, filename_without_directory);
	memcpy(output_filename_base,filename_without_directory, size_of_string(filename_without_directory)+1 );
	
	if(output_filename != NULL && *output_filename != '\0')
	{
		memcpy(output_filename_base,output_filename, size_of_string(output_filename)+1 );
	}

	if(output_vcf_file)
	{
		char * vcf_output_filename;
		vcf_output_filename = calloc(MAX_FILENAME_SIZE,sizeof(char));
		memcpy(vcf_output_filename, output_filename_base, (MAX_FILENAME_SIZE)*sizeof(char));
		if((output_vcf_file + output_phylip_file + output_multi_fasta_file) > 1 || (output_filename == NULL || *output_filename == '\0') )
		{
			char extension[5] = {".vcf"};
			concat_strings_created_with_malloc(vcf_output_filename,extension);
		}
		
	  create_vcf_file(vcf_output_filename, snp_locations, number_of_snps, bases_for_snps, sequence_names, number_of_samples, length_of_genome);
		free(vcf_output_filename);
  }

  if(output_phylip_file)
  {
		char *phylip_output_filename;
		phylip_output_filename = calloc(MAX_FILENAME_SIZE,sizeof(char));
		memcpy(phylip_output_filename, output_filename_base, (MAX_FILENAME_SIZE)*sizeof(char));
		if((output_vcf_file + output_phylip_file + output_multi_fasta_file) > 1 || (output_filename == NULL || *output_filename == '\0') )
		{
			char extension[10] = {".phylip"};
			concat_strings_created_with_malloc(phylip_output_filename,extension);
		}
	  create_phylib_of_snp_sites(phylip_output_filename, number_of_snps, bases_for_snps, sequence_names, number_of_samples);
		free(phylip_output_filename);
  }

  if((output_multi_fasta_file) || (output_vcf_file ==0 && output_phylip_file == 0 && output_multi_fasta_file == 0))
  {
		char *multi_fasta_output_filename;
		multi_fasta_output_filename = calloc(MAX_FILENAME_SIZE,sizeof(char));
		memcpy(multi_fasta_output_filename, output_filename_base, (MAX_FILENAME_SIZE)*sizeof(char));
		if((output_vcf_file + output_phylip_file + output_multi_fasta_file) > 1 || (output_filename == NULL || *output_filename == '\0') )
		{
			char extension[20] = {".snp_sites.aln"};
			concat_strings_created_with_malloc(multi_fasta_output_filename,extension);
		}
	  create_fasta_of_snp_sites(multi_fasta_output_filename, number_of_snps, bases_for_snps, sequence_names, number_of_samples);
	  free(multi_fasta_output_filename);
  }

  // free memory
	free(snp_locations);
	for(i = 0; i < number_of_samples; i++)
	{
		free(sequence_names[i]);
	}
	for(i = 0; i < number_of_snps; i++)
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
  for(i = 0; i< MAX_FILENAME_SIZE; i++)
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

// return new number of snps
int refilter_existing_snps(char * reference_bases, int number_of_snps, char ** column_names, int number_of_columns,int * snp_locations, int * filtered_snp_locations)
{
	// go through each snp column and check to see if there is still variation
	int i;
	int number_of_filtered_snps = number_of_snps;
	for(i = 0; i < number_of_snps; i++)
	{
		if( does_column_contain_snps(i, reference_bases[i]) == 0)
		{
			snp_locations[i] = -1;
			reference_bases[i] = '*';
			
			number_of_filtered_snps--;
		}
	}
	
	remove_filtered_snp_locations(filtered_snp_locations, snp_locations, number_of_snps);
	return number_of_filtered_snps;
}

void remove_filtered_snp_locations(int * filtered_snp_locations, int * snp_locations, int number_of_snps)
{
	int i;
	int filtered_counter=0;
	for(i = 0; i< number_of_snps; i++)
	{
		if(snp_locations[i] != -1)
		{
			filtered_snp_locations[filtered_counter] = snp_locations[i];
			filtered_counter++;
		}
	}
}













