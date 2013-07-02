/*
 *  Wellcome Trust Sanger Institute
 *  Copyright (C) 2011  Wellcome Trust Sanger Institute
 *  
 *  This program is free software; you can redistribute it and/or
 *  modify it under the terms of the GNU General Public License
 *  as published by the Free Software Foundation; either version 2
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
#include "parse_phylip.h"
#include "alignment_file.h"
#include "string_cat.h"

int num_samples;
int num_snps;
char ** sequences;
char ** phylip_sample_names;


void update_sequence_base(char new_sequence_base, int sequence_index, int base_index)
{
	sequences[sequence_index][base_index] = new_sequence_base;	
}

void get_sequence_for_sample_name(char * sequence_bases, char * sample_name)
{
	int sequence_index;
	sequence_index = find_sequence_index_from_sample_name( sample_name);
	if(sequence_index < 0)
	{
    printf("Couldnt find sequence name %s with index %d\n", sample_name,sequence_index);
	  exit(1);
  }

	memcpy(sequence_bases, sequences[sequence_index], size_of_string(sequences[sequence_index]) +1);
}



int does_column_contain_snps(int snp_column, char reference_base)
{
	int i;
	for(i = 0; i < num_samples; i++)
	{
		if(sequences[i][snp_column] == '\0' || sequences[i][snp_column] == '\n')
		{
			return 0;	
		}
		
		if(sequences[i][snp_column]  != '-' && toupper(sequences[i][snp_column])  != 'N' && sequences[i][snp_column] != reference_base)
		{
			return 1;
		}
	}
	return 0;
}

int number_of_samples_from_parse_phylip()
{
	return num_samples;
}

void get_sample_names_from_parse_phylip(char ** sample_names)
{
	int i;
	for(i = 0; i< num_samples; i++)
	{
		sample_names[i] =  (char *) calloc(MAX_SAMPLE_NAME_SIZE,sizeof(char));
		memcpy(sample_names[i], phylip_sample_names[i], size_of_string(phylip_sample_names[i]) +1);
	}
}

	

void filter_sequence_bases_and_rotate(char * reference_bases, char ** filtered_bases_for_snps, int number_of_filtered_snps)
{
	int i,j,reference_index;
	
	for(j = 0; j < number_of_filtered_snps; j++)
	{
		filtered_bases_for_snps[j] = (char *) calloc(num_samples,sizeof(char));
	}
		
	for(i = 0; i < num_samples; i++)
	{
		int filtered_base_counter = 0;
		
		for(reference_index = 0; reference_index < num_snps; reference_index++)
		{
			if(reference_bases[reference_index] == '\0')
			{
				break;	
			}
			
			if(reference_bases[reference_index] != '*' && sequences[i][reference_index] != '\0' && sequences[i][reference_index] != '\n')
			{
				filtered_bases_for_snps[filtered_base_counter][i] = sequences[i][reference_index];
				filtered_base_counter++;
			}
		}
	}
	for(j = 0; j < number_of_filtered_snps; j++)
	{
		 filtered_bases_for_snps[j][num_samples] = '\0';	
	}

}


int find_sequence_index_from_sample_name( char * sample_name)
{
	int i;
	
	for(i =0; i< num_samples; i++)
	{
		if(strcmp(sample_name,phylip_sample_names[i]) == 0)	
		{
			return i;
		}
	}
	return -1;
}


void load_sequences_from_phylib_file(char phylip_filename[])
{
	FILE * phylib_file_pointer;
	phylib_file_pointer=fopen(phylip_filename, "r");
	load_sequences_from_phylib(phylib_file_pointer);
}

void load_sequences_from_phylib(FILE * phylip_file_pointer)
{	
	rewind(phylip_file_pointer);
	char * line_buffer = (char *) calloc(MAX_READ_BUFFER,sizeof(char));
	int i;
	
	// The first line contains the number of samples and snps
	line_buffer[0] = '\0';
	line_buffer = read_line(line_buffer, phylip_file_pointer);
	
	num_samples = get_number_of_samples_from_phylip(line_buffer);
	num_snps = get_number_of_snps_from_phylip(line_buffer);
	
	sequences = (char **) calloc(num_samples,sizeof(char *));
	phylip_sample_names = (char **) calloc(num_samples,sizeof(char *));
	
	for(i = 0; i < num_samples; i++)
	{
		sequences[i] = (char *) calloc(num_snps,sizeof(char));
		phylip_sample_names[i] = (char *) calloc(MAX_SAMPLE_NAME_SIZE,sizeof(char));
	}
	
	int sample_counter = 0;

	do{
		line_buffer[0] = '\0';
		line_buffer = read_line(line_buffer, phylip_file_pointer);
		
		if(line_buffer[0] == '\0')
		{
			break;
		}
		
		int found_sequence = 0;
		int sequence_offset = 0;
		for(i = 0 ; i< (num_snps + MAX_SAMPLE_NAME_SIZE); i++)
		{
			if(line_buffer[i] == '\0' || line_buffer[i] == '\n')
			{
				sequences[sample_counter][i-sequence_offset] = '\0';
				break;	
			}
			else if(found_sequence == 1)
			{
				// Read in the sequence data of the sample
				sequences[sample_counter][i-sequence_offset] = line_buffer[i];
			}
			else
			{
				// Read in the name of the sample
				if(line_buffer[i] == '\t')
				{
					found_sequence = 1;
					sequence_offset = i+1;
					phylip_sample_names[sample_counter][i] = '\0';
				}
				else
				{
					phylip_sample_names[sample_counter][i] = line_buffer[i];
				}
			}
		}
		sample_counter++;

	}while(line_buffer[0] != '\0');
	
}

int get_number_of_samples_from_phylip(char * phylip_string)
{
	int i;
	char number_of_samples_string[20] = {0}; 

	for(i = 0; i< MAX_READ_BUFFER; i++)
	{
		if(phylip_string[i] == '\0' || phylip_string[i] == '\n' || phylip_string[i] == ' ')
		{
			break;	
		}
		else
		{
			number_of_samples_string[i] = phylip_string[i];
		}
	}
	
	return atoi(number_of_samples_string);
}

int get_number_of_snps_from_phylip(char * phylip_string)
{
	int found = 0;
	int i;
	int offset = 0;
	char number_of_snps_string[20] = {0}; 
	
	for(i = 0; i< MAX_READ_BUFFER; i++)
	{
		if(phylip_string[i] == '\0' || phylip_string[i] == '\n' )
		{
			break;
		}
		
		if(phylip_string[i] == ' ')
		{
			found = 1;
			offset = i+1;
		}
		else
		{
			if(found == 1)
			{
				number_of_snps_string[i-offset] = phylip_string[i];
			}
		}
	}
	
	return atoi(number_of_snps_string);
}