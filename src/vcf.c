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
#include <assert.h>

void create_vcf_file(char filename[], int snp_locations[],int number_of_snps, char ** bases_for_snps, char ** sequence_names, int number_of_samples, size_t length_of_genome)
{
	FILE *vcf_file_pointer;
	char * base_filename;
	base_filename = (char *) malloc(FILENAME_MAX*sizeof(char));
	strcpy(base_filename, filename);
	
	vcf_file_pointer=fopen(base_filename, "w");
	output_vcf_header(vcf_file_pointer,sequence_names, number_of_samples, length_of_genome);
	output_vcf_snps(vcf_file_pointer, bases_for_snps, snp_locations, number_of_snps, number_of_samples);
  fclose(vcf_file_pointer);
	free(base_filename);
}

void output_vcf_snps(FILE * vcf_file_pointer, char ** bases_for_snps, int * snp_locations, int number_of_snps, int number_of_samples)
{
	int i;
	for(i=0; i < number_of_snps; i++)
	{
		output_vcf_row(vcf_file_pointer, bases_for_snps[i], snp_locations[i], number_of_samples);
	}
}

void output_vcf_header( FILE * vcf_file_pointer, char ** sequence_names, int number_of_samples, size_t length_of_genome)
{
	int i;
	fprintf( vcf_file_pointer, "##fileformat=VCFv4.1\n" );	
	fprintf( vcf_file_pointer, "##contig=<ID=1,length=%i>\n", length_of_genome );	
	fprintf( vcf_file_pointer, "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n" );	
	fprintf( vcf_file_pointer, "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT" );
	
	for(i=0; i<number_of_samples; i++)
	{
		fprintf( vcf_file_pointer, "\t%s",  sequence_names[i]);
	}
	fprintf( vcf_file_pointer, "\n");
}

void output_vcf_row(FILE * vcf_file_pointer, char * bases_for_snp, int snp_location, int number_of_samples)
{
	char reference_base =  bases_for_snp[0];
	if(reference_base == '\0')
	{
		return;	
	}
	
	// Chromosome
	fprintf( vcf_file_pointer, "1\t");
	
	// Position
	fprintf( vcf_file_pointer, "%d\t", (int) snp_location +1 );	
	
	//ID
	fprintf( vcf_file_pointer, ".\t");
	
	// REF
	fprintf( vcf_file_pointer, "%c\t", (char) reference_base );
	
	// ALT
	// Need to look through list and find unique characters
	
	char * alt_bases = alternative_bases(reference_base, bases_for_snp, number_of_samples);
	char * alternative_bases_string = format_alternative_bases(alt_bases);
	fprintf( vcf_file_pointer, "%s\t", alternative_bases_string );
	free(alternative_bases_string);
	
	// QUAL
	fprintf( vcf_file_pointer, ".\t");
	
	// FILTER
	fprintf( vcf_file_pointer, ".\t");
	
	// INFO
	fprintf( vcf_file_pointer, ".\t");
	
	// FORMAT
	fprintf( vcf_file_pointer, "GT\t");
	
	// Bases for each sample
	output_vcf_row_samples_bases(vcf_file_pointer, reference_base, alt_bases, bases_for_snp, number_of_samples );
	free(alt_bases);
	
	fprintf( vcf_file_pointer, "\n");	
}


char * alternative_bases(char reference_base, char * bases_for_snp, int number_of_samples)
{
	int i;
	int num_alt_bases = 0;
	char * alt_bases = calloc(MAXIMUM_NUMBER_OF_ALT_BASES+1, sizeof(char));
	for(i=0; i< number_of_samples; i++ )
	{
		if((bases_for_snp[i] != reference_base) && (bases_for_snp[i] != '-') && (toupper(bases_for_snp[i]) != 'N') )
		{
			if(check_if_char_in_string(alt_bases, bases_for_snp[i], num_alt_bases) == 0)
			{
				if (num_alt_bases >= MAXIMUM_NUMBER_OF_ALT_BASES)
				{
					fprintf(stderr, "Unexpectedly large number of alternative bases found between sequences.  Please check input file is not corrupted\n\n");
					fflush(stderr);
					exit(EXIT_FAILURE);
				}
				alt_bases[num_alt_bases] = bases_for_snp[i];
				num_alt_bases++;
			}
		}
	}
	return alt_bases;
}

char * format_allele_index(char base, char reference_base, char * alt_bases)
{
	int length_of_alt_bases = strlen(alt_bases);
	assert(length_of_alt_bases < 100);
	char * result = calloc(3, sizeof(char));
	int index;
	if (reference_base == base || toupper(base) == 'N' || base == '-')
	{
		sprintf(result, "0");
	}
	else
	{
		sprintf(result, ".");
		for (index = 1; index <= length_of_alt_bases; index++)
		{
			if (alt_bases[index-1] == base)
			{
				sprintf(result, "%i", index);
				break;
			}
		}
	}
	return result;
}

char * format_alternative_bases(char * alt_bases)
{
	int number_of_alt_bases = strlen(alt_bases);
	assert( number_of_alt_bases < MAXIMUM_NUMBER_OF_ALT_BASES );
	char * formatted_alt_bases = calloc(number_of_alt_bases*2 + 1, sizeof(char));
	int i;
	formatted_alt_bases[0] = alt_bases[0];
	for (i = 1; i < number_of_alt_bases; i++)
	{
		formatted_alt_bases[i*2 - 1] = ',';
		formatted_alt_bases[i*2] = alt_bases[i];
	}
	return formatted_alt_bases;
}

int check_if_char_in_string(char search_string[], char target_char, int search_string_length)
{
	int i;
	for(i=0; i < search_string_length ; i++ )
	{
		if(search_string[i] == target_char)
		{
			return 1;
		}
	}
	return 0;
}

void output_vcf_row_samples_bases(FILE * vcf_file_pointer, char reference_base, char * alt_bases, char * bases_for_snp, int number_of_samples)
{
	int i;
	char * format;
	
	for(i=0; i < number_of_samples ; i++ )
	{
		format = format_allele_index(bases_for_snp[i], reference_base, alt_bases);
		fprintf( vcf_file_pointer, "%s", format);	
		free(format);
		if(i+1 != number_of_samples)
		{
			fprintf( vcf_file_pointer, "\t");
		}
	}
}




