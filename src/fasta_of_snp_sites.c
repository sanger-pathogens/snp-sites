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
#include <regex.h>
#include "fasta_of_snp_sites.h"

void create_fasta_of_snp_sites(char filename[], int number_of_snps, char ** bases_for_snps, char ** sequence_names, int number_of_samples)
{
	FILE *fasta_file_pointer;
	int sample_counter;
	int snp_counter; 
	char * base_filename;
	
	base_filename = (char *) calloc(MAX_FILENAME_SIZE,sizeof(char));
	memcpy(base_filename, filename, MAX_FILENAME_SIZE*sizeof(char));
	fasta_file_pointer = fopen(base_filename, "w");
	
	for(sample_counter=0; sample_counter< number_of_samples; sample_counter++)
	{
		fprintf( fasta_file_pointer, ">%s\n", sequence_names[sample_counter]);
		for(snp_counter=0; snp_counter< number_of_snps; snp_counter++)
		{
			fprintf( fasta_file_pointer, "%c", bases_for_snps[snp_counter][sample_counter]);
		}
		fprintf( fasta_file_pointer, "\n");
	}
  fclose(fasta_file_pointer);
}