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

#ifndef _ALIGNMENT_FILE_H_
#define _ALIGNMENT_FILE_H_

#include "kseq.h"

int detect_snps(char reference_sequence[],  char filename[], int length_of_genome);
int line_length(FILE * alignment_file_pointer);
int build_reference_sequence(char reference_sequence[], char filename[]);
void advance_to_sequence(FILE * alignment_file_pointer);
void advance_to_sequence_name(FILE * alignment_file_pointer);
int genome_length(char filename[]);
int read_line(char sequence[], FILE * pFilePtr);
int number_of_sequences_in_file(char filename[]);
void get_sample_names_for_header(char filename[], char ** sequence_names, int number_of_samples);
char filter_invalid_characters(char input_char);
void get_bases_for_each_snp(char filename[], int snp_locations[], char ** bases_for_snps, int length_of_genome, int number_of_snps);



#define MAX_READ_BUFFER 1048576
#define MAX_SAMPLE_NAME_SIZE 1024

#endif