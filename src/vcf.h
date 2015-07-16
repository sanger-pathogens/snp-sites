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


#ifndef _VCF_H_
#define _VCF_H_

void output_vcf_header( FILE * vcf_file_pointer, char ** sequence_names, int number_of_samples);
void create_vcf_file(char filename[], int snp_locations[], int number_of_snps, char ** bases_for_snps, char ** sequence_names, int number_of_samples);
void output_vcf_snps(FILE * vcf_file_pointer, char ** bases_for_snps, int * snp_locations, int number_of_snps, int number_of_samples);
void output_vcf_row(FILE * vcf_file_pointer, char * bases_for_snp, int snp_location, int number_of_samples);
void output_vcf_row_samples_bases(FILE * vcf_file_pointer, char reference_base, char * alt_bases, char * bases_for_snp, int number_of_samples);
void alternative_bases(char reference_base, char * bases_for_snp, char alt_bases[], int number_of_samples);
char * format_alternative_bases(char *);
char * format_allele_index(char, char, char *);
int check_if_char_in_string(char search_string[], char target_char, int search_string_length);
#define MAX_FILENAME_SIZE 250
#define MAXIMUM_NUMBER_OF_ALT_BASES 30

#endif
