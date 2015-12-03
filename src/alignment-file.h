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

#ifndef _ALIGNMENT_FILE_H_
#define _ALIGNMENT_FILE_H_

#include "kseq.h"

void detect_snps( char filename[]);
void get_bases_for_each_snp(char filename[], int snp_locations[], char ** bases_for_snps, size_t length_of_genome, int number_of_snps);
int is_unknown(char base);
int get_length_of_genome();
int get_number_of_samples();
int get_number_of_snps();
char ** get_sequence_names();
int * get_snp_locations();

#define MAX_SAMPLE_NAME_SIZE 2048
#define DEFAULT_NUM_SAMPLES 65536
#endif
