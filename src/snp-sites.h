/*
 *  Wellcome Trust Sanger Institute
 *  Copyright (C) 2013-2016  Wellcome Trust Sanger Institute
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


#ifndef _SNP_SITES_H_
#define _SNP_SITES_H_

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <errno.h>

int generate_snp_sites(char filename[],
                       int output_multi_fasta_file,
                       int output_vcf_file,
                       int output_phylip_file,
                       char output_filename[]);

int generate_snp_sites_with_ref(char filename[],
                                int output_multi_fasta_file,
                                int output_vcf_file,
                                int output_phylip_file,
                                char output_filename[]);

int generate_snp_sites_with_ref_pure_mono(char filename[],
                                          int output_multi_fasta_file,
                                          int output_vcf_file,
                                          int output_phylip_file,
                                          char output_filename[],
                                          int output_reference,
                                          int pure_mode,
                                          int output_monomorphic);

void count_constant_sites(char multi_fasta_filename[], char filename[]);

void strip_directory_from_filename(char *input_filename,
                                   char *output_filename);

#endif
