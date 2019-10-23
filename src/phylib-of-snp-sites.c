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
#include "phylib-of-snp-sites.h"


void create_phylib_of_snp_sites(char filename[], int number_of_snps, char **bases_for_snps, char **sequence_names,
                                int number_of_samples, int output_reference, char *pseudo_reference_sequence,
                                int snp_locations[])
{
    FILE *phylip_file_pointer;
    int sample_counter;
    int snp_counter;

    phylip_file_pointer = fopen(filename, "w");

    if (output_reference == 1) {
        fprintf(phylip_file_pointer, "%d %d\n", number_of_samples + 1, number_of_snps);
        fprintf(phylip_file_pointer, "pseudo_reference_sequence\t");
        for (snp_counter = 0; snp_counter < number_of_snps; snp_counter++) {
            fputc(pseudo_reference_sequence[snp_locations[snp_counter]], phylip_file_pointer);
        }
        fprintf(phylip_file_pointer, "\n");
    } else {
        fprintf(phylip_file_pointer, "%d %d\n", number_of_samples, number_of_snps);
    }

    for (sample_counter = 0; sample_counter < number_of_samples; sample_counter++) {
        // sequence_name can be more than 10 (relaxed phylib format) and contain [\w\s]
        //TODO check for illegal characters [^\w\s]
        fprintf(phylip_file_pointer, "%s\t", sequence_names[sample_counter]);

        for (snp_counter = 0; snp_counter < number_of_snps; snp_counter++) {
            fputc(bases_for_snps[snp_counter][sample_counter], phylip_file_pointer);
        }
        fprintf(phylip_file_pointer, "\n");
    }

    fclose(phylip_file_pointer);
}
