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
#include <ctype.h>
#include <unistd.h>
#include <getopt.h>
#include "snp-sites.h"
#include "config.h"

#define PROGRAM_NAME "snp-sites"
#define PROGRAM_VERSION PACKAGE_VERSION

static void print_usage()
{
	printf("Usage: snp-sites [-mvph] [-o output_filename] <file>\n");
	printf("This program finds snp sites from a multi FASTA alignment file.\n");
	printf(" -r		output internal pseudo reference sequence\n");
  printf(" -m		output a multi fasta alignment file (default)\n");
	printf(" -v		output a VCF file\n");
	printf(" -p		output a phylip file\n");
	printf(" -o		specify an output filename\n");
	printf(" -h		this help message\n");
	printf(" -V		print version and exit\n");
	printf(" <file>		input alignment file which can optionally be gzipped\n\n");

  printf("If you use this program, please cite:\n");
  printf("\"SNP-sites: rapid efficient extraction of SNPs from multi-FASTA alignments\",\n");
  printf("Andrew J. Page, Ben Taylor, Aidan J. Delaney, Jorge Soares, Torsten Seemann, Jacqueline A. Keane, Simon R. Harris (2016),\n");
  printf("bioRxiv doi: http://dx.doi.org/10.1101/038190\n");
  
}

static void print_version()
{
  printf("%s %s\n", PROGRAM_NAME, PROGRAM_VERSION);
}

int main (int argc, char **argv) {
	char multi_fasta_filename[FILENAME_MAX] = {""};
	char output_filename[FILENAME_MAX] = {""};

	int c;
	int index;
  int output_multi_fasta_file = 0;
  int output_vcf_file = 0;
  int output_phylip_file = 0;
  int output_reference = 0;
	
	 while ((c = getopt (argc, argv, "mvrpo:V")) != -1)
      switch (c)
        {
        case 'm':
          output_multi_fasta_file = 1;
          break;
        case 'v':
          output_vcf_file = 1;
          break;
        case 'V':
          print_version();
          return 0;
        case 'p':
          output_phylip_file = 1;
          break;
        case 'r':
          output_reference = 1;
          break;
	      case 'o':
          strncpy(output_filename, optarg, FILENAME_MAX);
	        break;
        case 'h':
          print_usage();
          exit(EXIT_SUCCESS);
      default:
        output_multi_fasta_file = 1;
      }
      
  
  if(optind < argc)
  {  
    // check to see if the input alignment file exists
    if( access( argv[optind], F_OK ) == -1 ) {
      fprintf(stderr,"ERROR: cannot access input alignment file '%s'\n", argv[optind]);
      fflush(stderr);
      exit(EXIT_FAILURE);
    }
    
    strncpy(multi_fasta_filename, argv[optind], FILENAME_MAX); 
    if (output_reference) {
      generate_snp_sites_with_ref(multi_fasta_filename,
                                  output_multi_fasta_file,
                                  output_vcf_file, output_phylip_file,
                                  output_filename);
    } else {
      generate_snp_sites(multi_fasta_filename, output_multi_fasta_file,
                         output_vcf_file, output_phylip_file,
                         output_filename);
    }
  }
  else
  {
    print_usage();
  }

	exit(EXIT_SUCCESS);
}
