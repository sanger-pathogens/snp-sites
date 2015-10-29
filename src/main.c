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
#include "snp-sites.h"
#include "string-cat.h"

#define PROGRAM_NAME "snp-sites"
#define PROGRAM_VERSION "2.0.2"

static void print_usage()
{
	printf("Usage: snp_sites [-mvph] [-o output_filename] <file>\n");
	printf("This program finds snp sites from a multi fasta alignment file.\n");
	printf(" -m		output a multi fasta alignment file (default)\n");
	printf(" -v		output a VCF file\n");
	printf(" -p		output a phylip file\n");
	printf(" -o		specify an output filename\n");
	printf(" -h		this help message\n");
	printf(" -V		print version and exit\n");
	printf(" <file>		input alignment file which can optionally be gzipped\n");
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
	
	 while ((c = getopt (argc, argv, "mvpo:V")) != -1)
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
	      case 'o':
          strncpy(output_filename, optarg, FILENAME_MAX);
	        break;
        case 'h':
          print_usage();
          return 0;
      default:
        output_multi_fasta_file = 1;
      }
  
  if(optind < argc)
  {
    strncpy(multi_fasta_filename, argv[optind], FILENAME_MAX); 
    generate_snp_sites(multi_fasta_filename, output_multi_fasta_file, output_vcf_file, output_phylip_file, output_filename);
  }
  else
  {
    print_usage();
  }
    
	return 0;
}


