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
#include "snp_sites.h"

#define MAX_FILENAME_SIZE 250

static void print_usage()
{
	printf("Find SNP sites from a multi fasta alignment file (which can be gzipped)\n");
	printf("./snp_sites file.aln\n");
	printf("./snp_sites file.aln.gz\n\n");
	
	printf("SNP sites are outputted in the following formats: Multi fasta alignment, phylip, VCF \n\n");
}

int main (int argc, const char * argv[]) {
	char multi_fasta_filename[MAX_FILENAME_SIZE];

  if(argc <=1)
  {
    print_usage();
  }
	else if(strcmp(argv[1], "--help") == 0)
  {
	  print_usage();
  }
  else
  {
    strcpy(multi_fasta_filename,argv[1]);
		generate_snp_sites(multi_fasta_filename);
  }

    
	return 0;
}


