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
#include <zlib.h>
#include <regex.h>
#include <sys/types.h>
#include "kseq.h"
#include "vcf.h"
#include "alignment-file.h"
#include "snp-sites.h"

KSEQ_INIT(gzFile, gzread)

int length_of_genome;
int number_of_samples;
int number_of_snps;
char ** sequence_names;
int * snp_locations;


int get_length_of_genome()
{
    return length_of_genome;
}

int get_number_of_samples()
{
    return number_of_samples;
}

int get_number_of_snps()
{
    return number_of_snps;
}

char ** get_sequence_names()
{
    return sequence_names;
}

int * get_snp_locations()
{
    return snp_locations;
}

void get_bases_for_each_snp(char filename[], int snp_locations[], char ** bases_for_snps, size_t length_of_genome, int number_of_snps)
{
  int l;
  int i = 0;
  int sequence_number = 0;
  size_t length_of_genome_found =0;
	
  gzFile fp;
  kseq_t *seq;
	
  fp = gzopen(filename, "r");
  seq = kseq_init(fp);

  while ((l = kseq_read(seq)) >= 0) 
    {
      if(sequence_number == 0)
	{
	  length_of_genome_found = seq->seq.l;
	}
      for(i = 0; i< number_of_snps; i++)
	{
	  bases_for_snps[i][sequence_number] = toupper(((char *) seq->seq.s)[snp_locations[i]]);
	}
      
      if(seq->seq.l != length_of_genome_found)
	{
	  fprintf(stderr, "Alignment %s contains sequences of unequal length. Expected length is %i but got %i in sequence %s\n\n",filename, (int) length_of_genome_found, (int) seq->seq.l,seq->name.s);
	  fflush(stderr);
	  exit(EXIT_FAILURE);
	}
		
      sequence_number++;
    }

  kseq_destroy(seq);
  gzclose(fp);
}

void detect_snps(char filename[])
{
  int i;
  int l;
  number_of_snps = 0;
  number_of_samples = 0; 
  length_of_genome = 0;
  
  gzFile fp;
  kseq_t *seq;
  kseq_t *first_sequence;
  
  fp = gzopen(filename, "r");
  seq_init(fp);  
  // First sequence is the reference sequence so skip it
  kseq_read(first_sequence);
  length_of_genome = first_sequence->seq.l;
  number_of_samples++; 
  sequence_names = calloc(number_of_samples, sizeof(char*));
  sequence_names[number_of_samples-1] = calloc(MAX_SAMPLE_NAME_SIZE,sizeof(char));
  strcpy(sequence_names[number_of_samples-1], first_sequence->name.s);
  
  while ((l = kseq_read(seq)) >= 0) {
    for(i = 0; i < length_of_genome; i++)
    {
	  // If there is an indel in the reference sequence, replace with the first proper base you find
	  if((first_sequence->seq.s[i] == '-' && seq->seq.s[i] != '-' ) || (toupper(first_sequence->seq.s[i]) == 'N' && seq->seq.s[i] != 'N' ))
	  {
	      first_sequence->seq.s[i] = toupper(seq->seq.s[i]);
	  }
	  
	  if(first_sequence->seq.s[i] != '*' && seq->seq.s[i] != '-' && toupper(first_sequence->seq.s[i]) != 'N' && toupper(first_sequence->seq.s[i]) != toupper(seq->seq.s[i]))
	  {
	      first_sequence->seq.s[i] = '*';
	      number_of_snps++;
	  }
   }
   sequence_names = realloc(sequence_names, number_of_samples * sizeof(char*));
   sequence_names[number_of_samples-1] = calloc(MAX_SAMPLE_NAME_SIZE,sizeof(char));
   strcpy(sequence_names[number_of_samples-1], seq->name.s);
   
   number_of_samples++;
  }
  
  int current_snp_index = 0;
  snp_locations = calloc(number_of_snps, sizeof(int));
  for(i = 0; i < length_of_genome; i++)
  {
      if(first_sequence->seq.s[i] == '*')
      {
          snp_locations[current_snp_index] = i;
          current_snp_index++;
      }
  }
  
  kseq_destroy(first_sequence);
  kseq_destroy(seq);
  gzclose(fp);

  return;
}

