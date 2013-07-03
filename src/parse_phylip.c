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
#include "parse_phylip.h"
#include "alignment_file.h"
#include "string_cat.h"

int num_samples;
int num_snps;
char ** sequences;
char ** phylip_sample_names;


int does_column_contain_snps(int snp_column, char reference_base)
{
	int i;
	for(i = 0; i < num_samples; i++)
	{
		if(sequences[i][snp_column] == '\0' || sequences[i][snp_column] == '\n')
		{
			return 0;	
		}
		
		if(sequences[i][snp_column]  != '-' && toupper(sequences[i][snp_column])  != 'N' && sequences[i][snp_column] != reference_base)
		{
			return 1;
		}
	}
	return 0;
}
