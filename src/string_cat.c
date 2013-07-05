/*
 *  Wellcome Trust Sanger Institute
 *  Copyright (C) 2013  Wellcome Trust Sanger Institute
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


int size_of_string(char *input_string)
{
	int i = 0;

	while( input_string[i] != '\0')
	{
		i++;
	}
	return i;
}

void concat_strings_created_with_malloc(char *input_string, char *string_to_concat)
{
	int input_str_size = size_of_string(input_string);
	int to_concat_str_size = size_of_string(string_to_concat);
	memcpy(input_string + input_str_size, string_to_concat, (to_concat_str_size+1)*sizeof(char));
}
