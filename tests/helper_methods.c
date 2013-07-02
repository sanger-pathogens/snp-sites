#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/stat.h>
#include <fcntl.h>
#include <unistd.h>
#include <errno.h>
#include "helper_methods.h"


int compare_files(char expected_output_filename[],char actual_output_filename[] )
{
  FILE *expected_output_fh;
  FILE *actual_output_fh;
  
  char    *expected_buffer;
  char    *actual_buffer;
  long    numbytes;
  
  expected_output_fh = fopen(expected_output_filename, "r");
  actual_output_fh = fopen(actual_output_filename, "r");
  
  fseek(expected_output_fh, 0L, SEEK_END);
  numbytes = ftell(expected_output_fh);
  fseek(expected_output_fh, 0L, SEEK_SET);	
  expected_buffer = (char*)calloc(numbytes, sizeof(char));	
  fread(expected_buffer, sizeof(char), numbytes, expected_output_fh);
  fclose(expected_output_fh);
  
  fseek(actual_output_fh, 0L, SEEK_END);
  numbytes = ftell(actual_output_fh);
  fseek(actual_output_fh, 0L, SEEK_SET);	
  actual_buffer = (char*)calloc(numbytes, sizeof(char));	
  fread(actual_buffer, sizeof(char), numbytes, actual_output_fh);
  fclose(actual_output_fh);
  
  if(strcmp(expected_buffer,actual_buffer) == 0)
  { 
    free(expected_buffer);
    free(actual_buffer);
    return 1;
  }

  free(expected_buffer);
  free(actual_buffer);
  
  return 0;
}
