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
#include <check.h>
#include <sys/stat.h>
#include <fcntl.h>
#include <unistd.h>
#include <errno.h>
#include "check-vcf.h"
#include "vcf.h"

void check_format_alternative_bases(char * test_case, char * expected_result)
{
  char * result;
  result = format_alternative_bases(test_case);
  ck_assert_str_eq(result, expected_result);
  free(result);
}

START_TEST (format_alternative_bases_test)
{
  check_format_alternative_bases("", "");
  check_format_alternative_bases("A", "A");
  check_format_alternative_bases("AC", "A,C");
  check_format_alternative_bases("ACT", "A,C,T");
}
END_TEST

void check_format_allele_index(char test_base, char reference_base, char alt_bases_array[], char expected_result_array[])
{
  char * result;
  char * alt_bases = malloc(50*sizeof(char));
  strcpy(alt_bases, alt_bases_array);
  char * expected_result = malloc(20*sizeof(char));
  strcpy(expected_result, expected_result_array);
  result = format_allele_index(test_base, reference_base, alt_bases);
  ck_assert_str_eq(result, expected_result);
  free(result);
  free(alt_bases);
  free(expected_result);
}

START_TEST (format_allele_index_test)
{
  check_format_allele_index('A', 'A', "", "0");
  check_format_allele_index('A', 'A', "C", "0");
  check_format_allele_index('A', 'A', "CA", "0");

  check_format_allele_index('A', 'C', "A", "1");
  check_format_allele_index('A', 'C', "GA", "2");

  check_format_allele_index('A', 'C', "", ".");
  check_format_allele_index('A', 'C', "G", ".");

  check_format_allele_index('A', 'B', "CDEFGHIJKLMNOPAQRST", "15");

  check_format_allele_index('-', 'A', "C", "0");
  check_format_allele_index('N', 'A', "C", "0");
  check_format_allele_index('n', 'A', "C", "0");
}
END_TEST

Suite * vcf_suite (void)
{
  Suite *s = suite_create ("Creating_VCF_file");

  TCase *tc_vcf_file = tcase_create ("vcf_file");
  tcase_add_test (tc_vcf_file, format_alternative_bases_test);
  tcase_add_test (tc_vcf_file, format_allele_index_test);
  suite_add_tcase (s, tc_vcf_file);

  return s;
}



