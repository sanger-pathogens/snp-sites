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
#include "check-snp-sites.h"
#include "snp-sites.h"
#include "alignment-file.h"
#include "helper-methods.h"

/*Defining constant values for the reference sequence lengths*/
#define actual_ref_seq_length1 2002 //Former was 2001
#define actual_ref_seq_length2 10  //Former was 9

START_TEST (valid_alignment_with_one_line_per_sequence)
{
  generate_snp_sites("../tests/data/alignment_file_one_line_per_sequence.aln",1,1,1,"");
  fail_unless( compare_files("../tests/data/alignment_file_one_line_per_sequence.aln.vcf", "alignment_file_one_line_per_sequence.aln.vcf" ) == 1, "Invalid VCF file for 1 line per seq" );
  fail_unless( compare_files("../tests/data/alignment_file_one_line_per_sequence.aln.phylip", "alignment_file_one_line_per_sequence.aln.phylip" ) == 1, "Invalid Phylip file for 1 line per seq" );
  fail_unless( compare_files("../tests/data/alignment_file_one_line_per_sequence.aln.snp_sites.aln","alignment_file_one_line_per_sequence.aln.snp_sites.aln" ) == 1 , "Invalid ALN file for 1 line per seq");
  remove("alignment_file_one_line_per_sequence.aln.vcf");
  remove("alignment_file_one_line_per_sequence.aln.phylip");
  remove("alignment_file_one_line_per_sequence.aln.snp_sites.aln");
}
END_TEST


START_TEST (valid_alignment_with_n_as_gap)
{
	  generate_snp_sites("../tests/data/alignment_file_with_n.aln",1,1,1, "");
	  fail_unless( compare_files("../tests/data/alignment_file_with_n.aln.vcf", "alignment_file_with_n.aln.vcf" ) == 1, "Invalid VCF file for 1 line per seq" );
	  fail_unless( compare_files("../tests/data/alignment_file_with_n.aln.phylip", "alignment_file_with_n.aln.phylip" ) == 1, "Invalid Phylip file for 1 line per seq" );
	  fail_unless( compare_files("../tests/data/alignment_file_with_n.aln.snp_sites.aln","alignment_file_with_n.aln.snp_sites.aln" ) == 1 , "Invalid ALN file for 1 line per seq");
	  remove("alignment_file_with_n.aln.vcf");
	  remove("alignment_file_with_n.aln.phylip");
	  remove("alignment_file_with_n.aln.snp_sites.aln");
}
END_TEST



START_TEST (valid_alignment_with_one_line_per_sequence_gzipped)
{
  generate_snp_sites("../tests/data/alignment_file_one_line_per_sequence.aln.gz",1,1,1, NULL);
  fail_unless( compare_files("../tests/data/alignment_file_one_line_per_sequence.aln.vcf", "alignment_file_one_line_per_sequence.aln.gz.vcf" ) == 1, "Invalid VCF file for 1 line per seq" );
  fail_unless( compare_files("../tests/data/alignment_file_one_line_per_sequence.aln.phylip", "alignment_file_one_line_per_sequence.aln.gz.phylip" ) == 1, "Invalid Phylip file for 1 line per seq" );
  fail_unless( compare_files("../tests/data/alignment_file_one_line_per_sequence.aln.snp_sites.aln","alignment_file_one_line_per_sequence.aln.gz.snp_sites.aln" ) == 1 , "Invalid ALN file for 1 line per seq");
  remove("alignment_file_one_line_per_sequence.aln.gz.vcf");
  remove("alignment_file_one_line_per_sequence.aln.gz.phylip");
  remove("alignment_file_one_line_per_sequence.aln.gz.snp_sites.aln");
}
END_TEST

START_TEST (valid_alignment_with_multiple_lines_per_sequence)
{
    generate_snp_sites("../tests/data/alignment_file_multiple_lines_per_sequence.aln",1,1,1,"");
    fail_unless( compare_files("../tests/data/alignment_file_one_line_per_sequence.aln.vcf", "alignment_file_multiple_lines_per_sequence.aln.vcf" ) == 1, "Invalid VCF file for multiple lines per seq" );
    fail_unless( compare_files("../tests/data/alignment_file_one_line_per_sequence.aln.phylip", "alignment_file_multiple_lines_per_sequence.aln.phylip" ) == 1, "Invalid Phylip file for multiple lines per seq" );
    fail_unless( compare_files("../tests/data/alignment_file_one_line_per_sequence.aln.snp_sites.aln","alignment_file_multiple_lines_per_sequence.aln.snp_sites.aln" ) == 1 ,"Invalid ALN file for multiple lines per seq");
    remove("alignment_file_multiple_lines_per_sequence.aln.vcf");
    remove("alignment_file_multiple_lines_per_sequence.aln.phylip");
    remove("alignment_file_multiple_lines_per_sequence.aln.snp_sites.aln");
}
END_TEST

START_TEST (valid_alignment_with_pure_mode)
{
    generate_snp_sites_with_ref_pure_mono("../tests/data/pure_mode_alignment.aln",1,1,1,"",0,1,0);
    fail_unless( compare_files("../tests/data/pure_mode_alignment.aln.vcf", "pure_mode_alignment.aln.vcf" ) == 1, "Invalid VCF file for multiple lines per seq" );
    fail_unless( compare_files("../tests/data/pure_mode_alignment.aln.phylip", "pure_mode_alignment.aln.phylip" ) == 1, "Invalid Phylip file for multiple lines per seq" );
    fail_unless( compare_files("../tests/data/pure_mode_alignment.aln.snp_sites.aln","pure_mode_alignment.aln.snp_sites.aln" ) == 1 ,"Invalid ALN file for multiple lines per seq");
    remove("pure_mode_alignment.aln.vcf");
    remove("pure_mode_alignment.aln.phylip");
    remove("pure_mode_alignment.aln.snp_sites.aln");
}
END_TEST
  
START_TEST (valid_alignment_with_monomorphic_sites)
{
    generate_snp_sites_with_ref_pure_mono("../tests/data/pure_mode_monomorphic_alignment.aln",1,1,1,"",0,1,1);
    fail_unless( compare_files("../tests/data/pure_mode_monomorphic_alignment.aln.vcf", "pure_mode_monomorphic_alignment.aln.vcf" ) == 1, "Invalid VCF file for multiple lines per seq" );
    fail_unless( compare_files("../tests/data/pure_mode_monomorphic_alignment.aln.phylip", "pure_mode_monomorphic_alignment.aln.phylip" ) == 1, "Invalid Phylip file for multiple lines per seq" );
    fail_unless( compare_files("../tests/data/pure_mode_monomorphic_alignment.aln.snp_sites.aln","pure_mode_monomorphic_alignment.aln.snp_sites.aln" ) == 1 ,"Invalid ALN file for multiple lines per seq");
    remove("pure_mode_monomorphic_alignment.aln.vcf");
    remove("pure_mode_monomorphic_alignment.aln.phylip");
    remove("pure_mode_monomorphic_alignment.aln.snp_sites.aln");
}
END_TEST
  

START_TEST (valid_with_only_aln_file_output_default)
{
	    generate_snp_sites("../tests/data/alignment_file_one_line_per_sequence.aln",0,0,0,"");
	    fail_unless( compare_files("../tests/data/alignment_file_one_line_per_sequence.aln.snp_sites.aln","alignment_file_one_line_per_sequence.aln.snp_sites.aln" ) == 1 ,"Invalid ALN file for multiple lines per seq");
	    remove("alignment_file_one_line_per_sequence.aln.snp_sites.aln");
}
END_TEST

START_TEST (valid_with_only_aln_file_output)
{
	    generate_snp_sites("../tests/data/alignment_file_one_line_per_sequence.aln",1,0,0,"");
	    fail_unless( compare_files("../tests/data/alignment_file_one_line_per_sequence.aln.snp_sites.aln","alignment_file_one_line_per_sequence.aln.snp_sites.aln" ) == 1 ,"Invalid ALN file for multiple lines per seq");
	    remove("alignment_file_one_line_per_sequence.aln.snp_sites.aln");
}
END_TEST

START_TEST (valid_with_only_aln_file_output_with_custom_name)
{
	char output_filename[100] = {"some_custom_name"};
	    generate_snp_sites("../tests/data/alignment_file_multiple_lines_per_sequence.aln",1,0,0,output_filename);
	    fail_unless( compare_files("../tests/data/alignment_file_one_line_per_sequence.aln.snp_sites.aln","some_custom_name" ) == 1 ,"Only aln with custom name should be outputted");
	    remove("some_custom_name");
}
END_TEST

START_TEST (valid_with_phylip_outputted_with_custom_name)
{
  generate_snp_sites("../tests/data/alignment_file_one_line_per_sequence.aln",0,0,1,"some_custom_name");
  fail_unless( compare_files("../tests/data/alignment_file_one_line_per_sequence.aln.phylip", "some_custom_name" ) == 1, "Custom name needs extra extension for phylip" );
  remove("some_custom_name");

}
END_TEST

START_TEST (invalid_with_uneven_file_lengths)
{
	    generate_snp_sites("../tests/data/uneven_alignment.aln",1,0,0,"");
}
END_TEST


START_TEST (valid_with_all_outputted_with_custom_name)
{
  generate_snp_sites("../tests/data/alignment_file_one_line_per_sequence.aln",1,1,1,"some_custom_name");
  fail_unless( compare_files("../tests/data/alignment_file_one_line_per_sequence.aln.vcf", "some_custom_name.vcf" ) == 1, "Custom name needs extra extension for VCF" );
  fail_unless( compare_files("../tests/data/alignment_file_one_line_per_sequence.aln.phylip", "some_custom_name.phylip" ) == 1, "Custom name needs extra extension for phylip" );
  fail_unless( compare_files("../tests/data/alignment_file_one_line_per_sequence.aln.snp_sites.aln","some_custom_name.snp_sites.aln" ) == 1 , "Custom name needs extra extension for ALN");
  remove("some_custom_name.vcf");
  remove("some_custom_name.phylip");
  remove("some_custom_name.snp_sites.aln");
}
END_TEST

  
START_TEST (valid_aln_plus_reference)
{
	    generate_snp_sites_with_ref("../tests/data/alignment_file_one_line_per_sequence.aln",1,0,0,"");
	    fail_unless( compare_files("../tests/data/alignment_file_one_line_per_sequence_reference.snp_sites.aln","alignment_file_one_line_per_sequence.aln.snp_sites.aln" ) == 1 ,"Invalid ALN file for multiple lines per seq");
	    remove("alignment_file_one_line_per_sequence.aln.snp_sites.aln");
}
END_TEST

START_TEST (valid_phylip_plus_reference)
{
  generate_snp_sites_with_ref("../tests/data/alignment_file_one_line_per_sequence.aln",0,0,1,"");
  fail_unless( compare_files("../tests/data/alignment_file_one_line_per_sequence_reference.aln.phylip", "alignment_file_one_line_per_sequence.aln.phylip" ) == 1, "invalid phylib with reference" );
  remove("alignment_file_one_line_per_sequence.aln.phylip");

}
END_TEST

START_TEST (valid_genome_length)
{
  detect_snps("../tests/data/alignment_file_one_line_per_sequence.aln",0,0);
  fail_unless( get_length_of_genome() == 2000 );
}
END_TEST

START_TEST (valid_genome_length_with_multiple_lines_per_sequence)
{
  detect_snps("../tests/data/alignment_file_multiple_lines_per_sequence.aln",0,0);
  fail_unless( get_length_of_genome() == 2000 );
}
END_TEST

START_TEST (valid_number_of_sequences_in_file)
{
  detect_snps("../tests/data/alignment_file_one_line_per_sequence.aln",0,0);
  fail_unless( get_number_of_samples() == 109 );
}
END_TEST

START_TEST (valid_number_of_sequences_in_file_with_multiple_lines_per_sequence)
{
  detect_snps("../tests/data/alignment_file_multiple_lines_per_sequence.aln",0,0);
  fail_unless( get_number_of_samples() == 109 );
}
END_TEST 

START_TEST (number_of_snps_detected)
{
  detect_snps("../tests/data/alignment_file_multiple_lines_per_sequence.aln",0,0);
  fail_unless( get_number_of_snps()  == 5);
}
END_TEST

START_TEST (number_of_snps_detected_small)
{
  detect_snps("../tests/data/small_alignment.aln",0,0);
  fail_unless(  get_number_of_snps()  == 1);
}
END_TEST
  
START_TEST (detect_snps_pure_mode)
{
  detect_snps("../tests/data/pure_mode_alignment.aln",1,0);
  fail_unless(  get_number_of_snps()  == 2);
}
END_TEST
  
  
START_TEST (detect_snps_pure_mode_monomorphic)
{
  detect_snps("../tests/data/pure_mode_monomorphic_alignment.aln",1,1);
  fail_unless(  get_number_of_snps()  == 3);
}
END_TEST

START_TEST (sample_names_from_alignment_file)
{
  detect_snps("../tests/data/small_alignment.aln",0,0);
  char ** current_sequence_names = get_sequence_names();

  fail_unless(strcmp(current_sequence_names[0],"reference_sequence") == 0);
  fail_unless(strcmp(current_sequence_names[1],"comparison_sequence") == 0);
  fail_unless(strcmp(current_sequence_names[2],"another_comparison_sequence") == 0);
}
END_TEST

START_TEST (check_strip_directory_from_filename_without_directory)
{
	char *input_filename_without_directory =  "my_file_name.aln";
	char output_filename[30];
	memset(output_filename,'\0',30);
	strip_directory_from_filename(input_filename_without_directory, output_filename);
	fail_unless( strcmp(input_filename_without_directory, output_filename) ==0 );
}
END_TEST

START_TEST (check_strip_directory_from_filename_with_directory)
{
	char *input_filename_without_directory =  "/some/directory/name/my_file_name.aln";
	char output_filename[30];
	memset(output_filename,'\0',30);
	strip_directory_from_filename(input_filename_without_directory, output_filename);
	fail_unless( strcmp("my_file_name.aln", output_filename) ==0 );
}
END_TEST

Suite * snp_sites_suite (void)
{
  Suite *s = suite_create ("Creating_SNP_Sites");

  TCase *tc_alignment_file = tcase_create ("alignment_file");
  tcase_add_test (tc_alignment_file, valid_genome_length);
  tcase_add_test (tc_alignment_file, valid_genome_length_with_multiple_lines_per_sequence);
  tcase_add_test (tc_alignment_file, valid_number_of_sequences_in_file);
  tcase_add_test (tc_alignment_file, valid_number_of_sequences_in_file_with_multiple_lines_per_sequence);
  tcase_add_test (tc_alignment_file, number_of_snps_detected_small);
  tcase_add_test (tc_alignment_file, number_of_snps_detected);
  tcase_add_test (tc_alignment_file, sample_names_from_alignment_file);
  tcase_add_test (tc_alignment_file, detect_snps_pure_mode);
  tcase_add_test (tc_alignment_file, check_strip_directory_from_filename_without_directory);
  tcase_add_test (tc_alignment_file, check_strip_directory_from_filename_with_directory);
  tcase_add_test (tc_alignment_file, detect_snps_pure_mode_monomorphic);
  
  suite_add_tcase (s, tc_alignment_file);

  TCase *tc_snp_sites = tcase_create ("snp_sites");
  tcase_add_test (tc_snp_sites, valid_alignment_with_one_line_per_sequence);
  tcase_add_test (tc_snp_sites, valid_alignment_with_multiple_lines_per_sequence);
  tcase_add_test (tc_snp_sites, valid_alignment_with_one_line_per_sequence_gzipped);
  tcase_add_test (tc_snp_sites, valid_alignment_with_n_as_gap);
  tcase_add_test (tc_snp_sites, valid_with_only_aln_file_output);
  tcase_add_test (tc_snp_sites, valid_with_only_aln_file_output_with_custom_name);
  tcase_add_test (tc_snp_sites, valid_with_all_outputted_with_custom_name);
  tcase_add_test (tc_snp_sites, valid_with_only_aln_file_output_default);
  tcase_add_test (tc_snp_sites, valid_with_phylip_outputted_with_custom_name);
  tcase_add_test (tc_snp_sites, valid_aln_plus_reference);
  tcase_add_test (tc_snp_sites, valid_phylip_plus_reference);
  tcase_add_test (tc_snp_sites, valid_alignment_with_pure_mode);
  tcase_add_test (tc_snp_sites, valid_alignment_with_monomorphic_sites);
  
  tcase_add_exit_test(tc_snp_sites, invalid_with_uneven_file_lengths,EXIT_FAILURE);
  remove("uneven_alignment.aln.snp_sites.aln");
  suite_add_tcase (s, tc_snp_sites);

  return s;
}

