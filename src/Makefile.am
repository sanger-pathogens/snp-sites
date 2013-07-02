# ship these headers in the "make dist" target
dist_data_DATA = alignment_file.h  vcf.h phylib_of_snp_sites.h snp_sites.h fasta_of_snp_sites.h parse_phylip.h string_cat.h 

# snp_sites is our top level progra
bin_PROGRAMS = snp_sites
snp_sites_SOURCES = main.c
snp_sites_LDADD=-lm libsnp_sites.la -lz

# libsnp_sites.so is our library
lib_LTLIBRARIES = libsnp_sites.la
libsnp_sites_la_SOURCES = alignment_file.c vcf.c phylib_of_snp_sites.c snp_sites.c fasta_of_snp_sites.c parse_phylip.c string_cat.c
libsnp_sites_la_LDFLAGS= -version-info 1:0

# "make check" target
TESTS = $(check_PROGRAMS)
check_PROGRAMS = run_all_tests
run_all_tests_SOURCES = \
	../tests/check_snp_sites.c \
	../tests/helper_methods.c \
	../tests/run_all_tests.c
run_all_tests_CFLAGS = -I../tests
run_all_tests_LDADD = -lm -lcheck libsnp_sites.la -lz