# ship these headers in the "make dist" target
dist_data_DATA = alignment-file.h  vcf.h phylib-of-snp-sites.h snp-sites.h fasta-of-snp-sites.h parse-phylip.h string-cat.h kseq.h

# snp_sites is our top level progra
bin_PROGRAMS = snp-sites
snp_sites_SOURCES = main.c
snp_sites_LDADD=-lm libsnp-sites.la -lz -lpthread -lrt -autoreconf

# libsnp_sites.so is our library
lib_LTLIBRARIES = libsnp-sites.la
libsnp_sites_la_SOURCES = alignment-file.c vcf.c phylib-of-snp-sites.c snp-sites.c fasta-of-snp-sites.c parse-phylip.c string-cat.c
libsnp_sites_la_LDFLAGS= -version-info 1:0

# "make check" target
TESTS = $(check_PROGRAMS)
check_PROGRAMS = run-all-tests
run_all_tests_SOURCES = \
	../tests/check-snp-sites.c \
	../tests/helper-methods.c \
	../tests/run-all-tests.c
run_all_tests_CFLAGS = -I../tests
run_all_tests_LDADD = -lm -lcheck libsnp-sites.la -lz -lpthread -lrt -autoreconf