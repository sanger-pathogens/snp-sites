snpsitesincludedir=$(includedir)/snp-sites
snpsitesinclude_HEADERS=alignment-file.h  vcf.h phylib-of-snp-sites.h snp-sites.h fasta-of-snp-sites.h kseq.h

bin_PROGRAMS = snp-sites
snp_sites_SOURCES = main.c
snp_sites_LDADD = libsnp-sites.la -lpthread

# libsnp_sites.so is our library
lib_LTLIBRARIES = libsnp-sites.la
libsnp_sites_la_SOURCES = alignment-file.c vcf.c phylib-of-snp-sites.c snp-sites.c fasta-of-snp-sites.c
libsnp_sites_la_LDFLAGS= -version-info 1:0

# "make check" target
TESTS = $(check_PROGRAMS)
check_PROGRAMS = run-all-tests
run_all_tests_SOURCES = \
	../tests/check-snp-sites.c \
	../tests/check-vcf.c \
	../tests/helper-methods.c \
	../tests/run-all-tests.c
run_all_tests_CFLAGS = $(CHECK_CFLAGS) $(CFLAGS)
run_all_tests_LDADD = libsnp-sites.la -lpthread $(CHECK_LIBS)
