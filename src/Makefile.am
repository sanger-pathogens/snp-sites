snpsitesincludedir=$(includedir)/snp-sites
snpsitesinclude_HEADERS=alignment-file.h  vcf.h phylib-of-snp-sites.h snp-sites.h fasta-of-snp-sites.h parse-phylip.h string-cat.h kseq.h

bin_PROGRAMS = snp-sites
snp_sites_SOURCES = main.c


snp_sites_LDADD =-lm libsnp-sites.la -lz -lpthread
if HOST_LINUX
snp_sites_LDADD += -lrt
endif

# libsnp_sites.so is our library
lib_LTLIBRARIES = libsnp-sites.la
libsnp_sites_la_SOURCES = alignment-file.c vcf.c phylib-of-snp-sites.c snp-sites.c fasta-of-snp-sites.c parse-phylip.c string-cat.c
libsnp_sites_la_LDFLAGS= -version-info 1:0

# "make check" target
TESTS = $(check_PROGRAMS)
check_PROGRAMS = run-all-tests
run_all_tests_SOURCES = \
	../tests/check-snp-sites.c \
	../tests/check-vcf.c \
	../tests/helper-methods.c \
	../tests/run-all-tests.c
run_all_tests_CFLAGS = -I../tests

run_all_tests_LDADD = -lm -lcheck libsnp-sites.la -lz -lpthread
if HOST_LINUX
run_all_tests_LDADD += -lsubunit -lrt
endif
