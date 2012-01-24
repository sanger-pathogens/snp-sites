test:
	gcc -c check_snp_sites.c  -I .
	gcc -o check_snp_sites check_snp_sites.o ../alignment_file.o ../vcf.o ../parse_phylip.o ../phylib_of_snp_sites.o ../fasta_of_snp_sites.o ../snp_sites.o -lm -lz -L /usr/local/lib/libcheck.dylib
	time ./check_snp_sites

