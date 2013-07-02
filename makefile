gubbins: alignment_file.o main.o  vcf.o phylib_of_snp_sites.o snp_sites.o fasta_of_snp_sites.o parse_phylip.o string_cat.o
	gcc -o snp_sites alignment_file.o main.o  vcf.o phylib_of_snp_sites.o snp_sites.o fasta_of_snp_sites.o parse_phylip.o string_cat.o -lm -lz

string_cat.o: string_cat.c
	gcc -c string_cat.c

alignment_file.o: alignment_file.c
	gcc -c alignment_file.c

fasta_of_snp_sites.o: 	fasta_of_snp_sites.c
	gcc -c fasta_of_snp_sites.c

main.o: main.c
	gcc -c main.c

phylib_of_snp_sites.o: 	phylib_of_snp_sites.c
	gcc -c phylib_of_snp_sites.c

snp_sites.o: snp_sites.c
	gcc -c snp_sites.c

parse_phylip.o: parse_phylip.c
	gcc -c parse_phylip.c

vcf.o: vcf.c
	gcc -c vcf.c

clean:
	-rm *.o

test:
	cd tests && make

