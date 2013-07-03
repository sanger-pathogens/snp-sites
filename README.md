Usage: snp_sites [-mvph] [-o output_filename] <file>
This program finds snp sites from a multi fasta alignment file.
 -m		output a multi fasta alignment file (default)
 -v		output a VCF file
 -p		output a phylip file
 -o		specify an output filename
 -h		this help message
 <file>		input alignment file which can optionally be gzipped


This application takes in a multi fasta alignment, finds all the SNP sites, then outputs the SNP sites in the following formats:

a multi fasta alignment,
VCF, 
relaxed phylip format.

Example usage:

snp_sites my_alignment.aln

snp_sites my_gzipped_alignment.aln.gz


Multi Fasta Alignment
=====
Similar to the input file but just containing the SNP sites.

VCF
=====
This contains the position of each SNP in the reference sequence, and the occurrence in each other sample. Can be loaded into Artemis for visualisation.

Relaxed Phylip format
=====
All the SNP sites in a format for RAxML and other tree building applications.
