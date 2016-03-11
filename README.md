[![Build Status](https://travis-ci.org/sanger-pathogens/snp-sites.png?branch=master)](https://travis-ci.org/sanger-pathogens/snp-sites)
# SNP Sites
Rapidly decreasing genome sequencing costs have led to a proportionate increase in the number of samples used in prokaryotic population studies. Extracting single nucleotide polymorphisms (SNPs) from a large whole genome alignment is now a routine task, but existing tools have failed to scale efficiently with the increased size of studies. These tools are slow, memory inefficient and are installed through non-standard procedures. We present SNP-sites which can rapidly extract SNPs from a multi-FASTA alignment using modest resources and can output results in multiple formats for downstream analysis. SNPs can be extracted from a 8.3 GB alignment file (1,842 taxa, 22,618 sites) in 267 seconds using 59 MB of RAM and 1 CPU core, making it feasible to run on modest computers. It is easy to install through the Debian and Homebrew package managers, and has been successfully tested on more than 20 operating systems. SNP-sites is implemented in C and is available under the open source license GNU GPL version 3.

The software can be cited as:

[Andrew J. Page, Ben Taylor, Aidan J. Delaney, Jorge Soares, Torsten Seemann, Jacqueline A. Keane, Simon R. Harris, "SNP-sites: rapid efficient extraction of SNPs from multi-FASTA alignments", bioRxiv, doi: http://dx.doi.org/10.1101/038190 (2016).](http://biorxiv.org/content/early/2016/01/29/038190)

```
Usage: snp-sites [-mvph] [-o output_filename] <file>
This program finds snp sites from a multi fasta alignment file.
 -r		    output internal pseudo reference sequence
 -m		    output a multi fasta alignment file (default)
 -v		    output a VCF file
 -p		    output a phylip file
 -o		    specify an output filename
 -h		    this help message
 -V		    print version and exit
 <file>		input alignment file which can optionally be gzipped
```

This application takes in a multi fasta alignment, finds all the SNP sites, then outputs the SNP sites in the following formats:

- a multi fasta alignment,
- VCF, 
- relaxed phylip format.

### Example input
For the given input file:
```
>sample1
AGACACAGTCAC
>sample1
AGACAC----AC
>sample1
AAACGCATTCAN
```
the output is:
```
>sample1
GAG
>sample1
GA-
>sample1
AGT
```

### Example usage

```
snp-sites my_alignment.aln
snp-sites my_gzipped_alignment.aln.gz
```

### Multi Fasta Alignment
Similar to the input file but just containing the SNP sites.

### VCF
This contains the position of each SNP in the reference sequence, and the occurrence in each other sample. Can be loaded into Artemis for visualisation.

### Relaxed Phylip format
All the SNP sites in a format for RAxML and other tree building applications.
