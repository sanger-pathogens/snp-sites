# SNP-sites
Rapidly extracts SNPs from a multi-FASTA alignment.

[![Build Status](https://travis-ci.org/sanger-pathogens/snp-sites.png?branch=master)](https://travis-ci.org/sanger-pathogens/snp-sites)   
[![License: GPL v3](https://img.shields.io/badge/License-GPL%20v3-brightgreen.svg)](https://github.com/sanger-pathogens/snp-sites/blob/master/LICENSE)   
[![status](https://img.shields.io/badge/MGEN-10.1099%2Fmgen.0.000056-brightgreen.svg)](http://mgen.microbiologyresearch.org/content/journal/mgen/10.1099/mgen.0.000056)   
[![install with bioconda](https://img.shields.io/badge/install%20with-bioconda-brightgreen.svg)](http://bioconda.github.io/recipes/gubbins/README.html)  
[![Container ready](https://img.shields.io/badge/container-ready-brightgreen.svg)](https://quay.io/repository/biocontainers/gubbins)  
[![Docker Build Status](https://img.shields.io/docker/build/sangerpathogens/gubbins.svg)](https://hub.docker.com/r/sangerpathogens/gubbins)  
[![Docker Pulls](https://img.shields.io/docker/pulls/sangerpathogens/gubbins.svg)](https://hub.docker.com/r/sangerpathogens/gubbins)  
[![codecov](https://codecov.io/gh/sanger-pathogens/snp-sites/branch/master/graph/badge.svg)](https://codecov.io/gh/sanger-pathogens/snp-sites)   


## Contents
* [Introduction](#introduction)
* [Installation](#installation)
  * [Linux \- Ubuntu/Debian](#linux---ubuntudebian)
  * [OSX/Linux \- using Bioconda](#osxlinux---using-bioconda)
  * [OSX/Linux \- from source](#osxlinux---from-source)
  * [OSX/Linux \- from a release tarball](#osxlinux---from-a-release-tarball)
  * [All platforms \- Docker](#all-platforms---docker)
  * [Running the tests](#running-the-tests)
* [Usage](#usage)
  * [Example input](#example-input)
  * [Example usage](#example-usage)
  * [Output](#output)
* [License](#license)
* [Feedback/Issues](#feedbackissues)
* [Citation](#citation)

## Introduction
Rapidly decreasing genome sequencing costs have led to a proportionate increase in the number of samples used in prokaryotic population studies. Extracting single nucleotide polymorphisms (SNPs) from a large whole genome alignment is now a routine task, but existing tools have failed to scale efficiently with the increased size of studies. These tools are slow, memory inefficient and are installed through non-standard procedures. We present SNP-sites which can rapidly extract SNPs from a multi-FASTA alignment using modest resources and can output results in multiple formats for downstream analysis. SNPs can be extracted from a 8.3 GB alignment file (1,842 taxa, 22,618 sites) in 267 seconds using 59 MB of RAM and 1 CPU core, making it feasible to run on modest computers. It is easy to install through the Debian and Homebrew package managers, and has been successfully tested on more than 20 operating systems. SNP-sites is implemented in C and is available under the open source license GNU GPL version 3.

## Installation
There are a few ways to install SNP-sites. The simpliest way is using apt (Debian/Ubuntu) or Conda. If you encounter an issue when installing SNP-sites please contact your local system administrator. If you encounter a bug please log it [here](https://github.com/sanger-pathogens/snp-sites/issues) or email us at snpsites-help@sanger.ac.uk

* Linux - Ubuntu/Debian
* OSX/Linux - using Bioconda
* OSX/Linux - from source
* OSX/Linux - from a release tarball

### Linux - Ubuntu/Debian
If you have a recent version of Ubuntu or Debian then you can install it using apt.
```
   apt-get install snp-sites
```

### OSX/Linux - using Bioconda
Install Conda and install the bioconda channels.
```
conda config --add channels conda-forge
conda config --add channels defaults
conda config --add channels r
conda config --add channels bioconda
conda install snp-sites
```

### OSX/Linux - from source
This is a difficult method and is only suitable for someone with advanced unix skills. No support is provided with this method, since you have advanced unix skills. Please consider using Conda instead. First install a standard development environment (e.g. gcc, automa\
ke, autoconf, libtool). Download the software from [GitHub](https://github.com/sanger-pathogens/snp-sites).

```
autoreconf -i -f
./configure
make
sudo make install
```

### OSX/Linux - from a release tarball
This is a difficult method and is only suitable for someone with advanced unix skills. No support is provided with this method, since you have advanced unix skills. Please consider using Conda instead. First install a standard development environment (e.g. gcc, automa\
ke, autoconf, libtool).

```
tar xzvf snp-sites-x.y.z.tar.gz
cd snp-sites-x.y.z
./configure
make
sudo make install
```

### All platforms - Docker
Bioconda produce a Docker container so you can use the software out of the box. Install Docker and then pull the container from Bioconda https://quay.io/repository/biocontainers/snp-sites

### Running the tests
The test can be run from the top level directory:  

```
autoreconf -i
./configure
make
make check
```

This requires libcheck (the `check` package in Ubuntu) to be installed.

## Usage

```
Usage: snp-sites [-mvph] [-o output_filename] <file>
This program finds snp sites from a multi fasta alignment file.
 -r     output internal pseudo reference sequence
 -m     output a multi fasta alignment file (default)
 -v     output a VCF file
 -p     output a phylip file
 -o STR specify an output filename [STDOUT]
 -c     only output columns containing exclusively ACGT
 -b     output monomorphic sites, used for BEAST
 -h     this help message
 -V     print version and exit
 <file> input alignment file which can optionally be gzipped
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
### Output
* Multi Fasta Alignment - Similar to the input file but just containing the SNP sites.

* VCF - This contains the position of each SNP in the reference sequence, and the occurrence in each other sample. Can be loaded into Artemis for visualisation.

* Relaxed Phylip format - All the SNP sites in a format for RAxML and other tree building applications.

## License
SNP-sites is free software, licensed under [GPLv3](https://github.com/sanger-pathogens/snp-sites/blob/master/LICENSE).

## Feedback/Issues
Please report any issues to the [issues page](https://github.com/sanger-pathogens/snp-sites/issues) or email snpsites-help@sanger.ac.uk

## Citation
If you use this software please cite:

"SNP-sites: rapid efficient extraction of SNPs from multi-FASTA alignments", Andrew J. Page, Ben Taylor, Aidan J. Delaney, Jorge Soares, Torsten Seemann, Jacqueline A. Keane, Simon R. Harris, [Microbial Genomics 2(4), (2016)](http://mgen.microbiologyresearch.org/content/journal/mgen/10.1099/mgen.0.000056)