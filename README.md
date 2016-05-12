[![Build Status](https://travis-ci.org/sanger-pathogens/snp-sites.png?branch=master)](https://travis-ci.org/sanger-pathogens/snp-sites)
# SNP-sites
Rapidly decreasing genome sequencing costs have led to a proportionate increase in the number of samples used in prokaryotic population studies. Extracting single nucleotide polymorphisms (SNPs) from a large whole genome alignment is now a routine task, but existing tools have failed to scale efficiently with the increased size of studies. These tools are slow, memory inefficient and are installed through non-standard procedures. We present SNP-sites which can rapidly extract SNPs from a multi-FASTA alignment using modest resources and can output results in multiple formats for downstream analysis. SNPs can be extracted from a 8.3 GB alignment file (1,842 taxa, 22,618 sites) in 267 seconds using 59 MB of RAM and 1 CPU core, making it feasible to run on modest computers. It is easy to install through the Debian and Homebrew package managers, and has been successfully tested on more than 20 operating systems. SNP-sites is implemented in C and is available under the open source license GNU GPL version 3.

The software can be cited as:

["SNP-sites: rapid efficient extraction of SNPs from multi-FASTA alignments", Andrew J. Page, Ben Taylor, Aidan J. Delaney, Jorge Soares, Torsten Seemann, Jacqueline A. Keane, Simon R. Harris, Microbial Genomics 2(4), (2016)](http://dx.doi.org/10.1099/mgen.0.000056)

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

# Installation
There are a few ways to install snp-sites. The simpliest way is using apt (Debian/Ubuntu), HomeBrew (OSX) or LinuxBrew.

* Linux - Ubuntu/Debian
* Linux - CentOS/RHEL
* Linux - using LinuxBrew
* OSX - using HomeBrew
* OSX/Linux - from source
* OSX/Linux - from a release tarball
* Windows/OSX/Linux - using a Virtual Machine

## Linux - Ubuntu/Debian
If you have a recent version of Ubuntu or Debian then you can install it using apt.
```
   apt-get install snp-sites
```

## Linux - CentOS/RHEL
The easiest way to install on CentOS/RHEL is to use [LinuxBrew](http://brew.sh/linuxbrew/). Enable EPEL and make sure compilers are installed.
```
sudo yum install epel-release gcc gcc-c++ automake ruby-irb
```
Install [LinuxBrew](http://brew.sh/linuxbrew/).
```
brew tap homebrew/science
ln -s $(which gcc) ~/.linuxbrew/bin/gcc-4.4
ln -s $(which g++) ~/.linuxbrew/bin/g++-4.4
brew install snp-sites
```

## Linux - using LinuxBrew
Install [LinuxBrew](http://brew.sh/linuxbrew/).
```
brew tap homebrew/science
brew install snp-sites
```

## OSX - using HomeBrew
Install [HomeBrew](http://brew.sh/). It requires a minimum of Xcode 5.1.1 (xcodebuild -version). Then run:
```
brew tap homebrew/science
brew install snp-sites
```
## OSX/Linux - from source
This is a difficult method and is only suitable for someone with advanced unix skills. No support is provided with this method, since you have advanced unix skills. Please consider using HomeBrew/LinuxBrew instead. First install a standard development environment (e.g. gcc, automake, autoconf, libtool). Download the software from [GitHub](https://github.com/sanger-pathogens/snp-sites).

```
autoreconf -i -f
./configure
make
sudo make install
```

## OSX/Linux - from a release tarball
This is a difficult method and is only suitable for someone with advanced unix skills. No support is provided with this method, since you have advanced unix skills. Please consider using HomeBrew/LinuxBrew instead. First install a standard development environment (e.g. gcc, automake, autoconf, libtool).

```
tar xzvf snp-sites-x.y.z.tar.gz
cd snp-sites-x.y.z
./configure
make
sudo make install
```

## Windows/OSX/Linux - using a Virtual Machine
A virtual machine (VM) is available containing the software. More details can be found on [the Pathogen VM page](http://sanger-pathogens.github.io/pathogens-vm/).
