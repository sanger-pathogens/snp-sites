#
#  From this base-image / starting-point
#
FROM debian:testing

#
#  Authorship
#
MAINTAINER path-help@sanger.ac.uk

#
# Pull in packages from testing
#
RUN apt-get update -qq

#
# Install SNP-sites
#
RUN apt-get -y install snp-sites
