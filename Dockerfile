################## BASE IMAGE ######################
FROM nfcore/base

################## METADATA ######################

LABEL base_image="nfcore/base"
LABEL version="1.0"
LABEL software="rnaseq-fusion-nf"
LABEL software.version="1.0"
LABEL about.summary="Container image containing all requirements for fusion gene discovery pipeline"
LABEL about.home="http://github.com/IARCbioinfo/rnaseq-fusion-nf"
LABEL about.documentation="http://github.com/IARCbioinfo/rnaseq-fusion-nf/README.md"
LABEL about.license_file="http://github.com/IARCbioinfo/rnaseq-fusion-nf/LICENSE.txt"
LABEL about.license="GNU-3.0"

################## MAINTAINER ######################
MAINTAINER **alcalan** <**alcalan@fellows.iarc.fr**>

################## INSTALLATION ######################
COPY environment.yml /
RUN conda env update -n root -f /environment.yml && conda clean -a
