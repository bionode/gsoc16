# Set the base image to Ubuntu
FROM biodckr/biodocker

################## BWA ######################

ENV ZIP=bwakit-0.7.12_x64-linux.tar.bz2
ENV URL=https://github.com/BioDocker/software-archive/releases/download/bwa
ENV FOLDER=bwa.kit
ENV DST=/home/biodocker/bin

RUN wget $URL/$ZIP -O $DST/$ZIP && \
  tar xvf $DST/$ZIP --strip-components=1 -C $DST && \
  rm $DST/$ZIP

# CHANGE WORKDIR TO /DATA
# WORKDIR /data/

################## samtools ###########################

# Change user to root
USER root

RUN apt-get update && \
  apt-get install -y \
  libncurses5-dev \
  libncursesw5-dev && \
  apt-get clean && \
  apt-get purge && \
  rm -rf /var/lib/apt/lists/* /tmp/* /var/tmp/*

ENV ZIP=samtools-1.2.tar.bz2
ENV URL=https://github.com/BioDocker/software-archive/releases/download/samtools
ENV FOLDER=samtools-1.2
ENV DST=/tmp

RUN wget $URL/$ZIP -O $DST/$ZIP && \
  tar xvf $DST/$ZIP -C $DST && \
  rm $DST/$ZIP && \
  cd $DST/$FOLDER && \
  make && \
  make install && \
  cd / && \
  rm -rf $DST/$FOLDER

# Change user to back to biodocker
USER biodocker

################## bcftools ###########################

# Change user to root
USER root

ENV ZIP=bcftools-1.2.tar.bz2
ENV URL=https://github.com/BioDocker/software-archive/releases/download/bcftools
ENV FOLDER=bcftools-1.2
ENV DST=/tmp

RUN wget $URL/$ZIP -O $DST/$ZIP && \
  tar xvf $DST/$ZIP -C $DST && \
  rm $DST/$ZIP && \
  cd $DST/$FOLDER && \
  make && \
  make install && \
  cd / && \
  rm -rf $DST/$FOLDER

# Change user to back to biodocker
USER biodocker

# CHANGE WORKDIR TO /DATA
WORKDIR /data/
