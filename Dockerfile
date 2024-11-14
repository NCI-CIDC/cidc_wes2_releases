# Base image 
FROM amd64/ubuntu:20.04
MAINTAINER Frankie Parks <frankie.parks@bioappdev.org>

## Run all root level commands 1st
# Create a user and group used to launch processes
RUN groupadd -r pipeline -g 1000 && useradd -u 1000 -r -g pipeline -m -d /home/pipeline -s /bin/bash -c "pipeline user" pipeline && chmod 755 /home/pipeline && chown pipeline:pipeline /home/pipeline

#Create mount point directory
RUN mkdir -p /media/analysis && chmod 777 /media/analysis

## Update to latest OS packages
RUN apt-get autoclean && apt-get update && apt-get -y dist-upgrade
RUN apt-get -y install wget vim curl bc less git

USER pipeline
WORKDIR /home/pipeline

## Download and install miniconda
RUN wget "https://github.com/conda-forge/miniforge/releases/download/23.11.0-0/Mambaforge-$(uname)-$(uname -m).sh" && bash Mambaforge-$(uname)-$(uname -m).sh -b -u
ENV PATH="/home/pipeline/mambaforge/bin:$PATH"

## Download and install gcloud utils
RUN curl -sSL https://sdk.cloud.google.com | bash
ENV PATH $PATH:/home/pipeline/google-cloud-sdk/bin


# Set the working directory to ubuntu' user home directory
RUN mkdir cidc_wes2
WORKDIR /home/pipeline/cidc_wes2

#Copy in all of the files/directories/etc.. from the project into the container
COPY . .
USER root
RUN chown -R pipeline:pipeline /home/pipeline/cidc_wes2
USER pipeline

## Install and initialize a mamba environment
#RUN conda install -y -n base -c conda-forge mamba
RUN mamba init

## Add conda channels and set priorities
RUN conda config --add channels defaults
RUN conda config --add channels bioconda
RUN conda config --add channels conda-forge
RUN conda config --set channel_priority flexible
RUN conda env create -f=environment.yaml
RUN conda create -n nci-cli -c bioconda python=3.11

#&& conda activate nci-cli && pip install nci-cidc-cli
#RUN conda activate cidc_wes2

#ENTRYPOINT ["conda", "activate", "cidc_wes2"]

CMD ["/bin/bash"]
