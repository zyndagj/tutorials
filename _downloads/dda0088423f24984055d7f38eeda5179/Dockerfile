FROM nvcr.io/nvidia/clara/clara-parabricks:4.3.1-1

ENV DEBIAN_FRONTEND=noninteractive 

RUN echo "#!/bin/bash\n\
apt-get clean\n\
[ -x \"\$(command -v conda)\" ] && conda clean -tipsy\n\
for dir in /tmp/* /var/tmp/* /home/jupyter/{.ccache,.cache/pip,conda-bld,.conda} /root/* /root/\.[^\.]* /var/lib/apt/lists/* /var/log/*; do\n\
	[ -e \$dir ] && rm -rf \$dir || true\n\
done" > /usr/bin/docker-clean && chmod a+rx /usr/bin/docker-clean && docker-clean

# Install dependencies
RUN apt-get update \
    && apt-get install -yq --no-install-recommends \
        unzip wget bzip2 ca-certificates git curl vim-nox less xterm \
    && docker-clean

# Add conda env variables
ENV MINICONDA_VERSION=py38_23.3.1-0 \
    CONDA_DIR=/opt/conda
ENV PATH=${CONDA_DIR}/bin:${PATH}
# Download and install miniconda
RUN mkdir $CONDA_DIR && chmod -R a+rX $CONDA_DIR \
    && wget --quiet https://repo.continuum.io/miniconda/Miniconda3-${MINICONDA_VERSION}-Linux-x86_64.sh \
    && bash Miniconda3-${MINICONDA_VERSION}-Linux-x86_64.sh -f -b -p $CONDA_DIR \
    && rm Miniconda3-${MINICONDA_VERSION}-Linux-x86_64.sh \
    && conda config --system --set auto_update_conda false \
    && conda config --system --set show_channel_urls true \
    && conda update -n base conda \
    && conda update --all --quiet --yes \
    && rm -rf ${CONDA_DIR}/pkgs/* \
    && docker-clean
# Activate conda on login
RUN ln -s /opt/conda/etc/profile.d/conda.sh /etc/profile.d/conda.sh

#https://docs.nvidia.com/clara/parabricks/4.2.0/documentation/tooldocs/compatiblecpusoftwareversions.html

RUN conda install -n base -c conda-forge mamba \
	&& docker-clean

RUN mamba install -y -n base -c conda-forge -c bioconda 'seqtk=1.3' 'samtools=1.16.1' 'gatk4=4.3.0.0' 'bwa=0.7.15' 'star=2.7.2a' 'star-fusion=1.7.0' \
    && docker-clean