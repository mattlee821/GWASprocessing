# GWASprocessing tools image: the environment in environment.yml plus bcftools 1.25 built with the BCFtools/liftover
# plugin, and a copy of the pipeline code. Built and published by .github/workflows/container.yml to
# ghcr.io/mattlee821/gwasprocessing. Nextflow runs every step inside it with -profile local,docker or slurm,apptainer
# (src/standardise.sh adds the profile when GWAS_CONTAINER is set). The references are not in the image: src/source.sh
# builds them on the host, running the build steps inside this image.
FROM mambaorg/micromamba:2.0.5

ARG BCFTOOLS_VERSION=1.25
ARG LIFTOVER_URL=https://raw.githubusercontent.com/freeseek/score/master/liftover.c

COPY --chown=$MAMBA_USER:$MAMBA_USER environment.yml /tmp/environment.yml
RUN micromamba install -y -n base -f /tmp/environment.yml && micromamba clean -a -y

# Build bcftools with the liftover plugin into the environment, with the environment's compiler and libraries.
ARG MAMBA_DOCKERFILE_ACTIVATE=1
RUN set -eux; \
    cd /tmp; \
    curl -fsSL -o bcftools.tar.bz2 \
      "https://github.com/samtools/bcftools/releases/download/${BCFTOOLS_VERSION}/bcftools-${BCFTOOLS_VERSION}.tar.bz2"; \
    tar -xjf bcftools.tar.bz2; \
    cd "bcftools-${BCFTOOLS_VERSION}"; \
    curl -fsSL -o plugins/liftover.c "${LIFTOVER_URL}"; \
    CC="${MAMBA_ROOT_PREFIX}/bin/cc" CPPFLAGS="-I${MAMBA_ROOT_PREFIX}/include" \
      LDFLAGS="-L${MAMBA_ROOT_PREFIX}/lib -Wl,-rpath,${MAMBA_ROOT_PREFIX}/lib" \
      ./configure --prefix="${MAMBA_ROOT_PREFIX}" --with-htslib=system --disable-libgsl --disable-perl-filters; \
    make -j"$(nproc)"; \
    make install; \
    cd /tmp; rm -rf "bcftools-${BCFTOOLS_VERSION}" bcftools.tar.bz2; \
    bcftools plugin -l | grep -qx liftover

USER root
COPY . /opt/GWASprocessing
RUN chmod -R a+rX /opt/GWASprocessing
USER $MAMBA_USER

# The scripts use the image's tools, not a host tool environment named in .env.
ENV GWAS_IN_CONTAINER=1 \
    PATH=/opt/conda/bin:$PATH \
    BCFTOOLS_PLUGINS=/opt/conda/libexec/bcftools
WORKDIR /opt/GWASprocessing
