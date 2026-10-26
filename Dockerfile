# Wanted to do the alpine version but that tag does not support arm64 (Mac M1) architecture
FROM julia:1.12
LABEL maintainer="Shaun Adkins (sadkins@som.umaryland.edu)"

WORKDIR /opt/FADU

COPY ./Project.toml ./Manifest.toml ./fadu.jl /opt/FADU/
COPY ./src /opt/FADU/src
# The tiny test dataset is used as the precompilation workload
COPY ./test /opt/FADU/test

# Install and precompile FADU and its dependencies into the image
ENV JULIA_DEPOT_PATH=/opt/julia_depot
ENV JULIA_PROJECT=/opt/FADU
RUN julia -e 'using Pkg; Pkg.instantiate(); Pkg.precompile()'

# Input and output files are expected to be mounted here
WORKDIR /data

ENTRYPOINT ["julia", "/opt/FADU/fadu.jl"]
