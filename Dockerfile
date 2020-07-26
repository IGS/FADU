FROM julia:1.5-alpine
LABEL maintainer="Shaun Adkins (sadkins@som.umaryland.edu)"


RUN mkdir -p /opt/FADU

WORKDIR /opt/FADU

COPY ./alignment_overlaps.jl /opt/FADU
COPY ./bam_record.jl /opt/FADU/
COPY ./fadu.jl /opt/FADU/
COPY ./feature_counts.jl /opt/FADU/
COPY ./gff_feature.jl /opt/FADU/

# force compilation of external, shared packages copied in from Manifest.toml./Project.toml
# For some reason these would not copy to the "root" home directory
RUN mkdir -p /tmp/.julia/environments/v1.5
COPY ./fadu/Manifest.toml /tmp/.julia/environments/v1.5/Manifest.toml
COPY ./fadu/Project.toml /tmp/.julia/environments/v1.5/Project.toml
ENV JULIA_DEPOT_PATH=/tmp/.julia
RUN julia -e 'using Pkg; Pkg.precompile(); Pkg.instantiate()' && julia -e 'using BGZFStreams, GenomicFeatures, GFF3, Indexes, XAM'


ENTRYPOINT ["julia", "/opt/FADU/fadu.jl"]