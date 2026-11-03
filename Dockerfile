# SEATREE container image
#
# Ubuntu 24.04 with the system GMT 6, GTK 4 python bindings, and the
# HC and ConMan codes built from the git submodules. The GUI needs an
# X display; run the container with X forwarding (for example
# docker run -e DISPLAY=$DISPLAY -v /tmp/.X11-unix:/tmp/.X11-unix ...)
# or use it for the command line tools only.
#
# This file mirrors the steps in README.md. It has not been tested on
# every Docker setup; please report problems.
FROM ubuntu:24.04

ENV DEBIAN_FRONTEND=noninteractive

RUN apt-get update && apt-get install -y --no-install-recommends \
    git make gcc gfortran gawk \
    python3 python3-numpy python3-scipy python3-matplotlib \
    python3-gi python3-gi-cairo gir1.2-gtk-4.0 \
    gmt gmt-dcw gmt-gshhg libgmt-dev \
    libnetcdf-dev netcdf-bin \
    mpich libnetcdff-dev \
    ghostscript imagemagick \
    ca-certificates \
    && rm -rf /var/lib/apt/lists/*

# the submodules (HC, ConMan, EQdyna) are needed for those modules
RUN git clone --recurse-submodules https://github.com/thwbecker/seatree.git /opt/seatree

WORKDIR /opt/seatree

ENV SEATREE=/opt/seatree

# build HC and ConMan, then run the python installer non-interactively
# (it builds the seismology modules and writes the configuration files)
RUN make all \
    && make optional \
    && make test \
    && yes '' | ./install/configure.python3.gtk4

CMD ["./seatree"]
