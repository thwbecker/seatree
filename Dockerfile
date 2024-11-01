# Use the official Ubuntu image as a base
FROM ubuntu:latest

# Prevent prompts during package installations
ENV DEBIAN_FRONTEND=noninteractive

# Update the package list and install required packages
RUN apt-get update && apt-get install -y \
    git \
    vim \
    wget \
    cmake \
    python3 \
    python3-pip \
    python3-matplotlib \
    gmt \
    python3-numpy \
    python3-gi \
    python3-gi-cairo \
    gir1.2-gtk-4.0 \
    libgtk-4-dev \
    build-essential \
    gfortran \
    x11-apps \
    libhdf5-dev

# Install the specific version of matplotlib using pip
RUN pip3 install --user matplotlib==3.9.2 --break-system-packages

# Clone the repository into /home
RUN git clone https://github.com/dunyuliu/seatree.dev.git /home/seatree

# Set working directory
WORKDIR /home/seatree

RUN bash install.gmt4.ubuntu22.sh

RUN bash install.netcdf.ubuntu22.sh

# Set environment variables
ENV SEATREE=/home/seatree \
    GMT4HOME=/home/seatree/gmt-4.5.18 \
    GMTHOME=/home/seatree/gmt-4.5.18 \
    GMT_GSHHG_DATA=/home/seatree/gmt-4.5.18/gshhg-gmt-2.3.7 \
    NETCDFHOME=/home/seatree/netcdf-c-4.9.3-rc1 \
    PATH=$PATH:/home/seatree/netcdf-c-4.9.3-rc1/lib:/home/seatree/gmt-4.5.18/gshhg-gmt-2.3.7 \
    LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/home/seatree/netcdf-c-4.9.3-rc1/lib

# Create symbolic link
RUN ln -s $NETCDFHOME/lib/libnetcdf.so.22 $NETCDFHOME/lib/libnetcdf.so.7



# Run the configuration script
RUN yes '' | ./configure.python3.gtk4

# Set the default command to run when starting the container
CMD ["bash", "-c", "cd /home/seatree && python3 run.seatree.python3.gtk4.py"]
