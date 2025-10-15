#!/bin/bash

# Function to check if XQuartz is installed
function check_xquartz_installed {
    if [ ! -d "/Applications/Utilities/XQuartz.app" ]; then
        echo "XQuartz is not installed. Please download and install XQuartz from https://www.xquartz.org"
        exit 1
    fi
}

# Function to get the IP address of the Mac
function get_mac_ip {
    IP=$(ipconfig getifaddr en0)
    if [ -z "$IP" ]; then
        IP=$(ipconfig getifaddr en1)
    fi
    if [ -z "$IP" ]; then
        echo "Unable to determine IP address. Ensure you're connected to a network."
        exit 1
    fi
    echo "Mac IP address: $IP"
}

# Function to allow XQuartz connections
function configure_xquartz {
    echo "Configuring XQuartz to allow connections..."
    open -a XQuartz
    sleep 2
    xhost + $IP
}

# Function to run the Docker container
function run_docker_container {
    echo "Running Docker container with X11 forwarding..."
    docker start seatree.dev && docker exec -e DISPLAY=$IP:0 -it seatree.dev bash && docker exec -it seatree.dev bash -c "su - st -c 'export DISPLAY=$IP:0'"
}

# Main script starts here
check_xquartz_installed
get_mac_ip
configure_xquartz
run_docker_container
