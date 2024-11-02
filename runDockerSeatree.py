import os
import subprocess
import sys

def main():
    if len(sys.argv) != 2:
        print("Usage: python3 script.py CONTAINER_NAME")
        sys.exit(1)
    
    container_name = sys.argv[1]

    # Check if XQuartz is installed
    if not os.path.isdir("/Applications/Utilities/XQuartz.app"):
        print("XQuartz is not installed. Please download and install XQuartz from https://www.xquartz.org")
        sys.exit(1)

    # Get the IP address of the Mac
    ip = get_mac_ip()
    if not ip:
        print("Unable to determine IP address. Ensure you're connected to a network.")
        sys.exit(1)
    print(f"Mac IP address: {ip}")

    # Configure XQuartz to allow connections
    configure_xquartz(ip)

    # Run the Docker container with X11 forwarding
    run_docker_container(container_name, ip)

def get_mac_ip():
    try:
        ip = subprocess.check_output(["ipconfig", "getifaddr", "en0"]).strip().decode('utf-8')
    except subprocess.CalledProcessError:
        try:
            ip = subprocess.check_output(["ipconfig", "getifaddr", "en1"]).strip().decode('utf-8')
        except subprocess.CalledProcessError:
            return None
    return ip

def configure_xquartz(ip):
    print("Configuring XQuartz to allow connections...")
    subprocess.call(["open", "-a", "XQuartz"])
    sleep(2)
    subprocess.call(["xhost", "+", ip])

def run_docker_container(container_name, ip):
    print("Running Docker container with X11 forwarding...")
    subprocess.call(["docker", "start", container_name])
    subprocess.call(["docker", "exec", "-e", f"DISPLAY={ip}:0", "-it", container_name, "bash", "-c", f"export DISPLAY={ip}:0; exec bash"])
    subprocess.call(["docker", "exec", "-it", container_name, "bash", "-c", f"su - st -c 'export DISPLAY={ip}:0'"])

def sleep(seconds):
    from time import sleep as real_sleep
    real_sleep(seconds)

if __name__ == "__main__":
    main()

