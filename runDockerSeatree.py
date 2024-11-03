import os
import subprocess
import sys

def run_command(command):
    try:
        return subprocess.check_output(command, shell=True, text=True).strip()
    except subprocess.CalledProcessError:
        return None

def check_x_server():
    x_display = os.getenv("DISPLAY")
    if x_display:
        print(f"X server found with DISPLAY: {x_display}")
        return x_display
    print("X server not found. Make sure X11 is running for GUI apps.")
    return None

def check_docker_image(image_name):
    result = run_command(f"docker images -q {image_name}")
    return bool(result)

def check_docker_container(container_name):
    result = run_command(f"docker ps -aq -f name={container_name}")
    return result if result else None

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

def main():
    container_name = "seatree"
    image_name = "dunyuliu/seatree:v1.0.0"

    x_display = check_x_server()
    if not x_display:
        return 
    
    os.environ["DISPLAY"] = x_display
    xhost_ip = get_mac_ip()

    container_id = check_docker_container(container_name)
    print('Check if container exists; '+str(container_id))

    if container_id:
        print("Container "+container_name+" exists. Starting container ...")
        run_command(f"docker start {container_id}")
    else:
        print("Container "+container_name+" doesn't exist. Check image ...")
        if not check_docker_image(image_name):

            print("Docker image "+image_name+ " not found. Pulling it from Docker Hub ...")
            run_command(f"docker pull {image_name}") 
        else:
            print("Docker image "+image_name+ " found.")
        print("Creating Docker container "+container_name+ " ...")
        run_command(f"docker create --name {container_name} {image_name} tail -f /dev/null")
    
    print("Running Docker container with X11 forwarding...")
    subprocess.call(["docker", "start", container_name])
    subprocess.call(["docker", "exec", "-e", f"DISPLAY={xhost_ip}:0", "-it", container_name, "bash", "-c", f"python3 run.seatree.python3.gtk4.py"])

    # Check if XQuartz is installed
    if not os.path.isdir("/Applications/Utilities/XQuartz.app"):
        print("XQuartz is not installed. Please download and install XQuartz from https://www.xquartz.org")
        sys.exit(1)

if __name__ == "__main__":
    main()

