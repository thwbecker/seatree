#!/usr/bin/env python3
"""
SEATREE Installation Verification Script

This script verifies that all components of SEATREE are installed correctly:
- Environment variables
- NetCDF libraries and executables
- GMT executables
- ConMan executable
- Python dependencies
- Required paths and symlinks
"""

import os
import sys
import subprocess
from pathlib import Path


# Source the install script to load environment variables
def load_seatree_env():
    """Load SEATREE environment variables by sourcing install.seatree.sh"""
    script_dir = Path(__file__).parent
    install_script = script_dir / "install.seatree.sh"

    if not install_script.exists():
        print(f"Warning: {install_script} not found, skipping environment setup")
        return

    # Source the script and export variables to current process
    command = f'source {install_script} && env'
    try:
        proc = subprocess.Popen(command, stdout=subprocess.PIPE, shell=True,
                               executable='/bin/bash', cwd=script_dir)
        output, _ = proc.communicate()

        # Parse environment variables from output
        for line in output.decode().split('\n'):
            if '=' in line:
                key, _, value = line.partition('=')
                # Only update SEATREE-related variables
                if key in ['SEATREE', 'GMT4HOME', 'GMTHOME', 'NETCDFHOME',
                          'GMT_GSHHG_DATA', 'ARCH', 'PATH', 'LD_LIBRARY_PATH']:
                    os.environ[key] = value
    except Exception as e:
        print(f"Warning: Failed to source environment: {e}")


class Colors:
    """ANSI color codes for terminal output"""
    GREEN = '\033[92m'
    RED = '\033[91m'
    YELLOW = '\033[93m'
    BLUE = '\033[94m'
    BOLD = '\033[1m'
    END = '\033[0m'


class VerificationResult:
    """Store verification results"""
    def __init__(self):
        self.passed = []
        self.failed = []
        self.warnings = []

    def add_pass(self, check_name, details=""):
        self.passed.append((check_name, details))

    def add_fail(self, check_name, details=""):
        self.failed.append((check_name, details))

    def add_warning(self, check_name, details=""):
        self.warnings.append((check_name, details))

    def is_success(self):
        return len(self.failed) == 0


def print_header(text):
    """Print a section header"""
    print(f"\n{Colors.BOLD}{Colors.BLUE}{'='*60}{Colors.END}")
    print(f"{Colors.BOLD}{Colors.BLUE}{text}{Colors.END}")
    print(f"{Colors.BOLD}{Colors.BLUE}{'='*60}{Colors.END}\n")


def print_check(status, message, details=""):
    """Print a check result"""
    if status == "PASS":
        symbol = f"{Colors.GREEN}✓{Colors.END}"
        status_text = f"{Colors.GREEN}PASS{Colors.END}"
    elif status == "FAIL":
        symbol = f"{Colors.RED}✗{Colors.END}"
        status_text = f"{Colors.RED}FAIL{Colors.END}"
    elif status == "WARN":
        symbol = f"{Colors.YELLOW}⚠{Colors.END}"
        status_text = f"{Colors.YELLOW}WARN{Colors.END}"
    else:
        symbol = "?"
        status_text = status

    print(f"{symbol} [{status_text}] {message}")
    if details:
        print(f"         {details}")


def check_env_var(var_name, result):
    """Check if an environment variable is set and valid"""
    value = os.environ.get(var_name)
    if not value:
        result.add_fail(f"Environment variable {var_name}", "Not set")
        print_check("FAIL", f"Environment variable: {var_name}", "Not set")
        return None

    path = Path(value)
    if not path.exists():
        result.add_fail(f"Environment variable {var_name}", f"Path does not exist: {value}")
        print_check("FAIL", f"Environment variable: {var_name}", f"Path does not exist: {value}")
        return None

    result.add_pass(f"Environment variable {var_name}", value)
    print_check("PASS", f"Environment variable: {var_name}", value)
    return value


def check_executable(exe_path, name, result):
    """Check if an executable exists and is executable"""
    path = Path(exe_path)
    if not path.exists():
        result.add_fail(f"Executable {name}", f"Not found: {exe_path}")
        print_check("FAIL", f"Executable: {name}", f"Not found: {exe_path}")
        return False

    if not os.access(path, os.X_OK):
        result.add_fail(f"Executable {name}", f"Not executable: {exe_path}")
        print_check("FAIL", f"Executable: {name}", f"Not executable: {exe_path}")
        return False

    result.add_pass(f"Executable {name}", exe_path)
    print_check("PASS", f"Executable: {name}", exe_path)
    return True


def check_library(lib_path, name, result):
    """Check if a library exists"""
    path = Path(lib_path)
    if not path.exists():
        result.add_fail(f"Library {name}", f"Not found: {lib_path}")
        print_check("FAIL", f"Library: {name}", f"Not found: {lib_path}")
        return False

    result.add_pass(f"Library {name}", lib_path)
    print_check("PASS", f"Library: {name}", lib_path)
    return True


def check_python_import(module_name, result):
    """Check if a Python module can be imported"""
    try:
        __import__(module_name)
        result.add_pass(f"Python module {module_name}", "Successfully imported")
        print_check("PASS", f"Python module: {module_name}", "Successfully imported")
        return True
    except ImportError as e:
        result.add_fail(f"Python module {module_name}", str(e))
        print_check("FAIL", f"Python module: {module_name}", str(e))
        return False


def check_command_version(command, name, result):
    """Check if a command exists and can report version"""
    try:
        output = subprocess.check_output([command, '--version'],
                                        stderr=subprocess.STDOUT,
                                        timeout=5,
                                        universal_newlines=True)
        version = output.strip().split('\n')[0][:80]  # First line, max 80 chars
        result.add_pass(f"Command {name}", version)
        print_check("PASS", f"Command: {name}", version)
        return True
    except (subprocess.CalledProcessError, FileNotFoundError, subprocess.TimeoutExpired) as e:
        result.add_fail(f"Command {name}", f"Failed to execute: {command}")
        print_check("FAIL", f"Command: {name}", f"Failed to execute: {command}")
        return False


def main():
    """Main verification routine"""
    # Load environment variables first
    load_seatree_env()

    print_header("SEATREE Installation Verification")

    result = VerificationResult()

    # Check environment variables
    print_header("Environment Variables")
    seatree_home = check_env_var("SEATREE", result)
    gmt4_home = check_env_var("GMT4HOME", result)
    gmthome = check_env_var("GMTHOME", result)
    netcdf_home = check_env_var("NETCDFHOME", result)

    # Check GMT_GSHHG_DATA (optional but recommended)
    gshhg_data = os.environ.get("GMT_GSHHG_DATA")
    if gshhg_data and Path(gshhg_data).exists():
        result.add_pass("Environment variable GMT_GSHHG_DATA", gshhg_data)
        print_check("PASS", "Environment variable: GMT_GSHHG_DATA", gshhg_data)
    else:
        result.add_warning("Environment variable GMT_GSHHG_DATA", "Not set or path doesn't exist")
        print_check("WARN", "Environment variable: GMT_GSHHG_DATA", "Not set or path doesn't exist")

    # Check ARCH
    arch = os.environ.get("ARCH")
    if arch:
        result.add_pass("Environment variable ARCH", arch)
        print_check("PASS", "Environment variable: ARCH", arch)
    else:
        result.add_warning("Environment variable ARCH", "Not set")
        print_check("WARN", "Environment variable: ARCH", "Not set")

    # Check NetCDF installation
    if netcdf_home:
        print_header("NetCDF Installation")
        check_library(f"{netcdf_home}/lib/libnetcdf.so", "libnetcdf.so", result)
        check_library(f"{netcdf_home}/lib/libnetcdf.so.22", "libnetcdf.so.22", result)
        check_library(f"{netcdf_home}/lib/libnetcdf.so.7", "libnetcdf.so.7 (symlink)", result)
        check_executable(f"{netcdf_home}/bin/nc-config", "nc-config", result)
        check_executable(f"{netcdf_home}/bin/ncdump", "ncdump", result)
        check_command_version(f"{netcdf_home}/bin/nc-config", "nc-config --version", result)

    # Check GMT installation
    if gmt4_home:
        print_header("GMT 4 Installation")
        # Note: GMT4 does not have unified 'gmt' command, only individual programs
        gmt_executables = ['pscoast', 'psxy', 'psbasemap', 'pstext', 'gmtset']
        for exe in gmt_executables:
            check_executable(f"{gmt4_home}/bin/{exe}", exe, result)

    # Check ConMan installation
    if seatree_home:
        print_header("ConMan Installation")
        # Check both possible ConMan locations
        conman_new = f"{seatree_home}/modules/mc/ConMan/conman"
        conman_legacy = f"{seatree_home}/modules/mc/ConMan.legacy/conman"

        if Path(conman_new).exists():
            check_executable(conman_new, "conman (new)", result)
        elif Path(conman_legacy).exists():
            check_executable(conman_legacy, "conman (legacy)", result)
        else:
            result.add_fail("Executable conman", "Not found in ConMan/ or ConMan.legacy/")
            print_check("FAIL", "Executable: conman", "Not found in ConMan/ or ConMan.legacy/")

    # Check HC (mantle convection) installation
    if seatree_home:
        print_header("HC (Mantle Convection) Installation")
        arch = os.environ.get("ARCH", "")

        # HC uses bin/ARCH directory when ARCH is set
        hc_bin_arch = f"{seatree_home}/modules/mc/hc/bin/{arch}" if arch else None
        hc_bin = f"{seatree_home}/modules/mc/hc/bin"

        hc_executables = ['hc', 'sh_ana', 'sh_syn', 'hc_extract_sh_layer']
        for exe in hc_executables:
            # Try ARCH-specific path first, then fall back to bin/
            exe_path_arch = f"{hc_bin_arch}/{exe}" if hc_bin_arch else None
            exe_path = f"{hc_bin}/{exe}"

            if exe_path_arch and Path(exe_path_arch).exists():
                check_executable(exe_path_arch, f"hc/{exe}", result)
            elif Path(exe_path).exists():
                check_executable(exe_path, f"hc/{exe}", result)
            else:
                result.add_warning(f"Executable hc/{exe}", f"Not found in bin/ or bin/{arch}/")
                print_check("WARN", f"Executable: hc/{exe}", f"Not found (optional)")

    # Check Seismo modules installation
    if seatree_home:
        print_header("Seismo Modules Installation")

        # Larry (local tomography)
        larry_invert = f"{seatree_home}/modules/seismo/larry/invert"
        if Path(larry_invert).exists():
            check_executable(larry_invert, "larry/invert", result)
        else:
            result.add_warning("Executable larry/invert", f"Not found (optional)")
            print_check("WARN", "Executable: larry/invert", "Not found (optional)")

        # Larry3D (3D tomography)
        larry3d_invert = f"{seatree_home}/modules/seismo/larry3d/invert"
        larry3d_extract = f"{seatree_home}/modules/seismo/larry3d/extract"
        if Path(larry3d_invert).exists():
            check_executable(larry3d_invert, "larry3d/invert", result)
        else:
            result.add_warning("Executable larry3d/invert", f"Not found (optional)")
            print_check("WARN", "Executable: larry3d/invert", "Not found (optional)")

        if Path(larry3d_extract).exists():
            check_executable(larry3d_extract, "larry3d/extract", result)
        else:
            result.add_warning("Executable larry3d/extract", f"Not found (optional)")
            print_check("WARN", "Executable: larry3d/extract", "Not found (optional)")

    # Check Python dependencies
    print_header("Python Dependencies")
    python_modules = ['numpy', 'matplotlib', 'gi']
    for module in python_modules:
        check_python_import(module, result)

    # Special check for GTK4
    try:
        import gi
        gi.require_version('Gtk', '4.0')
        from gi.repository import Gtk
        result.add_pass("GTK 4.0", "Successfully loaded")
        print_check("PASS", "GTK 4.0", "Successfully loaded")
    except (ImportError, ValueError) as e:
        result.add_fail("GTK 4.0", str(e))
        print_check("FAIL", "GTK 4.0", str(e))

    # Check PATH and LD_LIBRARY_PATH
    print_header("Environment Paths")
    path_env = os.environ.get("PATH", "")
    if netcdf_home and f"{netcdf_home}/bin" in path_env:
        result.add_pass("NetCDF in PATH", f"{netcdf_home}/bin")
        print_check("PASS", "NetCDF in PATH", f"{netcdf_home}/bin")
    else:
        result.add_warning("NetCDF in PATH", "NetCDF bin directory not in PATH")
        print_check("WARN", "NetCDF in PATH", "NetCDF bin directory not in PATH")

    ld_library_path = os.environ.get("LD_LIBRARY_PATH", "")
    if netcdf_home and f"{netcdf_home}/lib" in ld_library_path:
        result.add_pass("NetCDF in LD_LIBRARY_PATH", f"{netcdf_home}/lib")
        print_check("PASS", "NetCDF in LD_LIBRARY_PATH", f"{netcdf_home}/lib")
    else:
        result.add_warning("NetCDF in LD_LIBRARY_PATH", "NetCDF lib directory not in LD_LIBRARY_PATH")
        print_check("WARN", "NetCDF in LD_LIBRARY_PATH", "NetCDF lib directory not in LD_LIBRARY_PATH")

    # Print summary
    print_header("Verification Summary")

    total_checks = len(result.passed) + len(result.failed) + len(result.warnings)

    print(f"{Colors.GREEN}Passed:   {len(result.passed):3d}{Colors.END}")
    print(f"{Colors.RED}Failed:   {len(result.failed):3d}{Colors.END}")
    print(f"{Colors.YELLOW}Warnings: {len(result.warnings):3d}{Colors.END}")
    print(f"Total:    {total_checks:3d}")

    print()

    if result.is_success():
        print(f"{Colors.BOLD}{Colors.GREEN}{'='*60}{Colors.END}")
        print(f"{Colors.BOLD}{Colors.GREEN}✓ VERIFICATION SUCCESSFUL{Colors.END}")
        print(f"{Colors.BOLD}{Colors.GREEN}{'='*60}{Colors.END}")
        print(f"\n{Colors.GREEN}All critical components are properly installed!{Colors.END}")
        if result.warnings:
            print(f"{Colors.YELLOW}There are {len(result.warnings)} warnings that may need attention.{Colors.END}")
        return 0
    else:
        print(f"{Colors.BOLD}{Colors.RED}{'='*60}{Colors.END}")
        print(f"{Colors.BOLD}{Colors.RED}✗ VERIFICATION FAILED{Colors.END}")
        print(f"{Colors.BOLD}{Colors.RED}{'='*60}{Colors.END}")
        print(f"\n{Colors.RED}Installation is incomplete or has errors.{Colors.END}")
        print(f"\n{Colors.RED}Failed checks:{Colors.END}")
        for check_name, details in result.failed:
            print(f"  • {check_name}")
            if details:
                print(f"    {details}")
        print(f"\n{Colors.YELLOW}Please review the log file and fix the issues above.{Colors.END}")
        return 1


if __name__ == "__main__":
    sys.exit(main())
