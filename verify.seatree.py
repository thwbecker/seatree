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
import shutil
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
                          'GMT_GSHHG_DATA', 'GMTVERSION', 'ARCH', 'PATH', 'LD_LIBRARY_PATH']:
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


def check_env_var(var_name, result, required=True):
    """Check if an environment variable is set and valid"""
    value = os.environ.get(var_name)
    if not value:
        status = "FAIL" if required else "WARN"
        message = "Not set"
        if required:
            result.add_fail(f"Environment variable {var_name}", message)
        else:
            result.add_warning(f"Environment variable {var_name}", message)
        print_check(status, f"Environment variable: {var_name}", message)
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


def detect_netcdf_paths():
    """Discover NetCDF installation paths using env vars and nc-config"""
    info = {
        'home': None,
        'c_home': None,
        'libdir': None,
        'bindir': None,
        'nc_config': None,
        'ncdump': None,
        'lib_candidates': [],
    }

    env_home = os.environ.get("NETCDFHOME")
    env_c_home = os.environ.get("NETCDF_C_HOME")

    if env_home:
        home_path = Path(env_home)
        if home_path.exists():
            info['home'] = home_path
    if env_c_home:
        c_home_path = Path(env_c_home)
        if c_home_path.exists():
            info['c_home'] = c_home_path

    def add_lib_candidate(path, prefer=False):
        if not path or not path.exists():
            return
        if info['libdir'] is None or prefer:
            info['libdir'] = path
        if path not in info['lib_candidates']:
            info['lib_candidates'].append(path)

    def consider_prefix(prefix_path, prefer=False):
        if not prefix_path or not prefix_path.exists():
            return
        lib_candidate = prefix_path / "lib"
        add_lib_candidate(lib_candidate, prefer=prefer)
        bin_candidate = prefix_path / "bin"
        if bin_candidate.exists():
            if prefer or info['bindir'] is None:
                info['bindir'] = bin_candidate
            nc_candidate = bin_candidate / "nc-config"
            if nc_candidate.exists():
                if prefer or info['nc_config'] is None:
                    info['nc_config'] = nc_candidate
            ncd_candidate = bin_candidate / "ncdump"
            if ncd_candidate.exists():
                if prefer or info['ncdump'] is None:
                    info['ncdump'] = ncd_candidate

    consider_prefix(info.get('c_home'), prefer=True)
    consider_prefix(info.get('home'))

    if info['nc_config'] is None:
        which_nc = shutil.which("nc-config")
        if which_nc:
            info['nc_config'] = Path(which_nc)
            if info['bindir'] is None:
                info['bindir'] = Path(which_nc).parent

    if info['nc_config']:
        try:
            prefix_out = subprocess.check_output(
                [str(info['nc_config']), '--prefix'],
                stderr=subprocess.STDOUT,
                timeout=5,
                universal_newlines=True,
            ).strip().split('\n')[0]
            prefix_path = Path(prefix_out)
            if prefix_path.exists():
                if info['c_home'] is None:
                    info['c_home'] = prefix_path
                consider_prefix(prefix_path, prefer=True)
        except (subprocess.CalledProcessError, FileNotFoundError, subprocess.TimeoutExpired):
            pass

        try:
            libdir_out = subprocess.check_output(
                [str(info['nc_config']), '--libdir'],
                stderr=subprocess.STDOUT,
                timeout=5,
                universal_newlines=True,
            ).strip().split('\n')[0]
            libdir_path = Path(libdir_out)
            add_lib_candidate(libdir_path, prefer=True)
        except (subprocess.CalledProcessError, FileNotFoundError, subprocess.TimeoutExpired):
            pass

        try:
            bindir_out = subprocess.check_output(
                [str(info['nc_config']), '--bindir'],
                stderr=subprocess.STDOUT,
                timeout=5,
                universal_newlines=True,
            ).strip().split('\n')[0]
            bindir_path = Path(bindir_out)
            if bindir_path.exists():
                info['bindir'] = bindir_path
                ncd_candidate = bindir_path / "ncdump"
                if ncd_candidate.exists():
                    info['ncdump'] = ncd_candidate
        except (subprocess.CalledProcessError, FileNotFoundError, subprocess.TimeoutExpired):
            pass

        try:
            libs_out = subprocess.check_output(
                [str(info['nc_config']), '--libs'],
                stderr=subprocess.STDOUT,
                timeout=5,
                universal_newlines=True,
            )
            for token in libs_out.split():
                if token.startswith('-L'):
                    lib_path = Path(token[2:])
                    add_lib_candidate(lib_path)
        except (subprocess.CalledProcessError, FileNotFoundError, subprocess.TimeoutExpired):
            pass

    ld_library_path = os.environ.get("LD_LIBRARY_PATH", "")
    for entry in ld_library_path.split(os.pathsep):
        entry = entry.strip()
        if entry:
            add_lib_candidate(Path(entry))

    # Add a few common system locations as fallbacks
    for default_dir in [
        Path("/usr/lib/x86_64-linux-gnu"),
        Path("/usr/lib64"),
        Path("/usr/local/lib"),
        Path("/usr/lib"),
    ]:
        add_lib_candidate(default_dir)

    if info['ncdump'] is None:
        which_ncdump = shutil.which("ncdump")
        if which_ncdump:
            info['ncdump'] = Path(which_ncdump)

    return info


def locate_library(preferred_names, candidates):
    """Find the first existing library file in candidate directories."""
    seen = set()
    stem = preferred_names[0]
    for candidate in candidates:
        if not candidate:
            continue
        candidate = Path(candidate)
        if candidate in seen or not candidate.exists():
            continue
        seen.add(candidate)

        for name in preferred_names:
            lib_path = candidate / name
            if lib_path.exists():
                return lib_path

        for match in sorted(candidate.glob(f"{stem}*")):
            if match.is_file():
                return match
    return None


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
    netcdf_home = check_env_var("NETCDFHOME", result, required=False)

    # Check GMT_GSHHG_DATA (only needed for GMT4)
    gmtversion = os.environ.get("GMTVERSION", "6")
    gshhg_data = os.environ.get("GMT_GSHHG_DATA")
    if gmtversion == "4":
        if gshhg_data and Path(gshhg_data).exists():
            result.add_pass("Environment variable GMT_GSHHG_DATA", gshhg_data)
            print_check("PASS", "Environment variable: GMT_GSHHG_DATA", gshhg_data)
        else:
            result.add_warning("Environment variable GMT_GSHHG_DATA", "Not set or path doesn't exist (required for GMT4)")
            print_check("WARN", "Environment variable: GMT_GSHHG_DATA", "Not set or path doesn't exist (required for GMT4)")
    else:
        # GMT6 doesn't need GMT_GSHHG_DATA
        if gshhg_data:
            result.add_pass("Environment variable GMT_GSHHG_DATA", f"{gshhg_data} (GMT6 mode - not required)")
            print_check("PASS", "Environment variable: GMT_GSHHG_DATA", f"{gshhg_data} (GMT6 mode - not required)")
        else:
            result.add_pass("Environment variable GMT_GSHHG_DATA", "Not set (GMT6 mode - not required)")
            print_check("PASS", "Environment variable: GMT_GSHHG_DATA", "Not set (GMT6 mode - not required)")

    # Check ARCH
    arch = os.environ.get("ARCH")
    if arch:
        result.add_pass("Environment variable ARCH", arch)
        print_check("PASS", "Environment variable: ARCH", arch)
    else:
        result.add_warning("Environment variable ARCH", "Not set")
        print_check("WARN", "Environment variable: ARCH", "Not set")

    # Discover NetCDF paths
    netcdf_info = detect_netcdf_paths()

    # Check NetCDF installation
    print_header("NetCDF Installation")
    is_macos = sys.platform == 'darwin'

    lib_candidates = netcdf_info.get('lib_candidates', [])
    netcdf_lib_path = None
    netcdf_libdir = netcdf_info.get('libdir')
    nc_config_path = netcdf_info.get('nc_config')
    ncdump_path = netcdf_info.get('ncdump')

    lib_filename = "libnetcdf.dylib" if is_macos else "libnetcdf.so"
    netcdf_lib_path = locate_library([lib_filename], lib_candidates or ([netcdf_libdir] if netcdf_libdir else []))
    if netcdf_lib_path:
        netcdf_libdir = netcdf_lib_path.parent
        check_library(str(netcdf_lib_path), lib_filename, result)
    else:
        message = "Unable to locate NetCDF library directory (set NETCDFHOME or ensure nc-config is in PATH)"
        result.add_fail("NetCDF libraries", message)
        print_check("FAIL", "NetCDF libraries", message)

    if is_macos:
        fortran_lib_path = None
        if netcdf_info.get('home'):
            candidate = netcdf_info['home'] / "lib" / "libnetcdff.dylib"
            if candidate.exists():
                fortran_lib_path = candidate
        if not fortran_lib_path:
            fortran_lib_path = locate_library(["libnetcdff.dylib"], lib_candidates)
        if fortran_lib_path:
            check_library(str(fortran_lib_path), "libnetcdff.dylib (fortran)", result)
        else:
            result.add_fail("Library libnetcdff.dylib", "Not found; install netcdf-fortran or set NETCDFHOME")
            print_check("FAIL", "Library: libnetcdff.dylib", "Not found; install netcdf-fortran or set NETCDFHOME")

    if nc_config_path:
        check_executable(str(nc_config_path), "nc-config", result)
        check_command_version(str(nc_config_path), "nc-config --version", result)
    else:
        message = "nc-config not found; ensure NetCDF binaries are installed and in PATH"
        result.add_fail("Executable nc-config", message)
        print_check("FAIL", "Executable: nc-config", message)

    if ncdump_path:
        check_executable(str(ncdump_path), "ncdump", result)
    else:
        message = "ncdump not found; ensure NetCDF utilities are installed and in PATH"
        result.add_fail("Executable ncdump", message)
        print_check("FAIL", "Executable: ncdump", message)

    # Check GMT installation
    gmtversion = os.environ.get("GMTVERSION", "6")
    if gmtversion == "4":
        if gmt4_home:
            print_header("GMT 4 Installation")
            # Note: GMT4 does not have unified 'gmt' command, only individual programs
            gmt_executables = ['pscoast', 'psxy', 'psbasemap', 'pstext', 'gmtset']
            for exe in gmt_executables:
                check_executable(f"{gmt4_home}/bin/{exe}", exe, result)
    else:
        # GMT6 mode - check for system GMT
        print_header("GMT 6 Installation")
        if gmthome:
            # Check GMTHOME path exists
            gmthome_path = Path(gmthome)
            if gmthome_path.exists():
                result.add_pass("GMTHOME path", gmthome)
                print_check("PASS", "GMTHOME path", gmthome)
            else:
                result.add_fail("GMTHOME path", f"Path does not exist: {gmthome}")
                print_check("FAIL", "GMTHOME path", f"Path does not exist: {gmthome}")

        # Check for GMT6 executable in system PATH
        check_command_version("gmt", "GMT", result)

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

        # EQdyna (earthquake dynamics)
        eqdyna_bin = f"{seatree_home}/modules/seismo/EQdyna/bin/eqdyna"
        eqdyna_src = f"{seatree_home}/modules/seismo/EQdyna/src/eqdyna"
        if Path(eqdyna_bin).exists():
            check_executable(eqdyna_bin, "EQdyna/eqdyna", result)
        elif Path(eqdyna_src).exists():
            check_executable(eqdyna_src, "EQdyna/eqdyna (src)", result)
        else:
            result.add_warning("Executable EQdyna/eqdyna", f"Not found (optional)")
            print_check("WARN", "Executable: EQdyna/eqdyna", "Not found (optional)")

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

    # Check LD_LIBRARY_PATH (critical for runtime linking)
    print_header("Environment Paths")
    ld_library_path = os.environ.get("LD_LIBRARY_PATH", "")
    netcdf_libdir_str = str(netcdf_libdir) if netcdf_libdir else None
    ld_entries = [p for p in ld_library_path.split(os.pathsep) if p]
    if netcdf_libdir_str:
        if netcdf_libdir_str in ld_entries:
            result.add_pass("NetCDF in LD_LIBRARY_PATH", netcdf_libdir_str)
            print_check("PASS", "NetCDF in LD_LIBRARY_PATH", netcdf_libdir_str)
        elif netcdf_libdir_str.startswith(("/usr/lib", "/usr/local/lib", "/lib")):
            msg = f"{netcdf_libdir_str} (system default search path)"
            result.add_pass("NetCDF in LD_LIBRARY_PATH", msg)
            print_check("PASS", "NetCDF in LD_LIBRARY_PATH", msg)
        else:
            result.add_warning("NetCDF in LD_LIBRARY_PATH", f"{netcdf_libdir_str} not found in LD_LIBRARY_PATH")
            print_check("WARN", "NetCDF in LD_LIBRARY_PATH", f"{netcdf_libdir_str} not found in LD_LIBRARY_PATH")
    else:
        result.add_warning("NetCDF in LD_LIBRARY_PATH", "NetCDF lib directory unknown; skipping LD_LIBRARY_PATH check")
        print_check("WARN", "NetCDF in LD_LIBRARY_PATH", "NetCDF lib directory unknown; skipping check")

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
