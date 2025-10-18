import gi
gi.require_version('Gtk', '4.0')
from gi.repository import Gtk, GLib

from collections import OrderedDict
import os, sys, time, xml.dom.minidom, subprocess, shutil

from seatree.plotter.matPlotLib import matPlotLibPlotter

from seatree.modules.module import Module
from .eqdynaGUI import EQDYNAGUI
from seatree.util.scriptRunner import ScriptRunner

import seatree.plotter.matPlotLib.matPlotLibPlotter

class EQDYNA(Module):

    PLOT_FILE_SPECS = [
        ("Collage", "cRuptureDynamics.png"),
        ("Slip", "rupture_dynamics_slip.png"),
        ("Slip s", "rupture_dynamics_slip_s.png"),
        ("Slip d", "rupture_dynamics_slip_d.png"),
        ("Peak slip rate", "rupture_dynamics_peak_slip_rate.png"),
        ("Final slip rate", "rupture_dynamics_final_slip_rate.png"),
        ("Shear stress", "rupture_dynamics_shear_stress.png"),
        ("Dip shear stress", "rupture_dynamics_dip_shear_stress.png"),
        ("Normal stress", "rupture_dynamics_normal_stress.png"),
    ]

    def __init__(self):
        '''
        Constructor for EQDYNA module
        '''

        GLib.threads_init()
        # short name for the module
        shortName = "EQDYNA"

        # long, display name for the module
        longName = "Earthquake Dynamics (EQDYNA)"

        # version number
        version = 5.3

        # name of the directory that should be created inside of the users
        # home directory, inside of the .seatree folder
        storeName = "eqdyna"

        # base image
        baseImage = ""

        # call the Module constructor
        Module.__init__(self, shortName, longName, version, storeName, baseImage)

    def getPanel(self, mainWindow):
        '''
        Returns the GUI panel for the EQDYNA module
        '''
        self.gui = EQDYNAGUI(self)
        return self.gui

    def setDefaults(self, mainWindow):
        '''
        Initialize the module with default settings
        '''
        self.mainWindow = mainWindow
        self.plotter = matPlotLibPlotter.MatPlotLibPlotter(self, self.mainWindow, 300, 200, False)
        path = self.mainWindow.getPath()
        if not path.endswith(os.sep):
            path += os.sep
        self.loadConfFile()
        self.tempDir = self.mainWindow.getTempFileDir()
        if not self.tempDir.endswith(os.sep):
            self.tempDir += os.sep
        self.scriptRunner = ScriptRunner(workingDir=self.tempDir)
        self.commandString = ""  # Initialize command tracking
        self._log_handle = None  # Track streaming log handle
        self.generatedPlots = OrderedDict()

    def getPlotter(self):
        """
        Returns the plotter object for the module
        """
        return self.plotter

    def loadConfFile(self):
        '''
        Load configuration from eqdynaConf.xml
        '''
        confFile = os.path.join(self.seatreePath, "conf", "eqdyna", "eqdynaConf.xml")
        if os.path.exists(confFile):
            doc = xml.dom.minidom.parse(confFile)

            # load eqdyna path
            eqdynaNode = doc.getElementsByTagName("eqdynaPath")
            if eqdynaNode and eqdynaNode[0].firstChild:
                eqdynaPath = eqdynaNode[0].firstChild.nodeValue.strip()

                # Handle APPIMAGE_ROOT token replacement
                if eqdynaPath.startswith("APPIMAGE_ROOT"):
                    appimage_root = os.environ.get("APPIMAGE_ROOT", "")
                    if appimage_root:
                        eqdynaPath = eqdynaPath.replace("APPIMAGE_ROOT", appimage_root)

                # Convert relative paths to absolute paths
                if eqdynaPath and not os.path.isabs(eqdynaPath):
                    eqdynaPath = os.path.abspath(os.path.join(self.seatreePath, eqdynaPath))

                if not eqdynaPath:
                    eqdynaPath = ""
                elif not eqdynaPath.endswith(os.sep):
                    eqdynaPath = eqdynaPath + os.sep
            else:
                eqdynaPath = ""
        else:
            eqdynaPath = ""
        self.eqdynaPath = eqdynaPath
        print("EQDYNA path: " + self.eqdynaPath)

    def createCase(self, caseName, caseInputDir=None):
        '''
        Create a new EQDYNA case
        caseName: name of the case to create
        caseInputDir: path to directory containing user_defined_params.py (optional)
        '''
        if not self.eqdynaPath:
            print("ERROR: EQDYNA path not set!")
            return False

        # Run create.newcase script
        createScript = os.path.join(self.eqdynaPath, "scripts", "create.newcase")
        if not os.path.exists(createScript):
            print(f"ERROR: create.newcase script not found at {createScript}")
            return False

        # Create case directory path
        caseDir = os.path.join(self.tempDir, caseName)

        # Determine compset from case input directory
        compset = "test.tpv8"  # default
        if caseInputDir and os.path.exists(caseInputDir):
            # Get the basename of the case input directory (e.g., "test.tpv8")
            compset = os.path.basename(caseInputDir)
            print(f"Using compset: {compset} from {caseInputDir}")

        # Build command - create.newcase takes: <case_path> <compset>
        cmd = f"python3 {createScript} {caseDir} {compset}"

        self.commandString += f"# Creating EQDYNA case\n{cmd}\n"
        print(f"Creating case: {cmd}")

        # Set EQDYNAROOT environment variable required by create.newcase
        env = os.environ.copy()
        # Remove trailing separator for EQDYNAROOT
        eqdyna_root = self.eqdynaPath.rstrip(os.sep)
        env['EQDYNAROOT'] = eqdyna_root

        result = subprocess.run(cmd, shell=True, capture_output=True, text=True, env=env)
        if result.returncode != 0:
            print(f"ERROR creating case: {result.stderr}")
            print(f"stdout: {result.stdout}")
            return False

        print(f"Case created successfully at {caseDir}")
        print(result.stdout)

        self.currentCase = caseDir
        self.currentCaseInputDir = caseInputDir
        self.generatedPlots = OrderedDict()

        return True

    def setupCase(self):
        '''
        Setup the EQDYNA case
        '''
        if not hasattr(self, 'currentCase'):
            print("ERROR: No case created!")
            return False

        setupScript = os.path.join(self.currentCase, "case.setup")
        if not os.path.exists(setupScript):
            print(f"ERROR: case.setup script not found at {setupScript}")
            return False

        # Make sure script is executable
        os.chmod(setupScript, 0o755)

        cmd = f"python3 {setupScript}"
        self.commandString += f"# Setting up EQDYNA case\n{cmd}\n"
        print(f"Setting up case: {cmd}")

        result = subprocess.run(cmd, shell=True, cwd=self.currentCase, capture_output=True, text=True)
        if result.returncode != 0:
            print(f"ERROR setting up case: {result.stderr}")
            print(f"stdout: {result.stdout}")
            return False

        print("Case setup completed successfully")
        print(result.stdout)
        return True

    def runCase(self):
        '''
        Run the EQDYNA simulation
        '''
        if not hasattr(self, 'currentCase'):
            print("ERROR: No case created!")
            return False

        # Check if eqdyna executable exists
        eqdynaExe = os.path.join(self.eqdynaPath, "bin", "eqdyna")
        if not os.path.exists(eqdynaExe):
            print(f"ERROR: eqdyna executable not found at {eqdynaExe}")
            return False

        mpi_launcher = self._select_mpi_launcher()
        launcher_basename = os.path.basename(mpi_launcher)
        mpi_ranks = self._infer_mpi_ranks(self.currentCase)
        if mpi_ranks <= 0:
            mpi_ranks = 1

        # Build MPICH command directly (no run.sh)
        cmd = [mpi_launcher]
        if launcher_basename == "mpirun":
            cmd.append("--oversubscribe")
        cmd.extend(["-np", str(mpi_ranks), eqdynaExe])

        # Prepare log file (overwrite previous content)
        self.logFile = os.path.join(self.currentCase, "eqdyna.log")
        # Optional pre-clean step mirrors run.sh behaviour
        clean_script = os.path.join(self.currentCase, "clean.py")
        if os.path.exists(clean_script):
            clean_cmd = ["python3", "clean.py"]
            self.commandString += f"# Cleaning EQDYNA case\n{subprocess.list2cmdline(clean_cmd)}\n"
            with open(self.logFile, "w") as log_stream:
                log_stream.write("== Running clean.py ==\n")
                log_stream.flush()
                clean_result = subprocess.run(
                    clean_cmd,
                    cwd=self.currentCase,
                    stdout=log_stream,
                    stderr=subprocess.STDOUT,
                    text=True
                )
                log_stream.write(f"\n== clean.py exited with {clean_result.returncode} ==\n")
        else:
            with open(self.logFile, "w") as log_stream:
                log_stream.write("clean.py not found; skipping cleanup step.\n")

        command_str = subprocess.list2cmdline(cmd)
        self.commandString += f"# Running EQDYNA simulation\n{command_str}\n"

        print(f"Running simulation: {command_str}")
        print(f"Log file: {self.logFile}")

        # Launch EQDYNA asynchronously and stream output to log
        self._log_handle = open(self.logFile, "a")
        self._log_handle.write("\n== Launching EQDYNA ==\n")
        self._log_handle.flush()
        self.process = subprocess.Popen(
            cmd,
            cwd=self.currentCase,
            stdout=self._log_handle,
            stderr=subprocess.STDOUT,
            text=True
        )

        return True

    def _select_mpi_launcher(self):
        """
        Select the MPI launcher, preferring MPICH variants to match the bundled executable.
        """
        mpi_launcher = None
        for candidate in ("mpiexec.mpich", "mpirun.mpich", "mpiexec", "mpirun"):
            launcher_path = shutil.which(candidate)
            if launcher_path:
                mpi_launcher = launcher_path
                break
        if not mpi_launcher:
            print("WARNING: No MPI launcher found in PATH; defaulting to 'mpirun'")
            mpi_launcher = "mpirun"
        return mpi_launcher

    def _infer_mpi_ranks(self, case_dir):
        """
        Determine MPI ranks from bGlobal.txt (npx * npy * npz). Defaults to 1 on failure.
        """
        bglobal_path = os.path.join(case_dir, "bGlobal.txt")
        if not os.path.exists(bglobal_path):
            return 1

        try:
            with open(bglobal_path, "r") as f:
                nonempty = [line.strip() for line in f if line.strip()]
            if len(nonempty) <= 12:
                return 1
            parts = nonempty[12].split()
            if len(parts) < 3:
                return 1
            npx, npy, npz = (int(float(parts[i])) for i in range(3))
            total = npx * npy * npz
            return total if total > 0 else 1
        except Exception as exc:
            print(f"WARNING: Unable to infer MPI ranks from bGlobal.txt: {exc}")
            return 1

    def checkStatus(self):
        '''
        Check if simulation is still running
        '''
        if hasattr(self, 'process') and self.process:
            running = self.process.poll() is None
            if not running and self._log_handle:
                self._log_handle.close()
                self._log_handle = None
            return running
        return False

    def getLogFile(self):
        '''
        Get path to simulation log file
        '''
        if hasattr(self, 'logFile'):
            return self.logFile
        elif hasattr(self, 'currentCase'):
            return os.path.join(self.currentCase, "eqdyna.log")
        return None

    def readLog(self):
        '''
        Read the simulation log file
        '''
        logFile = self.getLogFile()
        if logFile and os.path.exists(logFile):
            with open(logFile, 'r') as f:
                return f.read()
        return "No log file available"

    def postProcess(self):
        '''
        Run post-processing scripts to generate plots
        '''
        if not hasattr(self, 'currentCase'):
            print("ERROR: No case to post-process!")
            return False

        # Check if simulation output files exist
        import glob
        frt_files = glob.glob(os.path.join(self.currentCase, "frt.txt*"))
        if not frt_files:
            print("ERROR: No simulation output files (frt.txt*) found!")
            print("Please make sure the simulation has completed successfully.")
            return False

        print(f"Found {len(frt_files)} output files to post-process")

        plotScript = os.path.join(self.currentCase, "plotRuptureDynamics")
        if not os.path.exists(plotScript):
            print(f"ERROR: plotRuptureDynamics script not found at {plotScript}")
            return False

        cmd = f"python3 {plotScript}"
        self.commandString += f"# Post-processing results\n{cmd}\n"
        print(f"Post-processing: {cmd}")

        result = subprocess.run(cmd, shell=True, cwd=self.currentCase, capture_output=True, text=True)
        if result.returncode != 0:
            print(f"ERROR post-processing: {result.stderr}")
            print(f"stdout: {result.stdout}")
            return False

        print("Post-processing completed successfully")
        print(result.stdout)
        self.refreshGeneratedPlots()
        return True

    def viewResults(self):
        '''
        Load and display results on the plotter canvas
        '''
        return self.viewPlot()

    def getKnownPlotLabels(self):
        """
        Return the ordered list of known plot labels for rupture dynamics outputs.
        """
        return [label for label, _ in self.PLOT_FILE_SPECS]

    def refreshGeneratedPlots(self):
        """
        Refresh the cached mapping of available plot labels to file paths.
        """
        plots = OrderedDict()
        if hasattr(self, 'currentCase'):
            for label, filename in self.PLOT_FILE_SPECS:
                plot_path = os.path.join(self.currentCase, filename)
                if os.path.exists(plot_path):
                    plots[label] = plot_path

            # Fallback: include any other PNG outputs in reverse chronological order
            import glob
            png_files = glob.glob(os.path.join(self.currentCase, "*.png"))
            png_files.sort(key=os.path.getmtime, reverse=True)
            for plot_path in png_files:
                if plot_path not in plots.values():
                    plots[os.path.basename(plot_path)] = plot_path

        self.generatedPlots = plots
        return plots

    def getGeneratedPlots(self, refresh=False):
        """
        Return the cached plot mapping, refreshing if requested or empty.
        """
        if refresh or not self.generatedPlots:
            return self.refreshGeneratedPlots()
        return self.generatedPlots

    def viewPlot(self, label=None):
        """
        Display a specific plot by label, or the first available plot if no label provided.
        """
        if not hasattr(self, 'currentCase'):
            print("ERROR: No case results to view!")
            return False

        plots = self.getGeneratedPlots(refresh=True)
        if not plots:
            print("ERROR: No plot files found in case directory")
            return False

        if label:
            plot_file = plots.get(label)
            if not plot_file:
                print(f"ERROR: Plot '{label}' not available")
                return False
        else:
            # Default to the first known plot if available, else the first cached entry
            plot_file = None
            for known_label, _ in self.PLOT_FILE_SPECS:
                if known_label in plots:
                    plot_file = plots[known_label]
                    break
            if plot_file is None:
                plot_file = next(iter(plots.values()))

        return self._displayPlotFile(plot_file)

    def _displayPlotFile(self, plot_file):
        """
        Load an image file and display it on the module plotter.
        """
        print(f"Loading plot: {plot_file}")
        try:
            from PIL import Image

            img = Image.open(plot_file)

            self.plotter.clearFigure()
            ax = self.plotter.getAxis()
            ax.imshow(img)
            ax.axis('off')
            self.plotter.drawFigure()

            return True
        except Exception as e:
            print(f"ERROR loading plot: {e}")
            return False

    def cleanup(self):
        """
        Cleanup when module is closed
        """
        if hasattr(self, 'process') and self.process and self.process.poll() is None:
            self.process.terminate()
        if self._log_handle:
            self._log_handle.close()
            self._log_handle = None
        return False

    def getOutput(self):
        """
        Returns command history
        """
        return self.commandString

    def clearOutput(self):
        """
        Clears command history
        """
        self.commandString = ""
