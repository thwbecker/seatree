import gi
gi.require_version('Gtk', '4.0')
from gi.repository import Gtk, GLib, Gio
import os
import numpy as np
from scipy.interpolate import griddata

class EQDYNAGUI(Gtk.Box):
    '''
    GUI panel for EQDYNA module
    '''

    def __init__(self, module):
        Gtk.Box.__init__(self, orientation=Gtk.Orientation.VERTICAL, spacing=6)
        self.module = module
        self.set_margin_top(10)
        self.set_margin_bottom(10)
        self.set_margin_start(10)
        self.set_margin_end(10)

        # Case input directory picker
        inputDirBox = Gtk.Box(orientation=Gtk.Orientation.HORIZONTAL, spacing=6)
        inputDirLabel = Gtk.Label(label="Case Input:")
        inputDirLabel.set_xalign(0)
        inputDirLabel.set_size_request(100, -1)

        self.inputDirEntry = Gtk.Entry()
        self.inputDirEntry.set_hexpand(True)
        self.inputDirEntry.set_placeholder_text("Select folder with user_defined_params.py")

        self.browseButton = Gtk.Button(label="Browse...")
        self.browseButton.connect("clicked", self.onBrowseCaseInput)

        self.editParamsButton = Gtk.Button(label="Edit Params")
        self.editParamsButton.connect("clicked", self.onEditParams)
        self.editParamsButton.set_sensitive(False)

        inputDirBox.append(inputDirLabel)
        inputDirBox.append(self.inputDirEntry)
        inputDirBox.append(self.browseButton)
        inputDirBox.append(self.editParamsButton)
        self.append(inputDirBox)

        # Case name input
        caseBox = Gtk.Box(orientation=Gtk.Orientation.HORIZONTAL, spacing=6)
        caseLabel = Gtk.Label(label="Case Name:")
        caseLabel.set_xalign(0)
        caseLabel.set_size_request(100, -1)
        self.caseEntry = Gtk.Entry()
        self.caseEntry.set_text("my_case")
        self.caseEntry.set_hexpand(True)
        caseBox.append(caseLabel)
        caseBox.append(self.caseEntry)
        self.append(caseBox)

        # Buttons in grid layout (2 per row)
        buttonGrid = Gtk.Grid()
        buttonGrid.set_row_spacing(6)
        buttonGrid.set_column_spacing(6)
        buttonGrid.set_column_homogeneous(True)

        # Create case button
        self.createButton = Gtk.Button(label="Create")
        self.createButton.connect("clicked", self.onCreateCase)
        buttonGrid.attach(self.createButton, 0, 0, 1, 1)

        # Setup case button
        self.setupButton = Gtk.Button(label="Setup")
        self.setupButton.connect("clicked", self.onSetupCase)
        self.setupButton.set_sensitive(False)
        buttonGrid.attach(self.setupButton, 1, 0, 1, 1)

        # Start simulation button
        self.startButton = Gtk.Button(label="Start")
        self.startButton.connect("clicked", self.onStartSimulation)
        self.startButton.set_sensitive(False)
        buttonGrid.attach(self.startButton, 0, 1, 1, 1)

        # Post-process button
        self.postProcessButton = Gtk.Button(label="Post-Process")
        self.postProcessButton.connect("clicked", self.onPostProcess)
        self.postProcessButton.set_sensitive(False)
        buttonGrid.attach(self.postProcessButton, 1, 1, 1, 1)

        # View results button (second row)
        self.viewButton = Gtk.Button(label="View Results")
        self.viewButton.connect("clicked", self.onViewResults)
        self.viewButton.set_sensitive(False)
        buttonGrid.attach(self.viewButton, 0, 2, 2, 1)

        # View log button (second row)
        self.viewLogButton = Gtk.Button(label="View Log")
        self.viewLogButton.connect("clicked", self.onViewLog)
        self.viewLogButton.set_sensitive(False)
        buttonGrid.attach(self.viewLogButton, 0, 3, 2, 1)

        self.append(buttonGrid)

        # Panel selection buttons for individual figures
        self.panelButtonsBox = Gtk.FlowBox()
        self.panelButtonsBox.set_selection_mode(Gtk.SelectionMode.NONE)
        self.panelButtonsBox.set_column_spacing(4)
        self.panelButtonsBox.set_row_spacing(4)
        self.panelButtonsBox.set_max_children_per_line(2)
        self.panelButtonsBox.set_min_children_per_line(1)
        self.panelButtonsBox.set_hexpand(True)
        self.panelButtonsBox.set_halign(Gtk.Align.FILL)
        self.panelButtons = {}
        for label in self.module.getKnownPlotLabels():
            button = Gtk.Button(label=label)
            button.connect("clicked", self.onViewPanel, label)
            button.set_sensitive(False)
            button.set_size_request(-1, 30)
            self.panelButtons[label] = button
            self.panelButtonsBox.append(button)
        self.append(self.panelButtonsBox)
        self.setPanelButtonsEnabled(set())

        # Add separator
        sep2 = Gtk.Separator(orientation=Gtk.Orientation.HORIZONTAL)
        self.append(sep2)

        # Status label
        self.statusLabel = Gtk.Label(label="Ready")
        self.statusLabel.set_xalign(0)
        self.statusLabel.set_margin_top(10)
        self.append(self.statusLabel)

        # Progress bar for long-running tasks
        self.progressBar = Gtk.ProgressBar()
        self.progressBar.set_hexpand(True)
        self.progressBar.set_margin_bottom(6)
        self.progressBar.set_visible(False)
        self.append(self.progressBar)
        self._progressPulseId = None

        # Ground motion monitoring state
        self._gm_monitor_id = None
        self._gm_last_file_size = 0  # Track file size to detect changes
        self._gm_num_snapshots = 0
        self._gm_selected_snapshot = -1  # -1 means latest
        self._gm_pending_snapshot_count = None  # Track pending snapshot to plot after delay
        self._gm_delay_timer_id = None  # Timer for delayed plotting
        self._gm_bytes_per_snapshot = None  # Cached snapshot size to avoid repeated file reads
        self._gm_surface_file_mtime = None  # Track surface file modification time
        self._gm_station_count = None  # Cached surface station count
        self._gm_surface_file_locked = False  # Indicates whether the station count has been locked
        self._gm_surface_file_warning_logged = False  # Avoid repeated warnings

        # Snapshot selection controls
        snapshotBox = Gtk.Box(orientation=Gtk.Orientation.HORIZONTAL, spacing=6)
        snapshotLabel = Gtk.Label(label="Snapshot:")
        snapshotLabel.set_xalign(0)
        snapshotLabel.set_size_request(100, -1)

        self.snapshotScale = Gtk.Scale.new_with_range(Gtk.Orientation.HORIZONTAL, 1, 12, 1)
        self.snapshotScale.set_hexpand(True)
        self.snapshotScale.set_value(12)  # Default to latest
        self.snapshotScale.set_draw_value(True)
        self.snapshotScale.set_value_pos(Gtk.PositionType.RIGHT)
        self.snapshotScale.connect("value-changed", self.onSnapshotChanged)
        self.snapshotScale.set_sensitive(False)

        self.latestCheckbox = Gtk.CheckButton(label="Latest")
        self.latestCheckbox.set_active(True)
        self.latestCheckbox.connect("toggled", self.onLatestToggled)

        snapshotBox.append(snapshotLabel)
        snapshotBox.append(self.snapshotScale)
        snapshotBox.append(self.latestCheckbox)
        self.append(snapshotBox)

        # Info text view
        scrolled = Gtk.ScrolledWindow()
        scrolled.set_vexpand(True)
        scrolled.set_policy(Gtk.PolicyType.AUTOMATIC, Gtk.PolicyType.AUTOMATIC)

        self.infoTextView = Gtk.TextView()
        self.infoTextView.set_editable(False)
        self.infoTextView.set_wrap_mode(Gtk.WrapMode.WORD)
        self.infoTextView.set_monospace(True)
        scrolled.set_child(self.infoTextView)

        self.append(scrolled)

    def onBrowseCaseInput(self, button):
        '''
        Handler for Browse button - opens folder chooser dialog
        '''
        chooser = Gtk.FileChooserDialog(
            title="Select Case Input Folder",
            action=Gtk.FileChooserAction.SELECT_FOLDER
        )

        # Default to EQDYNA case_input directory
        default_dir = ""
        if hasattr(self.module, "eqdynaPath") and self.module.eqdynaPath:
            case_input_dir = os.path.join(self.module.eqdynaPath.rstrip(os.sep), "case_input")
            if os.path.isdir(case_input_dir):
                default_dir = case_input_dir
        if default_dir and os.path.isdir(default_dir):
            try:
                chooser.set_current_folder(Gio.File.new_for_path(default_dir))
            except Exception:
                pass

        # Set parent window
        if hasattr(self.module, 'mainWindow'):
            root = self.module.mainWindow
            if isinstance(root, Gtk.Window):
                chooser.set_transient_for(root)
            elif hasattr(root, 'get_root'):
                root = root.get_root()
                if isinstance(root, Gtk.Window):
                    chooser.set_transient_for(root)

        chooser.add_button("_Cancel", Gtk.ResponseType.CANCEL)
        chooser.add_button("_Select", Gtk.ResponseType.OK)

        # GTK4: Use async pattern
        chooser.connect("response", self._on_folder_selected_response)
        chooser.show()

    def _on_folder_selected_response(self, dialog, response):
        '''
        Handle folder selection dialog response
        '''
        folder_path = None
        if response == Gtk.ResponseType.OK:
            gfile = dialog.get_file()
            if gfile:
                folder_path = gfile.get_path()
        dialog.destroy()

        if folder_path is not None:
            # Check if user_defined_params.py exists in selected folder
            params_file = os.path.join(folder_path, "user_defined_params.py")
            if os.path.exists(params_file):
                self.inputDirEntry.set_text(folder_path)
                self.updateStatus(f"Selected case input: {os.path.basename(folder_path)}")
                self.editParamsButton.set_sensitive(True)
                # Auto-populate case name from folder name
                case_name = os.path.basename(folder_path)
                self.caseEntry.set_text(case_name)
            else:
                self.updateStatus("ERROR: Selected folder does not contain user_defined_params.py")
                self.editParamsButton.set_sensitive(False)
        else:
            self.editParamsButton.set_sensitive(False)

    def onEditParams(self, button):
        '''
        Open dialog to edit user_defined_params.py in selected directory
        '''
        folder_path = self.inputDirEntry.get_text().strip()
        if not folder_path:
            self.updateStatus("ERROR: Please select a case input directory")
            return

        params_path = os.path.join(folder_path, "user_defined_params.py")
        if not os.path.exists(params_path):
            self.updateStatus("ERROR: user_defined_params.py not found in selected folder")
            self.editParamsButton.set_sensitive(False)
            return

        try:
            with open(params_path, "r", encoding="utf-8") as params_file:
                params_content = params_file.read()
        except OSError as exc:
            self.updateStatus(f"ERROR: Failed to read params file ({exc})")
            return

        dialog = Gtk.Dialog(title="Edit user_defined_params.py")
        dialog.set_default_size(700, 500)
        dialog.set_modal(True)

        # Set parent window when available
        parent_window = getattr(self.module, "mainWindow", None)
        if isinstance(parent_window, Gtk.Window):
            dialog.set_transient_for(parent_window)
        elif hasattr(parent_window, "get_root"):
            root = parent_window.get_root()
            if isinstance(root, Gtk.Window):
                dialog.set_transient_for(root)

        dialog.add_button("_Cancel", Gtk.ResponseType.CANCEL)
        dialog.add_button("_Save", Gtk.ResponseType.OK)

        content_area = dialog.get_content_area()
        scrolled = Gtk.ScrolledWindow()
        scrolled.set_hexpand(True)
        scrolled.set_vexpand(True)

        text_view = Gtk.TextView()
        text_view.set_wrap_mode(Gtk.WrapMode.NONE)
        text_view.set_monospace(True)
        buffer = text_view.get_buffer()
        buffer.set_text(params_content)
        scrolled.set_child(text_view)
        content_area.append(scrolled)

        def on_response(dlg, response_id):
            if response_id == Gtk.ResponseType.OK:
                start_iter = buffer.get_start_iter()
                end_iter = buffer.get_end_iter()
                updated_content = buffer.get_text(start_iter, end_iter, True)
                try:
                    with open(params_path, "w", encoding="utf-8") as params_file:
                        params_file.write(updated_content)
                    self.updateStatus(f"Saved {os.path.basename(params_path)}")
                except OSError as exc:
                    self.updateStatus(f"ERROR: Failed to save params file ({exc})")
            dlg.destroy()

        dialog.connect("response", on_response)
        dialog.present()

    def onCreateCase(self, button):
        '''
        Handler for Create Case button
        '''
        caseName = self.caseEntry.get_text().strip()
        if not caseName:
            self.updateStatus("ERROR: Please enter a case name")
            return

        # Get the case input directory from the entry field
        caseInputDir = self.inputDirEntry.get_text().strip()
        if not caseInputDir:
            self.updateStatus("ERROR: Please select a case input directory")
            return

        if not os.path.exists(caseInputDir):
            self.updateStatus("ERROR: Selected case input directory does not exist")
            return

        self.updateStatus(f"Creating case: {caseName}...")
        if self.module.createCase(caseName, caseInputDir):
            self.updateStatus(f"Case '{caseName}' created successfully")
            self.setupButton.set_sensitive(True)
            self.createButton.set_sensitive(False)
            self.setPanelButtonsEnabled(set())
        else:
            self.updateStatus("ERROR: Failed to create case")

    def onSetupCase(self, button):
        '''
        Handler for Setup Case button
        '''
        self.updateStatus("Setting up case...")
        if self.module.setupCase():
            self.updateStatus("Case setup completed")
            self.startButton.set_sensitive(True)
            self.setupButton.set_sensitive(False)
        else:
            self.updateStatus("ERROR: Failed to setup case")

    def onStartSimulation(self, button):
        '''
        Handler for Start Simulation button
        '''
        self.updateStatus("Starting simulation...")
        if self.module.runCase():
            logFile = self.module.getLogFile()
            self.updateStatus(f"Simulation running... (log: {logFile})")
            self.startProgressPulse()
            self.startButton.set_sensitive(False)
            self.viewLogButton.set_sensitive(True)
            # Start monitoring for ground motion files
            self.startGroundMotionMonitoring()
            # Start a timer to check simulation status
            GLib.timeout_add(1000, self.checkSimulationStatus)
        else:
            self.updateStatus("ERROR: Failed to start simulation")
            self.stopProgressPulse()

    def checkSimulationStatus(self):
        '''
        Check if simulation is still running
        '''
        if self.module.checkStatus():
            return True  # Continue checking
        else:
            # Simulation completed
            self.stopGroundMotionMonitoring()
            # Check if simulation generated output files
            import glob
            frt_files = glob.glob(os.path.join(self.module.currentCase, "frt.txt*"))
            if frt_files:
                self.updateStatus(f"Simulation completed ({len(frt_files)} output files)")
                self.postProcessButton.set_sensitive(True)
            else:
                self.updateStatus("Simulation completed (no output files found - check for errors)")
            self.stopProgressPulse()
            self.startButton.set_sensitive(True)
            return False  # Stop checking

    def onPostProcess(self, button):
        '''
        Handler for Post-Process button
        '''
        self.updateStatus("Post-processing results...")
        if self.module.postProcess():
            self.updateStatus("Post-processing completed")
            self.viewButton.set_sensitive(True)
            self.updatePanelButtons()
        else:
            self.updateStatus("ERROR: Failed to post-process results")
            self.setPanelButtonsEnabled(set())

    def onViewResults(self, button):
        '''
        Handler for View Results button
        '''
        self.updateStatus("Loading results...")
        if self.module.viewResults():
            self.updateStatus("Results displayed")
            self.updatePanelButtons()
        else:
            self.updateStatus("ERROR: Failed to load results")
            self.updatePanelButtons()

    def onViewLog(self, button):
        '''
        Handler for View Log button
        '''
        log_content = self.module.readLog()
        # Clear info text and show log
        buffer = self.infoTextView.get_buffer()
        buffer.set_text("")
        buffer.insert(buffer.get_end_iter(), "=== SIMULATION LOG ===\n\n")
        buffer.insert(buffer.get_end_iter(), log_content)
        self.updateStatus("Log displayed")

    def updateStatus(self, message):
        '''
        Update the status label
        '''
        self.statusLabel.set_text(message)
        print(message)

    def addInfoText(self, text):
        '''
        Add text to the info text view
        '''
        buffer = self.infoTextView.get_buffer()
        end_iter = buffer.get_end_iter()
        buffer.insert(end_iter, text)

    def setPanelButtonsEnabled(self, enabled_labels):
        '''
        Enable buttons matching labels in enabled_labels and disable the rest.
        '''
        if enabled_labels is None:
            enabled_labels = set()
        for label, button in self.panelButtons.items():
            button.set_sensitive(label in enabled_labels)

    def updatePanelButtons(self):
        '''
        Refresh panel buttons based on available plot outputs.
        '''
        try:
            plots = self.module.getGeneratedPlots(refresh=True)
        except Exception:
            plots = {}
        available = {label for label in self.panelButtons if label in plots}
        self.setPanelButtonsEnabled(available)

    def onViewPanel(self, button, label):
        '''
        Display a specific rupture dynamics panel.
        '''
        self.updateStatus(f"Loading {label} plot...")
        if self.module.viewPlot(label):
            self.updateStatus(f"{label} displayed")
        else:
            self.updateStatus(f"ERROR: Failed to load {label}")
        self.updatePanelButtons()

    def onSnapshotChanged(self, scale):
        '''
        Handler for snapshot slider change
        '''
        if not self.latestCheckbox.get_active():
            snapshot_num = int(scale.get_value())
            self._gm_selected_snapshot = snapshot_num - 1  # Convert to 0-based index
            self.plotGroundMotion()

    def onLatestToggled(self, checkbox):
        '''
        Handler for Latest checkbox toggle
        '''
        is_latest = checkbox.get_active()
        self.snapshotScale.set_sensitive(not is_latest)
        if is_latest:
            self._gm_selected_snapshot = -1
            self.plotGroundMotion()
        else:
            snapshot_num = int(self.snapshotScale.get_value())
            self._gm_selected_snapshot = snapshot_num - 1
            self.plotGroundMotion()

    def startProgressPulse(self):
        '''
        Begin pulsing the progress bar for long-running tasks.
        '''
        if self._progressPulseId is None:
            self.progressBar.set_fraction(0.0)
            self.progressBar.set_visible(True)
            self._progressPulseId = GLib.timeout_add(100, self._pulseProgress)

    def stopProgressPulse(self):
        '''
        Stop the progress bar pulse and hide it.
        '''
        if self._progressPulseId is not None:
            GLib.source_remove(self._progressPulseId)
            self._progressPulseId = None
        self.progressBar.set_fraction(0.0)
        self.progressBar.set_visible(False)

    def _pulseProgress(self):
        '''
        Internal helper to pulse the progress bar.
        '''
        if self._progressPulseId is None:
            return False
        self.progressBar.pulse()
        return True

    def startGroundMotionMonitoring(self):
        '''
        Start monitoring ground motion files and update plot in real-time
        Wait for surface_coor.txt to be fully written before starting
        '''
        if self._gm_monitor_id is None:
            self._gm_last_file_size = 0  # Track number of complete snapshots
            self._gm_bytes_per_snapshot = None  # Reset cached snapshot size
            self._gm_station_count = None  # Track how many stations were detected
            self._gm_surface_file_locked = False  # Lock once we trust the surface file contents
            self._gm_surface_file_warning_logged = False  # Avoid spamming warnings
            self._gm_surface_file_mtime = None  # Last known mtime of the surface file
            # Start by waiting for surface coordinates file to be ready
            self._gm_monitor_id = GLib.timeout_add(500, self.waitForSurfaceFile)
            print("Waiting for surface_coor.txt to be written...")

    def waitForSurfaceFile(self):
        '''
        Wait for surface_coor.txt file to be fully written before starting monitoring
        '''
        if not hasattr(self.module, 'currentCase'):
            return False

        import glob
        case_dir = self.module.currentCase

        # Look for any surface_coor.txt file (from any MPI rank)
        surface_files = glob.glob(os.path.join(case_dir, "surface_coor.txt*"))

        if not surface_files:
            return True  # Keep waiting

        # Check if the first surface file exists and has stable size
        surface_file = surface_files[0]
        if not os.path.exists(surface_file):
            return True  # Keep waiting

        # Get file size and modification time
        try:
            current_size = os.path.getsize(surface_file)
            current_mtime = os.path.getmtime(surface_file)

            # Check if we've seen this file before
            if not hasattr(self, '_surface_wait_size'):
                self._surface_wait_size = current_size
                self._surface_wait_mtime = current_mtime
                self._surface_wait_checks = 0
                print(f"Found surface file: {surface_file}, waiting for it to stabilize...")
                return True  # Wait one more cycle to check stability

            # File is stable if size and mtime haven't changed for at least 10 consecutive checks
            if current_size == self._surface_wait_size and current_mtime == self._surface_wait_mtime:
                self._surface_wait_checks += 1
                print(f"Surface file stable check {self._surface_wait_checks}/10 (size={current_size} bytes)")

                if self._surface_wait_checks >= 10:
                    print(f"✓ Surface file ready: {surface_file} ({current_size} bytes)")
                    # Clean up temporary tracking variables
                    del self._surface_wait_size
                    del self._surface_wait_mtime
                    del self._surface_wait_checks
                    # Switch to actual ground motion monitoring
                    self._gm_monitor_id = GLib.timeout_add(500, self.checkGroundMotionFiles)
                    print("Started ground motion monitoring")
                    return False  # Stop waiting, monitoring has taken over
                else:
                    return True  # Keep waiting for more stable checks
            else:
                # File is still changing, reset counter
                print(f"Surface file changed: size {self._surface_wait_size} → {current_size}")
                self._surface_wait_size = current_size
                self._surface_wait_mtime = current_mtime
                self._surface_wait_checks = 0
                return True  # Keep waiting

        except Exception as e:
            print(f"Error checking surface file: {e}")
            return True  # Keep waiting

    def stopGroundMotionMonitoring(self):
        '''
        Stop monitoring ground motion files
        '''
        if self._gm_monitor_id is not None:
            GLib.source_remove(self._gm_monitor_id)
            self._gm_monitor_id = None
        if self._gm_delay_timer_id is not None:
            GLib.source_remove(self._gm_delay_timer_id)
            self._gm_delay_timer_id = None
        self._gm_pending_snapshot_count = None
        # Clean up any temporary waiting state
        if hasattr(self, '_surface_wait_size'):
            del self._surface_wait_size
        if hasattr(self, '_surface_wait_mtime'):
            del self._surface_wait_mtime
        print("Stopped ground motion monitoring")

    def checkGroundMotionFiles(self):
        '''
        Check if ground motion files have new data and update plot
        Only plot when a complete snapshot is available (file size is exact multiple of snapshot size)
        '''
        if not hasattr(self.module, 'currentCase'):
            return False

        import glob
        case_dir = self.module.currentCase
        gm_files = sorted(glob.glob(os.path.join(case_dir, "gm[0-9]*")))

        if not gm_files:
            return True  # Continue checking

        # Check the first gm file to determine if we have complete snapshots
        first_gm = gm_files[0]
        if not os.path.exists(first_gm):
            return True

        # Get rank and corresponding surface file
        rank = int(os.path.basename(first_gm)[2:]) if len(os.path.basename(first_gm)) > 2 else 0
        surface_file = os.path.join(case_dir, f"surface_coor.txt{rank}")

        if not os.path.exists(surface_file):
            return True

        # Check if surface file has been modified - recalculate snapshot size if needed
        surface_mtime = os.path.getmtime(surface_file)
        if self._gm_bytes_per_snapshot is None or self._gm_bytes_per_snapshot == 0:
            # Count lines in surface file to get station count (it's a text file!)
            with open(surface_file, 'r') as f:
                n_stations = sum(1 for line in f if line.strip())

            if n_stations == 0:
                # File exists but has no stations yet - keep waiting.
                return True

            # Calculate expected bytes per snapshot: n_stations * 3 components * 8 bytes (float64)
            # The gm* files are BINARY files containing float64 data
            new_bytes_per_snapshot = n_stations * 3 * 8

            if self._gm_bytes_per_snapshot is None:
                print(f"Cached snapshot size: {new_bytes_per_snapshot} bytes ({n_stations} stations)")
            else:
                print(f"Surface file updated! Recalculated snapshot size: {new_bytes_per_snapshot} bytes ({n_stations} stations)")

            self._gm_bytes_per_snapshot = new_bytes_per_snapshot
            self._gm_surface_file_mtime = surface_mtime
            self._gm_station_count = n_stations
            self._gm_surface_file_locked = True
            self._gm_surface_file_warning_logged = False
        elif self._gm_surface_file_locked and surface_mtime != self._gm_surface_file_mtime:
            # The surface file changed after we locked the station count. This happens when
            # EQdyna appends diagnostic fault nodes at shutdown. Keep the original station count.
            if not self._gm_surface_file_warning_logged:
                print("surface_coor.txt was modified after monitoring started; keeping the original "
                      f"{self._gm_station_count} surface stations for snapshot sizing.")
                self._gm_surface_file_warning_logged = True
            self._gm_surface_file_mtime = surface_mtime

        # Get current file size
        current_size = os.path.getsize(first_gm)

        # Calculate how many complete snapshots we can read from the file
        # Integer division gives us the number of complete snapshots (ignore any partial remainder)
        complete_snapshots = current_size // self._gm_bytes_per_snapshot

        # Debug: print status when new data appears
        if complete_snapshots > self._gm_last_file_size:
            print(f"GM file check: size={current_size} bytes, complete_snapshots={complete_snapshots}, last_plotted={self._gm_last_file_size}")

        # Plot if we have a NEW complete snapshot available
        # We can safely read complete_snapshots worth of data, ignoring any partial remainder
        if (complete_snapshots > 0 and
            complete_snapshots > self._gm_last_file_size and
            complete_snapshots != self._gm_pending_snapshot_count):

            print(f"✓ Detected NEW complete snapshot {complete_snapshots} (file has {current_size} bytes, need {self._gm_bytes_per_snapshot} per snapshot)")
            print(f"  Waiting 200ms to ensure all MPI ranks finish writing...")

            # Cancel any existing delay timer
            if self._gm_delay_timer_id is not None:
                GLib.source_remove(self._gm_delay_timer_id)

            # Store the snapshot count we're about to plot
            self._gm_pending_snapshot_count = complete_snapshots

            # Schedule the plot after a small delay (200ms) to ensure all MPI ranks have written
            self._gm_delay_timer_id = GLib.timeout_add(200, self._delayedPlotGroundMotion)

        # Continue monitoring if timer is still active
        return self._gm_monitor_id is not None

    def _delayedPlotGroundMotion(self):
        '''
        Plot ground motion after delay to ensure all MPI ranks have finished writing
        '''
        if self._gm_pending_snapshot_count is not None:
            print(f"Plotting snapshot {self._gm_pending_snapshot_count} after delay")
            self._gm_last_file_size = self._gm_pending_snapshot_count
            self._gm_pending_snapshot_count = None
            self.plotGroundMotion()

        self._gm_delay_timer_id = None
        return False  # Don't repeat timer

    def plotGroundMotion(self):
        '''
        Read and plot ground motion data from gm* files
        '''
        try:
            import glob
            case_dir = self.module.currentCase

            # Find all gm files
            gm_files = sorted(glob.glob(os.path.join(case_dir, "gm[0-9]*")))
            if not gm_files:
                return

            # Read dt from user_defined_params.py
            dt = self._read_dt_from_params()

            # Use the number of complete snapshots already verified by checkGroundMotionFiles
            # This avoids recalculating and ensures consistency
            num_snapshots = self._gm_last_file_size

            if num_snapshots < 1:
                return  # Not enough data yet

            # Update snapshot slider range if needed
            if num_snapshots != self._gm_num_snapshots:
                self._gm_num_snapshots = num_snapshots
                self.snapshotScale.set_range(1, num_snapshots)
                if self.latestCheckbox.get_active():
                    self.snapshotScale.set_value(num_snapshots)
                self.snapshotScale.set_sensitive(True)

            # Determine which snapshot to plot
            if self._gm_selected_snapshot < 0 or self.latestCheckbox.get_active():
                # Plot latest snapshot (monitoring mode)
                # We can safely plot the latest because checkGroundMotionFiles ensures
                # we only plot complete snapshots (with 200ms delay for MPI synchronization)
                snapshot_index = num_snapshots - 1
                print(f"Plotting latest snapshot: {snapshot_index + 1} of {num_snapshots}")
            else:
                # Plot user-selected snapshot
                snapshot_index = min(self._gm_selected_snapshot, num_snapshots - 1)
                print(f"Plotting selected snapshot: {snapshot_index + 1} of {num_snapshots}")

            # Read the selected snapshot from all MPI ranks
            all_velocities = []
            all_coords = []

            for gm_file in gm_files:
                # Get MPI rank from filename
                rank = int(os.path.basename(gm_file)[2:]) if len(os.path.basename(gm_file)) > 2 else 0
                surface_file = os.path.join(case_dir, f"surface_coor.txt{rank}")

                if not os.path.exists(surface_file):
                    continue

                # Read surface coordinates
                coords = np.loadtxt(surface_file)
                if len(coords.shape) == 1:
                    coords = coords.reshape(1, -1)

                n_stations_rank = coords.shape[0]
                print(f"DEBUG: Loaded {n_stations_rank} stations from {surface_file}")

                # Read velocities for the specified snapshot
                # Binary format: all stations for each timestep
                # Layout: [t0: s0_vx,s0_vy,s0_vz, s1_vx,s1_vy,s1_vz, ..., t1: s0_vx,s0_vy,s0_vz, ...]
                values_per_snapshot_rank = n_stations_rank * 3
                bytes_per_snapshot_rank = values_per_snapshot_rank * 8

                # Seek to the start of this snapshot
                byte_offset = snapshot_index * bytes_per_snapshot_rank

                with open(gm_file, 'rb') as f:
                    f.seek(byte_offset)
                    snapshot_data = np.fromfile(f, dtype=np.float64, count=values_per_snapshot_rank)

                # Debug: check if we got all the data
                if len(snapshot_data) != values_per_snapshot_rank:
                    print(f"WARNING: Expected {values_per_snapshot_rank} values, got {len(snapshot_data)}")

                # Reshape to [n_stations_rank, 3]
                velocities = snapshot_data.reshape(n_stations_rank, 3)
                all_coords.append(coords)
                all_velocities.append(velocities)

            if not all_coords:
                return

            # Combine data from all MPI ranks
            coords = np.vstack(all_coords)
            velocities = np.vstack(all_velocities)

            print(f"DEBUG: Combined {len(coords)} total stations from {len(all_coords)} MPI rank(s)")

            # Extract velocity components
            vx = velocities[:, 0]
            vy = velocities[:, 1]
            vz = velocities[:, 2]

            # Calculate actual simulation time
            # Time = timestep * dt * 10 (GM written every 10 timesteps)
            if dt > 0:
                simulation_time = (snapshot_index + 1) * dt * 10
                title_time = f't = {simulation_time:.3f} s (snapshot {snapshot_index + 1})'
            else:
                title_time = f'Snapshot {snapshot_index + 1}'

            # Plot on the module's plotter with 3 subplots
            plotter = self.module.getPlotter()
            plotter.clearFigure()

            # Create 1x3 subplot layout
            fig = plotter.figure
            fig.clear()

            # Convert coordinates from meters to kilometers for display
            x_coords_km = coords[:, 0] / 1000.0
            y_coords_km = coords[:, 1] / 1000.0

            # Read dx from params to build regular grid
            dx = self._read_dx_from_params()
            if dx <= 0:
                dx = 100.0  # Default 100m
            dx_km = dx / 1000.0  # Convert to km
            dy_km = dx_km  # Assume dy = dx

            # Get map extent
            map_x_min, map_x_max = x_coords_km.min(), x_coords_km.max()
            map_y_min, map_y_max = y_coords_km.min(), y_coords_km.max()

            # Calculate number of grid points based on dx
            nx = round((map_x_max - map_x_min) / dx_km)
            ny = round((map_y_max - map_y_min) / dy_km)

            # Create regular grid using linspace
            xi = np.linspace(map_x_min, map_x_max, num=nx)
            yi = np.linspace(map_y_min, map_y_max, num=ny)
            xi_grid, yi_grid = np.meshgrid(xi, yi)

            # Interpolate velocities onto regular grid (griddata handles duplicates automatically)
            vx_grid = griddata((x_coords_km, y_coords_km), vx, (xi_grid, yi_grid), method='linear')
            vy_grid = griddata((x_coords_km, y_coords_km), vy, (xi_grid, yi_grid), method='linear')
            vz_grid = griddata((x_coords_km, y_coords_km), vz, (xi_grid, yi_grid), method='linear')

            # Set fixed color scale for all three components
            vmin = -0.2
            vmax = 0.2

            # Create subplots with space for shared colorbar at bottom
            from matplotlib.gridspec import GridSpec
            gs = GridSpec(2, 3, figure=fig, height_ratios=[30, 1], hspace=0.15, wspace=0.25,
                         top=0.95, bottom=0.08, left=0.05, right=0.98)

            # Add main title with time info
            fig.suptitle(f'{title_time}', fontsize=11, y=0.98)

            # Plot vx
            ax1 = fig.add_subplot(gs[0, 0])
            mesh1 = ax1.pcolormesh(xi_grid, yi_grid, vx_grid, cmap='seismic',
                                   vmin=vmin, vmax=vmax, shading='auto')
            ax1.set_xlabel('X (km)', fontsize=9)
            ax1.set_ylabel('Y (km)', fontsize=9)
            ax1.set_title('Vx', fontsize=10)
            ax1.set_aspect('equal')
            ax1.tick_params(labelsize=8)

            # Plot vy
            ax2 = fig.add_subplot(gs[0, 1])
            mesh2 = ax2.pcolormesh(xi_grid, yi_grid, vy_grid, cmap='seismic',
                                   vmin=vmin, vmax=vmax, shading='auto')
            ax2.set_xlabel('X (km)', fontsize=9)
            ax2.set_ylabel('Y (km)', fontsize=9)
            ax2.set_title('Vy', fontsize=10)
            ax2.set_aspect('equal')
            ax2.tick_params(labelsize=8)

            # Plot vz
            ax3 = fig.add_subplot(gs[0, 2])
            mesh3 = ax3.pcolormesh(xi_grid, yi_grid, vz_grid, cmap='seismic',
                                   vmin=vmin, vmax=vmax, shading='auto')
            ax3.set_xlabel('X (km)', fontsize=9)
            ax3.set_ylabel('Y (km)', fontsize=9)
            ax3.set_title('Vz', fontsize=10)
            ax3.set_aspect('equal')
            ax3.tick_params(labelsize=8)

            # Add single horizontal colorbar at bottom spanning all three plots
            cbar_ax = fig.add_subplot(gs[1, :])
            cbar = fig.colorbar(mesh1, cax=cbar_ax, orientation='horizontal')
            cbar.set_label('Velocity (m/s)', fontsize=9)
            cbar.ax.tick_params(labelsize=8)

            plotter.drawFigure()

            # Calculate max velocity for status
            max_vel = max(abs(vx).max(), abs(vy).max(), abs(vz).max())

            if dt > 0:
                print(f"Plotted ground motion: {len(coords)} nodes, t={simulation_time:.3f}s, max vel = {max_vel:.3e} m/s")
                self.updateStatus(f"Ground motion: t={simulation_time:.3f}s, {len(coords)} stations, max vel = {max_vel:.3e} m/s")
            else:
                print(f"Plotted ground motion: {len(coords)} nodes, snapshot {snapshot_index + 1}, max vel = {max_vel:.3e} m/s")
                self.updateStatus(f"Ground motion: snapshot {snapshot_index + 1}, {len(coords)} stations, max vel = {max_vel:.3e} m/s")

        except Exception as e:
            print(f"Error plotting ground motion: {e}")
            import traceback
            traceback.print_exc()

    def _read_dt_from_params(self):
        '''
        Read timestep (dt) from user_defined_params.py
        Returns dt in seconds, or 0 if not found
        '''
        try:
            if not hasattr(self.module, 'currentCase'):
                return 0.0

            params_file = os.path.join(self.module.currentCase, 'user_defined_params.py')
            if not os.path.exists(params_file):
                return 0.0

            # Read the file and look for dt calculation
            with open(params_file, 'r') as f:
                content = f.read()

            # Try to extract dx, vp, and calculate dt
            # dt is typically: dt = 0.5*dx/vp
            dx = None
            vp = None

            for line in content.split('\n'):
                line = line.strip()
                if line.startswith('par.dx') and '=' in line:
                    try:
                        dx = float(line.split('=')[1].strip().rstrip('.'))
                    except:
                        pass
                elif line.startswith('par.vp') and '=' in line:
                    try:
                        # Handle cases like: par.vp, par.vs, par.rou = 5716, 3300, 2700
                        rhs = line.split('=')[1].strip()
                        values = [v.strip() for v in rhs.split(',')]
                        if values:
                            vp = float(values[0])
                    except:
                        pass
                elif line.startswith('par.dt') and '=' in line and 'par.dx' not in line and 'par.vp' not in line:
                    # Direct dt assignment (not a formula)
                    try:
                        dt_str = line.split('=')[1].strip()
                        return float(dt_str)
                    except:
                        pass

            # If we found dx and vp, calculate dt using CFL condition
            if dx is not None and vp is not None and vp > 0:
                dt = 0.5 * dx / vp
                print(f"Calculated dt = {dt:.6f} s (dx={dx}, vp={vp})")
                return dt

            return 0.0

        except Exception as e:
            print(f"Error reading dt from params: {e}")
            return 0.0

    def _read_dx_from_params(self):
        '''
        Read grid spacing (dx) from user_defined_params.py
        Returns dx in meters, or 0 if not found
        '''
        try:
            if not hasattr(self.module, 'currentCase'):
                return 0.0

            params_file = os.path.join(self.module.currentCase, 'user_defined_params.py')
            if not os.path.exists(params_file):
                return 0.0

            # Read the file and look for dx
            with open(params_file, 'r') as f:
                content = f.read()

            for line in content.split('\n'):
                line = line.strip()
                if line.startswith('par.dx') and '=' in line:
                    try:
                        dx = float(line.split('=')[1].strip().rstrip('.'))
                        print(f"Found dx = {dx} m from user_defined_params.py")
                        return dx
                    except:
                        pass

            return 0.0

        except Exception as e:
            print(f"Error reading dx from params: {e}")
            return 0.0
