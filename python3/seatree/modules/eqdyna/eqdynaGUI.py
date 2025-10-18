import gi
gi.require_version('Gtk', '4.0')
from gi.repository import Gtk, GLib, Gio
import os

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

        # Title label
        titleLabel = Gtk.Label()
        titleLabel.set_markup("<b>Earthquake Dynamics Simulation (EQDYNA)</b>")
        titleLabel.set_xalign(0)
        self.append(titleLabel)

        # Add separator
        sep = Gtk.Separator(orientation=Gtk.Orientation.HORIZONTAL)
        self.append(sep)

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

        # Run simulation button
        self.runButton = Gtk.Button(label="Run")
        self.runButton.connect("clicked", self.onRunSimulation)
        self.runButton.set_sensitive(False)
        buttonGrid.attach(self.runButton, 0, 1, 1, 1)

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

        # Help text
        self.addInfoText("EQDYNA - Parallel Finite Element Software for Earthquake Dynamics\n\n")
        self.addInfoText("Instructions:\n")
        self.addInfoText("1. Click 'Browse...' to select a case input folder\n")
        self.addInfoText("   (folder must contain user_defined_params.py)\n")
        self.addInfoText("   Available test cases:\n")
        self.addInfoText("   - test.tpv8, test.tpv10, test.tpv36, test.tpv37\n")
        self.addInfoText("   - test.tpv104, test.tpv1053d, test.drv.a6\n")
        self.addInfoText("2. Enter a case name for your simulation\n")
        self.addInfoText("3. Click 'Create Case' to initialize the simulation\n")
        self.addInfoText("4. Click 'Setup Case' to configure the simulation\n")
        self.addInfoText("5. Click 'Run Simulation' to start the simulation\n\n")
        self.addInfoText("For more information, visit: https://github.com/EQDYNA/eqdyna\n")

    def onBrowseCaseInput(self, button):
        '''
        Handler for Browse button - opens folder chooser dialog
        '''
        chooser = Gtk.FileChooserDialog(
            title="Select Case Input Folder",
            action=Gtk.FileChooserAction.SELECT_FOLDER
        )

        # Default to EQDYNA case_input directory if available
        default_dir = ""
        if hasattr(self.module, "eqdynaPath") and self.module.eqdynaPath:
            default_dir = os.path.join(self.module.eqdynaPath.rstrip(os.sep), "case_input")
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
            self.runButton.set_sensitive(True)
            self.setupButton.set_sensitive(False)
        else:
            self.updateStatus("ERROR: Failed to setup case")

    def onRunSimulation(self, button):
        '''
        Handler for Run Simulation button
        '''
        self.updateStatus("Starting simulation...")
        if self.module.runCase():
            logFile = self.module.getLogFile()
            self.updateStatus(f"Simulation running... (log: {logFile})")
            self.runButton.set_sensitive(False)
            self.viewLogButton.set_sensitive(True)
            # Start a timer to check simulation status
            GLib.timeout_add(1000, self.checkSimulationStatus)
        else:
            self.updateStatus("ERROR: Failed to start simulation")

    def checkSimulationStatus(self):
        '''
        Check if simulation is still running
        '''
        if self.module.checkStatus():
            return True  # Continue checking
        else:
            # Check if simulation generated output files
            import glob
            frt_files = glob.glob(os.path.join(self.module.currentCase, "frt.txt*"))
            if frt_files:
                self.updateStatus(f"Simulation completed ({len(frt_files)} output files)")
                self.postProcessButton.set_sensitive(True)
            else:
                self.updateStatus("Simulation completed (no output files found - check for errors)")
            self.runButton.set_sensitive(True)
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
