#!/usr/bin/env python3
"""
Headless smoke test for the SEATREE GUI.

Starts the GTK application, selects a module in the start dialog the
way a user would, drives the module's compute and plot buttons, and
checks that a plot was produced. Meant to be run under xvfb-run, e.g.

    xvfb-run -a python3 test/gui_smoke_test.py hc
    xvfb-run -a python3 test/gui_smoke_test.py conman

The test exits with status 0 on success and 1 on failure, and prints
what it checked. A time limit (default 300 s) guards against hangs.
"""
import os, sys, time

root = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
pydir = os.path.join(root, "python3")
sys.path.insert(0, pydir)
os.environ["SEATREE"] = root
os.environ.setdefault("ARCH", os.uname().machine)
# the module-level code of SEATREE.py needs these for its temp directory
os.environ.setdefault("USER", "seatree")
os.environ.setdefault("HOST", "test")

import gi
gi.require_version("Gtk", "4.0")
from gi.repository import Gtk, GLib

import seatree.gui.SEATREE as S

which = sys.argv[1] if len(sys.argv) > 1 else "hc"
timeout = float(sys.argv[2]) if len(sys.argv) > 2 else 300.0

# class names as in python3/conf/conf.xml
wanted = {"hc": "FlowCalc", "conman": "ConMan"}
if which not in wanted:
    print("usage: gui_smoke_test.py hc|conman [timeout_s]")
    sys.exit(2)

state = {"stage": "start", "ok": False, "msgs": [], "t0": time.time()}

def finish(ok, msg):
    state["msgs"].append(msg)
    print(("PASS: " if ok else "FAIL: ") + msg)
    # do not go through the GTK shutdown path, which calls exit() in
    # several places; the process is done at this point
    sys.stdout.flush()
    os._exit(0 if ok else 1)

def note(msg):
    state["msgs"].append(msg)
    print("smoke: " + msg)
    sys.stdout.flush()

def module_index(app, className):
    for i, m in enumerate(app.modules):
        if m.classname == className:
            return i
    return -1

def drive_hc(app):
    mod = app.window.module
    gui = mod.gui
    note("hc: compute velocities")
    gui.computeVel(None)
    if not gui.plotVelButton.get_sensitive():
        finish(False, "hc: velocity computation did not succeed (plot button stays disabled)")
    note("hc: plot velocities")
    gui.plotVel(None)
    png = mod.gmtPlotterWidget.pngFile
    if not (png and os.path.isfile(png) and os.path.getsize(png) > 1000):
        finish(False, "hc: no velocity plot image produced (" + str(png) + ")")
    note("hc: velocity plot " + png + " " + str(os.path.getsize(png)) + " bytes")
    note("hc: plot geoid")
    gui.plotGeoid(None)
    png = mod.gmtPlotterWidget.pngFile
    if not (png and os.path.isfile(png) and os.path.getsize(png) > 1000):
        finish(False, "hc: no geoid plot image produced (" + str(png) + ")")
    finish(True, "hc: velocities and geoid computed and plotted, last image " + png)

def drive_conman(app):
    mod = app.window.module
    gui = mod.gui
    # keep the run short: smallest grid, fewest steps the GUI allows
    gui.nelzSelect.set_active(0) if hasattr(gui, "nelzSelect") else None
    gui.stepsEntry.setValue(1000)
    gui.stepsSavedEntry.setValue(10)
    note("conman: prepare input with nelz=" + str(gui.getNelz()) + " steps=" + str(gui.getSteps()))
    gui.GUIgawk()
    if not gui.startButton.get_sensitive():
        finish(False, "conman: input preparation failed")
    note("conman: start calculation")
    gui.startCalc()
    state["stage"] = "conman_running"

def poll_conman(app):
    mod = app.window.module
    gui = mod.gui
    if gui.calculating:
        return True
    data = mod.calcThread.getData() if mod.calcThread is not None else None
    nsteps = len(data) if data else 0
    if nsteps == 0 or not gui.hasData:
        finish(False, "conman: run finished without output data")
    fig_ok = mod.plotter is not None
    finish(fig_ok, "conman: run finished with " + str(nsteps) + " output steps, plotted step " + str(mod.getLastPlottedStepNum()))
    return False

def step(app):
    if time.time() - state["t0"] > timeout:
        finish(False, "timeout in stage " + state["stage"])
    try:
        if state["stage"] == "start":
            start = getattr(app, "start", None)
            if start is None:
                return True
            idx = module_index(app, wanted[which])
            if idx < 0:
                finish(False, "module " + wanted[which] + " not in the module list")
            note("selecting module index " + str(idx) + " (" + wanted[which] + ")")
            start.combo.set_active(idx)
            state["stage"] = "loading"
            start.dialog.response(Gtk.ResponseType.OK)
            return True
        if state["stage"] == "loading":
            win = getattr(app, "window", None)
            if win is None or win.module is None or win.module.classname != wanted[which]:
                return True
            if not hasattr(win.module, "gui"):
                return True
            note("module loaded: " + win.module.longname)
            state["stage"] = "drive"
            # give the window one more cycle to realize
            return True
        if state["stage"] == "drive":
            state["stage"] = "driving"
            if which == "hc":
                drive_hc(app)
            else:
                drive_conman(app)
            return True
        if state["stage"] == "conman_running":
            return poll_conman(app)
    except SystemExit:
        raise
    except Exception as e:
        import traceback
        traceback.print_exc()
        finish(False, "exception in stage " + state["stage"] + ": " + str(e))
    return True

def main():
    if not Gtk.init_check():
        finish(False, "GTK could not be initialized (no display?)")
    app = S.SEATREE(path=S.path, storeDir=S.storeDir)
    GLib.timeout_add(500, step, app)
    app.run(None)

if __name__ == "__main__":
    main()
