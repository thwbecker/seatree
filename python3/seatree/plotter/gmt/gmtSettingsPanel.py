import gi
gi.require_version('Gtk', '4.0')
from gi.repository import Gtk, Gdk

import seatree.gui.util.guiUtils as guiUtils
from seatree.gmt.gmtWrapper import *
from seatree.gui.util.rangeDialog import RangeDialog

class GMTSettingsPanel:
    
    def __init__(self, gmtPlotterWidget, gmtWrapper, width=200):
        self.width = width
        self.gmtPlotterWidget = gmtPlotterWidget

        self.eb = Gtk.Box()
        self.eb.set_size_request(self.width, 200)
        #self.eb.override_background_color(Gtk.StateFlags.NORMAL, Gdk.RGBA(0.5, 0.5, 0.5, 1))

        self.vBox = Gtk.Box(orientation=Gtk.Orientation.VERTICAL, spacing=1)
        self.row1 = Gtk.Box(orientation=Gtk.Orientation.HORIZONTAL, spacing=1)
        
        # Colormap Type
        self.colormapCombo = Gtk.ComboBoxText()
        self.colormapCombo.append_text("haxby")
        self.colormapCombo.append_text("gray")
        self.colormapCombo.append_text("wysiwyg")
        self.colormapCombo.append_text("polar")
        self.colormapCombo.append_text("seis")
        self.colormapCombo.append_text("spectral")
        self.colormapCombo.append_text("other...")
        self.colormapCombo.set_active(0)

        self.colormapCombo.connect("changed", self.checkColormap)

        self.colormapEntry = Gtk.Entry()
        self.colormapEntry.set_width_chars(5)
        self.colormapEntry.set_sensitive(False)

        self.colormapLabel = Gtk.Label(label="Colormap Type:")
        self.colormapBox = Gtk.Box(orientation=Gtk.Orientation.HORIZONTAL, spacing=5)
        self.colormapBox.append(self.colormapLabel)
        
        self.invertCheck = Gtk.CheckButton(label="Invert")
        self.invertCheckBox = Gtk.Box(orientation=Gtk.Orientation.HORIZONTAL, spacing=0)
        self.invertCheckBox.append(self.invertCheck)

        self.adjustCheck = Gtk.CheckButton(label="Auto scale")
        self.adjustCheck.set_active(True)
        self.adjustCheckBox = Gtk.Box(orientation=Gtk.Orientation.HORIZONTAL, spacing=0)
        self.adjustCheckBox.append(self.adjustCheck)

        self.labelCheck = Gtk.CheckButton(label="Labels")
        self.labelCheck.set_active(True)
        self.labelCheckBox = Gtk.Box(orientation=Gtk.Orientation.HORIZONTAL, spacing=0)
        self.labelCheckBox.append(self.labelCheck)

        self.plateBoundaryCheck = Gtk.CheckButton(label="Pbound")
        self.plateBoundaryCheckBox = Gtk.Box(orientation=Gtk.Orientation.HORIZONTAL, spacing=0)
        self.plateBoundaryCheckBox.append(self.plateBoundaryCheck)
        
        self.coastCheck = Gtk.CheckButton(label="Coastlines")
        self.coastCheckBox = Gtk.Box(orientation=Gtk.Orientation.HORIZONTAL, spacing=0)
        self.coastCheckBox.append(self.coastCheck)
        
        self.coastMaskCombo = Gtk.ComboBoxText()
        self.coastMaskCombo.append_text("No Mask")
        self.coastMaskCombo.append_text("Mask Sea")
        self.coastMaskCombo.append_text("Mask Land")
        self.coastMaskCombo.set_active(0)

        self.plotGridLines = Gtk.CheckButton(label="Grid lines")
        self.plotGridLinesBox = Gtk.Box(orientation=Gtk.Orientation.HORIZONTAL, spacing=0)
        self.plotGridLinesBox.append(self.plotGridLines)

        self.colorSymLabel = Gtk.Label(label="Colors:")

        self.coastColorButton = Gtk.Button(label="Coast")
        self.coastColorButton.connect("clicked", self.setCoastColor)
        
        self.coastMaskColorButton = Gtk.Button(label="Mask")
        self.coastMaskColorButton.connect("clicked", self.setCoastMaskColor)

        self.pbcolButton = Gtk.Button(label="Pbound")
        self.pbcolButton.connect("clicked", self.setPbcol)

        self.vecColButton = Gtk.Button(label="Vectors")
        self.vecColButton.connect("clicked", self.setVectColor)
        
        # grid resolution
        self.gridres = Gtk.Entry()
        self.gridres.set_width_chars(4)
        self.gridres.set_max_length(4)
        self.gridres.set_text("1.0")
        self.gridresLabel = Gtk.Label(label="resolution")

        # Range Settings
        self.rangeButton = Gtk.Button(label="Update Lat/Lon Range")
        self.rangeButton.connect("clicked", self.openRangeWindow)

        self.row1.append(self.colormapBox)
        self.row1.append(self.colormapCombo)
        self.row1.append(self.colormapEntry)
        self.row1.append(self.invertCheckBox)
        self.row1.append(self.adjustCheckBox)
        self.row1.append(self.labelCheckBox)
        self.row1.append(self.rangeButton)

        self.row2 = Gtk.Box(orientation=Gtk.Orientation.HORIZONTAL, spacing=5)

        self.applyChangesButton = Gtk.Button(label="Apply Changes")
        self.applyChangesButton.connect("clicked", self.applyChanges)

        self.row2.append(self.plotGridLinesBox )
        self.row2.append(self.coastCheckBox )
        self.row2.append(self.coastMaskCombo )
        self.row2.append(self.plateBoundaryCheckBox )
        self.row2.append(self.colorSymLabel )
        self.row2.append(self.pbcolButton )
        self.row2.append(self.coastColorButton )
        self.row2.append(self.vecColButton )
        self.row2.append(self.gridres )
        self.row2.append(self.gridresLabel )

        self.row3 = Gtk.Box(orientation=Gtk.Orientation.HORIZONTAL, spacing=5)

        # Projection
        self.projectionLabel = Gtk.Label(label="Projection:")
        self.projectionEntry = guiUtils.RangeSelectionBox(initial=170, min1=0, max1=359, digits=0, buttons=True)
        self.projectionEntryEnd = Gtk.Entry()
        self.projectionEntryEnd.set_width_chars(2)
        self.projectionEntryEnd.set_max_length(2)
        self.projectionEntryEnd.set_text("7")
        self.projectionLabel2 = Gtk.Label(label="PS width (inches)")

        self.projectionBox = Gtk.Box(orientation=Gtk.Orientation.HORIZONTAL, spacing=50)

        # projection type
        self.projectCombo = Gtk.ComboBoxText()
        self.projections = ['W Mollweide', 'H Hammer', 'N Robinson', 'Q Equi-Cyl', 'J Miller', 'M Mercator', 'Kf Eckert-IV', 'X Linear']

        for pt in self.projections:
            self.projectCombo.append_text(pt)
        self.projectCombo.set_active(0)

        self.projectionBox.append(self.projectionLabel )
        self.projectionBox.append(self.projectCombo )
        self.projectionBox.append(self.projectionEntry)
        self.projectionBox.append(self.projectionEntryEnd )
        self.projectionBox.append(self.projectionLabel2 )
        self.row3.append(self.projectionBox )

        self.row4 = Gtk.Box(orientation=Gtk.Orientation.HORIZONTAL, spacing=5)

        self.applyChangesButton = Gtk.Button(label="Apply Changes")
        self.applyChangesButton.connect("clicked", self.applyChanges)
        self.defaultsButton = Gtk.Button(label="Reset")
        self.defaultsButton.connect("clicked", self.resetToDefault)
        self.row4.append(self.applyChangesButton)
        self.row4.append(self.defaultsButton)

        self.vBox.append(self.row1)
        self.vBox.append(self.row2)
        self.vBox.append(self.row3)
        self.vBox.append(self.row4)

        self.eb.append(self.vBox)
        #
        # init with none, gets set in gmt plotter widget
        self.gmtPlotter = gmtWrapper
        self.loadDefaults()

    def setGMTPlotter(self, plotter):
        self.gmtPlotter = plotter

    def setPbcol(self, widget):
        rgb = self.get_color(self.gmtPlotter.pbColor, tlabel='plate boundaries')
        if rgb is not None:
            self.gmtPlotter.pbColor = rgb

    def setCoastColor(self, widget):
        rgb = self.get_color(self.gmtPlotter.pbColor, tlabel='coastlines')
        if rgb is not None:
            self.gmtPlotter.coastLineColor = rgb

    def setCoastMaskColor(self, widget):
        rgb = self.get_color(self.gmtPlotter.pbColor, tlabel='coastlines')
        if rgb is not None:
            self.gmtPlotter.coastLandColor = rgb
            self.gmtPlotter.coastSeaColor = rgb

    def setVectColor(self, widget):
        rgb = self.get_color(self.gmtPlotter.vectColor, tlabel='vectors')
        if rgb is not None:
            self.gmtPlotter.vectColor = rgb

    def get_color(self, rgb, tlabel=None):
        """ color selection dialog """
        def rgb_to_gdk_color(rgb):
            r, g, b = rgb
            return Gdk.RGBA(r, g, b, 1.0)

        def gdk_color_to_rgb(color):
            return int(color.red * 255), int(color.green * 255), int(color.blue * 255)

        title = 'Choose color'
        if tlabel:
            title += f' for {tlabel}'
        
        dialog = Gtk.ColorChooserDialog(title=title)
        dialog.set_rgba(rgb_to_gdk_color(rgb))
        response = dialog.run()
        if response == Gtk.ResponseType.OK:
            rgb = gdk_color_to_rgb(dialog.get_rgba())
        else:
            rgb = None
        dialog.destroy()
        return rgb

    def getPanel(self):
        return self.eb

    def checkColormap(self, widget):
        if self.colormapCombo.get_active_text() == "other...":
            self.colormapEntry.set_sensitive(True)
        else:
            self.colormapEntry.set_sensitive(False)

    def openRangeWindow(self, widget):
        plotter = self.gmtPlotter
        self.rangeWindow = RangeDialog(self, self.gmtPlotter.plotXmin, self.gmtPlotter.plotXmax, 
                                       self.gmtPlotter.plotYmin, self.gmtPlotter.plotYmax, self.gmtPlotterWidget)
        self.rangeWindow.show()
        self.rangeWindow.hide()

    def applyChanges(self, widget):
        cptType = str(self.colormapCombo.get_active_text())
        if cptType == "other...":
            cptType = self.colormapEntry.get_text()

        self.myPlotter = self.gmtPlotterWidget.getGMTPlotter()
        self.myPlotter.setColormapType(cptType)
        self.myPlotter.setColormapInvert(self.invertCheck.get_active())
        self.myPlotter.setDrawPlateBoundaries(self.plateBoundaryCheck.get_active())
        self.myPlotter.setDrawCoastlines(self.coastCheck.get_active())
        self.myPlotter.setMaskSea(self.coastMaskCombo.get_active() == 1)
        self.myPlotter.setMaskLand(self.coastMaskCombo.get_active() == 2)
        self.myPlotter.setAdjust(self.adjustCheck.get_active())
        self.myPlotter.setAddLabel(self.labelCheck.get_active())
        self.myPlotter.setGridLines(self.plotGridLines.get_active())
        
        projection_clon = str(self.projectionEntry.getValue())

        # plot width 
        pwidth = float(self.projectionEntryEnd.get_text())
        myprojection = GMTProjection(self.projectCombo.get_active_text()[0], projection_clon, "", pwidth, "")
        self.myPlotter.setMapProjection(myprojection)  # set the new projection
        
        # grid resolution
        gridres = float(self.gridres.get_text())
        self.myPlotter.setGridRes(gridres)

        # adjust the text projection to map width 
        myprojection = GMTProjection("X", "", "", pwidth, pwidth / 2.)
        self.myPlotter.setTextProjection(myprojection)  # set the new projection
        self.myPlotter.setColorbarPos(pwidth / 2, -pwidth * 0.042)
        self.myPlotter.setColorbarSize(pwidth / 2., pwidth * 0.036)

        self.gmtPlotterWidget.setGMTPlotter(self.myPlotter)
        self.gmtPlotterWidget.updatePlot()

    def resetToDefault(self, b):
        self.loadDefaults()
        self.gmtPlotter.plotXmin = 0
        self.gmtPlotter.plotYmin = -90
        self.gmtPlotter.plotXmax = 360
        self.gmtPlotter.plotYmax = 90

    def loadDefaults(self):
        if not self.gmtPlotter:
            print('error, no plotter instance')
            return
            
        self.invertCheck.set_active(self.gmtPlotter.colorbarInvert)
        self.projectionEntryEnd.set_text(str(self.gmtPlotter.projection.pwidth))
        self.plateBoundaryCheck.set_active(self.gmtPlotter.drawPlateBounds)
        self.coastCheck.set_active(self.gmtPlotter.drawCoastLines)
        if self.gmtPlotter.maskLand:
            self.coastMaskCombo.set_active(2)
        elif self.gmtPlotter.maskSea:
            self.coastMaskCombo.set_active(1)
        else:
            self.coastMaskCombo.set_active(0)
        self.adjustCheck.set_active(self.gmtPlotter.adjust)
        self.labelCheck.set_active(self.gmtPlotter.addLabel)
        
        if self.gmtPlotter.cptType == "haxby":
            self.colormapCombo.set_active(0)
        elif self.gmtPlotter.cptType == "gray":
            self.colormapCombo.set_active(1)
        elif self.gmtPlotter.cptType == "wysiwyg":
            self.colormapCombo.set_active(2)
        elif self.gmtPlotter.cptType == "polar":
            self.colormapCombo.set_active(3)
        elif self.gmtPlotter.cptType == "seis":
            self.colormapCombo.set_active(4)
        for i in range(len(self.projections)):
            proj = self.projections[i]
            if self.gmtPlotter.projection.type == proj[0]:
                self.projectCombo.set_active(i)
                break
        
        if self.projectCombo.get_active() != 7:
            self.projectionEntry.setValue(float(self.gmtPlotter.projection.clon))
