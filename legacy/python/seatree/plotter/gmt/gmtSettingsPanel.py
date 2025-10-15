import pygtk
pygtk.require('2.0')
import gtk
import seatree.gui.util.guiUtils as guiUtils
from seatree.gmt.gmtWrapper import *
from seatree.gui.util.rangeDialog import RangeDialog

class GMTSettingsPanel:
	
	def __init__(self, gmtPlotterWidget, gmtWrapper, width=200):


		self.tooltips = gtk.Tooltips()



		self.width = width
		
		self.gmtPlotterWidget = gmtPlotterWidget
		self.eb = gtk.EventBox()
		self.eb.modify_bg(gtk.STATE_NORMAL, gtk.gdk.color_parse("gray"))
		self.eb.set_size_request(self.width, 200)

		self.vBox = gtk.VBox(homogeneous=True, spacing=1)	
		self.row1 = gtk.HBox(homogeneous=False, spacing=1)
		
		# Colormap Type
		self.colormapCombo = gtk.combo_box_new_text()
		self.colormapCombo.append_text("roma")
		self.colormapCombo.append_text("grayC")
		self.colormapCombo.append_text("haxby")
		self.colormapCombo.append_text("vik")
		self.colormapCombo.append_text("seis")
		self.colormapCombo.append_text("spectral")
		self.colormapCombo.append_text("other...")
		self.colormapCombo.set_active(0)


		self.colormapCombo.connect("changed", self.checkColormap)
		self.tooltips.set_tip(self.colormapCombo, 'select some of the GMT colormaps, or type your preferred colormap in the text field', tip_private=None)


		self.colormapEntry = gtk.Entry()
		self.colormapEntry.set_width_chars(5)
		self.colormapEntry.set_sensitive(False)

		self.colormapLabel = gtk.Label("Colormap Type:")
		self.colormapBox = gtk.HBox(homogeneous=False, spacing=5)
		self.colormapBox.pack_start(self.colormapLabel)
		
		self.invertCheck = gtk.CheckButton(label="Invert")
		self.tooltips.set_tip(self.invertCheck, 'flip the colormap around, e.g. from red->blue to blue->red', tip_private=None)
		self.invertCheckBox = gtk.HBox(homogeneous=False, spacing=0)
		self.invertCheckBox.pack_start(self.invertCheck)

		self.adjustCheck = gtk.CheckButton(label="Auto scale")
		self.tooltips.set_tip(self.adjustCheck, 'adjust the colormap automatically to range of scalars plotted in background, else will used fixed range', tip_private=None)
		self.adjustCheck.set_active(True)
		self.adjustCheckBox = gtk.HBox(homogeneous = False, spacing = 0)
		self.adjustCheckBox.pack_start(self.adjustCheck)

		self.labelCheck = gtk.CheckButton(label="Labels")
		self.labelCheck.set_active(True)
		self.tooltips.set_tip(self.labelCheck, 'Add a colorbar and labels to the maps.', tip_private=None)
		self.labelCheckBox = gtk.HBox(homogeneous = False, spacing = 0)
		self.labelCheckBox.pack_start(self.labelCheck)

		self.plateBoundaryCheck = gtk.CheckButton(label="Pbound")
		self.tooltips.set_tip(self.plateBoundaryCheck, 'plot NUVEL1 major tectonic plate boundaries', tip_private=None)
		self.plateBoundaryCheckBox = gtk.HBox(homogeneous = False, spacing = 0)
		self.plateBoundaryCheckBox.pack_start(self.plateBoundaryCheck)
		
		self.coastCheck = gtk.CheckButton(label="Coastlines")
		self.tooltips.set_tip(self.coastCheck, 'Plot coastlines', tip_private=None)
		self.coastCheckBox = gtk.HBox(homogeneous = False, spacing = 0)
		self.coastCheckBox.pack_start(self.coastCheck)
		
		self.coastMaskCombo = gtk.combo_box_new_text()
		self.coastMaskCombo.append_text("No Mask")
		self.coastMaskCombo.append_text("Mask Sea")
		self.coastMaskCombo.append_text("Mask Land")
		self.coastMaskCombo.set_active(0)

		self.plotGridLines =    gtk.CheckButton(label="Grid lines")
		self.tooltips.set_tip(self.plotGridLines, 'plot grid lines on top of map', tip_private=None)
		self.plotGridLinesBox = gtk.HBox(homogeneous = False, spacing = 0)
		self.plotGridLinesBox.pack_start(self.plotGridLines)

		self.colorSymLabel = gtk.Label("Colors:")

		self.coastColorButton = gtk.Button("Coast")
		self.coastColorButton.connect("clicked", self.setCoastColor)
		self.tooltips.set_tip(self.coastColorButton, 'change color for coast lines', tip_private=None)
		
		self.coastMaskColorButton = gtk.Button("Mask")
		self.coastMaskColorButton.connect("clicked", self.setCoastMaskColor)
		self.tooltips.set_tip(self.coastMaskColorButton, 'change color for coast mask', tip_private=None)

		self.pbcolButton = gtk.Button("Pbound")
		self.pbcolButton.connect("clicked", self.setPbcol)
		self.tooltips.set_tip(self.pbcolButton, 'change color for plate boundaries', tip_private=None)

		self.vecColButton = gtk.Button("Vectors")
		self.vecColButton.connect("clicked", self.setVectColor)
		self.tooltips.set_tip(self.vecColButton, 'change color for vectors', tip_private=None)
		


         	# grid resolution
		self.gridres = gtk.Entry()
		self.gridres.set_width_chars(4)
		self.gridres.set_max_length(4)
		self.gridres.set_text("1.0")
		self.gridresLabel = gtk.Label("resolution")
		self.tooltips.set_tip(self.gridres, 'set the resolution (in degrees) of the background grid (if applicable). 0.25...2 are useful for global plots', tip_private=None)


		#Range Settings
#		self.rangeLabel = gtk.Label("Lat/Lon Range:")
#		self.range1Box = gtk.HBox(homogeneous = False, spacing = 0)
#		self.range1Box.pack_start(self.rangeLabel)
		
		self.rangeButton = gtk.Button("Update Lat/Lon Range")
		self.tooltips.set_tip(self.rangeButton, 
				      'change the geographic region to be mapped (not that useful for global plots)', tip_private=None)
		self.rangeButton.connect("clicked", self.openRangeWindow)


		self.row1.pack_start(self.colormapBox, expand=False)
		self.row1.pack_start(self.colormapCombo, expand=False)
		self.row1.pack_start(self.colormapEntry, expand=False)
		self.row1.pack_start(self.invertCheckBox, expand=False)
		self.row1.pack_start(self.adjustCheckBox,expand=False)
		self.row1.pack_start(self.labelCheckBox,expand=False)
		self.row1.pack_start(self.rangeButton, expand = False)
		

		self.row2 = gtk.HBox(homogeneous=False, spacing=5)

		self.applyChangesButton = gtk.Button("Apply Changes")
		self.applyChangesButton.connect("clicked", self.applyChanges)
		
#		self.row2.pack_start(self.range1Box, expand = False)
		self.row2.pack_start(self.plotGridLinesBox, expand = False)
		self.row2.pack_start(self.coastCheckBox, expand=False)
		self.row2.pack_start(self.coastMaskCombo, expand=False)
		self.row2.pack_start(self.plateBoundaryCheckBox, expand=False)
		self.row2.pack_start(self.colorSymLabel,expand=False)
		self.row2.pack_start(self.pbcolButton,expand=False)
		self.row2.pack_start(self.coastColorButton,expand=False)
#		self.row2.pack_start(self.coastMaskColorButton,expand=False)
		self.row2.pack_start(self.vecColButton,expand=False)
		
		self.row2.pack_start(self.gridres,expand=False)
		self.row2.pack_start(self.gridresLabel,expand=False)
		
		self.row3 = gtk.HBox(homogeneous=True, spacing=5)
		
		# Projection

		self.projectionLabel = gtk.Label("Projection:")
		self.projectionEntry = guiUtils.RangeSelectionBox(initial=170, min = 0, max = 359, digits = 0, buttons = True)
		# width of projection
		self.projectionEntryEnd = gtk.Entry()
		self.projectionEntryEnd.set_width_chars(2)
		self.projectionEntryEnd.set_max_length(2)
		self.projectionEntryEnd.set_text("7")
		self.projectionLabel2 = gtk.Label("PS width (inches)")

		self.projectionBox = gtk.HBox(homogeneous=False, spacing=50)

		# projection type
		self.projectCombo = gtk.combo_box_new_text()
		self.projections = [ 'W Mollweide',  'H Hammer', 'N Robinson' ,\
						'Q Equi-Cyl' , 'J Miller', 'M Mercator', 'Kf Eckert-IV', 'X Linear' ]

		for pt in self.projections:
			self.projectCombo.append_text(pt)
		self.projectCombo.set_active(0)

		self.projectionBox.pack_start(self.projectionLabel,expand = False)
		self.projectionBox.pack_start(self.projectCombo,expand = False)
		self.projectionBox.pack_start(self.projectionEntry,expand = True)
		self.projectionBox.pack_start(self.projectionEntryEnd,expand = False)
		self.projectionBox.pack_start(self.projectionLabel2,expand = False)
		self.row3.pack_start(self.projectionBox, expand=False)

	
		self.row4 = gtk.HBox(homogeneous = False, spacing = 5)

		self.applyChangesButton = gtk.Button("Apply Changes")
		self.applyChangesButton.connect("clicked", self.applyChanges)
		self.defaultsButton = gtk.Button("Reset")
		self.defaultsButton.connect("clicked", self.resetToDefault)
		self.row4.pack_start(self.applyChangesButton, expand = False)
		self.row4.pack_start(self.defaultsButton, expand = False)

		self.vBox.pack_start(self.row1, expand=False, fill=False)
		self.vBox.pack_start(self.row2, expand=False, fill=False)
		self.vBox.pack_start(self.row3, expand=True, fill=False)
		self.vBox.pack_start(self.row4, expand=False, fill=False)
		
		self.eb.add(self.vBox)
		#
		# init with none, gets set in gmt plotter widget
		self.gmtPlotter = gmtWrapper
		self.loadDefaults()

	def setGMTPlotter(self, plotter):
		self.gmtPlotter = plotter

	def setPbcol(self,widget):
		rgb = self.get_color(self.gmtPlotter.pbColor,tlabel = 'plate boundaries')
		if rgb != None:
			self.gmtPlotter.pbColor = rgb
		return
	
	def setCoastColor(self,widget):
		rgb = self.get_color(self.gmtPlotter.pbColor,tlabel = 'coastlines')
		if rgb != None:
			self.gmtPlotter.coastLineColor = rgb
		return
	
	def setCoastMaskColor(self,widget):
		rgb = self.get_color(self.gmtPlotter.pbColor,tlabel = 'coastlines')
		if rgb != None:
			self.gmtPlotter.coastLandColor = rgb
			self.gmtPlotter.coastSeaColor = rgb
		return

	def setVectColor(self,widgt):
		rgb = self.get_color(self.gmtPlotter.vectColor,tlabel = 'vectors')
		if rgb != None:
			self.gmtPlotter.vectColor = rgb
		return

	def get_color(self,rgb,tlabel = None):
		""" color selection dialog """
		def rgb_to_gdk_color(rgb):
			r,g,b = rgb
			color = gtk.gdk.Color(int(r*256), int(g*256), int(b*256))
			return color
		def gdk_color_to_rgb(color): # return ints
			return int(color.red/256.0), int(color.green/256.0), int(color.blue/256.0)
		title = 'Choose color'
		if tlabel != None:
			title += ' for ' + tlabel
		dialog = gtk.ColorSelectionDialog('Choose color')

		colorsel = dialog.colorsel
		color = rgb_to_gdk_color(rgb)
		colorsel.set_previous_color(color)
		colorsel.set_current_color(color)
		colorsel.set_has_palette(True)
		response = dialog.run()
		if response == gtk.RESPONSE_OK:
			rgb = gdk_color_to_rgb(colorsel.get_current_color())
		else:
			rgb = None
		dialog.destroy()
		return rgb

	def getPanel(self):
		return self.eb

	def checkColormap(self, widget):
		if(self.colormapCombo.get_active_text() == "other..."):
			self.colormapEntry.set_sensitive(True)
		else:
			self.colormapEntry.set_sensitive(False)

	def openRangeWindow(self, widget):
		plotter = self.gmtPlotter
		self.rangeWindow = RangeDialog(self, self.gmtPlotter.plotXmin, self.gmtPlotter.plotXmax, \
						self.gmtPlotter.plotYmin, self.gmtPlotter.plotYmax, self.gmtPlotterWidget)
		self.rangeWindow.show()
		self.rangeWindow.hide()

	def applyChanges(self, widget):
		
		cptType = str(self.colormapCombo.get_active_text())
		if(cptType == "other..."):
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
		myprojection = GMTProjection(self.projectCombo.get_active_text()[0],projection_clon,"",pwidth,"")
		self.myPlotter.setMapProjection(myprojection) # set the new projection
		
		# grid resolution
		gridres = float(self.gridres.get_text())
		self.myPlotter.setGridRes(gridres)

		# adjust the text projection to map width 
		myprojection = GMTProjection("X","","",pwidth,pwidth/2.)
		self.myPlotter.setTextProjection(myprojection) # set the new projection
		self.myPlotter.setColorbarPos(pwidth/2, -pwidth*0.042)
		self.myPlotter.setColorbarSize(pwidth/2., pwidth*0.036)

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
			print 'error, no plotter instance'
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
		
		if(self.gmtPlotter.cptType == "haxby"):
			self.colormapCombo.set_active(0)
		elif(self.gmtPlotter.cptType == "gray"):
			self.colormapCombo.set_active(1)
		elif(self.gmtPlotter.cptType == "wysiwyg"):
			self.colormapCombo.set_active(2)
		elif(self.gmtPlotter.cptType == "polar"):
			self.colormapCombo.set_active(3)
		elif(self.gmtPlotter.cptType == "seis"):
			self.colormapCombo.set_active(4)
		for i in xrange(len(self.projections)):
			proj = self.projections[i]
			if self.gmtPlotter.projection.type == proj[0]:
				self.projectCombo.set_active(i)
				break
#		if(self.gmtPlotter.projection.type == "H"):
#			self.projectCombo.set_active(0)
#		elif(self.gmtPlotter.projection.type == "N"):
#			self.projectCombo.set_active(1)
#		elif(self.gmtPlotter.projection.type == "Q"):
#			self.projectCombo.set_active(2)
#		elif(self.gmtPlotter.projection.type == "R"):
#			self.projectCombo.set_active(3)
#		elif(self.gmtPlotter.projection.type == "J"):
#			self.projectCombo.set_active(4)
#		elif(self.gmtPlotter.projection.type == "X"):
#			self.projectCombo.set_active(7)
		
		if (self.projectCombo.get_active() != 7):
			self.projectionEntry.setValue(float(self.gmtPlotter.projection.clon))
