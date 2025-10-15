#
# manipulate a HC type data profile
#
import math, copy

import sys,os

import seatree.plotter.matPlotLib.matPlotLibPlotter as matPlotLibPlotter

class DataPointChooser:
	"""
	"""
	
	# initialize class
	def __init__(self, figure, xy_diag, numX, numY, xtol = None, ytol = None):
		
		if xtol == None:
			xtol = float(numX) * 0.01
		
		if ytol == None:
			ytol = float(numY) * 0.01
		
		self.numX = numX
		self.numY = numY
		self.xy_diag = xy_diag
		
		self.xmin = 0
		self.xmax = self.numX
		
		self.ymin = 0
		self.ymax = self.numY
		
		self.figure = figure
		
		self.axis = self.figure.add_subplot(111)
		
		self.sources = []
		self.receivers = []
		
		self.origSources = []
		self.origReceivers = []
		
		self.xtol = xtol
		self.ytol = ytol
		
		#
		self.verbose = True	 # progress output
		
		# start a plot
		self.redraw_plot()
		self.figure.set_visible(True);

	def distance(self, x1, x2, y1, y2):
		"""
		return the distance between two points
		"""
		return(math.sqrt( (x1 - x2)**2 + (y1 - y2)**2 ))


	def __call__(self, event):
		#
		# get the x and y data coords
		#
		x, y = event.xdata, event.ydata

		if event.inaxes:
			print 'generic call?: data coords', x,y
	
	def getClosestPointInTolerance(self, x, y, list):
		minDist = -1
		i=0
		ret = (-1, None, minDist, list)
		for item in list:
			px,py = item[0], item[1]
			
			if abs(x - px) <= self.xtol and abs(y - py) <= self.ytol:
				dist = self.distance(x, px, y, py)
				if minDist == -1 or dist < minDist:
					minDist = dist
					ret = (i, item, minDist, list)
			
			i += 1
		
		return ret

	def on_click(self, event):
		# 
		# mouse click event
		# 
		if event.inaxes:		# within region
			# data coordinates
			x, y = event.xdata, event.ydata
			if self.axis == event.inaxes:
				#
				# look for close points
				#
				closest_source = self.getClosestPointInTolerance(x, y, self.sources)
				closest_reciever = self.getClosestPointInTolerance(x, y, self.receivers)
				
				closest = closest_source
				if closest[0] >= 0:
					if closest_reciever[0] >= 0:
						if closest_reciever[2] < closest_source[2]:
							closest = closest_reciever
				elif closest_reciever[0] >= 0:
					closest = closest_reciever
				else:
					closest = None
				button = event.button
				if button == 1 or button == 3: # center mouse click: add point to list
					if not closest:
						if self.verbose:
							print 'adding data point %7.2f, %7.2f ' % (x, y)
						if button == 3:
							self.sources.append((x, y))
						else:
							self.receivers.append((x, y))
						self.redraw_plot()
					else:
						if self.verbose:
							print 'there is already a point at %7.2f, %7.2f ' % (x, y)
				elif button == 2:
					if closest:
						val = closest[1]
						theList = closest[3]
						theList.remove(val)
						if self.verbose:
							print 'removing data point %7.2f, %7.2f ' % val
						self.redraw_plot()
					else:
						if self.verbose:
							print 'did not find data close to click  %7.2f, %7.2f' % (x,y)
	
	def plotScatterData(self, x, y, type=None, color='b', size=30, globalWidth=0.2, linewidths=None):
		if type == None:
			type = self.CIRCLE
		
		if globalWidth != None and not linewidths:
			linewidths = []
			for i in range(0, len(x)):
				linewidths.append(globalWidth)
		self.axis.scatter(x, y, s=size, c=color, marker=type, linewidths=linewidths)

	def redraw_plot(self):	  # refresh the plot
		"""

		redraw a plot

		"""
		# get the figure handle
		#p.axes(self.axis)
		self.axis.clear()
		
		size=30
		if len(self.sources) > 0:
			x = []
			y = []
			for source in self.sources:
				x.append(source[0])
				y.append(source[1])
			color = "r"
			type = matPlotLibPlotter.PENTAGRAM
			self.plotScatterData(x, y, type, color)
		
		if len(self.receivers) > 0:
			x = []
			y = []
			for reciever in self.receivers:
				x.append(reciever[0])
				y.append(reciever[1])
			color = "b"
			type = matPlotLibPlotter.TRIANGLE_DOWN
			self.plotScatterData(x, y, type, color)

		# fix the axes
		self.axis.set_ylim([self.ymin, self.ymax])
		self.axis.set_xlim([self.xmin,self.xmax])
		
		self.xy_diag.updateButtons(self.sources, self.receivers)
  
# what is the renderer?
#		self.axis.draw('GTKAgg')
		self.xy_diag.redraw()

	def reset_data(self):
		if self.verbose:
			print 'resetting to original data'
		self.sources = copy.copy(self.origSources)
		self.receivers = copy.copy(self.origReceivers)
		self.redraw_plot()
	
	def setOriginalData(self, sources=[], receivers=[]):
		self.origSources = sources
		self.origReceivers = receivers
	
	def getSources(self):
		return self.sources
	
	def getReceivers(self):
		return self.receivers