

class PGMImage:
	
	def __init__(self, fileName):
		self.fileName = fileName
		self.width = 0
		self.height = 0
		self.max = 255
		self.array = None
		self.loadImage()
		
	def loadImage(self):
		fp = open(self.fileName, 'r')
		line = fp.next() # P2 or comment
		while line.startswith('#'): # skip any leading comments
			line = fp.next()
		if (not line.startswith('P2')): # it's an unsupported image
			print "PGMImage:loadImage: ERROR: Unsupported Image Type!"
			return False
		line = fp.next() # skip the P2
		line = fp.next() # width height or comment
		while line.startswith('#'): # skip any comments
			line = fp.next()
		# line should be 'width height' here
		dimensions = line.split()
		self.width = int(dimensions[0]) # width of the image
		self.height = int(dimensions[1]) # height of the image
		line = fp.next() # depth of image or comment
		while line.startswith('#'): # skip any comments
			line = fp.next()
		# line should be 'depth' here
		self.max = int(line.strip())
		# the rest of the file should be the image itself!
		# iterate through it and load into array
		self.array = [ [ 0 for y in range(self.height)] for x in range(self.width) ]
		x = 0
		y = 0
		print "Loading PGM Image with Width=" + str(self.width) + ", Height=" + str(self.height) + ", and Depth=" + str(self.max)
		pixels = 0
		for line in fp:
			if (line.startswith('#')):
				continue
			tokens = line.split()
			for token in tokens:
				if (x == self.width): # we're at the end of a row
					x = 0 # move back to the start of the row
					y = y + 1 # and do the next row
				num = int(token)
				#print "x=" + str(x) + ", y=" + str(y) + ", num=" + str(num)
				self.array[x][y] = num
				x = x + 1
				pixels = pixels + 1
		print "Loaded " + str(pixels) + " pixels, x=" + str(x) + ", y=" + str(y)
	
	def getPixel(self, x, y, flip=False):
		"""
		Gets a pixel at the specified location. If flip is True, then (0,0) will be the lower left corner,
		otherwise it will be the upper leff.
		"""
		if (self.array == None):
			return None
		if (flip):
			return self.array[x][self.height - y - 1]
		else:
			return self.array[x][y]
	
	def getWidth(self):
		return self.width
	
	def getHeight(self):
		return self.height
	
	def getMax(self):
		return self.max
			
if (__name__ == '__main__'): # is being run from commmand line
	image = PGMImage("/home/kevin/workspace_seatree/SEATREE/modules/seismo/syn2d/makemodel/image2.pgm")
	print "(0,0): " + str(image.getPixel(0, 0, flip=True))