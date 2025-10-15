import xml.dom.minidom

class ReadXml:

	def __init__(self, fileName, Element = ""):
		if(not Element == ""):
			self.root = Element
		else:
			doc = xml.dom.minidom.parse(fileName)
			self.root = doc.firstChild

	def getRootChildren(self):
		i = 0
		self.rootChildren = []
		for e in self.root.childNodes:
			if (e.nodeType == e.ELEMENT_NODE):
				self.rootChildren.append(e)
				i = i + 1
		return self.rootChildren

	def getNumElements(self, theElement = ""):
		if(theElement == ""):
			theElement = self.root
		counter= 0
		for e in theElement.childNodes:
			counter = counter + 1
		return counter

	def getNodeLocalName(self, nodeNumber):
		return self.root.childNodes[nodeNumber].localName

	def getNode(self, nodeNumber):
		return self.root.childNodes[nodeNumber]

	def getNodeText(self, nodeNumber):
		string = self.root.childNodes[nodeNumber].childNodes[0].nodeValue.strip()
		return string.strip()
