import xml.dom.minidom

class ReadXml:

    def __init__(self, fileName, Element=""):
        if Element:
            self.root = Element
        else:
            doc = xml.dom.minidom.parse(fileName)
            self.root = doc.firstChild

    def getRootChildren(self):
        self.rootChildren = [e for e in self.root.childNodes if e.nodeType == e.ELEMENT_NODE]
        return self.rootChildren

    def getNumElements(self, theElement=""):
        if not theElement:
            theElement = self.root
        return sum(1 for e in theElement.childNodes)

    def getNodeLocalName(self, nodeNumber):
        return self.root.childNodes[nodeNumber].localName

    def getNode(self, nodeNumber):
        return self.root.childNodes[nodeNumber]

    def getNodeText(self, nodeNumber):
        return self.root.childNodes[nodeNumber].childNodes[0].nodeValue.strip()
