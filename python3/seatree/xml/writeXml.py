import xml.dom.minidom

class WriteXml:

    def __init__(self, name="root", fileLocation="doc.xml"):
        self.fileName = fileLocation
        self.xml_doc = xml.dom.minidom.Document()
        self.root = self.xml_doc.createElement(name)
        self.xml_doc.appendChild(self.root)

    def addToRoot(self, element):
        self.root.appendChild(element)

    def addNode(self, name):
        newNode = self.xml_doc.createElement(name)
        self.root.appendChild(newNode)
        return newNode
    
    def addSubNode(self, parent, name):
        newNode = self.xml_doc.createElement(name)
        parent.appendChild(newNode)
        return newNode

    def addAttribute(self, node, attrName, attr):
        node.setAttribute(attrName, attr)

    def addText(self, node, text):
        theText = self.xml_doc.createTextNode(text)
        node.appendChild(theText)
        return node

    def setFileName(self, name):
        self.fileName = name

    def getRoot(self):
        return self.root

    def writeToXml(self):
        with open(self.fileName, "w", encoding="utf-8") as file_object:
            self.xml_doc.writexml(file_object, indent="    ", addindent=" ", newl="\n", encoding="UTF-8")
