import gi
gi.require_version('Gtk', '4.0')
from gi.repository import Gtk, Gdk, GdkPixbuf
import os, sys
import shutil

from .plotter import *

class ImagePlotter(Plotter):
    
    def __init__(self, module, mainWindow, imageWidth, imageHeight, saveType=["png", "PNG Image"], bgColor="white"):
        """
        Constructor for an Image Plotter. This type of plotter just displays simple image plots
        
        module - the module using this plotter
        mainWindow - the MainWindow displaying this plotter
        imageWidth - the minimum width of the image
        imageHeight - the minimum height of the image
        saveType - the save type for this image (all images that are shown with this plotter should
                    be of this type). Default is a PNG image.
        bgColor - the name of the background color for the image widget (default is "white")
        """
        
        # call the superclass constructor
        Plotter.__init__(self, module, mainWindow)
        
        self.imageWidth = imageWidth
        self.imageHeight = imageHeight
        self.saveType = saveType
        self.bgColor = bgColor
        self.scaleFactor = 1.0 
 
        if self.module.baseimage:
            self.imageFile = os.path.join(self.module.directory, self.module.baseimage)
        else:
            self.imageFile = os.path.join(self.module.seatreePath, "img", "seatree.jpg")
        
    def getMainWidget(self):
        # return the Image Event Box for this image plotter
        
        self.image = Gtk.Image()
        #self.image.show()
        
        self.imageBuffer = 7
        self.imageEB = Gtk.Box()
        self.imageEB.append(self.image)
#        self.imageEB.override_background_color(Gtk.StateFlags.NORMAL, Gdk.RGBA(1, 1, 1, 1))
        self.imageEB.set_css_name("custom")
        
        # Load and apply CSS 
        css_provider = Gtk.CssProvider() 
        css_provider.load_from_data(""" 
            .custom { 
                background-color: rgba(255, 255, 255, 1); 
            } 
        """.encode('utf-8')) 
        Gtk.StyleContext.add_provider_for_display( 
            Gdk.Display.get_default(), 
            css_provider, 
            Gtk.STYLE_PROVIDER_PRIORITY_APPLICATION
        )
        self.imageEB.set_size_request(self.imageWidth*self.scaleFactor + self.imageBuffer, self.imageHeight*self.scaleFactor + self.imageBuffer)
        
        self.baseImage = self.imageFile
        
        self.displayImage(self.imageFile, default=True)
        #self.image.set_pixel_size(min(self.imageWidth, self.imageHeight))        
        print('Image pixel sizes are', self.imageWidth, self.imageHeight)
        self.imageEB.set_hexpand(True)
        self.imageEB.set_vexpand(True)
        self.image.set_hexpand(True)
        self.image.set_vexpand(True)
        #self.image.set_size_request(self.imageWidth*self.scaleFactor, self.imageHeight*self.scaleFactor)
        self.imageEB.show()
        
        return self.imageEB
    
    def getBottomPanel(self):
        # no bottom panel for default image plotter
        return None
    
    def getSaveTypes(self):
        saveTypes = []
        saveTypes.append(self.saveType)
        return saveTypes
    
    def savePlot(self, typeExtension, fileName):
        # don't check for save types since we only have 1
        if self.imageFile:
            shutil.copyfile(self.imageFile, fileName)
            return True
        else:
            return False
    
    def displayImage(self, imageFile, default=False):
        """
        This method will display an image file
        
        imageFile - the file name of the image to be displayed
        """
        self.imageFile = imageFile
        self.image.set_from_file(self.imageFile)
#        scaleFactor=1.3
#        pixbuf = GdkPixbuf.Pixbuf.new_from_file(imageFile)
#        scaled_pixbuf = pixbuf.scale_simple(self.imageWidth*scaleFactor, self.imageHeight*scaleFactor, GdkPixbuf.InterpType.BILINEAR)
#        self.image.set_from_pixbuf(scaled_pixbuf)
        if not default:
            self.mainWindow.setSaveActive(True)
