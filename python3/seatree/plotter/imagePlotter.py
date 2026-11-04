import gi
gi.require_version('Gtk', '4.0')
from gi.repository import Gtk, Gdk, GdkPixbuf
import os, sys
import shutil
import platform

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
        """
        Return the widget holding the plot image.

        The image is shown in a Gtk.Picture, which scales its content to
        the space it is given while keeping the aspect ratio, so the
        plot follows the size of the window and of the pane divider.
        The PNG files the plotters produce have a higher resolution than
        a typical window, so scaling them down keeps them sharp.
        imageWidth and imageHeight are used as the minimum size.
        """
        self.image = Gtk.Picture()
        self.image.set_can_shrink(True)
        if hasattr(self.image, "set_content_fit"):
            # GTK 4.8 and later
            self.image.set_content_fit(Gtk.ContentFit.CONTAIN)
        else:
            self.image.set_keep_aspect_ratio(True)
        self.image.set_hexpand(True)
        self.image.set_vexpand(True)
        self.image.set_size_request(self.imageWidth, self.imageHeight)

        self.imageBuffer = 7
        self.imageEB = Gtk.Box()
        self.imageEB.append(self.image)
        self.imageEB.set_css_name("custom")
        
        # white background behind the plot
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
        self.imageEB.set_hexpand(True)
        self.imageEB.set_vexpand(True)
        
        self.baseImage = self.imageFile
        self.displayImage(self.imageFile, default=True)
        
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
        print(f"displayImage called with: {imageFile}")

        if not imageFile:
            print("ERROR: imageFile is None or empty!")
            return

        if not os.path.isfile(imageFile):
            print(f"ERROR: Image file does not exist: {imageFile}")
            return

        file_size = os.path.getsize(imageFile)
        print(f"Image file exists: {imageFile}, size: {file_size} bytes")

        if file_size == 0:
            print(f"ERROR: Image file is empty: {imageFile}")
            return

        try:
            self.imageFile = imageFile
            self.image.set_filename(self.imageFile)
            print(f"Successfully loaded image: {imageFile}")
        except Exception as e:
            print(f"ERROR loading image {imageFile}: {e}")
            return

        if not default:
            self.mainWindow.setSaveActive(True)
