import os, sys
import math

from seatree.util.scriptRunner import ScriptRunner, ScriptResult

# find path to SEATREE root path (python)
mainpath = os.path.abspath(os.path.join(os.path.dirname(sys.argv[0]), '..'))

class GMTProjection:
    """
    GMT projection structure
    """
    def __init__(self, type, clon, clat, pwidth, pheight):
        self.type = type    # type of projection, e.g. H
        self.clon = clon    # central longitude
        self.clat = clat    # central latitude
        self.pwidth = pwidth    # width of plot
        self.pheight = pheight

class GMTWrapper:
    """
    GMT Wrapping object. Assumes that GMT is installed and is in the system path.
    """
    
    def __init__(self, verb=0, path="", tmpn=os.path.join(os.sep, "tmp", "gmtPlotter"), 
                 runDir="", awk="awk", grdfile=None):
        """
        Initializes HCWrapper and sets HC args
        
        Arguments:
        verb -- verbosity level from 0-3
        """
        #---------------
        # set defaults
        #---------------
        
        # awk command
        self.awk = awk
        
        # verbosity level
        self.verb = verb
        
        # gmt path
        self.gmtpath = path
        
        # temp path
        self.tmpn = tmpn
        
        # old PAPER_MEDIA selection
        self.oldMedia = ""
        
        # custom Colormap file
        self.cptFile = ""
    
        # Colormap type
        self.cptType = "haxby"

        self.last_grd_vals = []
            
        # post script file
        self.psFile = ""
        self.psState = 0 # ps state: 0 - closed or nonexistent, 1 - open
        
        # grid range
        self.gridXmin = 0.
        self.gridXmax = 360.
        self.gridYmin = -90.
        self.gridYmax = 90.
        
        # plotting range
        self.plotXmin = 0.
        self.plotXmax = 360.
        self.plotYmin = -90.
        self.plotYmax = 90.
        
        # plot offset
        self.plotXOffset = 0
        self.plotYOffset = 0
        
        # boundary annotation and tick intervals
        self.plotBoundAnnotation = ""
        
        # map projection
        self.projection = GMTProjection("H", 180, "", 7, "")
        
        # default to portrait mode
        self.portrait = 1
    
        # grid lines
        self.drawGridLines = False

        # coastline variables
        self.drawCoastLines = True
        self.maskLand = False
        self.maskSea = False
        self.coastMaskArea = 70000
        self.coastResolution = "c"
        self.coastWidth = 2
        self.coastLineColor = [100, 100, 100]
        self.coastLandColor = [128, 128, 128]
        self.coastSeaColor = [200, 200, 200]

        # colorbar variables
        self.colorbarN = 50
        self.colorbarXpos = "3.5"
        self.colorbarYpos = "-.3"
        self.colorbarLength = "3"
        self.colorbarWidth = ".25"
        self.colorbarHorizontal = 1
        self.colorbarTriangles = 1
        self.colorbarInterval = 50.0
        self.colorbarInvert = False
        self.colorbarUnits = ""
        
        # text variables
        self.textClip = 0
        self.textXmin = 0.
        self.textXmax = 1.
        self.textYmin = 0.
        self.textYmax = 1.
        self.textProjection = GMTProjection("X", "", "", 7, 3.5)
        
        # vector variables
        self.vectConvertToAngles = True # when using grdvector, adjust azimuth depending on map projection
        
        self.vectArrowWidth = "0.025i"
        self.vectHeadLength = "0.12i"
        self.vectHeadWidth = "0.045i"
        self.vectScaleShorterThanSize = .2
        self.vectScale = 3  # in cm/yr 
        self.vectScale0 = 4  # when fixed
        self.vectColor = [255, 165, 0]  # RGB style
        self.vectXmin = 0
        self.vectXmax = 350
        self.vectYmin = -85
        self.vectYmax = 85
        self.vectWidth = .5
        
        # adjust plot scales automatically?
        self.adjust = True

        # add a label to the plot
        self.addLabel = True

        # GMT style annotation, without the -B
        self.annotation = ""

        self.gridres = 1

        # plate boundary variables
        self.drawPlateBounds = False
        self.pbFile = os.path.join(mainpath, "..", "data", "common", "nuvel.360.xy")
        self.pbLinewidth = 5
        self.pbColor = [0, 0, 128]

        self.error = ""

        self.commandString = ""
        
        # move the bounding box to the front of the file
        self.BBFront = True
        # modify the bounding box, this also implies that you are moving it to the front
        self.modifyBB = False
        # user specified bounding box
        # defaults to HC's geoid BB
        self.bbLLX = 71
        self.bbLLY = 0
        self.bbURX = 575
        self.bbURY = 342
        
        # value for locations without data in xyz2grd
        self.noDataValue = float('nan')
        
        # force xyz2grd pixel registration
        self.forcePixelRegistration = False
        
        # create a script runner instance
        self.scriptRunner = ScriptRunner(workingDir=runDir)
        
        self.setup_GMT()

    def setVerbosity(self, verb):
        self.verb = verb
    
    def setGMTPath(self, gmtpath):
        self.gmtpath = gmtpath
    
    def setRunDir(self, runDir):
        '''
        Sets the directory that GMT commands should be run from.
        IMPORTANT: all file paths should be either absolute or relative to
        this run directory!
        '''
        self.scriptRunner.setWorkingDir(runDir)
    
    def runGMT(self, command):
        if self.verb > 2:
            print(f"Command: {command}")
        
        self.commandString += command + "\n"
        
        result = self.scriptRunner.runScript(command)
        out = result.getStandardOutput()
        err = result.getStandardError()
        ret = result.getReturnValue()
        
        if err:
            self.error += err
        if self.verb > 1 and out:
            print(out)
        if self.verb > 2 and err:
            print(err)
            print(f"return code: {ret}")
        return [ret, out, err]

    def clearCommandString(self):
        self.commandString = ""

    def getCommandString(self):
        return self.commandString
    
    def setup_GMT(self):
        """ Loads GMT configuration file and backs up Paper Media default
        before changing it """
        command = ""
        if self.gmtpath:
            command += os.path.join(self.gmtpath, "")
        command += "gmtdefaults -L"
        result = self.scriptRunner.runScript(command)
        out = result.getStandardOutput()
        self.gmt4 = True
        for line in out.splitlines():
            if 'PAPER_MEDIA' in line:
                line = line.split('=')[1].strip()
                self.oldMedia = line
                if self.verb > 2:
                    print(f"Old Paper Media Config: {line}")
            if 'GMT-SYSTEM' in line:
                if float(line.split()[2][0]) > 3:
                    if self.verb > 2:
                        print('detected GMT version > 4')
                    self.gmt4 = True
                else:
                    self.gmt4 = False
        
        command = ""
        if self.gmtpath:
            command += os.path.join(self.gmtpath, "")
        command += "gmtset PAPER_MEDIA letter+"
        self.runGMT(command)
    
    def cleanup(self):
        """ Restores old GMT Paper Media configuration as backed up by
        setup_GMT. """
        if "letter+" not in self.oldMedia:
            if self.verb > 2:
                print("Restoring old GMT paper media config")
            command = ""
            if self.gmtpath:
                command += os.path.join(self.gmtpath, "")
            command += f"gmtset PAPER_MEDIA {self.oldMedia}"
            self.runGMT(command)
    
    def setAwk(self, awk):
        self.awk = awk

    def makeCPT(self, z0, z1, dz, outFile, setBackgroundMax=False):
        """
        Makes a Colormap file with the given z values, color table, and out file
        
        Arguments:
        z0 -- minimum z-value
        z1 -- maximum z-value
        dz -- z step size
        colorTable -- selects master color table
        outFile -- selects the output cpt file
        setBackgroundMax -- selects the GMT -D option
        """
        command = ""
        if self.gmtpath:
            command += self.gmtpath + os.sep
        command += "makecpt "
        if self.colorbarInvert:
            command += "-I "
        if setBackgroundMax:
            command += "-D "
        # remove any comments for use of cpt file in larry 
        command += f"-T{z0}/{z1}/{dz} -C{self.cptType}"
        if not setBackgroundMax:
            command += f" | {self.awk} '{{if((substr($1,1,1) != \"#\") && (NF == 8))print($0)}}'"
        
        command += f' > {outFile}'
        self.runGMT(command)
        if os.path.exists(outFile):
            self.cptFile = outFile
        else:
            print("Error creating CPT file!")
    
    def setCPTFile(self, cptFile):
        """ Sets custom CPT file for use by GMT """
        self.cptFile = cptFile
    
    def getCPTFile(self):
        """ Returns custom CPT file in use by GMT """
        return self.cptFile

    def setColormapInvert(self, invert):
        self.colorbarInvert = invert
    
    def getColormapInvert(self):
        return self.colorbarInvert

    def setColormapType(self, type):
        self.cptType = type
    
    def setGridRes(self, gridres):
        self.gridres = gridres

    def getColormapType(self):
        return self.cptType
    
    def setGridRange(self, xmin, xmax, ymin, ymax):
        """
        Sets GMT plotting range
        
        Defaults:
        xmin = 0
        xmax = 360
        ymin = -90
        ymax = 90
        """
        self.gridXmin = xmin
        self.gridXmax = xmax
        self.gridYmin = ymin
        self.gridYmax = ymax
    
    def detectGridRange(self, inc, xyzFile):
        """
        Detects the grid range with minmax
        """
        command = ""
        if self.gmtpath:
            command += self.gmtpath + os.sep
        
        command += f"minmax -I{inc} {xyzFile}"
        result = self.runGMT(command)
        val = result[0]
        out = result[1]
        err = result[2]
        if val == 0:
            out = out.replace("\n", "").replace("-R", "")
            split = out.split("/")
            self.gridXmin = float(split[0])
            self.gridXmax = float(split[1])
            self.gridYmin = float(split[2])
            self.gridYmax = float(split[3])
            
            print(f"Detected region: xmin={self.gridXmin} xmax={self.gridXmax} ymin={self.gridYmin} ymax={self.gridYmax}")
        else:
            print("Error detecting region!")
    
    def getGridRange(self):
        """
        Returns the grid range in an array:
            [xmin, xmax, ymin, ymax]
        """
        return [self.gridXmin, self.gridXmax, self.gridYmin, self.gridYmax]
    
    def setNoDataValue(self, value):
        """
        Sets the -N flag in xyz2grd (for no flag, give it float('nan')):
        
        No data.  Set nodes with no input xyz triplet to this value [Default is NaN].  For z-tables, this  option
        is used to replace z-values that equal nodata with NaN.
        """
        self.noDataValue = value
    
    def setForcePixelRegistration(self, force):
        """
        Tells xyz2grd to force pixel registration (Default is grid registration).
        """
        self.forcePixelRegistration = force

    def spatialToNetCDF(self, inc, inpipe, outfile, interpolate, xmin=None, xmax=None, ymin=None, ymax=None, verbose=False):
        """
        Convert a file from spatial to NetCDF GRD
        
        Arguments:
        inc -- grid spacing in degrees
        inpipe -- shell command to pipe spherical harmonics file to sh_syn
        outFile -- selects the output cpt file
        interpolate -- 0: use xyz2grd assuming regularly spaced input
                       1: use surface to interpolate
        xmin, xmax, ymin, ymax: optional grid boundaries. if not set, will use gridXmin etc.
        """
        xmin = xmin if xmin is not None else self.gridXmin
        xmax = xmax if xmax is not None else self.gridXmax
        ymin = ymin if ymin is not None else self.gridYmin
        ymax = ymax if ymax is not None else self.gridYmax
        
        command = inpipe + " | "
        if self.gmtpath:
            command += self.gmtpath + os.sep
        is_global = (ymin == -90 and ymax == 90) and ((xmin == 0 and xmax == 360) or (xmin == -180 and xmax == 180))
        
        noDataString = f" -N{self.noDataValue} " if self.noDataValue != float('nan') else ""
        verboseString = " -V " if verbose else ""
        forcePixelRegistrationString = " -F " if self.forcePixelRegistration else ""
        gflags = "-Lg" if is_global and not self.gmt4 else "-fg" if is_global else ""
        
        if interpolate:
            command += f"surface {gflags} {self.getRangeString(xmin, xmax, ymin, ymax, self.projection, False)} -I{inc}{noDataString}{verboseString}{forcePixelRegistrationString} -G{outfile}"
        else:
            command += f"xyz2grd {gflags} {self.getRangeString(xmin, xmax, ymin, ymax, self.projection, False)} -I{inc}{noDataString}{verboseString}{forcePixelRegistrationString} -G{outfile}"
            
        self.runGMT(command)
        if os.path.exists(outfile):
            self.grdFile = outfile
        else:
            print("Error creating GRD file!")
    
    def getGRDFile(self):
        """ Returns custom CPT file in use by GMT """
        return self.grdFile
    
    def setPlotRange(self, xmin, xmax, ymin, ymax):
        """
        Sets GMT plotting range
        
        Defaults:
        xmin = 0
        xmax = 360
        ymin = -90
        ymax = 90
        """
        self.plotXmin = xmin
        self.plotXmax = xmax
        self.plotYmin = ymin
        self.plotYmax = ymax
    
    def setPortraitMode(self, portrait):
        """
        Sets plotter to portrait mode
        
        Arguments:
        portrait -- boolean, 1 for portrait, 0 for landscape
        """
        self.portrait = portrait
    
    def getRangeString(self, Xmin, Xmax, Ymin, Ymax, projection=None, adjust_for_plotting=False):
        """
        Will determine the -R range string

        If adjust_for_plotting is set, will check if Mercator and leave out poles if selected
        """
        if adjust_for_plotting and projection:
            if projection.type == 'M':
                if Ymin == -90 and Ymax == 90:
                    Ymin, Ymax = -75, 75

        if self.gmt4:
            if Ymin == -90 and Ymax == 90:
                if Xmin == 0 and Xmax == 360:
                    return '-Rg'
                if Xmin == -180 and Xmax == 180:
                    return '-Rd'
        return f"-R{Xmin}/{Xmax}/{Ymin}/{Ymax}"

    def getProjectionString(self, projection):
        """
        Form the -J projection string, adding arguments deciding on which one chosen

        Arguments:
        projection -- GMTProjection structure, init with 
        """
        if projection.type in ["H", "J", "R", "N", "W", "Q", "Kf"]:
            return f"-J{projection.type}{projection.clon}/{projection.pwidth}i"
        elif projection.type == "X":
            reg_string = f"-J{projection.type}{projection.pwidth}i"
            if projection.pheight:
                reg_string += f"/{projection.pheight}i"
            return reg_string
        elif projection.type == "M":
            return f"-J{projection.type}{projection.pwidth}i"
        else:
            print(f"Projection type {projection.type} undefined in getProjectionString")
            return f"-J{projection.type}{projection.clon}/{projection.pwidth}i"

    def setPSFile(self, outFile):
        self.psFile = outFile
    
    def getColorString(self, R, G, B):
        return f"{R}/{G}/{B}"

    def getLineString(self, linewidth, R, G, B):
        return f"-W{linewidth},{self.getColorString(R, G, B)}"  # new GMT
    
    def setModifyBoundingBox(self, modify):
        self.modifyBB = modify
    
    def setBoundingBox(self, llx, lly, urx, ury):
        self.modifyBB = True
        self.bbLLX = llx
        self.bbLLY = lly
        self.bbURX = urx
        self.bbURY = ury
    
    def initPSFile(self, outFile, xOff=0, yOff=1.25, basemap=False):
        self.psFile = outFile
        
        basemap = basemap and self.plotBoundAnnotation

        portraitString = " -P" if self.portrait else ""
        
        command = "" if basemap else "echo 1000 1000 | "
        if self.gmtpath:
            command += self.gmtpath + os.sep
        offset = ""
        if xOff != 0:
            offset += f" -X{xOff}i"
        if yOff != 0:
            offset += f" -Y{yOff}i"
        offset += " "
        
        if basemap:
            command += f"psbasemap {self.getRangeString(self.plotXmin, self.plotXmax, self.plotYmin, self.plotYmax, self.projection, True)} {self.getProjectionString(self.projection)} {offset}-B{self.plotBoundAnnotation}{portraitString} -K > {self.psFile}"
        else:
            command += f"psxy {self.getRangeString(self.plotXmin, self.plotXmax, self.plotYmin, self.plotYmax, self.projection, True)} {self.getProjectionString(self.projection)} {offset}{portraitString} -Sa1i -K > {self.psFile}"
        self.runGMT(command)
    
    def setPlotOffset(self, xOff, yOff):
        self.plotXOffset = xOff
        self.plotYOffset = yOff
    
    def closePSFile(self):
        """
        End a GMT postscript plot by adding a point outside the range 
        """
        command = "echo 1000 1000 | "
        if self.gmtpath:
            command += self.gmtpath + os.sep
        command += f"psxy {self.getRangeString(self.plotXmin, self.plotXmax, self.plotYmin, self.plotYmax, self.projection, True)} {self.getProjectionString(self.projection)} -Sa1i -O >> {self.psFile}"
        self.runGMT(command)
        
        if self.BBFront or self.modifyBB:
            with open(self.psFile, 'r', errors='ignore') as fp:
                lines = []
                atend = -1
                box = -1
                for line in fp:
                    if "%%BoundingBox:" in line and "atend" in line:
                        atend = len(lines)
                    if "%%BoundingBox:" in line and "atend" not in line:
                        box = len(lines)
                    lines.append(line)
            if self.modifyBB and box > -1 and atend < 0:  # we need to change it, but it's already at the top
                atend = box
            if box > -1 and atend > -1:  # we need to move it to the top
                if self.modifyBB:
                    lines[atend] = f"%%BoundingBox: {self.bbLLX} {self.bbLLY} {self.bbURX} {self.bbURY}\n"
                else:
                    lines[atend] = lines[box]
                lines.pop(box)

                with open(self.psFile, 'w', errors='ignore') as fp:
                    for line in lines:
                        fp.write(line)
    
    def setBoundaryAnnotation(self, annotation):
        """
        This sets the boundary annotation and tick flag, it should NOT contain the -B flag
        """
        self.plotBoundAnnotation = annotation
    
    def createImageFromGrid(self, grdfile):
        """
        Make a postscript plot from a GRD file
        """
        command = ""
        if self.gmtpath:
            command += self.gmtpath + os.sep
        command += f"grdimage {grdfile} "
        command += f"{self.getRangeString(self.plotXmin, self.plotXmax, self.plotYmin, self.plotYmax, self.projection, True)} "
        if self.plotXOffset:
            command += f"-X{self.plotXOffset} "
        if self.plotYOffset:
            command += f"-Y{self.plotYOffset} "
        if self.plotBoundAnnotation:
            command += f"-B{self.plotBoundAnnotation} "
        command += f"{self.getProjectionString(self.projection)} -C{self.cptFile} -O -K >> {self.psFile}"

        self.runGMT(command)
        
    def plotAnomalyBoxesWCPT(self, gmtFile):
        """
        Add main anomaly boxes to postscript file
        """
        command = ""
        if self.gmtpath:
            command += self.gmtpath + os.sep
        command += f"psxy {gmtFile} {self.getProjectionString(self.projection)} {self.getRangeString(self.plotXmin, self.plotXmax, self.plotYmin, self.plotYmax, self.projection, True)} -m -L -K -O -C{self.cptFile} >> {self.psFile}"
        self.runGMT(command)

    def plotAnomalyBoxes(self, gmtFile):
        """
        Add main anomaly boxes to postscript file
        """
        command = ""
        if self.gmtpath:
            command += self.gmtpath + os.sep
        command += f"psxy {gmtFile} {self.getProjectionString(self.projection)} {self.getRangeString(self.plotXmin, self.plotXmax, self.plotYmin, self.plotYmax, self.projection, True)} -A -K -O -M >> {self.psFile}"
        self.runGMT(command)

    def drawPlateBoundaries(self):
        """
        Plot the plate boundaries
        """
        self.plotPolygon(self.pbFile, self.pbLinewidth, self.pbColor[0], self.pbColor[1], self.pbColor[2])

    def setDrawPlateBoundaries(self, draw):
        self.drawPlateBounds = draw
        
    def isDrawPlateBoundaries(self):
        return self.drawPlateBounds
        
    def setDrawCoastlines(self, draw):
        self.drawCoastLines = draw
        
    def isDrawCoastLines(self):
        return self.drawCoastLines
    
    def setMaskLand(self, draw):
        self.maskLand = draw
    
    def isMaskLand(self):
        return self.maskLand
    
    def setMaskSea(self, draw):
        self.maskSea = draw
        
    def isMaskSea(self):
        return self.maskSea

    def setAdjust(self, adjust):
        self.adjust = adjust

    def setAddLabel(self, addLabel):
        self.addLabel = addLabel
    
    def setGridLines(self, draw):
        """ switch on grid lines, overriding other annotation """
        if draw:
            self.setAnnotation("g60/g30")
        else:
            self.setAnnotation("")

    def setAnnotation(self, ann_string):
        """ set the GMT annotation string, without the -B"""
        self.annotation = ann_string

    def setCoastlineMaskArea(self, area):
        """ Sets coastline masking area -Axxx type where xxx is the argument """
        self.coastMaskArea = area
    
    def setCoastlineResolution(self, resolution):
        """ Sets coastline resolution, give i, h, f, l, or c"""
        self.coastResolution = resolution
    
    def setCoastlineWidth(self, coastWidth):
        """ Sets coastline width """
        self.coastWidth = coastWidth

    def drawCoastline(self, drawCoasts=None, maskSea=None, maskLand=None):
        """
        Draws coastlines to a ps file and adds annotation
        drawCoasts - Boolean for drawing coastlines, or None for default
        maskSea - Boolean for masking oceans, or None for default
        maskLand - Boolean for masking land, or None for default
        """
        drawCoasts = drawCoasts if drawCoasts is not None else self.drawCoastLines
        maskSea = maskSea if maskSea is not None else self.maskSea
        maskLand = maskLand if maskLand is not None else self.maskLand

        if not (drawCoasts or maskSea or maskLand):
            return

        shape_string = ""
        if drawCoasts:
            shape_string = self.getLineString(self.coastWidth, self.coastLineColor[0], self.coastLineColor[1], self.coastLineColor[2]) + ' '
        if maskSea:
            shape_string += "-S" + self.getColorString(self.coastSeaColor[0], self.coastSeaColor[1], self.coastSeaColor[2]) + ' '
        if maskLand:
            shape_string += "-G" + self.getColorString(self.coastLandColor[0], self.coastLandColor[1], self.coastLandColor[2]) + ' '

        command = ""
        if self.gmtpath:
            command += self.gmtpath + os.sep
        add_string = f" -B{self.annotation} " if self.annotation else " "
        command += f"pscoast {self.getProjectionString(self.projection)} {self.getRangeString(self.plotXmin, self.plotXmax, self.plotYmin, self.plotYmax, self.projection, True)} -A{self.coastMaskArea} -D{self.coastResolution}{add_string} {shape_string} -O -K >> {self.psFile}"
        self.runGMT(command)
    
    def setColorbarN(self, n):
        """ Sets effective dots-per-inch for colorbar """
        self.colorbarN = n
    
    def setColorbarPos(self, xpos, ypos):
        """ Sets colorbar position """
        self.colorbarXpos = xpos
        self.colorbarYpos = ypos
    
    def setColorbarSize(self, length, width):
        """ Sets colorbar size """
        self.colorbarLength = length
        self.colorbarWidth = width
    
    def setColorbarHorizontal(self, horiz):
        """
        Sets colorbar to horizontal mode
        
        Arguments:
        horizontal -- boolean, 1 for horizontal, 0 for vertical
        """
        self.colorbarHorizontal = horiz
    
    def setColorbarTriangles(self, triangles):
        """ Adds triangles to ends of colorbar """
        self.colorbarTriangles = triangles
    
    def setColorbarInterval(self, interval):
        """ Sets colorbar interval (float) """
        self.colorbarInterval = interval

    def setColorbarUnits(self, units):
        """ Sets colorbar units (string) """
        self.colorbarUnits = units
    
    def drawColorbar(self):
        """
        Draw a colorbar to a ps file
        """
        dString = f"-D{self.colorbarXpos}i/{self.colorbarYpos}i/{self.colorbarLength}i/{self.colorbarWidth}i"
        if self.colorbarHorizontal:
            dString += "h"
        
        eString = " -E" if self.colorbarTriangles else ""
        
        command = ""
        if self.gmtpath:
            command += self.gmtpath + os.sep
        command += f"psscale -N{self.colorbarN} {dString}{eString} -C{self.cptFile} -B{self.colorbarInterval}/:\"{self.colorbarUnits}\": -O -K >> {self.psFile}"
        self.runGMT(command)
    
    def grdMath(self, infiles, operations, outfile):
        """
        Does math on NetCDF GRD files
        
        Arguments:
        infiles -- array of input files
        operations -- array of operations
        outfile -- output file
        """
        command = ""
        if self.gmtpath:
            command += self.gmtpath + os.sep
        command += "grdmath"
        
        for file in infiles:
            command += f" {file}"
        
        for op in operations:
            command += f" {op}"
        
        command += f" = {outfile}"
        self.runGMT(command)
    
    def grdNiceCmpRange(self, min, max, cutoff=1.0, ncol=21, nlabel=5, symmetric=True):
        """ 
        Determine the minimum and max of a grdfile, reduced it
        by fraction cutoff (full range for 1.0), and return
        bounds with ncol subdivisions for the actual colorbar,
        and with a "nice" spacing close to nlabel for the
        labeling. 

        If symmetric is set, will make sure the colorbar is
        symmetric around the absolute max

        i.e. returns xmin_red, xmax_red, dx, dx_nice
        """
        if max - min < 1e-5:
            min, max = -1, 1

        if min < 0 and max > 0 and symmetric:
            if -min > max:
                max = -min
            else:
                min = -max
        
        tr = [cutoff * min, cutoff * max]
        range = tr[1] - tr[0]
        tr.append(range / ncol)
        dx = range / nlabel
        dx = (10 ** int(math.log10(dx) + 0.5)) * 2
        while range / dx < 3.0:
            dx /= 2.0

        tr.append(dx)
        return tr

    def grdMinMaxMean(self, infile, geo=True):
        """ 
        Computes the min, max, and mean of a given NetCDF GRD file 

        Returns min, max, mean

        If geo = True, then take sphericity into account in even dx/dy
        spaced grids
        """
        if not infile or not os.path.exists(infile):
            print(f'grdMinMaxMean: file {infile} does not exist')
            return

        command = ""
        if self.gmtpath:
            command += self.gmtpath + os.sep
        command += f"grdinfo -C -L2 {infile} | {self.awk} '{{print($2,$3,$4,$5,$6,$7,$8,$9)}}'"
        if self.verb > 2:
            print(f"Command: {command}")
        result = self.scriptRunner.runScript(command)
        out = result.getStandardOutput().split()
        reg = self.getRangeString(float(out[0]), float(out[1]), float(out[2]), float(out[3]), self.projection, False)
        inc = f'-I{float(out[6])}/{float(out[7])}'
        min_val, max_val = float(out[4]), float(out[5])

        if geo:
            self.grdMath([reg + ' ' + inc], 'Y', infile + '.tmp.lat')
            self.grdMath([infile], [f'ISNAN 1 SUB ABS 0 NAN {infile}.tmp.lat COSD MUL '], infile + '.tmp.costheta')
            sumw, n = self.grdSum(infile + '.tmp.costheta')
            self.grdMath([infile], [f'{infile}.tmp.costheta MUL'], infile + '.tmp.scaled')
            sum_val, n = self.grdSum(infile + '.tmp.scaled')
            mean = sum_val / sumw
            rmstring = f'rm -f {infile}.tmp.lat {infile}.tmp.scaled {infile}.tmp.costheta'
            self.runGMT(rmstring)
        else:
            sum_val, n = self.grdSum(infile)
            mean = sum_val / n if n != 0 else None

        return min_val, max_val, mean

    def grdSum(self, infile):
        """ 
        Sum all non-NaN entries of a grdfile 

        Returns sum and number_of_entries
        """
        command = ""
        if self.gmtpath:
            command += self.gmtpath + os.sep
        command += f"grd2xyz -S -Z {infile}"
        with os.popen(command) as fp:
            out = fp.read()
        outs = out.split()
        total_sum = sum(float(val) for val in outs)
        return total_sum, len(outs)
    
    def plotPolygon(self, polygonfile, linewidth, R, G, B):
        """
        Plots GMT -M style polygons
        
        Arguments:
        polygonfile: filename of polygon file
        linewidth: linewidth
        R, G, B: RGB color for polygon
        """
        command = ""
        if self.gmtpath:
            command += self.gmtpath + os.sep
        command += f"psxy {polygonfile} -M {self.getProjectionString(self.projection)} {self.getRangeString(self.plotXmin, self.plotXmax, self.plotYmin, self.plotYmax, self.projection, True)} {self.getLineString(linewidth, R, G, B)} -O -K >> {self.psFile}"
        self.runGMT(command)
    
    def plotXY(self, xyFile, colorName="", colorR=-1, colorG=-1, colorB=-1, plotSymbols=False, symbol="", symbolSize=0):
        """
        Plots GMT XY plots
        
        Arguments:
        xyFile: filename of xy file
        colorName: name of the color to plot
        colorR: R component of color (0-255) if no colorName supplied
        colorG: G component of color (0-255) if no colorName supplied
        colorB: B component of color (0-255) if no colorName supplied
        plotSymbols: boolean to plot symbols
        symbol: symbol to plot (empty if specified in data)
        symbolSize: size of symbol to plot (empty if specified in data)
        """
        command = ""
        if self.gmtpath:
            command += self.gmtpath + os.sep
        command += f"psxy {xyFile}"
        if plotSymbols:
            command += " -S"
            if symbol:
                command += symbol
            if symbolSize > 0:
                command += f"{symbolSize}i"
        
        command += f" {self.getProjectionString(self.projection)} {self.getRangeString(self.plotXmin, self.plotXmax, self.plotYmin, self.plotYmax, self.projection, True)}"
        if colorName:
            command += f" -G{colorName}"
        elif colorR >= 0 and colorG >= 0 and colorB >= 0:
            command += f" -G{colorR}/{colorG}/{colorB}"
        command += f" -O -K >> {self.psFile}"
        self.runGMT(command)

    def setTextClipping(self, textClip):
        """
        Sets text plotter to clipping mode
        
        Arguments:
        textClip -- boolean, 1 for clipping, 0 for no clipping (default 0)
        """
        self.textClip = textClip
    
    def setTextRegion(self, xmin, xmax, ymin, ymax):
        """ Sets text plot region """
        self.textXmin = xmin
        self.textXmax = xmax
        self.textYmin = ymin
        self.textYmax = ymax
    
    def setMapProjection(self, projection):
        """ Sets custom GMT projection """
        self.projection = projection

    def setTextProjection(self, projection):
        """ Sets text plot projection """
        self.textProjection = projection

    def plotText(self, text):
        """ Plots text onto ps file """
        command = f"echo {text} | "
        if self.gmtpath:
            command += self.gmtpath + os.sep
        
        command += "pstext"
        if not self.textClip:
            command += " -N"
        command += f" {self.getRangeString(self.textXmin, self.textXmax, self.textYmin, self.textYmax, self.projection, True)} {self.getProjectionString(self.textProjection)} -O -K >> {self.psFile}"
        self.runGMT(command)
    
    def setVectConvertToAngles(self, convert):
        """
        Sets vector plotter to convert means to angles based off of projection
        
        Arguments:
        convert -- boolean, 1 for conversion, 0 for no conversion (default 1)
        """
        self.vectConvertToAngles = convert
    
    def setVectArrowSize(self, arrowWidth, headLength, headWidth):
        """ Sets arrow size to given arrow width, head length, and head width respectively """
        self.vectArrowWidth = arrowWidth
        self.vectHeadLength = headLength
        self.vectHeadWidth = headWidth
    
    def setVectScaleShorterThanSize(self, size):
        """ Sets size for which all vectors that are shorter are scaled by length / size """
        self.vectScaleShorterThanSize = size
    
    def setVectScale(self, scale):
        """ Sets vector scale factor """
        self.vectScale = scale
    
    def setVectColor(self, R, G, B):
        """ Sets vector to given RGB color """
        self.vectColor = [R, G, B]
    
    def setVectRegion(self, xmin, xmax, ymin, ymax):
        """ Sets vector plot to region """
        self.vectXmin = xmin
        self.vectXmax = xmax
        self.vectYmin = ymin
        self.vectYmax = ymax
    
    def setVectOutlineWidth(self, width):
        """ Sets vector width """
        self.vectWidth = width
    
    def plotVectors(self, inX, inY):
        """ 
        Plots the given X and Y NetCDF GRD files 
        
        Returns mean vector length
        """
        abs_grdfile = self.tmpn + ".abs.grd"
        self.grdMath([inX, inY], ["R2", "SQRT"], abs_grdfile)
        min_val, max_val, mean_vec_length = self.grdMinMaxMean(abs_grdfile, geo=True)
        
        command = ""
        if self.gmtpath:
            command += self.gmtpath + os.sep
        command += f"grdvector {inX} {inY}"
        if self.vectConvertToAngles:
            command += " -T"
        else:
            print('plotVectors: WARNING: not converting vectors to map projected directions')
        command += f" {self.getRangeString(self.plotXmin, self.plotXmax, self.plotYmin, self.plotYmax, self.projection, True)} {self.getProjectionString(self.projection)} -Q{self.vectArrowWidth}/{self.vectHeadLength}/{self.vectHeadWidth}"
        if self.vectScaleShorterThanSize:
            command += f"n{self.vectScaleShorterThanSize}i"

        if self.adjust:
            self.vectScale = mean_vec_length * 4.5
        else:
            self.vectScale = self.vectScale0 * 4.5

        if self.vectScale != 0:
            command += f" -S{self.vectScale}i"
        command += f" -G{self.getColorString(self.vectColor[0], self.vectColor[1], self.vectColor[2])}"
        command += f" {self.getLineString(self.vectWidth, 50, 50, 50)} -O -K >> {self.psFile}"
        self.runGMT(command)
        return mean_vec_length

