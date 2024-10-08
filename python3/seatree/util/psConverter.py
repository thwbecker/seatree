#!/usr/bin/env python3

import subprocess
import os

class PSConverter:
    
    def __init__(self, verb=0, convertPath="", psfile="", outfile="", maxWidth=650, maxHeight=400, density=90):
        self.convertPath = convertPath
        self.verb = verb
        self.maxWidth = maxWidth
        self.maxHeight = maxHeight
        self.density = density
        self.psfile = psfile
        self.outfile = outfile
        self.antialias = False
        self.use_eps2eps = False
        self.error = ""
        self.commandString = ""
    
    def calcDensity(self, maxWidth=0, maxHeight=0):
        if maxWidth != 0:
            self.maxWidth = maxWidth
        if maxHeight != 0:
            self.maxHeight = maxHeight
        if self.psfile is None:
            print('calcDensity: need psfile defined')
            return
        if not os.path.isfile(self.psfile):
            print(f'calcDensity: psfile {self.psfile} not found')
            return
        with open(self.psfile, 'r') as fp:
            for line in fp:
                index = line.find("%%BoundingBox:")
                if index > -1:
                    index += 14
                    vars = line[index:].split()
                    
                    llx = int(vars[0])
                    lly = int(vars[1])
                    urx = int(vars[2])
                    ury = int(vars[3])
                    
                    origWidth = urx - llx
                    origHeight = ury - lly
                    
                    # try scaling width first
                    density = float(72 * self.maxWidth) / float(origWidth)
                    newWidth = float(density * origWidth) / float(72)
                    newHeight = float(density * origHeight) / float(72)
                    
                    if newHeight > maxHeight:
                        # we need to scale with height
                        density = float(72 * self.maxHeight) / float(origHeight)
                    
                    self.density = int(density)
                    newWidth = float(self.density * origWidth) / float(72)
                    newHeight = float(self.density * origHeight) / float(72)

                    if self.verb > 2:
                        print(f"Calculated Density: {self.density}")
                    return self.density
        
        # if you got this far, you failed
        print("ERROR: no bounding box found!")
        return -1
    
    def command(self, command):
        if self.verb > 2:
            print(f"Command: {command}")
        proc = subprocess.Popen(command, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        self.commandString += command + "\n"
        output = proc.communicate()
        out = output[0].decode('utf-8')
        err = output[1].decode('utf-8')
        ret = proc.returncode
        if err:
            self.error += err
        if self.verb > 1 and out:
            print(out)
        if self.verb > 2 and err:
            print(err)
        return ret
    
    def convertPsToPng(self, psfile="", pngfile="", antialias=False):
        self.antialias = antialias
        if psfile:  # ps file specified here, use that
            self.psfile = psfile
        if self.psfile is None:
            print('convertPsToPng: error, ps file needs to be defined')
            print('convertPsToPng: possible cause: has GMT input already been created?')
            return

        if not pngfile:  # no png file specified, make one yourself
            if ".ps" in self.psfile:
                pngfile = self.psfile.replace(".ps", ".png")
            elif ".PS" in self.psfile:
                pngfile = self.psfile.replace(".PS", ".png")
            else:
                pngfile = self.psfile + ".png"
        
        # run a eps2eps to fix the bounding box
        if self.use_eps2eps:
            tmp_file = self.psfile.replace(".ps", ".tmp.ps")
            self.command(f"eps2eps {self.psfile} {tmp_file}")
            use_ps = tmp_file
        else:
            use_ps = self.psfile

        cmd = ""
        if self.convertPath:
            cmd += self.convertPath + os.sep
        cmd += f"convert -density {self.density} "
        if self.antialias:
            cmd += "-antialias "
        else:
            cmd += "+antialias "
        cmd += f"{use_ps} {pngfile}"
        self.command(cmd)
        
        return pngfile

    def getCommandString(self):
        temp = self.commandString
        self.commandString = ""
        return temp

if __name__ == '__main__':  # is being run from command line
    convert = PSConverter(verb=3, psfile="../py-hc/geoid.ps", maxWidth=650)
    convert.calcDensity()
    convert.convertPsToPng(pngfile="out.png")
