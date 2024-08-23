import os
import subprocess

class HCWrapper:
    
    def __init__(self, verb, dm, dt, dfac, dsf, use_dsf, tbc, platevelf, premfile, vf, spd, hcpath, computedir, tmpn, lkludge, awk="awk"):
        """
        Initializes HCWrapper and sets HC args
        
        Arguments:
        verb -- verbosity level from 0-3
        dm -- path to density model
        dt -- SH type
        dfac -- density scaling factor
        dsf -- density scaling factor file with depth-dependent density scaling, overrides dfac setting if use_dsf is set 
        use_dsf -- actually use the density scaling factor file
        tbc -- surface boundary condition: 0: free slip 1: no slip 2: plate velocities
        platevelf -- velocity boundary conditions file
        premfile -- prem model file
        vf -- viscosity file
        spd -- scale with prem density
        hcpath -- path to hc
        computedir -- directory to store computations
        tmpn -- path/prefix for temp files
        awk -- 
        lkludge -- minimum l from which to apply solver kludge 
        
        """
        self.verb = verb
        self.dm = os.path.abspath(dm)
        self.dt = dt
        self.dfac = dfac
        self.dsf = dsf
        self.use_dsf = use_dsf
        self.platevelf = platevelf
        self.premfile = premfile
        self.tbc = tbc
        self.lkludge = lkludge
        self.vf = os.path.abspath(vf)
        self.hcpath = hcpath
        self.tmpn = tmpn
        self.tmpdir = os.path.dirname(tmpn)
        self.commandString = ""
        self.spd = spd
        self.awk = awk
        
        self.setComputeDir(computedir)
        
        if verb == 0:
            self.hc_verb = ""
        elif verb == 1:
            self.hc_verb = "-v"
        elif verb == 2:
            self.hc_verb = "-vv"
        elif verb == 3:
            self.hc_verb = "-vvv"
        
        self.error = ""
    
    def setComputeDir(self, computedir):
        self.compdir = computedir
        if computedir:
            self.geoidFile = os.path.abspath(self.compdir + os.sep + "geoid.ab")
            self.velFile = os.path.abspath(self.compdir + os.sep + "vel.sol.bin")
            self.rtracFile = os.path.abspath(self.compdir + os.sep + "rtrac.sol.bin")
        else:
            self.geoidFile = os.path.abspath("geoid.ab")
            self.velFile = os.path.abspath("vel.sol.bin")
            self.rtracFile = os.path.abspath("rtrac.sol.bin")
    
    def getHCPath(self):
        return self.hcpath
    
    def run_hc(self, command):
        """ Run HC with given arguments """
        if self.hcpath:
            command = self.hcpath + os.sep + command
        if self.compdir:
            command = "cd " + self.compdir + "; " + command
        if self.verb > 2:
            print("Command: " + command)
        proc = subprocess.Popen(command, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        
        self.commandString += command + "\n"
        
        output = proc.communicate()
        out = output[0].decode()
        err = output[1].decode()
        ret = proc.returncode
        if err:
            self.error += err
        if self.verb > 1 and out:
            print(out)
        if self.verb > 2 and err:
            print(err)
        return ret

    def compute_vel(self):
        """ Computes velocities and geoid with hc with values from __init__ """
        if self.verb > 0:
            print("Computing Velocities...")
        
        if os.path.exists(self.velFile):
            if self.verb > 2:
                print("Velocity file exists, deleting...")
            os.remove(self.velFile)
        if os.path.exists(self.geoidFile):
            if self.verb > 2:
                print("Geoid file exists, deleting...")
            os.remove(self.geoidFile)
        ret = self.run_hc(self.get_common_hc_command())
        if not os.path.exists(self.velFile) or not os.path.exists(self.geoidFile):
            print("Error calculating Velocities/Geoid!")
        
        return ret == 0
    
    def compute_rtrac(self):
        """ Computes radial tractions with hc with values from __init__ """
        
        if self.verb > 0:
            print("Computing Radial Tractions...")
        if os.path.exists(self.rtracFile):
            if self.verb > 2:
                print("Radial Tractions file exists, deleting...")
            os.remove(self.rtracFile)
        ret = self.run_hc(self.get_common_hc_command() + " -rtrac ")
        if not os.path.exists(self.rtracFile):
            print("Error calculating Radial Tractions!")
        
        return ret == 0

    def get_common_hc_command(self):
        """ common hc command line strings """

        # top boundary conditions
        common_string = "hc "
        common_string += " -cbckl " + str(self.lkludge) + ' ' # l kludge 
        common_string += " -prem " + self.premfile
            
        if self.tbc == 0: # free slip
            common_string +=  " -fs "
        elif self.tbc == 1: # no slip
            common_string += " -ns "
        elif self.tbc == 2: # plates
            if self.platevelf == "":
                print('error: plate velocity file needs to be specified')
                exit()
            common_string += " -pvel " + self.platevelf + ' '	# specify plate velocity file
        else:
            print('error ', self.tbc, ' tbc code undefined ')
            exit()
        
        # density anomaly input
        if not self.use_dsf:	# fixed scaling factor
            common_string += " -ds " + str(self.dfac) + ' '
        else:		# depth-dependent in file 
            if self.dsf is None:
                print('hcWrapper: error, need dsf file')
            common_string += " -dsf " + str(self.dsf) + ' '
        
        # scale with PREM?
        if not self.spd:
            common_string += " -dnp "
        
        # density anomaly file and type
        common_string += " -dens " + self.dm + " -" + self.dt + ' '
        
        # viscosity file
        common_string += " -vf " + self.vf + ' '
        
        # verbosity
        common_string += self.hc_verb

        return common_string

    def extract_sh_layer(self, infile, mode, format, column, outfile):
        """ Extract given column of all depths to given file 

        mode 1: scalar, e.g. x_r
        mode 2: x_pol x_tor
        mode 3: x_r x_pol x_tor

        """
        command = "hc_extract_sh_layer " + infile + " " + str(mode) + " " + str(format) + \
            " | " + self.awk + " '{print($" + str(column) + ")}' > " + outfile
        self.run_hc(command)
    
    def extract_sh_layer_to_spherical(self, infile, layer, mode, format, inc, w, e, s, n, outfile):
        """ Extract given column of all depths to given file and convert to spatial

        mode 1: scalar, e.g. x_r
        mode 2: x_pol x_tor
        mode 3: x_r x_pol x_tor
        mode 4: x_pol x_tor, but use only poloidal component, toroidal set to zero 
        mode 5: x_pol x_tor, but use only toroidal component, poloidal set to zero 

        inc: spacing
        w,e,s,n: geographic range
        
        """
        if mode == 4:
            sstring = self.awk + " '{if(NR==1)print($0);else print($1,$2,0,0)}' | "
            mode = 2
        elif mode == 5:
            sstring = self.awk + " '{if(NR==1)print($0);else print(0,0,$3,$4)}' | "
            mode = 2
        else:
            sstring = ''
        command = "hc_extract_sh_layer " + infile + " " + str(layer) + " " + str(mode) + " " + \
            str(format) + " | " + sstring

        if self.hcpath:
            command += self.hcpath + os.sep
        command += self.get_sh_string(w, e, s, n, inc) + " > " + outfile
        self.run_hc(command)

    def get_sh_string(self, w, e, s, n, res):
        return "sh_syn 0 0 " + str(w) + " " + str(e) + " " + str(s) + " " + str(n) + " " + str(res) + " 2> /dev/null "
