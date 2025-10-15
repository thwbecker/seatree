#
# manipulate a HC type data profile
#
import math

import sys,os
# find path to SEATREE root path (python)
path = os.path.abspath(os.path.dirname(sys.argv[0]) + os.sep + '..' + os.sep)

class ManipulateXYData:
    """

    handles x y data that specifies HC-type profiles
    
    use like:

    >>>

    mp = ManipulateXYData(filename,mode)
    p.connect('button_press_event', mp.on_click)
    p.connect('button_release_event', mp.on_release)

    <<<

    INPUT

    filename: data name
    
    mode: 1: data are viscosity
          2: data are density scaling factors


    xtol amd ytol are relative tolerances

    inspired by http://www.scipy.org/Cookbook/Matplotlib/Interactive_Plotting
    
    """
    
    # initialize class
    def __init__(self, filename, mode, figure, xy_diag, tmpn="", \
                 xtol = None, ytol = None, figcount = 1, data_folder = None):
    	self.use = False;
    	self.saveToFile=False
    	self.tmpn = tmpn
    	self.outfile = ""
    	self.xy_diag = xy_diag
        if xtol == None:
            xtol = 0.1
        if ytol == None:
            ytol = 0.1
        self.xtol = xtol
        self.ytol = ytol
        
        if data_folder == None:
            data_folder = path + os.sep + '..' + os.sep + 'data' + os.sep + 'hc' + os.sep
        self.data_folder = data_folder

        self.visc_norm = 1.e21  # some constants
        self.radius_km = 6371.
        self.cmb_km = 2891.


        self.figure = figure
        self.axis = self.figure.add_subplot(111)

        #
        # 1: viscosity 
        # 2: density 
        self.plot_mode = mode
        if mode == 1:
            self.xmin,self.xmax = 1e-2,1e3
        else:
            self.xmin,self.xmax = -0.1,0.4
        #
        # read and convert data 
        self.datax, self.datay = self.read_data(filename,mode)

        self.datax0 = self.datax # copy for restore
        self.datay0 = self.datay

        self.moving = -1        # if point are being move


        self.zlabels = [300,660,1750] # for plot
        #
        self.verbose = True     # progress output

        #
        self.use_sh_bg = True   # show background best-fit models
        self.bg_init = False
        
        self.xl,self.bgxl,self.yl,self.bgyl =  [], [], [], []

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

    def on_click(self, event):
        # 
        # mouse click event
        # 
        if event.inaxes:        # within region
            # data coordinates
            x, y = event.xdata, event.ydata
            if self.axis == event.inaxes:
                #
                # look for close points
                #
                cps = []
                i=0
                data_cont = zip(self.datax,self.datay) # reformat
                for xd,yd in data_cont: # compute distances for
                    # those points that are
                    # within range compute
                    # tolerance
                    if xd == 0:         # compute tolerances
                        xts = 0.5
                    else:
                        xts = abs(xd)*self.xtol
                    if yd == 0:
                        yts = 0.5
                    else:
                        yts = abs(yd)*self.ytol
                    #
                    # if close, compute distance
                    if  (abs(x-xd) < xts) and  (abs(y-yd) < yts) :
                        cps.append( (self.distance(xd,x,yd,y),xd,yd,i) )
                    i=i+1
                if cps:             # if we found some point close enough, sort them and use the closest
                    cps.sort()
                    dist, xd, yd, i = cps[0] # closest

                if event.button == 2: # center mouse click: add point to list
                    if not cps or dist > 1:
                        if self.verbose:
                            print 'adding data point %7.2f, %7.2f ' % (x, y)
                        self.datax.append(x)
                        self.datay.append(y)
                        self.redraw_plot()
                    else:
                        if self.verbose:
                            print 'there is already a point at %7.2f, %7.2f ' % (x, y)
                else:
                    # 
                    # left or right
                    # 
                    if cps:
                        if event.button == 1: 
                            # left mouse button, move this data point
                            self.moving = i
                            if self.verbose:
                                print 'moving data point %5i ' % i, 'from %7.2f, %7.2f ' % (xd, yd)
                        elif event.button == 3: 
                    # right click: removing this data point
                            if self.verbose:
                                print 'removing data point %7.2f, %7.2f ' % (self.datax[i],self.datay[i])
                            del self.datax[i]
                            del self.datay[i]
                            self.redraw_plot()
                    else:
                        if self.verbose:
                            print 'did not find data close to click  %7.2f, %7.2f' % (x,y)


    def on_release(self, event):
        # mouse has been released
        if event.inaxes:
            xd, yd = event.xdata, event.ydata
            if self.axis == event.inaxes:
                if self.moving > -1: # are we actually moving a point?
                    if self.verbose:
                        print 'assigning %7.2f, %7.2f to data point %5i' % (xd, yd, self.moving)
                    i=0;xn,yn=[],[]
                    data_cont=zip(self.datax,self.datay)
                    # this could be dealt with smarter
                    self.datax, self.datay = [], []
                    for x,y in data_cont: # replace the self.moving-th point with the current location
                        if i==self.moving:
                            self.datax.append(xd);self.datay.append(yd)
                        else:
                            self.datax.append(x);self.datay.append(y)
                        i+=1
                    self.redraw_plot()
                    self.moving = -1 # reset

    def redraw_plot(self):      # refresh the plot
        """

        redraw a plot

        """
        self.datax, self.datay = self.sortlevels(self.datax,self.datay) # sort data and get layer entries, adjust max and min

        # get the figure handle
        #p.axes(self.axis)
        self.axis.clear()

  # adjust range
        if self.plot_mode == 1:  # viscosity 
            if min(self.datax) < self.xmin:
                self.xmin /= 10.
            if max(self.datax) > self.xmax:
                self.xmax *= 10.
        else:                   # density or 
            if min(self.datax) < self.xmin:
                self.xmin = self.datax *0.8
            if max(self.datax) > self.xmax:
                self.xmax = self.datay *1.2
        
        # make layer plot
        self.xl, self.yl = self.get_layer_from_xy_data(self.datax,self.datay)


        if self.plot_mode == 1:  # viscosity 
            self.add_pmantle_ornaments()
            self.axis.semilogx(self.xl,self.yl,linewidth=3,color='red') # plot layers
            self.axis.semilogx(self.datax,self.datay,'o') # plot actual profile
            self.axis.set_xlabel('viscosity [1e21 Pas]')

        elif self.plot_mode == 2: # density scaling factor
            self.add_pmantle_ornaments()
            self.axis.plot(self.xl,self.yl,linewidth=3,color='blue')
            self.axis.plot(self.datax,self.datay,'o') # plot actual profile

            self.axis.set_title('left mouse: move center: add right: remove point')
            self.axis.set_xlabel('scale factor')

        # fix the axes
        self.axis.set_ylim([-self.cmb_km, 0])
        self.axis.set_xlim([self.xmin,self.xmax])
  
# what is the renderer?
#        self.axis.draw('GTKAgg')
        self.xy_diag.redraw()
  

    def add_pmantle_ornaments(self):
        """
        add ornaments typical for the earth's mantle to the plot
        """
        self.axis.grid(True)
        self.axis.set_title('left mouse: move center: add right: remove point')
        self.axis.set_ylabel('depth [km]')
        if self.plot_mode == 1:
            xoff = self.xmin*2.5
            uselogx = True
        else:
            xoff = 0.025*(self.xmax-self.xmin)
            uselogx = False
        x = [self.xmin, self.xmax]
        for z in self.zlabels: # add a few labels
            y=[-z,-z];
            self.axis.text(self.xmin+xoff,-z+10.,str(z)+' km',fontstyle='italic')
            if uselogx:
                self.axis.semilogx(x,y,linewidth=2,color='black',linestyle='dashed')
            else:
                self.axis.plot(x,y,linewidth=2,color='black',linestyle='dashed')
        # add some best fit models?
        if self.use_sh_bg:
            if not self.bg_init:
                # read in file name 
                if self.plot_mode == 1:
                    filename = self.data_folder + os.sep + 'viscosity'+ os.sep + 'visc.sh08'
                else:
                    filename = self.data_folder + os.sep + 'dscale'+ os.sep + 'dscale_sh08.dat'
                if not os.path.exists(filename):
                    print 'add_pmantle_ornaments: reference ', filename , ' not found'
                else:
                    # initialize files
                    tmpx, tmpy = self.read_data(filename,self.plot_mode)
                    self.bgxl, self.bgyl = self.get_layer_from_xy_data(tmpx,tmpy)
                    self.bg_init = True
            if self.bg_init:
                if uselogx:
                    self.axis.semilogx(self.bgxl,self.bgyl,linewidth=1,color='black',linestyle='-')
                else:
                    self.axis.plot(self.bgxl,self.bgyl,linewidth=1,color='black',linestyle='-')
  

    def reset_data(self):
        if self.verbose:
            print 'resetting to original data'
        self.datax = self.datax0
        self.datay = self.datay0
        self.redraw_plot()

    def sortlevels(self,datax,datay):  
        """
        
        sort through a list of weirdly formatted viscosity file
        values and add data point to make a plot look nice

        also, assign layer plot data 

        """
        # sort the z and eta vectors
        data = []
        for zl, el in zip(datay, datax):
            data.append((zl, el))
        data.sort();
        z,eta = [], []
        for zl,el in data:
            z.append(zl); eta.append(el)

        zn, en =[], []
        n = len(z)
        if n:
            if z[0] > -self.cmb_km:
                zn.append(-self.cmb_km)
                en.append(eta[0])
            i=0
            while i < n:
                zn.append(z[i])
                en.append(eta[i])
                i += 1
            if n and z[n-1] < 0:
                zn.append(0)
                en.append(eta[n-1])
        datax = en
        datay = zn

        return datax,datay
       
    def get_layer_from_xy_data(self,datax,datay):

        """
        convert the point-based data to one that can be plotted as  layers
        """
        data = []
        i=0; n = len(datax);
        while i < n:
            if i > 0 and datay[i] != datay[i-1]:
               data.append((datay[i],datax[i-1]))
            data.append((datay[i],datax[i]))
            i += 1

        xl,yl = [],[]
        for yloc,xloc in data:
            xl.append(xloc)
            yl.append(yloc)

        return xl,yl
      

    def read_data(self,filename,mode):
        """
        
        read HC profile data from file, convert into our x-y format,
        and return datax, datay

        (to convert back, call convert_data with invert set to True)
        
        """
        datax,datay = [],[]
        f=open(filename,'r')
        for line in f:
            val = line.split()
            if len(val) != 2:
                print 'error file ', filename, ' appears to be in wrong format'
                print 'expecting'
                if self.mode == 1:
                    print 'radius[non_dim] viscosity[Pas]'
                elif self.mode == 2:
                    print 'radius[non_dim] density_scale'
                else:
                    print 'unknown'
                exit()
            datax.append(val[0])
            datay.append(val[1])
        f.close()
        # convert to plotting format
        return self.convert_data(mode,False,datax,datay)


    def convert_data(self,mode,reverse,datax,datay):
        """ convert input data to plotting format and vice versa """
        tmpx, tmpy = [],[]
        i=0
        for i in range(0,len(datax)):
            if mode == 1:       # viscosity
                if not reverse:
                    tmpx.append(float(datay[i])/self.visc_norm)
                    tmpy.append(-(1.-float(datax[i]))*self.radius_km)
                else:
                    tmpy.append(float(datax[i])*self.visc_norm)
                    tmpx.append((self.radius_km+float(datay[i]))/self.radius_km)
            elif mode == 2:     # density 
                if not reverse:
                    tmpx.append(float(datay[i]))
                    tmpy.append(-(1.-float(datax[i]))*self.radius_km)
                else:
                    tmpy.append(float(datax[i]))
                    tmpx.append((self.radius_km+float(datay[i]))/self.radius_km)
        return tmpx, tmpy

    def save_to_file(self):
        if self.verbose:
            print 'saving modified data'
        if self.plot_mode == 1:
        	self.outfile = self.tmpn + 'visc.dat'
        	print 'saving modified viscosity profile data to ', self.outfile
        elif self.plot_mode == 2:
        	self.outfile = self.tmpn + 'dens.dat'
        	print 'saving modified density profile data to ', self.outfile
            
        #
        # convert data back
        self.datax, self.datay = self.convert_data(self.plot_mode,True,self.datax,self.datay)
        f=open(self.outfile,'w')
        i=0
        for i in range(0,len(self.datax)):
            ostring = "%8.5f\t%12.7e\n" % (self.datax[i], self.datay[i])
            f.writelines(ostring)
        f.close()
        self.use = True;
        #p.close(self.figure)
        #self.figure.set_visible(False)
