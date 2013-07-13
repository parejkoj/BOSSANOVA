"""
Defines properties of a BOSS plate, including harness positions and fiber reach.

example:
  import bossplate
  bossplate = bossplate.BossPlate()
  # make some plots:
  bossplate.plot_fiber_reach(1.0)
  bossplate.plot_fiber_reach(0.75)
"""
import numpy as np
from scipy.interpolate import interp1d

import matplotlib
import matplotlib.pyplot as plt
from mpl_toolkits.basemap import Basemap
from matplotlib import cm
plt.ion()

def inside_circle(x,y,radius):
    """Boolean array of the points that are inside circle centered at 0."""
    return np.sqrt(x**2+y**2) <= radius

class BossPlate(object):
    def __init__(self):
        self.read_block_par_file()
        self.nblocks = len(self.blocks)
        self.fibers_per_block = len(self.blocks[1])
        self.reach = self.fiber_reach()
        self.colors = [cm.hsv(x) for x in np.linspace(0,1,len(self.blocks))]
    
    def _setup_plot(self):
        """Sets up a single equal-axis plot with a plate circle (r=1.5deg)."""
        fig = plt.figure()
        ax = fig.add_subplot(111)
        ax.axis('equal')
        ax.axis((-2.,2.,-2.,2.))
        self.fig = fig
        self.ax = ax
        self._draw_circle()
    
    def _draw_circle(self,radius=1.5,lw=4):
        """Draws an r degree radius circle at 0,0 on the current ax instance."""
        circ = matplotlib.patches.Circle((0,0),radius,facecolor='none',edgecolor='black',lw=lw)
        patch = self.ax.add_patch(circ)
        patch.set_zorder(4) # force the circle patch to draw on top

    def plot_fiber_reach(self,radius=1.5,save=False):
        """Plot the area reachable by each fiber."""
        self._setup_plot()
        for b,fibers in self.blocks.items():
            for f in fibers:
                x = f[0] + self.reach[:,0]
                y = f[1] + self.reach[:,1]
                inside = inside_circle(x,y,radius)
                self.ax.plot(x[inside],y[inside],color=self.colors[b-1],lw=1)
        if radius != 1.5:
            self._draw_circle(radius,2)
        self.ax.set_title('fiber reach, r=%4.2f$^\circ$'%radius)
        if save:
            plt.savefig('../plots/fiber_reach_%4.2f.png'%radius,dpi=72,bbox_inches='tight')
    #...
    
    def plot_fiber_noreach(self,radius=1.5,save=False):
        """
        Plots the fibers that do not reach into the given radius.
        
        Rather stupid, as if r is small, the harness inside r will register as
        not reaching because none of their reach positions are within r.
        """
        self._setup_plot()
        self._draw_circle(radius,2)
        reaches = np.zeros(len(self.fibers),dtype=bool)
        for b,fibers in self.blocks.items():
            for i,f in enumerate(fibers):
                x = f[0] + self.reach[:,0]
                y = f[1] + self.reach[:,1]
                inside = inside_circle(x,y,radius)
                fiberid = (b-1)*self.fibers_per_block + i
                reaches[fiberid] = np.any(inside)
                if not reaches[fiberid]:
                    self.ax.scatter(f[0],f[1],facecolor=self.colors[b-1],edgecolor='none')
        missing = sum(~reaches)
        frac_missing = missing/float(len(self.fibers))
        self.ax.set_title("%d (%3f%%) fibers didn't reach"%(missing,frac_missing))
        if save:
            plt.savefig('../plots/fiber_noreach_%4.2f.png'%radius,dpi=72,bbox_inches='tight')        
    
    def read_block_par_file(self):
        """
        Reads fiberBlocksBOSS.par and returns a numpy array.
        Should really use the .parfile stuff from yanny, but...
        """
        import csv
        types = ['i4','f8','f8']
        names = ['blockid','x','y']
        infile = file('../data/fiberBlocksBOSS.par')
        incsv = csv.reader(infile,delimiter=' ')
        i = 0
        # remove the header, etc.
        while i < 18:
            line = infile.readline()
            i += 1
        fibers = []
        for line in incsv:
            fibers.append((line[1],line[2],line[3]))
        fibers = np.array(fibers,dtype=np.dtype(zip(names,types)))
        blocks = {}
        for f in fibers:
            block = f['blockid']
            if block  not in blocks:
                blocks[block] = np.zeros(20,dtype=zip(names[1:],types[1:]))
                i = 0
            blocks[block][i] = (f['x'],f['y'])
            i += 1
        self.fibers = fibers
        self.blocks = blocks
    
    def fiber_reach(self,interp='linear',npoints=1e2):
        """Returns an array of points representing the reach of a fiber in xy (degrees)."""
        reach = self._fiber_reach_values()
        r = np.sqrt((reach[:]**2).sum(axis=1))
        theta = np.arctan2(reach[:,1],reach[:,0])
        # interp1d needs monotonically increasing values
        idx = np.argsort(theta)
        theta = theta[idx]
        r = r[idx]
        new_theta = np.arange(-1,1,1/npoints)*np.pi
        # 3rd order spline interpolator
        f = interp1d(theta,r,kind=interp,bounds_error=False)
        new_r = f(new_theta)
        return np.array(zip(new_r*np.cos(new_theta),new_r*np.sin(new_theta)))
    
    def _fiber_reach_values(self):
        """
        Returns an array of the fiber reach values for a boss cart.
        In order to determine the actual reach of a given fiber,
        one should interpolate between these values in polar coordinates.
        
        Values are given in degrees relative to each individual fiber position.
        
        Taken from platedesign/pro/plate/boss_reachvalues.pro
        """
        platescale = 217.7358
        xcm= [25., 20., 15., 10., 5., 0., -5., -10., -15., -17., -17., 
            -17., -15., -10., -5., 0., 5., 10., 15., 20.]
        ycm= [0., 15., 19., 21., 22., 22., 22., 20., 15., 2., 0.,
            -2., -15., -20., -22., -22., -22., -21., -19., -15.]
        return 10*np.array(zip(xcm,ycm))/platescale
#...
