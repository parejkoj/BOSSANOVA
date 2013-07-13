"""
Simple and stupid "plate design" for the BOSSANOVA MilkyWay analog SDSS3 2014
proposal.
"""
import os.path

import numpy as np
from scipy.interpolate import interp1d

import matplotlib
import matplotlib.pyplot as plt
from mpl_toolkits.basemap import Basemap
from matplotlib import cm
plt.ion()

import pyfits
import pairs.mr_spherematch

coll_r = 62./3600. # fiber collision radius

def one_plate(radec,mask=None):
    """
    Make a list of targets for one "plate", by checking fiber collision
    radii of everything in radec.
    mask is an index array of targets to include.
    """
    # make a list of collided objects
    sp = pairs.mr_spherematch.Spherematch(radec,radec.copy())
    match = sp(coll_r,maxmatch=len(radec))
    dictmatch = {}
    for x in match:
        # ignore self matches and objects that weren't already picked.
        if (x['idx1'] != x['idx2']):# and (x['idx1'] in mask):
            if (x['idx1'] not in dictmatch):
                dictmatch[x['idx1']] = []
            dictmatch[x['idx1']].append(x['idx2'])
    
    print len(dictmatch.keys()),'collision groups'
    # for each collision group, pick one object, and remove the others
    # from consideration
    intargs = np.ones(len(radec),dtype=bool)
    counter = 0
    # have to make a copy of the keys, in order to shuffle them.
    randtargs = list(dictmatch.keys())
    np.random.shuffle(randtargs)
    for targ in randtargs:
        collided = dictmatch[targ]
        if intargs[targ]:
            intargs[collided] = False
            counter += 1
    
    # select a random 900 elements to actually target:
    targs = np.nonzero(intargs)[0]
    np.random.shuffle(targs)
    return targs[:900]
#...

def read_block_par_file():
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
    return fibers,blocks


def plot_fiber_reach(blocks):
    """Plot the area reachable by each fiber."""
    fig = plt.figure()
    ax = fig.add_subplot(111)
    ax.axis('equal')
    reach = fiber_reach()
    colors = [cm.hsv(x) for x in np.linspace(0,1,len(blocks))]
    for b,fibers in blocks.items():
        for f in fibers:
            ax.plot(f[0]+reach[:,0],f[1]+reach[:,1],color=colors[b-1],lw=1)
        #if b > 10: break
    circ = matplotlib.patches.Circle((0,0),1.5,facecolor='none',edgecolor='black',lw=4)
    patch = ax.add_patch(circ)
    patch.set_zorder(4) # force the circle patch to draw on top
    ax.axis((-2.,2.,-2.,2.))
    ax.set_title('fiber reach, colored by harness')
    plt.savefig('../plots/fiber_reach.png',dpi=72,bbox_inches='tight')
            
def fiber_reach(interp='linear',npoints=1e2):
    """Returns an array of points representing the reach of a fiber in xy (degrees)."""
    reach = fiber_reach_values()
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

def fiber_reach_values():
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

def check_pluggability():
    """
    Needs fiberBlocksBOSS.par from PLATEDESIGN
    Taken from boss/bosstile/pro/chunk/plugtest.pro
    88	lambda= boss_lambdaeff(obj)
    89	ad2xyfocal, obj.ra, obj.dec, xf, yf, racen=tiles.racen,
    deccen=tiles.deccen, $
    90	  lambda=lambda
    
    101	blockfile=getenv('PLATEDESIGN_DIR')+'/data/boss/fiberBlocksBOSS.par'
    102	sdss_plugprob, xf, yf, fiberid, maxinblock=maxinblock,
    mininblock=mininblock, $
    103	  blockfile=blockfile, reachfunc='boss_reachcheck'
    104	
    105	;; check if all assigned
    106	ibad= where(fiberid le 0L, nbad)
    107	pluggable= nbad eq 0L
    """
#...

def plot_plate(radec,center,hexbin=False,filename=None):
    """Make a plot of the targets on their 'plate'."""
    fig = plt.figure()
    ax = fig.add_subplot(111)
    ax.axis('equal')
    height = 3
    width = 3
    #map = Basemap(projection='laea',lon_0=center[0],lat_0=center[1],
    #              celestial=True,height=height,width=width,rsphere=180./np.pi)
    circ = plt.Circle(center,radius=1.5,facecolor='none',edgecolor='red')
    ax.add_patch(circ)
    #x,y = map(radec[:,0],radec[:,1])
    x,y = radec[:,0],radec[:,1]
    if hexbin:
        ax.hexbin(x,y,mincnt=1)
    else:
        ax.scatter(x,y,s=6,marker='x')
    ax.text(0.84,0.93,'N points=%d'%len(radec),transform=ax.transAxes)
    # this somehow breaks the plot!
    #map.tissot(center[0],center[1],1.5,1000,edgecolor='red',facecolor='none')
    if filename:
        plt.savefig(os.path.join('../plots/',filename),bbox_inches='tight')
#...

def main(argv=None):
    data = pyfits.open('../data/ngc3254-r21.5_90arcmin.fits.gz')[1].data
    print "total targs:",len(data)
    # use this to not restrict targets at all.
    #test = np.ones(len(data),dtype=bool)
    test = data.r < 20.5
    radec = np.array(zip(data.ra,data.dec))[test]
    # the above is too much to actually use. We have to cut it down somehow.
    center = np.array(((157.33313,29.49180),))
    sp_cut = pairs.mr_spherematch.Spherematch(radec,center)
    match_cut = sp_cut(45./60,maxmatch=len(radec))
    data_cut = data[test][match_cut['idx1']]
    radec_cut = np.array(zip(data_cut.ra,data_cut.dec))
    print "reduced targs:",len(radec_cut)
    
    plot_plate(radec_cut,center[0],hexbin=True,filename='ngc3254-alltargs.pdf')
    
    unpicked = np.ones(len(radec_cut),dtype=bool)
    remaining = len(radec_cut[unpicked])
    while remaining > 0:
        targs = one_plate(radec_cut[unpicked])#,np.nonzero(unpicked)[0])
        plate = [a[targs] for a in np.where(unpicked)][0]
        unpicked[plate] = False
        remaining = len(radec_cut[unpicked])
    plot_plate(radec_cut[plate],center[0])
    print 'Remaining:',remaining
    import pdb
    pdb.set_trace()
    
    raw_input('blah')
#...

if __name__ == '__main__':
    sys.exit(main())
