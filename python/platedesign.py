"""
Simple and stupid "plate design" for the BOSSANOVA MilkyWay analog SDSS3 2014
proposal.
"""
import os.path

import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.basemap import Basemap
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
