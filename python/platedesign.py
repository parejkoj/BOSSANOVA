"""
Simple and stupid "plate design" for the BOSSANOVA MilkyWay analog SDSS3 2014
proposal.
"""

import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.basemap import Basemap
plt.ion()

import csv_parse
import pairs.mr_spherematch

coll_r = 62./3600. # fiber collision radius

def one_plate(radec,mask):
    """Make a list of targets for one "plate", by checking fiber collision
    radii of everything in radec.
    mask is an index array of targets to include."""
    # make a list of collided objects
    sp = pairs.mr_spherematch.Spherematch(radec,radec.copy())
    match = sp(coll_r,maxmatch=len(radec))
    dictmatch = {}
    for x in match:
        if (x['idx1'] != x['idx2']) and (x['idx1'] in mask):
            if (x['idx1'] not in dictmatch):
                dictmatch[x['idx1']] = []
            dictmatch[x['idx1']].append(x['idx2'])
    
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
    return np.random.choice(np.nonzero(intargs)[0],900)
#...

def plot_plate(radec,center,hexbin=False):
    """Make a plot of the targets on their "plate"."""
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
    # this somehow breaks the plot!
    #map.tissot(center[0],center[1],1.5,1000,edgecolor='red',facecolor='none')
#...

head,data = csv_parse.read('ngc3254-r21.5_90arcmin.csv')
print "total targs:",len(data)
#test = np.ones(len(data),dtype=bool)
test = data.r < 21.1
radec = np.array(zip(data[test].ra,data[test].dec))
# the above is too much to actually use. We have to cut it down somehow.
center = np.array(((157.33313,29.49180),))
sp_cut = pairs.mr_spherematch.Spherematch(radec,center)
match_cut = sp_cut(45./60,maxmatch=len(radec))
data_cut = data[match_cut['idx1']]
radec_cut = np.array(zip(data_cut.ra,data_cut.dec))
print "reduced targs:",len(radec_cut)

plot_plate(radec_cut,center[0],hexbin=True)

unpicked = np.ones(len(radec_cut),dtype=bool)
while len(radec_cut[unpicked]) > 0:
    print 'Remaining:',len(radec_cut[unpicked])
    targs = one_plate(radec_cut,np.nonzero(unpicked)[0])
    unpicked[targs] = False
    plot_plate(radec_cut[targs],center[0])
    import pdb
    pdb.set_trace()

raw_input('blah')
