"""
Simple and stupid "plate design" for the BOSSANOVA MilkyWay analog SDSS3 2014
proposal.
"""
import os.path
import sys

import numpy as np

import matplotlib
import matplotlib.pyplot as plt
from matplotlib import cm
plt.ion()

import pyfits
import pairs.mr_spherematch

coll_r = 62./3600. # fiber collision radius

def draw_circle(ax,center=(0,0),radius=1.5,lw=4,color='black'):
    """Draws an r degree radius circle at center on the current ax instance."""
    circ = matplotlib.patches.Circle(center,radius,facecolor='none',edgecolor=color,lw=lw)
    patch = ax.add_patch(circ)
    patch.set_zorder(4) # force the circle patch to draw on top

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
    targs = targs[:900]
    sp = pairs.mr_spherematch.Spherematch(radec[targs],radec[targs].copy())
    match = sp(coll_r,maxmatch=len(targs))
    print sum(np.nonzero(match['dist'] > 0)), 'still colliding'
    return targs
#...

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

def get_xy(radec,center):
    """Returns x,y to be plotted in ra/dec degrees."""
    return ((radec[:,0]-center[0])*np.cos(np.radians(radec[:,1])))+center[0],radec[:,1]

def plot_plate(radec,center,hexbin=False,filename=None):
    """Make a plot of the targets on their 'plate'."""
    fig = plt.figure()
    ax = fig.add_subplot(111)
    ax.axis('equal')
    draw_circle(ax,center=center,radius=1.5)
    draw_circle(ax,center=center,radius=1.,lw=2)
    ax.set_xlabel('RA (deg)')
    ax.set_ylabel('Dec (deg)')
    x,y = get_xy(radec,center)
    if hexbin:
        ax.hexbin(x,y,mincnt=1)
    else:
        ax.scatter(x,y,s=6,marker='x')
    ax.text(0.84,0.93,'N points=%d'%len(radec),transform=ax.transAxes)
    if filename:
        plt.savefig(os.path.join('../plots/',filename),bbox_inches='tight')
#...

def plot_one_plate(fig,ax,radec,center,color=None,label=''):
    """Make a plot of the targets on a single plate."""
    x,y = get_xy(radec,center)
    ax.scatter(x,y,s=6,marker='x',c=color,label=label)

def do_dat(targfile):
    """Generate a plate list for .dat files."""
    infile = file(targfile)
    header = infile.readline()
    name,ra,dec = header.split(' ')
    center = np.array((ra,dec),dtype=float)
    radec = np.loadtxt(infile)
    print 'targs: ',len(radec)
    plot_plate(radec,center,hexbin=True,filename='%s-alltargs.png'%name)
    plates = make_plates(radec,center,name)
    total = sum([len(p) for p in plates])
    remain = len(radec) - total
    complete = total/float(len(radec))*100
    
    fig = plt.figure(figsize=(8,8))
    ax = fig.add_subplot(111)
    ax.axis('equal')
    r = 1.5
    ax.axis((center[0]-r,center[0]+r,center[1]-r,center[1]+r))
    ax.set_xlabel('RA (deg)')
    ax.set_ylabel('Dec (deg)')
    colors = ['red','green','blue','cyan','brown','magenta','yellow']
    for i,(c,p) in enumerate(zip(colors,plates)):
        plot_one_plate(fig,ax,radec[p],center,color=c,label='%d: %d'%(i+1,len(p)))
    ax.set_title('NSA%s: %d targets remain (%d%% complete)'%(name,remain,complete))
    ax.legend()
    draw_circle(ax,center,1.5)
    draw_circle(ax,center,1.,lw=2)
    plt.savefig('../plots/%s-allplates.png'%name,dpi=100,bbox_inches='tight')
#...

def make_plates(radec,center,name):
    """Turn an radec list of targets into pluggable plates."""
    unpicked = np.ones(len(radec),dtype=bool)
    remaining = len(radec[unpicked])
    prev = len(radec)
    plates = []
    while remaining > 500:
        prev = remaining
        targs = one_plate(radec[unpicked])#,np.nonzero(unpicked)[0])
        plate = [a[targs] for a in np.where(unpicked)][0]
        unpicked[plate] = False
        remaining = len(radec[unpicked])
        plates.append(plate)
        #plot_plate(radec[plate],center)
        print 'Remaining:',remaining
        #raw_input('blah')
    return plates

def do_fits(targfile):
    data = pyfits.open(targfile)[1].data
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
    make_plates(radec_cut,center[0])
#...

def main(argv=None):
    if argv is None: argv = sys.argv[1:]
    from optparse import OptionParser, OptionGroup
    usage = '%prog [OPTS] TARGETFILE'
    parser = OptionParser(usage)

    i = 0
    while True:
        try:
            targfile = argv[i]
        except IndexError:
            if i == 0:
                print 'Must provide TARGETFILE.'
                return -1
            else:
                return
        if '.fits' in targfile:
            do_fits(targfile)
        elif '.dat' in targfile:
            do_dat(targfile)
        i+=1
#...

if __name__ == '__main__':
    sys.exit(main())
