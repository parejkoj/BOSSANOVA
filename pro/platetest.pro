; Read in a sample plate.
openr,1,"~/16559-sample.dat"
data=dblarr(2,900)
readf,1,data
close,1


nfibers = 900
nrepeat = 0
nskyperblock = 2
mininblock= (long(n_elements(nfibers)/50L)-nskyperblock+1L)>0L
maxinblock= (20L-1L)>mininblock


center = [154.887,58.2058]
ad2xyfocal,data[0,*],data[1,*],xf,yf,racen=center[0],deccen=center[1]

blockfile=getenv('PLATEDESIGN_DIR')+'/data/boss/fiberBlocksBOSS.par'
sdss_plugprob, xf, yf, fiberid, maxinblock=maxinblock, mininblock=mininblock, $
  blockfile=blockfile, reachfunc='boss_reachcheck'

