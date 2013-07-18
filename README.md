BOSSANOVA
=========

For work on the BOSSANOVA proposal for SDSS3 spring 2014 time. Put general
useful notes in this readme, to be incorporated into other things later.

Redshift completness of BOSS
====
This is described in Dawson et al. 2013, section 5.1. Taking SN^2 >20 for
red and >10 for blue gives 85% redshift completeness for CMASS targets with
21.5 < ifiber2 < 21.7, and 63% for fainter CMASS targets.

http://adsabs.harvard.edu/abs/2013AJ....145...10D

We would definitely want to use the version of idlspec2d that uses galaxy
priors, to reduce catastrophic redshift failures, though we may have to
modify the galaxy templates. See Bolton et al. 2012, section 3.2 and 4.1,
and also Figures 7 and 8 for a S/N vs. redshift failure rate plot. Takeaway
point: SN_r^2 ~> 9 should get us what we want.

http://adsabs.harvard.edu/abs/2012AJ....144..144B

Sky fibers, etc.
====
Need room for sky fibers (20), calibration stars (>=80), guide stars (16). 
We need to ensure that there is room for these fibers to distribute in
the same region as the object fibers. At least one sky fiber per
20-fiber bundle. Dawson et al. Section 3.2.

APOGEE bulge plates plug 300 fibers quite tightly into the middle.

Tightly packing the fibers makes it more likey that a fiber will fall
out. If we're willing to have a larger fraction of fibers fall out,
it will be fine. Re-plugging fallen fibers is hard if they are very
tight.

Bright Stars
====
Drill holes around stars brighter than m_V=7.5 of 52.5" diameter to
prevent scattered light, chosen from Tycho-2 catalog.
Dawson et al. Section 3.2

Plate design/plugging
====
The relevant code seems to be spread between the bosstile and
platedesign products.

  bosstile/pro/chunk/bosstile_plugtest.pro
runs through the tiles and calls:
  bosstile/pro/chunk/plugtest.pro
which converts positions in sky-coordinates to focal plane positions given a
plate center, reads in the fiber harness block descriptions and then calls:
  platedesign/pro/plate/sdss_plugprob.pro
which actually assigns fibers to targets. 
sdss_plugprob.pro uses (via the parameter reachfunc):
  platedesign/pro/plate/boss_reachcheck.pro
to determine how far each fiber can reach, which takes the values from:
  platedesign/pro/plate/boss_reachvalues.pro
which should probably actually be in a .par file, but whatever.
