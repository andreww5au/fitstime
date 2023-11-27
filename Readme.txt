
FITS header time analysis

Written by Andrew Williams (at the time, at Perth Observatory, now Andrew.Williams@curtin.edu.au)

This code was written in the late 1990's to support on-site data reduction for the PLANET (Probing Lensing Anomolies
NETwork) microlensing group. We used a diverse set of optical telescopes around the world, dispersed in longitude
so that we could follow microlensing events continuously, 24 hours a day.

These telescopes had CCD cameras and software from the late 80s and early 90s, custom written. Some of the sites
(eg the European Southern Observatory) left observers without admin permissions on the computer, or even the ability
to install Python libraries into the system (and this was in Python 1.5 days, before local python environments
existed).

This ruled out the use of any Python library not in the standard library set, including FITS file IO and astrometry,
because they relied on C libraries which would need to be compiled (CFITSIO and EPHEM). So, I implemented my own FITS
library in pure Python (fits.py), and my own time conversion and astrometry code (coords.py).

The big problem was that each telescope and camera wrote FITS files, but with entirely different ways to specify
the observation date and time. The 'fitstime.py' and 'parseing.py' code reads in every date, time and coordinate field
used by any of the telescope sites, rated each value with a confidence level, and compared every date/time combination
with every other, including a few possible offsets (beginning vs midpoint of observation, 0.5 days for JD/MJD/etc
offsets, and heliocentric barycentre arrival time if an RA/DEC was specified).

The fitstime.py module itself, if run from the command line, will look for all of these headers:

  DATE     = geth(h,'DATE')
  DATEmOBS = geth(h,'DATE-OBS')
  DEC      = gethe(h,'DEC')
  DEC_OBJ  = gethf(h,'DEC_OBJ')
  DEC_OBS  = gethf(h,'DEC_OBS')
  EPOCH    = gethe(h,'EPOCH')
  EQUINOX  = gethe(h,'EQUINOX')
  EXPT     = gethf(h,'EXPT')
  EXPTIME  = gethf(h,'EXPTIME')
  EXPOSURE = gethf(h,'EXPOSURE')
  HJD      = gethf(h,'HJD')
  ITIME    = gethf(h,'ITIME')
  JD       = gethf(h,'JD')
  JDSTART  = gethf(h,'JDSTART')
  LJD      = gethf(h,'LJD')
  MJD      = gethf(h,'MJD')
  MJDmOBS  = gethf(h,'MJD-OBS')
  OBJECT   = geth(h,'OBJECT')
  OBSERVAT = geth(h,'OBSERVAT')
  OBSERVER = geth(h,'OBSERVER')
  RA       = gethe(h,'RA')
  RA_OBJ   = gethf(h,'RA_OBJ')
  RA_OBS   = gethf(h,'RA_OBS')
  TELESCOP = geth(h,'TELESCOP')
  TIME     = geth(h,'TIME')
  TIMEmOBS = geth(h,'TIME-OBS')
  TM_END   = gethf(h,'TM_END')
  TM_START = gethf(h,'TM_START')
  TMmSTART = gethf(h,'TM-START')
  UT       = geth(h,'UT')
  UTCmOBS  = geth(h,'UTC-OBS')
  UTDATE   = geth(h,'UTDATE')
  UTmDATE  = geth(h,'UT-DATE')
  UTmSTART = geth(h,'UT-START')
  UTmTIME  = geth(h,'UT-TIME')
  UTMIDDLE = geth(h,'UTMIDDLE')
  UTSHUT   = geth(h,'UTSHUT')

(and a few extra are created, for example, DATE-OBS, TM_START, etc, might just be a date, or just a time, or might
 be something like YYYY-MM-DD,HH:mm:SS, in which case the field is split into a date part and a time part).

For each date-like field, it attemps to determine the date order (YY-MM-DD, YYYY-MM-DD, DD-MM-YY or DD-MM-YYYY),
and how confident it is about that order (if what it thinks is the day value is > 12, for example, and what it thinks
is the year value is > 31, then it can be completely confident, while '01/02/03' could be in any order.

Then it combines all combinations of dates and times to create JD values, and compares those with all the JD-like
fields (JD, LJD, MJD, MJD-OBS, HJD), and checks to see if the differences between any two of them are equal to some
combination of:
   - Half the exposure duration
   - 0.5 days (for JD/MJD/etc offsets)
   - The heliocentric correction, for arrival time from the given RA/DEC to the solar system barycentre.

Once it has the best guess at the JD value for the HJD at the observation midpoint (HJD_Calc), it prints out either
just that value, or, if run with the -s option, the full analysis - eg:

------------------------------
ZOB03208R020.fits
Finding time in: ZOB03208R020.fits
No equinox or epoch of coordinates in FITS header, assuming 2000
Calculated HJD_Calc from (DATE-OBS & TM_START + Hel.Corr. + Exptime/2)
Hel.Corr = 0.005778 (499.26 sec)  Exptime/2 = 0.001736 (150.00 sec)
HJD_Calc = 2452812.855546943
         = MJD-OBS +Hel.Corr. +Exptime/2   (0.0 sec error)
MJD-OBS = 2452812.8480324
2452812.85555
------------------------------

Note that the analysis is likely to fail on modern FITS images - in particular, MJD is now usually JD-2450000.5, and
not 2400000.5, which will break the JD to date conversion.