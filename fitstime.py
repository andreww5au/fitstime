#!/usr/bin/python

version = "$Revision$"

import sys

import parseing
import fits
import coords

basefield = 'HJD_Calc'       #The base julian day field to use for output times
                             #The default, HJD_Calc, is derived from the best date
                             #and time fields, plus half the exptime, plus the 
                             #heliocentric offset.
#These can be changed to either +1 or -1 by command line arguments
hcorr=0   #Don't add or subtract the heliocentric correction to the base field 
ecorr=0   #Don't add or subtract half the exptime to the base field result
mcorr=0   #Any extra modifier, generally +/- 0.5 for broken JD/MJD conversions



usage="""FITS header time analysis - Andrew Williams
usage:  fitstime [-h|-help|--help]  OR
        fitstime [options] [filename] [filename] ...

where options can be:
-s, -S, --show   Give a full analysis of all time fields in the header,
                 and how they relate. Note that if this is specified, all
                 of the time-offset parameters are ignored, and the result
                 will be an analysis of the times IN the FITS header, not
                 corrected and output for use.

=[fieldname]     The base header field to use for output time, subject to the
                 offsets described below. The fieldname can be HJD_Calc, HJD,
                 JD, MJD, or MJD-OBS, and all have previously been converted
                 to full Julian day numbers.
                 The default, HJD_Calc, is not from the FITS header directly,
                 but calculated from the best available date and time fields
                 in the header, plus the heliocentric offset, and half the
                 exposure time.

[+|-][hel|exp]   A plus or minus, then 'hel' or 'exp' (or similar). Add or
                 subtract either the heliocentric correction, or half the
                 exposure time, to get the final output time for each image.
                 Any word containin 'hel' or 'exp' will do,
                 eg "+Hel.Corr. -Exptime/2"

[+|-][number]    An arbitrary offset to add to get the final output time for
                 each image. Must have a leading plus or minus sign. Typically
                 +0.5 or -0.5 to correct for unknown original MJD and MJD-OBS
                 offsets in the FITS header.

-[dmy|ymd]       If two-digit years, combined with a year after '00', make the
                 order ambiguous, the default bahaviour is to issue a warning,
		 and guess. If '-dmy' or '-ymd' is given on the command line, 
		 the specified order is used as a last resort, instead of 
		 warning or guessing. The (broken) 'month/day/year' order isn't
		 supported.

Typically you would do 'fitstime --show *.fits' to examine the time fields
present in your images, and check them for consistency. Normally, you could
then 'fitstime *.fits > times.txt' to produce a sorted HJD-midpoint time
list. 

If you see from the analysis output that a different field (eg JD) produces 
midpoint-HJD times, or you wish to apply an offset, then specify this on the 
command line. For example, if there is no date/time to get HJD_Calc, you might
want to specify 'fitstime =MJD-OBS +hel +exp -0.5 *.fits > times.txt'.
"""

def cmpconf(a,b):
  """compare a and b based on the confidence level of the tuple, used in sorting"""
  return -cmp(a[2],b[2])


def getra(ras, verbose=1):
  """Extract the best RA from the values in the FITS header. Guess at 18h if none provided"""
  outstring = ''
  if ras:
    ras.sort(cmpconf)
    raval=ras[0][0]
    rafield=ras[0][1]
    if verbose and (len(ras)>1):
      if (raval-ras[1][0]) > 1e-2:         #We only care about RAs for HJD calculations
        outstring += rafield+" and "+ras[1][1]+" disagree, using "+rafield+"="+coords.sexstring(raval,':') + '\n'
    return raval,rafield,outstring
  else:
    if verbose:
      outstring += "No RA field, assume 18h for HJD_Calc\n"
    return 18.0, "guess", outstring


def getdec(decs, verbose=1):
  """Extract the best DEC from the values in the FITS header. Guess at -28d if none provided"""
  outstring = ''
  if decs:
    decs.sort(cmpconf)
    decval=decs[0][0]
    decfield=decs[0][1]
    if verbose and (len(decs)>1):
      if (decval-decs[1][0]) > 1e-2:         #We only care about DECs for HJD calculations
        outstring += decfield+" and "+decs[1][1]+" disagree, using "+decfield+"="+coords.sexstring(decval,':') + '\n'
    return decval,decfield,outstring
  else:
    if verbose:
      outstring += "No DEC field, assume -28d for HJD_Calc\n"
    return -28.0,"guess",outstring


def getdate(dates, verbose=1):
  """Extract the best date from the values in the FITS header. Return None if none provided"""
  outstring = ''
  if dates:
    dates.sort(cmpconf)
    dateval=dates[0][0]
    datefield=dates[0][1]
    datecon=dates[0][2]
    if verbose and (len(dates)>1):
      jd1 = coords.juldate(data=(dateval[0], dateval[1], int(dateval[2]), 0,0,0, 0,0,0))
      d2 = dates[1][0]
      jd2 = coords.juldate(data=(d2[0], d2[1], int(d2[2]), 0,0,0, 0,0,0))
      if ((jd1-jd2) > 1e-2) and (dates[1][2]>0):     #If dates differ by >1day, and confidence>0    
        outstring += datefield+" and "+dates[1][1]+" disagree, using "+datefield+"="+parseing.ptuple(dateval,':') +'\n'
    return dateval,datefield,outstring
  else:
    if verbose:
      outstring += "No detected date field in header\n"
    return None,"None",outstring


def gettime(times, dates, verbose=1):
  """Extract the best UT time from the values in the FITS header. Return None if none provided"""
  outstring = ''
  if times:
    times.sort(cmpconf)
    timeval=times[0][0]
    timefield=times[0][1]
    timecon=times[0][2]
    if verbose and (len(times)>1):
      t1 = timeval[0]*3600 + timeval[1]*60 + timeval[2]     #in seconds
      t2 = times[1][0][0]*3600 + times[1][0][1]*60 + times[1][0][2]  #in seconds
      if ((t1-t2) > 10) and (times[1][2]>0):       #If times differ by >10sec, and confidence>0
        outstring += timefield+" and "+times[1][1]+" disagree, using "+timefield+"="+parseing.ptuple(timeval,':') + '\n'
    return timeval,timefield,outstring
  else:
    if verbose:
      outstring += "No detected time field in header\n"
    return None,"None",outstring


def getjd(jds, verbose=1):
  """Extract the best JD from the values in the FITS header. Return None if none provided"""
  if jds:
    jds.sort(cmpconf)
    jdval=jds[0][0]
    jdfield=jds[0][1]
    return jdval,jdfield,''
  else:
    return None,"None",'No JD fields found\n'


def getexptime(exptimes, verbose=1):
  """Extract the best JD from the values in the FITS header. Return None if none provided"""
  if exptimes:
    exptimes.sort(cmpconf)
    exptimeval=exptimes[0][0]
    exptimefield=exptimes[0][1]
    return exptimeval,exptimefield,''
  else:
    return None,"None",'No exposure time field found\n'


def getequinox(equinoxes, verbose=1):
  """Extract the best equinox of coordinates from the values in the FITS header. Return 2000 if
     none provided
  """
  outstring = ''
  if equinoxes:
    equinoxes.sort(cmpconf)
    equinoxval = equinoxes[0][0]
    equinoxfield = equinoxes[0][1]
    return equinoxval,equinoxfield,outstring
  else:
    if verbose:
      outstring += "No equinox or epoch of coordinates in FITS header, assuming 2000\n"
    return 2000.0,"guess",outstring


def yearfromheaders(hdic=None):
  """Parse other header fields (JD's, object name, whatever) to try and guess the year, so
     that broken 2-digit years in dates can be distinguished from day numbers.
  """
  jdl = []
  for field in ['JD', 'MJD', 'MJD-OBS', 'HJD', 'LJD']:
    try:
      val = float(hdic[field])
      if val < 2000000:
        val = val + 2400000
      jdl.append(val)
    except:
      continue
  if len(jdl) == 0:
    return None

  jdl.sort()
  jd = jdl[(len(jdl)-1) // 2]     #Take the middle value

  y,m,d = coords.caldate(jd)
  return y
    
  

class HeaderFields:
  pass          #An instance of this is used to store the header fields sorted by group.


def findtime(fname='', fimage=None, verbose=1, allfields=0):
  try:
    if not fimage:
      f = fits.FITS(fname,'h')
    else:
      f = fimage

    yearguess = yearfromheaders(f.headers)    #Parse other header fields for year to break 2-digit-year degeneracy
    hf = HeaderFields()
    (hf.dates,hf.times,hf.jds,hf.hjds,
     hf.ras,hf.decs,hf.equinoxes,hf.exptimes, outstring) = parseing.parseheader(f.headers, f.comments, yearguess)

    #parseheader returns lists of all values in each category (all dates, all ras, etc). Each list
    #is composed of tuples, being (value, field, confidence), where value is the number, field is 
    #a string containing the name of the FITS header field it was derived from, and confidence
    #is a number from 0 to ~200 containing the 'confidence' that that field is valid'. 

  except AssertionError:
    print "#Error opening or parseing FITS headers in file: "+fname
    sys.excepthook(*sys.exc_info())
    print
    if allfields:
      return None,"#Error opening or parseing FITS headers in file: "+fname+"\n", None
    else:
      return None,"#Error opening or parseing FITS headers in file: "+fname+"\n"

  if verbose and fname:
    outstring += "\nFinding time in: "+fname+"\n"

  ftime,ftimefield,os2 = gettime(hf.times, hf.dates, verbose=verbose)
  fdate,fdatefield,os1 = getdate(hf.dates, verbose=verbose)
  fjd,fjdfield,os3 = getjd(hf.jds+hf.hjds, verbose=verbose)    
  fra,frafield,os4 = getra(hf.ras, verbose=verbose)
  fdec,fdecfield,os5 = getdec(hf.decs, verbose=verbose)
  fequinox,fequinoxfield,os6 = getequinox(hf.equinoxes, verbose=verbose)
  fexptime,fexptimefield,os7 = getexptime(hf.exptimes, verbose=verbose)
  outstring += "".join([os1,os2,os3,os4,os5,os6,os7])

  #Above calls extract the 'best' value in each category, by sorting based on confidence. If there is
  #more than one value for a category, compare the best and second best to check consistency. If there's
  #no value in a category, make up a reasonable default if appropriate (eg equinox of J2000).

  if fexptime:
    edelta = (fexptime/2.0)/86400       #Half the exposure time, in days
  else:
    if verbose:
      outstring += "No exposure time information, can't verify or calculate time offsets\n"
    edelta = 0.0


  #Calculate HJD_Calc from the best date and time field, and find the heliocentric offset

  if fdate and ftime:
    cjd = coords.juldate(data=(fdate[0], fdate[1], int(fdate[2]),
                                 ftime[0], ftime[1], ftime[2], 0,0,0))
    chjd = coords.hjd(jd=cjd, ra=fra*15.0, dec=fdec)
    hdelta = chjd - cjd
    chjd = chjd + edelta
    hf.hjds.insert(0,(chjd,"HJD_Calc",200))
    chjdfield = "(" + fdatefield + " & " + ftimefield + " + Hel.Corr. + Exptime/2)"
  elif fjd:
    ghjd = coords.hjd(jd=fjd, ra=fra*15.0, dec=fdec)
    hdelta = ghjd - fjd
  else:
    if verbose:
      outstring += "No date and time, or any form of JD field. No idea how to find the time for this image...\n"
    chjdfield = "(no data)"
    if allfields:
      return None, outstring, hf
    else:
      return None, outstring


  if (abs(hdelta<1e-4) or (edelta<1e-4)) and verbose:
    outstring += "Heliocentric or exptime/2 offset is less than 0.0001 days, can't reliably compare JD offsets\n"

  if verbose:
    outstring += "Calculated HJD_Calc from "+chjdfield+"\n"
    outstring += "Hel.Corr = %8.6f (%6.2f sec)  Exptime/2 = %8.6f (%6.2f sec) \n" % (hdelta,
                                                                                     hdelta*86400,
                                                                                     edelta,
                                                                                     edelta*86400)

  hdl = [(-hdelta," -Hel.Corr."), (0.0,""), (+hdelta," +Hel.Corr.")]
  edl = [(-edelta," -Exptime/2"), (0.0,""), (+edelta," +Exptime/2")]
  mdl = [(-0.5," -0.5"), (0.0,""), (+0.5," +0.5")]

  jdict = {}
  for item in hf.jds + hf.hjds:
    jdict[item[1]] = item[0]

  if not verbose:
    try:
      out = jdict[basefield]
    except KeyError:         #The base field isn't available
      return None, outstring + "Invalid base field specified\n"
    out = out + hcorr*hdelta + ecorr*edelta + mcorr
    if allfields:
      return out, outstring, hf
    else:
      return out, outstring

  jlist = jdict.keys()
  jlist.sort()

  #Now compare every field containing a JD or HJD value with every other, checking for close matches,
  #for all possible combinations of heliocentric, half-exptime and half-day offsets. For each match, 
  #store the match key names, a string specifying which offsets were used, and and error, in days.

  matches = {}
  clashes = {}
  for akey in jlist:
    matches[akey] = {}
    clashes[akey] = []
    for bkey in jlist[jlist.index(akey)+1:]:
      ajd = jdict[akey]
      bjd = jdict[bkey]
      if abs(ajd-bjd)>1:     #Two JD values differ by more than a day
        clashes[akey].append(bkey)
      for hd in hdl:
        for ed in edl:
          for md in mdl:
            offset = hd[0] + ed[0] + md[0]
            offsetstring = hd[1] + ed[1] + md[1]
            if abs(ajd - (bjd+offset)) < 2e-4:       #About 17 seconds
              matches[akey][bkey] = ( offsetstring, ajd - (bjd+offset) )

  for akey in jlist:
    outstring += akey + " = " + `jdict[akey]` + '\n'
    for bkey in matches[akey].keys():
      outstring += " "*len(akey) + " = " + bkey + (matches[akey][bkey][0] + 
                   "   (" + str(round(matches[akey][bkey][1]*86400,2)) +" sec error)\n" )

  for akey in jlist:
    for bkey in clashes[akey]:
      outstring += "Warning: %s=%9.5f and %s=%9.5s differ by more than one day!" % (akey, jdict[akey], bkey, jdict[bkey])
      print "Warning: %s=%9.5f and %s=%9.5s differ by more than one day!" % (akey, jdict[akey], bkey, jdict[bkey])

  try:
    out = jdict[basefield]
  except KeyError:         #The base field isn't available
    return None, outstring + "Invalid base field specified\n"
  out = out + hcorr*hdelta + ecorr*edelta + mcorr
  if allfields:
    return out, outstring, hf
  else:
    return out, outstring



####################################################################

#Main program

if __name__ == '__main__':
  args=sys.argv[1:]
  if not args:
    print usage
    sys.exit()

  verbose=0     #Don't verbosely analyse the file, just print "filename time"
  parseing.dateorder=None      #Don't override best guess at date order - can also be 'DMY' or 'YMD'

  signs={'-':-1, '+':+1}

  files=[]
  for ar in args:
    if ar=='-s' or ar=='-S' or ar=='-show' or ar=='--show':
      verbose=1
    elif ar=='-h' or ar=='-help' or ar=='--help':
      print usage
      sys.exit()
    elif ar=='-DMY' or ar=='-dmy' or ar=='--DMY' or ar=='--dmy':
      parseing.dateorder='DMY'
    elif ar=='-YMD' or ar=='-ymd' or ar=='--YMD' or ar=='--ymd':
      parseing.dateorder='YMD'
    elif ar[0]=='=':
      ac=ar[1:]
      if not ac:
        sys.exit("Invalid option '=', must specify a field name after the '='")
      ac = ac.upper()
      if ac=='HJD_CALC':
        basefield='HJD_Calc'
      elif ac=='HJD' or ac=='JD' or ac=='MJD' or ac=='MJD-OBS':
        basefield=ac
      else:
        sys.exit("Invalid base field '" + ac + "'")
    elif ar[0]=='-' or ar[0]=='+':
      scorr=signs[ar[0]]
      ac=ar[1:]
      if not ac:
        break             #If we're passed a lone + or - arg, treat rest of command line as files
      ac = ac.upper()
      if ac.find('HEL') > -1:
        hcorr=scorr
      elif ac.find('EXP') > -1:
        ecorr=scorr
      else:
        try:
          mcorr=float(ar)
        except:
          mcorr=0            #Whatever it is after the +/- it's not a number or a valid modifier
          sys.exit("Unknown uption '" + ar + "'")
    else:
      files.append(ar)

  for f in files:
    if verbose:
      print '\n',f,
    else:
      print f,
    t,comments=findtime(fname=f,verbose=verbose)
    if t:
      print comments,t
    else:
      print "***No Data***"
  
