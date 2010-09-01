
"""FITS header field parsing functions

   Written by Andrew Williams, Perth Observatory
   <andrew@physics.uwa.edu.au>
"""



version = "$Revision$"

import coords

import re
import string

resf=r"[-+]?(?:\d+(?:\.\d*)?|\d*\.\d+)(?:[eE][-+]?\d+)?"
reuf=r"(?:\d+(?:\.\d*)?|\d*\.\d+)"
mlen=[31,29,31,30,31,30,31,31,30,31,30,31]
months=['JAN','FEB','MAR','APR','MAY','JUN','JUL','AUG','SEP','OCT','NOV','DEC']

dateorder=None        #Set to 'YMD' or 'DMY' to override guessing in getdate


def ptuple(t,sp=' '):
  "Prints a tuple of integers with the given seperator string"
  return string.join(map(str,map(int,t)),sp)


def getdate(s=""):
  """Attempt to parse the given string as a date. The format is unknown,
     but it's assumed that the broken American m/d/y isn't a possibility.
     It attempts to distinguish between d/m/y and y/m/d using the values,
     and handles any seperator. 

     if dateorder (global) is 'DMY' or 'YMD', it overrides the last-resort
     guessing used if the order is ambiguous, but otherwise has no effect. 
     For example, "2003-11-12" is still parsed correctly if dateorder is 'DMY'.

     returns y,m,d,confidence where confidence is 1 if the date order is
     definitely correct, and 0 if it's ambiguous. Throws an AssertionError
     exception if the input is definitely not a date triple.
  """
  global parseoutput, YearGuess
  nums=re.findall(reuf,s)
  nums=map(float,nums)     #Convert from strings to floats
  if len(nums)==2:
    month=0
    #May be broken Canopus UTdate with month name
    s=string.upper(s)
    for i in range(12):
      if string.find(s, months[i])>=0:
        month=i+1
    assert month, "Only two numbers, and month name not found in '"+s+"'"
    nums=[nums[0],month,nums[1]]
    monthname=1
  else:
    monthname=0

  assert len(nums)==3, "Too many/few numbers in '"+s+"'"

  #Now either nums[0] is the year and nums[2] is the day, or vice-versa
  #either way, nums[1] is the month

  month=int(nums[1])
  assert (month>=1) and (month<=12), "Month invalid in '"+s+"'"

  if nums[0]>100:
    #If first number is >100, it must be a 4-digit-year, so the third must be a day
    assert (nums[2]>=1) and (nums[2]<=mlen[month-1]) and (nums[0]>1980) and (nums[0]<2050), "YMD, Day invalid in '"+s+"'"
    day=nums[2]
    year=nums[0]
    if YearGuess:
      if year <> YearGuess:
        parseoutput += "Best guess at year is "+`YearGuess`+", clashes with "+`year`+" from "+s
    else:
      YearGuess = int(year)
    return (year,month,day),1       #confident it's YMD

  if nums[2]>100:
    #If third number is >100, it must be a 4-digit-year, so the first must be a day
    assert (nums[0]>=1) and (nums[0]<=mlen[month-1]) and (nums[2]>1980) and (nums[2]<2050), "DMY, Day invalid in '"+s+"'"
    day=nums[0]
    year=nums[2]
    if YearGuess:
      if year <> YearGuess:
        parseoutput += "Best guess at year is "+`YearGuess`+", clashes with "+`year`+" from "+s
    else:
      YearGuess = int(year)
    return (year,month,day),1       #confident it's DMY

  #OK, at this point all three numbers are less than or equal to 100

  if nums[0]>50:
    #If first number is >50, it must be a pre-2000 2-digit-year, so the third must be a day
    assert (nums[2] >= 1) and (nums[2] <= mlen[month-1]), "YMD, Day invalid in '"+s+"'"
    day=nums[2]
    if nums[0]<100:
      year=nums[0]+1900
    else:
      year=nums[0]
    if YearGuess:
      if year <> YearGuess:
        parseoutput += "Best guess at year is "+`YearGuess`+", clashes with "+`year`+" from "+s
    else:
      YearGuess = int(year)
    return (year,month,day),1       #confident it's YMD

  if nums[2]>50:
    #If third number is >50, it must be a pre-2000 2-digit-year, so the first must be a day
    assert (nums[0] >= 1) and (nums[0] <= mlen[month-1]), "DMY, Day invalid in '"+s+"'"
    day=nums[0]
    if nums[2]<100:
      year=nums[2]+1900
    else:
      year=nums[2]
    if YearGuess:
      if year <> YearGuess:
        parseoutput += "Best guess at year is "+`YearGuess`+", clashes with "+`year`+" from "+s
    else:
      YearGuess = int(year)
    return (year,month,day),1       #confident it's DMY

  #At this point, all numbers are <=50, so could conceivably be in either order. Try yearguess first

  if YearGuess:
    yg = int(str(int(YearGuess))[-2:])
    if (nums[0] == YearGuess) or (nums[0] == yg):
      assert (nums[2] >= 1) and (nums[2] <= mlen[month-1]), "YMD, Day invalid in '"+s+"'"
      year = YearGuess
      day = nums[2]   
      return (year,month,day),1       #confident it's YMD
    elif (nums[2] == YearGuess) or (nums[2] == yg):
      assert (nums[0] >= 1) and (nums[0] <= mlen[month-1]), "DMY, Day invalid in '"+s+"'"
      year = YearGuess
      day = nums[0]      
      return (year,month,day),1       #confident it's DMY

  #First and last number are both valid days, so we hope the LARGER one is the day

  if ( (nums[2]>nums[0]) and (not monthname) ) or (dateorder=='YMD'):
    assert (nums[2] >= 1) and (nums[2] <= mlen[month-1]), "guess YMD, Day invalid in '"+s+"'"
    day=nums[2]
    year=nums[0]+2000
    if dateorder == 'YMD':
      return (year,month,day),0.5       #Relatively sure since we have a specified date order
    else:
      print "Warning - guessing at YMD order for '"+s+"'"
      parseoutput += "Warning - guessing at YMD order for '"+s+"'\n"
      return (year,month,day),0       #Only guess it's YMD
  else:
    assert (nums[0] >= 1) and (nums[0] <= mlen[month-1]), "guess DMY, Day invalid in '"+s+"'"
    day=nums[0]
    year=nums[2]+2000
    if dateorder == 'DMY':
      return (year,month,day),0.5       #Relatively sure since we have a specified date order
    else:
      print "Warning - guessing at DMY order for '"+s+"'"
      parseoutput += "Warning - guessing at DMY order for '"+s+"'\n"
      return (year,month,day),0       #Only guess it's DMY


def gettimestring(s="", angle=None):
  """Attempt to parse the given string as a time. The format is three 
     numbers, with any seperator/s, H,M,S order.

     returns h,m,s as floats
     Throws an AssertionError exception if the input is definitely not a time triple.

     if angle=="dms", allow angles in degrees from -90 to +90
     if angle=="udms" allow angles in degrees from 0 to 360
     if angle=="hms" allow angles in hours from 0 to 24

     for 'dms' case, if angle is negative, _all_ components returned are negative,
     to handle cases -01:00:00 < v < 00:00:00
  """
  assert len(s)>5,"Too short to be a time in '"+s+"'"

  if angle=="dms" and s[0]=='-':
    sign=-1
  else:
    sign=1
  nums=re.findall(reuf,s)
  nums=map(float,nums)     #Convert from strings to floats

  assert len(nums)==3, "Too many/few numbers in '"+s+"'"
  hour,minute,second = tuple(nums)

  if angle=="dms":
    assert (hour>=0) and (hour<=90), "Angle outside DMS range in '"+s+"'"
  elif angle=="udms":
    assert (hour>=0) and (hour<360), "Angle outside UDMS range in '"+s+"'"
  else:
    assert (hour>=0) and (hour<24), "Hour outside HMS range in '"+s+"'"
  assert (minute>=0) and (minute<=59), "Minute outside range in '"+s+"'"
  assert (second>=0) and (second<=60), "Second outside range in '"+s+"'"  #Leap seconds possible
  
  return (hour*sign, minute*sign, second*sign)


def gettime(v=None):
  """Accepts either a number or a string. If it's a number, interprets that as a time either 
     as an integer number of seconds since midnight, or decimal hours. If a string, calls 
     gettimestring to parse it as an H,M,S triple. Either way, returns h,m,s,confidence
     (where confidence is 0 if it's impossible to distinguish between seconds and hours)
     or None,0 if all attempts fail.
  """
  if not v:
    return None,0
  if type(v)==type(""):
    try:
      h,m,s = gettimestring(v)
      return (h,m,s),1           #Confident it's a time tuple
    except AssertionError:
      return None,0            #Its a string, but not a time tuple
  else:
    if v>86400:            #Probably a unix time stamp
      print "Warning - UNIX timestamp in '"+`v`+"' not parsed."
      parseoutput += "Warning - UNIX timestamp in '"+`v`+"' not parsed\n"
      return None,0
    elif v>24.0:           #definitely time in seconds since midnight
      tmp=coords.sexstring(v/3600.0, ' ')     #Convert from seconds into decimal hours
      h,m,s = tuple(map(float, string.split(tmp)))
      return (h,m,s),1          #Confident it's seconds since midnight
    elif int(v) <> v:           #<24 and it has a fractional part, almost certainly decimal hours
      tmp=coords.sexstring(v, ' ') 
      h,m,s = tuple(map(float, string.split(tmp)))
      return (h,m,s),1         
    else:                   #<24 but an integer, probably decimal hours but not sure
      tmp=coords.sexstring(v, ' ') 
      h,m,s = tuple(map(float, string.split(tmp)))
      return (h,m,s),0          #Hard to be sure, it could be <24 seconds after midnight UT



def getnumber(s=None, signed=1):
  """Accepts a string. If there is ONE number anywhere in that string, return it as a float,
     otherwise return None. If 'signed' is true (default), accept negative numbers. Note that 
     this returns None if there are multiple distinct numbers in the string.
  """
  if signed:
    nums=re.findall(resf,s)
  else:
    nums=re.findall(reuf,s)
  nums=map(float,nums)     #Convert from strings to floats
  if len(nums)==1:
    return nums[0]
  else:
    return None


def geth(h,name):
  "Grabs a FITS string field and returns the stripped string minus the quotes"
  try:
    v = string.strip(h[name])
    if len(v)>=2:
      if v[0]=="'" or v[0]=='"':
        v=string.strip(v[1:-1])
    return v
  except KeyError:
    return None


def gethf(h,name):
  """Assumes the field is a numeric value (not a string) and if so returns a float. If
     it doesn't exist, or if it's not a valid number, return None
  """
  try:
    s=h[name]
  except KeyError:
    return None
  try:
    return float(s)
  except ValueError:
    return None


def gethe(h,name):
  """Firstly try interpreting the field as a float, and return that if it works, or None
     if the key doesn't exist. If it's not parseable as a float, return the stripped
     string, minus the quotes.
  """
  v=gethf(h,name)
  if v is None:
    v=geth(h,name)
  return v


def parseheader(h=None,comments=None, yearguess=None):
  """
     parseheader returns lists of all values in each category (all dates, all ras, etc). Each list
     is composed of tuples, being (value, field, confidence), where value is the number, field is 
     a string containing the name of the FITS header field it was derived from, and confidence
     is a number from 0 to ~200 containing the 'confidence' that that field is valid'. All field
     values are converted to the same units:

     dates: (y,m,d) tuples
     times: (h,m,s) tuples
     jds:   full julian days, including fractional part
     hjds:  full julian days, with heliocentric correction
     ras:   fractional hours
     decs:  fractional degrees
     equinoxes:  fractional years
     exptimes:   seconds

  """

  global parseoutput, YearGuess
  parseoutput = ''
  YearGuess = yearguess

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


  DATEmOBSd = None        #d and t suffixes refer to components, eg '2003-06-24T06:39:12.152'
  DATEmOBSt = None
  if DATEmOBS:
    tmp=string.split(DATEmOBS, 'T')
    if len(tmp)==2:
      DATEmOBSd = tmp[0]
      DATEmOBSt = tmp[1]

  DATEd = None
  DATEt = None
  if DATE:
    tmp=string.split(DATE, 'T')
    if len(tmp)==2:
      DATEd = tmp[0]
      DATEt = tmp[1]
  
  dates=[]       #A list of valid ((y,m,d),"header fields",confidence) tuples
  times=[]       #A list of valid ((h,m,s),"header fields",confidence) tuples
  jds=[]         #A list of valid (jd,"header fields",confidence) tuples (non-HJD)
  hjds=[]        #A list of valid (hjd,"header fields",confidence) tuples
  ras=[]         #A list of valid (ra,"header fields",confidence) tuples (RA in hours)
  decs=[]        #A list of valid (dec,"header fields",confidence) tuples (DEC in degrees)
  equinoxes=[]   #A list of valid (equinox,"header fields",confidence) tuples (eg 1950, 2000, etc)
  exptimes=[]    #A list of valid (exptime,"header fields",confidence) tuples, exptime in seconds

#Dates, and time fields merged with date values

  if DATEmOBSd:
    date,c = getdate(DATEmOBSd)
    dates.append((date,"DATE-OBSd",c*100))
  elif DATEmOBS:
    date,c = getdate(DATEmOBS)
    dates.append((date,"DATE-OBS",c*100))
    if (date[2]-int(date[2]))>1e-9:       #Fractional day number
      tmp=coords.sexstring((date[2]-int(date[2])*24), ' ')    #decimal days -> hours
      h,m,s = tuple(map(float, string.split(tmp)))
      times.append( ((h,m,s),"DATE-OBS",c*100) )


#Note that DATE and TIME fields are often only the file creation date/time, not image acquisition
  if UTDATE:
    date,c = getdate(UTDATE)
    dates.append((date,"UTDATE",c*120))
    if (date[2]-int(date[2]))>1e-9:       #Fractional day number
      tmp=coords.sexstring((date[2]-int(date[2])*24), ' ')    #decimal days -> hours
      h,m,s = tuple(map(float, string.split(tmp)))
      times.append( ((h,m,s),"UTDATE",c*100) )

  if UTmDATE:
    date,c = getdate(UTmDATE)
    dates.append((date,"UT-DATE",c*120))
    if (date[2]-int(date[2]))>1e-9:       #Fractional day number
      tmp=coords.sexstring((date[2]-int(date[2])*24), ' ')    #decimal days -> hours
      h,m,s = tuple(map(float, string.split(tmp)))
      times.append( ((h,m,s),"UT-DATE",c*100) )

  if EPOCH:
    epy=None
    if (type(EPOCH)==type("")):
      try:
        date,c = getdate(EPOCH)
        dates.append((date,"EPOCH",c*50))
        if (date[2]-int(date[2]))>1e-9:       #Fractional day number
          tmp=coords.sexstring((date[2]-int(date[2])*24), ' ')    #decimal days -> hours
          h,m,s = tuple(map(float, string.split(tmp)))
          times.append( ((h,m,s),"EPOCH",c*100) )
      except AssertionError:
        #If EPOCH isn't a valid date, it's probably a decimal year for RA/DEC coordinates
        epy=getnumber(EPOCH)
    else:
      epy=EPOCH
    if epy is not None:
      if (abs(epy-1950.0)<1e-5) or (abs(epy-2000)<1e-5):  #If it's 1950 or 2000 exactly
        equinoxes.append((epy,"EPOCH",100))
      else:         #It's probably a julian day number with an unknown offset, ignore it
        pass

  if DATEd:
    date,c = getdate(DATEd)
    dates.append((date,"DATEd",c*20))
  elif DATE:
    date,c = getdate(DATE)
    dates.append((date,"DATE",c*20))
    if (date[2]-int(date[2]))>1e-9:       #Fractional day number
      tmp=coords.sexstring((date[2]-int(date[2])*24), ' ')    #decimal days -> hours
      h,m,s = tuple(map(float, string.split(tmp)))
      times.append( ((h,m,s),"DATE",c*20) )

#Times in the header

  if TIMEmOBS:
    tm,c=gettime(TIMEmOBS)
    if tm:
      times.append((tm,"TIME-OBS",c*100))

  if DATEmOBSt:
    tm,c=gettime(DATEmOBSt)
    if tm:
      times.append((tm,"DATE-OBSt",c*100))

  if UTCmOBS:
    tm,c=gettime(UTCmOBS)
    if tm:
      times.append((tm,"UTC-OBS",c*100))

  if TM_START:
    tm,c=gettime(TM_START)
    if tm:
      times.append((tm,"TM_START",c*150))

  if TMmSTART:
    tm,c=gettime(TMmSTART)
    if tm:
      times.append((tm,"TM-START",c*150))

  if UTmSTART:
    tm,c=gettime(UTmSTART)
    if tm:
      times.append((tm,"UT-START",c*150))

  if UTmTIME:
    tm,c=gettime(UTmTIME)
    if tm:
      times.append((tm,"UT-TIME",c*150))

  if UTSHUT:
    tm,c=gettime(UTSHUT)
    if tm:
      times.append((tm,"UTSHUT",c*50))

#Note that DATE and TIME fields are often only the file creation date/time, not image acquisition
  if TIME:
    tm,c=gettime(TIME)
    if tm:
      times.append((tm,"TIME",c*20))

  if DATEt:
    tm,c=gettime(DATEt)
    if tm:
      times.append((tm,"DATEt",c*20))

  if UT:
    tm,c=gettime(UT)
    if tm:
      times.append((tm,"UT",c*50))

#JDs in the header, from JD, MJD, MJD-OBS fields. Non-heliocentric, image start times

  if JD:
    jds.append((JD,"JD",100))

  if MJDmOBS:
    jds.append((MJDmOBS+2400000.5, "MJD-OBS", 90))

  if MJD:
    jds.append((MJD+2400000.5, "MJD", 80))

  if JDSTART:
    jds.append((JDSTART, "JDSTART", 100))

#Heliocentric JD values in the header

  if HJD:
    hjds.append((HJD, "HJD", 100))

#RA's in the header

  if RA:
    if (type(RA)==type("")):
      ra=gettimestring(RA, angle="hms")
      if ra:
        ra=ra[0]+ra[1]/60.0+ra[2]/3600.0
        ras.append((ra,"RAstr",100))
    else:
      try:
        racomm=comments['RA']
      except KeyError:
        racomm=''
      if (RA>=24.0) or string.find(racomm,"deg")>=0:
        ras.append((RA/15.0,"RAdeg",100))
      else:
        ras.append((RA,"RAhour",50))
      
  if RA_OBJ:
    if (type(RA_OBJ)==type("")):
      ra=gettimestring(RA_OBJ, angle="hms")
      if ra:
        ra=ra[0]+ra[1]/60.0+ra[2]/3600.0
        ras.append((ra,"RA_OBJstr",100))
    else:
      try:
        racomm=comments['RA_OBJ']
      except KeyError:
        racomm=''
      if (RA_OBJ>=24.0) or string.find(racomm,"deg")>=0:
        ras.append((RA_OBJ/15.0,"RA_OBJdeg",100))
      else:
        ras.append((RA_OBJ,"RA_OBJhour",50))
      
  if RA_OBS:
    if (type(RA_OBS)==type("")):
      ra=gettimestring(RA_OBS, angle="hms")
      if ra:
        ra=ra[0]+ra[1]/60.0+ra[2]/3600.0
        ras.append((ra,"RA_OBSstr",100))
    else:
      try:
        racomm=comments['RA_OBJ']
      except KeyError:
        racomm=''
      if (RA_OBS>=24.0) or string.find(racomm,"deg")>=0:
        ras.append((RA_OBS/15.0,"RA_OBSdeg",100))
      else:
        ras.append((RA_OBS,"RA_OBShour",50))
      
#DECs in the header

  if DEC:
    if (type(DEC)==type("")):
      dec=gettimestring(DEC, angle="dms")
      if dec:
        dec=dec[0]+dec[1]/60.0+dec[2]/3600.0
        decs.append((dec,"DECstr",100))
    else:
      decs.append((DEC,"DECdeg",100))
   
  if DEC_OBJ:
    if (type(DEC_OBJ)==type("")):
      dec=gettimestring(DEC_OBJ, angle="dms")
      if dec:
        dec=dec[0]+dec[1]/60.0+dec[2]/3600.0
        decs.append((dec,"DEC_OBJstr",100))
    else:
      decs.append((DEC_OBJ,"DEC_OBJdeg",100))
   
  if DEC_OBS:
    if (type(DEC_OBS)==type("")):
      dec=gettimestring(DEC_OBS, angle="dms")
      if dec:
        dec=dec[0]+dec[1]/60.0+dec[2]/3600.0
        decs.append((dec,"DEC_OBSstr",100))
    else:
      decs.append((DEC_OBS,"DEC_OBSdeg",100))

#Equinoxes in the header. For EPOCH field, if  equal to 1950 or 2000, high confidence, if
# 1995<eq<2005, assume its a fractional year with medium confidence, otherwise assume
#the epoch is 2000 with low confidence. EPOCH handled above with dates

  if EQUINOX:
    if type(EQUINOX)==type(""):
      eq=getnumber(EQUINOX)
      if eq:
        equinoxes.append((eq,"EQUINOX",100))
    else:
      equinoxes.append((EQUINOX,"EQUINOX",100))

  if EXPTIME:
    if type(EXPTIME)==type(""):
      et=getnumber(EXPTIME)
      if et:
        exptimes.append((et,"EXPTIME",100))
    else:
      exptimes.append((EXPTIME,"EXPTIME",100))

  if EXPOSURE:
    if type(EXPOSURE)==type(""):
      et=getnumber(EXPOSURE)
      if et:
        exptimes.append((et,"EXPOSURE",100))
    else:
      exptimes.append((EXPOSURE,"EXPOSURE",100))

  if ITIME:
    if type(ITIME)==type(""):
      et=getnumber(ITIME)
      if et:
        exptimes.append((et,"ITIME",90))
    else:
      exptimes.append((ITIME,"ITIME",90))


  return dates,times,jds,hjds,ras,decs,equinoxes,exptimes,parseoutput


