
"""PLANET Event name and date handling routines (Julian day, PJD, etc)
   Andrew Williams
"""

import math
import string
import time

def dsin(x):
  return math.sin(float(x)/180.0*math.pi)

def dcos(x):
  return math.cos(float(x)/180.0*math.pi)

def dtan(x):
  return dsin(x)/dcos(x)

def dasin(x):
  return math.asin(x)*180/math.pi

def datan2(y,x):
  return math.atan2(y,x)*180/math.pi

J2000=2451544.5


sites=['Canopus', 'Perth', 'SAAO', 'ESO', 'Danish', 'CTIO']


def sexstring(value=0,sp=':'):
  """Convert the floating point 'value' into a sexagecimal string.
     The character in 'sp' is used as a spacer between components. Useful for
     within functions, not on its own.
     eg: sexstring(status.TJ.ObjRA,' ')
  """
  try:
    aval=abs(value)
    error=0
  except:
    aval=0.0
    error=1
  if value<0:
    outs='-'
  else:
    outs=''
  D=int(aval)
  M=int((aval-float(D))*60)
  S=float(int((aval-float(D)-float(M)/60)*36000))/10
  outs=outs+str(D)+sp+str(M)+sp+str(S)
  if error:
    return ''
  else:
    return outs


def juldate(data=None):
  "Return full Julian Day for a given time tuple. Use current date/time if no arg given"
  if data:
    year,month,day,hour,minute,second,wd,dnum,dst = data
  else:
    year,month,day,hour,minute,second,wd,dnum,dst = time.gmtime(time.time())

  if (month == 1) or (month == 2):
    year = year - 1
    month = month + 12

  A = math.floor(year/100.0);
  B = 2 - A + math.floor(A/4.0);
  jd = math.floor(365.25 * year) + math.floor(30.6001 * (month + 1))
  jd = jd + day + (hour + (minute/60.0) + (second/3600.0)) / 24.0
  jd = jd + 1720994 + B + 0.5;
  return jd
  

def pjd(data=None):
  return juldate(data) - 2450000.0


def caldate(JD=0):
  "Return tuple (year,month,day) for full Julian Day. Use current date/time if no arg given"
  if not JD:
    JD=juldate()
  Z = int(JD + 0.5)
  F = (JD + 0.5) - int(JD + 0.5)
  if Z<2299161:
    A = Z
  else:
    alpha = int( (Z - 1867216.25)/36524.25 )
    A = Z + 1 + alpha - int(alpha/4)
  B = A + 1524
  C = int( (B - 122.1)/365.25 )
  D = int( 365.25*C )
  E = int( (B - D)/30.6001 )
  day = B - D - int(30.6001*E) + F
  if E<13.5:
    month = E - 1
  else:
    month = E - 13
  if month>2.5:
    year = C - 4716
  else:
    year = C - 4715
  return (year,month,day)



def formatpjd(lpjd=None):
  """Returns a formatted string for a given PJD. Use the current date if arg=None
     if the argument can not be converted to a float, return it unchanged.
  """
  if lpjd==None:
    lpjd=pjd()
  try:
    lpjd=round(float(lpjd), 1)
  except:
    return `lpjd`
  year,month,day=caldate(lpjd+2450000)
  month=month-1  #List elements numbered from 0 not 1
  months=['Jan','Feb','Mar','Apr','May','Jun','Jul','Aug','Sep','Oct',
          'Nov','Dec']
  return "%6.1f (%s %3.1f)" % (lpjd, months[month], day)


def cname(n='HB-2K-060'):
  return n[:2]+n[3:5]+n[6:]


def lname(n='HB-2K-060'):
  t=n[0]
  if t=='M':
    team='MACHO'
  elif t=='O':
    team='OGLE'
  elif t=='E':
    team='EROS'
  elif t=='K':
    team='MOA'
  else:
    team=t

  s=n[1]
  if s=='B':
    source='BLG'
  elif s=='L':
    source='LMC'
  elif s=='S':
    source='SMC'
  else:
    source=s

  y=n[3:5]
  if y[0]=='9':
    year='19'+y
  elif y[0]=='0':
    year='20'+y
  elif y=='2K':
    year='2000'
  else:
    year=y
  return team+'-'+year+'-'+source+'-'+n[6:]



def _pnparse(obj):
  """
     If obj is a PLANET name, then code,yr,num,filt and imgnum are returned as
     a tuple.  PLANET names have:

  -the letter O, M, K, or E.
  -the letter B
  -two digit year, equal to '98', '99', '2K', '01', '02', or '03'
  -a two or three digit event number (a leading 0 is added for two digit cases)
  -an optional filter letter, assumed to be I if missing. If an image 
     number is present, the filter is _not_ optional.
  -an optional image number, returned as '' if missing.

  If the argument is not a valid PLANET name, then None is returned. Check
  that the result is true before attempting to resolve the tuple.

  """
  pl=0    #Start by assuming the argument is _not_ valid
  obj=string.upper(obj)

  if len(obj)>0:             #Second letter must be O, M, or E
    if obj[0]=='O' or obj[0]=='M' or obj[0]=='E' or obj[0]=='K':
      pl=1                   #Seems to be valid so far...

  if len(obj) < 6:
    pl=0                #If it's too short, it can't be valid
  else:
    if obj[1]<>'B':
      pl=0              #if the scond letter isn't B it's not valid
    yr=obj[2:4]         #Test the 3rd and 4th chars for a valid year
    if (yr<>'98') and (yr<>'99') and (yr<>'2K') and (yr<>'01') and (yr<>'02') and (yr<>'03'):
      pl=0 

    tail=obj[4:]       #Take the rest of the name
    num=''
    im=''
    filt=''
    gf=0

    #Assign all digits before the first non-digit to 'num', and all
    #digits after the first non-digit to 'im'. The non-digit/s are
    #stored in 'filt'.
    for a in tail:
      if a in string.digits and not gf:
        num=num+a
      else:
        gf=1
        if a in string.digits:
          im=im+a
        else:
          filt=filt+a

    if len(filt)<>0 and len(filt)<>1:
      pl=0       #If filter field is present, but >1 letter, it's invalid
      filt=''

    if (tail<>num+filt+im) or (len(num) > 3) or (len(num) < 2) or (len(im)>3):
      pl=0   #Invalid if any extra junk, or if num/im are too long/short

  if pl==1:
    code=obj[:2]    #eg 'OB'
    num=string.zfill(num,3)     #Add leading zeroes with zfill
    if im<>'':                  #If there's an image number, pad it out
      im=string.zfill(im,3)
    return (code,yr,num,filt,im)
  else:
    return None


def zstrip(s=''):
  "Strip leading zeroes from the given string"
  l=list(s)
  while l[0]=='0':
    l=l[1:]
  return string.join(l,'')


def jd_year(jd=None):
  "Turn a JD into a decimal year"
  y,m,d = caldate(jd)
  if y == -1:
    y = -2
  e0 = juldate(data=(y,1,1,0,0,0,0,0,0))
  e1 = juldate(data=(y+1,1,1,0,0,0,0,0,0))
  return y + (jd-e0)/(e1-e0)


def JDtoT(jd=None):
  "Return T value, decimal century, from JD"
  return (float(jd)-2415020.0)/36525

                                                                                
def range(v=0, r=360.0):
  "ensures that v is in the range 0 <= v < r"
  return v - r*math.floor(v/r);


def precess(jd1=None, jd2=None, ra=None, dec=None):
  "Precess coords (in degrees) from epoch jd1 to jd2"
  alpha_in = ra
  delta_in = dec
  
  #precession progresses about 1 arc second in .047 years */
  #From from_equinox to 2000.0 */
                                                                                  
  if (abs(jd1 - J2000)/365.25) > .04:
    T = JDtoT(jd1) - 1
    zeta_A  = 0.6406161* T + 0.0000839* T*T + 0.0000050* T*T*T
    z_A     = 0.6406161* T + 0.0003041* T*T + 0.0000051* T*T*T
    theta_A = 0.5567530* T - 0.0001185* T*T - 0.0000116* T*T*T
                                                                                
    A = dsin(alpha_in - z_A) * dcos(delta_in)
    B = ( dcos(alpha_in - z_A) * dcos(theta_A) * dcos(delta_in)
                               + dsin(theta_A) * dsin(delta_in) )
    C = (-dcos(alpha_in - z_A) * dsin(theta_A) * dcos(delta_in)
                              + dcos(theta_A) * dsin(delta_in) )
                                                                                
    alpha2000 = datan2(A,B) - zeta_A
    alpha2000 = range (alpha2000, 360.0)
    delta2000 = dasin(C)
  else:
    alpha2000 = alpha_in
    delta2000 = delta_in

                                                                            
  #From 2000.0 to to_equinox */
                                                                                
  if (abs(jd2 - J2000)/365.25) > .04:
    T = JDtoT(jd2) - 1
    zeta_A  = 0.6406161* T + 0.0000839* T*T + 0.0000050* T*T*T
    z_A     = 0.6406161* T + 0.0003041* T*T + 0.0000051* T*T*T
    theta_A = 0.5567530* T - 0.0001185* T*T - 0.0000116* T*T*T
                                                                                
    A = dsin(alpha2000 + zeta_A) * dcos(delta2000)
    B = ( dcos(alpha2000 + zeta_A) * dcos(theta_A) * dcos(delta2000)
         - dsin(theta_A) * dsin(delta2000) )
    C = ( dcos(alpha2000 + zeta_A) * dsin(theta_A) * dcos(delta2000)
         + dcos(theta_A) * dsin(delta2000) )
                                                                                   
    alpha = datan2(A,B) + z_A
    alpha = range(alpha, 360.0)
    delta = dasin(C)
  else:
    alpha = alpha2000
    delta = delta2000
                                                                                
  return alpha, delta


def hjd(jd=None, ra=None, dec=None):
  """ Calculate heliocentric correction - input angles in degrees
      double e;            /* obliquity of ecliptic */
      double n;            /* day number */
      double g;            /* solar mean anomaly */
      double L;            /* solar ecliptic longitude */
      double l;            /* mean solar ecliptic longitude */
      double R;            /* sun distance, AU */
      double X, Y;         /* equatorial rectangular solar coords */
      double cdec, sdec;
      double cra, sra;
      double deltajd;      /* HJD = JD - deltajd */
  """
                                                                                
  #precess from J2000 (input equinox) to the observation equinox,
  #where times are given in MJD
                                                                                
  ra, dec = precess(J2000, jd, ra, dec)
                                                                                
  # do it to it
                                           
  cdec = dcos(dec)
  sdec = dsin(dec)
  cra = dcos(ra)
  sra = dsin(ra)
                                                                                
  n = jd - 2451545.0                    #use epoch 2000
  e = 23.439 - 0.0000004*n
  g = 357.528 + 0.9856003*n
  L = 280.461 + 0.9856474*n
  l = L + 1.915*dsin(g) + 0.02*dsin(2.0*g)
  R = 1.00014 - 0.01671*dcos(g) - 0.00014*dcos(2.0*g)
  X = R*dcos(l)
  Y = R*dcos(e)*dsin(l)
   
  deltajd = 0.0057755 * (cdec*cra*X + (cdec*sra + dtan(e)*sdec)*Y)
 
  return (jd-deltajd)



