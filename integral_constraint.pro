;+
;  NAME:
;    integral_constraint
;  PURPOSE:
;    Calculate the integral constraint correction for an angular
;    autocorrelation (see Equation 8 of Roche & Eales, 1999, MNRAS,
;    307, 703).  Need to fit a power law and multiply by amplitude
;    to get correction factor to add to w(theta)
;
;  USE:
;    C=integral_constraint(rr=rr,theta=theta,slope=slope)
;
;  INPUT:
;    
;  Optional Inputs:
;    rr - Random-Random counts in theta bins.  Can be array or file
;         name.  Defaults to RR.txt (written by ang_cluster.pro)
;    theta - angular scales of rr values (degrees).  If not set, will
;            generate them using define_bins.pro and 4 bins per dex
;            out to 2 degrees.
;    slope - the assumed slope of the angular autocorrelation.
;            Defaults to -1.
;
;  OUTPUT:
;    C - the correction factor (subtract this times the amplitude from
;        the measured autocorrelation at each scale)
;
;  HISTORY:
;    1-4-16 - Written - MAD (Dartmouth)
;-
FUNCTION integral_constraint,rr=rr,theta=theta,slope=slope

  ;MAD Set some defaults, if not supplied
  IF ~keyword_set(rr) THEN rr='RR.txt'
  IF ~keyword_set(slope) THEN slope=-1.
  IF ~keyword_set(theta) THEN BEGIN
     bininfo=define_bins(4.,2.)
     theta=(bininfo.cents/60.)
  ENDIF

  ;MAD Read in RR counts from file, or use supplied array
  IF (size(rr,/type) EQ 7) THEN readcol,rr,rrs,format='D' ELSE rrs=rr
  
  ;MAD Calculate correction
  num=total(rrs*(theta^slope))
  den=total(rrs)
  C=num/den

  return, C
END

