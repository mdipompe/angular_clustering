;+
;  NAME:
;    ang_corr_pix
;  PURPOSE:
;    Use pixelization scheme to calculate angular auto/cross
;    correlation (Scranton et al. 2002 equations 7 & 8)
;
;  USE:
;    ang_corr_pix,delta1,clong,clat,binedges,wth,delta2=delta2
;
;  INPUT:
;    delta1 - delta (value relative to mean) in each subregion.
;    clong - longitude of subregion centers (degrees)
;    clat - latitude of subregion centers (degrees)
;    binedges - array with edges of bins defined
;
;  OPTIONAL INPUT:
;    delta2 - delta (value relative to mean) in each subregion for
;             second parameter.  If set, does a cross-correlation with
;             delta1, if not an autocorrelation of delta1 is done.
;
;  KEYWORDS:        
;
;  OUTPUT:
;    wth - correlation function
;
;  NOTES:
;
;  HISTORY:
;    6-1-15 - Written - MAD (UWyo)
;-
PRO ang_corr_pix,delta1,clong,clat,binedges,wth,delta2=delta2

IF ~keyword_set(delta2) THEN delta2=delta1

FOR i=0L,n_elements(binedges)-2 DO BEGIN
   counter,i,n_elements(binedges)
   spherematch,clong,clat,clong,clat,binedges[i+1],m1,m2,sep,maxmatch=0
   xx=where(sep GT binedges[i])
   tmp=total((delta1[m1[xx]]*delta2[m2[xx]]))*(1./n_elements(m1[xx]))
   IF (n_elements(wth) EQ 0) THEN wth=tmp ELSE wth=[wth,tmp]
ENDFOR

return
END
