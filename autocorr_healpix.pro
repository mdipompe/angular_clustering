;+
;  NAME:
;    autocorr_healpix
;
;  PURPOSE:
;    Use healpix at various resolutions to measure clustering instead of Landy & Szalay
;    (1993) estimator. 
;
;  USE:
;    autocorr_healpix,data,max_nside,bins,w_theta,doerrs=doerrs,outfile='out.txt',poly=manglepolys
;
;  INPUTS:
;    data - Structure with your data.  Should have tags ra and dec.
;    max_nside - maximium nside to calculate w_theta for (sets minimum scale)
;
;  KEYWORDS:
;    doerrs - if set, jackknife estimation of errors    
;
;  OPTIONAL INPUTS:
;    outfile - if set writes a text file with outputs (3 cols: bin
;              center, w_theta, error)
;    poly - only use pixels with centers within this region (defined
;           by mangle polygons)
;
;  OUTPUTS:
;    bins - bin centers
;    w_theta - correlation function
;
;  HISTORY:
;    2013 - Written - MAD (UWyo) 
;    5-27-15 - Cleaned up and generalized - MAD (UWyo)
;    7-20-15 - Added doerrs keyword, edits for clarity - MAD (UWyo)
;-
PRO autocorr_healpix,data,max_nside,bins,w_theta,outfile=outfile,poly=poly

;MAD Get start time
st=systime(1)

;MAD Convert RAs and DEC into phi and theta for HEALpix
phi=data.ra*(!dpi/180.)
theta=(90.-data.dec)*(!dpi/180.)

;MAD Get number in each sample in double format
n_data=double(n_elements(data))

;MAD Set full area of sphere in sq deg
fullarea=(4.*!dpi*(180./!dpi)^2.)

maxpow=round(ALOG10(max_nside)/ALOG10(2.))
w_theta=fltarr(maxpow+1)
w_theta_err=fltarr(maxpow+1)
bins=fltarr(maxpow+1)

;MAD Loop over each HEALpix level
FOR i=3L,maxpow+0.01 DO BEGIN
   print,'ANG_CLUSTER_PIX - Working on nside = ' + $
         strtrim(2.^i,2) + '...'
   Nside=2.^i
   Npix=12.*Nside^2.
   pixarea=fullarea/Npix
   pixsize=SQRT(pixarea)
   bins[i]=pixsize

   ang2pix_ring,nside,theta,phi,pix

   h=histogram(pix,bin=1,min=0,max=Npix-1)

   IF keyword_set(poly) THEN BEGIN
      healpix_coords,Nside,ra,dec,coord='Q'
      in=is_in_window(ra=ra,dec=dec,poly)
      usepix=where(in EQ 1)
   ENDIF

   delta=(h-mean(h[usepix]))/mean(h[usepix])
   delta=delta[usepix]

   w_theta[i]=total(delta^2.)/n_elements(usepix)

   IF keyword_set(doerrs) THEN BEGIN
      IF (n_elements(h) LT 50.) THEN BEGIN
         N=n_elements(h)
         flag=0
      ENDIF ELSE BEGIN
         N=50.
         flag=1
      ENDELSE
      w_theta_jack=fltarr(N)
      FOR j=0L,N-1. DO BEGIN
         IF (flag EQ 0.) THEN k=j ELSE k=round(randomu(systime_seed)*n_elements(h)-1.)
         new=h
         remove,k,new
         delta_new=(new-mean(new))/mean(new)
         w_theta_jack[j]=total(delta_new^2.)/(Npix-1.)
      ENDFOR 
      w_theta_err[i]=SQRT(((N-1.)/N)*total((mean(w_theta_jack)-w_theta_jack)^2.))
   ENDIF
ENDFOR

IF keyword_set(outfile) THEN BEGIN
   openw,1,outfile
   IF keyword_set(doerrs) THEN BEGIN
      FOR i=0L,n_elements(w_theta)-1 DO BEGIN
         printf,1,bins[i],w_theta[i],w_theta_err[i]
      ENDFOR
   ENDIF ELSE BEGIN
      FOR i=0L,n_elements(w_theta)-1 DO BEGIN
         printf,1,bins[i],w_theta[i]
      ENDFOR
   ENDELSE
   close,1
ENDIF

return
END
