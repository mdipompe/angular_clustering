;+
;  NAME:
;    model_autocorr
;  PURPOSE:
;    Generate a model projected angular autocorrelation of matter
;
;  USE:
;    model_autocorr,theta,mod_w,power_spec='power_spec.fits',dndz='dndz.txt',$
;                      zarray=zarray,$
;                      omega_m=omega_m,omega_l=omega_l,h0=h0,omega_b=omega_b,$
;                      outfile='model_autocorr.txt',paramfile='params.ini'
;
;  INPUT:
;    theta - Array of angular scales you want the model at (degrees)
;
;  Optional Inputs:
;    power_spec - structure or string name of fits file 
;                 containing the matter power
;                 spectrum.  Needs tags pk (the power spectrum), k
;                 (wavenumber), z (redshift).
;                 Can be made from CAMB data using
;                 matter_power_spec.pro/camb4idl.pro and
;                 combine_camb.pro. If not set, will call necessary
;                 procedures to make it.
;    dndz - either a string name of file with your real z values,
;           which will be fit with a spline function and written out
;           to dndz_fit.txt, or an array of dndz values that
;           correspond to zarray values (can generate with
;           fit_dndz.pro). If not set, defaults to look for file
;           called dndz.txt with real z values.
;    zarray - array of z values power spectrum and dndz are calculated
;             for.  If not sets, defaults to 0.01 through 4.0 in steps
;             of 0.01
;    omega_m - Omega_matter, defaults to 0.273.  Be sure it is
;              consistent with what you used when you calculated the
;              power spectrum if you supply one!
;    omega_l - Omega_lambda, defaults to 0.727.  See above...
;    h0 - little h (H0/100), defaults to 0.702.  See above...
;    omega_b - omega_baryon, defaults to 0.046
;    outfile - if supplied, writes model power out to text file
;    plotfile - if supplied, makes a plot of the model
;    bz - coefficients of b(z) model, of the form b(z) = bz[0] + bz[1](1+z)^2
;
;  KEYWORDS:
;    
;  OUTPUT:
;    mod_w - the model autocorrelation at each input theta
;
;  NOTES:
;    Make sure your cosmology agrees with how you calculated the power
;    spectrum, if you supply your own. Also check the dndz_fit.txt file to make sure nothing
;    went wrong if you don't supply your own fit!
;
;  HISTORY:
;    8-12-15 - Written - MAD (UWyo)
;    8-21-15 - If power spec not supplied, calls CAMB4IDL to get it -
;              MAD (UWyo)
;    11-6-15 - Added b(z) functionality - MAD (Dartmouth)
;   11-19-15 - Added option to supply own dndz rather than auto
;              calling fit_dndz.txt (more flexible) - MAD (Dartmouth)
;-
PRO model_autocorr,theta,mod_w,power_spec=power_spec,dndz=dndz,$
                   zarray=zarray,paramfile=paramfile,$
                   omega_m=omega_m,omega_l=omega_l,h0=h0,omega_b=omega_b,$
                   outfile=outfile,plotfile=plotfile,$
                   bz=bz

;MAD If output file already exists, don't run just read it in
IF ~keyword_set(outfile) THEN check='' ELSE check=file_search(outfile)
IF (check NE '') THEN BEGIN
   print,'MODEL_AUTOCORR: Output file already exists, reading in and returning...'
   readcol,outfile,mod_theta,mod_w,format='F,D'
   return
ENDIF

;MAD Set constants/defaults
IF ~keyword_set(paramfile) THEN paramfile='default_params.ini'
IF ~keyword_set(zarray) THEN zarray=(findgen(400)/100.)+0.01
IF ~keyword_set(dndz) THEN dndz='dndz.txt'

c=2.99792458e5
IF ~keyword_set(h0) THEN h0=0.702
IF ~keyword_set(omega_m) THEN omega_m=0.275
IF ~keyword_set(omega_l) THEN omega_l=0.725
IF ~keyword_set(omega_b) THEN omega_b=0.046

;MAD Convert z array to comoving distance array
chi=fltarr(n_elements(zarray))
FOR i=0L,n_elements(zarray)-1 DO BEGIN
   d=cosmocalc(zarray[i],h=h0,om=omega_m,lambda=omega_l)
   chi[i]=d.d_c
ENDFOR

IF ~keyword_set(power_spec) THEN BEGIN
   root='camb'
   revz=reverse(zarray)
   numloop=floor((n_elements(revz)/150))
   IF (numloop NE 0) THEN BEGIN
      numrem=(n_elements(revz) MOD 150)
      FOR i=0L,numloop-1 DO BEGIN
         indx=indgen(150)+(i*150)
         matter_power_spec,paramfile,zarray[indx],h0=h0,omega_b=omega_b,$
                           omega_dm=omega_m-omega_b,omega_l=omega_l,maxk=50
      ENDFOR
      IF (numrem NE 0) THEN BEGIN
         matter_power_spec,paramfile,zarray[(n_elements(zarray)-numrem):n_elements(zarray)-1],$
                           h0=h0,omega_b=omega_b,omega_dm=omega_m-omega_b,omega_l=omega_l,maxk=50
      ENDIF
   ENDIF ELSE BEGIN
      matter_power_spec,paramfile,zarray,h0=h0,omega_b=omega_b,omega_dm=omega_m-omega_b,omega_l=omega_l,maxk=50
   ENDELSE
   combine_camb,'./',zarray,pspec,outfile='power_spec_camb.fits'
ENDIF ELSE BEGIN
;MAD Read in power spectrum if given file name instead of structure
   check=size(power_spec)
   IF ((check[0] EQ 0) AND (check[1] EQ 7)) THEN $
      pspec=mrdfits(power_spec,1) ELSE $
         pspec=power_spec
ENDELSE


;MAD Fit dndz if needed
IF (size(dndz,/type) EQ 7) THEN BEGIN
   readcol,dndz,zdist,format='D'
   fit_dndz,zdist,zarray,dndz
ENDIF ELSE BEGIN
   IF (n_elements(zarray) NE n_elements(dndz)) THEN $
      message,'Length of dndz and zarray don''t match!'
ENDELSE


   
;MAD convert to dimensionless power spectrum, put in factors of h
delsq=pspec.pk*(1./(2.*!dpi^2.))*pspec.k^3.
pkk=pspec.k*h0
pkz=pspec.z
;MAD Pull out unique k values
kvals=pkk[where(pkz EQ min(pkz))]

;MAD Convert theta to radians
thetarad=theta*(!dpi/180)

;MAD Define function of z array that will be integrated over
fz=fltarr(n_elements(where(pkk EQ pkk[0])))

;MAD Set dz/d\chi array
dzdchi=((h0*100.)/c)*(((omega_m*(1+zarray)^3.)+omega_l)^0.5)

;MAD Set output array
mod_w=fltarr(n_elements(theta))

;MAD Loop over theta values, integrate k then z
print,'MODEL_AUTOCORR: Integrating at each theta value...'
FOR i=0L,n_elements(thetarad)-1 DO BEGIN
   counter,i,n_elements(thetarad)
   FOR j=0L,n_elements(fz)-1 DO BEGIN
      ;MAD Get power spectrum at appropriate z
      xx=where(round(pkz*100.)/100. EQ round(zarray[j]*100.)/100.)
      ;MAD Interpolate on new grid of k, deltsq values for more accurate integral
      exponents=cgScaleVector(Findgen(10000),-5,max(alog10(kvals)))
      newk=10.^exponents
      newdelsq=interpol(delsq[xx],kvals,newk)
      ;MAD Don't keep interpolated values that go negative
      newk=newk[where(newdelsq GT 0)]
      newdelsq=newdelsq[where(newdelsq GT 0)]
      ;MAD Define Bessel function at each k
      J0=beselj(newk*thetarad[i]*chi[j],0,/double)
      ;MAD Integrate over k, put into array for integration over z
      fk=(newdelsq/newk)*J0*(dndz[j]^2)*dzdchi[j]*(1./newk)
      fz[j]=int_tabulated(newk,fk,/double)
   ENDFOR      
   ;MAD Integrate over z, multiply by pi
   IF ~keyword_set(bz) THEN mod_w[i]=!dpi*int_tabulated(zarray,fz,/double) ELSE $
      mod_w[i]=!dpi*int_tabulated(zarray,fz*(bz[0] + bz[1]*(1+zarray)^2),/double)
ENDFOR

IF keyword_set(plotfile) THEN BEGIN
   PS_start,filename=plotfile,xsize=8,ysize=6
   plot,theta,mod_w,psym=1,/xlog,/ylog,xra=[0.002,1],$
        yra=[0.0001,3],xsty=1,ysty=1,$
        xtit=textoidl('\theta'),ytit=textoidl('\omega(\theta)')
   PS_end,/png
ENDIF
   
IF keyword_set(outfile) THEN BEGIN
   openw,1,outfile
   printf,1,';theta (deg)    w'
   FOR i=0L,n_elements(mod_w)-1 DO BEGIN
      printf,1,theta[i],mod_w[i]
   ENDFOR
   close,1
ENDIF

return
END


