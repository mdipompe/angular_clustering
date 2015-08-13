;+
;  NAME:
;    model_autocorr
;  PURPOSE:
;    Generate a model projected angular autocorrelation of matter
;
;  USE:
;    model_cross_corr,theta,mod_w,power_spec='power_spec.fits',dndz='dndz.txt',zarray='z_chi.txt',$
;                      omega_m=omega_m,omega_l=omega_l,h0=h0,$
;                      outfile='model_autocorr.txt'
;
;  INPUT:
;    theta - Array of angular scales you want the model at (degrees)
;
;  Optional Inputs:
;    power_spec - structure or string name of fits file 
;                 containing the matter power
;                 spectrum.  Needs tags pk (the power spectrum), k
;                 (wavenumber), z (redshift).
;                 Can be made from CAMB data using camb4idl and
;                 combine_camb.pro. Defaults to 'power_spec_camb.fits'
;    dndz - string name of text file with dndz.  Will fit a spline
;           function and write it out to dndz_fit.txt.  You should
;           check this!!  Defaults to dndz.txt
;    zsample - string name of file containing z and chi values, as
;              made by chi_list.pro.  Defaults to z_chi.txt
;    omega_m - Omega_matter, defaults to 0.273.  Be sure it is
;              consistent with what you used when you calculated the
;              power spectrum!
;    omega_l - Omega_lambda, defaults to 0.727.  See above...
;    h0 - little h (H0/100), defaults to 0.702.  See above...
;    outfile - if supplied, writes model power out to text file
;    plotfile - if supplied, makes a plot of the model
;
;  KEYWORDS:
;    
;  OUTPUT:
;    mod_w - the model autocorrelation
;    dndz_fit.txt - a file containing the normalized dndz fit as a check
;
;  NOTES:
;    Make sure your cosmology agrees with how you calculated the power
;    spectrum, if you supply your own. Also check the dndz_fit.txt file to make sure nothing
;    went wrong!
;
;  HISTORY:
;    8-12-15 - Written - MAD (UWyo)
;-
PRO model_autocorr,theta,mod_w,power_spec=power_spec,dndz=dndz,$
                   zsample=zsample,omega_m=omega_m,omega_l=omega_l,h0=h0,$
                   outfile=outfile,plotfile=plotfile

;MAD If output file already exists, don't run just read it in
IF ~keyword_set(outfile) THEN check='' ELSE check=file_search(outfile)
IF (check NE '') THEN BEGIN
   print,'MODEL_AUTOCORR: Output file already exists, reading in and returning...'
   readcol,outfile,mod_theta,mod_w,format='F,D'
   return
ENDIF

;MAD Set constants/defaults
c=2.99792458e5
IF ~keyword_set(h0) THEN h0=0.702
IF ~keyword_set(omega_m) THEN omega_m=0.273
IF ~keyword_set(omega_l) THEN omega_l=0.727

IF ~keyword_set(power_spec) THEN power_spec='power_spec_camb.fits'
IF ~keyword_set(dndz) THEN dndz='dndz.txt'
IF ~keyword_set(zsample) THEN zsample='z_chi.txt'

;MAD Read in power spectrum if given file name instead of structure
check=size(power_spec)
IF ((check[0] EQ 0) AND (check[1] EQ 7)) THEN $
   pspec=mrdfits(power_spec,1) ELSE $
      pspec=power_spec

;MAD Read in z, comoving distances 
readcol,zsample,z,chi,format='D'

;MAD Fit dndz
readcol,dndz,zdist,format='D'
fit_dndz,zdist,z,dndz

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
dzdchi=((h0*100.)/c)*(((omega_m*(1+z)^3.)+omega_l)^0.5)

;MAD Set output array
mod_w=fltarr(n_elements(theta))

;MAD Loop over theta values, integrate k then z
print,'MODEL_AUTOCORR: Integrating at each theta value...'
FOR i=0L,n_elements(thetarad)-1 DO BEGIN
   counter,i,n_elements(thetarad)
   FOR j=0L,n_elements(fz)-1 DO BEGIN
      ;MAD Get power spectrum at appropriate z
      xx=where(round(pkz*100.)/100. EQ round(z[j]*100.)/100.)
      ;MAD Interpolate on new grid of k, deltsq values for more accurate integral
      exponents=cgScaleVector(Findgen(10000),-5,max(alog10(kvals)))
      newk=10.^exponents
      newdelsq=interpol(delsq[xx],kvals,newk)
      ;MAD Define Bessel function at each k
      J0=beselj(newk*thetarad[i]*chi[j],0,/double)
      ;MAD Integrate over k, put into array for integration over z
      fk=(newdelsq/newk)*J0*(dndz[j]^2)*dzdchi[j]*(1./newk)
      fz[j]=int_tabulated(newk,fk,/double)
   ENDFOR      
   ;MAD Integrate over z, multiply by pi
   mod_w[i]=!dpi*int_tabulated(z,fz,/double)
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

END


