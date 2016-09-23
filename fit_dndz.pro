;+
;  NAME:
;    fit_dndz
;  PURPOSE:
;    Take a z distribution and fit a model dndz
;
;  USE:
;    fit_dndz,realz,sampledz,dndz,binsize=binsize,outfile=outfile
;
;  INPUT:
;    realz - the z values of the real data set to fit
;    sampledz - the z values you want to know the fit at
;
;  Optional Inputs:
;    binsize - the bin size of the histogram to fit (default 0.2)
;    outfile - name of file with z and the fit (for testing/checks)
;    sigma - "tension" in the spline fit.  0-10, defaults to 1
;
;  OUTPUT:
;    dndz - the values of the fit at sampledz
;    dndz_fit.txt - a text file with z and the fit, so you can check
;                   that things worked if you want.
;
;  HISTORY:
;    11-7-14 - Written - MAD (UWyo)
;    11-18-15 - Fixed accuracy bug in dndz normalization for small
;               bins - MAD (Dartmouth)
;     9-15-16 - Added output file option - MAD (Dartmouth)
;-
PRO fit_dndz,realz,sampledz,dndz,binsize=binsize,outfile=outfile,sigma=sigma

;MAD Set defaults
IF ~keyword_set(binsize) THEN binsize=0.2
IF ~keyword_set(outfile) THEN outfile='dndz_fit.txt'
IF ~keyword_set(sigma) THEN sigma=1

;MAD Bin the real data
h=histogram(realz,binsize=binsize,min=0,max=max(realz))

;MAD Make array of bin centers
x=fltarr(ceil((max(realz)/binsize))+1)
x[0]=binsize/2.
i=1
WHILE (max(x) LT max(realz)) DO BEGIN
 x[i]=x[i-1]+binsize
 i=i+1
ENDWHILE

;MAD Smooth over gaps in the histogram
xx=where(h NE 0)
h=h[xx]
x=x[xx]

;MAD Force fit to N=0 at z=0
x=[0,x]
h=[0,h]

;MAD Fit cubic spline, normalize
fit=spline(x,h,sampledz,sigma)
fit[where((sampledz GT max(realz)) OR (sampledz LT min(realz)))]=0
area=int_tabulated(sampledz,fit,/double)
dndz=fit/area

openw,lun,outfile,/get_lun
printf,lun,';z     fit    dndz (normalized fit)'
FOR i=0,n_elements(sampledz)-1 DO BEGIN
  printf,lun,sampledz[i],fit[i],dndz[i],format='(F,1x,F,1x,F)'
ENDFOR
close,lun

return
END
