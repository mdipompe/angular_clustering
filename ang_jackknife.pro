;+
;  NAME:
;    ang_jackknife
;
;  PURPOSE:
;    This procedure will take the outputs of ang_cluster.pro (with
;    jackknife flag set), and jackknife the errors using the covariance matrix.
;
;  USE:
;    ang_jackknife,data,rand,theta_full,w_theta_full,errors $
;                [,data2=data2,maxscale=maxscale,outfile='results_with_errs.txt',bins=bins]
;
;  INPUTS:
;    data - The real data structure, with the tags ra, dec, and
;           reg. (Made by split_regions.pro)
;    rand - The random structure, with the tags ra, dec,
;           and reg. (Made by split_regions.pro)
;    theta_full - the angular bin centers from the full measurement
;    w_theta_full - the correlation funciton using the full sample          
;
;  OPTIONAL INPUTS:
;    data2 - second data set for cross-correlation, with tags ra,dec,
;            and reg.  (Made by split regions.pro)
;    maxscale - The maximum angular scale of interest, in degrees.  Should be
;               the same as what was used with ang_cluster.pro.  Defaults to 2
;    bins - options for binning include 3, 4 or 5 (refers to bins/dex)
;           default is 5.  Can also set to 0 to make two large bins.
;
;    outfile - string name of output, will have columns of theta,
;              w_theta, error
;   
;  OUTPUTS:
;    errors - the 1-sigma errors (diagonal of covariance matrix)
;    jackknife_results.txt - Three columns; bin center (degrees), w_theta,
;                            and pixel/region excluded from analysis
;    covariance.txt - A file with the covariance matrix.
;
;  HISTORY:
;    2013 - Written - MAD (UWyo) 
;    12-1-13 - First combined version - MAD (UWyo)
;    3-27-15 - Reworked and cleaned - MAD (UWyo)
;    7-20-15 - Added cross-correlation capability - MAD (UWyo)
;-
PRO ang_jackknife,data,rand,theta_full,w_theta_full,errors,data2=data2,maxscale=maxscale,outfile=outfile,bins=bins

IF (n_elements(data) EQ 0) THEN message,'Syntax - ang_jackknife,data,rand,theta_full,w_theta_full,errors[,data2=data2,maxscale=maxscale,outfile=''results_with_errs.txt'',bins=bins]'

;MAD Get the start time
st=systime(1)

;MAD Set default max scale
IF ~keyword_set(maxscale) THEN maxscale=2.

;MAD Set dataset sizes
n_data_tot=double(n_elements(data))
n_rand_tot=double(n_elements(rand))
IF keyword_set(data2) THEN n_data2_tot=double(n_elements(data2))

;MAD Set binning (5dex is the default), make bin edges, centers
IF ~keyword_set(bins) THEN bins=5
IF (bins NE 0) THEN BEGIN
   bininfo=define_bins(bins,maxscale)
   bin_edge=bininfo.edges
   bin_cent=bininfo.cents
ENDIF

;MAD Define big bins, if needed
IF (bins EQ 0) THEN BEGIN
   bin_edge=[0.2,6.,maxscale*60.]
   bin_cent=[3.1,((6.+(maxscale*60.))/2.),(maxscale*60.)+30.]
ENDIF


;MAD Open file to write jackknife results
openw,1,'jackknife_results.txt'

;MAD Read in RR data from ang_cluster.pro 
print,'Ang_jackknife - Reading counts/region files...'
readcol,'RR.txt',rr_tot,format='D'
readcol,'rr_reg.txt',rr1,rr2,rr3,rr4,rr5,rr6,rr7,rr8,rr9,rr10,rr11,$
        rr12,rr13,rr14,rr15,rr16,format='D'
;MAD Multiply original measurements by total to un-normalize the values
rr_tot=rr_tot*n_rand_tot*n_rand_tot

;MAD Initialize an array to store the number of random points used in
;each iteration, for use in calculating covariance later
n_rand_pix=dblarr(max(rand.reg))

;MAD Read in DD/DR data from ang_cluster.pro
readcol,'DD_DR.txt',dd_tot,dr_tot,format='D,D'
readcol,'dd_reg.txt',dd1,dd2,dd3,dd4,dd5,dd6,dd7,dd8,dd9,dd10,dd11,$
                     dd12,dd13,dd14,dd15,dd16,format='D'
readcol,'dr_reg.txt',dr1,dr2,dr3,dr4,dr5,dr6,dr7,dr8,dr9,dr10,dr11,$
                     dr12,dr13,dr14,dr15,dr16,format='D'

;MAD Multiply original measurements by total to un-normalize the values
IF ~keyword_set(data2) THEN dd_tot=dd_tot*n_data_tot*n_data_tot ELSE $
   dd_tot=dd_tot*n_data_tot*n_data2_tot   
dr_tot=dr_tot*n_data_tot*n_rand_tot


;MAD Start loop over included regions
print,'Ang_jackknife - Looping over regions...'
FOR i=1L,max(rand.reg) DO BEGIN
   use_data=data[where(data.reg NE i)]
   use_rand=rand[where(rand.reg NE i)]
   IF keyword_set(data2) THEN use_data2=data2[where(data2.reg NE i)]

   ;MAD Get number in each sample in double format
   n_data=double(n_elements(use_data))
   n_rand=double(n_elements(use_rand))
   n_rand_pix[i-1]=n_rand
   IF keyword_set(data2) THEN n_data2=double(n_elements(use_data2))

   ;MAD Use info from ang_cluster.pro to do DD/DR/RR quickly 
   rr_use='rr'+strtrim(i,2)
   cmd='h_rr = rr_tot - '+rr_use
   R=execute(cmd)
   h_rr=h_rr*(1./(n_rand*n_rand))
   
   dd_use='dd'+strtrim(i,2)
   cmd='h_dd = dd_tot - '+dd_use
   R=execute(cmd)
   IF ~keyword_set(data2) THEN h_dd=h_dd*(1./(n_data*n_data)) ELSE $
      h_dd=h_dd*(1./(n_data*n_data2))

   dr_use='dr'+strtrim(i,2)
   cmd='h_dr = dr_tot - '+dr_use
   R=execute(cmd)
   h_dr=h_dr*(1./(n_data*n_rand))

   ;MAD Calculate autocorrelation
   IF ~keyword_set(data2) THEN $
      w_theta=(1./h_rr)*(h_dd-(2.*h_dr)+h_rr) ELSE $
         w_theta=(h_dd/h_dr)-1.

   ;If the scale of interest is slightly larger than the maximum edge of
   ;the bins, the last value is nonsense
   IF (maxscale*60. GT max(bin_edge)) THEN BEGIN
      bin_cent=bin_cent[0:n_elements(w_theta)-2]
      w_theta=w_theta[0:n_elements(w_theta)-2]
   ENDIF
 
   xx=where(finite(w_theta) EQ 0,cnt)
   IF (cnt NE 0) THEN BEGIN
      w_theta[xx] = -9999
      bin_cent[xx] = -9999
   ENDIF

   FOR k=0L,n_elements(w_theta)-1 DO BEGIN
      printf,1,bin_cent[k],w_theta[k],i
   ENDFOR

   print,'Ang_jackknife - finished region ',strtrim(i,2)
ENDFOR
close,1


;MAD Calculate covariance matrix 
print,'Ang_jackknife - calculating covariance matrix...'
readcol,'jackknife_results.txt',theta_pix,w_theta_pix,pix,format='D,D,F'

num=n_elements(theta_full)
readcol,'rr_reg.txt',rr1,rr2,rr3,rr4,rr5,rr6,rr7,rr8,rr9,rr10,rr11,$
        rr12,rr13,rr14,rr15,rr16,format='D',numline=num
readcol,'RR.txt',rr_tot,format='D',numline=num
rr_tot_norm=rr_tot
rr_tot=rr_tot*n_rand_tot*n_rand_tot

C=dblarr(n_elements(theta_full),n_elements(theta_full))
diag=dblarr(n_elements(theta_full))
FOR k=1,max(pix) DO BEGIN
   string='rr_pix = rr'+strtrim(k,2)
   R=Execute(string)
   rr_pix=rr_tot-rr_pix
   w_theta=w_theta_pix[where(pix EQ k)]
   temp1=(SQRT(rr_pix*(1./rr_tot)))*(w_theta-w_theta_full)
   temp2=(SQRT(rr_pix*(1./rr_tot)))*(w_theta-w_theta_full)
   C= C + (temp1 # temp2)
ENDFOR

;MAD Write out covariance matrix file for fitting
print,'Ang_jackknife - writing out covariance matrix...'
openw,1,'covariance.txt'
formatstring=strarr(n_elements(theta_full))+'D,1x,'
formatstring[0]='(D,1x,'
formatstring[n_elements(formatstring)-1]='D)'
formatstring=strjoin(formatstring)
FOR i=0,n_elements(theta_full)-1 DO BEGIN
   printf,1,C[i,*],format=formatstring
ENDFOR
close,1

;MAD Pull out the diagonals of the covariance matrix
FOR i=0L,n_elements(theta_full)-1 DO BEGIN
   FOR j=0L,n_elements(theta_full)-1 DO BEGIN
      IF (i EQ j) THEN diag[i] = C[i,j]
   ENDFOR
ENDFOR

;MAD Errors are diagonal of covariance matrix
errors=SQRT(diag)

;MAD Calculate regression matrix (in case you want it as a sanity check)
reg=dblarr(n_elements(theta_full),n_elements(theta_full))
FOR i=0,n_elements(theta_full)-1. DO BEGIN
   FOR j=0,n_elements(theta_full)-1. DO BEGIN
      reg[i,j]=c[i,j]/(SQRT(diag[i])*SQRT(diag[j]))
   ENDFOR
ENDFOR

;MAD Write out final file with 1-sigma errors from diagonals
IF keyword_set(outfile) THEN BEGIN
   print,'Ang_jackknife - Writing final results file with standard errors...'
   openw,1,outfile
   FOR i=0,n_elements(diag)-1 DO BEGIN
      printf,1,theta_full[i],w_theta_full[i],errors[i],format='(F,1x,D,1x,D)'
   ENDFOR
   close,1
ENDIF


;MAD Get finish, elapsed time
et=systime(1)
print,'Ang_jackknife - Elapsed time= ',strtrim((et-st)/60,2),' min'

return
END


