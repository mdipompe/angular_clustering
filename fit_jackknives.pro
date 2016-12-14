;+
;  NAME:
;    fit_jackknives
;
;  PURPOSE:
;    Randomly sample clustering amplitudes from jackknife iterations
;    and fit the bias, to build a distribuiton of bias values.
;    Cross-between usual and Bayesian fitting, without relying on a
;    prior. Relies on several files from ang_cluster.pro and
;    ang_jackknife.pro, which can be names explicitly if desired.
;
;  USE:
;    fit_jackknives,theta,w_theta,b,[optional inputs]

;  INPUTS:
;    theta - bin centers in angular scale
;    w_theta - autocorrelation values
;
;  OPTIONAL INPUTS:
;    n - number of random draws/fits (default 10000)
;    minscale - set lower limit for fitting range (default min(theta))
;    maxscale - set upper limit for fitting range (default max(theta)
;    datafile - fits file of the data (needed to get N_data for
;               Poisson noise; default data_reg.fits)
;    data2file - fits file of second data set (in cross-correlation)
;    randfile - fits file of the randoms (needed to get N_rand for
;               width of jackknife draws; default rand_reg.fits)
;    dmfile - name of file with DM model autocorrelation. 
;             Two columns, angular scale (in deg) and w_theta  
;    jackkfile - text file with jackknife results (3 columns, theta,
;                w, iteration; default jackknife_results.txt)
;    dcountfile - text file with DD counts (default DD_DR.txt)
;    rcountfile - text file with RR counts (default RR.txt)
;    rregfile - text file with RR counts per jackknife iteration
;               (default rr_reg.txt')
;    path - path to directory that has files (makes easier to use
;           default names; defaults to current directory)
;    acc - accuracy with which to do the bias fitting in
;          fit_angular_cluster (more accurate=slower; default 1e-3)
;    outfile - if set, will print the bias and error of each fit
;              iteration to this file
;    seed - the random seed for the random jackknife iteration draw
;           (default 818)
;
;  OUTPUTS:
;    bias - distribution of bias values
;
;  HISTORY:
;    12-8-16 - Written - MAD (Dartmouth)
;-
PRO fit_jackknives,theta,w_theta,b,n=n,minscale=minscale,maxscale=maxscale,$
                   datafile=datafile,$
                   data2file=data2file,$
                   dmfile=dmfile,$
                   randfile=randfile,jackfile=jackfile,$
                   covfile=covfile,$
                   dcountfile=dcountfile,$
                   dregfile=dregfile,$
                   rcountfile=rcountfile,rregfile=rregfile,$
                   path=path,acc=acc,outfile=outfile,seed=seed,$
                   poisson=poisson
  
  ;MAD Start timer
  st=timer()

  ;MAD Set defaults
  IF (n_elements(n) EQ 0) THEN n=1000
  IF (n_elements(n) EQ 0) THEN acc=1e-3
  IF (n_elements(seed) EQ 0) THEN seed=818
  IF (n_elements(path) EQ 0) THEN path='./'
  IF (n_elements(datafile) EQ 0) THEN datafile=path+'data_reg.fits'
  IF (n_elements(randfile) EQ 0) THEN datafile=path+'rand_reg.fits'
  IF (n_elements(jackfile) EQ 0) THEN jackfile=path+'jackknife_results.txt'
  IF (n_elements(dcountfile) EQ 0) THEN dcountfile=path+'DD_DR.txt'
  IF (n_elements(dregfile) EQ 0) THEN dregfile=path+'dd_reg.txt'
  IF (n_elements(rcountfile) EQ 0) THEN rcountfile=path+'RR.txt'
  IF (n_elements(rregfile) EQ 0) THEN rregfile=path+'rr_reg.txt'
  IF (~keyword_set(poisson) AND n_elements(covfile) EQ 0) THEN covfile=path+'covariance.txt'

  ;MAD read in necessary data, convert as necessary
  data=mrdfits(datafile,1)
  IF (n_elements(data2file) NE 0) THEN data2=mrdfits(data2file,1)
  rand=mrdfits(randfile,1)
  readcol,jackfile,thall,wall,pix,format='D'
  readcol,dcountfile,dd,format='D',numline=n_elements(theta)
  IF (n_elements(data2) NE 0) THEN dd=dd*n_elements(data)*n_elements(data2) ELSE $
     dd=dd*(n_elements(data)^2.)
  ddreg=read_matrix(dregfile)
  readcol,rcountfile,rr,format='D',numline=n_elements(theta)
  rr=rr*(n_elements(rand)^2.)
  rrreg=read_matrix(rregfile)

  IF ~keyword_set(poisson) THEN C=read_matrix(covfile)
  
  ;MAD Set default scales (full range)
  IF (n_elements(minscale) EQ 0) THEN minscale=min(theta)/60.
  IF (n_elements(maxscale) EQ 0) THEN minscale=max(theta)/60.

  ;MAD Find unique angular bins
  th=thall[rem_dup(thall)]

  ;Initialize outputs
  b=fltarr(n)
  b_err=fltarr(n)

  ;Start n iterations
  FOR i=0,n-1 DO BEGIN
     counter,i,n
     w=dblarr(n_elements(th))
     err=dblarr(n_elements(th))
     ;MAD pick a w for each bin randomly from jackknives
     FOR j=0,n_elements(th)-1 DO BEGIN
        xx=where(thall EQ th[j])
        ;Have to shift to 0, widen distribution, shift back
        ws=((wall[xx]-w_theta[j])*sqrt(total((rr[j]-rrreg[j,*])/rr[j])))+w_theta[j]
        indx=floor(randomu(seed)*max(pix))
        w[j]=ws[indx]
        ;MAD Get Poisson error for this point
        IF (keyword_set(poisson)) THEN BEGIN
           IF (n_elements(data2) EQ 0) THEN factor=2. ELSE factor=1.
           err[j]=sqrt((factor*(1+w[j])^2.)/(dd[j]-ddreg[j,indx]))
        ENDIF
     ENDFOR
     ;MAD Fit random iteration
     IF (keyword_set(poisson)) THEN $
        fit_ang_cluster,th,w,err,bias,bias_err,minscale=minscale,maxscale=maxscale,/fitbias,$
                        dmfile=dmfile,acc=acc,/silent,/nocovar ELSE $
                           fit_ang_cluster,th,w,err,bias,bias_err,minscale=minscale,maxscale=maxscale,$
                                           /fitbias,dmfile=dmfile,filepath=path,acc=acc,/silent
     b[i]=bias
     b_err[i]=bias_err
     
     ;MAD Hack, because readcol will die after many
     ;MAD iterations because it doesn't close luns for some reason
     close,/all
  ENDFOR

  ;MAD Write out, if needed
  IF (n_elements(outfile) NE 0) THEN BEGIN
     openw,lun,outfile,/get_lun
     printf,lun,';b       b_err',format='(A)'
     FOR i=0L,n-1 DO $
        printf,strtrim(b[i],2) + '     ' + strtrim(b_err[i],2),format='(A)'
     free_lun,lun
  ENDIF
  
  et=timer(st=st,/fin)
  return
END
