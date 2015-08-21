;+
;  NAME:
;    fig_ang_cluster
;
;  PURPOSE:
;    This procedure fits a fixed slope and full power-law to angular
;    clustering measurements, first using mpfitfun (requires 
;    'powerlaw_fixedslope.pro' and 'powerlaw.pro' functions), and then using the
;    covariance matrix.  Plots to screen or file, prints fit parameters to
;    screen.  Expects there to be a file called
;    'jackknife_results.txt' and 'covariance.txt'
;
;  USE:
;    fit_ang_cluster,theta,w_theta,errors,minscale=minscale,maxscale=maxscale,$
;                   /fitplaws,/fitbias,dmfile='model_file.txt',$
;                   plawplot='plawplot.png',biasplot='biasplot.png'
;
;  INPUTS:
;    theta - bin centers in angular scale
;    w_theta - autocorrelation values
;    errors - 1sigma errors on wtheta
;
;  KEYWORDS:
;    fitplaws - set to fit power laws (free slope, fixed slope)
;    fitbias - set to fit model bias (need to set dmfile to model file)
;
;  OPTIONAL INPUTS:
;    minscale - set lower limit for fitting range (default 0.0027 deg)
;    maxscale - set upper limit for fitting range (default 2 deg)
;    plawplot - name of plot file for power law fits
;    biasplot - name of plot file for bias fits
;    dmfile - name of file with DM model autocorrelation. 
;             Two columns, angular scale (in deg) and w_theta
;    filepath - path to location of files like, covariance file.
;               Defaults to current directory.
;
;  OUTPUTS:
;    bias 
;    biaserr
;
;  HISTORY:
;    2013 - Written - MAD (UWyo) 
;    12-1-13 - First combined version - MAD (UWyo)
;    3-27-15 - Reworked and cleaned - MAD (UWyo)
;-
PRO fit_ang_cluster,theta,w_theta,errors,bias,biaserr,minscale=minscale,maxscale=maxscale,fitplaws=fitplaws,fitbias=fitbias,dmfile=dmfile,plawplot=plawplot,biasplot=biasplot,filepath=filepath

IF (n_elements(theta) EQ 0) THEN message,'Syntax - fit_ang_cluster,theta,w_theta,errors[,minscale=minscale,maxscale=maxscale,/fitplaws,/fitbias,plawplot=''plawplot.png'',biasplot=''biasplot.png'']'

;MAD Get start time
st=systime(1)

;MAD Set default min/max scale
IF ~keyword_set(minscale) THEN minscale=0.0027
IF ~keyword_set(maxscale) THEN maxscale=2.

;MAD Set default file path
IF ~keyword_set(filepath) THEN filepath='./'

;MAD Define filled circle for plot
circsym

;MAD Read in covariance matrix
print,'fit_ang_cluster - reading in covariance matrix...'
read_square_matrix,n_elements(theta),filepath+'covariance.txt',C

;MAD Limit to regions of interest (set by min/max scale keywords)
inscale=where((theta GT minscale*60.) AND (theta LT maxscale*60.))
bin=theta[inscale]
wtheta=w_theta[inscale]
errors=errors[inscale]
C=C[inscale,*]
C=C[*,inscale]

;MAD Look at jackknife iterations, don't fit bins where these might go negative
fitflag=0
w_theta_fit=wtheta
readcol,filepath+'jackknife_results.txt',jackbin,jackw_theta,jackpix,format='D,D,F'
FOR i=1,max(jackpix) DO BEGIN
   xx=where(jackpix EQ i)
   tempbin=jackbin[xx]
   tempw=jackw_theta[xx]
   xx=where((tempbin GT minscale*60.) AND (tempbin LT maxscale*60.))
   tempbin=tempbin[xx]
   tempw=tempw[xx]
   yy=where(tempw LE 0)
   IF (yy[0] NE -1) THEN BEGIN
      w_theta_fit[yy] = -9999
      fitflag=1
   ENDIF
ENDFOR

;MAD limit fitting to values of w_theta > 0 in all jackknives
okfit=where(w_theta_fit GT 0.)


;MAD Fit power laws if keyword set
IF keyword_set(fitplaws) THEN BEGIN
   ;MAD Fit power-law with mpfitfun to generate an initial guess
   print,'fit_ang_cluster - getting intial power-law parameter guesses...'
   guess=[10.,-1.]
   fit_plaw=mpfitfun('powerlaw',bin[okfit]/60.,wtheta[okfit],errors[okfit],guess,perror=fit_plaw_errs,/quiet)
   x_plaw=findgen(100)+0.00001
   y_plaw=fit_plaw[0]*x_plaw^(fit_plaw[1])

   ;MAD Generate array of A and delta guesses around initial power-law points
   min_A_plaw=fit_plaw[0]-(1.5*fit_plaw[0])
   max_A_plaw=fit_plaw[0]+(1.5*fit_plaw[0])
   A_guess_plaw=(findgen(10000.)*((max_A_plaw-min_A_plaw)/10000.))+min_A_plaw
   min_d_plaw=fit_plaw[1]-(1.5*fit_plaw[1])
   max_d_plaw=fit_plaw[1]+(1.5*fit_plaw[1])
   d_guess_plaw=(findgen(10000.)*((max_d_plaw-min_d_plaw)/10000.))+min_d_plaw

   ;MAD Fit fixed slope power-law with mpfitfun for initial guess
   print,'fit_ang_cluster - getting initial fixed-slope power-law guesses...'
   guess=[10.]
   fit_fixed=mpfitfun('powerlaw_fixedslope',bin[okfit]/60.,wtheta[okfit],errors[okfit],guess,perror=fit_fixed_errs,/quiet)
   x_fixed=findgen(100)+0.00001
   y_fixed=fit_fixed[0]*x_fixed^(-1.)

   ;MAD Generate array of A guesses around initial fixed-slope point
   min_A_fixed=fit_fixed[0]-(1.5*fit_fixed[0])
   max_A_fixed=fit_fixed[0]+(1.5*fit_fixed[0])
   A_guess_fixed=(findgen(20000.)*((max_A_fixed-min_A_fixed)/20000.))+min_A_fixed

   ;MAD Initialize arrays of Chi^2 values to fill; invert covariance matrix
   chisq_plaw=dblarr(n_elements(A_guess_plaw),n_elements(d_guess_plaw))
   chisq_fixed=dblarr(n_elements(A_guess_fixed))
   C_inv=invert(C,/double)
   C_inv_use=C_inv[*,okfit]
   C_inv_use=C_inv_use[okfit,*]

   ;;MAD Fit using just the variance terms of the covariance matrix
   ;FOR i=0L,n_elements(A_guess_plaw)-1 DO BEGIN
   ; FOR j=0L,n_elements(d_guess_plaw)-1 DO BEGIN
   ;  FOR k=0L,n_elements(bin[okfit])-1 DO BEGIN
   ;   val=A_guess_plaw[i]*((bin[okfit[k]]/60.)^d_guess_plaw[j])
   ;   temp=(wtheta[okfit[k]]-val) * C_inv_use[k,k] * (wtheta[okfit[k]]-val)
   ;   chisq_plaw[i,j]=chisq_plaw[i,j]+temp
   ;  ENDFOR
   ; ENDFOR
   ;ENDFOR

   ;MAD Loop over power-law guesses to fill chi^2 matrix to fit using covariance
   print,'fit_ang_cluster - Building power-law chi^2 matrix...'
   FOR i=0L,n_elements(A_guess_plaw)-1 DO BEGIN
      FOR j=0L,n_elements(d_guess_plaw)-1 DO BEGIN
         val=A_guess_plaw[i]*((bin[okfit]/60.)^d_guess_plaw[j])
         temp=(wtheta[okfit]-val) # C_inv_use # (wtheta[okfit]-val)
         chisq_plaw[i,j]=chisq_plaw[i,j]+temp
      ENDFOR
   ENDFOR

   ;MAD Plot the chi^2 countours (at intervals of min(chi^2) + 10%)
   ;contour,chisq_plaw,A_guess_plaw,d_guess_plaw,charsize=2,xtit='A',ytit='delta',xra=[0.,0.005],yra=[-1.6,-0.5],levels=[min(chisq_plaw),min(chisq_plaw)+0.1*min(chisq_plaw),min(chisq_plaw)+0.2*min(chisq_plaw),min(chisq_plaw)+0.3*min(chisq_plaw),min(chisq_plaw)+0.4*min(chisq_plaw),min(chisq_plaw)+0.5*min(chisq_plaw),min(chisq_plaw)+0.6*min(chisq_plaw),min(chisq_plaw)+0.7*min(chisq_plaw),min(chisq_plaw)+0.8*min(chisq_plaw),min(chisq_plaw)+0.9*min(chisq_plaw)]

   ;MAD Find where the power-law chi^2 value is minumum
   print,'fit_ang_cluster - finding min chi^2 and best fit power-law params...'
   minchi_plaw=where(chisq_plaw EQ min(chisq_plaw))
   s=size(chisq_plaw)
   ncol=s[1]
   col=minchi_plaw MOD ncol
   row=minchi_plaw/ncol

   ;MAD get best A and g values, define lines for plots
   best_A_plaw=A_guess_plaw[col]
   best_d_plaw=d_guess_plaw[row]
   x_covar_plaw=findgen(100)+0.00001
   y_covar_plaw=best_A_plaw[0]*x_covar_plaw^(best_d_plaw[0])

   ;MAD Find errors on A and g (using delta chi^2 = 2.3, or closest to it)
   chi_col=chisq_plaw[*,row]-min(chisq_plaw)
   chi_col1=chi_col[0:col]
   chi_col2=chi_col[col:n_elements(chi_col)-1]
   A_vals1=A_guess_plaw[0:col]
   A_vals2=A_guess_plaw[col:n_elements(A_guess_plaw)-1]
   xx=closest(chi_col1,2.3)
   A_err1=A_vals1[xx]
   yy=closest(chi_col2,2.3)
   A_err2=A_vals2[yy]
   A_plaw_err=(abs(best_A_plaw[0]-A_err1)+abs(best_A_plaw[0]-A_err2))*(0.5)
   dchi_A_plaw=(chi_col1[xx]+chi_col2[yy])*0.5

   chi_row=chisq_plaw[col,*]-min(chisq_plaw)
   chi_row1=chi_row[0:row]
   chi_row2=chi_row[row:n_elements(chi_row)-1]
   d_vals1=d_guess_plaw[0:row]
   d_vals2=d_guess_plaw[row:n_elements(d_guess_plaw)-1]
   xx=closest(chi_row1,2.3)
   d_err1=d_vals1[xx]
   yy=closest(chi_row2,2.3)
   d_err2=d_vals2[yy]
   d_plaw_err=(abs(best_d_plaw[0]-d_err1)+abs(best_d_plaw[0]-d_err2))*(0.5)
   dchi_d_plaw=(chi_row1[xx]+chi_row2[yy])*0.5


   ;MAD Loop over fixed-splope power-law guesses to fill chi^2 array
   print,'fit_ang_cluster - Building power-law fixed chi^2 array...'
   FOR i=0L,n_elements(A_guess_fixed)-1 DO BEGIN
      val=A_guess_fixed[i]*((bin[okfit]/60.)^(-1.))
      temp=(wtheta[okfit]-val) # C_inv_use # (wtheta[okfit]-val)
      chisq_fixed[i]=chisq_fixed[i]+temp
   ENDFOR

   ;MAD Plot the chi^2 values vs A
   ;plot,A_guess_fixed,chisq_fixed,linestyle=1,xra=[0.,0.004],yra=[0,10000],xtit='A',ytit='chi^2'

   ;MAD Find where the fixed-slope chi^2 values is minumum, get best A value
   print,'fit_ang_cluster - finding min chi^2 and best fit fixed power-law amplitude...'
   minchi_fixed=where(chisq_fixed EQ min(chisq_fixed))
   best_A_fixed=A_guess_fixed[minchi_fixed]
   x_covar_fixed=findgen(100)+0.00001
   y_covar_fixed=(best_A_fixed[0])*((x_covar_fixed)^(-1.))

   ;MAD find errors on best fixed A, using delta chi^2 = 1
   chi1=chisq_fixed[0:minchi_fixed]-min(chisq_fixed)
   chi2=chisq_fixed[minchi_fixed:n_elements(chisq_fixed)-1]-min(chisq_fixed)
   A_vals1=A_guess_fixed[0:minchi_fixed]
   A_vals2=A_guess_fixed[minchi_fixed:n_elements(A_guess_fixed)-1]
   xx=closest(chi1,1.)
   A_err1=A_vals1[xx]
   yy=closest(chi2,1.)
   A_err2=A_vals2[yy]
   A_fixed_err=(abs(best_A_fixed[0]-A_err1)+abs(best_A_fixed[0]-A_err2))*(0.5)
   dchi_A_fixed=(chi1[xx]+chi2[yy])*0.5

   ;MAD Set any w_theta values < 0 to a very small number (so it can be
   ;plotted in log scale)
   IF (total(where(wtheta LT 0)) NE -1) THEN wtheta[where(wtheta LT 0)]=1e-6


   ;MAD Plot results (if keyword set)
   IF (keyword_set(plawplot)) THEN PS_start,filename=plawplot
   xtit=textoidl('\theta (deg)')
   ytit=textoidl('\omega_{\theta}')
   temp=textoidl('\delta')
   legendstring1='mpfit fixed '+textoidl('\delta=-1')
   legendstring2='covar fixed '+textoidl('\delta=-1')
   legendstring3='mpfit full power-law'
   legendstring4='covar full power-law'
   nice_plot,0.002,5.0,0.0005,10.,xtit=xtit,ytit=ytit,/xlog,/ylog
   oplot,bin/60.,wtheta,psym=-8,color=cgcolor('blue')
   oploterror,bin/60.,wtheta,errors,psym=3,color=cgcolor('blue')
   oplot,x_fixed,y_fixed,linestyle=1
   oplot,x_plaw,y_plaw,linestyle=2
   oplot,x_covar_plaw,y_covar_plaw,linestyle=2,color=cgcolor('red')
   oplot,x_covar_fixed,y_covar_fixed,linestyle=1,color=cgcolor('red')
   legend,[legendstring1,legendstring2,legendstring3,legendstring4],$
          linestyle=[1,1,2,2],$
          color=[cgcolor('black'),cgcolor('red'),cgcolor('black'),cgcolor('red')],$
          /top,/right,box=0,charsize=1.2,charthick=1.2
   IF (keyword_set(plawplot)) THEN PS_end,/png


   ;MAD Print results to screen
   print,' '

   IF (fitflag EQ 1) THEN print,'***************************************************************************'
   IF (fitflag EQ 1) THEN print,'***SOME POINTS WERE LEFT OUT OF FIT DUE TO NEGATIVE W_THETA IN JACKKNIFE***'
   IF (fitflag EQ 1) THEN print,'***************************************************************************'

   print,'Fixed slope covar fit errors are based on a delta chi^2 of ',strtrim(dchi_A_fixed,2),' (for A)'
   print,'Power-law covar fit errors are based on a delta chi^2 of ',strtrim(dchi_A_plaw,2),' (for A)'
   print,'Power-law covar fit errors are based on a delta chi^2 of ',strtrim(dchi_d_plaw,2),' (for delta)'

   print,' '

   print,'Fits from mpfitfun ('+strtrim(minscale,2)+' - '+strtrim(maxscale,2)+' degrees):'
   print,'Fixed slope power-law amplitude A = ',strtrim(fit_fixed[0],2),' +/- ',strtrim(fit_fixed_errs[0],2)
   print,'Full power-law amplitude A = ',strtrim(fit_plaw[0],2),' +/- ',strtrim(fit_plaw_errs[0],2)
   print,'Full power-law exponent delta = ',strtrim(fit_plaw[1],2),' +/- ',strtrim(fit_plaw_errs[1],2)

   print,' '

   print,'Fits using covariance ('+strtrim(minscale,2)+' - '+strtrim(maxscale,2)+' degrees):'
   print,'Fixed slope power-law amplitude A = ',strtrim(best_A_fixed[0],2),' +/- ',strtrim(A_fixed_err,2)
   print,'Full power-law amplitude A = ',strtrim(best_A_plaw[0],2),' +/- ',strtrim(A_plaw_err,2)
   print,'Full power-law exponent delta = ',strtrim(best_d_plaw[0],2),' +/- ',strtrim(d_plaw_err,2)

   print,' '

   ;MAD If chi^2 minimization went all the way to the edge of the input
   ;grids, give a warning
   IF ((best_A_fixed[0] EQ min(A_guess_fixed)) OR (best_A_fixed[0] EQ max(A_guess_fixed))) THEN $
      print,'WARNING: Maybe didn''t explore enough parameter space for fixed-slope A!!'
   IF ((best_A_plaw[0] EQ min(A_guess_plaw)) OR (best_A_plaw[0] EQ max(A_guess_plaw))) THEN $
      print,'WARNING: Maybe didn''t explore enough parameter space for power-law A!!'
   IF ((best_d_plaw[0] EQ min(d_guess_plaw)) OR (best_d_plaw[0] EQ max(d_guess_plaw))) THEN $
      print,'WARNING: Maybe didn''t explore enough parameter space for power-law delta!!'
ENDIF



IF keyword_set(fitbias) THEN BEGIN
   ;MAD read in model data
   readcol,dmfile,dmbin,dmw,format='F,D'
   dmbin=dmbin*60.
   inscale=where((dmbin GT minscale*60.) AND (dmbin LT maxscale*60.))
   dmbin=dmbin[inscale]
   dmw=dmw[inscale]

   ;MAD Generate array of bias values
   b_guess=(findgen(500000.)*(5./500000.)+1.)
   ;MAD Initialize arrays of Chi^2 values to fill; invert covariance matrix
   chisq=dblarr(n_elements(b_guess))
   C_inv=invert(C,/double)
   C_inv_use=C_inv[*,okfit]
   C_inv_use=C_inv_use[okfit,*]

   ;;MAD Loop over bias values to get chi-sq using just variance...
   ;print,'fit_bias - Building bias chi^2 array...'
   ;FOR i=0L,n_elements(b_guess)-1 DO BEGIN
   ; FOR k=0L,n_elements(bin[okfit])-1 DO BEGIN
   ;  val=dmw[okfit[k]]*(b_guess[i]^2.)
   ;  temp=(wtheta[okfit[k]]-val) * C_inv_use[k,k] * (wtheta[okfit[k]]-val)
   ;  chisq[i]=chisq[i]+temp
   ; ENDFOR
   ;ENDFOR

   ;MAD Loop over bias values to get chi-sq using covariance...
   print,'fit_bias - Building bias chi^2 array...'
   FOR i=0L,n_elements(b_guess)-1 DO BEGIN
      val=dmw[okfit]*(b_guess[i]^2.)
      temp=(wtheta[okfit]-val) # C_inv_use # (wtheta[okfit]-val)
      chisq[i]=chisq[i]+temp
   ENDFOR

   ;MAD Plot the chi^2 values vs A
   ;plot,b_guess,chisq,linestyle=1,xra=[0.8,7],yra=[0,500],xtit='bias',ytit='chi^2'
   ;stop

   ;MAD Find where the fixed-slope chi^2 values is minumum, get best bias value
   print,'fit_bias - finding min chi^2 and best fit...'
   minchi=where(chisq EQ min(chisq))
   best_b=b_guess[minchi]
   y=dmw[okfit]*(best_b[0]^2.)

   ;MAD find errors on best bias, using delta chi^2 = 1
   chi1=chisq[0:minchi]-min(chisq)
   chi2=chisq[minchi:n_elements(chisq)-1]-min(chisq)
   b_vals1=b_guess[0:minchi]
   b_vals2=b_guess[minchi:n_elements(b_guess)-1]
   xx=closest(chi1,1.)
   b_err1=b_vals1[xx]
   yy=closest(chi2,1.)
   b_err2=b_vals2[yy]
   b_err=(abs(best_b[0]-b_err1)+abs(best_b[0]-b_err2))*(0.5)
   dchi_b=(chi1[xx]+chi2[yy])*0.5

   ;MAD Plot results (to hardcopy of keyword set)
   IF (keyword_set(biasplot)) THEN PS_start,filename=biasplot
   xtit=textoidl('\theta (deg)')
   ytit=textoidl('\omega_{\theta}')
   temp=textoidl('\delta')
   nice_plot,0.002,5.0,0.0005,10.,xtit=xtit,ytit=ytit,/xlog,/ylog
   oplot,bin/60.,wtheta,psym=-8,color=cgcolor('blue')
   oploterror,bin/60.,wtheta,errors,psym=3,color=cgcolor('blue')
   oplot,dmbin[okfit]/60.,y,linestyle=1
   IF (keyword_set(biasplot)) THEN PS_end,/png

   bias=best_b[0]
   biaserr=b_err

   ;MAD Print results to screen
   print,' '

   IF (fitflag EQ 1) THEN print,'***************************************************************************'
   IF (fitflag EQ 1) THEN print,'***SOME POINTS WERE LEFT OUT OF FIT DUE TO NEGATIVE W_THETA IN JACKKNIFE***'
   IF (fitflag EQ 1) THEN print,'***************************************************************************'

   print,'Bias fit errors are based on a delta chi^2 of ',strtrim(dchi_b,2)
   print,' '
   print,'Fits using covariance ('+strtrim(minscale,2)+' - '+strtrim(maxscale,2)+' degrees):'
   print,'Bias = ',strtrim(best_b[0],2),' +/- ',strtrim(b_err,2)
   print,' '

   ;MAD If chi^2 minimization went all the way to the edge of the input
   ;grids, give a warning
   IF ((best_b[0] EQ min(b_guess)) OR (best_b[0] EQ max(b_guess))) THEN $
      print,'WARNING: Maybe didn''t explore enough parameter space for bias!!'
ENDIF


;MAD Get finish, elapsed time
et=systime(1)
print,'fit_ang_cluster - Elapsed time= ',strtrim((et-st)/60,2),' min'

return
END
