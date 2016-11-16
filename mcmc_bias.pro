;+
;  NAME:
;    mcmc_bias
;
;  PURPOSE:
;    Generate posterior probability distribution of bias values from
;    angular clustering results.  Uses Metropolis-Hastings MCMC, with
;    a Gaussian proposal mechanism.  Also assumes Gaussian prior on bias.
;
;  USE:
;    mcmc_bias,theta,omega_theta,omega_error,mod_theta,mod_omega,bias,$
;              [optional inputs]
;
;  INPUT:
;    theta - angular size bins (degrees, for plotting)
;    omega - angular clustering power at each theta
;    error - errors on omega
;    modtheta - angular bins of model (should be same as theta)
;    modomega - model DM clustering
;
;  OPTIONAL INPUT:
;    covariance - file name or array of covariance matrix.  If not
;                 supplied, uses just variance terms (supplied errors)
;    minscale - minimum scale of interest (degrees)
;    maxscale - maximim scale of interest (degrees)
;    b_prior_mean - mean of the Gaussian prior on bias (default = 2.5)
;    b_prior_sd - standard deviation of Gaussian prior on bias
;                 (default = 2.)
;    prop_width - standard deviation of distribution from which
;                 proposed bias steps are drawn.  Should tune until
;                 accepted fraction is ~0.5
;    samples - number of MCMC chain steps (default = 10000)
;    propseed - seed for random number generator for proposed bias
;    acceptseed - seed for random number generator to determine
;                 acceptence of proposed bias
;    figname - string name for an output figure of results
;    trim - fraction of samples to trim off the start of the chain
;           (default 0.05)
;    
;  OUTPUT:
;    bias - the posterior distribution of bias values
;    f_accept - set to variable to return accepted fraction of
;               proposal values.  Should be ~0.5, tune with prop_width keyword
;
;  NOTES:
;
;  HISTORY:
;    9-12-16 - Written - MAD (Dartmouth)
;-
PRO mcmc_bias,theta,omega,error,modtheta,modomega,bias,$
              covariance=covariance,$
              b_prior_mean=b_prior_mean,b_prior_sd=b_prior_sd,$
              prop_width=prop_width,samples=samples,$
              propseed=propseed,acceptseed=acceptseed,$
              minscale=minscale,maxscale=maxscale,f_accept=f_accept,$
              figname=figname,trim=trim

  ;MAD Start timer
  st=timer()
  
  ;MAD Set defaults
  IF ~keyword_set(b_prior_mean) THEN b_prior_mean=25
  IF ~keyword_set(b_prior_sd) THEN b_prior_sd=1.5
  IF ~keyword_set(samples) THEN samples=10000.
  IF ~keyword_set(propseed) THEN propseed=154
  IF ~keyword_set(acceptseed) THEN acceptseed=918
  IF ~keyword_set(prop_width) THEN prop_width=0.5
  IF ~keyword_set(trim) THEN trim=0.05
  
  ;MAD Read in covariance if given file
  IF keyword_set(covariance) THEN BEGIN
     IF (size(covariance,/type) EQ 7) THEN C=read_matrix(covariance) ELSE $
        C=covariance
  ENDIF
  
  th=theta
  w=omega
  modth=modtheta
  modw=modomega
  err=error

  ;MAD Limit to scales of interest (if set)
  IF keyword_set(maxscale) THEN BEGIN
     use=where(th LE maxscale)
     th=th[use]
     w=w[use]
     err=err[use]
     IF keyword_set(covariance) THEN BEGIN
        C=C[use,*]
        C=C[*,use]
     ENDIF
     use=where(modth LE maxscale)
     modth=modth[use]
     modw=modw[use]
  ENDIF
  IF keyword_set(minscale) THEN BEGIN
     use=where(th GE minscale)
     th=th[use]
     w=w[use]
     err=err[use]
     IF keyword_set(covariance) THEN BEGIN
        C=C[use,*]
        C=C[*,use]
     ENDIF
     use=where(modth GE minscale)
     modth=modth[use]
     modw=modw[use]
  ENDIF

  ;MAD Starting b, initialize n_accept
  b_current=b_prior_mean
  n_accept=0.

  ;MAD Start MCMC chain (Metropolis-Hastings chain)
  print,'MCMC_BIAS - Starting MCMC chain...'
  FOR i=0L,samples-1 DO BEGIN
     counter,i,samples
     ;MAD Get proposed new bias
     b_proposal=(randomn(propseed)*prop_width)+b_current

     ;MAD Calculate likelihoods
     IF keyword_set(covariance) THEN BEGIN
        term1=1/sqrt(DETERM(C,/double))
        ;MAD Current L
        term2=EXP((-0.5)*((w-((b_current^2.)*modw)) # INVERT(C,/double) # (w-((b_current^2.)*modw))))
        like_current=(1./sqrt(2.*!dpi))*product(term1*term2)

        ;MAD Proposed L
        term2=EXP((-0.5)*((w-((b_proposal^2.)*modw)) # INVERT(C,/double) # (w-((b_proposal^2.)*modw))))
        like_proposal=(1./sqrt(2.*!dpi))*product(term1*term2)
     ENDIF ELSE BEGIN
        ;MAD Current L
        term1=1./SQRT(2.*!dpi*err^2.)
        term2=EXP(-((w-((b_current^2.)*modw))^2.)/(2.*err^2.))
        like_current=product(term1*term2)

        ;MAD Proposed L
        term1=1./SQRT(2.*!dpi*err^2.)
        term2=EXP(-((w-((b_proposal^2.)*modw))^2.)/(2.*err^2.))
        like_proposal=product(term1*term2)
     ENDELSE
        
     ;MAD Get prior probability for current and proposed values
     prior_current=norm_pdf(b_current,mu=b_prior_mean,sd=b_prior_sd)
     prior_proposal=norm_pdf(b_proposal,mu=b_prior_mean,sd=b_prior_sd)

     ;MAD posterior probability (likelihood * prior)
     p_current=like_current*prior_current
     p_proposal=like_proposal*prior_proposal

     ;MAD Decide whether to accept or reject proposed bias
     IF ((randomu(acceptseed) LT (p_proposal/p_current)) AND (b_proposal GT 0)) THEN BEGIN
        n_accept=n_accept+1
        b_current=b_proposal
     ENDIF

     ;MAD Save either current or proposed bias
     IF (n_elements(bias) EQ 0) THEN bias=b_current ELSE $
        bias=[bias,b_current]
  ENDFOR
  
  bias=bias[(n_elements(bias)*trim):(n_elements(bias)-1)]

  ;MAD determine accepted fraction
  f_accept=n_accept/samples
  bias_err=conf_int(bias,median(bias))

  ;MAD Make some output plots, if file specified
  IF keyword_set(figname) THEN BEGIN
     PS_start,filename=figname,xsize=6,ysize=10,/encapsul
     ;MAD For some reason, if samples > 13333, without
     ;MAD setting font=1, everything becomes greek?
     !P.Font=1
     !P.Multi=[0,1,3]
     plot,findgen(samples*(1.-trim)),bias,psym=3,xtit='Sample',ytit='bias',charsize=2.,thick=2,$
          yra=[min(bias)-0.5,max(bias)+1.],ystyle=0,xra=[0.,max(samples)*(1.-trim)],xsty=1
     legend,['Accepted fraction = '+strtrim(n_accept/samples,2)],/top,/right,box=0,charsize=1.5

     plothist,bias,xvals,yvals,bin=0.05,charsize=2,$
              thick=3,xtit='b',ytit='N'

     legend,['Mean = '+strtrim(mean(bias),2),'Median = '+strtrim(median(bias),2),$
             'Sigma = '+strtrim(stddev(bias),2),'67% conf = '+strtrim(bias_err,2)]$
            ,/top,/left,box=0,charsize=1.5
     oplot,[mean(bias),mean(bias)],[0,1e6],linestyle=1,thick=1
     oplot,[median(bias),median(bias)],[0,1e6],linestyle=3,thick=1

     circsym,/fill
     plot,th,w,psym=8,xtit=textoidl('\theta [deg]'),ytit=textoidl('\omega'),$
          xra=[0.009,3],xsty=1,yra=[1e-4,2],ysty=1,/xlog,/ylog,charsize=2.,symsize=1.5
     oploterror,th,w,err,psym=3,thick=1.5
     oplot,modth,modw,linestyle=1,thick=1.5
     oplot,modth,modw*(median(bias)^2.),linestyle=2,thick=1.5

     PS_end
  ENDIF

  print,'MCMC_BIAS - accepted fraction is '+strtrim(f_accept,2)
  
  ;MAD End timer
  print,'MCMC_BIAS - '
  et=timer(st=st,/fin,unit='s')
  
  return
END
