;+
;  NAME:
;    ang_cluster
;
;  PURPOSE:
;    Calculate the angular auto/cross-correlation function using the Landy &
;    Szalay (1993) estimator.  Smallest possible bin edge is 0.94" (2e-4 deg).
;    Can also call several other procedures to jackknife errors, fit
;    power laws, fit bias, etc.
;
;    !!!!IMPORTANT!!!! Code automatically outputs some files with
;    standard names.  It will look for these and then read them in if
;    they exist to skip steps.  If you change things (like your
;    sample) make sure you delete these files!
;    !!!!!!!!!!!!!!!!!
;
;  USE:
;    ang_cluster,data,rand,theta,wtheta[,data2=data2,errs=errs,maxscale=maxscale, $
;                minscale=minscale,outfile=''outfile.txt'',bins=bins, $
;                /fitplaws,/jackknife,/fitbias,dmfile=''model_dm.txt'']
;
;  INPUTS:
;    data - Structure with your data.  Should have tags ra and dec.
;    rand - Structure with random sample for normalization.  Should
;           have tags ra and dec
;
;  OPTIONAL INPUT:
;    data2_in - Structure with second data set for cross-correlation.
;            Always assumed that randoms follow regular data, not
;            data2, and DR values calculated as such.
;
;    rand2_in - structure with randoms mimicking second data set. If
;               not set, uses same random for both in cross-correlation.
;
;    maxscale - the maximum scale of interest (in degrees).  Defaults to
;            2.
;    minscale - the minimum scale of interest (in degrees). Only
;               applies to fitting.  Defaults to 0.0027 (~10")
;    outfile - set to a string name of output file.  Will contain
;              scale, w_theta, and error (if jackknife set)
;    bins - options for binning include 3, 4 or 5 (refers to bins/dex)
;           default is 5.  Can also set to 0 to make two large bins.
;    dmfile - must be set if /bias is set.  String name of file with
;             DM model, two columns: scale, w_theta)
;    n - number of pixels per side to use if jackknifing errors
;        (default = 9).
;    frac - fraction of rands to use when splitting up regions for jackknife
;
;  KEYWORDS:
;    fitplaws - Set to fit powerlaws (2 parameter and fixed slope=-1) to
;          result
;    jackknife - set to break up data into regions with equal
;                areas, and measure the clustering neglecting one
;                region at a time to estimate errors
;    fitbias - set to fit a DM model to the data (requires dmcluster to
;           be set)
;
;  OUTPUT:
;    theta - centers of scale bins
;    wtheta - autocorrelation values
;    errs - error bars (optional, must be set if /jackknife set)
;    RR.txt - If the random-random counts are being counted for the first
;             time, will output them for each bin to the file RR.txt.
;             These can be read in next time to save computation time- 
;             BE SURE TO erase RR.txt if you change the random catalog at
;             all!!
;    DD_DR.txt - same as above, for DD and DR counts.  Make sure to delete
;                if you change the data sample at all!
;    DR2.txt - same as above for second data set if doing a cross correlation.
;    rr_reg.txt - If jackknife keyword is set, writes a file containing the number
;                 of RR counts in each bin and pixel.  Has n_rows=n_bins for clustering
;                 measurement, and n_columns=n_pix (currently only works
;                 for 16 pixels)
;    dd_reg.txt - same as above, for DD counts.
;    dr_reg.txt - same as above, for DR counts.
;
;  HISTORY:
;    2013 - Written - MAD (UWyo) 
;    12-1-13 - First combined version - MAD (UWyo)
;    3-27-15 - Reworked and cleaned - MAD (UWyo)
;    7-17-15 - Added data2 keyword for cross-correlations - MAD (UWyo)
;    9-17-15 - Fixed some cross correlation bugs - MAD (Dartmouth)
;     9-5-16 - Generalized to arbitrary number or jackknife pixels -
;              MAD (Dartmouth)
;    3-15-17 - Generalized cross-corr, added second random catalog
;              option - MAD (Dartmouth)
;-
PRO ang_cluster,data_in,rand_in,theta,w_theta,data2_in=data2_in,rand2_in=rand2_in,$
                errs=errs,maxscale=maxscale,minscale=minscale,outfile=outfile,bins=bins, $
                fitplaws=fitplaws,jackknife=jackknife,n=n,frac=frac,fitbias=fitbias,dmfile=dmfile

IF (n_elements(data_in) EQ 0) THEN message,'Syntax - ang_cluster,data,rand,theta,w_theta[,errs=errs,maxscale=maxscale,minscale=minscale,outfile=''outfile.txt'',bins=bins,/fitplaws,/jackknife,n=n,frac=frac,/fitbias,dmfile=''model_dm.txt'']'

;MAD Get start time
st=systime(1)

;MAD Set default min/max scale
IF ~keyword_set(minscale) THEN minscale=0.0027
IF ~keyword_set(maxscale) THEN maxscale=2.

;MAD Set default number of pixels per side if doing jackknife errors
IF (~keyword_set(n) AND keyword_set(jackkinfe)) THEN n=9
IF (~keyword_set(frac) AND keyword_set(jackknife)) THEN frac=0.5

;MAD If jackknife keyword set (and they haven't been already), split data into 
;16 pixels for jackknife calculations
IF (keyword_set(jackknife)) THEN BEGIN
   tagflag_data=tag_exist(data_in,'reg',/quiet)
   tagflag_rand=tag_exist(rand_in,'reg',/quiet)
   IF (keyword_set(data2_in) AND keyword_set(rand2_in)) THEN BEGIN
      tagflag_data2=tag_exist(data2_in,'reg',/quiet)
      tagflag_rand2=tag_exist(rand2_in,'reg',/quiet)
      IF ((tagflag_data EQ 0) OR (tagflag_rand EQ 0) OR (tagflag_data2 EQ 0) OR (tagflag_rand2 EQ 0)) THEN BEGIN
         split_regions_gen,data_in,rand_in,data,rand,n=n, frac=frac, $
                           data_fileout='data_reg.fits',rand_fileout='rand_reg.fits', $
                           data2_in=data2_in,data2_fileout='data2_reg.fits',$
                           rand2_in=rand2_in,rand2_fileout='rand2_reg.fits',$
                           split_file='jack_splits.txt',/figures
         data2=mrdfits('data2_reg.fits',1)
         rand2=mrdfits('rand2_reg.fits',1)
      ENDIF ELSE BEGIN
         data=mrdfits('data_reg.fits',1)
         rand=mrdfits('rand_reg.fits',1)
         data2=mrdfits('data2_reg.fits',1)
         rand2=mrdfits('rand2_reg.fits',1)
      ENDELSE
   ENDIF ELSE IF (keyword_set(data2_in) AND ~keyword_set(rand2_in)) THEN BEGIN
      tagflag_data2=tag_exist(data2_in,'reg',/quiet)
      IF ((tagflag_data EQ 0) OR (tagflag_rand EQ 0) OR (tagflag_data2 EQ 0)) THEN BEGIN
         split_regions_gen,data_in,rand_in,data,rand,n=n, frac=frac, $
                           data_fileout='data_reg.fits',rand_fileout='rand_reg.fits', $
                           data2_in=data2_in,data2_fileout='data2_reg.fits',$
                           split_file='jack_splits.txt',/figures
         data2=mrdfits('data2_reg.fits',1)
      ENDIF ELSE BEGIN
         data=mrdfits('data_reg.fits',1)
         rand=mrdfits('rand_reg.fits',1)
         data2=mrdfits('data2_reg.fits',1)
      ENDELSE
   ENDIF ELSE BEGIN
      IF ((tagflag_data EQ 0) OR (tagflag_rand EQ 0)) THEN BEGIN
         split_regions_gen,data_in,rand_in,data,rand, n=n, frac=frac, $
                           data_fileout='data_reg.fits',rand_fileout='rand_reg.fits',$
                           split_file='jack_splits.txt',/figures
      ENDIF ELSE BEGIN
         data=mrdfits('data_reg.fits',1)
         rand=mrdfits('rand_reg.fits',1)
      ENDELSE
   ENDELSE
ENDIF ELSE BEGIN
   data=data_in
   rand=rand_in
   IF keyword_set(data2_in) THEN data2=data2_in
   IF keyword_set(rand2_in) THEN rand2=rand2_in
ENDELSE

;MAD Get number in each sample in double format
n_data=double(n_elements(data))
n_rand=double(n_elements(rand))
IF keyword_set(data2_in) THEN n_data2=double(n_elements(data2))
IF keyword_set(rand2_in) THEN n_rand2=double(n_elements(rand2))

;MAD Set binning (5dex is the default), make bin edges, centers
IF (n_elements(bins) EQ 0) THEN bins=5
IF (bins NE 0) THEN BEGIN
   bininfo=define_bins(bins,maxscale)
   bin_edge=bininfo.edges
   bin_cent=bininfo.cents
ENDIF

;MAD Define big bins, if needed
IF (bins EQ 0) THEN BEGIN
;   bin_edge=[0.2,6.,maxscale*60.]
;   bin_cent=[3.1,((6.+(maxscale*60.))/2.),(maxscale*60.)+30.]
   bin_edge=[1,10,60.]
   bin_cent=[5.5,35,90]
ENDIF

;MAD If including outputs information of RR for jackknife errors,
;intialize array to fill
IF (keyword_set(jackknife)) THEN BEGIN
   rr_regions=dblarr(max(rand.reg),n_elements(bin_cent))
   dr_regions=dblarr(max(rand.reg),n_elements(bin_cent))
   IF (keyword_set(data2_in) AND keyword_set(rand2_in)) THEN $
      dr2_regions=dblarr(max(rand.reg),n_elements(bin_cent))
   dd_regions=dblarr(max(rand.reg),n_elements(bin_cent))
ENDIF

;MAD Look for file of RR counts done previously
rr_file=file_search('RR.txt')
rr_reg_file=file_search('rr_reg.txt')
;MAD If RR counts haven't been done, do them
;MAD (This all looks very messy, but most is just careful bookkeeping
;MAD to 1: Break it up into steps to speed things up, 2: to make it
;MAD easier to rerun if needed, and 3: to make the jackknife error
;MAD SUPER fast!!)
IF (rr_file EQ '') THEN BEGIN
   ;MAD Loop over random sample in pieces to calculate RR
   print,'Ang_cluster - starting loop for RR counts...'
   h_rr=0.D
   k=long(0)
   step=10000.
   WHILE (k LT n_elements(rand)) DO BEGIN
      temprand=rand[k:k+step-1]
      IF ~keyword_set(rand2_in) THEN $
         spherematch,rand.ra,rand.dec,temprand.ra,temprand.dec,maxscale,m2_rr,m1_rr,sep_rr,maxmatch=0 ELSE $
            spherematch,rand2.ra,rand2.dec,temprand.ra,temprand.dec,maxscale,m2_rr,m1_rr,sep_rr,maxmatch=0
      sep_rr=sep_rr*60.
      xx=where(sep_rr GE min(bin_edge))
      ;Only continue if there are separations > the minimum bin limit
      IF (xx[0] NE -1) THEN BEGIN
         sep_rr=sep_rr[xx]
         m1_rr=m1_rr[xx]
         m2_rr=m2_rr[xx]
         rebinned=value_locate(bin_edge,sep_rr)
         temp_rr=double(histogram(rebinned,binsize=1,min=0,max=n_elements(bin_edge)-1))
         IF ~keyword_set(rand2_in) THEN $
            temp_rr=temp_rr*(1./(n_rand*n_rand)) ELSE $
               temp_rr=temp_rr*(1./(n_rand*n_rand2))
         h_rr=h_rr+temp_rr
         ;If error keyword set, find pixel each pair member lives in
         IF (keyword_set(jackknife)) THEN BEGIN
            FOR v=0L,n_elements(bin_edge)-1 DO BEGIN
               xx=where(rebinned EQ v)
               IF (xx[0] NE -1) THEN BEGIN
                  trand1=[temprand[m1_rr[xx]].reg]
                  IF ~keyword_set(rand2_in) THEN $
                     trand2=[rand[m2_rr[xx]].reg] ELSE $
                        trand2=[rand2[m2_rr[xx]].reg]
                  yy=where(trand1 EQ trand2)
                  zz=where(trand1 NE trand2)
                  IF ((zz[0] NE -1) AND (yy[0] NE -1)) THEN BEGIN
                     h1=double(histogram(trand1[zz],binsize=1,min=1,max=max(rand.reg)))
                     h2=double(histogram(trand2[zz],binsize=1,min=1,max=max(rand.reg)))
                     h3=double(histogram(trand1[yy],binsize=1,min=1,max=max(rand.reg)))
                     rr_regions[*,v]=rr_regions[*,v]+h1+h2+h3
                  ENDIF 
                  IF ((yy[0] NE -1) AND (zz[0] EQ -1)) THEN BEGIN
                     h3=double(histogram(trand1[yy],binsize=1,min=1,max=max(rand.reg)))
                     rr_regions[*,v]=rr_regions[*,v]+h3
                  ENDIF
                  IF ((zz[0] NE -1) AND (yy[0] EQ -1)) THEN BEGIN
                     h1=double(histogram(trand1[zz],binsize=1,min=1,max=max(rand.reg)))
                     h2=double(histogram(trand2[zz],binsize=1,min=1,max=max(rand.reg)))
                     rr_regions[*,v]=rr_regions[*,v]+h1+h2
                  ENDIF
               ENDIF
            ENDFOR
         ENDIF
      ENDIF
      print, format='("RR iteration done, random points ",i8," - ",i8," (",f8.4,"%)",a1,$)',$
             k,k+step,(k+step)*(1./n_elements(rand))*100.,string(13B)
      IF (k+step GE n_elements(rand)) THEN print,'ANG_CLUSTER - RR counts done...     '
      k=k+step
      IF (k+10000 LT n_elements(rand)) THEN BEGIN
         step=10000
      ENDIF ELSE BEGIN
         IF (k+1000 LT n_elements(rand)) THEN BEGIN
            step=1000
         ENDIF ELSE BEGIN
            IF (k+100 LT n_elements(rand)) THEN BEGIN
               step=100
            ENDIF ELSE BEGIN
               IF (k+10 LT n_elements(rand)) THEN BEGIN
                  step=10
               ENDIF ELSE BEGIN
                  step=1
               ENDELSE
            ENDELSE
         ENDELSE
      ENDELSE
   ENDWHILE

   ;MAD Write out RR counts text file
   openw,lun,'RR.txt',/get_lun
   FOR i=0L,n_elements(h_rr)-1 DO BEGIN
      printf,lun,h_rr[i],format='(D)'
   ENDFOR
   close,lun
   
;MAD If RR counts have already been done, read them in
ENDIF ELSE BEGIN
   print,'Ang_cluster - reading in previous RR counts file...'
   readcol,rr_file[0],h_rr,format='D'
ENDELSE

;MAD If doing errors, then write out RR file for those
IF (keyword_set(jackknife)) THEN BEGIN
   IF (rr_reg_file EQ '') THEN BEGIN
      print,'Ang_cluster - writing out RR counts per region...'
      openw,lun,'rr_reg.txt',/get_lun
      FOR w=0L,max(rand.reg)-1 DO BEGIN
         xx=where(finite(rr_regions[w,*]) EQ 0,cnt)
         IF (cnt NE 0) THEN rr_regions[w,xx] = -9999
         fmtstring='('
         FOR l=0L,n_elements(bin_cent)-2 DO BEGIN
            fmtstring=fmtstring+'E,1x,'
         ENDFOR
         fmtstring=fmtstring+'E)'
         printf,lun,rr_regions[w,*],format=fmtstring
      ENDFOR
      close,lun
   ENDIF
ENDIF

;MAD Look for file of DD/DR counts done previously
d_file=file_search('DD_DR.txt')
dd_reg_file=file_search('dd_reg.txt')
dr_reg_file=file_search('dr_reg.txt')
;MAD If DD/DR counts haven't been done, do them
IF (d_file EQ '') THEN BEGIN
   ;MAD Loop over chunks of data to find DD and DR
   print,'Ang_cluster - starting loop for DD/DR counts...'
   h_dd=0.D
   h_dr=0.D
   k=long(0)
   IF (n_elements(data) LT  10000.-1.) THEN step=n_elements(data)-1 ELSE step=10000.
   WHILE (k LT n_elements(data)) DO BEGIN
   ;Data-data counts
      tempdata=data[k:k+step-1]
      IF ~keyword_set(data2) THEN $
         spherematch,data.ra,data.dec,tempdata.ra,tempdata.dec,maxscale,m2_dd,m1_dd,sep_dd,maxmatch=0 ELSE $
            spherematch,data2.ra,data2.dec,tempdata.ra,tempdata.dec,maxscale,m2_dd,m1_dd,sep_dd,maxmatch=0
      sep_dd=sep_dd*60.
      xx=where(sep_dd GE min(bin_edge))
      IF (xx[0] NE -1) THEN BEGIN
         sep_dd=sep_dd[xx]
         m1_dd=m1_dd[xx]
         m2_dd=m2_dd[xx]
         rebinned=value_locate(bin_edge,sep_dd)
         temp_dd=double(histogram(rebinned,binsize=1,min=0,max=n_elements(bin_edge)-1))
         IF ~keyword_set(data2) THEN temp_dd=temp_dd*(1./(n_data*n_data)) ELSE $
            temp_dd=temp_dd*(1./(n_data*n_data2))
         h_dd=h_dd+temp_dd 
         ;If error keyword set, find pixel each pair member lives in
         IF (keyword_set(jackknife)) THEN BEGIN
            FOR v=0L,n_elements(bin_edge)-1 DO BEGIN
               xx=where(rebinned EQ v)
               IF (xx[0] NE -1) THEN BEGIN
                  tdata1=[tempdata[m1_dd[xx]].reg]
                  IF ~keyword_set(data2) THEN tdata2=[data[m2_dd[xx]].reg] ELSE tdata2=[data2[m2_dd[xx]].reg]
                  yy=where(tdata1 EQ tdata2)
                  zz=where(tdata1 NE tdata2)
                  IF ((zz[0] NE -1) AND (yy[0] NE -1)) THEN BEGIN
                     h1=double(histogram(tdata1[zz],binsize=1,min=1,max=max(data.reg)))
                     h2=double(histogram(tdata2[zz],binsize=1,min=1,max=max(data.reg)))
                     h3=double(histogram(tdata1[yy],binsize=1,min=1,max=max(data.reg)))
                     dd_regions[*,v]=dd_regions[*,v]+h1+h2+h3
                  ENDIF 
                  IF ((yy[0] NE -1) AND (zz[0] EQ -1)) THEN BEGIN
                     h3=double(histogram(tdata1[yy],binsize=1,min=1,max=max(rand.reg)))
                     dd_regions[*,v]=dd_regions[*,v]+h3
                  ENDIF
                  IF ((zz[0] NE -1) AND (yy[0] EQ -1)) THEN BEGIN
                     h1=double(histogram(tdata1[zz],binsize=1,min=1,max=max(rand.reg)))
                     h2=double(histogram(tdata2[zz],binsize=1,min=1,max=max(rand.reg)))
                     dd_regions[*,v]=dd_regions[*,v]+h1+h2
                  ENDIF
               ENDIF
            ENDFOR
         ENDIF
      ENDIF

      ;Data-random counts
      tempdata=data[k:k+step-1]
      IF ~keyword_set(rand2_in) THEN $
         spherematch,rand.ra,rand.dec,tempdata.ra,tempdata.dec,maxscale,m2_dr,m1_dr,sep_dr,maxmatch=0 ELSE $
            spherematch,rand2.ra,rand2.dec,tempdata.ra,tempdata.dec,maxscale,m2_dr,m1_dr,sep_dr,maxmatch=0
      sep_dr=sep_dr*60.
      xx=where(sep_dr GE min(bin_edge))
      IF (xx[0] NE -1) THEN BEGIN
         sep_dr=sep_dr[xx]
         m1_dr=m1_dr[xx]
         m2_dr=m2_dr[xx]
         rebinned=value_locate(bin_edge,sep_dr)
         temp_dr=double(histogram(rebinned,binsize=1,min=0,max=n_elements(bin_edge)-1))
         IF ~keyword_set(rand2_in) THEN $
            temp_dr=temp_dr*(1./(n_data*n_rand)) ELSE $
               temp_dr=temp_dr*(1./(n_data*n_rand2))
         h_dr=h_dr+temp_dr
         ;If error keyword set, find pixel each pair member lives in
         IF (keyword_set(jackknife)) THEN BEGIN
            FOR v=0L,n_elements(bin_edge)-1 DO BEGIN
               xx=where(rebinned EQ v)
               IF (xx[0] NE -1) THEN BEGIN
                  tdata=[tempdata[m1_dr[xx]].reg]
                  IF ~keyword_set(rand2_in) THEN $
                     trand=[rand[m2_dr[xx]].reg] ELSE $
                        trand=[rand2[m2_dr[xx]].reg]
                  yy=where(tdata EQ trand)
                  zz=where(tdata NE trand)
                  IF ((zz[0] NE -1) AND (yy[0] NE -1)) THEN BEGIN
                     h1=double(histogram(tdata[zz],binsize=1,min=1,max=max(data.reg)))
                     h2=double(histogram(trand[zz],binsize=1,min=1,max=max(data.reg)))
                     h3=double(histogram(tdata[yy],binsize=1,min=1,max=max(data.reg)))
                     dr_regions[*,v]=dr_regions[*,v]+h1+h2+h3
                  ENDIF
                  IF ((yy[0] NE -1) AND (zz[0] EQ -1)) THEN BEGIN
                     h3=double(histogram(tdata[yy],binsize=1,min=1,max=max(rand.reg)))
                     dr_regions[*,v]=dr_regions[*,v]+h3
                  ENDIF
                  IF ((zz[0] NE -1) AND (yy[0] EQ -1)) THEN BEGIN
                     h1=double(histogram(tdata[zz],binsize=1,min=1,max=max(rand.reg)))
                     h2=double(histogram(trand[zz],binsize=1,min=1,max=max(rand.reg)))
                     dr_regions[*,v]=dr_regions[*,v]+h1+h2
                  ENDIF
               ENDIF
            ENDFOR
         ENDIF
      ENDIF
      print, format='("DD/DR iteration done, data points ",i8," - ",i8," (",f8.4,"%)",a1,$)',$
                k,k+step,(k+step)*(1./n_elements(data))*100.,string(13B)
      IF (k+step GE n_elements(data)) THEN print,'ANG_CLUSTER - DD/DR counts done...   '
      k=k+step
      IF (k+10000 LT n_elements(data)) THEN BEGIN
         step=10000
      ENDIF ELSE BEGIN
         IF (k+1000 LT n_elements(data)) THEN BEGIN
            step=1000
         ENDIF ELSE BEGIN
            IF (k+100 LT n_elements(data)) THEN BEGIN
               step=100
            ENDIF ELSE BEGIN
               IF (k+10 LT n_elements(data)) THEN BEGIN
                  step=10
               ENDIF ELSE BEGIN
                  step=1
               ENDELSE
            ENDELSE
         ENDELSE
      ENDELSE
   ENDWHILE
   ;MAD Write out DD counts text file
   openw,lun,'DD_DR.txt',/get_lun
   FOR i=0L,n_elements(h_dd)-1 DO BEGIN
      printf,lun,h_dd[i],h_dr[i],format='(D,D)'
   ENDFOR
   close,lun
;MAD If DD/DR counts have already been done, read them in
ENDIF ELSE BEGIN
   print,'Ang_cluster - reading in previous DD/DR counts...'
   readcol,d_file,h_dd,h_dr,format='D,D'
ENDELSE

;MAD If doing errors, then write out DD file for those
IF (keyword_set(jackknife)) THEN BEGIN
   IF (dd_reg_file EQ '') THEN BEGIN
      openw,lun,'dd_reg.txt',/get_lun
      FOR w=0L,max(rand.reg)-1 DO BEGIN
         xx=where(finite(dd_regions[w,*]) EQ 0,cnt)
         IF (cnt NE 0) THEN dd_regions[w,xx] = -9999
         fmtstring='('
         FOR l=0L,n_elements(bin_cent)-2 DO BEGIN
            fmtstring=fmtstring+'E,1x,'
         ENDFOR
         fmtstring=fmtstring+'E)'
         printf,lun,dd_regions[w,*],format=fmtstring
      ENDFOR      
      close,lun
   ENDIF
ENDIF

;MAD If doing errors, then write out DR file for those
IF (keyword_set(jackknife)) THEN BEGIN
   IF (dr_reg_file EQ '') THEN BEGIN
      print,'Ang_cluster - writing out DD/DR counts per region...'
      openw,lun,'dr_reg.txt',/get_lun
      FOR w=0L,max(rand.reg)-1 DO BEGIN
         xx=where(finite(dr_regions[w,*]) EQ 0,cnt)
         IF (cnt NE 0) THEN dr_regions[w,xx] = -9999
         fmtstring='('
         FOR l=0L,n_elements(bin_cent)-2 DO BEGIN
            fmtstring=fmtstring+'E,1x,'
         ENDFOR
         fmtstring=fmtstring+'E)'
         printf,lun,dr_regions[w,*],format=fmtstring
      ENDFOR
      close,lun
   ENDIF
ENDIF




IF (keyword_set(data2_in) AND keyword_set(rand2_in)) THEN BEGIN
   ;MAD Look for file of D2R counts done previously
   d2_file=file_search('DR2.txt')
   dr2_reg_file=file_search('dr2_reg.txt')
   ;MAD If DR2 counts haven't been done, do them
   IF (d2_file EQ '') THEN BEGIN
   ;MAD Loop over chunks of data to find DD and DR
      print,'Ang_cluster - starting loop for D2R counts...'
      h_dr2=0.D
      k=long(0)
      IF (n_elements(data2) LT  10000.-1.) THEN step=n_elements(data2)-1 ELSE step=10000.
      WHILE (k LT n_elements(data2)) DO BEGIN
         ;Data-random counts
         tempdata=data2[k:k+step-1]
         spherematch,rand.ra,rand.dec,tempdata.ra,tempdata.dec,maxscale,m2_dr,m1_dr,sep_dr,maxmatch=0
         sep_dr=sep_dr*60.
         xx=where(sep_dr GE min(bin_edge))
         IF (xx[0] NE -1) THEN BEGIN
            sep_dr=sep_dr[xx]
            m1_dr=m1_dr[xx]
            m2_dr=m2_dr[xx]
            rebinned=value_locate(bin_edge,sep_dr)
            temp_dr=double(histogram(rebinned,binsize=1,min=0,max=n_elements(bin_edge)-1))
            temp_dr=temp_dr*(1./(n_data2*n_rand))
            h_dr2=h_dr2+temp_dr
            ;If error keyword set, find pixel each pair member lives in
            IF (keyword_set(jackknife)) THEN BEGIN
               FOR v=0L,n_elements(bin_edge)-1 DO BEGIN
                  xx=where(rebinned EQ v)
                  IF (xx[0] NE -1) THEN BEGIN
                     tdata=[tempdata[m1_dr[xx]].reg]
                     trand=[rand[m2_dr[xx]].reg]
                     yy=where(tdata EQ trand)
                     zz=where(tdata NE trand)
                     IF ((zz[0] NE -1) AND (yy[0] NE -1)) THEN BEGIN
                        h1=double(histogram(tdata[zz],binsize=1,min=1,max=max(data2.reg)))
                        h2=double(histogram(trand[zz],binsize=1,min=1,max=max(data2.reg)))
                        h3=double(histogram(tdata[yy],binsize=1,min=1,max=max(data2.reg)))
                        dr2_regions[*,v]=dr2_regions[*,v]+h1+h2+h3
                     ENDIF
                     IF ((yy[0] NE -1) AND (zz[0] EQ -1)) THEN BEGIN
                        h3=double(histogram(tdata[yy],binsize=1,min=1,max=max(rand2.reg)))
                        dr2_regions[*,v]=dr2_regions[*,v]+h3
                     ENDIF
                     IF ((zz[0] NE -1) AND (yy[0] EQ -1)) THEN BEGIN
                        h1=double(histogram(tdata[zz],binsize=1,min=1,max=max(rand2.reg)))
                        h2=double(histogram(trand[zz],binsize=1,min=1,max=max(rand2.reg)))
                        dr2_regions[*,v]=dr2_regions[*,v]+h1+h2
                     ENDIF
                  ENDIF
               ENDFOR
            ENDIF
         ENDIF
         print, format='("DR2 iteration done, data points ",i8," - ",i8," (",f8.4,"%)",a1,$)',$
                k,k+step,(k+step)*(1./n_elements(data2))*100.,string(13B)
         IF (k+step GE n_elements(data2)) THEN print,'ANG_CLUSTER - DD/DR counts done...   '
         k=k+step
         IF (k+10000 LT n_elements(data2)) THEN BEGIN
            step=10000
         ENDIF ELSE BEGIN
            IF (k+1000 LT n_elements(data2)) THEN BEGIN
               step=1000
            ENDIF ELSE BEGIN
               IF (k+100 LT n_elements(data2)) THEN BEGIN
                  step=100
               ENDIF ELSE BEGIN
                  IF (k+10 LT n_elements(data2)) THEN BEGIN
                     step=10
                  ENDIF ELSE BEGIN
                     step=1
                  ENDELSE
               ENDELSE
            ENDELSE
         ENDELSE
      ENDWHILE
      ;MAD Write out DR2 counts text file
      openw,lun,'DR2.txt',/get_lun
      FOR i=0L,n_elements(h_dr2)-1 DO BEGIN
         printf,lun,h_dr2[i],format='(D,D)'
      ENDFOR
      close,lun
   ;MAD If DR2 counts have already been done, read them in
   ENDIF ELSE BEGIN
      print,'Ang_cluster - reading in previous DR2 counts...'
      readcol,d2_file,h_dr2,format='D'
   ENDELSE

   ;MAD If doing errors, then write out DR2 file for those
   IF (keyword_set(jackknife)) THEN BEGIN
      IF (dr2_reg_file EQ '') THEN BEGIN
         openw,lun,'dr2_reg.txt',/get_lun
         FOR w=0L,max(rand2.reg)-1 DO BEGIN
            xx=where(finite(dr2_regions[w,*]) EQ 0,cnt)
            IF (cnt NE 0) THEN dr2_regions[w,xx] = -9999
            fmtstring='('
            FOR l=0L,n_elements(bin_cent)-2 DO BEGIN
               fmtstring=fmtstring+'E,1x,'
            ENDFOR
            fmtstring=fmtstring+'E)'
            printf,lun,dr2_regions[w,*],format=fmtstring
         ENDFOR      
         close,lun
      ENDIF
   ENDIF
ENDIF





print,'Ang_cluster - calculating auto/cross-correlation, trimming bad data...'
;MAD Calculate auto/cross-correlation (Landy & Szalay 1993)
IF (~keyword_set(data2_in) AND ~keyword_set(rand2_in)) THEN $
   w_theta=(1./h_rr)*(h_dd-(2.*h_dr)+h_rr) ELSE $
      w_theta=(1./h_rr)*(h_dd-h_dr-h_dr2+h_rr)
theta=bin_cent

;If the scale of interest is larger than the maximum edge of
;the bins, the last value is nonsense
IF (maxscale*60. GT max(bin_edge)) THEN BEGIN
   theta=bin_cent[0:n_elements(w_theta)-2]
   w_theta=w_theta[0:n_elements(w_theta)-2]
ENDIF

;MAD Limit only to values that aren't INF or NaN
xx=where(finite(w_theta) EQ 0,cnt)
IF (cnt NE 0) THEN BEGIN
   w_theta[xx]=-9999
   theta[xx]=-9999
ENDIF

;MAD Print main results to file (if errors not being done)
IF ~keyword_set(jackknife) THEN BEGIN
   print,'Ang_cluster - writing out auto/cross-correlation results...'
   openw,lun,outfile,/get_lun
   FOR i=0L,n_elements(w_theta)-1 DO BEGIN
      printf,lun,theta[i],w_theta[i]
   ENDFOR
   close,lun
ENDIF


;MAD Do jackknife error calculation
IF (keyword_set(jackknife)) THEN BEGIN
   print,'Ang_cluster - calling ang_jackknife.pro to do errors...'
   IF ~keyword_set(data2) THEN data2=0
   IF ~keyword_set(rand2) THEN rand2=0
;   IF ~keyword_set(data2) THEN $
;      ang_jackknife,data,rand,theta,w_theta,errs,maxscale=maxscale,outfile=outfile,bins=bins ELSE $
   ang_jackknife,data,rand,theta,w_theta,errs,data2=data2,rand2=rand2,maxscale=maxscale,outfile=outfile,bins=bins
ENDIF

;MAD Fit the results
IF (keyword_set(fitplaws)) THEN BEGIN
   print,'Ang_cluster - calling fit_ang_cluster to do fits, make plots...'
   IF (keyword_set(fitbias)) THEN BEGIN
      fit_ang_cluster,theta,w_theta,errs,bias,biaserr,minscale=minscale,maxscale=maxscale,dmfile=dmfile,$
                      /fitplaws,/fitbias,plawplot='ang_cluster_plaw.png',biasplot='ang_cluster_bias.png'
   ENDIF ELSE BEGIN
      fit_ang_cluster,theta,w_theta,errs,bias,biaserr,minscale=minscale,maxscale=maxscale,$
                      /fitplaws,plawplot='ang_cluster_plaw.png'
   ENDELSE
ENDIF



;MAD Get finish, elapsed time
et=systime(1)
print,'Ang_cluster - Elapsed time= ',strtrim((et-st)/60,2),' min'

return

END
