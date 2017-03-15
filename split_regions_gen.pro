;+
;  NAME:
;    split_regions_gen
;  PURPOSE:
;    Split a set of data into an N by N grid of equal area points.
;    Works in the presence of masked data.  Can't be N by M.
;
;  USE:
;    split_regions_gen,data_in,rand_in,data_out,rand_out,$
;                      N=n,frac=frac,data_fileout='data.fits',$
;                      rand_fileout='rand.fits',data2_in=data2,$
;                      data2_fileout='data2.fits',$
;                      split_file='splits_file.txt',$
;                      /figures
;
;  INPUT:
;    data_in - structure or fits file name of data.  Must have tags RA
;              and DEC
;    rand_in - structure or fits file of random catalog.  Must have
;              tags RA and DEC
;    
;  OPTIONAL INPUT:
;    N - number of regions per side (Defaults to 4)
;    frac - fraction of random catalog to use when finding
;           boundaries. Can be slow for large catalogs, so this is an
;           option to speed at the cost of area accuracy.
;    data_fileout - name of fits file to write data to, with
;                   additional "region" tag
;    rand_fileout - name of fits file to write randoms to, with
;                   additional "region" tag
;    data2_in - structure or fits filename of second data set (if
;               doing a cross-correlation)
;    data2_fileout - name of fits file to write second data set to,
;                    with additional region tag
;    rand2_in - structure or fits filename of second random catalog
;    rand2_filout - filename to write second random catalog to
;    split_file - name of file for location of splits. If already
;                 exists, will read in and apply.  If not, will
;                 write results to file.  Format is pixel/region,lower
;                 dec limit, upper dec limit, lower ra limit, upper ra limit.
;
;
;  KEYWORDS:
;    figures - make plots highlighting pixelization of data/randoms
;
;  OUTPUT:
;    data_out - structure of data with additional region tag
;    rand_out - structure of rands with additional region tag
;
;  NOTES:
;
;  HISTORY:
;    8-31-16 - Written - MAD (Dartmouth)
;    9-13-16 - Added ra and dec split in/out option - MAD (Dartmouth)
;    3-15-17 - Added support for second random catalog for cross-corr
;              - MAD (Dartmouth)
;-
PRO split_regions_gen,data_in,rand_in,data_out,rand_out,N=N,frac=frac,$
                      data_fileout=data_fileout,rand_fileout=rand_fileout,$
                      data2_in=data2_in,data2_fileout=data2_fileout,$
                      rand2_in=rand2_in,rand2_fileout=rand2_fileout,$
                      figures=figures,$
                      split_file=split_file

  ;MAD Set defaults
  IF ~keyword_set(N) THEN N=4

  ;MAD if given file names, read in
  IF (size(data_in,/type) EQ 7) THEN data=mrdfits(data_in,1) ELSE data=data_in
  IF (size(rand_in,/type) EQ 7) THEN rand=mrdfits(rand_in,1) ELSE rand=rand_in
  IF keyword_set(data2_in) THEN BEGIN
     IF (size(data2_in,/type) EQ 7) THEN data2=mrdfits(data2_in,1) ELSE data2=data2_in
  ENDIF
  IF keyword_set(rand2_in) THEN BEGIN
     IF (size(rand2_in,/type) EQ 7) THEN rand2=mrdfits(rand2_in,1) ELSE rand2=rand2_in
  ENDIF

  ;MAD Add region tag
  add={reg:0}
  adddata=replicate(add,n_elements(data))
  data_out=struct_addtags(data,adddata)
  addrand=replicate(add,n_elements(rand))
  rand_out=struct_addtags(rand,addrand)
  IF keyword_set(data2_in) THEN BEGIN
     adddata2=replicate(add,n_elements(data2))
     data2_out=struct_addtags(data2,adddata2)
  ENDIF
  IF keyword_set(rand2_in) THEN BEGIN
     addrand2=replicate(add,n_elements(rand2))
     rand2_out=struct_addtags(rand2,addrand2)
  ENDIF
  
  ;MAD Set total number of regions, number of cuts
  ntot=n^2.
  n_dec_cuts=n-1.
  n_ra_cuts=n-1.

  ;MAD if file supplied, use to apply splits, otherwise find them
  IF keyword_set(split_file) THEN check=file_search(split_file) ELSE check=''
  IF (check NE '') THEN BEGIN
     print,'SPLIT_REGIONS - using supplied RA and Dec cuts...'
     readcol,split_file,pix,dec_cuts_low,dec_cuts_high,ra_cuts_low,ra_cuts_high,format='I,D,D,D,D'
     FOR i=0,n_elements(pix)-1 DO BEGIN
        data_out[where((data_out.dec GE dec_cuts_low[i]) AND $
                       (data_out.dec LT dec_cuts_high[i]) AND $
                       (data_out.ra GE ra_cuts_low[i]) AND $
                       (data_out.ra LT ra_cuts_high[i]))].reg = pix[i]
        rand_out[where((rand_out.dec GE dec_cuts_low[i]) AND $
                       (rand_out.dec LT dec_cuts_high[i]) AND $
                       (rand_out.ra GE ra_cuts_low[i]) AND $
                       (rand_out.ra LT ra_cuts_high[i]))].reg = pix[i]
        IF keyword_set(data2_in) THEN BEGIN
           data2_out[where((data2_out.dec GE dec_cuts_low[i]) AND $
                           (data2_out.dec LT dec_cuts_high[i]) AND $
                           (data2_out.ra GE ra_cuts_low[i]) AND $
                           (data2_out.ra LT ra_cuts_high[i]))].reg = pix[i]
        ENDIF
        IF keyword_set(rand2_in) THEN BEGIN
           rand2_out[where((rand2_out.dec GE dec_cuts_low[i]) AND $
                           (rand2_out.dec LT dec_cuts_high[i]) AND $
                           (rand2_out.ra GE ra_cuts_low[i]) AND $
                           (rand2_out.ra LT ra_cuts_high[i]))].reg = pix[i]
        ENDIF   
     ENDFOR
  ENDIF ELSE BEGIN  
     ;MAD open output file if set
     IF keyword_set(split_file) THEN openw,lun,split_file,/get_lun
     ;MAD sort declinations to find splits
     print,'SPLIT_REGIONS - finding splits in DEC...'
     decs=rand_out[bsort(rand_out.dec)].dec
     ;MAD Only use some of the random sample if frac is set
     IF keyword_set(frac) THEN BEGIN
        indx=randomu(718,n_elements(decs)*frac)*(n_elements(decs)-1)
        decs=decs[indx]
        decs=decs[bsort(decs)]
     ENDIF
     ;MAD Loop over sorted decs to find splits
     j=0
     i=0L
     WHILE (n_elements(dec_cuts) LT n_dec_cuts) DO BEGIN
        counter,i,n_elements(decs)
        tmp=n_elements(where(decs LT decs[i]))*1./n_elements(decs)
        IF (tmp GE ((j+1.)*(ntot/n))/ntot) THEN BEGIN
           IF (n_elements(dec_cuts) EQ 0) THEN dec_cuts=decs[i] ELSE $
              dec_cuts=[dec_cuts,decs[i]]
           j=j+1
        ENDIF
        i=i+1
     ENDWHILE  

     ;MAD Add min and max decs to cuts
     dec_cuts=[min(rand_out.dec)-0.5,dec_cuts,max(rand_out.dec)+0.5]
     
     ;MAD Apply dec splits
     FOR i=0L,n_elements(dec_cuts)-2 DO BEGIN
        xx=where(rand_out.dec GE dec_cuts[i] AND rand_out.dec LT dec_cuts[i+1])
        rand_out[xx].reg=i+1
        xx=where(data_out.dec GE dec_cuts[i] AND data_out.dec LT dec_cuts[i+1])
        data_out[xx].reg=i+1
        IF keyword_set(data2_in) THEN BEGIN
           xx=where(data2_out.dec GE dec_cuts[i] AND data2_out.dec LT dec_cuts[i+1])
           data2_out[xx].reg=i+1
        ENDIF
        IF keyword_set(rand2_in) THEN BEGIN
           xx=where(rand2_out.dec GE dec_cuts[i] AND rand2_out.dec LT dec_cuts[i+1])
           rand2_out[xx].reg=i+1
        ENDIF
     ENDFOR

     ;MAD Loop over splits by DEC to split in ra.
     ;MAD The "w" index counts final regions
     w=1.
     count=1.
     print,'SPLIT_REGIONS - finding splits in RA, for each DEC strip...'
     FOR i=0L,n-1 DO BEGIN
        counter,i,n
        used=where(data_out.reg EQ i+1)
        usedata=data_out[used]
        IF keyword_set(data2_in) THEN BEGIN
           used2=where(data2_out.reg EQ i+1)
           usedata2=data2_out[used2]
        ENDIF
        user=where(rand_out.reg EQ i+1)
        userand=rand_out[user]
        ras=userand[bsort(userand.ra)].ra
        IF keyword_set(rand2_in) THEN BEGIN
           user2=where(rand2_out.reg EQ i+1)
           userand2=rand2_out[user2]
        ENDIF        
        ;MAD Only use subset if frac is set
        IF keyword_set(frac) THEN BEGIN
           indx=randomu(123,n_elements(ras)*frac)*(n_elements(ras)-1)
           ras=ras[indx]
           ras=ras[bsort(ras)]
        ENDIF

        k=0
        j=0L
        WHILE (n_elements(ra_cuts) LT n_ra_cuts) DO BEGIN
           tmp=n_elements(where(ras LT ras[j]))*(1./n_elements(ras))
           IF (tmp GE ((k+1.)*(ntot/n))/ntot) THEN BEGIN
              IF (n_elements(ra_cuts) EQ 0) THEN ra_cuts=ras[j] ELSE $
                 ra_cuts=[ra_cuts,ras[j]]
              k=k+1
           ENDIF
           j=j+1
        ENDWHILE

        ;MAD Add min and max RA, write out if needed
        ra_cuts=[min(userand.ra)-0.5,ra_cuts,max(userand.ra)+0.5]
        IF keyword_set(split_file) THEN BEGIN
           FOR j=0L,n_elements(ra_cuts)-2 DO BEGIN
              printf,lun,count,dec_cuts[i],dec_cuts[i+1],ra_cuts[j],ra_cuts[j+1],$
                     format='(I,1x,D,1x,D,1x,D,1x,D)'
              count=count+1
           ENDFOR
        ENDIF

        ;MAD Loop and give unique region names
        ;MAD that can't be the same as the ones used for the first
        ;MAD declination regions (which is why they are negative)
        FOR v=0L,n_elements(ra_cuts)-2 DO BEGIN
           xx=where(rand_out[user].ra GE ra_cuts[v] AND rand_out[user].ra LT ra_cuts[v+1])
           rand_out[user[xx]].reg=w*(-1)
           xx=where(data_out[used].ra GE ra_cuts[v] AND data_out[used].ra LT ra_cuts[v+1])
           data_out[used[xx]].reg=w*(-1)
           IF keyword_set(data2_in) THEN BEGIN
              xx=where(data2_out[used2].ra GE ra_cuts[v] AND data2_out[used2].ra LT ra_cuts[v+1])
              data2_out[used2[xx]].reg=w*(-1)
           ENDIF
           IF keyword_set(rand2_in) THEN BEGIN
              xx=where(rand2_out[user2].ra GE ra_cuts[v] AND rand2_out[user2].ra LT ra_cuts[v+1])
              rand2_out[user2[xx]].reg=w*(-1)
           ENDIF
           w=w+1
        ENDFOR
       ;MAD Delete the RA cuts so the WHILE loop above works again
        IF (n_elements(ra_cuts_save) EQ 0) THEN ra_cuts_save=ra_cuts ELSE $
           ra_cuts_save=[ra_cuts_save,ra_cuts]
        undefine,ra_cuts
     ENDFOR
     close,lun
     
     ;MAD Fix the region indices so they are positive
     rand_out.reg=rand_out.reg*(-1)
     data_out.reg=data_out.reg*(-1)
     IF keyword_set(data2_in) THEN data2_out.reg=data2_out.reg*(-1)
     IF keyword_set(rand2_in) THEN rand2_out.reg=rand2_out.reg*(-1)
  ENDELSE

  ;MAD Write out files if needed
  IF keyword_set(data_fileout) THEN mwrfits,data_out,data_fileout,/create
  IF keyword_set(rand_fileout) THEN mwrfits,rand_out,rand_fileout,/create
  IF keyword_set(data2_fileout) THEN mwrfits,data2_out,data2_fileout,/create
  IF keyword_set(rand2_fileout) THEN mwrfits,rand2_out,rand2_fileout,/create
  
  ;MAD Make some plots as checks
  IF keyword_set(figures) THEN BEGIN
     colseed=123
     indices=rand_out.reg
     indices=indices[rem_dup(indices)]
     PS_start,filename='rand_regions.png'
     plot,rand_out.ra,rand_out.dec,xtit='RA',ytit='Dec',psym=3
     loadct,2
     FOR i=0L,n_elements(indices)-1 DO BEGIN
        xx=where(rand_out.reg EQ indices[i])
        oplot,rand_out[xx].ra,rand_out[xx].dec,psym=3,color=randomu(colseed)*!D.TABLE_SIZE
     ENDFOR
     PS_end,/png

     PS_start,filename='data_regions.png'
     plot,data_out.ra,data_out.dec,xtit='RA',ytit='Dec',psym=3
     loadct,2
     FOR i=0L,n_elements(indices)-1 DO BEGIN
        xx=where(data_out.reg EQ indices[i])
        oplot,data_out[xx].ra,data_out[xx].dec,psym=3,color=randomu(colseed)*!D.TABLE_SIZE
     ENDFOR
     PS_end,/png

     IF keyword_set(data2) THEN BEGIN
        PS_start,filename='data2_regions.png'
        plot,data2_out.ra,data2_out.dec,xtit='RA',ytit='Dec',psym=3
        loadct,2
        FOR i=0L,n_elements(indices)-1 DO BEGIN
           xx=where(data2_out.reg EQ indices[i])
           oplot,data2_out[xx].ra,data2_out[xx].dec,psym=3,color=randomu(colseed)*!D.TABLE_SIZE
        ENDFOR
        PS_end,/png
     ENDIF
     IF keyword_set(rand2_in) THEN BEGIN
        PS_start,filename='rand2_regions.png'
        plot,rand2_out.ra,rand2_out.dec,xtit='RA',ytit='Dec',psym=3
        loadct,2
        FOR i=0L,n_elements(indices)-1 DO BEGIN
           xx=where(rand2_out.reg EQ indices[i])
           oplot,rand2_out[xx].ra,rand2_out[xx].dec,psym=3,color=randomu(colseed)*!D.TABLE_SIZE
        ENDFOR
        PS_end,/png
     ENDIF
  ENDIF
  return
END
