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
;                      data2_fileout='data2.fits',/figures
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
;-
PRO split_regions_gen,data_in,rand_in,data_out,rand_out,N=N,frac=frac,$
                      data_fileout=data_fileout,rand_fileout=rand_fileout,$
                      data2_in=data2_in,data2_fileout=data2_fileout,figures=figures

  ;MAD Set defaults
  IF ~keyword_set(N) THEN N=16

  ;MAD if given file names, read in
  IF (size(data_in,/type) EQ 7) THEN data=mrdfits(data_in,1) ELSE data=data_in
  IF (size(rand_in,/type) EQ 7) THEN rand=mrdfits(rand_in,1) ELSE rand=rand_in
  IF keyword_set(data2_in) THEN BEGIN
     IF (size(data2_in,/type) EQ 7) THEN data2=mrdfits(data2_in,1) ELSE data2=data2_in
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
  
  ;MAD Set total number of regions, number of cuts
  ntot=n^2.
  n_dec_cuts=n-1.
  n_ra_cuts=n-1.

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

  ;MAD Add min and max decs to cuts, split up by region
  dec_cuts=[min(rand_out.dec),dec_cuts,max(rand_out.dec)+0.1]
  FOR i=0L,n_elements(dec_cuts)-2 DO BEGIN
     xx=where(rand_out.dec GE dec_cuts[i] AND rand_out.dec LT dec_cuts[i+1])
     rand_out[xx].reg=i+1
     xx=where(data_out.dec GE dec_cuts[i] AND data_out.dec LT dec_cuts[i+1])
     data_out[xx].reg=i+1
     IF keyword_set(data2_in) THEN BEGIN
        xx=where(data2_out.dec GE dec_cuts[i] AND data2_out.dec LT dec_cuts[i+1])
        data2_out[xx].reg=i+1
     ENDIF
  ENDFOR

  ;MAD Loop over splits by DEC to split in ra.
  ;MAD The "w" index counts final regions
  w=1.
  FOR i=0L,n-1 DO BEGIN
     print,'SPLIT_REGIONS - finding splits in RA, DEC strip '+strtrim(i+1,2)+'...'
     used=where(data_out.reg EQ i+1)
     usedata=data_out[used]
     IF keyword_set(data2_in) THEN BEGIN
        used2=where(data2_out.reg EQ i+1)
        usedata2=data2_out[used2]
     ENDIF
     user=where(rand_out.reg EQ i+1)
     userand=rand_out[user]
     ras=userand[bsort(userand.ra)].ra
     ;MAD Only use subset if frac is set
     IF keyword_set(frac) THEN BEGIN
        indx=randomu(123,n_elements(ras)*frac)*(n_elements(ras)-1)
        ras=ras[indx]
        ras=ras[bsort(ras)]
     ENDIF

     k=0
     j=0L
     WHILE (n_elements(ra_cuts) LT n_ra_cuts) DO BEGIN
        counter,j,n_elements(ras)
        tmp=n_elements(where(ras LT ras[j]))*(1./n_elements(ras))
        IF (tmp GE ((k+1.)*(ntot/n))/ntot) THEN BEGIN
           IF (n_elements(ra_cuts) EQ 0) THEN ra_cuts=ras[j] ELSE $
              ra_cuts=[ra_cuts,ras[j]]
           k=k+1
        ENDIF
        j=j+1
     ENDWHILE

     ;MAD Add min and max RA, loop and give unique region names
     ;MAD That can't be the same as the ones used for the first
     ;MAD declination regions (which is why they are negative)
     ra_cuts=[min(userand.ra),ra_cuts,max(userand.ra)+0.1]

     FOR v=0L,n_elements(ra_cuts)-2 DO BEGIN
        xx=where(rand_out[user].ra GE ra_cuts[v] AND rand_out[user].ra LT ra_cuts[v+1])
        rand_out[user[xx]].reg=w*(-1)
        xx=where(data_out[used].ra GE ra_cuts[v] AND data_out[used].ra LT ra_cuts[v+1])
        data_out[used[xx]].reg=w*(-1)
        IF keyword_set(data2_in) THEN BEGIN
           xx=where(data2_out[used2].ra GE ra_cuts[v] AND data2_out[used2].ra LT ra_cuts[v+1])
           data2_out[used2[xx]].reg=w*(-1)
        ENDIF
        w=w+1
     ENDFOR
     ;MAD Delete the RA cuts so the WHILE loop above works again
     IF (n_elements(ra_cuts_save) EQ 0) THEN ra_cuts_save=ra_cuts ELSE $
        ra_cuts_save=[ra_cuts_save,ra_cuts]
     undefine,ra_cuts
  ENDFOR

  ;MAD Fix the region indices so they are positive
  rand_out.reg=rand_out.reg*(-1)
  data_out.reg=data_out.reg*(-1)
  IF keyword_set(data2_in) THEN data2_out.reg=data2_out.reg*(-1)
  
  ;MAD Write out files if needed
  IF keyword_set(data_fileout) THEN mwrfits,data_out,data_fileout,/create
  IF keyword_set(rand_fileout) THEN mwrfits,rand_out,rand_fileout,/create
  IF keyword_set(data2_fileout) THEN mwrfits,data2_out,data2_fileout,/create
  
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
  ENDIF
  return
END
