PRO split_regions,data_in,rand_in,data_fileout,rand_fileout,data2_in=data2_in,data2_fileout=data2_fileout,figures=figures
;-----------------------------------------------------------------------------
;This procedure will split up a set of positional data and
;accompanying random catalog into 16 regions with equal numbers of
;random catalog points (and thus roughly equal area).  Will output new
;fits structures with the original ra, dec, and an additional 'reg'
;tag that identifies which region each point is in.
;
;INPUTS
;data_filein - string name of the fits file containing the real data.  Should
;              be a structure with at least ra and dec tags.
;
;rand_filein - string tame of the fits file containing the random
;              data. Should be a structure with at least ra,dec tags.
;
;data_fileout - the name of the file to output the new data structure with
;                ra, dec, and reg tags
;
;rand_fileout - the name of the file to output the new radnom structure with
;               ra, dec, and reg tags
;
;OPTIONAL INPUTS
;data2 - second data set if doing a cross-correlation
;data2_fileout - filename for output of second data set
;
;OPTIONAL KEYWORDS
;figures - set this flag to output PNG format figures showing all of
;          the data and random points color-coded by region
;
;HISTORY
;   10-1-2012 - Written - MAD (UWyo)
;   7-20-2015 - Added data2 for cross-correlation, removed unneeded loops - MAD (UWyo)
;-----------------------------------------------------------------------------
IF (n_elements(rand_fileout) EQ 0) THEN message,'Syntax - split_regions,''data.fits'',''rand.fits'',''data_out.fits'',''rand_out.fits''[,data2=data2,data2_fileout=''data2_out.fits'',/figures]'

;MAD Get start time
st=systime(1)

;MAD Set input data
data=data_in
rand=rand_in
IF keyword_set(data2_in) THEN data2=data2_in

;MAD Make new structures with region flags
randregion={ra:0., dec:0., reg:0}
randregion=replicate(randregion,n_elements(rand))
randregion.ra=rand.ra
randregion.dec=rand.dec

dataregion={ra:0.,dec:0.,reg:0}
dataregion=replicate(dataregion,n_elements(data))
dataregion.ra=data.ra
dataregion.dec=data.dec

IF keyword_set(data2) THEN BEGIN
   data2region={ra:0.,dec:0.,reg:0}
   data2region=replicate(data2region,n_elements(data2))
   data2region.ra=data2.ra
   data2region.dec=data2.dec
ENDIF

print,'Split_regions - finding RA and DEC cuts...'
;MAD Determine first split in DEC
dec_split1=median(rand.dec)
xx=[where(rand.dec GE dec_split1)]
yy=[where(rand.dec LT dec_split1)]
;MAD Determine first split in RA
ra_split1_1=median(rand[xx].ra)
ra_split1_2=median(rand[yy].ra)

;MAD Determine second DEC splits
xx=where((rand.dec GE dec_split1) AND (rand.ra LT ra_split1_1))
dec_split2_1=median(rand[xx].dec)
xx=where((rand.dec GE dec_split1) AND (rand.ra GE ra_split1_1))
dec_split2_2=median(rand[xx].dec)
xx=where((rand.dec LT dec_split1) AND (rand.ra LT ra_split1_1))
dec_split2_3=median(rand[xx].dec)
xx=where((rand.dec LT dec_split1) AND (rand.ra GE ra_split1_1))
dec_split2_4=median(rand[xx].dec)

;MAD Determine second RA splits
xx=where((rand.dec GE dec_split2_1) AND (rand.ra LT ra_split1_1))
ra_split2_1=median(rand[xx].ra)
xx=where((rand.dec LT dec_split2_1) AND (rand.dec GE dec_split1) AND (rand.ra LT ra_split1_1))
ra_split2_2=median(rand[xx].ra)
xx=where((rand.dec LT dec_split1) AND (rand.dec GE dec_split2_3) AND (rand.ra LT ra_split1_2))
ra_split2_3=median(rand[xx].ra)
xx=where((rand.dec LT dec_split2_3) AND (rand.ra LT ra_split1_2))
ra_split2_4=median(rand[xx].ra)
xx=where((rand.dec GE dec_split2_2) AND (rand.ra GE ra_split1_1))
ra_split2_5=median(rand[xx].ra)
xx=where((rand.dec LT dec_split2_2) AND (rand.dec GE dec_split1) AND (rand.ra GE ra_split1_1))
ra_split2_6=median(rand[xx].ra)
xx=WHERE((rand.dec LT dec_split1) AND (rand.dec GE dec_split2_4) AND (rand.ra GE ra_split1_2))
ra_split2_7=median(rand[xx].ra)
xx=where((rand.dec LT dec_split2_4) AND (rand.ra GE ra_split1_2))
ra_split2_8=median(rand[xx].ra)


;MAD Assign region to each data point
print,'Split_regions - Assigning region flags to each data point...'
xx=where((dataregion.ra LT ra_split2_1) AND (dataregion.dec GE dec_split2_1),cnt)
IF (cnt NE 0) THEN dataregion[xx].reg=1
xx=where((dataregion.ra GE ra_split2_1) AND (dataregion.ra LT ra_split1_1) AND (dataregion.dec GE dec_split2_1),cnt)
IF (cnt NE 0) THEN dataregion[xx].reg=2
xx=where((dataregion.ra GE ra_split1_1) AND (dataregion.ra LT ra_split2_5) AND (dataregion.dec GE dec_split2_2),cnt)
IF (cnt NE 0) THEN dataregion[xx].reg=3
xx=where((dataregion.ra GE ra_split2_5) AND (dataregion.dec GE dec_split2_2),cnt)
IF (cnt NE 0) THEN dataregion[xx].reg=4
xx=where((dataregion.ra LT ra_split2_2) AND (dataregion.dec GE dec_split1) AND (dataregion.dec LT dec_split2_1),cnt)
IF (cnt NE 0) THEN dataregion[xx].reg=5
xx=where((dataregion.ra GE ra_split2_2) AND (dataregion.ra LT ra_split1_1) AND (dataregion.dec GE dec_split1) AND (dataregion.dec LT dec_split2_1),cnt)
IF (cnt NE 0) THEN dataregion[xx].reg=6
xx=where((dataregion.ra GE ra_split1_1) AND (dataregion.ra LT ra_split2_6) AND (dataregion.dec GE dec_split1) AND (dataregion.dec LT dec_split2_2),cnt)
IF (cnt NE 0) THEN dataregion[xx].reg=7
xx=where((dataregion.ra GE ra_split2_6) AND (dataregion.dec GE dec_split1) AND (dataregion.dec LT dec_split2_2),cnt)
IF (cnt NE 0) THEN dataregion[xx].reg=8
xx=where((dataregion.ra LT ra_split2_3) AND (dataregion.dec GE dec_split2_3) AND (dataregion.dec LT dec_split1),cnt)
IF (cnt NE 0) THEN dataregion[xx].reg=9
xx=where((dataregion.ra GE ra_split2_3) AND (dataregion.ra LT ra_split1_2) AND (dataregion.dec GE dec_split2_3) AND (dataregion.dec LT dec_split1),cnt)
IF (cnt NE 0) THEN dataregion[xx].reg=10
xx=where((dataregion.ra GE ra_split1_2) AND (dataregion.ra LT ra_split2_7) AND (dataregion.dec GE dec_split2_4) AND (dataregion.dec LT dec_split1),cnt)
IF (cnt NE 0) THEN dataregion[xx].reg=11
xx=where((dataregion.ra GE ra_split2_7) AND (dataregion.dec GE dec_split2_4) AND (dataregion.dec LT dec_split1),cnt)
IF (cnt NE 0) THEN dataregion[xx].reg=12
xx=where((dataregion.ra LT ra_split2_4) AND (dataregion.dec LT dec_split2_3),cnt)
IF (cnt NE 0) THEN dataregion[xx].reg=13
xx=where((dataregion.ra GE ra_split2_4) AND (dataregion.ra LT ra_split1_2) AND (dataregion.dec LT dec_split2_3),cnt)
IF (cnt NE 0) THEN dataregion[xx].reg=14
xx=where((dataregion.ra GE ra_split1_2) AND (dataregion.ra LT ra_split2_8) AND (dataregion.dec LT dec_split2_4),cnt)
IF (cnt NE 0) THEN dataregion[xx].reg=15
xx=where((dataregion.ra GE ra_split2_8) AND (dataregion.dec LT dec_split2_4),cnt)
IF (cnt NE 0) THEN dataregion[xx].reg=16

;MAD Write out data file with region data
mwrfits,dataregion,data_fileout,/create


;MAD Assign region to each random point
print,'Split_regions - Assigning region flags to each random point...'
xx=where((randregion.ra LT ra_split2_1) AND (randregion.dec GE dec_split2_1),cnt)
IF (cnt NE 0) THEN randregion[xx].reg=1
xx=where((randregion.ra GE ra_split2_1) AND (randregion.ra LT ra_split1_1) AND (randregion.dec GE dec_split2_1),cnt)
IF (cnt NE 0) THEN randregion[xx].reg=2
xx=where((randregion.ra GE ra_split1_1) AND (randregion.ra LT ra_split2_5) AND (randregion.dec GE dec_split2_2),cnt)
IF (cnt NE 0) THEN randregion[xx].reg=3
xx=where((randregion.ra GE ra_split2_5) AND (randregion.dec GE dec_split2_2),cnt)
IF (cnt NE 0) THEN randregion[xx].reg=4
xx=where((randregion.ra LT ra_split2_2) AND (randregion.dec GE dec_split1) AND (randregion.dec LT dec_split2_1),cnt)
IF (cnt NE 0) THEN randregion[xx].reg=5
xx=where((randregion.ra GE ra_split2_2) AND (randregion.ra LT ra_split1_1) AND (randregion.dec GE dec_split1) AND (randregion.dec LT dec_split2_1),cnt)
IF (cnt NE 0) THEN randregion[xx].reg=6
xx=where((randregion.ra GE ra_split1_1) AND (randregion.ra LT ra_split2_6) AND (randregion.dec GE dec_split1) AND (randregion.dec LT dec_split2_2),cnt)
IF (cnt NE 0) THEN randregion[xx].reg=7
xx=where((randregion.ra GE ra_split2_6) AND (randregion.dec GE dec_split1) AND (randregion.dec LT dec_split2_2),cnt)
IF (cnt NE 0) THEN randregion[xx].reg=8
xx=where((randregion.ra LT ra_split2_3) AND (randregion.dec GE dec_split2_3) AND (randregion.dec LT dec_split1),cnt)
IF (cnt NE 0) THEN randregion[xx].reg=9
xx=where((randregion.ra GE ra_split2_3) AND (randregion.ra LT ra_split1_2) AND (randregion.dec GE dec_split2_3) AND (randregion.dec LT dec_split1),cnt)
IF (cnt NE 0) THEN randregion[xx].reg=10
xx=where((randregion.ra GE ra_split1_2) AND (randregion.ra LT ra_split2_7) AND (randregion.dec GE dec_split2_4) AND (randregion.dec LT dec_split1),cnt)
IF (cnt NE 0) THEN randregion[xx].reg=11
xx=where((randregion.ra GE ra_split2_7) AND (randregion.dec GE dec_split2_4) AND (randregion.dec LT dec_split1),cnt)
IF (cnt NE 0) THEN randregion[xx].reg=12
xx=where((randregion.ra LT ra_split2_4) AND (randregion.dec LT dec_split2_3),cnt)
IF (cnt NE 0) THEN randregion[xx].reg=13
xx=where((randregion.ra GE ra_split2_4) AND (randregion.ra LT ra_split1_2) AND (randregion.dec LT dec_split2_3),cnt)
IF (cnt NE 0) THEN randregion[xx].reg=14
xx=where((randregion.ra GE ra_split1_2) AND (randregion.ra LT ra_split2_8) AND (randregion.dec LT dec_split2_4),cnt)
IF (cnt NE 0) THEN randregion[xx].reg=15
xx=where((randregion.ra GE ra_split2_8) AND (randregion.dec LT dec_split2_4),cnt)
IF (cnt NE 0) THEN randregion[xx].reg=16

;MAD Write out random file with region info
mwrfits,randregion,rand_fileout,/create



;MAD Do second data set if doing cross-corr
IF keyword_set(data2) THEN BEGIN
   xx=where((data2region.ra LT ra_split2_1) AND (data2region.dec GE dec_split2_1),cnt)
   IF (cnt NE 0) THEN data2region[xx].reg=1
   xx=where((data2region.ra GE ra_split2_1) AND (data2region.ra LT ra_split1_1) AND (data2region.dec GE dec_split2_1),cnt)
   IF (cnt NE 0) THEN data2region[xx].reg=2
   xx=where((data2region.ra GE ra_split1_1) AND (data2region.ra LT ra_split2_5) AND (data2region.dec GE dec_split2_2),cnt)
   IF (cnt NE 0) THEN data2region[xx].reg=3
   xx=where((data2region.ra GE ra_split2_5) AND (data2region.dec GE dec_split2_2),cnt)
   IF (cnt NE 0) THEN data2region[xx].reg=4
   xx=where((data2region.ra LT ra_split2_2) AND (data2region.dec GE dec_split1) AND (data2region.dec LT dec_split2_1),cnt)
   IF (cnt NE 0) THEN data2region[xx].reg=5
   xx=where((data2region.ra GE ra_split2_2) AND (data2region.ra LT ra_split1_1) AND (data2region.dec GE dec_split1) AND (data2region.dec LT dec_split2_1),cnt)
   IF (cnt NE 0) THEN data2region[xx].reg=6
   xx=where((data2region.ra GE ra_split1_1) AND (data2region.ra LT ra_split2_6) AND (data2region.dec GE dec_split1) AND (data2region.dec LT dec_split2_2),cnt)
   IF (cnt NE 0) THEN data2region[xx].reg=7
   xx=where((data2region.ra GE ra_split2_6) AND (data2region.dec GE dec_split1) AND (data2region.dec LT dec_split2_2),cnt)
   IF (cnt NE 0) THEN data2region[xx].reg=8
   xx=where((data2region.ra LT ra_split2_3) AND (data2region.dec GE dec_split2_3) AND (data2region.dec LT dec_split1),cnt)
   IF (cnt NE 0) THEN data2region[xx].reg=9
   xx=where((data2region.ra GE ra_split2_3) AND (data2region.ra LT ra_split1_2) AND (data2region.dec GE dec_split2_3) AND (data2region.dec LT dec_split1),cnt)
   IF (cnt NE 0) THEN data2region[xx].reg=10
   xx=where((data2region.ra GE ra_split1_2) AND (data2region.ra LT ra_split2_7) AND (data2region.dec GE dec_split2_4) AND (data2region.dec LT dec_split1),cnt)
   IF (cnt NE 0) THEN data2region[xx].reg=11
   xx=where((data2region.ra GE ra_split2_7) AND (data2region.dec GE dec_split2_4) AND (data2region.dec LT dec_split1),cnt)
   IF (cnt NE 0) THEN data2region[xx].reg=12
   xx=where((data2region.ra LT ra_split2_4) AND (data2region.dec LT dec_split2_3),cnt)
   IF (cnt NE 0) THEN data2region[xx].reg=13
   xx=where((data2region.ra GE ra_split2_4) AND (data2region.ra LT ra_split1_2) AND (data2region.dec LT dec_split2_3),cnt)
   IF (cnt NE 0) THEN data2region[xx].reg=14
   xx=where((data2region.ra GE ra_split1_2) AND (data2region.ra LT ra_split2_8) AND (data2region.dec LT dec_split2_4),cnt)
   IF (cnt NE 0) THEN data2region[xx].reg=15
   xx=where((data2region.ra GE ra_split2_8) AND (data2region.dec LT dec_split2_4),cnt)
   IF (cnt NE 0) THEN data2region[xx].reg=16

   mwrfits,data2region,data2_fileout,/create
ENDIF








;MAD Make figures if keyword is set
IF (keyword_set(figures)) THEN BEGIN
   print,'Split_regions - making and writing out figures...'
   PS_start,filename='data_regions.png'
   nice_plot,min(randregion.ra)-5.,max(randregion.ra)+5.,min(randregion.dec)-5.,max(randregion.dec)+5.,xtit='RA (deg)',ytit='DEC (deg)'
   cols=['grey','blue','yellow','red','orchid','sienna','green','aquamarine','hot pink','cyan','olive','charcoal','maroon','khaki','orange','violet']
   FOR i=1L,max(dataregion.reg) DO BEGIN
      oplot,data[WHERE(dataregion.reg EQ i)].ra,data[WHERE(dataregion.reg EQ i)].dec,psym=3,color=cgcolor(cols[i-1])
   ENDFOR
   PS_End,/png

   PS_start,filename='rand_regions.png'
   nice_plot,min(randregion.ra)-5.,max(randregion.ra)+5.,min(randregion.dec)-5.,max(randregion.dec)+5.,xtit='RA (deg)',ytit='DEC (deg)'
   FOR i=1L,max(randregion.reg) DO BEGIN
      oplot,randregion[WHERE(randregion.reg EQ i)].ra,randregion[WHERE(randregion.reg EQ i)].dec,psym=3,color=cgcolor(cols[i-1])
   ENDFOR
   PS_End,/png

   IF keyword_set(data2) THEN BEGIN
      PS_start,filename='data2_regions.png'
      nice_plot,min(randregion.ra)-5.,max(randregion.ra)+5.,min(randregion.dec)-5.,max(randregion.dec)+5.,xtit='RA (deg)',ytit='DEC (deg)'
      cols=['grey','blue','yellow','red','orchid','sienna','green','aquamarine','hot pink','cyan','olive','charcoal','maroon','khaki','orange','violet']
      FOR i=1L,max(data2region.reg) DO BEGIN
         oplot,data2[WHERE(data2region.reg EQ i)].ra,data2[WHERE(data2region.reg EQ i)].dec,psym=3,color=cgcolor(cols[i-1])
      ENDFOR
      PS_End,/png
   ENDIF
ENDIF

;MAD Get finish, elapsed time
et=systime(1)
print,'Split regions - Elapsed time= ',strtrim((et-st)/60,2),' min'

END
