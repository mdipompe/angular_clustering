PRO split_regions,data_in,rand_in,data_fileout,rand_fileout,figures=figures

;-----------------------------------------------------------------------------
;Written by MAD sometime in 2013, finalized 12/13
;
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
;OPTIONAL KEYWORDS
;figures - set this flag to output PNG format figures showing all of
;          the data and random points color-coded by region
;-----------------------------------------------------------------------------
IF (n_elements(rand_fileout) EQ 0) THEN message,'Syntax - split_regions,''data.fits'',''rand.fits'',''data_out.fits'',''rand_out.fits''[,/figures]'

;MAD Get start time
st=systime(1)

;MAD Set input data
data=data_in
rand=rand_in

;MAD Make new structures with region flags
randregion={ra:0., dec:0., reg:0}
randregion=replicate(randregion,n_elements(rand))
randregion.ra=rand.ra
randregion.dec=rand.dec

dataregion={ra:0.,dec:0.,reg:0}
dataregion=replicate(dataregion,n_elements(data))
dataregion.ra=data.ra
dataregion.dec=data.dec

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
FOR i=0L,n_elements(data.ra)-1 DO BEGIN
 IF ((dataregion[i].ra LT ra_split2_1) AND (dataregion[i].dec GE dec_split2_1)) THEN dataregion[i].reg=1
 IF ((dataregion[i].ra GE ra_split2_1) AND (dataregion[i].ra LT ra_split1_1) AND (dataregion[i].dec GE dec_split2_1)) THEN dataregion[i].reg=2
 IF ((dataregion[i].ra GE ra_split1_1) AND (dataregion[i].ra LT ra_split2_5) AND (dataregion[i].dec GE dec_split2_2)) THEN dataregion[i].reg=3
 IF ((dataregion[i].ra GE ra_split2_5) AND (dataregion[i].dec GE dec_split2_2)) THEN dataregion[i].reg=4
 IF ((dataregion[i].ra LT ra_split2_2) AND (dataregion[i].dec GE dec_split1) AND (dataregion[i].dec LT dec_split2_1)) THEN dataregion[i].reg=5
 IF ((dataregion[i].ra GE ra_split2_2) AND (dataregion[i].ra LT ra_split1_1) AND (dataregion[i].dec GE dec_split1) AND (dataregion[i].dec LT dec_split2_1)) THEN dataregion[i].reg=6 
 IF ((dataregion[i].ra GE ra_split1_1) AND (dataregion[i].ra LT ra_split2_6) AND (dataregion[i].dec GE dec_split1) AND (dataregion[i].dec LT dec_split2_2)) THEN dataregion[i].reg=7 
 IF ((dataregion[i].ra GE ra_split2_6) AND (dataregion[i].dec GE dec_split1) AND (dataregion[i].dec LT dec_split2_2)) THEN dataregion[i].reg=8 
 IF ((dataregion[i].ra LT ra_split2_3) AND (dataregion[i].dec GE dec_split2_3) AND (dataregion[i].dec LT dec_split1)) THEN dataregion[i].reg=9  
 IF ((dataregion[i].ra GE ra_split2_3) AND (dataregion[i].ra LT ra_split1_2) AND (dataregion[i].dec GE dec_split2_3) AND (dataregion[i].dec LT dec_split1)) THEN dataregion[i].reg=10
 IF ((dataregion[i].ra GE ra_split1_2) AND (dataregion[i].ra LT ra_split2_7) AND (dataregion[i].dec GE dec_split2_4) AND (dataregion[i].dec LT dec_split1)) THEN dataregion[i].reg=11
 IF ((dataregion[i].ra GE ra_split2_7) AND (dataregion[i].dec GE dec_split2_4) AND (dataregion[i].dec LT dec_split1)) THEN dataregion[i].reg=12
 IF ((dataregion[i].ra LT ra_split2_4) AND (dataregion[i].dec LT dec_split2_3)) THEN dataregion[i].reg=13
 IF ((dataregion[i].ra GE ra_split2_4) AND (dataregion[i].ra LT ra_split1_2) AND (dataregion[i].dec LT dec_split2_3)) THEN dataregion[i].reg=14
 IF ((dataregion[i].ra GE ra_split1_2) AND (dataregion[i].ra LT ra_split2_8) AND (dataregion[i].dec LT dec_split2_4)) THEN dataregion[i].reg=15
 IF ((dataregion[i].ra GE ra_split2_8) AND (dataregion[i].dec LT dec_split2_4)) THEN dataregion[i].reg=16
ENDFOR

;MAD Write out data file with region data
mwrfits,dataregion,data_fileout,/create


;MAD Assign region to each random point
print,'Split_regions - Assigning region flags to each random point...'
FOR i=0L,n_elements(randregion.ra)-1 DO BEGIN
 IF ((randregion[i].ra LT ra_split2_1) AND (randregion[i].dec GE dec_split2_1)) THEN randregion[i].reg=1
 IF ((randregion[i].ra GE ra_split2_1) AND (randregion[i].ra LT ra_split1_1) AND (randregion[i].dec GE dec_split2_1)) THEN randregion[i].reg=2
 IF ((randregion[i].ra GE ra_split1_1) AND (randregion[i].ra LT ra_split2_5) AND (randregion[i].dec GE dec_split2_2)) THEN randregion[i].reg=3
 IF ((randregion[i].ra GE ra_split2_5) AND (randregion[i].dec GE dec_split2_2)) THEN randregion[i].reg=4
 IF ((randregion[i].ra LT ra_split2_2) AND (randregion[i].dec GE dec_split1) AND (randregion[i].dec LT dec_split2_1)) THEN randregion[i].reg=5
 IF ((randregion[i].ra GE ra_split2_2) AND (randregion[i].ra LT ra_split1_1) AND (randregion[i].dec GE dec_split1) AND (randregion[i].dec LT dec_split2_1)) THEN randregion[i].reg=6 
 IF ((randregion[i].ra GE ra_split1_1) AND (randregion[i].ra LT ra_split2_6) AND (randregion[i].dec GE dec_split1) AND (randregion[i].dec LT dec_split2_2)) THEN randregion[i].reg=7 
 IF ((randregion[i].ra GE ra_split2_6) AND (randregion[i].dec GE dec_split1) AND (randregion[i].dec LT dec_split2_2)) THEN randregion[i].reg=8 
 IF ((randregion[i].ra LT ra_split2_3) AND (randregion[i].dec GE dec_split2_3) AND (randregion[i].dec LT dec_split1)) THEN randregion[i].reg=9  
 IF ((randregion[i].ra GE ra_split2_3) AND (randregion[i].ra LT ra_split1_2) AND (randregion[i].dec GE dec_split2_3) AND (randregion[i].dec LT dec_split1)) THEN randregion[i].reg=10
 IF ((randregion[i].ra GE ra_split1_2) AND (randregion[i].ra LT ra_split2_7) AND (randregion[i].dec GE dec_split2_4) AND (randregion[i].dec LT dec_split1)) THEN randregion[i].reg=11
 IF ((randregion[i].ra GE ra_split2_7) AND (randregion[i].dec GE dec_split2_4) AND (randregion[i].dec LT dec_split1)) THEN randregion[i].reg=12
 IF ((randregion[i].ra LT ra_split2_4) AND (randregion[i].dec LT dec_split2_3)) THEN randregion[i].reg=13
 IF ((randregion[i].ra GE ra_split2_4) AND (randregion[i].ra LT ra_split1_2) AND (randregion[i].dec LT dec_split2_3)) THEN randregion[i].reg=14
 IF ((randregion[i].ra GE ra_split1_2) AND (randregion[i].ra LT ra_split2_8) AND (randregion[i].dec LT dec_split2_4)) THEN randregion[i].reg=15
 IF ((randregion[i].ra GE ra_split2_8) AND (randregion[i].dec LT dec_split2_4)) THEN randregion[i].reg=16
ENDFOR

;MAD Write out random file with region info
mwrfits,randregion,rand_fileout,/create

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
ENDIF

;MAD Get finish, elapsed time
et=systime(1)
print,'Split regions - Elapsed time= ',strtrim((et-st)/60,2),' min'

END
