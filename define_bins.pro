FUNCTION define_bins,nperdex,max

;Define step size of exponent
step=1./nperdex

;Calculate bin edges
power=-1.8
bin_edge=0.
WHILE (max(bin_edge) LT max*60.) DO BEGIN
   bin_edge=[bin_edge,10.^(power)]
   power=power+step
ENDWHILE

;Goes one too far, get rid of last one
bin_edge=bin_edge[1:n_elements(bin_edge)-2]

;Define output structure
bins={edges:0.,cents:0.}
bins=replicate(bins,n_elements(bin_edge))
bins.edges=bin_edge

;Calculate bin centers (yes, its odd that it doesnt have 
;one less element then edges, this is a relic that must
;remain for now for other codes...)
bin_cent=fltarr(n_elements(bin_edge))
power=(-1.8)+(step/2.)
FOR i=0L,n_elements(bin_cent)-1 DO BEGIN
   bin_cent[i]=10.^(power)
   power=power+step
ENDFOR

bins.cents=bin_cent

return,bins
END
