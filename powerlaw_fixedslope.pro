FUNCTION powerlaw_fixedslope,x,p

;-----------------------------------------
;Created by MAD 1/21/13
;This function generates an exponential model
;to be used by mpfitfunc for fitting decaying
;expononential lines.  Two parameters define the model;
;A the normalization and B the exponential factor.
;------------------------------------------

A=p[0]
B=-1.
ymod=A*(x^(B))

return,ymod

END
