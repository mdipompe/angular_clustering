;+
;Created by MAD 1/21/13
; Generates an exponential model
; to be used by mpfitfunc for fitting a power law.
; Two parameters define the model,
; A the normalization and B the power.
;-
FUNCTION powerlaw,x,p

A=p[0]
B=p[1]
ymod=A*x^(B)

return,ymod

END
