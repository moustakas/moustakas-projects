
; weight is 1/error^2 and calculated outside code

FUNCTION GAPPY_ages, weight, flux, espec

if (size(espec))[0] eq 1 then message,'GAPPY needs more than one espectrum'

nrecon = (size(espec,/dim))[1]
nbin = (size(espec,/dim))[0]

tmp = (size(flux,/dim))
if n_elements(tmp) eq 2 then ngal = tmp[1] else if n_elements(tmp) eq 1 then ngal=1

pcs = dblarr(nrecon,ngal)

for j=0,ngal-1 do begin


    M = fltarr(nrecon,nrecon)
    M2 = fltarr(nrecon,nrecon)

    ;    for l=0,nbin-1 do $
    ;      M = M + weight[l]*MATRIX_MULTIPLY(espec[l,*],espec[l,*], /atranspose)
    ;this is the same
    M = MATRIX_MULTIPLY(espec*rebin(weight[*,j],nbin,nrecon),espec, /atranspose)
    
    F = MATRIX_MULTIPLY((weight[*,j]*(flux[*,j])),espec)
        
    M2 = invert(TRANSPOSE(M),status)
        
;if status = 1 or 2 there's a definite problem
    if status ne 0 then begin
        print, 'invert not worked', j,status
        continue                ;pcs will be 0.0 for this file
    endif

    pcs[*,j] = reform(MATRIX_MULTIPLY(F,M2,/BTRANSPOSE))



endfor

return,pcs

END
