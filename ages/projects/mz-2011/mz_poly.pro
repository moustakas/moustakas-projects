function mz_poly, logmass, params
; return a simple polynomial
    norm = 10D
    model = poly(logmass-norm,params)
return, model
end

