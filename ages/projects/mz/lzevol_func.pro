function lzevol_func, absmag, params, z=z, qz0=qz0, pivotmag=pivotmag
; simple linear evolution model (see MLFIT_MZLZEVOL)
    c0 = params[0]
    c1 = params[1]
    ohevol = params[2]
    mbevol = params[3]
    model = c0 + c1*(absmag-pivotmag) + (ohevol+c1*mbevol)*(z-qz0)
return, model
end
