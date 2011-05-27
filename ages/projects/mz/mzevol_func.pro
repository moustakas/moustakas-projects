function mzevol_func, logmass, params, z=z1, qz0=qz0
; model: evolutionary model (see MLFIT_MZLZEVOL)
    npar = n_elements(params)
    ngal = n_elements(logmass)
    if (n_elements(z1) eq 1) then z = replicate(z1,ngal) else z = z1

    mratio = (params[1]*10^(params[3]*(z-qz0)))/(10^(logmass-9D))
    ohmodel = params[0] + params[4]*(z-qz0) - alog10(1+mratio^params[2])
;   plot, logmass, ohmodel, psym=6, ysty=3 & cc = get_kbrd(1)    

return, ohmodel
end
