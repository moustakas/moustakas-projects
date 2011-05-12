function mzevol_func, logmass, params, z=z1, qz0=qz0
; model: evolutionary model (see MLFIT_MZLZEVOL)
    npar = n_elements(params)
    ngal = n_elements(logmass)
    if (n_elements(z1) eq 1) then z = replicate(z1,ngal) else z = z1
    ohmodel = fltarr(ngal)
    for ii = 0, ngal-1 do begin
       newparams = params[0:2]
;      for bb = 4, npar-1 do newparams[0] += params[bb]*(z[ii]-qz0)^(bb-4) ; O/H*
;      for bb = 4, npar-1 do newparams[0] += params[bb]*(z[ii]-qz0)^(bb-3) ; O/H*
       newparams[0] = params[0] + params[4]*(z[ii]-qz0)    ; O/H*
       newparams[1] = params[1]*10^(params[3]*(z[ii]-qz0)) ; M*
       ohmodel[ii] = mz_closedbox(logmass[ii],newparams)
    endfor
return, ohmodel
end
