function mzevol_func, logmass, params, z=z, qz0=qz0
; model: evolutionary model (see MLFIT_MZLZEVOL)

    newparams = params[0:2]
    newparams[0] = params[0] + params[3]*(z-qz0) ; O/H*
    newparams[1] = params[1] + params[4]*(z-qz0) ; M*
    
    ohmodel = mz_closedbox(logmass,newparams)
    
return, ohmodel
end
