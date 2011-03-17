function cosmicimf_init_mf_parinfo, alpha=alpha, fitslope=fitslope
; jm10mar23ucsd - initialize the MF fitting parameters

    parinfo = {value: 0.0D, limited: [0,0], $
      limits: [0.0D,0.0D], fixed: 0, step: 0.0D, $
      relstep: 0.0D}
    parinfo = replicate(parinfo,3)
    
    parinfo[0].value = 4.0D-3
;   parinfo[0].limited[0] = 1
;   parinfo[0].limits[0] = 1D-8
    
    parinfo[1].value = 10.5D
    parinfo[1].limited = 1
    parinfo[1].limits = [9.5D,11.5D]

    if (n_elements(alpha) eq 0) then parinfo[2].value = -1.1D else $
      parinfo[2].value = alpha
    parinfo[2].fixed = (keyword_set(fitslope) eq 0)

return, parinfo
end
