function init_lf24_parinfo, alpha=alpha, beta=beta, phistar=phistar
; jm10mar25ucsd - initialize the L(24) LF fitting parameters

    if (n_elements(alpha) eq 0) then alpha = 0.37D
    if (n_elements(beta) eq 0) then beta = 0.37D
    if (n_elements(phistar) eq 0) then phistar = 1.2D-3
    
    parinfo = {value: 0.0D, limited: [0,0], $
      limits: [0.0D,0.0D], fixed: 0}
    parinfo = replicate(parinfo,4)
    
    parinfo[0].value = phistar
    parinfo[1].value = alog10(4D9)
    parinfo[2].value = alpha
    parinfo[3].value = beta
    parinfo[[0,2,3]].fixed = 1 ; by default just fit for L*

return, parinfo
end
