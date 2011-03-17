function init_lf24lit, nz, nimf=nimf
    lit = {$
      color1: '', color2: '', psym: 0.0, symsize: 0.0, $
      z:       0.0, zerr:        0.0, $
      alpha:   0.0, alpha_err:   0.0, $
      beta:    0.0, beta_err:    0.0, $
      lstar:   0.0, lstar_err:   0.0, $
      phistar: 0.0, phistar_err: 0.0, phistar_lss_err: 0.0}
    if (n_elements(nimf) ne 0) then lit = struct_addtags(lit,$
      {rhosfr: fltarr(nimf), rhosfr_err: fltarr(nimf)})
    lit = replicate(lit,nz)
return, lit
end
    
