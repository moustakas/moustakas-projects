function init_mflit, nz
    lit = {$
      color1: '', color2: '', psym: 0.0, symsize: 0.0, $
      z:       0.0, zerr:       0.0, $
      alpha:   0.0, alpha_err:   0.0, $
      mstar:   0.0, mstar_err:   0.0, $
      phistar: 0.0, phistar_err: 0.0, $
      rho:     0.0, rho_err: 0.0}
    lit = replicate(lit,nz)
return, lit
end
