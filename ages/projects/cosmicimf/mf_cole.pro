function mf_cole
; Cole+01 [SDSS+2MASS]; Omega_0=0.3, Omega_lamba=0.7, h=0.7;
; Salpeter IMF
    data = init_mflit(1)
    data.color1      = 'orchid'
    data.color2      = 'orange'
    data.psym        = 11 ; 34
    data.symsize     = 3.5
    data.z           = 0.05
    data.zerr        = 0.0
    data.alpha       = -1.18
    data.alpha_err   = 0.034
    data.mstar       = alog10(7.07D10)
    data.mstar_err   = 0.18/7.07/alog(10)
    data.phistar     = 0.009
    data.phistar_err = 0.0014
    data.rho         = 1.36D11*0.0029 ; Omega*-->rho*
    data.rho_err     = 1.36D11*0.00039
;   data.rho         = alog10(1.36D11*0.0029/0.7) ; Omega*-->rho*
;   data.rho_err     = 0.00039/0.0026/alog(10)
return, data
end

