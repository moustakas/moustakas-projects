function mf_gallazzi
; Gallazzi+08
    data = init_mflit(1)
;   data.color1      = 'orchid'
;   data.color2      = 'orange'
;   data.psym        = 11 ; 34
;   data.symsize     = 3.5
;   data.z           = 0.05
;   data.zerr        = 0.0
;   data.alpha       = -1.18
;   data.alpha_err   = 0.034
;   data.mstar       = alog10(7.07D10/0.7^2) ; h=1-->0.7
;   data.mstar_err   = 0.18/7.07/alog(10)
;   data.phistar     = 0.009*0.7^3   ; h=1-->0.7
;   data.phistar_err = 0.0014
    data.rho         = 3.413D8*10^0.26 ; Chabrier --> Salpeter
    data.rho_err     = 0.562D8*10^0.26
return, data
end

