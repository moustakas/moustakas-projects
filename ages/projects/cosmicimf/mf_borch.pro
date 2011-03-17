function mf_borch
; Borch+06 [COMBO-17]; Omega_0=0.3, Omega_lamba=0.7, h=0.7; Kroupa+93 IMF
    data = init_mflit(4)
    data.color1      = 'dodger blue'
    data.color2      = 'dodger blue'
    data.psym        = 9 ; 16
    data.symsize     = 3.2 ; 2.8
    data.z           = [0.3,0.5,0.7,0.9]
    data.zerr        = 0.1
    data.alpha       = -1.1
    data.alpha_err   = 0.02
    data.mstar       = [11.03,11.02,11.09,11.08]+0.14 ; Kroupa-->Salpeter
    data.mstar_err   = [0.08,0.08,0.15,0.11]
    data.phistar     = [19,18,16,12]*1E-4
    data.phistar_err = [9,6,5,4]*1E-4
    data.rho         = 10D^([8.34,8.32,8.33,8.17]+0.14) ; Kroupa-->Salpeter
    data.rho_err     = alog(10)*[0.15,0.11,0.10,0.18]*10D^[8.34,8.32,8.33,8.17]
return, data
end

