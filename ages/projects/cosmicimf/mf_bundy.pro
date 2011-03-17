function mf_bundy
; Bundy+06 [DEEP2]; Omega_0=0.3, Omega_lamba=0.7, h=0.7; Chabrier+03
; IMF 
    data = init_mflit(3)
    data.color1      = 'red'
    data.color2      = 'red'
    data.psym        = 5
    data.symsize     = 3.8
    data.z           = [0.55,0.875,1.2]
    data.zerr        = [0.15,0.125,0.2]
    data.alpha       = [-0.81,-0.59,-0.51]
    data.alpha_err   = [0.2,0.4,0.5]
    data.mstar       = [10.94,10.87,10.97]+0.26 ; Chabrier-->Salpeter
    data.mstar_err   = [0.08,0.10,0.25]
    data.phistar     = [0.0027,0.0031,0.0012]
    data.phistar_err = [0.0004,0.0010,0.0008]
    data.rho         = 10D^([8.31,8.30,8.15]+0.26) ; Chabrier-->Salpeter
    data.rho_err     = alog(10)*[0.07,0.1,0.1]*10D^[8.31,8.30,8.15]
return, data
end

