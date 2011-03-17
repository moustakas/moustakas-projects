function mf_bell
; Bell+03 [SDSS+2MASS]; Omega_0=0.3, Omega_lamba=0.7, h=1;
; diet-Salpeter IMF
    data = init_mflit(1)
    data.color1      = 'orange'
    data.color2      = 'navy'
    data.psym        = 4 ; 14
    data.symsize     = 3.2
    data.z           = 0.1
    data.zerr        = 0.0
    data.alpha       = -1.10
    data.alpha_err   = 0.02
    data.mstar       = 10.7-alog10(0.7^2)+0.14 ; h=1-->0.7, diet-Salpeter-->Salpeter
    data.mstar_err   = sqrt(0.02^2+(0.05/alog(10))^2) ; include 5% systematic error
i think 0.7 should be divided, not multiplied
    data.phistar     = 0.0102*0.7^3 ; h=1-->0.7
    data.phistar_err = sqrt(0.0005^2 + (0.0102*0.1)^2) ; include 10% systematic
    data.rho         = 5.47E8*0.7*10D^0.14 ; h=1-->0.7, diet-Salpeter-->Salpeter
    data.rho_err     = sqrt(0.11E8^2 + (0.3*5.47E8)^2) ; include 30% systematic error
;   data.rho         = 5.47E8*0.7*10D^0.14 ; h=1-->0.7, diet-Salpeter-->Salpeter
;   data.rho_err     = sqrt((0.11E8/5.47E8/alog(10))^2+(0.3/alog(10))^2) ; include 30% systematic error
return, data
end

