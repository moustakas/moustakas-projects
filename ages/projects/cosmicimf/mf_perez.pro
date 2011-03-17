function mf_perez
; Perez-Gonzalez+08 [HDF,CDFS,LockmanHole]; Omega_0=0.3,
; Omega_lamba=0.7, h=0.7; Salpeter IMF
    data = init_mflit(4)
    data.color1      = 'forest green'
    data.color2      = 'tomato'
    data.psym        = 5
    data.symsize     = 3.5
    data.z           = [0.1,0.3,0.5,0.7]
    data.zerr        = 0.1
    data.alpha       = [-1.18,-1.19,-1.22,-1.26]
    data.alpha_err   = [0.12,0.08,0.07,0.08]
    data.mstar       = [11.16,11.20,11.26,11.25]
    data.mstar_err   = [0.25,0.10,0.11,0.08]
    data.phistar     = 10D^[-2.47,-2.65,-2.76,-2.82]
    data.phistar_err = alog(10)*[0.22,0.15,0.13,0.12]*10D^[-2.47,-2.65,-2.76,-2.82]
    data.rho         = 10D^[8.75,8.61,8.57,8.52]
    data.rho_err     = alog(10)*[0.12,0.06,0.04,0.05]*10D^[8.75,8.61,8.57,8.52]
return, data
end

