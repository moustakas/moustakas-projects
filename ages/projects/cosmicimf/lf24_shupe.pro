function lf24_shupe, nimf=nimf
; Shupe+98 (from Rujopakarn+10)
    data = init_lf24lit(1,nimf=nimf)
    data.color1      = 'forest green'
    data.color2      = 'blue'
    data.psym        = 11
    data.symsize     = 3.5
    data.z           = 0.0
    data.zerr        = 0.0
    data.alpha       = 0.437
    data.alpha_err   = 0.032
    data.beta        = 1.749
    data.beta_err    = 0.067
    data.lstar       = alog10(4.47D9)
    data.lstar_err   = 0.4/4.47/alog(10)
    data.phistar     = 1.0D-3
    data.phistar_err = 0.6D-3
return, data
end

