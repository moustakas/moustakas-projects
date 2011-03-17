function lf24_wiphu, nimf=nimf
; Rujopakarn+10
    data = init_lf24lit(7,nimf=nimf)
    data.color1      = 'cyan'
    data.color2      = 'firebrick'
    data.psym        = 6
    data.symsize     = 3.5
    data.z = [0.0,mean([0.05,0.2]),mean([0.2,0.35]),mean([0.35,0.5]),$
      mean([0.5,0.65]),mean([0.7,1.0]),mean([1.0,1.3])]
    data.zerr = [0.0,stddev([0.05,0.2]),stddev([0.2,0.35]),stddev([0.35,0.5]),$
      stddev([0.5,0.65]),stddev([0.7,1.0]),stddev([1.0,1.3])]
    data.alpha       = 0.37
    data.alpha_err   = 0.04
    data.beta        = 2.36
    data.beta_err    = 0.41
    data.lstar       = alog10([4.27,7.25,10.6,15.5,23.6,30.3,62.4]*1D9)
    data.lstar_err   = [0.71,1.12,1.8,2.4,3.6,7.8,18.8]/[4.27,7.25,10.6,15.5,23.6,30.3,62.4]/alog(10)
    data.phistar     = 1.2D-3
    data.phistar_err = 0.8D-3
; uncertainty due to large-scale structure on larger than ~1 deg scales
    data.phistar_lss_err = data.phistar*[0.0,0.15,0.11,0.1,0.09,0.0,0.0]
return, data
end

