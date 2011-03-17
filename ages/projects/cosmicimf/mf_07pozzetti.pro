function mf_07pozzetti
; Pozzetti+07 [VVDS]; Omega_0=0.3, Omega_lamba=0.7, h=0.7; Chabrier+03 
; IMF 
    data = init_mflit(4)
    data.color1      = 'forest green'
    data.color2      = 'forest green'
    data.psym        = 4
    data.symsize     = 3.8
    data.z           = [0.225,0.55,0.8,1.05] ; [0.27,0.58,0.8,1.05]
    data.zerr        = [0.175,0.15,0.1,0.15] ; [0.15,0.1,0.1,0.15] ; these are approximate!
; average all the values in Table 2!
    data.alpha       = [$
      mean([-1.26,-1.28,-1.38,-1.39]),$
      mean([-1.23,-1.22,-1.14,-1.16]),$
      mean([-1.23,-1.04,-1.01,-1.16]),$
      mean([-1.09,-1.16,-1.10,-1.20])]
    data.alpha_err       = [$
      stddev([-1.26,-1.28,-1.38,-1.39]),$
      stddev([-1.23,-1.22,-1.14,-1.16]),$
      stddev([-1.23,-1.04,-1.01,-1.16]),$
      stddev([-1.09,-1.16,-1.10,-1.20])]
    data.mstar       = [$
      mean([11.00,11.15,10.93,11.12]),$
      mean([11.00,11.15,10.93,11.12]),$
      mean([10.88,10.83,10.67,10.98]),$
      mean([10.85,10.89,10.78,11.07])]+0.26 ; Chabrier-->Salpeter
    data.mstar_err       = [$
      stddev([11.00,11.15,10.93,11.12]),$
      stddev([11.00,11.15,10.93,11.12]),$
      stddev([10.88,10.83,10.67,10.98]),$
      stddev([10.85,10.89,10.78,11.07])]
    data.phistar       = [$
      mean([1.90,1.75,1.29,1.17]),$
      mean([1.72,1.58,1.83,1.58]),$
      mean([1.60,3.02,2.60,1.74]),$
      mean([1.30,1.80,1.83,1.34])]*1E-3
    data.phistar_err       = [$
      stddev([1.90,1.75,1.29,1.17]),$
      stddev([1.72,1.58,1.83,1.58]),$
      stddev([1.60,3.02,2.60,1.74]),$
      stddev([1.30,1.80,1.83,1.34])]*1E-3
    data.rho     = 10D^([8.45,8.34,8.22,8.14]+0.26) ; Chabrier-->Salpeter
    data.rho_err = alog(10)*[0.09,0.08,0.11,0.12]*10D^([8.45,8.34,8.22,8.14]+0.26)
return, data
end

