pro ilbert_mass_function
; jm09apr14nyu - build the mass function for *all* star-forming
;   galaxies (intermediate- and high-activity) from Ilbert+10 
; jm10mar29ucsd - added the quiescent points to form the "total" MF 

    zz = [0.3,0.5,0.7,0.9,1.1,1.35,1.75]
    nzz = n_elements(zz)

; quiescent
    alpha_q = [-0.91,-0.56,-0.25,0.04,0.25,0.5,0.5]
    mstar_q = [11.13,10.97,10.83,10.77,10.70,10.64,10.67]
    phistar_q = 1D-3*[1.12,0.87,1.15,1.43,0.55,0.26,0.10]
    
; intermediate activity    
    alpha_inter = [-1.26,-1.09,-0.94,-0.59,-0.37,-0.54,-0.43]
    mstar_inter = [11.01,10.95,10.89,10.78,10.73,10.76,10.86]
    phistar_inter = 1D-3*[1.08,0.78,0.82,1.28,1.10,0.65,0.30]

; high-activity    
    alpha_high = [-1.28,-1.29,-1.30,-1.21,-1.16,-1.17,-1.30]
    mstar_high = [10.26,10.30,10.39,10.42,10.42,10.49,10.86]
    phistar_high = 1D-3*[0.62,0.77,1.04,1.40,1.36,1.13,0.30]

    mass = im_array(9.2,12.5,0.02)
    nmass = n_elements(mass)
    quiescent = fltarr(nmass,nzz)
    inter = fltarr(nmass,nzz)
    high = fltarr(nmass,nzz)

    for jj = 0, nzz-1 do quiescent[*,jj] = mf_schechter(mass,$
      phistar_q[jj],mstar_q[jj],alpha_q[jj])
    for jj = 0, nzz-1 do inter[*,jj] = mf_schechter(mass,$
      phistar_inter[jj],mstar_inter[jj],alpha_inter[jj])
    for jj = 0, nzz-1 do high[*,jj] = mf_schechter(mass,$
      phistar_high[jj],mstar_high[jj],alpha_high[jj])

; all galaxies
    for jj = 0, nzz-1 do begin
       totmf = quiescent[*,jj]+inter[*,jj]+high[*,jj]
       mf_fit_schechter, mass, totmf, 0.1*totmf, ss1, /quiet
       if (jj eq 0) then ss = ss1 else ss = [ss,ss1]
       djs_plot, mass, totmf, /ylog, xsty=3, ysty=3, psym=6, sym=0.2
       djs_oplot, mass, quiescent, color='red'
       djs_oplot, mass, inter, color='green'
       djs_oplot, mass, high, color='blue'
       djs_oplot, mass, mf_schechter(mass,ss1), $
         color='orange', thick=2, line=5
       cc = get_kbrd(1)
    endfor
    struct_print, ss

stop    
    
; all star-forming galaxies    
    for jj = 0, nzz-1 do begin
       mf_fit_schechter, mass, inter[*,jj]+high[*,jj],$
         0.1*(inter[*,jj]+high[*,jj]), ss1
       if (jj eq 0) then ss = ss1 else ss = [ss,ss1]
       djs_plot, alog10(mass), inter[*,jj]+high[*,jj], $
         /ylog, xsty=3, ysty=3
       djs_oplot, alog10(mass), mf_schechter(mass,ss1), $
         color='red', thick=2
;      cc = get_kbrd(1)
    endfor
    struct_print, ss

    niceprint, zz, alog10(ss.mstar)


stop    
    
return
end
    
