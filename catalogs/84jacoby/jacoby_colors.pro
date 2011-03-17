function jacoby_colors, doplot=doplot
; jm03sep4uofa
; measure synthetic (spectral) colors for the jacoby atlas

    forward_function readfilt, filtmag
    
; read the atlas

    jhc = read_jacoby_atlas()
    njhc = n_elements(jhc)
    
    result = {$
      star:      '', $
      sp_type:   '', $
      U_B_JHC:  0.0, $
      B_V_JHC:  0.0, $
      U_B:      0.0, $
      B_V:      0.0, $
      B_R:      0.0}
    result = replicate(result,njhc)

    result.star = jhc.star
    result.sp_type = jhc.sp_type
    result.U_B_JHC = jhc.U_B
    result.B_V_JHC = jhc.B_V
    
; load the filters
    
    loadnmz

    u = readfilt('JohnU.dat')
    b = readfilt('JohnB.dat')
    v = readfilt('JohnV.dat')
    r = readfilt('HARRR.dat')

    for i = 0L, njhc-1L do begin

       print, format='("Computing colors for star ",I3,"/",I3,".",A1,$)', i+1, njhc, string(13b)

       wave = jhc[i].wave
       flux = jhc[i].flux
       sed = create_struct('wave', wave, 'flux', flux)

       result[i].U_B = filtmag(sed,u,/vega)-filtmag(sed,b,/vega)
       result[i].B_V = filtmag(sed,b,/vega)-filtmag(sed,v,/vega)
       result[i].B_R = filtmag(sed,b,/vega)-filtmag(sed,r,/vega)

       if keyword_set(doplot) then begin
       
          plot, sed.wave, sed.flux/(0.5*max(sed.flux)), xsty=3, ysty=3
          oplot, u.filtw, u.filtf, line=0
          oplot, b.filtw, b.filtf, line=0
          oplot, v.filtw, v.filtf, line=2
          oplot, r.filtw, r.filtf, line=3
          cc = get_kbrd(1)
          
       endif
       
    endfor
    print

; generate postscript output

    x = result.B_V
    y = result.B_V_JHC
    xtitle = '(B-V)  [Synthetic]'
    ytitle = '(B-V)  [Jacoby]'

    xrange = [-1.0,2.0]
    yrange = xrange

    plotsym, 0, 1.3
    djs_plot, x, y, ps=8, xsty=3, ysty=3, xtitle=xtitle, ytitle=ytitle, $
      charsize=2.0, charthick=2.0, xthick=2.0, ythick=2.0, $
      xrange=xrange, yrange=yrange
    oplot, [-10,10], [-10,10], line=0, thick=2.0

return, result
end    
