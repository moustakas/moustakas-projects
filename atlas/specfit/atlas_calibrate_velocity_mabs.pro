pro atlas_calibrate_velocity_mabs, postscript=postscript, encapsulated=encapsulated
; jm05may16uofa
; calibrate a relationship between maximum rotational velocity and M_B
; using average values in Hyperleda; remember, velocity dispersion is
; v_c/sqrt(2) 

    datapath = atlas_path(/ned)
    data = rsex(datapath+'leda_velocity_mabs.dat')

    im_openclose, datapath+'atlas2d_calibrate_velocity_mabs', $
      postscript=postscript, encapsulated=encapsulated, xsize=8.5, ysize=8.5
    
    if keyword_set(postscript) then begin
       postthick = 8.0
    endif else begin
       im_window, 0, xratio=0.4, /square
       postthick = 2.0
    endelse

; maximum rotational velocity    
    
    pagemaker, nx=1, ny=1, xspace=0, yspace=0.0, width=6.8, height=6.8, $
      xmargin=[1.3,0.4], ymargin=[0.4,1.3], xpage=8.5, ypage=8.5, $
      position=pos, /normal

    djs_plot, data.mb, data.logavmm/sqrt(2), ps=3, xsty=3, ysty=3, xrange=[-12,-25], $
      yrange=[0.5,2.3], xthick=postthick, ythick=postthick, charsize=2.0, $
      charthick=postthick, xtitle='M_{B}', ytitle='log \sigma (km s^{-1})', $
      position=pos

    fit = where((data.mb lt -16.0) and (data.mb gt -23.0))
;   sixlin, data.mb, data.logavmm, a, siga, b, sigb
;   coeff = [a[2],b[2]]
    coeff = linfit(data[fit].mb,data[fit].logavmm/sqrt(2))

    MBaxis = findgen(((max(data[fit].mb))-(min(data[fit].mb)))/0.01+1)*0.01+min(data[fit].mb)
    yfit = poly(MBaxis,coeff)
    djs_oplot, MBaxis, yfit, thick=postthick, line=0, color='red'

    label = 'log \sigma = '+strtrim(string(coeff[0],format='(F12.3)'),2)+''+$
      strtrim(string(coeff[1],format='(F12.3)'),2)+' M_{B}'
    legend, textoidl(label), /left, /top, box=0, charsize=1.8, charthick=postthick
    
; central velocity dispersion

;   good = where(data.logs gt -900)
;   djs_plot, data[good].mb, data[good].logs, ps=3, xsty=3, ysty=3, xrange=[-11,-24], $
;     yrange=[0.1,4.0], xthick=postthick, ythick=postthick, charsize=2.0, $
;     charthick=postthick, xtitle='M_{B}', ytitle='log \sigma (km s^{-1})'
;
;   sixlin, data.mb, data.logavmm, a, siga, b, sigb
;   coeff = [a[2],b[2]]
;
;   yfit = poly(data.mb,coeff)
;   djs_oplot, data.mb, yfit, thick=postthick, line=0, color='dark green'

    if keyword_set(postscript) then dfpsclose

stop    
    
return
end
    
