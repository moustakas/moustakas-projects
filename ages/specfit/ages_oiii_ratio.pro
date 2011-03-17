pro ages_oiii_ratio, ages, postscript=postscript
; jm05aug22uofa - plot the [O III] doublet ratio 
; jm06mar07uofa - updated

    if (n_elements(ages) eq 0L) then ages = read_ages_mz_sample()

    pspath = ages_path(/plots)
    if keyword_set(postscript) then begin
       postthick = 8.0
       dfpsplot, pspath+'oiii_doublet.ps', /color, /square
    endif else begin
       im_window, 0, xratio=0.5, /square
       postthick = 2.0
    endelse

    good = where((ages.oiii_5007[1] gt 0.0) and (ages.oiii_5007_ew[1] gt 0.0) and $
      (ages.oiii_5007[0] gt 0.0) and (ages.oiii_5007_ew[0] gt 0.0) and $
      (ages.oiii_4959_orig[1] gt 0.0))
    ages = ages[good]

    intr = alog10(2.984)

    x = alog10(ages.oiii_5007_ew[0])
;   x = alog10(ages.oiii_4959_ew[0])
;   x = ages.oiii_5007[0]/ages.oiii_5007[1]
    y = alog10(ages.oiii_5007[0]/ages.oiii_4959_orig[0]) - intr

    yrange = [-1.5,+1.5]

;   ages_lineplot, x, y, x*0.0, y*0.0, xsty=3, ysty=3, xthick=postthick, $
;     ythick=postthick, charsize=1.5, charthick=postthick, /bin2d, $
;     yrange=yrange, yminor=5, xtitle='EW([O III] \lambda4959)', $
;     ytitle='log ([O III] \lambda5007 / [O III] \lambda4959)'
    
;   plotsym, 8, 0.5, /fill
    djs_plot, x, y, ps=3, xsty=3, ysty=3, xthick=postthick, ythick=postthick, $
      charsize=1.5, charthick=postthick, xtitle='EW([O III] \lambda5007)', $
      ytitle='log ([O III] \lambda5007 / [O III] \lambda4959)', $
      yrange=yrange, yminor=5;, /xlog
;   oplot, 10^!x.crange, [0,0], line=0, thick=2.0
    oplot, !x.crange, [0,0], line=0, thick=2.0

    res = im_medxbin(x,y,0.1,minpts=5)
    djs_oplot, res.binctr, res.medy, color='red', thick=postthick
    djs_oplot, res.binctr, res.medy+res.sigy, color='red', thick=postthick, line=2
    djs_oplot, res.binctr, res.medy-res.sigy, color='red', thick=postthick, line=2

    niceprint, 10^res.binctr, 100*(10^res.medy-1), 100*(10^res.sigy-1)
    print
    niceprint, interpol(10^res.binctr,100*(10^res.sigy-1),[20.0,10.0,5.0])

    if keyword_set(postscript) then dfpsclose

return
end
    
