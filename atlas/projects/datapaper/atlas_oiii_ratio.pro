pro atlas_oiii_ratio, atlas, postscript=postscript
; jm04apr27uofa
; plot the [O III] doublet ratio

    if (n_elements(atlas) eq 0L) then atlas = read_integrated(/silent)

    pspath = atlas_path(/papers)+'atlas/FIG_ATLAS/'
    if keyword_set(postscript) then begin
       postthick = 8.0
       dfpsplot, pspath+'oiii_doublet.ps', /color, /square
    endif else begin
       im_window, 0, xratio=0.5, /square
       postthick = 2.0
    endelse

    good = where((atlas.oiii_5007[1] gt 0.0) and (atlas.oiii_4959_orig[1] gt 0.0))
    atlas = atlas[good]

    intr = alog10(2.984)

    x = alog10(atlas.oiii_4959_ew[0])
;   x = atlas.oiii_5007[0]/atlas.oiii_5007[1]
    y = alog10(atlas.oiii_5007[0]/atlas.oiii_4959_orig[0]) - intr

    yrange = [-0.9,+0.9]

    plotsym, 8, 1.5;, /fill
    djs_plot, x, y, ps=8, xsty=3, ysty=3, xthick=postthick, ythick=postthick, $
      charsize=1.5, charthick=postthick, xtitle='EW([O III] \lambda4959)', $
      ytitle='log ([O III] \lambda5007 / [O III] \lambda4959)', $
      yrange=yrange, yminor=5;, /xlog
    oplot, 10^!x.crange, [0,0], line=0, thick=2.0
;   djs_oplot, [5,5], !y.crange, line=1, color='navy', thick=3.0

    res = im_medxbin(x,y,0.3,minpts=5)
    djs_oplot, res.binctr, res.medy, color='red', thick=postthick
    djs_oplot, res.binctr, res.medy+res.sigy, color='red', thick=postthick, line=2
    djs_oplot, res.binctr, res.medy-res.sigy, color='red', thick=postthick, line=2

    niceprint, 10^res.binctr, 100*(10^res.medy-1), 100*(10^res.sigy-1)
    
stop    
    
    if keyword_set(postscript) then dfpsclose

return
end
    
