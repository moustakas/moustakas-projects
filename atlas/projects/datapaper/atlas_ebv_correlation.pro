pro atlas_ebv_correlation, atlas, atlasnodust, nuclear, $
  nuclearnodust, postscript=postscript
; jm05aug03uofa - correlation between the continuum and the nebular
;                 reddening  

    if (n_elements(atlas) eq 0L) then $
      atlas = read_integrated(atlasnodust=atlasnodust)
    if (n_elements(nuclear) eq 0L) then $
      nuclear = read_nuclear(atlasnodust=nuclearnodust)

    pspath = atlas_path(/papers)+'atlas/FIG_ATLAS/'
    if keyword_set(postscript) then begin
       postthick = 8.0
       dfpsplot, pspath+'ebv_correlation.eps', /color, xsize=8.5, ysize=7.5, /encapsulated;/square
    endif else begin
       im_window, 0, xratio=0.5, /square
       postthick = 2.0
    endelse

    good = where((atlasnodust.ebv_hahb_err gt 0.0),ngood)
    good_nuclear = where((nuclearnodust.ebv_hahb_err gt 0.0),ngood)

    x = atlasnodust[good].ebv_hahb
    y = atlasnodust[good].continuum_ebv[0]

    xnuc = nuclearnodust[good_nuclear].ebv_hahb
    ynuc = nuclearnodust[good_nuclear].continuum_ebv[0]

    rcor = r_correlate(x,y,zd=zd,probd=probd)
    splog, 'Integrated Spectra: ', rcor, probd, zd

    rcor = r_correlate([x,xnuc],[y,ynuc],zd=zd,probd=probd)
    splog, 'Integrated+Nuclear Spectra: ', rcor, probd, zd

    intstats = im_stats(y/x,/verbose,sigrej=3.0)
    nucstats = im_stats(ynuc/xnuc,/verbose,sigrej=3.0)
    
    xrange = [-0.05,1.6]
    yrange = [-0.03,0.7]

    xebv = findgen(18.0)/10.0
    
    plotsym, 8, 1.0, /fill
    djs_plot, x, y, ps=8, xsty=3, ysty=3, xthick=postthick, ythick=postthick, $
      charsize=2.0, charthick=postthick, xtitle='Nebular E(B-V) [mag]', $
      ytitle='Continuum E(B-V) [mag]', yrange=yrange, xrange=xrange

    plotsym, 0, 1.0, thick=postthick;, /fill
    djs_oplot, xnuc, ynuc, ps=8

    oplot, xebv, 0.5*xebv, line=2, thick=postthick
    oplot, xebv, 0.25*xebv, line=2, thick=postthick
    oplot, xebv, 0.1*xebv, line=2, thick=postthick

    xyouts, 0.93, 0.55, textoidl('\alpha = 0.5'), /data, charsize=1.5, $
      charthick=postthick, align=0.5
    xyouts, 1.3, 0.39, textoidl('\alpha = 0.25'), /data, charsize=1.5, charthick=postthick, align=0.5
    xyouts, 1.4, 0.17, textoidl('\alpha = 0.1'), /data, charsize=1.5, charthick=postthick, align=0.5

    if keyword_set(postscript) then dfpsclose

stop    
    
return
end
    
