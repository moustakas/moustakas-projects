pro order_oplot, xx, yy, allorder
; plot each order using different colors

    names = fsc_color(/all,/name)
    keep = where($
      (strmatch(names,'*light*',/fold) eq 0) and $
      (strmatch(names,'*white*',/fold) eq 0) and $
      (strmatch(names,'*snow*',/fold) eq 0) and $
      (strmatch(names,'*linen*',/fold) eq 0) and $
      (strmatch(names,'*beige*',/fold) eq 0) and $
      (strmatch(names,'*cornsilk*',/fold) eq 0) and $
      (strmatch(names,'*seashell*',/fold) eq 0) and $
      (strmatch(names,'*almond*',/fold) eq 0) and $
      (strmatch(names,'*papaya*',/fold) eq 0) and $
      (strmatch(names,'*bisque*',/fold) eq 0) and $
      (strmatch(names,'*moccasin*',/fold) eq 0) and $
      (strmatch(names,'*ivory*',/fold) eq 0))
    allcolor = names[keep]
    
    order = allorder[uniq(allorder,sort(allorder))]
    norder = n_elements(order)
    for ii = 0, norder-1 do begin
       these = where(order[ii] eq allorder)
       djs_oplot, xx[these], yy[these], psym=symcat(16), $
         symsize=0.7, color=fsc_color(allcolor[ii],100)
;      print, allcolor[ii]
    endfor

return
end    

pro qaplot_alpha_arcfit, arcfile, qafile=qafile
; jm10jan04ucsd - make QAplots

; read the output from GET_ALPHA_ARCFIT
    splog, 'Reading '+arcfile
    info = mrdfits(arcfile,1)

    statsfile = repstr(arcfile,'.fits','_stats.fits')
    splog, 'Reading '+statsfile
    stats = mrdfits(statsfile,1)

    allarc = info.night+info.arc
    narc = n_elements(uniq(allarc,sort(allarc)))
    
; finally make the plots
    light = 2.99792458D5 ; speed of light [km/s]
    symsize = 0.7

    im_plotconfig, 0, pos, psfile=qafile, xmargin=[1.6,0.3], $
      width=6.6, height=6.6
    
    djs_plot, [0], [0], /nodata, position=pos, /ylog, $
      xsty=3, ysty=3, xrange=minmax(stats.wave), yrange=[0.1,1E3], $
      xtitle='Wavelength (\AA)', ytitle='Dispersion (m\AA)'
    order_oplot, stats.wave, 1E3*stats.sigma, stats.order
    
    djs_plot, [0], [0], /nodata, position=pos, /ylog, $
      xsty=3, ysty=3, xrange=minmax(stats.wave), yrange=[0.02,40], $
      xtitle='Wavelength (\AA)', ytitle='Fractional Dispersion (10^{-6})'
    order_oplot, stats.wave, 1E6*stats.sigma/stats.wave, stats.order
    
    djs_plot, [0], [0], /nodata, position=pos, $
      xsty=3, ysty=3, xrange=minmax(stats.wave), yrange=[-30.0,30.0], $
      xtitle='Wavelength (\AA)', ytitle='Mean Offset (m\AA)'
    order_oplot, stats.wave, 1E3*stats.mean, stats.order

    djs_plot, [0], [0], /nodata, position=pos, /ylog, $
      xsty=3, ysty=3, xrange=minmax(stats.wave), yrange=[2E-3,500.0], $
      xtitle='Wavelength (\AA)', ytitle='Mean Offset (absolute value, m\AA)'
    order_oplot, stats.wave, abs(1E3*stats.mean), stats.order

    djs_plot, [0], [0], /nodata, position=pos, $
      xsty=3, ysty=3, xrange=minmax(stats.wave), yrange=[2E-4,10.0], /ylog, $
      xtitle='Wavelength (\AA)', ytitle='Mean Offset (absolute value, km s^{-1})'
    order_oplot, stats.wave, abs(light*stats.mean/stats.wave), $
      stats.order

    djs_plot, [0], [0], /nodata, position=pos, $
      xsty=3, ysty=3, xrange=minmax(stats.wave), yrange=[0.0,1.0], $
      xtitle='Wavelength (\AA)', ytitle='Fraction of Time Line Rejected'
    order_oplot, stats.wave, stats.nbad/float(stats.nused), $
      stats.order

    im_plotconfig, psfile=qafile, /psclose, /gzip
    
return
end
