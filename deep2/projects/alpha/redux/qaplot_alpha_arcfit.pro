pro order_oplot, xx, yy, allorder, symsize=symsize, psym=psym, color=color
; plot each order using different colors

    order = allorder[uniq(allorder,sort(allorder))]
    norder = n_elements(order)
    
    if (n_elements(symsize) eq 0) then symsize = 1.0
    if (n_elements(psym) eq 0) then psym = 16
    if (n_elements(color) eq 0) then begin
       names = fsc_color(/all,/name)
       keep = where($
         (strmatch(names,'*light*',/fold) eq 0) and $
         (strmatch(names,'*white*',/fold) eq 0) and $
         (strmatch(names,'*snow*',/fold) eq 0) and $
         (strmatch(names,'*linen*',/fold) eq 0) and $
         (strmatch(names,'*beige*',/fold) eq 0) and $
         (strmatch(names,'*ivory*',/fold) eq 0) and $
         (strmatch(names,'*rose*',/fold) eq 0) and $
         (strmatch(names,'*thistle*',/fold) eq 0) and $
         (strmatch(names,'*wheat*',/fold) eq 0) and $
         (strmatch(names,'*honey*',/fold) eq 0) and $
         (strmatch(names,'*lavender*',/fold) eq 0) and $
         (strmatch(names,'*cornsilk*',/fold) eq 0) and $
         (strmatch(names,'*burlywood*',/fold) eq 0) and $
         (strmatch(names,'*medium*',/fold) eq 0) and $
         (strmatch(names,'*seashell*',/fold) eq 0) and $
         (strmatch(names,'*almond*',/fold) eq 0) and $
         (strmatch(names,'*papaya*',/fold) eq 0) and $
         (strmatch(names,'*bisque*',/fold) eq 0) and $
         (strmatch(names,'*moccasin*',/fold) eq 0),nkeep)
       allcolor = names[keep]
    endif else allcolor = replicate(color,norder)
    
    for ii = 0, norder-1 do begin
       these = where(order[ii] eq allorder,nthese)
       if (nthese eq 1) then begin
          plots, xx[these], yy[these], psym=symcat(psym,thick=4), $
            symsize=symsize, color=fsc_color(allcolor[ii],101)
       endif else begin
          djs_oplot, xx[these], yy[these], psym=symcat(psym,thick=4), $
            symsize=symsize, color=fsc_color(allcolor[ii],101)
       endelse
;      print, allcolor[ii]
    endfor

return
end    

pro qaplot_alpha_arcfit, arcfile, iter=iter, rejcut=rejcut, qafile=qafile
; jm10jan04ucsd - make QAplots

    title = 'Iteration '+string(iter,format='(I2.2)')
    
; read the output from GET_ALPHA_ARCFIT
    splog, 'Reading '+arcfile
    info = mrdfits(arcfile,1)

    statsfile = repstr(arcfile,'.fits','_stats.fits')
    splog, 'Reading '+statsfile
    stats = mrdfits(statsfile,1)
    
    allarc = info.night+info.arc
    narc = n_elements(uniq(allarc,sort(allarc)))
    
    frac = stats.nbad/float(stats.nused)
    rej = where(frac gt rejcut,nrej,ncomp=ngood) ; rejected lines

; finally make the plots
    light = 2.99792458D5 ; speed of light [km/s]

    im_plotconfig, 0, pos, psfile=qafile, xmargin=[1.8,0.3], $
      width=6.4, height=5.8, ymargin=[0.6,1.1]
    
    djs_plot, [0], [0], /nodata, position=pos, /ylog, $
      xsty=3, ysty=1, xrange=minmax(stats.wave), yrange=[0.1,500], $
      xtitle='Wavelength (\AA)', ytitle='Dispersion (m\AA)', title=title
    im_legend, ['N_{good}='+strtrim(ngood,2),'N_{rej}='+strtrim(nrej,2)], $
      /left, /top, box=0, charsize=1.5
    order_oplot, stats.wave, 1E3*stats.sigma, stats.order
    if (nrej ne 0) then order_oplot, stats[rej].wave, 1E3*stats[rej].sigma, $
      stats[rej].order, psym=6, symsize=1.2, color='red'
    
    djs_plot, [0], [0], /nodata, position=pos, /ylog, $
      xsty=3, ysty=1, xrange=minmax(stats.wave), yrange=[0.05,50], $
      xtitle='Wavelength (\AA)', ytitle='Fractional Dispersion (10^{-6})', title=title
    im_legend, ['N_{good}='+strtrim(ngood,2),'N_{rej}='+strtrim(nrej,2)], $
      /left, /top, box=0, charsize=1.5
    order_oplot, stats.wave, 1E6*stats.sigma/stats.wave, stats.order
    if (nrej ne 0) then order_oplot, stats[rej].wave, 1E6*stats[rej].sigma/stats[rej].wave, $
      stats[rej].order, psym=6, symsize=1.2, color='red'
    
    djs_plot, [0], [0], /nodata, position=pos, $
      xsty=3, ysty=1, xrange=minmax(stats.wave), yrange=40*[-1,1], $
      xtitle='Wavelength (\AA)', ytitle='Mean Offset (m\AA)', title=title
    im_legend, ['N_{good}='+strtrim(ngood,2),'N_{rej}='+strtrim(nrej,2)], $
      /left, /top, box=0, charsize=1.5
    order_oplot, stats.wave, 1E3*stats.mean, stats.order
    if (nrej ne 0) then order_oplot, stats[rej].wave, 1E3*stats[rej].mean, $
      stats[rej].order, psym=6, symsize=1.2, color='red'

    djs_plot, [0], [0], /nodata, position=pos, /ylog, $
      xsty=3, ysty=1, xrange=minmax(stats.wave), yrange=[5E-3,1.3E3], $
      xtitle='Wavelength (\AA)', ytitle='Mean Offset (absolute value, m\AA)', title=title
    im_legend, ['N_{good}='+strtrim(ngood,2),'N_{rej}='+strtrim(nrej,2)], $
      /left, /top, box=0, charsize=1.5
    order_oplot, stats.wave, abs(1E3*stats.mean), stats.order
    if (nrej ne 0) then order_oplot, stats[rej].wave, abs(1E3*stats[rej].mean), $
      stats[rej].order, psym=6, symsize=1.2, color='red'

    djs_plot, [0], [0], /nodata, position=pos, $
      xsty=3, ysty=1, xrange=minmax(stats.wave), yrange=[3E-4,90.0], /ylog, $
      xtitle='Wavelength (\AA)', ytitle='Mean Offset (absolute value, km s^{-1})', title=title
    im_legend, ['N_{good}='+strtrim(ngood,2),'N_{rej}='+strtrim(nrej,2)], $
      /left, /top, box=0, charsize=1.5
    order_oplot, stats.wave, abs(light*stats.mean/stats.wave), stats.order
    if (nrej ne 0) then order_oplot, stats[rej].wave, abs(light*stats[rej].mean/stats[rej].wave), $
      stats[rej].order, psym=6, symsize=1.2, color='red'

    djs_plot, [0], [0], /nodata, position=pos, /ylog, $
      xsty=3, ysty=3, xrange=minmax(stats.wave), yrange=[0.8,max(stats.nused)*2], $
      xtitle='Wavelength (\AA)', ytitle='Number of Times Line Used', title=title
    im_legend, ['N_{good}='+strtrim(ngood,2),'N_{rej}='+strtrim(nrej,2)], $
      /left, /top, box=0, charsize=1.5
    order_oplot, stats.wave, stats.nused, stats.order
    if (nrej ne 0) then order_oplot, stats[rej].wave, stats[rej].nused, $
      stats[rej].order, psym=6, symsize=1.2, color='red'

    djs_plot, [0], [0], /nodata, position=pos, $
      xsty=3, ysty=1, xrange=minmax(stats.wave), yrange=[-0.05,1.05], $
      xtitle='Wavelength (\AA)', ytitle='Fraction of Time Line Rejected', title=title
    im_legend, ['N_{good}='+strtrim(ngood,2),'N_{rej}='+strtrim(nrej,2)], $
      /left, /top, box=0, charsize=1.5
    order_oplot, stats.wave, frac, stats.order
    if (nrej ne 0) then order_oplot, stats[rej].wave, frac[rej], $
      stats[rej].order, psym=6, symsize=1.2, color='red'
    djs_oplot, !x.crange, rejcut*[1,1], line=0, thick=3

    im_plotconfig, psfile=qafile, /psclose, /gzip
    
return
end
