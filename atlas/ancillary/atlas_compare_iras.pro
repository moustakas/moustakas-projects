pro atlas_compare_iras, postscript=postscript
; jm05nov08uofa - compare the three IRAS flux measurements

    datapath = atlas_path(/analysis)
    photoname = 'atlas_ned_photo.fits.gz'
    
    photo = mrdfits(datapath+photoname,1,/silent)

; set up some plotting variables    

    plotsym, 0, 1.8, /fill

    if keyword_set(postscript) then begin
       postthick = 8.0
    endif else begin
       im_window, 0, xratio=0.5, /square
       postthick = 2.0
    endelse

    pagemaker, nx=1, ny=2, yspace=0, xspace=0, xmargin=[1.5,0.2], $
      ymargin=[0.5,1.5], position=pos, height=[5.0,2.5], /normal, $
      xpage=8.5, ypage=9.0, width=6.8

    charsize = 2.0

; ---------------------------------------------------------------------------    
; SOIFER vs MOSHIR    
; ---------------------------------------------------------------------------    

; 12 micron

    psname = 'soifer_vs_moshir_12'
    im_openclose, psname, xsize=8.5, ysize=9.0, postscript=postscript

    xtitle = 'Soifer et al. (1989)'
    ytitle = 'Moshir et al. (1990)'
    title = '12 micron'
    
    good = where((photo.soifer_iras_12 gt 0.0) and (photo.moshir_iras_12 gt 0.0))

    x = photo[good].soifer_iras_12
    y = photo[good].moshir_iras_12
    resid = 100.0*(y/x-1)

    print, 'Soifer vs Moshir - 12'
    niceprint, photo[good].galaxy, x, y, resid

    xrange = [min(x)<min(y),max(x)>max(y)]
    yrange = xrange
    yrange2 = max(abs(resid))*[-1.1,1.1]
    
    djs_plot, x, y, ps=8, xsty=3, ysty=3, xrange=xrange, yrange=yrange, $
      xthick=postthick, ythick=postthick, charsize=charsize, charthick=postthick, $
      position=pos[*,0], xtickname=replicate(' ',10), xtitle='', $
      ytitle=ytitle, /xlog, /ylog, title=title
    djs_oplot, 10^!x.crange, 10^!y.crange, line=0, thick=2.0

    djs_plot, x, resid, ps=8, xsty=3, ysty=3, xrange=xrange, yrange=yrange2, $
      xthick=postthick, ythick=postthick, charsize=charsize, charthick=postthick, $
      position=pos[*,1], /noerase, xtitle=xtitle, ytitle='Residuals [%]', /xlog
    djs_oplot, 10^!x.crange, [0,0], line=0, thick=2.0

    im_openclose, postscript=postscript, /close

; 25 micron

    psname = 'soifer_vs_moshir_25'
    im_openclose, psname, xsize=8.5, ysize=9.0, postscript=postscript

    xtitle = 'Soifer et al. (1989)'
    ytitle = 'Moshir et al. (1990)'
    title = '25 micron'
    
    good = where((photo.soifer_iras_25 gt 0.0) and (photo.moshir_iras_25 gt 0.0))

    x = photo[good].soifer_iras_25
    y = photo[good].moshir_iras_25
    resid = 100.0*(y/x-1)
    
    print, 'Soifer vs Moshir - 25'
    niceprint, photo[good].galaxy, x, y, resid

    xrange = [min(x)<min(y),max(x)>max(y)]
    yrange = xrange
    yrange2 = max(abs(resid))*[-1.1,1.1]
    
    djs_plot, x, y, ps=8, xsty=3, ysty=3, xrange=xrange, yrange=yrange, $
      xthick=postthick, ythick=postthick, charsize=charsize, charthick=postthick, $
      position=pos[*,0], xtickname=replicate(' ',10), xtitle='', $
      ytitle=ytitle, /xlog, /ylog, title=title
    djs_oplot, 10^!x.crange, 10^!y.crange, line=0, thick=2.0

    djs_plot, x, resid, ps=8, xsty=3, ysty=3, xrange=xrange, yrange=yrange2, $
      xthick=postthick, ythick=postthick, charsize=charsize, charthick=postthick, $
      position=pos[*,1], /noerase, xtitle=xtitle, ytitle='Residuals [%]', /xlog
    djs_oplot, 10^!x.crange, [0,0], line=0, thick=2.0

    im_openclose, postscript=postscript, /close

; 60 micron

    psname = 'soifer_vs_moshir_60'
    im_openclose, psname, xsize=8.5, ysize=9.0, postscript=postscript

    xtitle = 'Soifer et al. (1989)'
    ytitle = 'Moshir et al. (1990)'
    title = '60 micron'
    
    good = where((photo.soifer_iras_60 gt 0.0) and (photo.moshir_iras_60 gt 0.0))

    x = photo[good].soifer_iras_60
    y = photo[good].moshir_iras_60
    resid = 100.0*(y/x-1)
    
    print, 'Soifer vs Moshir - 60'
    niceprint, photo[good].galaxy, x, y, resid

    xrange = [min(x)<min(y),max(x)>max(y)]
    yrange = xrange
    yrange2 = max(abs(resid))*[-1.1,1.1]
    
    djs_plot, x, y, ps=8, xsty=3, ysty=3, xrange=xrange, yrange=yrange, $
      xthick=postthick, ythick=postthick, charsize=charsize, charthick=postthick, $
      position=pos[*,0], xtickname=replicate(' ',10), xtitle='', $
      ytitle=ytitle, /xlog, /ylog, title=title
    djs_oplot, 10^!x.crange, 10^!y.crange, line=0, thick=2.0

    djs_plot, x, resid, ps=8, xsty=3, ysty=3, xrange=xrange, yrange=yrange2, $
      xthick=postthick, ythick=postthick, charsize=charsize, charthick=postthick, $
      position=pos[*,1], /noerase, xtitle=xtitle, ytitle='Residuals [%]', /xlog
    djs_oplot, 10^!x.crange, [0,0], line=0, thick=2.0

    im_openclose, postscript=postscript, /close

; 100 micron

    psname = 'soifer_vs_moshir_100'
    im_openclose, psname, xsize=8.5, ysize=9.0, postscript=postscript

    xtitle = 'Soifer et al. (1989)'
    ytitle = 'Moshir et al. (1990)'
    title = '100 micron'
    
    good = where((photo.soifer_iras_100 gt 0.0) and (photo.moshir_iras_100 gt 0.0))

    x = photo[good].soifer_iras_100
    y = photo[good].moshir_iras_100
    resid = 100.0*(y/x-1)
    
    print, 'Soifer vs Moshir - 100'
    niceprint, photo[good].galaxy, x, y, resid

    xrange = [min(x)<min(y),max(x)>max(y)]
    yrange = xrange
    yrange2 = max(abs(resid))*[-1.1,1.1]
    
    djs_plot, x, y, ps=8, xsty=3, ysty=3, xrange=xrange, yrange=yrange, $
      xthick=postthick, ythick=postthick, charsize=charsize, charthick=postthick, $
      position=pos[*,0], xtickname=replicate(' ',10), xtitle='', $
      ytitle=ytitle, /xlog, /ylog, title=title
    djs_oplot, 10^!x.crange, 10^!y.crange, line=0, thick=2.0

    djs_plot, x, resid, ps=8, xsty=3, ysty=3, xrange=xrange, yrange=yrange2, $
      xthick=postthick, ythick=postthick, charsize=charsize, charthick=postthick, $
      position=pos[*,1], /noerase, xtitle=xtitle, ytitle='Residuals [%]', /xlog
    djs_oplot, 10^!x.crange, [0,0], line=0, thick=2.0

    im_openclose, postscript=postscript, /close

stop    
    
; ---------------------------------------------------------------------------    
; RICE vs MOSHIR    
; ---------------------------------------------------------------------------    

; 12 micron

    psname = 'rice_vs_moshir_12'
    im_openclose, psname, xsize=8.5, ysize=9.0, postscript=postscript

    xtitle = 'Rice et al. (1988)'
    ytitle = 'Moshir et al. (1990)'
    title = '12 micron'
    
    good = where((photo.rice_iras_12 gt 0.0) and (photo.moshir_iras_12 gt 0.0))

    x = photo[good].rice_iras_12
    y = photo[good].moshir_iras_12
    resid = 100.0*(y/x-1)

    print, 'Rice vs Moshir - 12'
    niceprint, photo[good].galaxy, x, y, resid

    xrange = [min(x)<min(y),max(x)>max(y)]
    yrange = xrange
    yrange2 = max(abs(resid))*[-1.1,1.1]
    
    djs_plot, x, y, ps=8, xsty=3, ysty=3, xrange=xrange, yrange=yrange, $
      xthick=postthick, ythick=postthick, charsize=charsize, charthick=postthick, $
      position=pos[*,0], xtickname=replicate(' ',10), xtitle='', $
      ytitle=ytitle, /xlog, /ylog, title=title
    djs_oplot, 10^!x.crange, 10^!y.crange, line=0, thick=2.0

    djs_plot, x, resid, ps=8, xsty=3, ysty=3, xrange=xrange, yrange=yrange2, $
      xthick=postthick, ythick=postthick, charsize=charsize, charthick=postthick, $
      position=pos[*,1], /noerase, xtitle=xtitle, ytitle='Residuals [%]', /xlog
    djs_oplot, 10^!x.crange, [0,0], line=0, thick=2.0

    im_openclose, postscript=postscript, /close

; 25 micron

    psname = 'rice_vs_moshir_25'
    im_openclose, psname, xsize=8.5, ysize=9.0, postscript=postscript

    xtitle = 'Rice et al. (1988)'
    ytitle = 'Moshir et al. (1990)'
    title = '25 micron'
    
    good = where((photo.rice_iras_25 gt 0.0) and (photo.moshir_iras_25 gt 0.0))

    x = photo[good].rice_iras_25
    y = photo[good].moshir_iras_25
    resid = 100.0*(y/x-1)
    
    print, 'Rice vs Moshir - 25'
    niceprint, photo[good].galaxy, x, y, resid

    xrange = [min(x)<min(y),max(x)>max(y)]
    yrange = xrange
    yrange2 = max(abs(resid))*[-1.1,1.1]
    
    djs_plot, x, y, ps=8, xsty=3, ysty=3, xrange=xrange, yrange=yrange, $
      xthick=postthick, ythick=postthick, charsize=charsize, charthick=postthick, $
      position=pos[*,0], xtickname=replicate(' ',10), xtitle='', $
      ytitle=ytitle, /xlog, /ylog, title=title
    djs_oplot, 10^!x.crange, 10^!y.crange, line=0, thick=2.0

    djs_plot, x, resid, ps=8, xsty=3, ysty=3, xrange=xrange, yrange=yrange2, $
      xthick=postthick, ythick=postthick, charsize=charsize, charthick=postthick, $
      position=pos[*,1], /noerase, xtitle=xtitle, ytitle='Residuals [%]', /xlog
    djs_oplot, 10^!x.crange, [0,0], line=0, thick=2.0

    im_openclose, postscript=postscript, /close

; 60 micron

    psname = 'rice_vs_moshir_60'
    im_openclose, psname, xsize=8.5, ysize=9.0, postscript=postscript

    xtitle = 'Rice et al. (1988)'
    ytitle = 'Moshir et al. (1990)'
    title = '60 micron'
    
    good = where((photo.rice_iras_60 gt 0.0) and (photo.moshir_iras_60 gt 0.0))

    x = photo[good].rice_iras_60
    y = photo[good].moshir_iras_60
    resid = 100.0*(y/x-1)
    
    print, 'Rice vs Moshir - 60'
    niceprint, photo[good].galaxy, x, y, resid

    xrange = [min(x)<min(y),max(x)>max(y)]
    yrange = xrange
    yrange2 = max(abs(resid))*[-1.1,1.1]
    
    djs_plot, x, y, ps=8, xsty=3, ysty=3, xrange=xrange, yrange=yrange, $
      xthick=postthick, ythick=postthick, charsize=charsize, charthick=postthick, $
      position=pos[*,0], xtickname=replicate(' ',10), xtitle='', $
      ytitle=ytitle, /xlog, /ylog, title=title
    djs_oplot, 10^!x.crange, 10^!y.crange, line=0, thick=2.0

    djs_plot, x, resid, ps=8, xsty=3, ysty=3, xrange=xrange, yrange=yrange2, $
      xthick=postthick, ythick=postthick, charsize=charsize, charthick=postthick, $
      position=pos[*,1], /noerase, xtitle=xtitle, ytitle='Residuals [%]', /xlog
    djs_oplot, 10^!x.crange, [0,0], line=0, thick=2.0

    im_openclose, postscript=postscript, /close

; 100 micron

    psname = 'rice_vs_moshir_100'
    im_openclose, psname, xsize=8.5, ysize=9.0, postscript=postscript

    xtitle = 'Rice et al. (1988)'
    ytitle = 'Moshir et al. (1990)'
    title = '100 micron'
    
    good = where((photo.rice_iras_100 gt 0.0) and (photo.moshir_iras_100 gt 0.0))

    x = photo[good].rice_iras_100
    y = photo[good].moshir_iras_100
    resid = 100.0*(y/x-1)
    
    print, 'Rice vs Moshir - 100'
    niceprint, photo[good].galaxy, x, y, resid

    xrange = [min(x)<min(y),max(x)>max(y)]
    yrange = xrange
    yrange2 = max(abs(resid))*[-1.1,1.1]
    
    djs_plot, x, y, ps=8, xsty=3, ysty=3, xrange=xrange, yrange=yrange, $
      xthick=postthick, ythick=postthick, charsize=charsize, charthick=postthick, $
      position=pos[*,0], xtickname=replicate(' ',10), xtitle='', $
      ytitle=ytitle, /xlog, /ylog, title=title
    djs_oplot, 10^!x.crange, 10^!y.crange, line=0, thick=2.0

    djs_plot, x, resid, ps=8, xsty=3, ysty=3, xrange=xrange, yrange=yrange2, $
      xthick=postthick, ythick=postthick, charsize=charsize, charthick=postthick, $
      position=pos[*,1], /noerase, xtitle=xtitle, ytitle='Residuals [%]', /xlog
    djs_oplot, 10^!x.crange, [0,0], line=0, thick=2.0

    im_openclose, postscript=postscript, /close

; ---------------------------------------------------------------------------    
; RICE vs SOIFER    
; ---------------------------------------------------------------------------    

; 12 micron

    psname = 'rice_vs_soifer_12'
    im_openclose, psname, xsize=8.5, ysize=9.0, postscript=postscript

    xtitle = 'Rice et al. (1989)'
    ytitle = 'Soifer et al. (1990)'
    title = '12 micron'
    
    good = where((photo.rice_iras_12 gt 0.0) and (photo.soifer_iras_12 gt 0.0))

    x = photo[good].rice_iras_12
    y = photo[good].soifer_iras_12
    resid = 100.0*(y/x-1)
    
    print, 'Rice vs Soifer - 12'
    niceprint, photo[good].galaxy, x, y, resid

    xrange = [min(x)<min(y),max(x)>max(y)]
    yrange = xrange
    yrange2 = max(abs(resid))*[-1.1,1.1]
    
    djs_plot, x, y, ps=8, xsty=3, ysty=3, xrange=xrange, yrange=yrange, $
      xthick=postthick, ythick=postthick, charsize=charsize, charthick=postthick, $
      position=pos[*,0], xtickname=replicate(' ',10), xtitle='', $
      ytitle=ytitle, /xlog, /ylog, title=title
    djs_oplot, 10^!x.crange, 10^!y.crange, line=0, thick=2.0

    djs_plot, x, resid, ps=8, xsty=3, ysty=3, xrange=xrange, yrange=yrange2, $
      xthick=postthick, ythick=postthick, charsize=charsize, charthick=postthick, $
      position=pos[*,1], /noerase, xtitle=xtitle, ytitle='Residuals [%]', /xlog
    djs_oplot, 10^!x.crange, [0,0], line=0, thick=2.0

    im_openclose, postscript=postscript, /close

; 25 micron

    psname = 'rice_vs_soifer_25'
    im_openclose, psname, xsize=8.5, ysize=9.0, postscript=postscript

    xtitle = 'Rice et al. (1989)'
    ytitle = 'Soifer et al. (1990)'
    title = '25 micron'
    
    good = where((photo.rice_iras_25 gt 0.0) and (photo.soifer_iras_25 gt 0.0))

    x = photo[good].rice_iras_25
    y = photo[good].soifer_iras_25
    resid = 100.0*(y/x-1)
    
    print, 'Rice vs Soifer - 25'
    niceprint, photo[good].galaxy, x, y, resid

    xrange = [min(x)<min(y),max(x)>max(y)]
    yrange = xrange
    yrange2 = max(abs(resid))*[-1.1,1.1]
    
    djs_plot, x, y, ps=8, xsty=3, ysty=3, xrange=xrange, yrange=yrange, $
      xthick=postthick, ythick=postthick, charsize=charsize, charthick=postthick, $
      position=pos[*,0], xtickname=replicate(' ',10), xtitle='', $
      ytitle=ytitle, /xlog, /ylog, title=title
    djs_oplot, 10^!x.crange, 10^!y.crange, line=0, thick=2.0

    djs_plot, x, resid, ps=8, xsty=3, ysty=3, xrange=xrange, yrange=yrange2, $
      xthick=postthick, ythick=postthick, charsize=charsize, charthick=postthick, $
      position=pos[*,1], /noerase, xtitle=xtitle, ytitle='Residuals [%]', /xlog
    djs_oplot, 10^!x.crange, [0,0], line=0, thick=2.0

    im_openclose, postscript=postscript, /close

; 60 micron

    psname = 'rice_vs_soifer_60'
    im_openclose, psname, xsize=8.5, ysize=9.0, postscript=postscript

    xtitle = 'Rice et al. (1989)'
    ytitle = 'Soifer et al. (1990)'
    title = '60 micron'
    
    good = where((photo.rice_iras_60 gt 0.0) and (photo.soifer_iras_60 gt 0.0))

    x = photo[good].rice_iras_60
    y = photo[good].soifer_iras_60
    resid = 100.0*(y/x-1)
    
    print, 'Rice vs Soifer - 60'
    niceprint, photo[good].galaxy, x, y, resid

    xrange = [min(x)<min(y),max(x)>max(y)]
    yrange = xrange
    yrange2 = max(abs(resid))*[-1.1,1.1]
    
    djs_plot, x, y, ps=8, xsty=3, ysty=3, xrange=xrange, yrange=yrange, $
      xthick=postthick, ythick=postthick, charsize=charsize, charthick=postthick, $
      position=pos[*,0], xtickname=replicate(' ',10), xtitle='', $
      ytitle=ytitle, /xlog, /ylog, title=title
    djs_oplot, 10^!x.crange, 10^!y.crange, line=0, thick=2.0

    djs_plot, x, resid, ps=8, xsty=3, ysty=3, xrange=xrange, yrange=yrange2, $
      xthick=postthick, ythick=postthick, charsize=charsize, charthick=postthick, $
      position=pos[*,1], /noerase, xtitle=xtitle, ytitle='Residuals [%]', /xlog
    djs_oplot, 10^!x.crange, [0,0], line=0, thick=2.0

    im_openclose, postscript=postscript, /close

; 100 micron

    psname = 'rice_vs_soifer_100'
    im_openclose, psname, xsize=8.5, ysize=9.0, postscript=postscript

    xtitle = 'Rice et al. (1989)'
    ytitle = 'Soifer et al. (1990)'
    title = '100 micron'
    
    good = where((photo.rice_iras_100 gt 0.0) and (photo.soifer_iras_100 gt 0.0))

    x = photo[good].rice_iras_100
    y = photo[good].soifer_iras_100
    resid = 100.0*(y/x-1)
    
    print, 'Rice vs Soifer - 100'
    niceprint, photo[good].galaxy, x, y, resid

    xrange = [min(x)<min(y),max(x)>max(y)]
    yrange = xrange
    yrange2 = max(abs(resid))*[-1.1,1.1]
    
    djs_plot, x, y, ps=8, xsty=3, ysty=3, xrange=xrange, yrange=yrange, $
      xthick=postthick, ythick=postthick, charsize=charsize, charthick=postthick, $
      position=pos[*,0], xtickname=replicate(' ',10), xtitle='', $
      ytitle=ytitle, /xlog, /ylog, title=title
    djs_oplot, 10^!x.crange, 10^!y.crange, line=0, thick=2.0

    djs_plot, x, resid, ps=8, xsty=3, ysty=3, xrange=xrange, yrange=yrange2, $
      xthick=postthick, ythick=postthick, charsize=charsize, charthick=postthick, $
      position=pos[*,1], /noerase, xtitle=xtitle, ytitle='Residuals [%]', /xlog
    djs_oplot, 10^!x.crange, [0,0], line=0, thick=2.0

    im_openclose, postscript=postscript, /close

stop    


return
end
    
