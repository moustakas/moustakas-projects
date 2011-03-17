pro nfgs_compare, atlas, nfgs, jansen, postscript=postscript
; jm04mar21uofa
; jm04dec10uofa

    if keyword_set(postscript) then begin
       dfpsplot, 'compare_ispec1d_nfgs.ps', /square, xsize=8.5, ysize=6.5
       postthick = 8.0 
    endif else begin
       im_window, 0, xratio=0.5, /square
       postthick = 2.0
    endelse

    pagemaker, nx=2, ny=2, width=[2.9,2.9], height=[2.8,2.2], $
      position=pos, /normal, xspace=1.2, yspace=0.0, $
      xmargin=[1.1,0.4], ymargin=[0.4,1.1], xpage=8.5, ypage=6.5
    
;   pspath = atlas_path(/papers)+'atlas/'
    
    if (n_elements(atlas) eq 0L) then atlas = read_integrated()
    if (n_elements(nfgs) eq 0L) then nfgs = read_nfgs()
    if (n_elements(jansen) eq 0L) then jansen = read_00jansen()
    
; compare my measurements of the NFGS spectra with Rolf's measurements 

    match, strtrim(nfgs.ned_galaxy,2), strtrim(jansen.nedgalaxy,2), nfgsmatch, jansenmatch
;   match, nfgs.nfgs_id, jansen.nfgs_id, nfgsmatch, jansenmatch
    niceprint, nfgs[nfgsmatch].ned_galaxy, jansen[jansenmatch].nedgalaxy
    
    plotsym, 0, 0.8, /fill
    scale = 1E16
    ytitle2 = 'Residuals ['+angstrom()+']'
    
; --------------------------------------------------    
; [O II]
; --------------------------------------------------    

; equivalent widths    
    
    indx = where((nfgs[nfgsmatch].oii_3727_ew[1] gt 0.0) and $
      (jansen[jansenmatch].oii_3727_ew[1] gt 0.0),nindx)

    x = nfgs[nfgsmatch[indx]].oii_3727_ew[0]
    y = jansen[jansenmatch[indx]].oii_3727_ew[0]
    resid = y-x
;   resid = 100*(y-x)/x
    stats = im_stats(resid,sigrej=3.0,/verbose)

    xstr = strtrim(string(stats.median_rej,format='(F12.2)'),2)+' ('+$
      strtrim(string(stats.mean_rej,format='(F12.2)'),2)+$
      '\pm'+strtrim(string(stats.sigma_rej,format='(F12.2)'),2)+')'

    xtitle = 'EW([O II] \lambda3727) [iSPEC1d]'
    ytitle = 'EW([O II] \lambda3727) [Jansen]'

    xrange = [min(x)<min(y),max(x)>max(y)]
    yrange = xrange
    yrange2 = (abs(min(resid))>abs(max(resid)))*[-1.1,1.1]
    
    djs_plot, x, y, xtitle='', ytitle=ytitle, ps=8, position=pos[*,0], $
      xrange=xrange, yrange=yrange, xstyle=3, ystyle=3, yminor=3, $
      xthick=postthick, ythick=postthick, charsize=1.2, charthick=postthick, $
      xtickname=replicate(' ',10), /xlog, /ylog
    djs_oplot, 10^!x.crange, 10^!y.crange, line=0, thick=postthick

    djs_plot, x, resid, /noerase, xtitle=xtitle, ytitle=ytitle2, ps=8, $
      position=pos[*,2], xrange=xrange, yrange=yrange2, xstyle=3, ystyle=3, yminor=3, $
      xthick=postthick, ythick=postthick, charsize=1.2, charthick=postthick, /xlog
    djs_oplot, 10^!x.crange, [0,0], line=0, thick=postthick

    legend, textoidl(xstr), /right, /top, box=0, charsize=1.0, $
      charthick=postthick, clear=keyword_set(postscript)

; fluxes
    
    indx = where((nfgs[nfgsmatch].oii_3727[1] gt 0.0) and $
      (jansen[jansenmatch].oii_3727[1] gt 0.0),nindx)

    x = nfgs[nfgsmatch[indx]].oii_3727[0]
    y = jansen[jansenmatch[indx]].oii_3727[0]

    scale = median(y/x)
    x = scale*x
    
    resid = 100*(y-x)/x
    stats = im_stats(resid,sigrej=3.0,/verbose,/no_head)
    
    xstr = strtrim(string(stats.median_rej,format='(F12.2)'),2)+' ('+$
      strtrim(string(stats.mean_rej,format='(F12.2)'),2)+$
      '\pm'+strtrim(string(stats.sigma_rej,format='(F12.2)'),2)+')'

    xtitle = 'F([O II] \lambda3727) [iSPEC1d]'
    ytitle = 'F([O II] \lambda3727) [Jansen]'

    xrange = [min(x)<min(y),max(x)>max(y)]
    yrange = xrange
    yrange2 = (abs(min(resid))>abs(max(resid)))*[-1.1,1.1]
    
    djs_plot, x, y, xtitle='', ytitle=ytitle, ps=8, position=pos[*,1], $
      xrange=xrange, yrange=yrange, xstyle=3, ystyle=3, yminor=3, $
      xthick=postthick, ythick=postthick, charsize=1.2, charthick=postthick, $
      xtickname=replicate(' ',10), /xlog, /ylog, /noerase
    djs_oplot, 10^!x.crange, 10^!y.crange, line=0, thick=postthick

    djs_plot, x, resid, /noerase, xtitle=xtitle, ytitle='Residuals [%]', ps=8, $
      position=pos[*,3], xrange=xrange, yrange=yrange2, xstyle=3, ystyle=3, yminor=3, $
      xthick=postthick, ythick=postthick, charsize=1.2, charthick=postthick, /xlog
    djs_oplot, 10^!x.crange, [0,0], line=0, thick=postthick

    legend, textoidl(xstr), /right, /top, box=0, charsize=1.0, $
      charthick=postthick, clear=keyword_set(postscript)

    if not keyword_set(postscript) then cc = get_kbrd(1)

; --------------------------------------------------    
; H-beta    
; --------------------------------------------------    

    indx = where((nfgs[nfgsmatch].h_beta_ew[1] gt 0.0) and $
      (jansen[jansenmatch].h_beta_ew[1] gt 0.0),nindx)

    x = nfgs[nfgsmatch[indx]].h_beta_ew[0]
    y = jansen[jansenmatch[indx]].h_beta_ew[0]
    resid = y-x
;   resid = 100*(y-x)/x
    stats = im_stats(resid,sigrej=3.0,/verbose,/no_head)

    xstr = strtrim(string(stats.median_rej,format='(F12.2)'),2)+' ('+$
      strtrim(string(stats.mean_rej,format='(F12.2)'),2)+$
      '\pm'+strtrim(string(stats.sigma_rej,format='(F12.2)'),2)+')'

    xtitle = 'EW(H\beta) [iSPEC1d]'
    ytitle = 'EW(H\beta) [Jansen]'

    xrange = [min(x)<min(y),max(x)>max(y)]
    yrange = xrange
    yrange2 = (abs(min(resid))>abs(max(resid)))*[-1.1,1.1]
    
    djs_plot, x, y, xtitle='', ytitle=ytitle, ps=8, position=pos[*,0], $
      xrange=xrange, yrange=yrange, xstyle=3, ystyle=3, yminor=3, $
      xthick=postthick, ythick=postthick, charsize=1.2, charthick=postthick, $
      xtickname=replicate(' ',10), /xlog, /ylog
    djs_oplot, 10^!x.crange, 10^!y.crange, line=0, thick=postthick

    djs_plot, x, resid, /noerase, xtitle=xtitle, ytitle=ytitle2, ps=8, $
      position=pos[*,2], xrange=xrange, yrange=yrange2, xstyle=3, ystyle=3, yminor=3, $
      xthick=postthick, ythick=postthick, charsize=1.2, charthick=postthick, /xlog
    djs_oplot, 10^!x.crange, [0,0], line=0, thick=postthick
    
    legend, textoidl(xstr), /right, /top, box=0, charsize=1.0, $
      charthick=postthick, clear=keyword_set(postscript)

; fluxes
    
    indx = where((nfgs[nfgsmatch].h_beta[1] gt 0.0) and $
      (jansen[jansenmatch].h_beta[1] gt 0.0),nindx)

    x = nfgs[nfgsmatch[indx]].h_beta[0]
    y = jansen[jansenmatch[indx]].h_beta[0]

    scale = median(y/x)
    x = scale*x
    
    resid = 100*(y-x)/x
    stats = im_stats(resid,sigrej=3.0,/verbose,/no_head)
    
    xstr = strtrim(string(stats.median_rej,format='(F12.2)'),2)+' ('+$
      strtrim(string(stats.mean_rej,format='(F12.2)'),2)+$
      '\pm'+strtrim(string(stats.sigma_rej,format='(F12.2)'),2)+')'

    xtitle = 'F(H\beta) [iSPEC1d]'
    ytitle = 'F(H\beta) [Jansen]'

    xrange = [min(x)<min(y),max(x)>max(y)]
    yrange = xrange
    yrange2 = (abs(min(resid))>abs(max(resid)))*[-1.1,1.1]
    
    djs_plot, x, y, xtitle='', ytitle=ytitle, ps=8, position=pos[*,1], $
      xrange=xrange, yrange=yrange, xstyle=3, ystyle=3, yminor=3, $
      xthick=postthick, ythick=postthick, charsize=1.2, charthick=postthick, $
      xtickname=replicate(' ',10), /xlog, /ylog, /noerase
    djs_oplot, 10^!x.crange, 10^!y.crange, line=0, thick=postthick

    djs_plot, x, resid, /noerase, xtitle=xtitle, ytitle='Residuals [%]', ps=8, $
      position=pos[*,3], xrange=xrange, yrange=yrange2, xstyle=3, ystyle=3, yminor=3, $
      xthick=postthick, ythick=postthick, charsize=1.2, charthick=postthick, /xlog
    djs_oplot, 10^!x.crange, [0,0], line=0, thick=postthick

    legend, textoidl(xstr), /right, /top, box=0, charsize=1.0, $
      charthick=postthick, clear=keyword_set(postscript)

    if not keyword_set(postscript) then cc = get_kbrd(1)

; --------------------------------------------------    
; [O III] 4959
; --------------------------------------------------    

    indx = where((nfgs[nfgsmatch].oiii_4959_ew_orig[1] gt 0.0) and $
      (jansen[jansenmatch].oiii_4959_ew[1] gt 0.0),nindx)

    x = nfgs[nfgsmatch[indx]].oiii_4959_ew_orig[0]
    y = jansen[jansenmatch[indx]].oiii_4959_ew[0]
    resid = y-x
;   resid = 100*(y-x)/x
    stats = im_stats(resid,sigrej=3.0,/verbose,/no_head)

    xstr = strtrim(string(stats.median_rej,format='(F12.2)'),2)+' ('+$
      strtrim(string(stats.mean_rej,format='(F12.2)'),2)+$
      '\pm'+strtrim(string(stats.sigma_rej,format='(F12.2)'),2)+')'

    xtitle = 'EW([O III] \lambda4959) [iSPEC1d]'
    ytitle = 'EW([O III] \lambda4959) [Jansen]'

    xrange = [min(x)<min(y),max(x)>max(y)]
    yrange = xrange
    yrange2 = (abs(min(resid))>abs(max(resid)))*[-1.1,1.1]
    
    djs_plot, x, y, xtitle='', ytitle=ytitle, ps=8, position=pos[*,0], $
      xrange=xrange, yrange=yrange, xstyle=3, ystyle=3, yminor=3, $
      xthick=postthick, ythick=postthick, charsize=1.2, charthick=postthick, $
      xtickname=replicate(' ',10), /xlog, /ylog
    djs_oplot, 10^!x.crange, 10^!y.crange, line=0, thick=postthick

    djs_plot, x, resid, /noerase, xtitle=xtitle, ytitle=ytitle2, ps=8, $
      position=pos[*,2], xrange=xrange, yrange=yrange2, xstyle=3, ystyle=3, yminor=3, $
      xthick=postthick, ythick=postthick, charsize=1.2, charthick=postthick, /xlog
    djs_oplot, 10^!x.crange, [0,0], line=0, thick=postthick

    legend, textoidl(xstr), /right, /top, box=0, charsize=1.0, $
      charthick=postthick, clear=keyword_set(postscript)

; fluxes
    
    indx = where((nfgs[nfgsmatch].oiii_4959[1] gt 0.0) and $
      (jansen[jansenmatch].oiii_4959[1] gt 0.0),nindx)

    x = nfgs[nfgsmatch[indx]].oiii_4959[0]
    y = jansen[jansenmatch[indx]].oiii_4959[0]

    scale = median(y/x)
    x = scale*x
    
    resid = 100*(y-x)/x
    stats = im_stats(resid,sigrej=3.0,/verbose,/no_head)
    
    xstr = strtrim(string(stats.median_rej,format='(F12.2)'),2)+' ('+$
      strtrim(string(stats.mean_rej,format='(F12.2)'),2)+$
      '\pm'+strtrim(string(stats.sigma_rej,format='(F12.2)'),2)+')'

    xtitle = 'F([O III] \lambda4959) [iSPEC1d]'
    ytitle = 'F([O III] \lambda4959) [Jansen]'

    xrange = [min(x)<min(y),max(x)>max(y)]
    yrange = xrange
    yrange2 = (abs(min(resid))>abs(max(resid)))*[-1.1,1.1]
    
    djs_plot, x, y, xtitle='', ytitle=ytitle, ps=8, position=pos[*,1], $
      xrange=xrange, yrange=yrange, xstyle=3, ystyle=3, yminor=3, $
      xthick=postthick, ythick=postthick, charsize=1.2, charthick=postthick, $
      xtickname=replicate(' ',10), /xlog, /ylog, /noerase
    djs_oplot, 10^!x.crange, 10^!y.crange, line=0, thick=postthick

    djs_plot, x, resid, /noerase, xtitle=xtitle, ytitle='Residuals [%]', ps=8, $
      position=pos[*,3], xrange=xrange, yrange=yrange2, xstyle=3, ystyle=3, yminor=3, $
      xthick=postthick, ythick=postthick, charsize=1.2, charthick=postthick, /xlog
    djs_oplot, 10^!x.crange, [0,0], line=0, thick=postthick
    
    legend, textoidl(xstr), /right, /top, box=0, charsize=1.0, $
      charthick=postthick, clear=keyword_set(postscript)

    if not keyword_set(postscript) then cc = get_kbrd(1)

; --------------------------------------------------    
; [O III] 5007
; --------------------------------------------------    

    indx = where((nfgs[nfgsmatch].oiii_5007_ew[1] gt 0.0) and $
      (jansen[jansenmatch].oiii_5007_ew[1] gt 0.0),nindx)

    x = nfgs[nfgsmatch[indx]].oiii_5007_ew[0]
    y = jansen[jansenmatch[indx]].oiii_5007_ew[0]
    resid = y-x
;   resid = 100*(y-x)/x
    stats = im_stats(resid,sigrej=3.0,/verbose,/no_head)

    xstr = strtrim(string(stats.median_rej,format='(F12.2)'),2)+' ('+$
      strtrim(string(stats.mean_rej,format='(F12.2)'),2)+$
      '\pm'+strtrim(string(stats.sigma_rej,format='(F12.2)'),2)+')'

    xtitle = 'EW([O III] \lambda5007) [iSPEC1d]'
    ytitle = 'EW([O III] \lambda5007) [Jansen]'

    xrange = [min(x)<min(y),max(x)>max(y)]
    yrange = xrange
    yrange2 = (abs(min(resid))>abs(max(resid)))*[-1.1,1.1]
    
    djs_plot, x, y, xtitle='', ytitle=ytitle, ps=8, position=pos[*,0], $
      xrange=xrange, yrange=yrange, xstyle=3, ystyle=3, yminor=3, $
      xthick=postthick, ythick=postthick, charsize=1.2, charthick=postthick, $
      xtickname=replicate(' ',10), /xlog, /ylog
    djs_oplot, 10^!x.crange, 10^!y.crange, line=0, thick=postthick

    djs_plot, x, resid, /noerase, xtitle=xtitle, ytitle='Residuals ['+angstrom()+']', ps=8, $
      position=pos[*,2], xrange=xrange, yrange=yrange2, xstyle=3, ystyle=3, yminor=3, $
      xthick=postthick, ythick=postthick, charsize=1.2, charthick=postthick, /xlog
    djs_oplot, 10^!x.crange, [0,0], line=0, thick=postthick
    
    legend, textoidl(xstr), /right, /top, box=0, charsize=1.0, $
      charthick=postthick, clear=keyword_set(postscript)

; fluxes
    
    indx = where((nfgs[nfgsmatch].oiii_5007[1] gt 0.0) and $
      (jansen[jansenmatch].oiii_5007[1] gt 0.0),nindx)

    x = nfgs[nfgsmatch[indx]].oiii_5007[0]
    y = jansen[jansenmatch[indx]].oiii_5007[0]

    scale = median(y/x)
    x = scale*x
    
    resid = 100*(y-x)/x
    stats = im_stats(resid,sigrej=3.0,/verbose,/no_head)
    
    xstr = strtrim(string(stats.median_rej,format='(F12.2)'),2)+' ('+$
      strtrim(string(stats.mean_rej,format='(F12.2)'),2)+$
      '\pm'+strtrim(string(stats.sigma_rej,format='(F12.2)'),2)+')'

    xtitle = 'F([O III] \lambda5007) [iSPEC1d]'
    ytitle = 'F([O III] \lambda5007) [Jansen]'

    xrange = [min(x)<min(y),max(x)>max(y)]
    yrange = xrange
    yrange2 = (abs(min(resid))>abs(max(resid)))*[-1.1,1.1]
    
    djs_plot, x, y, xtitle='', ytitle=ytitle, ps=8, position=pos[*,1], $
      xrange=xrange, yrange=yrange, xstyle=3, ystyle=3, yminor=3, $
      xthick=postthick, ythick=postthick, charsize=1.2, charthick=postthick, $
      xtickname=replicate(' ',10), /xlog, /ylog, /noerase
    djs_oplot, 10^!x.crange, 10^!y.crange, line=0, thick=postthick

    djs_plot, x, resid, /noerase, xtitle=xtitle, ytitle=ytitle2, ps=8, $
      position=pos[*,3], xrange=xrange, yrange=yrange2, xstyle=3, ystyle=3, yminor=3, $
      xthick=postthick, ythick=postthick, charsize=1.2, charthick=postthick, /xlog
    djs_oplot, 10^!x.crange, [0,0], line=0, thick=postthick
    
    legend, textoidl(xstr), /right, /top, box=0, charsize=1.0, $
      charthick=postthick, clear=keyword_set(postscript)

    if not keyword_set(postscript) then cc = get_kbrd(1)

; --------------------------------------------------    
; [N II] 6548
; --------------------------------------------------    

    indx = where((nfgs[nfgsmatch].nii_6548_ew[1] gt 0.0) and $
      (jansen[jansenmatch].nii_6548_ew[1] gt 0.0),nindx)

    x = nfgs[nfgsmatch[indx]].nii_6548_ew[0]
    y = jansen[jansenmatch[indx]].nii_6548_ew[0]
    resid = y-x
;   resid = 100*(y-x)/x
    stats = im_stats(resid,sigrej=3.0,/verbose,/no_head)

    xstr = strtrim(string(stats.median_rej,format='(F12.2)'),2)+' ('+$
      strtrim(string(stats.mean_rej,format='(F12.2)'),2)+$
      '\pm'+strtrim(string(stats.sigma_rej,format='(F12.2)'),2)+')'

    xtitle = 'EW([N II] \lambda6548) [iSPEC1d]'
    ytitle = 'EW([N II] \lambda6548) [Jansen]'

    xrange = [min(x)<min(y),max(x)>max(y)]
    yrange = xrange
    yrange2 = (abs(min(resid))>abs(max(resid)))*[-1.1,1.1]
    
    djs_plot, x, y, xtitle='', ytitle=ytitle, ps=8, position=pos[*,0], $
      xrange=xrange, yrange=yrange, xstyle=3, ystyle=3, yminor=3, $
      xthick=postthick, ythick=postthick, charsize=1.2, charthick=postthick, $
      xtickname=replicate(' ',10), /xlog, /ylog
    djs_oplot, 10^!x.crange, 10^!y.crange, line=0, thick=postthick

    djs_plot, x, resid, /noerase, xtitle=xtitle, ytitle='Residuals ['+angstrom()+']', ps=8, $
      position=pos[*,2], xrange=xrange, yrange=yrange2, xstyle=3, ystyle=3, yminor=3, $
      xthick=postthick, ythick=postthick, charsize=1.2, charthick=postthick, /xlog
    djs_oplot, 10^!x.crange, [0,0], line=0, thick=postthick
    
    legend, textoidl(xstr), /right, /top, box=0, charsize=1.0, $
      charthick=postthick, clear=keyword_set(postscript)

; fluxes
    
    indx = where((nfgs[nfgsmatch].nii_6548[1] gt 0.0) and $
      (jansen[jansenmatch].nii_6548[1] gt 0.0),nindx)

    x = nfgs[nfgsmatch[indx]].nii_6548[0]
    y = jansen[jansenmatch[indx]].nii_6548[0]

    scale = median(y/x)
    x = scale*x
    
    resid = 100*(y-x)/x
    stats = im_stats(resid,sigrej=3.0,/verbose,/no_head)
    
    xstr = strtrim(string(stats.median_rej,format='(F12.2)'),2)+' ('+$
      strtrim(string(stats.mean_rej,format='(F12.2)'),2)+$
      '\pm'+strtrim(string(stats.sigma_rej,format='(F12.2)'),2)+')'

    xtitle = 'F([N II] \lambda6548) [iSPEC1d]'
    ytitle = 'F([N II] \lambda6548) [Jansen]'

    xrange = [min(x)<min(y),max(x)>max(y)]
    yrange = xrange
    yrange2 = (abs(min(resid))>abs(max(resid)))*[-1.1,1.1]
    
    djs_plot, x, y, xtitle='', ytitle=ytitle, ps=8, position=pos[*,1], $
      xrange=xrange, yrange=yrange, xstyle=3, ystyle=3, yminor=3, $
      xthick=postthick, ythick=postthick, charsize=1.2, charthick=postthick, $
      xtickname=replicate(' ',10), /xlog, /ylog, /noerase
    djs_oplot, 10^!x.crange, 10^!y.crange, line=0, thick=postthick

    djs_plot, x, resid, /noerase, xtitle=xtitle, ytitle='Residuals [%]', ps=8, $
      position=pos[*,3], xrange=xrange, yrange=yrange2, xstyle=3, ystyle=3, yminor=3, $
      xthick=postthick, ythick=postthick, charsize=1.2, charthick=postthick, /xlog
    djs_oplot, 10^!x.crange, [0,0], line=0, thick=postthick
    
    legend, textoidl(xstr), /right, /top, box=0, charsize=1.0, $
      charthick=postthick, clear=keyword_set(postscript)

    if not keyword_set(postscript) then cc = get_kbrd(1)

; --------------------------------------------------    
; H-alpha    
; --------------------------------------------------    

    indx = where((nfgs[nfgsmatch].h_alpha_ew[1] gt 0.0) and $
      (jansen[jansenmatch].h_alpha_ew[1] gt 0.0),nindx)

    x = nfgs[nfgsmatch[indx]].h_alpha_ew[0]
    y = jansen[jansenmatch[indx]].h_alpha_ew[0]
    resid = y-x
;   resid = 100*(y-x)/x
    stats = im_stats(resid,sigrej=3.0,/verbose,/no_head)

    xstr = strtrim(string(stats.median_rej,format='(F12.2)'),2)+' ('+$
      strtrim(string(stats.mean_rej,format='(F12.2)'),2)+$
      '\pm'+strtrim(string(stats.sigma_rej,format='(F12.2)'),2)+')'

    xtitle = 'EW(H\alpha) [iSPEC1d]'
    ytitle = 'EW(H\alpha) [Jansen]'

    xrange = [min(x)<min(y),max(x)>max(y)]
    yrange = xrange
    yrange2 = (abs(min(resid))>abs(max(resid)))*[-1.1,1.1]
    
    djs_plot, x, y, xtitle='', ytitle=ytitle, ps=8, position=pos[*,0], $
      xrange=xrange, yrange=yrange, xstyle=3, ystyle=3, yminor=3, $
      xthick=postthick, ythick=postthick, charsize=1.2, charthick=postthick, $
      xtickname=replicate(' ',10), /xlog, /ylog
    djs_oplot, 10^!x.crange, 10^!y.crange, line=0, thick=postthick

    djs_plot, x, resid, /noerase, xtitle=xtitle, ytitle='Residuals ['+angstrom()+']', ps=8, $
      position=pos[*,2], xrange=xrange, yrange=yrange2, xstyle=3, ystyle=3, yminor=3, $
      xthick=postthick, ythick=postthick, charsize=1.2, charthick=postthick, /xlog
    djs_oplot, 10^!x.crange, [0,0], line=0, thick=postthick
    
    legend, textoidl(xstr), /right, /top, box=0, charsize=1.0, $
      charthick=postthick, clear=keyword_set(postscript)

; fluxes
    
    indx = where((nfgs[nfgsmatch].h_alpha[1] gt 0.0) and $
      (jansen[jansenmatch].h_alpha[1] gt 0.0),nindx)

    x = nfgs[nfgsmatch[indx]].h_alpha[0]
    y = jansen[jansenmatch[indx]].h_alpha[0]

    scale = median(y/x)
    x = scale*x
    
    resid = 100*(y-x)/x
    stats = im_stats(resid,sigrej=3.0,/verbose,/no_head)
    
    xstr = strtrim(string(stats.median_rej,format='(F12.2)'),2)+' ('+$
      strtrim(string(stats.mean_rej,format='(F12.2)'),2)+$
      '\pm'+strtrim(string(stats.sigma_rej,format='(F12.2)'),2)+')'

    xtitle = 'F(H\alpha) [iSPEC1d]'
    ytitle = 'F(H\alpha) [Jansen]'

    xrange = [min(x)<min(y),max(x)>max(y)]
    yrange = xrange
    yrange2 = (abs(min(resid))>abs(max(resid)))*[-1.1,1.1]
    
    djs_plot, x, y, xtitle='', ytitle=ytitle, ps=8, position=pos[*,1], $
      xrange=xrange, yrange=yrange, xstyle=3, ystyle=3, yminor=3, $
      xthick=postthick, ythick=postthick, charsize=1.2, charthick=postthick, $
      xtickname=replicate(' ',10), /xlog, /ylog, /noerase
    djs_oplot, 10^!x.crange, 10^!y.crange, line=0, thick=postthick

    djs_plot, x, resid, /noerase, xtitle=xtitle, ytitle=ytitle2, ps=8, $
      position=pos[*,3], xrange=xrange, yrange=yrange2, xstyle=3, ystyle=3, yminor=3, $
      xthick=postthick, ythick=postthick, charsize=1.2, charthick=postthick, /xlog
    djs_oplot, 10^!x.crange, [0,0], line=0, thick=postthick

    legend, textoidl(xstr), /right, /top, box=0, charsize=1.0, $
      charthick=postthick, clear=keyword_set(postscript)

    if not keyword_set(postscript) then cc = get_kbrd(1)

; --------------------------------------------------    
; [N II] 6584
; --------------------------------------------------    

    indx = where((nfgs[nfgsmatch].nii_6584_ew[1] gt 0.0) and $
      (jansen[jansenmatch].nii_6584_ew[1] gt 0.0),nindx)

    x = nfgs[nfgsmatch[indx]].nii_6584_ew[0]
    y = jansen[jansenmatch[indx]].nii_6584_ew[0]
    resid = y-x
;   resid = 100*(y-x)/x
    stats = im_stats(resid,sigrej=3.0,/verbose,/no_head)

    xstr = strtrim(string(stats.median_rej,format='(F12.2)'),2)+' ('+$
      strtrim(string(stats.mean_rej,format='(F12.2)'),2)+$
      '\pm'+strtrim(string(stats.sigma_rej,format='(F12.2)'),2)+')'

    xtitle = 'EW([N II] \lambda6584) [iSPEC1d]'
    ytitle = 'EW([N II] \lambda6584) [Jansen]'

    xrange = [min(x)<min(y),max(x)>max(y)]
    yrange = xrange
    yrange2 = (abs(min(resid))>abs(max(resid)))*[-1.1,1.1]
    
    djs_plot, x, y, xtitle='', ytitle=ytitle, ps=8, position=pos[*,0], $
      xrange=xrange, yrange=yrange, xstyle=3, ystyle=3, yminor=3, $
      xthick=postthick, ythick=postthick, charsize=1.2, charthick=postthick, $
      xtickname=replicate(' ',10), /xlog, /ylog
    djs_oplot, 10^!x.crange, 10^!y.crange, line=0, thick=postthick

    djs_plot, x, resid, /noerase, xtitle=xtitle, ytitle='Residuals ['+angstrom()+']', ps=8, $
      position=pos[*,2], xrange=xrange, yrange=yrange2, xstyle=3, ystyle=3, yminor=3, $
      xthick=postthick, ythick=postthick, charsize=1.2, charthick=postthick, /xlog
    djs_oplot, 10^!x.crange, [0,0], line=0, thick=postthick

    legend, textoidl(xstr), /right, /top, box=0, charsize=1.0, $
      charthick=postthick, clear=keyword_set(postscript)

; fluxes
    
    indx = where((nfgs[nfgsmatch].nii_6584[1] gt 0.0) and $
      (jansen[jansenmatch].nii_6584[1] gt 0.0),nindx)

    x = nfgs[nfgsmatch[indx]].nii_6584[0]
    y = jansen[jansenmatch[indx]].nii_6584[0]

    scale = median(y/x)
    x = scale*x
    
    resid = 100*(y-x)/x
    stats = im_stats(resid,sigrej=3.0,/verbose,/no_head)
    
    xstr = strtrim(string(stats.median_rej,format='(F12.2)'),2)+' ('+$
      strtrim(string(stats.mean_rej,format='(F12.2)'),2)+$
      '\pm'+strtrim(string(stats.sigma_rej,format='(F12.2)'),2)+')'

    xtitle = 'F([N II] \lambda6584) [iSPEC1d]'
    ytitle = 'F([N II] \lambda6584) [Jansen]'

    xrange = [min(x)<min(y),max(x)>max(y)]
    yrange = xrange
    yrange2 = (abs(min(resid))>abs(max(resid)))*[-1.1,1.1]
    
    djs_plot, x, y, xtitle='', ytitle=ytitle, ps=8, position=pos[*,1], $
      xrange=xrange, yrange=yrange, xstyle=3, ystyle=3, yminor=3, $
      xthick=postthick, ythick=postthick, charsize=1.2, charthick=postthick, $
      xtickname=replicate(' ',10), /xlog, /ylog, /noerase
    djs_oplot, 10^!x.crange, 10^!y.crange, line=0, thick=postthick

    djs_plot, x, resid, /noerase, xtitle=xtitle, ytitle='Residuals [%]', ps=8, $
      position=pos[*,3], xrange=xrange, yrange=yrange2, xstyle=3, ystyle=3, yminor=3, $
      xthick=postthick, ythick=postthick, charsize=1.2, charthick=postthick, /xlog
    djs_oplot, 10^!x.crange, [0,0], line=0, thick=postthick
    
    legend, textoidl(xstr), /right, /top, box=0, charsize=1.0, $
      charthick=postthick, clear=keyword_set(postscript)

    if not keyword_set(postscript) then cc = get_kbrd(1)

; --------------------------------------------------    
; [S II] 6716
; --------------------------------------------------    

    indx = where((nfgs[nfgsmatch].sii_6716_ew[1] gt 0.0) and $
      (jansen[jansenmatch].sii_6716_ew[1] gt 0.0),nindx)

    x = nfgs[nfgsmatch[indx]].sii_6716_ew[0]
    y = jansen[jansenmatch[indx]].sii_6716_ew[0]
    resid = y-x
;   resid = 100*(y-x)/x
    stats = im_stats(resid,sigrej=3.0,/verbose,/no_head)

    xstr = strtrim(string(stats.median_rej,format='(F12.2)'),2)+' ('+$
      strtrim(string(stats.mean_rej,format='(F12.2)'),2)+$
      '\pm'+strtrim(string(stats.sigma_rej,format='(F12.2)'),2)+')'

    xtitle = 'EW([S II] \lambda6716) [iSPEC1d]'
    ytitle = 'EW([S II] \lambda6716) [Jansen]'

    xrange = [min(x)<min(y),max(x)>max(y)]
    yrange = xrange
    yrange2 = (abs(min(resid))>abs(max(resid)))*[-1.1,1.1]
    
    djs_plot, x, y, xtitle='', ytitle=ytitle, ps=8, position=pos[*,0], $
      xrange=xrange, yrange=yrange, xstyle=3, ystyle=3, yminor=3, $
      xthick=postthick, ythick=postthick, charsize=1.2, charthick=postthick, $
      xtickname=replicate(' ',10), /xlog, /ylog
    djs_oplot, 10^!x.crange, 10^!y.crange, line=0, thick=postthick

    djs_plot, x, resid, /noerase, xtitle=xtitle, ytitle='Residuals ['+angstrom()+']', ps=8, $
      position=pos[*,2], xrange=xrange, yrange=yrange2, xstyle=3, ystyle=3, yminor=3, $
      xthick=postthick, ythick=postthick, charsize=1.2, charthick=postthick, /xlog
    djs_oplot, 10^!x.crange, [0,0], line=0, thick=postthick
    
    legend, textoidl(xstr), /right, /top, box=0, charsize=1.0, $
      charthick=postthick, clear=keyword_set(postscript)

; fluxes
    
    indx = where((nfgs[nfgsmatch].sii_6716[1] gt 0.0) and $
      (jansen[jansenmatch].sii_6716[1] gt 0.0),nindx)

    x = nfgs[nfgsmatch[indx]].sii_6716[0]
    y = jansen[jansenmatch[indx]].sii_6716[0]

    scale = median(y/x)
    x = scale*x
    
    resid = 100*(y-x)/x
    stats = im_stats(resid,sigrej=3.0,/verbose,/no_head)
    
    xstr = strtrim(string(stats.median_rej,format='(F12.2)'),2)+' ('+$
      strtrim(string(stats.mean_rej,format='(F12.2)'),2)+$
      '\pm'+strtrim(string(stats.sigma_rej,format='(F12.2)'),2)+')'

    xtitle = 'F([S II] \lambda6716) [iSPEC1d]'
    ytitle = 'F([S II] \lambda6716) [Jansen]'

    xrange = [min(x)<min(y),max(x)>max(y)]
    yrange = xrange
    yrange2 = (abs(min(resid))>abs(max(resid)))*[-1.1,1.1]
    
    djs_plot, x, y, xtitle='', ytitle=ytitle, ps=8, position=pos[*,1], $
      xrange=xrange, yrange=yrange, xstyle=3, ystyle=3, yminor=3, $
      xthick=postthick, ythick=postthick, charsize=1.2, charthick=postthick, $
      xtickname=replicate(' ',10), /xlog, /ylog, /noerase
    djs_oplot, 10^!x.crange, 10^!y.crange, line=0, thick=postthick

    djs_plot, x, resid, /noerase, xtitle=xtitle, ytitle=ytitle2, ps=8, $
      position=pos[*,3], xrange=xrange, yrange=yrange2, xstyle=3, ystyle=3, yminor=3, $
      xthick=postthick, ythick=postthick, charsize=1.2, charthick=postthick, /xlog
    djs_oplot, 10^!x.crange, [0,0], line=0, thick=postthick
    
    legend, textoidl(xstr), /right, /top, box=0, charsize=1.0, $
      charthick=postthick, clear=keyword_set(postscript)

    if not keyword_set(postscript) then cc = get_kbrd(1)

; --------------------------------------------------    
; [S II] 6731
; --------------------------------------------------    

    indx = where((nfgs[nfgsmatch].sii_6731_ew[1] gt 0.0) and $
      (jansen[jansenmatch].sii_6731_ew[1] gt 0.0),nindx)

    x = nfgs[nfgsmatch[indx]].sii_6731_ew[0]
    y = jansen[jansenmatch[indx]].sii_6731_ew[0]
    resid = y-x
;   resid = 100*(y-x)/x
    stats = im_stats(resid,sigrej=3.0,/verbose,/no_head)

    xstr = strtrim(string(stats.median_rej,format='(F12.2)'),2)+' ('+$
      strtrim(string(stats.mean_rej,format='(F12.2)'),2)+$
      '\pm'+strtrim(string(stats.sigma_rej,format='(F12.2)'),2)+')'

    xtitle = 'EW([S II] \lambda6731) [iSPEC1d]'
    ytitle = 'EW([S II] \lambda6731) [Jansen]'

    xrange = [min(x)<min(y),max(x)>max(y)]
    yrange = xrange
    yrange2 = (abs(min(resid))>abs(max(resid)))*[-1.1,1.1]
    
    djs_plot, x, y, xtitle='', ytitle=ytitle, ps=8, position=pos[*,0], $
      xrange=xrange, yrange=yrange, xstyle=3, ystyle=3, yminor=3, $
      xthick=postthick, ythick=postthick, charsize=1.2, charthick=postthick, $
      xtickname=replicate(' ',10), /xlog, /ylog
    djs_oplot, 10^!x.crange, 10^!y.crange, line=0, thick=postthick

    djs_plot, x, resid, /noerase, xtitle=xtitle, ytitle='Residuals ['+angstrom()+']', ps=8, $
      position=pos[*,2], xrange=xrange, yrange=yrange2, xstyle=3, ystyle=3, yminor=3, $
      xthick=postthick, ythick=postthick, charsize=1.2, charthick=postthick, /xlog
    djs_oplot, 10^!x.crange, [0,0], line=0, thick=postthick
    
    legend, textoidl(xstr), /right, /top, box=0, charsize=1.0, $
      charthick=postthick, clear=keyword_set(postscript)

    if not keyword_set(postscript) then cc = get_kbrd(1)

; fluxes
    
    indx = where((nfgs[nfgsmatch].sii_6731[1] gt 0.0) and $
      (jansen[jansenmatch].sii_6731[1] gt 0.0),nindx)

    x = nfgs[nfgsmatch[indx]].sii_6731[0]
    y = jansen[jansenmatch[indx]].sii_6731[0]

    scale = median(y/x)
    x = scale*x
    
    resid = 100*(y-x)/x
    stats = im_stats(resid,sigrej=3.0,/verbose,/no_head)
    
    xstr = strtrim(string(stats.median_rej,format='(F12.2)'),2)+' ('+$
      strtrim(string(stats.mean_rej,format='(F12.2)'),2)+$
      '\pm'+strtrim(string(stats.sigma_rej,format='(F12.2)'),2)+')'

    xtitle = 'F([S II] \lambda6731) [iSPEC1d]'
    ytitle = 'F([S II] \lambda6731) [Jansen]'

    xrange = [min(x)<min(y),max(x)>max(y)]
    yrange = xrange
    yrange2 = (abs(min(resid))>abs(max(resid)))*[-1.1,1.1]
    
    djs_plot, x, y, xtitle='', ytitle=ytitle, ps=8, position=pos[*,1], $
      xrange=xrange, yrange=yrange, xstyle=3, ystyle=3, yminor=3, $
      xthick=postthick, ythick=postthick, charsize=1.2, charthick=postthick, $
      xtickname=replicate(' ',10), /xlog, /ylog, /noerase
    djs_oplot, 10^!x.crange, 10^!y.crange, line=0, thick=postthick

    djs_plot, x, resid, /noerase, xtitle=xtitle, ytitle='Residuals [%]', ps=8, $
      position=pos[*,3], xrange=xrange, yrange=yrange2, xstyle=3, ystyle=3, yminor=3, $
      xthick=postthick, ythick=postthick, charsize=1.2, charthick=postthick, /xlog
    djs_oplot, 10^!x.crange, [0,0], line=0, thick=postthick
    
    legend, textoidl(xstr), /right, /top, box=0, charsize=1.0, $
      charthick=postthick, clear=keyword_set(postscript)

    if keyword_set(postscript) then dfpsclose

stop    
    
return
end
