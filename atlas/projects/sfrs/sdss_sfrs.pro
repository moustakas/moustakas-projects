;+
; NAME:
;       SDSS_PLOTS
;
; PURPOSE:
;       Investigate nebular star formation rate indicators in the
;       Tremonti SDSS sample.
;
; INPUTS:
;
; KEYWORD PARAMETERS:
;
; OUTPUTS:
;
; COMMENTS:
;
; EXAMPLES:
;
; MODIFICATION HISTORY:
;       J. Moustakas, 2004 Nov 04, U of A
;-

pro sdss_plots, sdssdust, sdssnodust, postscript=postscript, paper=paper, $
  cleanpng=cleanpng, hiiregions=hiiregions, atlas=atlas, nfgs=nfgs, $
  kewley_grids=kewley_grids, bptmodels=bptmodels, _extra=extra

; sdss_plots, sdssdust, sdssnodust, /atlas, /nfgs, /hii, /bpt, /paper, /post, /clean

    if keyword_set(paper) then postscript = 1L
    if keyword_set(postscript) then postthick = 8.0 else postthick = 2.0

    if n_elements(snrcut) eq 0L then snrcut = 3.0

    light = 2.99792458D18       ; speed of light [A/s]
    lsun = 3.826D33             ; [erg/s]

    htmlbase = 'sdss'

    html_path = atlas_path(/web)+'analysis/'
    pspath = atlas_path(/web)+'analysis/'+htmlbase+'/'
    paperpath = atlas_path(/papers)+'sfrs/FIG_SFRS/'
    paperpath2 = atlas_path(/papers)+'abundances/FIG_ABUNDANCES/'

; read the data and the models    
    
    if (n_elements(sdssdust) eq 0L) then sdssdust = read_sdss(sdssnodust=sdssnodust,/hiionly)

    sbgrids = read_kewley_grids(model=3,Z=Z,U=U)  ; model=1
    hiigrids = read_kewley_grids(model=8,Z=Z,U=U) ; model=6

    if keyword_set(atlas) then atlasdust = read_integrated(linefitnodust=atlasnodust,/snrcuts,/hiionly,/silent)
    if keyword_set(nfgs) then nfgsdust = read_nfgs(linefitnodust=nfgsnodust,/snrcuts,/hiionly,/silent)

    if keyword_set(hiiregions) then hii = read_hii_regions(/silent)
    if keyword_set(bptmodels) then models = kewley_bpt_lines(/kauffmann,_extra=extra)

    HaHb = 2.86
;   HaHb = return_tbalmer(/HaHb)
    
;   sfhgrid = read_sfhgrid()

; initialize plotting variables

    @'../plotting_ranges.idl'

    oiiharange = [-1.5,0.8] ; overwrite the default range
    
    Zsun_old = 8.9
    Zsun_new = 8.69

    atlaspsize = 0.5
    nfgspsize = 0.5
    
; determine upper and lower branches    

    divsep = 0.05
    div = findgen((8.5-7.5)/divsep+1)*divsep+7.5
    ndiv = n_elements(div)
    resid = fltarr(ndiv)

; clean up old PS and PNG files

    if keyword_set(cleanpng) then begin
       splog, 'Deleting all PNG and PS files in '+pspath
       spawn, ['rm -f '+pspath+'*.png*'], /sh
       spawn, ['rm -f '+pspath+'*.ps'], /sh
    endif

    if (not keyword_set(postscript)) then window, 0, xs=700, ys=700

; ---------------------------------------------------------------------------    
; compare the Kauffmann and the Bell masses
; ---------------------------------------------------------------------------    

    psname = 'sdss_mass_gr_r_vs_mass_kauffmann.ps'
    im_openclose, pspath+psname, postscript=postscript

    indx = where(sdssdust.mass_gr_r gt -900,nindx)

    x = sdssdust[indx].mass_gr_r
    xerr = sdssdust[indx].mass_gr_r_err
    
    y = sdssdust[indx].kauffmann_mass
    yerr = y*0.0

    residuals = x - y
    residuals_err = sqrt(xerr^2 + yerr^2)

    stats = im_stats(residuals,/verbose)
    xstr = strtrim(string(stats.median,format='(F12.2)'),2)+'\pm'+strtrim(string(stats.sigma_rej,format='(F12.2)'),2)

    xtitle = 'log (M / M'+sunsymbol()+') [(g-r) & r]'
    ytitle = 'log (M / M'+sunsymbol()+') [Kauffmann]'
    ytitle2 = 'Residuals'

    xrange = [8.0,12.5]
    yrange = xrange

    xrange2 = xrange
    yrange2 = stats.median + stats.sigma_rej*[-5,5] ; max(abs(residuals))*[-1.1,1.1]
    
    pagemaker, nx=1, ny=2, yspace=0, xspace=0, xmargin=[1.3,0.3], $
      ymargin=[0.3,1.3], position=pos, height=[6.0,3.4], /normal

    sdss_lineplot, x, y, xerr, yerr, plottype=1, postscript=postscript, $
      xtitle='', ytitle=ytitle, xrange=xrange, yrange=yrange, legendtype=0, $
      /left, /top, position=pos[*,0], xtickname=replicate(' ',10)
    djs_oplot, !x.crange, !y.crange, line=0, thick=postthick

    sdss_lineplot, x, residuals, xerr, residuals_err, plottype=1, postscript=postscript, $
      xtitle=xtitle, ytitle=ytitle2, xrange=xrange2, yrange=yrange2, legendtype=0, $
      /left, /top, position=pos[*,1], /noerase
    djs_oplot, !x.crange, [0,0], line=0, thick=postthick
    legend, textoidl(xstr), /right, /bottom, box=0, charsize=1.8, charthick=postthick

    im_openclose, postscript=postscript, /close    

; ------------------------------------------------------------
; 12+log(O/H) [M91_O32] vs 12+log(O/H) [N2]
; ------------------------------------------------------------

    psname = 'sdss_12OH_M91_O32_vs_12OH_N2.ps'
    im_openclose, pspath+psname, postscript=postscript

    offset = 0.20-0.08 ; the 0.08 shifts the N2 calibration to O3N2+N2 and the 0.20 shifts O3N2+N2 to M91
    stroffset = string(offset,format='(F4.2)')

    indx = where((sdssnodust.z_12oh_m91_o32 gt -900) and (sdssnodust.z_12oh_n2_denicolo gt -900),nindx)

    x = sdssnodust[indx].z_12oh_m91_o32
    xerr = sdssnodust[indx].z_12oh_m91_o32_err

    y = sdssnodust[indx].z_12oh_n2_denicolo+offset
    yerr = sdssnodust[indx].z_12oh_n2_denicolo_err
    
    residuals = x - y
    residuals_err = sqrt(xerr^2 + yerr^2)

    stats = im_stats(residuals,/verbose)
    xstr = strtrim(string(stats.median,format='(F12.2)'),2)+'\pm'+strtrim(string(stats.sigma_rej,format='(F12.2)'),2)

    xtitle = '12 + log (O/H) [M91 O_{32}]'
    ytitle = '12 + log (O/H) + '+stroffset+' [N_{2}]'
    ytitle2 = 'Residuals'

    xrange = [8.0,9.5]
    yrange = [8.0,9.5]

    xrange2 = xrange
    yrange2 = max(abs(residuals))*[-1.1,1.1]

    pagemaker, nx=1, ny=2, yspace=0, xspace=0, xmargin=[1.3,0.3], $
      ymargin=[0.3,1.3], position=pos, height=[6.0,3.4], /normal

    sdss_lineplot, x, y, xerr, yerr, plottype=1, postscript=postscript, $
      xtitle='', ytitle=ytitle, xrange=xrange, yrange=yrange, legendtype=0, $
      /left, /top, position=pos[*,0], xtickname=replicate(' ',10)
    djs_oplot, !x.crange, !y.crange, line=0, thick=postthick

    sdss_lineplot, x, residuals, xerr, residuals_err, plottype=1, postscript=postscript, $
      xtitle=xtitle, ytitle=ytitle2, xrange=xrange2, yrange=yrange2, legendtype=0, $
      /left, /top, position=pos[*,1], /noerase
    djs_oplot, !x.crange, [0,0], line=0, thick=postthick
    legend, textoidl(xstr), /right, /bottom, box=0, charsize=1.8, charthick=postthick

    im_openclose, postscript=postscript, /close    
    
; ------------------------------------------------------------
; 12+log(O/H) [O3N2+N2] vs 12+log(O/H) [Strong]
; ------------------------------------------------------------

    psname = 'sdss_12OH_O3N2_N2_vs_12OH_strong.ps'
    im_openclose, pspath+psname, postscript=postscript

    pagemaker, nx=2, ny=2, position=pos, /normal, xspace=0.0, $
      xmargin=[1.3,1.3], ymargin=[0.2,1.3], yspace=0.0

    xtitle = '12 + log (O/H) [O_{3} N_{2} + N_{2}]'
    ytitle = 'Residuals'

    xrange = [8.0,9.1]
    yrange = [-0.9,0.9]
;   yrange = [7.2,9.8]
    multicharsize = 1.6

; --------------------------------------------------    
; ZKH94
; --------------------------------------------------    

    indx = where((sdssnodust.z_12oh_o3n2_n2 gt 8.1) and (sdssnodust.Z_12OH_zkh94 gt -900),nindx)

    x = sdssnodust[indx].Z_12OH_O3N2_N2
    xerr = sdssnodust[indx].Z_12OH_O3N2_N2_err

    y = sdssnodust[indx].z_12oh_zkh94
    yerr = sdssnodust[indx].z_12oh_zkh94_err

    y = x-y ; NOTE!
    
    xbig = x
    ybig = y

;   ytitle = '12 + log (O/H) [N_{2} O_{2}]'

    stats = im_stats(ybig,/verbose)
;   stats = im_stats(xbig-ybig,/verbose)
    xstr = strtrim(string(stats.median,format='(F12.2)'),2)+'\pm'+strtrim(string(stats.sigma_rej,format='(F12.2)'),2)
    
    rcor = r_correlate(xbig,ybig,zd=zd,probd=probd)
    splog, 'Spearman rank test: ', rcor, zd

    sdss_lineplot, x, y, xerr, yerr, plottype=1, postscript=postscript, $
      xtitle='', ytitle=ytitle, xrange=xrange, yrange=yrange, legendtype=0, $
      xstyle=3, ystyle=3, /right, /top, position=pos[*,0], charsize=multicharsize, $
      xtickname=replicate(' ',10);, xatlas=xatlas, yatlas=yatlas, $
;     xerratlas=xerratlas, yerratlas=yerratlas, xnfgs=xnfgs, ynfgs=ynfgs, $
;     xerrnfgs=xerrnfgs, yerrnfgs=yerrnfgs
    djs_oplot, !x.crange, [0,0], line=2, thick=postthick
;   djs_oplot, !x.crange, !y.crange, thick=postthick

    legend, textoidl(xstr), /right, /bottom, box=0, charsize=multicharsize, charthick=postthick
    legend, '(a)', /left, /top, box=0, charsize=multicharsize, charthick=postthick
    legend, textoidl('ZKH94'), /right, /top, box=0, charsize=multicharsize, charthick=postthick

; --------------------------------------------------    
; M91
; --------------------------------------------------    

    indx = where((sdssnodust.z_12oh_o3n2_n2 gt 8.1) and (sdssnodust.Z_12OH_m91_upper gt -900.0),nindx)

    x = sdssnodust[indx].z_12oh_o3n2_n2
    xerr = sdssnodust[indx].z_12oh_o3n2_n2_err
    
    y = sdssnodust[indx].z_12oh_m91_upper
    yerr = sdssnodust[indx].z_12oh_m91_upper_err

    y = x-y ; NOTE!
    
    xbig = x
    ybig = y

;   ytitle = '12 + log (O/H) [M91]'

    stats = im_stats(ybig,/verbose,/no_head)
;   stats = im_stats(xbig-ybig,/verbose,/no_head)
    xstr = strtrim(string(stats.median,format='(F12.2)'),2)+'\pm'+strtrim(string(stats.sigma_rej,format='(F12.2)'),2)
    
    rcor = r_correlate(xbig,ybig,zd=zd,probd=probd)
    splog, 'Spearman rank test: ', rcor, zd

    sdss_lineplot, x, y, xerr, yerr, plottype=1, postscript=postscript, $
      xtitle='', ytitle='', xrange=xrange, yrange=yrange, legendtype=0, $
      xstyle=3, ystyle=11, /right, /top, position=pos[*,1], /noerase, charsize=multicharsize, $
      xtickname=replicate(' ',10), ytickname=replicate(' ',10)
    axis, /yaxis, ystyle=3, ythick=postthick, charsize=multicharsize, charthick=postthick, $
      yrange=yrange, ytitle=ytitle
    djs_oplot, !x.crange, [0,0], line=2, thick=postthick
;   djs_oplot, !x.crange, !y.crange, thick=postthick

    legend, textoidl(xstr), /right, /bottom, box=0, charsize=multicharsize, charthick=postthick
    legend, '(b)', /left, /top, box=0, charsize=multicharsize, charthick=postthick
    legend, textoidl('M91'), /right, /top, box=0, charsize=multicharsize, charthick=postthick
    
; --------------------------------------------------    
; KD02-Combined
; --------------------------------------------------    

    indx = where((sdssnodust.Z_12OH_O3N2_N2 gt 8.1) and (sdssnodust.z_12oh_kd02_combined gt -900.0),nindx)
 
    x = sdssnodust[indx].z_12oh_o3n2_n2
    xerr = sdssnodust[indx].z_12oh_o3n2_n2_err

    y = sdssnodust[indx].z_12oh_kd02_combined
    yerr = sdssnodust[indx].z_12oh_kd02_combined_err

    y = x-y ; NOTE!

    xbig = x
    ybig = y

;   ytitle = '12 + log (O/H) [KD02]'

    stats = im_stats(ybig,/verbose,/no_head)
;   stats = im_stats(xbig-ybig,/verbose,/no_head)
    xstr = strtrim(string(stats.median,format='(F12.2)'),2)+'\pm'+strtrim(string(stats.sigma_rej,format='(F12.2)'),2)
    
    rcor = r_correlate(xbig,ybig,zd=zd,probd=probd)
    splog, 'Spearman rank test: ', rcor, zd

    sdss_lineplot, x, y, xerr, yerr, plottype=1, postscript=postscript, $
      xtitle='', ytitle=ytitle, xrange=xrange, yrange=yrange, legendtype=0, charsize=multicharsize, $
      xstyle=3, ystyle=3, /right, /top, position=pos[*,2], /noerase;, xnfgs=xnfgs, ynfgs=ynfgs, $
;     xerrnfgs=xerrnfgs, yerrnfgs=yerrnfgs
    djs_oplot, !x.crange, [0,0], line=2, thick=postthick
;   djs_oplot, !x.crange, !y.crange, thick=postthick
 
    legend, textoidl(xstr), /right, /bottom, box=0, charsize=multicharsize, charthick=postthick
    legend, '(c)', /left, /top, box=0, charsize=multicharsize, charthick=postthick
    legend, textoidl('KD02'), /right, /top, box=0, charsize=multicharsize, charthick=postthick


; --------------------------------------------------    
; T04
; --------------------------------------------------    

    indx = where((sdssnodust.z_12oh_o3n2_n2 gt 8.1) and (sdssnodust.Z_12OH_t04 gt -900.0),nindx)

    x = sdssnodust[indx].z_12oh_o3n2_n2
    xerr = sdssnodust[indx].z_12oh_o3n2_n2_err
    
    y = sdssnodust[indx].z_12oh_t04
    yerr = sdssnodust[indx].z_12oh_t04_err

    y = x-y ; NOTE!
    
    xbig = x
    ybig = y

;   ytitle = '12 + log (O/H) [T04]'

    stats = im_stats(ybig,/verbose,/no_head)
;   stats = im_stats(xbig-ybig,/verbose,/no_head)
    xstr = strtrim(string(stats.median,format='(F12.2)'),2)+'\pm'+strtrim(string(stats.sigma_rej,format='(F12.2)'),2)
    
    rcor = r_correlate(xbig,ybig,zd=zd,probd=probd)
    splog, 'Spearman rank test: ', rcor, zd

    sdss_lineplot, x, y, xerr, yerr, plottype=1, postscript=postscript, $
      xtitle='', ytitle='', xrange=xrange, yrange=yrange, legendtype=0, $
      xstyle=3, ystyle=11, /right, /top, position=pos[*,3], /noerase, charsize=multicharsize, $
      ytickname=replicate(' ',10)
    axis, /yaxis, ystyle=3, ythick=postthick, charsize=multicharsize, charthick=postthick, $
      yrange=yrange, ytitle=ytitle
    djs_oplot, !x.crange, [0,0], line=2, thick=postthick

    legend, textoidl(xstr), /right, /bottom, box=0, charsize=multicharsize, charthick=postthick
    legend, '(d)', /left, /top, box=0, charsize=multicharsize, charthick=postthick
    legend, textoidl('T04'), /right, /top, box=0, charsize=multicharsize, charthick=postthick

; title    
    
    xyouts, 0.5*(pos[2,1]-pos[0,0])+pos[0,0], pos[1,2]*0.2, textoidl(xtitle), align=0.5, $
      charsize=multicharsize, charthick=postthick, /normal

    im_openclose, postscript=postscript, /close    

; ------------------------------------------------------------
; inter-compare empirical and strong-line calibrations
; ------------------------------------------------------------

    psname = 'sdss_12OH_strong_vs_12OH_empirical.ps'
    im_openclose, pspath+psname, postscript=postscript

    pagemaker, nx=3, ny=3, position=pos, /normal, xspace=0.0, $
      xmargin=[1.3,0.2], ymargin=[0.2,1.1], yspace=0.0

    xrange = [6.9,9.5]
    yrange = xrange
    multicharsize = 1.3

; --------------------------------------------------    
; Panel 1 - O3N2+N2 versus KD02-Combined
; --------------------------------------------------    

    indx = where((sdssnodust.Z_12OH_O3N2_N2 gt -900.0) and (sdssnodust.z_12oh_kd02_combined gt -900.0),nindx)
 
    x = sdssnodust[indx].z_12oh_o3n2_n2
    xerr = sdssnodust[indx].z_12oh_o3n2_n2_err

    y = sdssnodust[indx].z_12oh_kd02_combined
    yerr = sdssnodust[indx].z_12oh_kd02_combined_err

    xtitle = '12 + log (O/H) [O_{3} N_{2} + N_{2}]'
    ytitle = '12 + log (O/H) [KD02]'

    stats = im_stats(x-y,/verbose)
    xstr = strtrim(string(stats.median,format='(F12.2)'),2)+'\pm'+strtrim(string(stats.sigma_rej,format='(F12.2)'),2)
    
    sdss_lineplot, x, y, xerr, yerr, plottype=1, postscript=postscript, $
      xtitle='', ytitle=ytitle, xrange=xrange, yrange=yrange, legendtype=0, charsize=multicharsize, $
      xstyle=3, ystyle=3, /right, /top, position=pos[*,0], xtickname=replicate(' ',10)
    djs_oplot, !x.crange, !y.crange, thick=postthick
 
    legend, textoidl(xstr), /right, /bottom, box=0, charsize=multicharsize, charthick=postthick
    legend, '(a)', /left, /top, box=0, charsize=multicharsize, charthick=postthick
    
; --------------------------------------------------    
; Panel 2 - N2O2 versus KD02-Combined
; --------------------------------------------------    

    indx = where((sdssnodust.Z_12OH_kd02_nii_oii gt -900.0) and (sdssnodust.z_12oh_kd02_combined gt -900.0),nindx)

    x = sdssnodust[indx].z_12oh_kd02_nii_oii
    xerr = sdssnodust[indx].z_12oh_kd02_nii_oii_err

    y = sdssnodust[indx].z_12oh_kd02_combined
    yerr = sdssnodust[indx].z_12oh_kd02_combined_err

    xtitle = '12 + log (O/H) [N_{2} O_{2}]'
    ytitle = '12 + log (O/H) [KD02]'

    stats = im_stats(x-y,/verbose,/no_head)
    xstr = strtrim(string(stats.median,format='(F12.2)'),2)+'\pm'+strtrim(string(stats.sigma_rej,format='(F12.2)'),2)
    
    sdss_lineplot, x, y, xerr, yerr, plottype=1, postscript=postscript, $
      xtitle='', ytitle='', xrange=xrange, yrange=yrange, legendtype=0, charsize=multicharsize, $
      xstyle=3, ystyle=3, /right, /top, position=pos[*,1], /noerase, ytickname=replicate(' ',10), $
      xtickname=replicate(' ',10)
    djs_oplot, !x.crange, !y.crange, thick=postthick

    legend, textoidl(xstr), /right, /bottom, box=0, charsize=multicharsize, charthick=postthick
    legend, '(b)', /left, /top, box=0, charsize=multicharsize, charthick=postthick

; --------------------------------------------------    
; Panel 3 - M91 versus KD02-Combined
; --------------------------------------------------    
 
    indx = where((sdssnodust.Z_12OH_m91_upper gt -900.0) and (sdssnodust.Z_12OH_m91_lower gt -900.0) and $
      (sdssnodust.z_12oh_kd02_combined gt -900.0),nindx)

    y = sdssnodust[indx].z_12oh_kd02_combined
    yerr = sdssnodust[indx].z_12oh_kd02_combined_err

    x = y*0.0
    xerr = yerr*0.0
    
    for idiv = 0L, ndiv-1L do begin

       up = where((sdssnodust[indx].z_12oh_kd02_combined gt div[idiv]),nup,comp=lo,ncomp=nlo)
       if (nup ne 0L) then begin
          x[up] = sdssnodust[indx[up]].z_12oh_m91_upper
          xerr[up] = sdssnodust[indx[up]].z_12oh_m91_upper_err
       endif

       if (nlo ne 0L) then begin
          x[lo] = sdssnodust[indx[lo]].z_12oh_m91_lower
          xerr[lo] = sdssnodust[indx[lo]].z_12oh_m91_lower_err
       endif

       resid[idiv] = stddev(x-y)
;      resid[idiv] = djsig(x-y,sigrej=3.0)
;      plot, x, x-y, ps=4, xrange=xrange, yrange=[-1,1]
;      print, div[idiv], resid[idiv]
;      cc = get_kbrd(1)
       
    endfor

;   mindiv = find_nminima(resid,div,nfind=1,width=0.2)
    minresid = min(resid,minindx)
    mindiv = div[minindx]
    splog, 'M91 vs KD02-Combined: '+string(mindiv,format='(F4.2)')

    up = where((sdssnodust[indx].z_12oh_kd02_combined gt mindiv),nup)
    if (nup ne 0L) then begin
       x = sdssnodust[indx[up]].z_12oh_m91_upper
       xerr = sdssnodust[indx[up]].z_12oh_m91_upper_err

       y = sdssnodust[indx[up]].z_12oh_kd02_combined
       yerr = sdssnodust[indx[up]].z_12oh_kd02_combined_err
    endif

    lo = where((sdssnodust[indx].z_12oh_kd02_combined lt mindiv),nlo)
    if (nlo ne 0L) then begin
       x = [x,sdssnodust[indx[lo]].z_12oh_m91_lower]
       xerr = [xerr,sdssnodust[indx[lo]].z_12oh_m91_lower_err]

       y = [y,sdssnodust[indx[lo]].z_12oh_kd02_combined]
       yerr = [yerr,sdssnodust[indx[lo]].z_12oh_kd02_combined_err]
    endif

    xtitle = '12 + log (O/H) [M91]'
    ytitle = '12 + log (O/H) [KD02]'

    stats = im_stats(x-y,/verbose,/no_head)
    xstr = strtrim(string(stats.median,format='(F12.2)'),2)+'\pm'+strtrim(string(stats.sigma_rej,format='(F12.2)'),2)
    
    sdss_lineplot, x, y, xerr, yerr, plottype=1, postscript=postscript, $
      xtitle=xtitle, ytitle='', xrange=xrange, yrange=yrange, legendtype=0, $
      xstyle=3, ystyle=3, /right, /top, position=pos[*,2], /noerase, charsize=multicharsize, $
      ytickname=replicate(' ',10)
    djs_oplot, !x.crange, !y.crange, thick=postthick
 
    legend, textoidl(xstr), /right, /bottom, box=0, charsize=multicharsize, charthick=postthick
    legend, '(c)', /left, /top, box=0, charsize=multicharsize, charthick=postthick

; --------------------------------------------------    
; Panel 4 - O3N2+N2 versus M91
; --------------------------------------------------    

    indx = where((sdssnodust.z_12oh_o3n2_n2 gt -900) and (sdssnodust.Z_12OH_m91_upper gt -900.0) and $
      (sdssnodust.Z_12OH_m91_lower gt -900.0),nindx)

    x = sdssnodust[indx].z_12oh_o3n2_n2
    xerr = sdssnodust[indx].z_12oh_o3n2_n2_err
    
    y = x*0.0
    yerr = xerr*0.0
    
    for idiv = 0L, ndiv-1L do begin

       up = where((x gt div[idiv]),nup,comp=lo,ncomp=nlo)
       if (nup ne 0L) then begin
          y[up] = sdssnodust[indx[up]].z_12oh_m91_upper
          yerr[up] = sdssnodust[indx[up]].z_12oh_m91_upper_err
       endif

       if (nlo ne 0L) then begin
          y[lo] = sdssnodust[indx[lo]].z_12oh_m91_lower
          yerr[lo] = sdssnodust[indx[lo]].z_12oh_m91_lower_err
       endif

       resid[idiv] = stddev(x-y)
;      resid[idiv] = djsig(x-y,sigrej=3.0)
;      plot, x, x-y, ps=4, xrange=xrange, yrange=[-1,1]
;      print, div[idiv], resid[idiv]
;      cc = get_kbrd(1)
       
    endfor

;   mindiv = find_nminima(resid,div,nfind=1,width=0.2)
    minresid = min(resid,minindx)
    mindiv = div[minindx]
    splog, 'O3N2+N2 vs M91: '+string(mindiv,format='(F4.2)')

    up = where((sdssnodust[indx].z_12oh_o3n2_n2 gt mindiv),nup)
    if (nup ne 0L) then begin
       x = sdssnodust[indx[up]].z_12oh_o3n2_n2
       xerr = sdssnodust[indx[up]].z_12oh_o3n2_n2_err

       y = sdssnodust[indx[up]].z_12oh_m91_upper
       yerr = sdssnodust[indx[up]].z_12oh_m91_upper_err
    endif

    lo = where((sdssnodust[indx].z_12oh_o3n2_n2 lt mindiv),nlo)
    if (nlo ne 0L) then begin
       x = [x,sdssnodust[indx[lo]].z_12oh_o3n2_n2]
       xerr = [xerr,sdssnodust[indx[lo]].z_12oh_o3n2_n2_err]

       y = [y,sdssnodust[indx[lo]].z_12oh_m91_lower]
       yerr = [yerr,sdssnodust[indx[lo]].z_12oh_m91_lower_err]
    endif

    xtitle = '12 + log (O/H) [O_{3} N_{2} + N_{2}]'
    ytitle = '12 + log (O/H) [M91]'

    stats = im_stats(x-y,/verbose,/no_head)
    xstr = strtrim(string(stats.median,format='(F12.2)'),2)+'\pm'+strtrim(string(stats.sigma_rej,format='(F12.2)'),2)
    
    sdss_lineplot, x, y, xerr, yerr, plottype=1, postscript=postscript, $
      xtitle='', ytitle=ytitle, xrange=xrange, yrange=yrange, legendtype=0, $
      xstyle=3, ystyle=3, /right, /top, position=pos[*,3], /noerase, charsize=multicharsize, $
      xtickname=replicate(' ',10)
    djs_oplot, !x.crange, !y.crange, thick=postthick

    legend, textoidl(xstr), /right, /bottom, box=0, charsize=multicharsize, charthick=postthick
    legend, '(d)', /left, /top, box=0, charsize=multicharsize, charthick=postthick
    
; --------------------------------------------------    
; Panel 5 - N2O2 versus M91
; --------------------------------------------------    

    indx = where((sdssnodust.Z_12OH_kd02_nii_oii gt -900.0) and (sdssnodust.z_12oh_m91_upper gt -900.0) and $
      (sdssnodust.z_12oh_m91_lower gt -900.0),nindx)
    
    x = sdssnodust[indx].z_12oh_kd02_nii_oii
    xerr = sdssnodust[indx].z_12oh_kd02_nii_oii_err
    
    y = x*0.0
    yerr = xerr*0.0
    
    for idiv = 0L, ndiv-1L do begin

       up = where((x gt div[idiv]),nup,comp=lo,ncomp=nlo)
       if (nup ne 0L) then begin
          y[up] = sdssnodust[indx[up]].z_12oh_m91_upper
          yerr[up] = sdssnodust[indx[up]].z_12oh_m91_upper_err
       endif

       if (nlo ne 0L) then begin
          y[lo] = sdssnodust[indx[lo]].z_12oh_m91_lower
          yerr[lo] = sdssnodust[indx[lo]].z_12oh_m91_lower_err
       endif

       resid[idiv] = stddev(x-y)
;      resid[idiv] = djsig(x-y,sigrej=3.0)
;      plot, x, x-y, ps=4, xrange=xrange, yrange=[-1,1]
;      print, div[idiv], resid[idiv]
;      cc = get_kbrd(1)
       
    endfor

;   mindiv = find_nminima(resid,div,nfind=1,width=0.2)
    minresid = min(resid,minindx)
    mindiv = div[minindx]
    splog, 'N2O2 vs M91: '+string(mindiv,format='(F4.2)')

    up = where((sdssnodust[indx].z_12oh_kd02_nii_oii gt mindiv),nup)
    if (nup ne 0L) then begin
       x = sdssnodust[indx[up]].z_12oh_kd02_nii_oii
       xerr = sdssnodust[indx[up]].z_12oh_kd02_nii_oii_err

       y = sdssnodust[indx[up]].z_12oh_m91_upper
       yerr = sdssnodust[indx[up]].z_12oh_m91_upper_err
    endif

    lo = where((sdssnodust[indx].z_12oh_kd02_nii_oii lt mindiv),nlo)
    if (nlo ne 0L) then begin
       x = [x,sdssnodust[indx[lo]].z_12oh_kd02_nii_oii]
       xerr = [xerr,sdssnodust[indx[lo]].z_12oh_kd02_nii_oii_err]

       y = [y,sdssnodust[indx[lo]].z_12oh_m91_lower]
       yerr = [yerr,sdssnodust[indx[lo]].z_12oh_m91_lower_err]
    endif

    xtitle = '12 + log (O/H) [N_{2} O_{2}]'
    ytitle = '12 + log (O/H) [M91]'

    stats = im_stats(x-y,/verbose,/no_head)
    xstr = strtrim(string(stats.median,format='(F12.2)'),2)+'\pm'+strtrim(string(stats.sigma_rej,format='(F12.2)'),2)
    
    sdss_lineplot, x, y, xerr, yerr, plottype=1, postscript=postscript, $
      xtitle=xtitle, ytitle='', xrange=xrange, yrange=yrange, legendtype=0, charsize=multicharsize, $
      xstyle=3, ystyle=3, /right, /top, position=pos[*,4], /noerase, $
      ytickname=replicate(' ',10)
    djs_oplot, !x.crange, !y.crange, thick=postthick

    legend, textoidl(xstr), /right, /bottom, box=0, charsize=multicharsize, charthick=postthick
    legend, '(e)', /left, /top, box=0, charsize=multicharsize, charthick=postthick
    
; --------------------------------------------------    
; Panel 6 - No Data
; --------------------------------------------------    

; --------------------------------------------------    
; Panel 7 - O3N2+N2 versus N2O2
; --------------------------------------------------    

    indx = where((sdssnodust.z_12oh_o3n2_n2 gt -900) and (sdssnodust.Z_12OH_kd02_nii_oii gt -900.0),nindx)

    x = sdssnodust[indx].Z_12OH_O3N2_N2
    xerr = sdssnodust[indx].Z_12OH_O3N2_N2_err

    y = sdssnodust[indx].z_12oh_kd02_nii_oii
    yerr = sdssnodust[indx].z_12oh_kd02_nii_oii_err

    xtitle = '12 + log (O/H) [O_{3} N_{2} + N_{2}]'
    ytitle = '12 + log (O/H) [N_{2} O_{2}]'

    stats = im_stats(x-y,/verbose,/no_head)
    xstr = strtrim(string(stats.median,format='(F12.2)'),2)+'\pm'+strtrim(string(stats.sigma_rej,format='(F12.2)'),2)
    
    sdss_lineplot, x, y, xerr, yerr, plottype=1, postscript=postscript, $
      xtitle=xtitle, ytitle=ytitle, xrange=xrange, yrange=yrange, legendtype=0, $
      xstyle=3, ystyle=3, /right, /top, position=pos[*,6], /noerase, charsize=multicharsize
    djs_oplot, !x.crange, !y.crange, thick=postthick

    legend, textoidl(xstr), /right, /bottom, box=0, charsize=multicharsize, charthick=postthick
    legend, '(f)', /left, /top, box=0, charsize=multicharsize, charthick=postthick
    
; --------------------------------------------------    
; Panel 8 - No Data
; --------------------------------------------------    

; --------------------------------------------------    
; Panel 9 - No Data
; --------------------------------------------------    

    im_openclose, postscript=postscript, /close    

; ------------------------------------------------------------
; BPT diagnostic plot 1 [N II]/Ha versus [O III]/Hb
; ------------------------------------------------------------

    psname = 'sdss_niiha_vs_oiiihb.ps'
    im_openclose, pspath+psname, postscript=postscript

    lineratio, sdssdust, 'NII_6584', 'H_ALPHA', 'OIII_5007', 'H_BETA', $
      x, xerr, y, yerr, index=indx, nindex=nindx, snrcut=snrcut
    
    xtitle = 'log ([N II] \lambda6584 / H\alpha)_{obs}'
    ytitle = 'log ([O III] \lambda5007 / H\beta)_{obs}'

    xrange = niiharange
    yrange = oiiihbrange

    if keyword_set(hiiregions) then begin
       good = where((hii.nii_6584_h_alpha gt -900.0) and (hii.oiii_5007_h_beta gt -900.0))
       xregion = hii[good].nii_6584_h_alpha & xerrregion = hii[good].nii_6584_h_alpha_err
       yregion = hii[good].oiii_5007_h_beta & yerrregion = hii[good].oiii_5007_h_beta_err
    endif

    sdss_lineplot, x, y, xerr, yerr, plottype=3, postscript=postscript, $
      xtitle=xtitle, ytitle=ytitle, xrange=xrange, yrange=yrange, legendtype=0, $
      /right, /top;, xregion=xregion, yregion=yregion, xerrregion=xerrregion, yerrregion=yerrregion
    
; overplot L. Kewley's starburst mixing line

    if keyword_set(bptmodels) then begin
    
       oplot, models.x_nii, models.y_nii, line=0, thick=2.0
;      oplot, models.x_nii, models.y_nii_upper, line=2, thick=2.0
;      oplot, models.x_nii, models.y_nii_lower, line=2, thick=2.0

    endif
    
    im_openclose, postscript=postscript, /close

; ------------------------------------------------------------
; 2-panel f(3727) - E(B-V)
; ------------------------------------------------------------
    
    psname = 'sdss_f3727_ehbha_multipanel.ps'
    im_openclose, pspath+psname, postscript=postscript

    indx = where((sdssnodust.h_alpha[0]/sdssnodust.h_alpha[1] gt snrcut) and $   
      (sdssdust.oii_3727_continuum[0]/sdssdust.oii_3727_continuum[1] gt snrcut) and $
      (sdssnodust.oii_3727_continuum[0]/sdssnodust.oii_3727_continuum[1] gt snrcut) and $
      (sdssnodust.ebv_hahb_err gt 0.0),nindx)

    x = sdssnodust[indx].ehbha
    xerr = sdssnodust[indx].ehbha_err
    xebv = x/(k_lambda(4861,/odonnell)-k_lambda(6563,/odonnell))
    
    ha = sdssnodust[indx].h_alpha[0]
    ha_err = sdssnodust[indx].h_alpha[1]
    
    xrange = ehbharange
    yrange = f3727harange

    xtitle = 'E(H\beta-H\alpha)'

    pagemaker, nx=1, ny=2, yspace=0, xspace=0, xmargin=[1.3,0.4], $
      ymargin=[1.1,1.1], position=pos, /normal

; ##########################################################
; Panel 1: No reddening correction
; ##########################################################

    f3727 = sdssdust[indx].oii_3727_continuum[0]
    f3727_err = sdssdust[indx].oii_3727_continuum[1]
    
    y = alog10(f3727/ha)
    yerr = im_compute_error(f3727,f3727_err,ha,ha_err,/log)

    junk = im_stats(y,/verbose)

    ytitle = 'log (F_{\lambda3727}_{obs} / H\alpha_{cor})'

    sdss_lineplot, x, y, xerr, yerr, plottype=3, postscript=postscript, $
      xtitle='', ytitle=ytitle, xrange=xrange, yrange=yrange, legendtype=0, $
      charsize=1.5, position=pos[*,0], xtickname=replicate(' ',10), $
      xstyle=11, ymargin=[4,3]
    axis, /xaxis, xrange=interpol(xebv,x,!x.crange), xthick=postthick, $
      charthick=postthick, charsize=1.5, xtitle = 'E(B-V)', xsty=1

    legend, 'a)', /left, /top, box=0, charsize=1.3, charthick=postthick

    label = ['No Reddening Correction']
    legend, textoidl(label), /left, /bottom, box=0, charsize=1.5, $
      charthick=postthick
    
; ############################################################
; Panel 2: Good Reddening Correction
; ############################################################

    f3727 = sdssnodust[indx].oii_3727_continuum[0]
    f3727_err = sdssnodust[indx].oii_3727_continuum[1]

    y = alog10(f3727/ha)
    yerr = im_compute_error(f3727,f3727_err,ha,ha_err,/log)

    junk = im_stats(y,/verbose,/no_head)

    ytitle = 'log (F_{\lambda3727}_{cor} / H\alpha_{cor})'

    sdss_lineplot, x, y, xerr, yerr, plottype=3, postscript=postscript, $
      xtitle=xtitle, ytitle=ytitle, xrange=xrange, yrange=yrange, legendtype=0, $
      charsize=1.5, position=pos[*,1], /noerase

    legend, 'c)', /left, /top, box=0, charsize=1.3, charthick=postthick

    label = ['Individual Reddening Correction']
    legend, textoidl(label), /left, /bottom, box=0, charsize=1.5, $
      charthick=postthick

    im_openclose, postscript=postscript, /close

; ------------------------------------------------------------
; [O III]/[O II] vs E(B-V)
; ------------------------------------------------------------

    psname = 'sdss_oiiioii_vs_ehbha.ps'
    im_openclose, pspath+psname, postscript=postscript

    lineratio, sdssnodust, 'OIII_5007', 'OII_3727', '', '', x, xerr, $
      dum1, dum2, index=indx, nindex=nindx, snrcut=snrcut

    y = sdssnodust[indx].ehbha
    yerr = sdssnodust[indx].ehbha_err
    yebv = y/(k_lambda(4861,/odonnell)-k_lambda(6563,/odonnell))

    xtitle = 'log ([O III] \lambda5007 / [O II] \lambda3727)_{cor}'
    ytitle = 'E(H\beta-H\alpha)'

    xrange = oiiioiirange
    yrange = ehbharange

    sdss_lineplot, x, y, xerr, yerr, plottype=3, postscript=postscript, $
      xtitle=xtitle, ytitle=ytitle, xrange=xrange, yrange=yrange, legendtype=0, $
      ystyle=11, xmargin=[8,6], ymargin=[4,3], /right, /top
    axis, /yaxis, yrange=interpol(yebv,y,!y.crange), ythick=postthick, $
      charthick=postthick, charsize=2.0, ytitle = 'E(B-V)', ysty=1

    im_openclose, postscript=postscript, /close        

; ------------------------------------------------------------
; 4-panel H-beta plot - L(B) vs hbha_cor
; ------------------------------------------------------------
    
    psname = 'sdss_LB_vs_hbha_cor_multipanel.ps'
    im_openclose, pspath+psname, postscript=postscript

    indx = where((sdssdust.h_beta_uncor[0]/sdssdust.h_beta_uncor[1] gt snrcut) and $      
      (sdssdust.B_lum gt -900) and (sdssnodust.h_alpha[0]/sdssnodust.h_alpha[1] gt snrcut) and $
      (sdssnodust.ebv_hahb_err gt 0.0)$
      ,nindx)

    x = sdssdust[indx].B_lum
    xabs = sdssdust[indx].M_B
    xerr = sdssdust[indx].B_lum_err
    
    ha = sdssnodust[indx].h_alpha[0]
    ha_err = sdssnodust[indx].h_alpha[1]

    hahb_true = sdssdust[indx].h_alpha[0]/sdssdust[indx].h_beta[0]
    hahb_uncor = sdssdust[indx].h_alpha_uncor[0]/sdssdust[indx].h_beta_uncor[0]

    hb_ew_abs     = 3.5 ; EW [mean absorption correction]
    hb_ew_abs_err = 1.0

    xrange = LBrange
    yrange = [-2.8,0.45]

    xtitle = 'log [L(B) / L(B'+sunsymbol()+')]'
    ytitle = 'log (H\beta / H\alpha_{cor})'

    pagemaker, nx=2, ny=2, yspace=0, xspace=0, xmargin=[1.3,0.4], $
      ymargin=[0.9,1.1], position=pos, /normal

; ##########################################################
; Panel 1: No absorption correction, no reddening correction
; ##########################################################

    hb = sdssdust[indx].h_beta_uncor[0]
    hb_err = sdssdust[indx].h_beta_uncor[1]
    
    y = alog10(hb/ha)
    yerr = im_compute_error(hb,hb_err,ha,ha_err,/log)

    junk = im_stats(y,/verbose)

    sdss_lineplot, x, y, xerr, yerr, plottype=3, postscript=postscript, $
      xtitle='', ytitle=ytitle, xrange=xrange, yrange=yrange, legendtype=0, $
      charsize=1.4, xstyle=11, xtickname=replicate(' ',10), position=pos[*,0]
    axis, /xaxis, xrange=interpol(xabs,x,!x.crange), xthick=postthick, charthick=postthick, $
      charsize=1.4, xtitle = textoidl('M_{B}'), xsty=1
    djs_oplot, !x.crange, alog10(1.0/HaHb)*[1,1], line=2, thick=postthick

    legend, 'a)', /left, /top, box=0, charsize=1.3, charthick=postthick

    label = ['No Abs. Correction','No Dust Correction']
    legend, textoidl(label), /left, /bottom, box=0, charsize=1.3, $
      spacing=0, charthick=postthick

; ############################################################
; Panel 2: Mean absorption correction, no reddening correction
; ############################################################

    hb = sdssdust[indx].h_beta_uncor[0] + hb_ew_abs*sdssdust[indx].babs_h_beta_continuum[0]
    hb_err = sqrt(sdssdust[indx].h_beta[1]^2 + (hb_ew_abs_err*sdssdust[indx].babs_h_beta_continuum[0])^2)
    
    y = alog10(hb/ha)
    yerr = im_compute_error(hb,hb_err,ha,ha_err,/log)

    junk = im_stats(y,/verbose,/no_head)

    sdss_lineplot, x, y, xerr, yerr, plottype=3, postscript=postscript, $
      xtitle='', ytitle='', xrange=xrange, yrange=yrange, legendtype=0, $
      charsize=1.4, xstyle=11, xtickname=replicate(' ',10), ytickname=replicate(' ',10), $
      position=pos[*,1], /noerase
    axis, /xaxis, xrange=interpol(xabs,x,!x.crange), xthick=postthick, charthick=postthick, $
      charsize=1.4, xtitle = textoidl('M_{B}'), xsty=1
    djs_oplot, !x.crange, alog10(1.0/HaHb)*[1,1], line=2, thick=postthick

    legend, 'b)', /left, /top, box=0, charsize=1.3, charthick=postthick

    label = ['Mean Abs. Correction','No Dust Correction']
    legend, textoidl(label), /left, /bottom, box=0, charsize=1.3, $
      spacing=0, charthick=postthick
    
; ############################################################
; Panel 3: Good absorption correction, no reddening correction
; ############################################################

    hb = sdssdust[indx].h_beta[0]
    hb_err = sdssdust[indx].h_beta[1]
    
    y = alog10(hb/ha)
    yerr = im_compute_error(hb,hb_err,ha,ha_err,/log)

    junk = im_stats(y,/verbose,/no_head)
    
    sdss_lineplot, x, y, xerr, yerr, plottype=3, postscript=postscript, $
      xtitle=xtitle, ytitle=ytitle, xrange=xrange, yrange=yrange, legendtype=0, $
      charsize=1.4, position=pos[*,2], /noerase
    djs_oplot, !x.crange, alog10(1.0/HaHb)*[1,1], line=2, thick=postthick

    legend, 'c)', /left, /top, box=0, charsize=1.3, charthick=postthick

    label = ['Individual Abs. Correction','No Dust Correction']
    legend, textoidl(label), /left, /bottom, box=0, charsize=1.3, $
      spacing=0, charthick=postthick
    
; ###############################################################
; Panel 4: Mean absorption correction, Mean reddening correction
; ###############################################################

    Aha = 1.0
    Aha_err = 0.1

    yha = sdssdust[indx].h_alpha[0]*10^(0.4*Aha)
    yha_err = sdssdust[indx].h_alpha[1]*10^(0.4*Aha) + $
      sdssdust[indx].h_alpha[0]*alog(10.0)*10^(0.4*Aha)*Aha_err

    Ahb = Aha*k_lambda(4861,/odonnell)/k_lambda(6563,/odonnell)
    Ahb_err = 0.1

    mhb = sdssdust[indx].h_beta_uncor[0] + hb_ew_abs*sdssdust[indx].babs_h_beta_continuum[0]
    mhb_err = sqrt(sdssdust[indx].h_beta[1]^2 + (hb_ew_abs_err*sdssdust[indx].babs_h_beta_continuum[0])^2)

    hb = mhb*10^(0.4*Ahb)
    hb_err = mhb_err*10^(0.4*Ahb) + mhb*alog(10.0)*10^(0.4*Ahb)*Ahb_err
    
    y = alog10(hb/yha)
    yerr = im_compute_error(hb,hb_err,yha,yha_err,/log)
       
    junk = im_stats(y,/verbose,/no_head)

    sdss_lineplot, x, y, xerr, yerr, plottype=3, postscript=postscript, $
      xtitle=xtitle, ytitle='', xrange=xrange, yrange=yrange, legendtype=0, $
      charsize=1.4, ytickname=replicate(' ',10), position=pos[*,3], /noerase
    djs_oplot, !x.crange, alog10(1.0/HaHb)*[1,1], line=2, thick=postthick

    legend, 'd)', /left, /top, box=0, charsize=1.3, charthick=postthick

    label = ['Mean Abs. Correction','Mean Dust Correction']
    legend, textoidl(label), /left, /bottom, box=0, charsize=1.3, $
      spacing=0, charthick=postthick

    im_openclose, postscript=postscript, /close

; ------------------------------------------------------------
; 4-panel H-beta plot - [O II] & Ha_cor
; ------------------------------------------------------------
    
    psname = 'sdss_hbha_cor_ewoii_multipanel.ps'
    im_openclose, pspath+psname, postscript=postscript

    indx = where((sdssdust.h_beta_uncor[0]/sdssdust.h_beta_uncor[1] gt snrcut) and $      
      (sdssdust.oii_3727_ew[0]/sdssdust.oii_3727_ew[1] gt snrcut) and $
      (sdssnodust.h_alpha[0]/sdssnodust.h_alpha[1] gt snrcut) and $
      (sdssnodust.ebv_hahb_err gt 0.0)$
      ,nindx)

    x = alog10(sdssdust[indx].oii_3727_ew[0])
    xerr = sdssdust[indx].oii_3727_ew[1]/sdssdust[indx].oii_3727_ew[0]/alog(10)
    
    ha = sdssnodust[indx].h_alpha[0]
    ha_err = sdssnodust[indx].h_alpha[1]

    hahb_true = sdssdust[indx].h_alpha[0]/sdssdust[indx].h_beta[0]
    hahb_uncor = sdssdust[indx].h_alpha_uncor[0]/sdssdust[indx].h_beta_uncor[0]

    xrange = ewoiirange
    yrange = [-3.2,0.4]

    xtitle = 'log EW([O II] \lambda3727)'
    ytitle = 'log (H\beta / H\alpha_{cor})'

    pagemaker, nx=2, ny=2, yspace=0, xspace=0, xmargin=[1.3,0.4], $
      ymargin=[0.3,1.1], position=pos, /normal

; ##########################################################
; Panel 1: No absorption correction, no reddening correction
; ##########################################################

    hb = sdssdust[indx].h_beta_uncor[0]
    hb_err = sdssdust[indx].h_beta_uncor[1]
    
    y = alog10(hb/ha)
    yerr = im_compute_error(hb,hb_err,ha,ha_err,/log)

    junk = im_stats(y,/verbose)

    sdss_lineplot, x, y, xerr, yerr, plottype=3, postscript=postscript, $
      xtitle='', ytitle=ytitle, xrange=xrange, yrange=yrange, legendtype=0, $
      charsize=1.4, xtickname=replicate(' ',10), position=pos[*,0]
    djs_oplot, !x.crange, alog10(1.0/HaHb)*[1,1], line=2, thick=postthick

    legend, 'a)', /left, /top, box=0, charsize=1.3, charthick=postthick

    label = ['No Absorption','No Reddening']
    legend, textoidl(label), /left, /bottom, box=0, charsize=1.3, $
      spacing=1.7, charthick=postthick

; ############################################################
; Panel 2: Mean absorption correction, no reddening correction
; ############################################################

    hb_ew_abs     = 3.5 ; EW [mean absorption]
    hb_ew_abs_err = 1.0

    hb = sdssdust[indx].h_beta_uncor[0] + hb_ew_abs*sdssdust[indx].babs_h_beta_continuum[0]
    hb_err = sqrt(sdssdust[indx].h_beta[1]^2 + (hb_ew_abs_err*sdssdust[indx].babs_h_beta_continuum[0])^2)
    
    y = alog10(hb/ha)
    yerr = im_compute_error(hb,hb_err,ha,ha_err,/log)

    junk = im_stats(y,/verbose,/no_head)

    sdss_lineplot, x, y, xerr, yerr, plottype=3, postscript=postscript, $
      xtitle='', ytitle='', xrange=xrange, yrange=yrange, legendtype=0, $
      charsize=1.4, xtickname=replicate(' ',10), ytickname=replicate(' ',10), $
      position=pos[*,1], /noerase
    djs_oplot, !x.crange, alog10(1.0/HaHb)*[1,1], line=2, thick=postthick

    legend, 'b)', /left, /top, box=0, charsize=1.3, charthick=postthick

    label = ['Mean Absorption','No Reddening']
    legend, textoidl(label), /left, /bottom, box=0, charsize=1.3, $
      spacing=1.7, charthick=postthick
    
; ############################################################
; Panel 3: Good absorption correction, no reddening correction
; ############################################################

    hb = sdssdust[indx].h_beta[0]
    hb_err = sdssdust[indx].h_beta[1]
    
    y = alog10(hb/ha)
    yerr = im_compute_error(hb,hb_err,ha,ha_err,/log)

    junk = im_stats(y,/verbose,/no_head)
    
    sdss_lineplot, x, y, xerr, yerr, plottype=3, postscript=postscript, $
      xtitle=xtitle, ytitle=ytitle, xrange=xrange, yrange=yrange, legendtype=0, $
      charsize=1.4, position=pos[*,2], /noerase
    djs_oplot, !x.crange, alog10(1.0/HaHb)*[1,1], line=2, thick=postthick

    legend, 'c)', /left, /top, box=0, charsize=1.3, charthick=postthick

    label = ['Good Absorption','No Reddening']
    legend, textoidl(label), /left, /bottom, box=0, charsize=1.3, $
      spacing=1.7, charthick=postthick
    
; ###############################################################
; Panel 4: Good absorption correction, Mean reddening correction
; ###############################################################

    Aha = 1.0
    Aha_err = 0.1

    yha = sdssdust[indx].h_alpha[0]*10^(0.4*Aha)
    yha_err = sdssdust[indx].h_alpha[1]*10^(0.4*Aha) + $
      sdssdust[indx].h_alpha[0]*alog(10.0)*10^(0.4*Aha)*Aha_err

    Ahb = Aha*k_lambda(4861)/k_lambda(6563)
    Ahb_err = 0.1

    hb = sdssdust[indx].h_beta[0]*10^(0.4*Ahb)
    hb_err = sdssdust[indx].h_beta[1]*10^(0.4*Ahb) + $
      sdssdust[indx].h_beta[0]*alog(10.0)*10^(0.4*Ahb)*Ahb_err
    
    y = alog10(hb/yha)
    yerr = im_compute_error(hb,hb_err,yha,yha_err,/log)
       
    junk = im_stats(y,/verbose,/no_head)

    sdss_lineplot, x, y, xerr, yerr, plottype=3, postscript=postscript, $
      xtitle=xtitle, ytitle='', xrange=xrange, yrange=yrange, legendtype=0, $
      charsize=1.4, ytickname=replicate(' ',10), position=pos[*,3], /noerase
    djs_oplot, !x.crange, alog10(1.0/HaHb)*[1,1], line=2, thick=postthick

    legend, 'd)', /left, /top, box=0, charsize=1.3, charthick=postthick

    label = ['Good Absorption','Mean Reddening']
    legend, textoidl(label), /left, /bottom, box=0, charsize=1.3, $
      spacing=1.7, charthick=postthick

    im_openclose, postscript=postscript, /close

; --------------------------------------------------    
; GENERATE THE WEB PAGE
; --------------------------------------------------    
    
    if keyword_set(postscript) then begin

       im_ps2html, htmlbase, html_path=html_path, cleanpng=0, npscols=2, _extra=extra
       
    endif

; --------------------------------------------------    
; SELECT PLOTS FOR THE PAPER HERE
; --------------------------------------------------    

    if keyword_set(paper) then begin

       splog, 'Writing paper plots to '+paperpath+'.'
       paperplots = [$
         'sdss_LB_vs_oiiha_multipanel.ps','sdss_ehbha_vs_oiiha_multipanel.ps',$
         'sdss_niiha_vs_oiiha_cor.ps','sdss_O3N2_vs_oiiha_cor.ps',$
         'sdss_LB_vs_LU_LHa_multipanel.ps','sdss_ewoii_vs_ehbha.ps','sdss_LB_vs_ehbha.ps',$
         'sdss_LOII_vs_ebv.ps']
         
       for k = 0L, n_elements(paperplots)-1L do $
         spawn, ['/bin/cp -f '+html_path+htmlbase+'/'+paperplots[k]+' '+paperpath], /sh

       splog, 'Writing paper plots to '+paperpath2+'.'
       paperplots = [$
         'sdss_12OH_M91_O32_vs_12OH_N2.ps',$         
         'sdss_12OH_O3N2_N2_vs_12OH_strong.ps']
       for k = 0L, n_elements(paperplots)-1L do $
         spawn, ['/bin/cp -f '+html_path+htmlbase+'/'+paperplots[k]+' '+paperpath2], /sh
       
    endif

stop    
       
; --------------------------------------------------    
; [N II]/[O II]
; --------------------------------------------------    

    indx = where((sdssnodust.z_12oh_o3n2_n2 gt -900) and (sdssnodust.Z_12OH_kd02_nii_oii gt -900.0),nindx)

    x = sdssnodust[indx].Z_12OH_O3N2_N2
    xerr = sdssnodust[indx].Z_12OH_O3N2_N2_err

    y = sdssnodust[indx].z_12oh_kd02_nii_oii
    yerr = sdssnodust[indx].z_12oh_kd02_nii_oii_err

    y = x-y ; NOTE!
    
    xbig = x
    ybig = y

    if keyword_set(atlas) and keyword_set(nfgs) then begin
       
       indxatlas = where((atlasnodust.z_12oh_o3n2_n2 gt -900) and (atlasnodust.Z_12OH_kd02_nii_oii gt -900.0),nindxatlas)

       xatlas = atlasnodust[indxatlas].Z_12OH_O3N2_N2
       xerratlas = atlasnodust[indxatlas].Z_12OH_O3N2_N2_err

       yatlas = atlasnodust[indxatlas].z_12oh_kd02_nii_oii
       yerratlas = atlasnodust[indxatlas].z_12oh_kd02_nii_oii_err

       xbig = [xbig,xatlas]
       ybig = [ybig,yatlas]
       
       indxnfgs = where((nfgsnodust.z_12oh_o3n2_n2 gt -900) and (nfgsnodust.Z_12OH_kd02_nii_oii gt -900.0),nindxnfgs)

       xnfgs = nfgsnodust[indxnfgs].Z_12OH_O3N2_N2
       xerrnfgs = nfgsnodust[indxnfgs].Z_12OH_O3N2_N2_err

       ynfgs = nfgsnodust[indxnfgs].z_12oh_kd02_nii_oii
       yerrnfgs = nfgsnodust[indxnfgs].z_12oh_kd02_nii_oii_err

       xbig = [xbig,xnfgs]
       ybig = [ybig,ynfgs]
       
    endif
    
;   ytitle = '12 + log (O/H) [N_{2} O_{2}]'

    stats = im_stats(ybig,/verbose)
;   stats = im_stats(xbig-ybig,/verbose)
    xstr = strtrim(string(stats.median,format='(F12.2)'),2)+'\pm'+strtrim(string(stats.sigma_rej,format='(F12.2)'),2)
    
    rcor = r_correlate(xbig,ybig,zd=zd,probd=probd)
    splog, 'Spearman rank test: ', rcor, zd

    sdss_lineplot, x, y, xerr, yerr, plottype=1, postscript=postscript, $
      xtitle='', ytitle=ytitle, xrange=xrange, yrange=yrange, legendtype=0, $
      xstyle=3, ystyle=3, /right, /top, position=pos[*,0], charsize=multicharsize, $
      xtickname=replicate(' ',10);, xatlas=xatlas, yatlas=yatlas, $
;     xerratlas=xerratlas, yerratlas=yerratlas, xnfgs=xnfgs, ynfgs=ynfgs, $
;     xerrnfgs=xerrnfgs, yerrnfgs=yerrnfgs
    djs_oplot, !x.crange, [0,0], line=2, thick=postthick
;   djs_oplot, !x.crange, !y.crange, thick=postthick

    legend, textoidl(xstr), /right, /bottom, box=0, charsize=multicharsize, charthick=postthick
    legend, '(a)', /left, /top, box=0, charsize=multicharsize, charthick=postthick
    legend, textoidl('N_{2} O_{2}'), /right, /top, box=0, charsize=multicharsize, charthick=postthick

; --------------------------------------------------    
; M91
; --------------------------------------------------    

    indx = where((sdssnodust.z_12oh_o3n2_n2 gt -900) and (sdssnodust.Z_12OH_m91_upper gt -900.0) and $
      (sdssnodust.Z_12OH_m91_lower gt -900.0),nindx)

    x = sdssnodust[indx].z_12oh_o3n2_n2
    xerr = sdssnodust[indx].z_12oh_o3n2_n2_err
    
    y = x*0.0
    yerr = xerr*0.0
    
    for idiv = 0L, ndiv-1L do begin

       up = where((x gt div[idiv]),nup,comp=lo,ncomp=nlo)
       if (nup ne 0L) then begin
          y[up] = sdssnodust[indx[up]].z_12oh_m91_upper
          yerr[up] = sdssnodust[indx[up]].z_12oh_m91_upper_err
       endif

       if (nlo ne 0L) then begin
          y[lo] = sdssnodust[indx[lo]].z_12oh_m91_lower
          yerr[lo] = sdssnodust[indx[lo]].z_12oh_m91_lower_err
       endif

       resid[idiv] = stddev(x-y)
;      resid[idiv] = djsig(x-y,sigrej=3.0)
;      plot, x, x-y, ps=4, xrange=xrange, yrange=[-1,1]
;      print, div[idiv], resid[idiv]
;      cc = get_kbrd(1)
       
    endfor

;   mindiv = find_nminima(resid,div,nfind=1,width=0.2)
    minresid = min(resid,minindx)
    mindiv = div[minindx]
    splog, 'O3N2+N2 vs M91: '+string(mindiv,format='(F4.2)')

    up = where((sdssnodust[indx].z_12oh_o3n2_n2 gt mindiv),nup)
    if (nup ne 0L) then begin
       x = sdssnodust[indx[up]].z_12oh_o3n2_n2
       xerr = sdssnodust[indx[up]].z_12oh_o3n2_n2_err

       y = sdssnodust[indx[up]].z_12oh_m91_upper
       yerr = sdssnodust[indx[up]].z_12oh_m91_upper_err
    endif

    lo = where((sdssnodust[indx].z_12oh_o3n2_n2 lt mindiv),nlo)
    if (nlo ne 0L) then begin
       x = [x,sdssnodust[indx[lo]].z_12oh_o3n2_n2]
       xerr = [xerr,sdssnodust[indx[lo]].z_12oh_o3n2_n2_err]

       y = [y,sdssnodust[indx[lo]].z_12oh_m91_lower]
       yerr = [yerr,sdssnodust[indx[lo]].z_12oh_m91_lower_err]
    endif

    y = x-y ; NOTE!
    
    xbig = x
    ybig = y

;   if keyword_set(nfgs) then begin
;
;      indxnfgs = where((nfgsnodust.z_12oh_o3n2_n2 gt -900) and (nfgsnodust.Z_12OH_m91_upper gt -900.0) and $
;        (nfgsnodust.Z_12OH_m91_lower gt -900.0),nindxnfgs)
;
;      xnfgs = nfgsnodust[indxnfgs].z_12oh_o3n2_n2
;      xnfgserr = nfgsnodust[indxnfgs].z_12oh_o3n2_n2_err
;      
;      ynfgs = xnfgs*0.0
;      ynfgserr = xnfgserr*0.0
;      
;      for idiv = 0L, ndiv-1L do begin
;
;         up = where((xnfgs gt div[idiv]),nup,comp=lo,ncomp=nlo)
;         if (nup ne 0L) then begin
;            ynfgs[up] = nfgsnodust[indxnfgs[up]].z_12oh_m91_upper
;            ynfgserr[up] = nfgsnodust[indxnfgs[up]].z_12oh_m91_upper_err
;         endif
;
;         if (nlo ne 0L) then begin
;            ynfgs[lo] = nfgsnodust[indxnfgs[lo]].z_12oh_m91_lower
;            ynfgserr[lo] = nfgsnodust[indxnfgs[lo]].z_12oh_m91_lower_err
;         endif
;
;         resid[idiv] = stddev(xnfgs-ynfgs)
;
;      endfor
;
;      minresid = min(resid,minindxnfgs)
;      mindiv = div[minindxnfgs]
;      splog, 'O3N2+N2 vs M91 [NFGS]: '+string(mindiv,format='(F4.2)')
;
;      up = where((nfgsnodust[indxnfgs].z_12oh_o3n2_n2 gt mindiv),nup)
;      if (nup ne 0L) then begin
;         xnfgs = nfgsnodust[indxnfgs[up]].z_12oh_o3n2_n2
;         xnfgserr = nfgsnodust[indxnfgs[up]].z_12oh_o3n2_n2_err
;
;         ynfgs = nfgsnodust[indxnfgs[up]].z_12oh_m91_upper
;         ynfgserr = nfgsnodust[indxnfgs[up]].z_12oh_m91_upper_err
;      endif
;
;      lo = where((nfgsnodust[indxnfgs].z_12oh_o3n2_n2 lt mindiv),nlo)
;      if (nlo ne 0L) then begin
;         xnfgs = [xnfgs,nfgsnodust[indxnfgs[lo]].z_12oh_o3n2_n2]
;         xnfgserr = [xnfgserr,nfgsnodust[indxnfgs[lo]].z_12oh_o3n2_n2_err]
;
;         ynfgs = [ynfgs,nfgsnodust[indxnfgs[lo]].z_12oh_m91_lower]
;         ynfgserr = [ynfgserr,nfgsnodust[indxnfgs[lo]].z_12oh_m91_lower_err]
;      endif
;
;      xbig = [xbig,xnfgs]
;      ybig = [ybig,ynfgs]
;      
;   endif

;   ytitle = '12 + log (O/H) [M91]'

    stats = im_stats(ybig,/verbose,/no_head)
;   stats = im_stats(xbig-ybig,/verbose,/no_head)
    xstr = strtrim(string(stats.median,format='(F12.2)'),2)+'\pm'+strtrim(string(stats.sigma_rej,format='(F12.2)'),2)
    
    rcor = r_correlate(xbig,ybig,zd=zd,probd=probd)
    splog, 'Spearman rank test: ', rcor, zd

    sdss_lineplot, x, y, xerr, yerr, plottype=1, postscript=postscript, $
      xtitle='', ytitle=ytitle, xrange=xrange, yrange=yrange, legendtype=0, $
      xstyle=3, ystyle=3, /right, /top, position=pos[*,1], /noerase, charsize=multicharsize, $
      xtickname=replicate(' ',10);, xnfgs=xnfgs, ynfgs=ynfgs, $
;     xerrnfgs=xerrnfgs, yerrnfgs=yerrnfgs
    djs_oplot, !x.crange, [0,0], line=2, thick=postthick
;   djs_oplot, !x.crange, !y.crange, thick=postthick

    legend, textoidl(xstr), /right, /bottom, box=0, charsize=multicharsize, charthick=postthick
    legend, '(b)', /left, /top, box=0, charsize=multicharsize, charthick=postthick
    legend, textoidl('M91'), /right, /top, box=0, charsize=multicharsize, charthick=postthick
    
; --------------------------------------------------    
; KD02-Combined
; --------------------------------------------------    

    indx = where((sdssnodust.Z_12OH_O3N2_N2 gt -900.0) and (sdssnodust.z_12oh_kd02_combined gt -900.0),nindx)
 
    x = sdssnodust[indx].z_12oh_o3n2_n2
    xerr = sdssnodust[indx].z_12oh_o3n2_n2_err

    y = sdssnodust[indx].z_12oh_kd02_combined
    yerr = sdssnodust[indx].z_12oh_kd02_combined_err

    y = x-y ; NOTE!

    xbig = x
    ybig = y

;   if keyword_set(nfgs) then begin
;
;      indxnfgs = where((nfgsnodust.Z_12OH_O3N2_N2 gt -900.0) and (nfgsnodust.z_12oh_kd02_combined gt -900.0),nindxnfgs)
;      
;      xnfgs = nfgsdust[indxnfgs].z_12oh_o3n2_n2
;      xerrnfgs = nfgsdust[indxnfgs].z_12oh_o3n2_n2_err
;
;      ynfgs = nfgsdust[indxnfgs].z_12oh_kd02_combined
;      yerrnfgs = nfgsdust[indxnfgs].z_12oh_kd02_combined_err
;
;      xbig = [xbig,xnfgs]
;      ybig = [ybig,ynfgs]
;      
;   endif
       
;   ytitle = '12 + log (O/H) [KD02]'

    stats = im_stats(ybig,/verbose,/no_head)
;   stats = im_stats(xbig-ybig,/verbose,/no_head)
    xstr = strtrim(string(stats.median,format='(F12.2)'),2)+'\pm'+strtrim(string(stats.sigma_rej,format='(F12.2)'),2)
    
    rcor = r_correlate(xbig,ybig,zd=zd,probd=probd)
    splog, 'Spearman rank test: ', rcor, zd

    sdss_lineplot, x, y, xerr, yerr, plottype=1, postscript=postscript, $
      xtitle=xtitle, ytitle=ytitle, xrange=xrange, yrange=yrange, legendtype=0, charsize=multicharsize, $
      xstyle=3, ystyle=3, /right, /top, position=pos[*,2], /noerase;, xnfgs=xnfgs, ynfgs=ynfgs, $
;     xerrnfgs=xerrnfgs, yerrnfgs=yerrnfgs
    djs_oplot, !x.crange, [0,0], line=2, thick=postthick
;   djs_oplot, !x.crange, !y.crange, thick=postthick
 
    legend, textoidl(xstr), /right, /bottom, box=0, charsize=multicharsize, charthick=postthick
    legend, '(c)', /left, /top, box=0, charsize=multicharsize, charthick=postthick
    legend, textoidl('KD02'), /right, /top, box=0, charsize=multicharsize, charthick=postthick

stop    
    
    im_openclose, postscript=postscript, /close    

return
end    
