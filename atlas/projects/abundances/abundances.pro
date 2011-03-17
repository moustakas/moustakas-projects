;+
; NAME:
;       ABUNDANCES
;
; CALLING SEQUENCE:
;
; PURPOSE:
;       Investigate the nebular line ratio sequences of galaxies.
;       Derive integrated nebular abundances.
;
; COMMENTS:
;       Currently, black-background images are not supported with
;       SDSS_LINEPLOT. 
;
; MODIFICATION HISTORY:
;       J. Moustakas, 2004-2005, U of A
;-

pro abundances, hii, intdust, intnodust, sdssdust, sdssnodust, sdssancillary, $
  postscript=postscript, paper=paper, cleanpng=cleanpng, encapsulated=encapsulated, $
  _extra=extra

; abundances,intdust,intnodust,sdssdust,sdssnodust,sdssancillary,hii
    
; read the data and the models    
    
    if (n_elements(hii) eq 0L) then hii = read_hii_regions(/nosdss)
;   if (n_elements(intdust) eq 0L) then intdust = read_integrated_abundances_sample(intnodust=intnodust)
;   if (n_elements(sdssdust) eq 0L) then sdssdust = read_sdss_abundances_sample(sdssnodust=sdssnodust,sdssancillary=sdssancillary)

    models = kewley_bpt_lines(/kauffmann,_extra=extra)

    htmlbase = 'abundances'
    html_path = atlas_path(/web)+'analysis/'
    pspath = atlas_path(/web)+'analysis/'+htmlbase+'/'
    paperpath = atlas_path(/papers)+'abundances/FIG_ABUNDANCES/'

; clean up old PS and PNG files

    if keyword_set(cleanpng) then begin
       splog, 'Deleting all PNG and PS files in '+pspath
       spawn, ['rm -f '+pspath+'*.png*'], /sh
       spawn, ['rm -f '+pspath+'*.*ps'], /sh
    endif

; initialize plotting variables
    
    if keyword_set(paper) then begin
       postscript = 1L
       encapsulated = 1L
    endif

    if keyword_set(postscript) then begin
       postthick = 5.0 
    endif else begin
       postthick = 2.0
       im_window, 0, xratio=0.6, /square
    endelse

    intcolor = 'grey'
    intpsize = 0.5
    
    @'xyrange_abundances'

; ###########################################################################
; Paper Plots    
; ###########################################################################

    
    
    
    
; ------------------------------------------------------------
; inter-compare empirical and strong-line calibrations
; ------------------------------------------------------------

    psname = 'hii_12oh_strong_vs_12oh_empirical'
    im_openclose, pspath+psname, postscript=postscript, xsize=8.5, ysize=8.5, encapsulated=encapsulated

    pagemaker, nx=3, ny=3, yspace=0, xspace=0, width=2.4*[1,1,1], height=2.4*[1,1,1], $
      xmargin=[1.0,0.3], ymargin=[0.3,1.0], xpage=8.5, ypage=8.5, position=pos, /normal

    xrange = ohrange
    yrange = xrange

; --------------------------------------------------    
; Panel 1 - KK04 vs T04
; --------------------------------------------------    

    indx = where((intnodust.zstrong_12oh_t04 gt -900.0) and (intnodust.zstrong_12oh_kk04_r23_upper gt -900.0),nindx)
 
    x = intnodust[indx].zstrong_12oh_kk04_r23_upper
    xerr = intnodust[indx].zstrong_12oh_kk04_r23_upper_err

    y = intnodust[indx].zstrong_12oh_t04
    yerr = intnodust[indx].zstrong_12oh_t04_err

    ytitle = '12 + log (O/H) [T04]'

    stats = im_stats(y-x,/verbose)
    xstr = strtrim(string(stats.median,format='(F12.2)'),2)+'\pm'+strtrim(string(stats.sigma_rej,format='(F12.2)'),2)
    
    atlas1d_lineplot, x, y, xerr, yerr, plottype=1, postscript=postscript, $
      xtitle='', ytitle=ytitle, xrange=xrange, yrange=yrange, legendtype=0, charsize=charsize_2, $
      xstyle=3, ystyle=3, /right, /top, position=pos[*,0], xtickname=replicate(' ',10), $
      atlascolor=intcolor, atlaspsize=intpsize
    djs_oplot, !x.crange, !y.crange, thick=postthick
 
    legend, textoidl(xstr), /right, /bottom, box=0, charsize=charsize_2, charthick=postthick
    legend, '(a)', /left, /top, box=0, charsize=charsize_2, charthick=postthick

; --------------------------------------------------    
; Panel 2 - ZKH94 vs T04
; --------------------------------------------------    

    indx = where((intnodust.zstrong_12oh_zkh94 gt -900.0) and (intnodust.zstrong_12oh_t04 gt -900.0),nindx)

    x = intnodust[indx].zstrong_12oh_zkh94
    xerr = intnodust[indx].zstrong_12oh_zkh94_err

    y = intnodust[indx].zstrong_12oh_t04
    yerr = intnodust[indx].zstrong_12oh_t04_err

    stats = im_stats(y-x,/verbose,/no_head)
    xstr = strtrim(string(stats.median,format='(F12.2)'),2)+'\pm'+strtrim(string(stats.sigma_rej,format='(F12.2)'),2)
    
    atlas1d_lineplot, x, y, xerr, yerr, plottype=1, postscript=postscript, $
      xtitle='', ytitle='', xrange=xrange, yrange=yrange, legendtype=0, charsize=charsize_2, $
      xstyle=3, ystyle=3, /right, /top, position=pos[*,1], /noerase, ytickname=replicate(' ',10), $
      xtickname=replicate(' ',10), atlascolor=intcolor, atlaspsize=intpsize
    djs_oplot, !x.crange, !y.crange, thick=postthick, color=djs_icolor(talkcolor)

    legend, textoidl(xstr), /right, /bottom, box=0, charsize=charsize_2, charthick=postthick
    legend, '(b)', /left, /top, box=0, charsize=charsize_2, charthick=postthick, textcolor=djs_icolor(talkcolor)

; --------------------------------------------------    
; Panel 3 - M91 vs T04
; --------------------------------------------------    
 
    indx = where((intnodust.zstrong_12oh_m91_upper gt -900.0) and (intnodust.zstrong_12oh_t04 gt -900.0),nindx)

    x = intnodust[indx].zstrong_12oh_m91_upper
    xerr = intnodust[indx].zstrong_12oh_m91_upper_err

    y = intnodust[indx].zstrong_12oh_t04
    yerr = intnodust[indx].zstrong_12oh_t04_err

    xtitle = '12 + log (O/H) M91'

    stats = im_stats(y-x,/verbose,/no_head)
    xstr = strtrim(string(stats.median,format='(F12.2)'),2)+'\pm'+strtrim(string(stats.sigma_rej,format='(F12.2)'),2)
    
    atlas1d_lineplot, x, y, xerr, yerr, plottype=1, postscript=postscript, $
      xtitle=xtitle, ytitle='', xrange=xrange, yrange=yrange, legendtype=0, $
      xstyle=3, ystyle=3, /right, /top, position=pos[*,2], /noerase, charsize=charsize_2, $
      ytickname=replicate(' ',10), atlascolor=intcolor, atlaspsize=intpsize
    djs_oplot, !x.crange, !y.crange, thick=postthick, color=djs_icolor(talkcolor)
 
    legend, textoidl(xstr), /right, /bottom, box=0, charsize=charsize_2, charthick=postthick, $
      textcolor=djs_icolor(talkcolor)
    legend, '(c)', /left, /top, box=0, charsize=charsize_2, charthick=postthick, textcolor=djs_icolor(talkcolor)

; --------------------------------------------------    
; Panel 4 - KK04 vs M91
; --------------------------------------------------    

    indx = where((intnodust.zstrong_12oh_kk04_r23_upper gt -900.0) and (intnodust.zstrong_12oh_m91_upper gt -900.0),nindx)

    
    
    x = intnodust[indx].zstrong_12oh_kk04
    xerr = intnodust[indx].zstrong_12oh_kk04_err

    y = intnodust[indx].zstrong_12oh_m91
    yerr = intnodust[indx].zstrong_12oh_m91_err

    ytitle = '12 + log (O/H) [M91]'

    stats = im_stats(y-x,/verbose,/no_head)
    xstr = strtrim(string(stats.median,format='(F12.2)'),2)+'\pm'+strtrim(string(stats.sigma_rej,format='(F12.2)'),2)

stop    
    
    indx = where((intnodust.zstrong_12oh_oiiinii_pettini gt -900) and (intnodust.zstrong_12oh_m91_upper gt -900.0) and $
      (intnodust.zstrong_12oh_m91_lower gt -900.0),nindx)

    x = intnodust[indx].zstrong_12oh_oiiinii_pettini
    xerr = intnodust[indx].zstrong_12oh_oiiinii_pettini_err
    
    y = x*0.0
    yerr = xerr*0.0
    
    for idiv = 0L, ndiv-1L do begin

       up = where((x gt div[idiv]),nup,comp=lo,ncomp=nlo)
       if (nup ne 0L) then begin
          y[up] = intnodust[indx[up]].zstrong_12oh_m91_upper
          yerr[up] = intnodust[indx[up]].zstrong_12oh_m91_upper_err
       endif

       if (nlo ne 0L) then begin
          y[lo] = intnodust[indx[lo]].zstrong_12oh_m91_lower
          yerr[lo] = intnodust[indx[lo]].zstrong_12oh_m91_lower_err
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
    splog, 'oiiinii+niiha vs M91: '+string(mindiv,format='(F4.2)')

    up = where((intnodust[indx].zstrong_12oh_oiiinii_pettini gt mindiv),nup)
    if (nup ne 0L) then begin
       x = intnodust[indx[up]].zstrong_12oh_oiiinii_pettini
       xerr = intnodust[indx[up]].zstrong_12oh_oiiinii_pettini_err

       y = intnodust[indx[up]].zstrong_12oh_m91_upper
       yerr = intnodust[indx[up]].zstrong_12oh_m91_upper_err
    endif

    lo = where((intnodust[indx].zstrong_12oh_oiiinii_pettini lt mindiv),nlo)
    if (nlo ne 0L) then begin
       x = [x,intnodust[indx[lo]].zstrong_12oh_oiiinii_pettini]
       xerr = [xerr,intnodust[indx[lo]].zstrong_12oh_oiiinii_pettini_err]

       y = [y,intnodust[indx[lo]].zstrong_12oh_m91_lower]
       yerr = [yerr,intnodust[indx[lo]].zstrong_12oh_m91_lower_err]
    endif

    xtitle = '12 + log (O/H) ([O III]/H\beta)/([N II]/H\alpha)'
    ytitle = '12 + log (O/H) M91'

    stats = im_stats(y-x,/verbose,/no_head)
    xstr = strtrim(string(stats.median,format='(F12.2)'),2)+'\pm'+strtrim(string(stats.sigma_rej,format='(F12.2)'),2)
    
    atlas1d_lineplot, x, y, xerr, yerr, plottype=1, postscript=postscript, $
      xtitle='', ytitle=ytitle, xrange=xrange, yrange=yrange, legendtype=0, $
      xstyle=3, ystyle=3, /right, /top, position=pos[*,3], /noerase, charsize=charsize_2, $
      xtickname=replicate(' ',10)
    djs_oplot, !x.crange, !y.crange, thick=postthick, color=djs_icolor(talkcolor)

    legend, textoidl(xstr), /right, /bottom, box=0, charsize=charsize_2, charthick=postthick, $
      textcolor=djs_icolor(talkcolor)
    legend, '(d)', /left, /top, box=0, charsize=charsize_2, charthick=postthick, textcolor=djs_icolor(talkcolor)
    
; --------------------------------------------------    
; Panel 5 - niihaO2 versus M91
; --------------------------------------------------    

    indx = where((intnodust.zstrong_12oh_kd02_niioii gt -900.0) and (intnodust.zstrong_12oh_m91_upper gt -900.0) and $
      (intnodust.zstrong_12oh_m91_lower gt -900.0),nindx)
    
    x = intnodust[indx].zstrong_12oh_kd02_niioii
    xerr = intnodust[indx].zstrong_12oh_kd02_niioii_err
    
    y = x*0.0
    yerr = xerr*0.0
    
    for idiv = 0L, ndiv-1L do begin

       up = where((x gt div[idiv]),nup,comp=lo,ncomp=nlo)
       if (nup ne 0L) then begin
          y[up] = intnodust[indx[up]].zstrong_12oh_m91_upper
          yerr[up] = intnodust[indx[up]].zstrong_12oh_m91_upper_err
       endif

       if (nlo ne 0L) then begin
          y[lo] = intnodust[indx[lo]].zstrong_12oh_m91_lower
          yerr[lo] = intnodust[indx[lo]].zstrong_12oh_m91_lower_err
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
    splog, 'niihaO2 vs M91: '+string(mindiv,format='(F4.2)')

    up = where((intnodust[indx].zstrong_12oh_kd02_niioii gt mindiv),nup)
    if (nup ne 0L) then begin
       x = intnodust[indx[up]].zstrong_12oh_kd02_niioii
       xerr = intnodust[indx[up]].zstrong_12oh_kd02_niioii_err

       y = intnodust[indx[up]].zstrong_12oh_m91_upper
       yerr = intnodust[indx[up]].zstrong_12oh_m91_upper_err
    endif

    lo = where((intnodust[indx].zstrong_12oh_kd02_niioii lt mindiv),nlo)
    if (nlo ne 0L) then begin
       x = [x,intnodust[indx[lo]].zstrong_12oh_kd02_niioii]
       xerr = [xerr,intnodust[indx[lo]].zstrong_12oh_kd02_niioii_err]

       y = [y,intnodust[indx[lo]].zstrong_12oh_m91_lower]
       yerr = [yerr,intnodust[indx[lo]].zstrong_12oh_m91_lower_err]
    endif

    xtitle = '12 + log (O/H) [N II]/[O II]'
    ytitle = '12 + log (O/H) M91'

    stats = im_stats(y-x,/verbose,/no_head)
    xstr = strtrim(string(stats.median,format='(F12.2)'),2)+'\pm'+strtrim(string(stats.sigma_rej,format='(F12.2)'),2)
    
    atlas1d_lineplot, x, y, xerr, yerr, plottype=1, postscript=postscript, $
      xtitle=xtitle, ytitle='', xrange=xrange, yrange=yrange, legendtype=0, charsize=charsize_2, $
      xstyle=3, ystyle=3, /right, /top, position=pos[*,4], /noerase, $
      ytickname=replicate(' ',10)
    djs_oplot, !x.crange, !y.crange, thick=postthick, color=djs_icolor(talkcolor)

    legend, textoidl(xstr), /right, /bottom, box=0, charsize=charsize_2, charthick=postthick, $
      textcolor=djs_icolor(talkcolor)
    legend, '(e)', /left, /top, box=0, charsize=charsize_2, charthick=postthick, textcolor=djs_icolor(talkcolor)
    
; --------------------------------------------------    
; Panel 6 - No Data
; --------------------------------------------------    

; --------------------------------------------------    
; Panel 7 - oiiinii+niiha versus niihaO2
; --------------------------------------------------    

    indx = where((intnodust.zstrong_12oh_oiiinii_pettini gt -900) and (intnodust.zstrong_12oh_kd02_niioii gt -900.0),nindx)

    x = intnodust[indx].zstrong_12oh_oiiinii_pettini
    xerr = intnodust[indx].zstrong_12oh_oiiinii_pettini_err

    y = intnodust[indx].zstrong_12oh_kd02_niioii
    yerr = intnodust[indx].zstrong_12oh_kd02_niioii_err

    xtitle = '12 + log (O/H) ([O III]/H\beta)/([N II]/H\alpha)'
    ytitle = '12 + log (O/H) [N II]/[O II]'

    stats = im_stats(y-x,/verbose,/no_head)
    xstr = strtrim(string(stats.median,format='(F12.2)'),2)+'\pm'+strtrim(string(stats.sigma_rej,format='(F12.2)'),2)
    
    atlas1d_lineplot, x, y, xerr, yerr, plottype=1, postscript=postscript, $
      xtitle=xtitle, ytitle=ytitle, xrange=xrange, yrange=yrange, legendtype=0, $
      xstyle=3, ystyle=3, /right, /top, position=pos[*,6], /noerase, charsize=charsize_2
    djs_oplot, !x.crange, !y.crange, thick=postthick, color=djs_icolor(talkcolor)

    legend, textoidl(xstr), /right, /bottom, box=0, charsize=charsize_2, charthick=postthick, $
      textcolor=djs_icolor(talkcolor)
    legend, '(f)', /left, /top, box=0, charsize=charsize_2, charthick=postthick, textcolor=djs_icolor(talkcolor)
    
; --------------------------------------------------    
; Panel 8 - No Data
; --------------------------------------------------    

; --------------------------------------------------    
; Panel 9 - No Data
; --------------------------------------------------    

    im_openclose, postscript=postscript, /close    

stop    
    
;;    y = intnodust[indx].zstrong_12oh_kd02_combined
;;    yerr = intnodust[indx].zstrong_12oh_kd02_combined_err
;;
;;    x = y*0.0
;;    xerr = yerr*0.0
;;    
;;    for idiv = 0L, ndiv-1L do begin
;;
;;       up = where((intnodust[indx].zstrong_12oh_kd02_combined gt div[idiv]),nup,comp=lo,ncomp=nlo)
;;       if (nup ne 0L) then begin
;;          x[up] = intnodust[indx[up]].zstrong_12oh_m91_upper
;;          xerr[up] = intnodust[indx[up]].zstrong_12oh_m91_upper_err
;;       endif
;;
;;       if (nlo ne 0L) then begin
;;          x[lo] = intnodust[indx[lo]].zstrong_12oh_m91_lower
;;          xerr[lo] = intnodust[indx[lo]].zstrong_12oh_m91_lower_err
;;       endif
;;
;;       resid[idiv] = stddev(x-y)
;;;      resid[idiv] = djsig(x-y,sigrej=3.0)
;;;      plot, x, x-y, ps=4, xrange=xrange, yrange=[-1,1]
;;;      print, div[idiv], resid[idiv]
;;;      cc = get_kbrd(1)
;;       
;;    endfor
;;
;;;   mindiv = find_nminima(resid,div,nfind=1,width=0.2)
;;    minresid = min(resid,minindx)
;;    mindiv = div[minindx]
;;    splog, 'M91 vs KD02-Combined: '+string(mindiv,format='(F4.2)')
;;
;;    up = where((intnodust[indx].zstrong_12oh_kd02_combined gt mindiv),nup)
;;    if (nup ne 0L) then begin
;;       x = intnodust[indx[up]].zstrong_12oh_m91_upper
;;       xerr = intnodust[indx[up]].zstrong_12oh_m91_upper_err
;;
;;       y = intnodust[indx[up]].zstrong_12oh_kd02_combined
;;       yerr = intnodust[indx[up]].zstrong_12oh_kd02_combined_err
;;    endif
;;
;;    lo = where((intnodust[indx].zstrong_12oh_kd02_combined lt mindiv),nlo)
;;    if (nlo ne 0L) then begin
;;       x = [x,intnodust[indx[lo]].zstrong_12oh_m91_lower]
;;       xerr = [xerr,intnodust[indx[lo]].zstrong_12oh_m91_lower_err]
;;
;;       y = [y,intnodust[indx[lo]].zstrong_12oh_kd02_combined]
;;       yerr = [yerr,intnodust[indx[lo]].zstrong_12oh_kd02_combined_err]
;;    endif

; ------------------------------------------------------------
; 12+log(O/H) [KK04] vs 12+log(O/H) [Empirical] - Integrated
; ------------------------------------------------------------

    psname = 'int_12oh_kk04_vs_12oh_empirical'
    im_openclose, pspath+psname, postscript=postscript, encapsulated=encapsulated, xsize=8.5, ysize=10.5

    pagemaker, nx=1, ny=3, height=3.0*[1,1,1], width=6.8, xmargin=[1.3,0.4], $
      ymargin=[0.4,1.1], xspace=0, yspace=0, xpage=8.5, ypage=10.5, $
      position=pos, /normal

    xtitle = '12+log(O/H) [KK04]'
    ytitle = '12+log(O/H)'

    xrange = ohrange
    yrange = xrange

; --------------------------------------------------    
; [N II]/Ha
; --------------------------------------------------    

    indx = where((intnodust.zstrong_12oh_niiha_pettini gt -900) and $
      (intnodust.zstrong_12oh_kk04 gt -900),nindx)

    x = intnodust[indx].zstrong_12oh_kk04
    xerr = intnodust[indx].zstrong_12oh_kk04_err

    y = intnodust[indx].zstrong_12oh_niiha_pettini
    yerr = intnodust[indx].zstrong_12oh_niiha_pettini_err

    residuals = y-x
    stats = im_stats(residuals)
    xstr = strtrim(string(stats.median_rej,format='(F12.2)'),2)+' ('+$
      strtrim(string(stats.mean_rej,format='(F12.2)'),2)+$
      '\pm'+strtrim(string(stats.sigma_rej,format='(F12.2)'),2)+')'

    atlas1d_lineplot, x, y, xerr, yerr, plottype=1, postscript=postscript, $
      xtitle='', ytitle=ytitle, xrange=xrange, yrange=yrange, legendtype=0, $
      /right, /top, position=pos[*,0], charsize=charsize_8, xtickname=replicate(' ',10)
    djs_oplot, !x.crange, !y.crange, thick=postthick, color=djs_icolor(talkcolor), line=2

    legend, textoidl(xstr), /right, /bottom, box=0, charsize=charsize_5, $
      charthick=postthick, textcolor=djs_icolor(talkcolor)
    legend, textoidl('[N II]/H\alpha'), /left, /top, box=0, $
      charsize=charsize_8, charthick=postthick, textcolor=djs_icolor(talkcolor)

; --------------------------------------------------    
; ([O III]/Hb)/([N II]/Ha)
; --------------------------------------------------    

    indx = where((intnodust.zstrong_12oh_oiiinii_pettini gt -900) and $
      (intnodust.zstrong_12oh_kk04 gt -900),nindx)

    x = intnodust[indx].zstrong_12oh_kk04
    xerr = intnodust[indx].zstrong_12oh_kk04_err

    y = intnodust[indx].zstrong_12oh_oiiinii_pettini
    yerr = intnodust[indx].zstrong_12oh_oiiinii_pettini_err

    residuals = y-x
    stats = im_stats(residuals)
    xstr = strtrim(string(stats.median_rej,format='(F12.2)'),2)+' ('+$
      strtrim(string(stats.mean_rej,format='(F12.2)'),2)+$
      '\pm'+strtrim(string(stats.sigma_rej,format='(F12.2)'),2)+')'

    atlas1d_lineplot, x, y, xerr, yerr, plottype=1, postscript=postscript, $
      xtitle='', ytitle=ytitle, xrange=xrange, yrange=yrange, legendtype=0, /noerase, $
      /right, /top, position=pos[*,1], charsize=charsize_8, xtickname=replicate(' ',10)
    djs_oplot, !x.crange, !y.crange, thick=postthick, color=djs_icolor(talkcolor), line=2

    legend, textoidl(xstr), /right, /bottom, box=0, charsize=charsize_5, $
      charthick=postthick, textcolor=djs_icolor(talkcolor)
    legend, 'O3N2', /left, /top, box=0, $
;   legend, textoidl('([O III]/H\beta)/([N II]/H\alpha)'), /left, /top, box=0, $
      charsize=charsize_8, charthick=postthick, textcolor=djs_icolor(talkcolor)
;   legend, textoidl('(a) ([O III]/H\beta)/([N II]/H\alpha)'), /left, /top, box=0, $
;     charsize=charsize_8, charthick=postthick, textcolor=djs_icolor(talkcolor)

; --------------------------------------------------    
; P-method
; --------------------------------------------------    

    indx = where((intnodust.zstrong_12oh_kk04 gt -900) and (intnodust.zstrong_12oh_pt05 gt -900.0),nindx)

    x = intnodust[indx].zstrong_12oh_kk04
    xerr = intnodust[indx].zstrong_12oh_kk04_err

    y = intnodust[indx].zstrong_12oh_pt05
    yerr = intnodust[indx].zstrong_12oh_pt05_err

    residuals = y-x
    stats = im_stats(residuals)
    xstr = strtrim(string(stats.median_rej,format='(F12.2)'),2)+' ('+$
      strtrim(string(stats.mean_rej,format='(F12.2)'),2)+$
      '\pm'+strtrim(string(stats.sigma_rej,format='(F12.2)'),2)+')'

    atlas1d_lineplot, x, y, xerr, yerr, plottype=1, postscript=postscript, $
      xtitle=xtitle, ytitle=ytitle, xrange=xrange, yrange=yrange, legendtype=0, $
      /right, /top, position=pos[*,2], charsize=charsize_8, /noerase
    djs_oplot, !x.crange, !y.crange, thick=postthick, color=djs_icolor(talkcolor), line=2

    legend, textoidl(xstr), /right, /bottom, box=0, charsize=charsize_5, $
      charthick=postthick, textcolor=djs_icolor(talkcolor)
    legend, textoidl('PT05'), /left, /top, box=0, $
      charsize=charsize_8, charthick=postthick, textcolor=djs_icolor(talkcolor)
;   legend, textoidl('(c) P-method'), /left, /top, box=0, $
;     charsize=charsize_8, charthick=postthick, textcolor=djs_icolor(talkcolor)

    im_openclose, postscript=postscript, /close    

stop    
    
; ------------------------------------------------------------
; 12+log(O/H) [M91] vs 12+log(O/H) [T04]
; ------------------------------------------------------------

    psname = 'sdss_12oh_m91_vs_12oh_t04'
    im_openclose, pspath+psname, postscript=postscript, xsize=8.5, ysize=9.0, encapsulated=encapsulated

    pagemaker, nx=1, ny=2, xspace=0, yspace=0.0, width=6.5, height=[5.0,2.5], $
      xmargin=[1.5,0.5], ymargin=[0.4,1.1], xpage=8.5, ypage=9.0, $
      position=pos, /normal

    indx = where((sdssnodust.zstrong_12oh_m91 gt -900) and (sdssnodust.zstrong_12oh_t04 gt -900.0),nindx)
;   indx = where((sdssnodust.zstrong_12oh_m91 gt -900) and (sdssancillary.tremonti_oh gt -900.0),nindx)

    x = sdssnodust[indx].zstrong_12oh_m91
    xerr = sdssnodust[indx].zstrong_12oh_m91_err

    y = sdssnodust[indx].zstrong_12oh_t04
    yerr = sdssnodust[indx].zstrong_12oh_t04_err
;   y = sdssancillary[indx].tremonti_oh
;   yerr = sdssancillary[indx].tremonti_oh_err

    yresid = y-x
    stats = im_stats(yresid,/verbose)
    xstr = strtrim(string(stats.median_rej,format='(F12.2)'),2)+' ('+$
      strtrim(string(stats.mean_rej,format='(F12.2)'),2)+$
      '\pm'+strtrim(string(stats.sigma_rej,format='(F12.2)'),2)+')'
    
    xtitle = '12 + log (O/H) [M91]'
    ytitle = '12 + log (O/H) [T04]'
    yresidtitle = '\Delta'+'log(O/H)'

    xrange = ohrange11
    yrange = xrange
    yresidrange = 0.9*[-1,1]

    sdss_lineplot, x, y, xerr, yerr, plottype=1, postscript=postscript, $
      xtitle='', ytitle=ytitle, xrange=xrange, yrange=yrange, legendtype=0, $
      /right, /top, position=pos[*,0], xtickname=replicate(' ',10), charsize=charsize_8
    djs_oplot, !x.crange, !y.crange, line=0, thick=postthick

    sdss_lineplot, x, yresid, xerr, yerr, plottype=1, postscript=postscript, $
      xtitle=xtitle, ytitle=yresidtitle, xrange=xrange, yrange=yresidrange, $
      legendtype=0, /right, /top, /noerase, position=pos[*,1], charsize=charsize_8
    djs_oplot, !x.crange, [0,0], line=0, thick=postthick

    legend, textoidl(xstr), /right, /top, box=0, charsize=charsize_5, charthick=postthick

    im_openclose, postscript=postscript, /close    

; ------------------------------------------------------------
; 12+log(O/H) [T04] vs 12+log(O/H) [Strong] - SDSS
; ------------------------------------------------------------

    psname = 'sdss_12oh_t04_vs_12oh_strong'
    im_openclose, pspath+psname, postscript=postscript, encapsulated=encapsulated, xsize=8.5, ysize=10.5

    pagemaker, nx=1, ny=3, height=3.0*[1,1,1], width=6.8, xmargin=[1.3,0.4], $
      ymargin=[0.4,1.1], xspace=0, yspace=0, xpage=8.5, ypage=10.5, $
      position=pos, /normal

;   if keyword_set(dotalk) then polyfill, [1,1,0,0,1], [1,0,0,1,1], /normal, $
;     color=djs_icolor('black')

    xtitle = '12+log(O/H) [T04]'
    ytitle = '12+log(O/H)'

    xrange = ohrange11
    yrange = xrange

; --------------------------------------------------    
; McGaugh (1991)
; --------------------------------------------------    

    indx = where((sdssnodust.zstrong_12oh_m91 gt -900) and (sdssnodust.zstrong_12oh_t04 gt -900.0),nindx)
;   indx = where((sdssnodust.zstrong_12oh_m91 gt -900) and (sdssancillary.tremonti_oh gt -900.0),nindx)

    x = sdssnodust[indx].zstrong_12oh_t04
    xerr = sdssnodust[indx].zstrong_12oh_t04_err
;   x = sdssancillary[indx].tremonti_oh
;   xerr = sdssancillary[indx].tremonti_oh_err
 
    y = sdssnodust[indx].zstrong_12oh_m91
    yerr = sdssnodust[indx].zstrong_12oh_m91_err

    residuals = y-x
    stats = im_stats(residuals)
    xstr = strtrim(string(stats.median_rej,format='(F12.2)'),2)+' ('+$
      strtrim(string(stats.mean_rej,format='(F12.2)'),2)+$
      '\pm'+strtrim(string(stats.sigma_rej,format='(F12.2)'),2)+')'
    xstr2 = strtrim(string(stats.mean_rej,format='(F12.2)'),2)+$
      '\pm'+strtrim(string(stats.sigma_rej,format='(F12.2)'),2)

    sdss_lineplot, x, y, xerr, yerr, plottype=1, postscript=postscript, $
      xtitle='', ytitle=ytitle, xrange=xrange, yrange=yrange, legendtype=0, $
      /right, /top, position=pos[*,0], charsize=charsize_8, xtickname=replicate(' ',10)
    djs_oplot, !x.crange, !y.crange, thick=postthick, color=djs_icolor(talkcolor), line=2

    legend, textoidl(xstr2), /right, /bottom, box=0, charsize=charsize_5, $
;   legend, textoidl(xstr), /right, /bottom, box=0, charsize=charsize_5, $
      charthick=postthick, textcolor=djs_icolor(talkcolor)
    legend, textoidl('M91'), /left, /top, box=0, $
;   legend, textoidl('McGaugh (1991, Upper Branch)'), /left, /top, box=0, $
      charsize=charsize_8, charthick=postthick, textcolor=djs_icolor(talkcolor)

; --------------------------------------------------    
; Zaritsky et al. (1994)
; --------------------------------------------------    

    indx = where((sdssnodust.zstrong_12oh_zkh94 gt -900) and (sdssnodust.zstrong_12oh_t04 gt -900.0),nindx)
;   indx = where((sdssnodust.zstrong_12oh_zkh94 gt -900) and (sdssancillary.tremonti_oh gt -900.0),nindx)

    x = sdssnodust[indx].zstrong_12oh_t04
    xerr = sdssnodust[indx].zstrong_12oh_t04_err
;   x = sdssancillary[indx].tremonti_oh
;   xerr = sdssancillary[indx].tremonti_oh_err
 
    y = sdssnodust[indx].zstrong_12oh_zkh94
    yerr = sdssnodust[indx].zstrong_12oh_zkh94_err

    residuals = y-x
    stats = im_stats(residuals)
    xstr = strtrim(string(stats.median_rej,format='(F12.2)'),2)+' ('+$
      strtrim(string(stats.mean_rej,format='(F12.2)'),2)+$
      '\pm'+strtrim(string(stats.sigma_rej,format='(F12.2)'),2)+')'
    xstr2 = strtrim(string(stats.mean_rej,format='(F12.2)'),2)+$
      '\pm'+strtrim(string(stats.sigma_rej,format='(F12.2)'),2)

    sdss_lineplot, x, y, xerr, yerr, plottype=1, postscript=postscript, $
      xtitle='', ytitle=ytitle, xrange=xrange, yrange=yrange, legendtype=0, $
      /right, /top, position=pos[*,1], charsize=charsize_8, xtickname=replicate(' ',10), /noerase
    djs_oplot, !x.crange, !y.crange, thick=postthick, color=djs_icolor(talkcolor), line=2

    legend, textoidl(xstr2), /right, /bottom, box=0, charsize=charsize_5, $
;   legend, textoidl(xstr), /right, /bottom, box=0, charsize=charsize_5, $
      charthick=postthick, textcolor=djs_icolor(talkcolor)
    legend, textoidl('ZKH94'), /left, /top, box=0, $
;   legend, textoidl('Zaritsky et al. (1994)'), /left, /top, box=0, $
      charsize=charsize_8, charthick=postthick, textcolor=djs_icolor(talkcolor)

; --------------------------------------------------    
; Kobulnicky & Kewley (2004)
; --------------------------------------------------    

    indx = where((sdssnodust.zstrong_12oh_kk04 gt -900) and (sdssnodust.zstrong_12oh_t04 gt -900.0),nindx)
;   indx = where((sdssnodust.zstrong_12oh_kk04 gt -900) and (sdssancillary.tremonti_oh gt -900.0),nindx)

    x = sdssnodust[indx].zstrong_12oh_t04
    xerr = sdssnodust[indx].zstrong_12oh_t04_err
;   x = sdssancillary[indx].tremonti_oh
;   xerr = sdssancillary[indx].tremonti_oh_err
 
    y = sdssnodust[indx].zstrong_12oh_kk04
    yerr = sdssnodust[indx].zstrong_12oh_kk04_err

    residuals = y-x
    stats = im_stats(residuals)
    xstr = strtrim(string(stats.median_rej,format='(F12.2)'),2)+' ('+$
      strtrim(string(stats.mean_rej,format='(F12.2)'),2)+$
      '\pm'+strtrim(string(stats.sigma_rej,format='(F12.2)'),2)+')'
    xstr2 = strtrim(string(stats.mean_rej,format='(F12.2)'),2)+$
      '\pm'+strtrim(string(stats.sigma_rej,format='(F12.2)'),2)

    sdss_lineplot, x, y, xerr, yerr, plottype=1, postscript=postscript, $
      xtitle=xtitle, ytitle=ytitle, xrange=xrange, yrange=yrange, legendtype=0, $
      /right, /top, position=pos[*,2], charsize=charsize_8, /noerase
    djs_oplot, !x.crange, !y.crange, thick=postthick, color=djs_icolor(talkcolor), line=2

    legend, textoidl(xstr2), /right, /bottom, box=0, charsize=charsize_5, $
;   legend, textoidl(xstr), /right, /bottom, box=0, charsize=charsize_5, $
      charthick=postthick, textcolor=djs_icolor(talkcolor)
    legend, textoidl('KK04'), /left, /top, box=0, $
;   legend, textoidl('Kobulnicky & Kewley (2004, Upper Branch)'), /left, /top, box=0, $
      charsize=charsize_8, charthick=postthick, textcolor=djs_icolor(talkcolor)

    im_openclose, postscript=postscript, /close    

; ------------------------------------------------------------
; 12+log(O/H) [KK04] vs 12+log(O/H) [Empirical] - SDSS
; ------------------------------------------------------------

    psname = 'sdss_12oh_kk04_vs_12oh_empirical'
    im_openclose, pspath+psname, postscript=postscript, encapsulated=encapsulated, xsize=8.5, ysize=10.5

    pagemaker, nx=1, ny=3, height=3.0*[1,1,1], width=6.8, xmargin=[1.3,0.4], $
      ymargin=[0.4,1.1], xspace=0, yspace=0, xpage=8.5, ypage=10.5, $
      position=pos, /normal

    xtitle = '12+log(O/H) [KK04]'
    ytitle = '12+log(O/H)'

    xrange = ohrange11
    yrange = xrange

; --------------------------------------------------    
; [N II]/Ha
; --------------------------------------------------    

    indx = where((sdssnodust.zstrong_12oh_niiha_pettini gt -900) and $
      (sdssnodust.zstrong_12oh_kk04 gt -900),nindx)

    x = sdssnodust[indx].zstrong_12oh_kk04
    xerr = sdssnodust[indx].zstrong_12oh_kk04_err

    y = sdssnodust[indx].zstrong_12oh_niiha_pettini
    yerr = sdssnodust[indx].zstrong_12oh_niiha_pettini_err

    residuals = y-x
    stats = im_stats(residuals)
    xstr = strtrim(string(stats.median_rej,format='(F12.2)'),2)+' ('+$
      strtrim(string(stats.mean_rej,format='(F12.2)'),2)+$
      '\pm'+strtrim(string(stats.sigma_rej,format='(F12.2)'),2)+')'

    sdss_lineplot, x, y, xerr, yerr, plottype=1, postscript=postscript, $
      xtitle='', ytitle=ytitle, xrange=xrange, yrange=yrange, legendtype=0, $
      /right, /top, position=pos[*,0], charsize=charsize_8, xtickname=replicate(' ',10)
    djs_oplot, !x.crange, !y.crange, thick=postthick, color=djs_icolor(talkcolor), line=2

    legend, textoidl(xstr), /right, /bottom, box=0, charsize=charsize_5, $
      charthick=postthick, textcolor=djs_icolor(talkcolor)
    legend, 'N2', /left, /top, box=0, $
;   legend, textoidl('[N II]/H\alpha'), /left, /top, box=0, $
      charsize=charsize_8, charthick=postthick, textcolor=djs_icolor(talkcolor)
;   legend, textoidl('(b) [N II]/H\alpha'), /left, /top, box=0, $
;     charsize=charsize_8, charthick=postthick, textcolor=djs_icolor(talkcolor)

; --------------------------------------------------    
; ([O III]/Hb)/([N II]/Ha)
; --------------------------------------------------    

    indx = where((sdssnodust.zstrong_12oh_oiiinii_pettini gt -900) and $
      (sdssnodust.zstrong_12oh_kk04 gt -900),nindx)

    x = sdssnodust[indx].zstrong_12oh_kk04
    xerr = sdssnodust[indx].zstrong_12oh_kk04_err

    y = sdssnodust[indx].zstrong_12oh_oiiinii_pettini
    yerr = sdssnodust[indx].zstrong_12oh_oiiinii_pettini_err

    residuals = y-x
    stats = im_stats(residuals)
    xstr = strtrim(string(stats.median_rej,format='(F12.2)'),2)+' ('+$
      strtrim(string(stats.mean_rej,format='(F12.2)'),2)+$
      '\pm'+strtrim(string(stats.sigma_rej,format='(F12.2)'),2)+')'

    sdss_lineplot, x, y, xerr, yerr, plottype=1, postscript=postscript, $
      xtitle='', ytitle=ytitle, xrange=xrange, yrange=yrange, legendtype=0, /noerase, $
      /right, /top, position=pos[*,1], charsize=charsize_8, xtickname=replicate(' ',10)
    djs_oplot, !x.crange, !y.crange, thick=postthick, color=djs_icolor(talkcolor), line=2

    legend, textoidl(xstr), /right, /bottom, box=0, charsize=charsize_5, $
      charthick=postthick, textcolor=djs_icolor(talkcolor)
    legend, 'O3N2', /left, /top, box=0, $
;   legend, textoidl('([O III]/H\beta)/([N II]/H\alpha)'), /left, /top, box=0, $
      charsize=charsize_8, charthick=postthick, textcolor=djs_icolor(talkcolor)
;   legend, textoidl('(a) ([O III]/H\beta)/([N II]/H\alpha)'), /left, /top, box=0, $
;     charsize=charsize_8, charthick=postthick, textcolor=djs_icolor(talkcolor)

; --------------------------------------------------    
; P-method
; --------------------------------------------------    

    indx = where((sdssnodust.zstrong_12oh_kk04 gt -900) and (sdssnodust.zstrong_12oh_pt05 gt -900.0),nindx)

    x = sdssnodust[indx].zstrong_12oh_kk04
    xerr = sdssnodust[indx].zstrong_12oh_kk04_err

    y = sdssnodust[indx].zstrong_12oh_pt05
    yerr = sdssnodust[indx].zstrong_12oh_pt05_err

    residuals = y-x
    stats = im_stats(residuals)
    xstr = strtrim(string(stats.median_rej,format='(F12.2)'),2)+' ('+$
      strtrim(string(stats.mean_rej,format='(F12.2)'),2)+$
      '\pm'+strtrim(string(stats.sigma_rej,format='(F12.2)'),2)+')'

    sdss_lineplot, x, y, xerr, yerr, plottype=1, postscript=postscript, $
      xtitle=xtitle, ytitle=ytitle, xrange=xrange, yrange=yrange, legendtype=0, $
      /right, /top, position=pos[*,2], charsize=charsize_8, /noerase
    djs_oplot, !x.crange, !y.crange, thick=postthick, color=djs_icolor(talkcolor), line=2

    legend, textoidl(xstr), /right, /bottom, box=0, charsize=charsize_5, $
      charthick=postthick, textcolor=djs_icolor(talkcolor)
    legend, textoidl('PT05'), /left, /top, box=0, $
      charsize=charsize_8, charthick=postthick, textcolor=djs_icolor(talkcolor)
;   legend, textoidl('(c) P-method'), /left, /top, box=0, $
;     charsize=charsize_8, charthick=postthick, textcolor=djs_icolor(talkcolor)

    im_openclose, postscript=postscript, /close    

; ------------------------------------------------------------
; 12+log(O/H) [M91] vs 12+log(O/H) [Empirical] - SDSS
; ------------------------------------------------------------

    psname = 'sdss_12oh_m91_vs_12oh_empirical'
    im_openclose, pspath+psname, postscript=postscript, encapsulated=encapsulated, xsize=8.5, ysize=10.5

    pagemaker, nx=1, ny=3, height=3.0*[1,1,1], width=6.8, xmargin=[1.3,0.4], $
      ymargin=[0.4,1.1], xspace=0, yspace=0, xpage=8.5, ypage=10.5, $
      position=pos, /normal

;   if keyword_set(dotalk) then polyfill, [1,1,0,0,1], [1,0,0,1,1], /normal, $
;     color=djs_icolor('black')

    xtitle = '12+log(O/H) [M91]'
    ytitle = '12+log(O/H)'

    xrange = ohrange11
    yrange = xrange

; --------------------------------------------------    
; [N II]/Ha
; --------------------------------------------------    

    indx = where((sdssnodust.zstrong_12oh_niiha_pettini gt -900) and $
      (sdssnodust.zstrong_12oh_m91 gt -900),nindx)

    x = sdssnodust[indx].zstrong_12oh_m91
    xerr = sdssnodust[indx].zstrong_12oh_m91_err

    y = sdssnodust[indx].zstrong_12oh_niiha_pettini
    yerr = sdssnodust[indx].zstrong_12oh_niiha_pettini_err

    residuals = y-x
    stats = im_stats(residuals)
    xstr = strtrim(string(stats.median_rej,format='(F12.2)'),2)+' ('+$
      strtrim(string(stats.mean_rej,format='(F12.2)'),2)+$
      '\pm'+strtrim(string(stats.sigma_rej,format='(F12.2)'),2)+')'

    sdss_lineplot, x, y, xerr, yerr, plottype=1, postscript=postscript, $
      xtitle='', ytitle=ytitle, xrange=xrange, yrange=yrange, legendtype=0, $
      /right, /top, position=pos[*,0], charsize=charsize_8, xtickname=replicate(' ',10)
    djs_oplot, !x.crange, !y.crange, thick=postthick, color=djs_icolor(talkcolor), line=2

    legend, textoidl(xstr), /right, /bottom, box=0, charsize=charsize_5, $
      charthick=postthick, textcolor=djs_icolor(talkcolor)
    legend, 'N2', /left, /top, box=0, $
;   legend, textoidl('[N II]/H\alpha'), /left, /top, box=0, $
      charsize=charsize_8, charthick=postthick, textcolor=djs_icolor(talkcolor)
;   legend, textoidl('(b) [N II]/H\alpha'), /left, /top, box=0, $
;     charsize=charsize_8, charthick=postthick, textcolor=djs_icolor(talkcolor)

; --------------------------------------------------    
; ([O III]/Hb)/([N II]/Ha)
; --------------------------------------------------    

    indx = where((sdssnodust.zstrong_12oh_oiiinii_pettini gt -900) and $
      (sdssnodust.zstrong_12oh_m91 gt -900),nindx)

    x = sdssnodust[indx].zstrong_12oh_m91
    xerr = sdssnodust[indx].zstrong_12oh_m91_err

    y = sdssnodust[indx].zstrong_12oh_oiiinii_pettini
    yerr = sdssnodust[indx].zstrong_12oh_oiiinii_pettini_err

    residuals = y-x
    stats = im_stats(residuals)
    xstr = strtrim(string(stats.median_rej,format='(F12.2)'),2)+' ('+$
      strtrim(string(stats.mean_rej,format='(F12.2)'),2)+$
      '\pm'+strtrim(string(stats.sigma_rej,format='(F12.2)'),2)+')'

    sdss_lineplot, x, y, xerr, yerr, plottype=1, postscript=postscript, $
      xtitle='', ytitle=ytitle, xrange=xrange, yrange=yrange, legendtype=0, /noerase, $
      /right, /top, position=pos[*,1], charsize=charsize_8, xtickname=replicate(' ',10)
    djs_oplot, !x.crange, !y.crange, thick=postthick, color=djs_icolor(talkcolor), line=2

    legend, textoidl(xstr), /right, /bottom, box=0, charsize=charsize_5, $
      charthick=postthick, textcolor=djs_icolor(talkcolor)
    legend, 'O3N2', /left, /top, box=0, $
;   legend, textoidl('([O III]/H\beta)/([N II]/H\alpha)'), /left, /top, box=0, $
      charsize=charsize_8, charthick=postthick, textcolor=djs_icolor(talkcolor)
;   legend, textoidl('(a) ([O III]/H\beta)/([N II]/H\alpha)'), /left, /top, box=0, $
;     charsize=charsize_8, charthick=postthick, textcolor=djs_icolor(talkcolor)

; --------------------------------------------------    
; P-method
; --------------------------------------------------    

    indx = where((sdssnodust.zstrong_12oh_m91 gt -900) and (sdssnodust.zstrong_12oh_pt05 gt -900.0),nindx)

    x = sdssnodust[indx].zstrong_12oh_m91
    xerr = sdssnodust[indx].zstrong_12oh_m91_err

    y = sdssnodust[indx].zstrong_12oh_pt05
    yerr = sdssnodust[indx].zstrong_12oh_pt05_err

    residuals = y-x
    stats = im_stats(residuals)
    xstr = strtrim(string(stats.median_rej,format='(F12.2)'),2)+' ('+$
      strtrim(string(stats.mean_rej,format='(F12.2)'),2)+$
      '\pm'+strtrim(string(stats.sigma_rej,format='(F12.2)'),2)+')'

    sdss_lineplot, x, y, xerr, yerr, plottype=1, postscript=postscript, $
      xtitle=xtitle, ytitle=ytitle, xrange=xrange, yrange=yrange, legendtype=0, $
      /right, /top, position=pos[*,2], charsize=charsize_8, /noerase
    djs_oplot, !x.crange, !y.crange, thick=postthick, color=djs_icolor(talkcolor), line=2

    legend, textoidl(xstr), /right, /bottom, box=0, charsize=charsize_5, $
      charthick=postthick, textcolor=djs_icolor(talkcolor)
    legend, textoidl('PT05'), /left, /top, box=0, $
      charsize=charsize_8, charthick=postthick, textcolor=djs_icolor(talkcolor)
;   legend, textoidl('(c) P-method'), /left, /top, box=0, $
;     charsize=charsize_8, charthick=postthick, textcolor=djs_icolor(talkcolor)

    im_openclose, postscript=postscript, /close    

; ------------------------------------------------------------
; 12+log(O/H) [N2] versus 12+log(O/H) [T04]
; ------------------------------------------------------------

    psname = 'sdss_12oh_n2_vs_12oh_t04'
    im_openclose, pspath+psname, postscript=postscript, xsize=8.5, ysize=9.0, encapsulated=encapsulated

    pagemaker, nx=1, ny=2, xspace=0, yspace=0.0, width=6.5, height=[5.0,2.5], $
      xmargin=[1.5,0.5], ymargin=[0.4,1.1], xpage=8.5, ypage=9.0, $
      position=pos, /normal

    indx = where((sdssdust.zstrong_12oh_niiha_pettini gt -900) and (sdssnodust.zstrong_12oh_t04 gt -900.0),nindx)
;   indx = where((sdssdust.zstrong_12oh_niiha_pettini gt -900) and (sdssancillary.tremonti_oh gt -900.0),nindx)

    x = sdssdust[indx].zstrong_12oh_niiha_pettini
    xerr = sdssdust[indx].zstrong_12oh_niiha_pettini_err

    y = sdssnodust[indx].zstrong_12oh_t04
    yerr = sdssnodust[indx].zstrong_12oh_t04_err
;   y = sdssancillary[indx].tremonti_oh
;   yerr = sdssancillary[indx].tremonti_oh_err

    yresid = y-x
;   w = where(x lt 8.55)
    stats = im_stats(yresid,/verbose)
    xstr = strtrim(string(stats.median_rej,format='(F12.2)'),2)+' ('+$
      strtrim(string(stats.mean_rej,format='(F12.2)'),2)+$
      '\pm'+strtrim(string(stats.sigma_rej,format='(F12.2)'),2)+')'
    
    xtitle = '12 + log (O/H) [N2]'
    ytitle = '12 + log (O/H) [T04]'
    yresidtitle = '\Delta'+'log(O/H)'

    xrange = ohrange11
    yrange = xrange
    yresidrange = 0.9*[-1,1]

    sdss_lineplot, x, y, xerr, yerr, plottype=1, postscript=postscript, $
      xtitle='', ytitle=ytitle, xrange=xrange, yrange=yrange, legendtype=0, $
      /right, /top, position=pos[*,0], xtickname=replicate(' ',10), charsize=charsize_8
    djs_oplot, !x.crange, !y.crange, line=0, thick=postthick
    djs_oplot, 8.55*[1,1], !y.crange, line=2, thick=postthick

    sdss_lineplot, x, yresid, xerr, yerr, plottype=1, postscript=postscript, $
      xtitle=xtitle, ytitle=yresidtitle, xrange=xrange, yrange=yresidrange, $
      legendtype=0, /right, /top, /noerase, position=pos[*,1], charsize=charsize_8
    djs_oplot, !x.crange, [0,0], line=0, thick=postthick
    djs_oplot, 8.55*[1,1], !y.crange, line=2, thick=postthick

    legend, textoidl(xstr), /right, /top, box=0, charsize=charsize_5, charthick=postthick

    im_openclose, postscript=postscript, /close    

; ------------------------------------------------------------
; 12+log(O/H) [O3N2] vs 12+log(O/H) [T04]
; ------------------------------------------------------------

    psname = 'sdss_12oh_o3n2_vs_12oh_t04'
    im_openclose, pspath+psname, postscript=postscript, xsize=8.5, ysize=9.0, encapsulated=encapsulated

    pagemaker, nx=1, ny=2, xspace=0, yspace=0.0, width=6.5, height=[5.0,2.5], $
      xmargin=[1.5,0.5], ymargin=[0.4,1.1], xpage=8.5, ypage=9.0, $
      position=pos, /normal

    indx = where((sdssdust.zstrong_12oh_oiiinii_pettini gt -900) and (sdssnodust.zstrong_12oh_t04 gt -900.0),nindx)
;   indx = where((sdssdust.zstrong_12oh_oiiinii_pettini gt -900) and (sdssancillary.tremonti_oh gt -900.0),nindx)

    x = sdssdust[indx].zstrong_12oh_oiiinii_pettini
    xerr = sdssdust[indx].zstrong_12oh_oiiinii_pettini_err

    y = sdssnodust[indx].zstrong_12oh_t04
    yerr = sdssnodust[indx].zstrong_12oh_t04_err
;   y = sdssancillary[indx].tremonti_oh
;   yerr = sdssancillary[indx].tremonti_oh_err

    yresid = y-x
    stats = im_stats(yresid,/verbose)
    xstr = strtrim(string(stats.median_rej,format='(F12.2)'),2)+' ('+$
      strtrim(string(stats.mean_rej,format='(F12.2)'),2)+$
      '\pm'+strtrim(string(stats.sigma_rej,format='(F12.2)'),2)+')'
    
    xtitle = '12 + log (O/H) [O3N2]'
    ytitle = '12 + log (O/H) [T04]'
    yresidtitle = '\Delta'+'log(O/H)'

    xrange = ohrange11
    yrange = xrange
    yresidrange = 0.9*[-1,1]

    sdss_lineplot, x, y, xerr, yerr, plottype=1, postscript=postscript, $
      xtitle='', ytitle=ytitle, xrange=xrange, yrange=yrange, legendtype=0, $
      /right, /top, position=pos[*,0], xtickname=replicate(' ',10), charsize=charsize_8
    djs_oplot, !x.crange, !y.crange, line=0, thick=postthick

    sdss_lineplot, x, yresid, xerr, yerr, plottype=1, postscript=postscript, $
      xtitle=xtitle, ytitle=yresidtitle, xrange=xrange, yrange=yresidrange, $
      legendtype=0, /right, /top, /noerase, position=pos[*,1], charsize=charsize_8
    djs_oplot, !x.crange, [0,0], line=0, thick=postthick

    legend, textoidl(xstr), /right, /top, box=0, charsize=charsize_5, charthick=postthick

    im_openclose, postscript=postscript, /close    

; ------------------------------------------------------------
; 12+log(O/H) [O3N2] versus 12+log(O/H) [N2]
; ------------------------------------------------------------

    psname = 'sdss_12oh_o3n2_vs_12oh_n2'
    im_openclose, pspath+psname, postscript=postscript, xsize=8.5, ysize=9.0, encapsulated=encapsulated

    pagemaker, nx=1, ny=2, xspace=0, yspace=0.0, width=6.5, height=[5.0,2.5], $
      xmargin=[1.5,0.5], ymargin=[0.4,1.1], xpage=8.5, ypage=9.0, $
      position=pos, /normal

    indx = where((sdssdust.zstrong_12oh_oiiinii_pettini gt -900) and (sdssdust.zstrong_12oh_niiha_pettini gt -900.0),nindx)

    x = sdssdust[indx].zstrong_12oh_oiiinii_pettini
    xerr = sdssdust[indx].zstrong_12oh_oiiinii_pettini_err

    y = sdssdust[indx].zstrong_12oh_niiha_pettini
    yerr = sdssdust[indx].zstrong_12oh_niiha_pettini_err

    yresid = y-x
;   w = where(x lt 8.55)
    stats = im_stats(yresid,/verbose)
    xstr = strtrim(string(stats.median_rej,format='(F12.2)'),2)+' ('+$
      strtrim(string(stats.mean_rej,format='(F12.2)'),2)+$
      '\pm'+strtrim(string(stats.sigma_rej,format='(F12.2)'),2)+')'
    
    xtitle = '12 + log (O/H) [O3N2]'
    ytitle = '12 + log (O/H) [N2]'
    yresidtitle = '\Delta'+'log(O/H)'

    xrange = ohrange11
    yrange = xrange
    yresidrange = 0.9*[-1,1]

    sdss_lineplot, x, y, xerr, yerr, plottype=1, postscript=postscript, $
      xtitle='', ytitle=ytitle, xrange=xrange, yrange=yrange, legendtype=0, $
      /right, /top, position=pos[*,0], xtickname=replicate(' ',10), charsize=charsize_8
    djs_oplot, !x.crange, !y.crange, line=0, thick=postthick
    djs_oplot, 8.55*[1,1], !y.crange, line=2, thick=postthick

    sdss_lineplot, x, yresid, xerr, yerr, plottype=1, postscript=postscript, $
      xtitle=xtitle, ytitle=yresidtitle, xrange=xrange, yrange=yresidrange, $
      legendtype=0, /right, /top, /noerase, position=pos[*,1], charsize=charsize_8
    djs_oplot, !x.crange, [0,0], line=0, thick=postthick
    djs_oplot, 8.55*[1,1], !y.crange, line=2, thick=postthick

    legend, textoidl(xstr), /right, /top, box=0, charsize=charsize_5, charthick=postthick

    im_openclose, postscript=postscript, /close    
    
; --------------------------------------------------    
; GENERATE THE WEB PAGE
; --------------------------------------------------    
    
    if keyword_set(postscript) and keyword_set(encapsulated) then $
      im_ps2html, htmlbase, html_path=html_path, cleanpng=0, npscols=3, _extra=extra

stop
    
; # WORK IN PROGRESS

; ------------------------------------------------------------
; 3-panel 12+log(O/H) KK04 vs 12+log(O/H) P-method [HII,Integrated,SDSS]
; ------------------------------------------------------------

    psname = '12oh_kk04_vs_12oh_pt05_3panel'
    im_openclose, pspath+psname, postscript=postscript, encapsulated=encapsulated, $
      xsize=8.5, ysize=3.9;, /landscape

    pagemaker, nx=3, ny=1, height=2.5, width=2.5*[1,1,1], xmargin=[0.8,0.2], $
      ymargin=[0.4,1.0], xspace=0, yspace=0, xpage=8.5, ypage=3.9, $
      position=pos, /normal

;   if keyword_set(dotalk) then polyfill, [1,1,0,0,1], [1,0,0,1,1], /normal, $
;     color=djs_icolor('black')

    xtitle = '12+log(O/H) [KK04]'
    ytitle = '12+log(O/H) [PT05]'

    xrange = ohrange4
    yrange = ohrange4

; --------------------------------------------------    
; HII Regions
; --------------------------------------------------    

    indx = where((hii.zstrong_12oh_kk04 gt -900) and $
      (hii.zstrong_12oh_kk04 gt -900) and $
      (hii.zstrong_12oh_pt05 gt -900) and $
      (hii.zstrong_12oh_pt05 gt -900) and $
      (hii.zstrong_niiha gt -900),nindx)

    up = where((hii[indx].zstrong_12oh_pt05 gt 8.2) and (hii[indx].zstrong_niiha gt -1.0),nup)
;   up = where((hii[indx].zstrong_12oh_niiha_pettini gt 8.2),nup)
    if (nup ne 0L) then begin
       xup = hii[indx[up]].zstrong_12oh_kk04
       xerrup = hii[indx[up]].zstrong_12oh_kk04_err

       yup = hii[indx[up]].zstrong_12oh_pt05
       yerrup = hii[indx[up]].zstrong_12oh_pt05_err

       x = xup & xerr = xerrup
       y = yup & yerr = yerrup
    endif

    lo = where((hii[indx].zstrong_12oh_pt05 lt 7.95) and (hii[indx].zstrong_niiha lt -1.0),nlo)
;   lo = where((hii[indx].zstrong_12oh_niiha_pettini lt 7.95),nlo)
    if (nlo ne 0L) then begin
       xlo = hii[indx[lo]].zstrong_12oh_kk04
       xerrlo = hii[indx[lo]].zstrong_12oh_kk04_err

       ylo = hii[indx[lo]].zstrong_12oh_pt05
       yerrlo = hii[indx[lo]].zstrong_12oh_pt05_err

       x = [xup,xlo] & xerr = [xerrup,xerrlo]
       y = [yup,ylo] & yerr = [yerrup,yerrlo]
    endif

    upstats = im_stats(yup-xup)
    upxstr = strtrim(string(upstats.median_rej,format='(F12.2)'),2)+' ('+$
      strtrim(string(upstats.mean_rej,format='(F12.2)'),2)+$
      '\pm'+strtrim(string(upstats.sigma_rej,format='(F12.2)'),2)+')'

    lostats = im_stats(ylo-xlo)
    loxstr = strtrim(string(lostats.median_rej,format='(F12.2)'),2)+' ('+$
      strtrim(string(lostats.mean_rej,format='(F12.2)'),2)+$
      '\pm'+strtrim(string(lostats.sigma_rej,format='(F12.2)'),2)+')'

    atlas1d_lineplot, x, y, xerr, yerr, plottype=1, postscript=postscript, $
      xtitle='', ytitle=ytitle, xrange=xrange, yrange=yrange, legendtype=0, $
      /right, /top, position=pos[*,0], charsize=charsize_2, hiicolor='purple'
    djs_oplot, !x.crange, !y.crange, thick=postthick, color=djs_icolor(talkcolor), line=2

    legend, textoidl([upxstr+' (Upper)',loxstr+' (Lower)']), /left, /top, box=0, $
      charsize=charsize_1, charthick=postthick, textcolor=djs_icolor(talkcolor), margin=0
;   legend, textoidl('(a)'), /left, /top, box=0, $
;     charsize=charsize_2, charthick=postthick, textcolor=djs_icolor(talkcolor)

; --------------------------------------------------    
; Integrated
; --------------------------------------------------    

; Atlas    
    
    indx = where((intnodust.zstrong_12oh_kk04 gt -900) and $
      (intnodust.zstrong_12oh_kk04 gt -900) and $
      (intnodust.zstrong_12oh_pt05 gt -900) and $
      (intnodust.zstrong_12oh_pt05 gt -900) and $
      (intnodust.zstrong_niiha gt -900),nindx)

    up = where((intnodust[indx].zstrong_12oh_pt05 gt 8.2) and $
      (intnodust[indx].zstrong_niiha gt -1.0),nup)
;   up = where((intnodust[indx].zstrong_12oh_niiha_pettini gt 8.2),nup)
    if (nup ne 0L) then begin
       xup = intnodust[indx[up]].zstrong_12oh_kk04
       xerrup = intnodust[indx[up]].zstrong_12oh_kk04_err

       yup = intnodust[indx[up]].zstrong_12oh_pt05
       yerrup = intnodust[indx[up]].zstrong_12oh_pt05_err

       x = xup & xerr = xerrup
       y = yup & yerr = yerrup
    endif

    lo = where((intnodust[indx].zstrong_12oh_pt05 lt 7.95) and $
      (intnodust[indx].zstrong_niiha lt -1.0),nlo)
;   lo = where((intnodust[indx].zstrong_12oh_niiha_pettini lt 7.95),nlo)
    if (nlo ne 0L) then begin
       xlo = intnodust[indx[lo]].zstrong_12oh_kk04
       xerrlo = intnodust[indx[lo]].zstrong_12oh_kk04_err

       ylo = intnodust[indx[lo]].zstrong_12oh_pt05
       yerrlo = intnodust[indx[lo]].zstrong_12oh_pt05_err

       x = [xup,xlo] & xerr = [xerrup,xerrlo]
       y = [yup,ylo] & yerr = [yerrup,yerrlo]
    endif

    atlas1d_lineplot, x, y, xerr, yerr, plottype=1, postscript=postscript, $
      xtitle=xtitle, ytitle='', xrange=xrange, yrange=yrange, legendtype=0, /noerase, $
      /right, /top, position=pos[*,1], charsize=charsize_2, ytickname=replicate(' ',10)
    djs_oplot, !x.crange, !y.crange, thick=postthick, color=djs_icolor(talkcolor), line=2

    legend, textoidl([upxstr+' (Upper)',loxstr+' (Lower)']), /left, /top, box=0, $
      charsize=charsize_1, charthick=postthick, textcolor=djs_icolor(talkcolor), margin=0
;   legend, textoidl('(b)'), /left, /top, box=0, $
;     charsize=charsize_2, charthick=postthick, textcolor=djs_icolor(talkcolor)

; --------------------------------------------------    
; SDSS
; --------------------------------------------------    

    indx = where((sdssnodust.zstrong_12oh_kk04 gt -900) and $
      (sdssnodust.zstrong_12oh_kk04 gt -900) and $
      (sdssnodust.zstrong_12oh_pt05 gt -900) and $
      (sdssnodust.zstrong_12oh_pt05 gt -900) and $
      (sdssnodust.zstrong_niiha gt -900),nindx)

    up = where((sdssnodust[indx].zstrong_12oh_pt05 gt 8.2) and $
      (sdssnodust[indx].zstrong_niiha gt -1.0),nup)
;   up = where((sdssnodust[indx].zstrong_12oh_niiha_pettini gt 8.2),nup)
    if (nup ne 0L) then begin
       xup = sdssnodust[indx[up]].zstrong_12oh_kk04
       xerrup = sdssnodust[indx[up]].zstrong_12oh_kk04_err

       yup = sdssnodust[indx[up]].zstrong_12oh_pt05
       yerrup = sdssnodust[indx[up]].zstrong_12oh_pt05_err

       x = xup & xerr = xerrup
       y = yup & yerr = yerrup
    endif

    lo = where((sdssnodust[indx].zstrong_12oh_pt05 lt 7.95) and $
      (sdssnodust[indx].zstrong_niiha lt -1.0),nlo)
;   lo = where((sdssnodust[indx].zstrong_12oh_niiha_pettini lt 7.95),nlo)
    if (nlo ne 0L) then begin
       xlo = sdssnodust[indx[lo]].zstrong_12oh_kk04
       xerrlo = sdssnodust[indx[lo]].zstrong_12oh_kk04_err

       ylo = sdssnodust[indx[lo]].zstrong_12oh_pt05
       yerrlo = sdssnodust[indx[lo]].zstrong_12oh_pt05_err

       x = [xup,xlo] & xerr = [xerrup,xerrlo]
       y = [yup,ylo] & yerr = [yerrup,yerrlo]
    endif

    upstats = im_stats(yup-xup)
    upxstr = strtrim(string(upstats.median_rej,format='(F12.2)'),2)+' ('+$
      strtrim(string(upstats.mean_rej,format='(F12.2)'),2)+$
      '\pm'+strtrim(string(upstats.sigma_rej,format='(F12.2)'),2)+')'

    lostats = im_stats(ylo-xlo)
    loxstr = strtrim(string(lostats.median_rej,format='(F12.2)'),2)+' ('+$
      strtrim(string(lostats.mean_rej,format='(F12.2)'),2)+$
      '\pm'+strtrim(string(lostats.sigma_rej,format='(F12.2)'),2)+')'

    sdss_lineplot, x, y, xerr, yerr, plottype=1, postscript=postscript, $
      xtitle='', ytitle='', xrange=xrange, yrange=yrange, legendtype=0, $
      /right, /top, position=pos[*,2], charsize=charsize_2, ytickname=replicate(' ',10), $
      /noerase
    djs_oplot, !x.crange, !y.crange, thick=postthick, color=djs_icolor(talkcolor), line=2

    legend, textoidl(upxstr+' (Upper)'), /left, /top, box=0, $
      charsize=charsize_1, charthick=postthick, textcolor=djs_icolor(talkcolor), margin=0
;   legend, textoidl([upxstr+' (Upper)',loxstr+' (Lower)']), /right, /bottom, box=0, $
;     charsize=charsize_1, charthick=postthick, textcolor=djs_icolor(talkcolor)
;   legend, textoidl('(c)'), /left, /top, box=0, $
;     charsize=charsize_2, charthick=postthick, textcolor=djs_icolor(talkcolor)

    im_openclose, postscript=postscript, /close    

; ------------------------------------------------------------
; 3-panel BPT diagrams [HII,Integrated,SDSS]
; ------------------------------------------------------------

    psname = 'niiha_vs_oiiihb_3panel'
    im_openclose, pspath+psname, postscript=postscript, encapsulated=encapsulated, xsize=8.5, ysize=3.9

    pagemaker, nx=3, ny=1, height=2.5, width=2.5*[1,1,1], xmargin=[0.8,0.2], $
      ymargin=[0.4,1.0], xspace=0, yspace=0, xpage=8.5, ypage=3.9, $
      position=pos, /normal

;   if keyword_set(dotalk) then polyfill, [1,1,0,0,1], [1,0,0,1,1], /normal, $
;     color=djs_icolor('black')

    xtitle = 'log ([N II]/H\alpha)'
    ytitle = 'log ([O III]/H\beta)'

    xrange = niiharange
    yrange = oiiihbrange

; --------------------------------------------------    
; HII Regions
; --------------------------------------------------    

    indx = where((hii.zstrong_niiha gt -900) and (hii.zstrong_oiiihbeta gt -900),nindx)
    
    x = hii[indx].zstrong_niiha
    xerr = hii[indx].zstrong_niiha_err
 
    y = hii[indx].zstrong_oiiihbeta
    yerr = hii[indx].zstrong_oiiihbeta_err
 
    atlas1d_lineplot, x, y, xerr, yerr, plottype=1, postscript=postscript, $
      xtitle=xtitle, ytitle=ytitle, xrange=xrange, yrange=yrange, legendtype=0, $
      /right, /top, position=pos[*,0], charsize=charsize_1, hiicolor='purple'
;   legend, textoidl('(a)'), /right, /top, box=0, $
;     charsize=charsize_1, charthick=postthick, textcolor=djs_icolor(talkcolor)

; overplot L. Kewley's starburst mixing line

    oplot, models.x_nii, models.y_nii, line=0, thick=2.0

; --------------------------------------------------    
; Integrated
; --------------------------------------------------    

; Atlas    
    
    indx = where((intdust.zstrong_niiha gt -900) and (intdust.zstrong_oiiihbeta gt -900),nindx)
    
    x = intdust[indx].zstrong_niiha
    xerr = intdust[indx].zstrong_niiha_err
 
    y = intdust[indx].zstrong_oiiihbeta
    yerr = intdust[indx].zstrong_oiiihbeta_err
 
    atlas1d_lineplot, x, y, xerr, yerr, plottype=1, postscript=postscript, $
      xtitle=xtitle, ytitle='', xrange=xrange, yrange=yrange, legendtype=0, /noerase, $
      /right, /top, position=pos[*,1], charsize=charsize_1, ytickname=replicate(' ',10)
;   legend, textoidl('(b)'), /right, /top, box=0, $
;     charsize=charsize_1, charthick=postthick, textcolor=djs_icolor(talkcolor)

; overplot L. Kewley's starburst mixing line

    oplot, models.x_nii, models.y_nii, line=0, thick=2.0

; --------------------------------------------------    
; SDSS
; --------------------------------------------------    

    indx = where((sdssdust.zstrong_niiha gt -900) and (sdssdust.zstrong_oiiihbeta gt -900),nindx)
    
    x = sdssdust[indx].zstrong_niiha
    xerr = sdssdust[indx].zstrong_niiha_err
 
    y = sdssdust[indx].zstrong_oiiihbeta
    yerr = sdssdust[indx].zstrong_oiiihbeta_err
 
    sdss_lineplot, x, y, xerr, yerr, plottype=1, postscript=postscript, $
      xtitle=xtitle, ytitle='', xrange=xrange, yrange=yrange, legendtype=0, $
      /right, /top, position=pos[*,2], charsize=charsize_1, ytickname=replicate(' ',10), $
      /noerase
;   legend, textoidl('(c)'), /right, /top, box=0, $
;     charsize=charsize_1, charthick=postthick, textcolor=djs_icolor(talkcolor)

; overplot L. Kewley's starburst mixing line

    oplot, models.x_nii, models.y_nii, line=0, thick=2.0

    im_openclose, postscript=postscript, /close    

; ------------------------------------------------------------
; 3-panel 12+log(O/H) M91 vs 12+log(O/H) P-method [HII,Integrated,SDSS]
; ------------------------------------------------------------

    psname = '12oh_m91_vs_12oh_pt05_3panel'
    im_openclose, pspath+psname, postscript=postscript, encapsulated=encapsulated, xsize=8.5, ysize=3.9

    pagemaker, nx=3, ny=1, height=2.5, width=2.5*[1,1,1], xmargin=[0.8,0.2], $
      ymargin=[0.4,1.0], xspace=0, yspace=0, xpage=8.5, ypage=3.9, $
      position=pos, /normal

;   if keyword_set(dotalk) then polyfill, [1,1,0,0,1], [1,0,0,1,1], /normal, $
;     color=djs_icolor('black')

    xtitle = '12+log(O/H) [M91]'
    ytitle = '12+log(O/H) [PT05]'

    xrange = ohrange4
    yrange = ohrange4

; --------------------------------------------------    
; HII Regions
; --------------------------------------------------    

    indx = where((hii.zstrong_12oh_m91 gt -900) and (hii.zstrong_12oh_pt05 gt -900) and $
      (hii.zstrong_12oh_pt05 gt -900) and (hii.zstrong_12oh_niiha_pettini gt -900),nindx)

    up = where((hii[indx].zstrong_12oh_pt05 gt 8.2) and (hii[indx].zstrong_niiha gt -1.0),nup)
;   up = where((hii[indx].zstrong_12oh_niiha_pettini gt 8.2),nup)
    if (nup ne 0L) then begin
       xup = hii[indx[up]].zstrong_12oh_m91
       xerrup = hii[indx[up]].zstrong_12oh_m91_err

       yup = hii[indx[up]].zstrong_12oh_pt05
       yerrup = hii[indx[up]].zstrong_12oh_pt05_err

       x = xup & xerr = xerrup
       y = yup & yerr = yerrup
    endif

    lo = where((hii[indx].zstrong_12oh_pt05 lt 7.95) and (hii[indx].zstrong_niiha lt -1.0),nlo)
;   lo = where((hii[indx].zstrong_12oh_niiha_pettini lt 7.95),nlo)
    if (nlo ne 0L) then begin
       xlo = hii[indx[lo]].zstrong_12oh_m91
       xerrlo = hii[indx[lo]].zstrong_12oh_m91_err

       ylo = hii[indx[lo]].zstrong_12oh_pt05
       yerrlo = hii[indx[lo]].zstrong_12oh_pt05_err

       x = [xup,xlo] & xerr = [xerrup,xerrlo]
       y = [yup,ylo] & yerr = [yerrup,yerrlo]
    endif

    upstats = im_stats(yup-xup)
    upxstr = strtrim(string(upstats.median_rej,format='(F12.2)'),2)+' ('+$
      strtrim(string(upstats.mean_rej,format='(F12.2)'),2)+$
      '\pm'+strtrim(string(upstats.sigma_rej,format='(F12.2)'),2)+')'

    lostats = im_stats(ylo-xlo)
    loxstr = strtrim(string(lostats.median_rej,format='(F12.2)'),2)+' ('+$
      strtrim(string(lostats.mean_rej,format='(F12.2)'),2)+$
      '\pm'+strtrim(string(lostats.sigma_rej,format='(F12.2)'),2)+')'

    atlas1d_lineplot, x, y, xerr, yerr, plottype=1, postscript=postscript, $
      xtitle='', ytitle=ytitle, xrange=xrange, yrange=yrange, legendtype=0, $
      /right, /top, position=pos[*,0], charsize=charsize_2, hiicolor='purple'
    djs_oplot, !x.crange, !y.crange, thick=postthick, color=djs_icolor(talkcolor), line=2

    legend, textoidl([upxstr+' (Upper)',loxstr+' (Lower)']), /top, /left, box=0, $
      charsize=charsize_1, charthick=postthick, textcolor=djs_icolor(talkcolor)
;   legend, textoidl('(a)'), /left, /top, box=0, $
;     charsize=charsize_2, charthick=postthick, textcolor=djs_icolor(talkcolor)

; --------------------------------------------------    
; Integrated
; --------------------------------------------------    

; Atlas    
    
    indx = where((intnodust.zstrong_12oh_m91 gt -900) and (intnodust.zstrong_12oh_pt05 gt -900) and $
      (intnodust.zstrong_12oh_pt05 gt -900) and (intnodust.zstrong_12oh_niiha_pettini gt -900),nindx)

    up = where((intnodust[indx].zstrong_12oh_pt05 gt 8.2) and $
      (intnodust[indx].zstrong_niiha gt -1.0),nup)
;   up = where((intnodust[indx].zstrong_12oh_niiha_pettini gt 8.2),nup)
    if (nup ne 0L) then begin
       xup = intnodust[indx[up]].zstrong_12oh_m91
       xerrup = intnodust[indx[up]].zstrong_12oh_m91_err

       yup = intnodust[indx[up]].zstrong_12oh_pt05
       yerrup = intnodust[indx[up]].zstrong_12oh_pt05_err

       x = xup & xerr = xerrup
       y = yup & yerr = yerrup
    endif

    lo = where((intnodust[indx].zstrong_12oh_pt05 lt 7.95) and $
      (intnodust[indx].zstrong_niiha lt -1.0),nlo)
;   lo = where((intnodust[indx].zstrong_12oh_niiha_pettini lt 7.95),nlo)
    if (nlo ne 0L) then begin
       xlo = intnodust[indx[lo]].zstrong_12oh_m91
       xerrlo = intnodust[indx[lo]].zstrong_12oh_m91_err

       ylo = intnodust[indx[lo]].zstrong_12oh_pt05
       yerrlo = intnodust[indx[lo]].zstrong_12oh_pt05_err

       x = [xup,xlo] & xerr = [xerrup,xerrlo]
       y = [yup,ylo] & yerr = [yerrup,yerrlo]
    endif
    
    atlas1d_lineplot, x, y, xerr, yerr, plottype=1, postscript=postscript, $
      xtitle=xtitle, ytitle='', xrange=xrange, yrange=yrange, legendtype=0, /noerase, $
      /right, /top, position=pos[*,1], charsize=charsize_2, ytickname=replicate(' ',10)
    djs_oplot, !x.crange, !y.crange, thick=postthick, color=djs_icolor(talkcolor), line=2

    legend, textoidl([upxstr+' (Upper)',loxstr+' (Lower)']), /left, /top, box=0, $
      charsize=charsize_1, charthick=postthick, textcolor=djs_icolor(talkcolor)
;   legend, textoidl('(b)'), /left, /top, box=0, $
;     charsize=charsize_2, charthick=postthick, textcolor=djs_icolor(talkcolor)

; --------------------------------------------------    
; SDSS
; --------------------------------------------------    

    indx = where((sdssnodust.zstrong_12oh_m91 gt -900) and (sdssnodust.zstrong_12oh_pt05 gt -900) and $
      (sdssnodust.zstrong_12oh_pt05 gt -900) and (sdssnodust.zstrong_12oh_niiha_pettini gt -900),nindx)

    up = where((sdssnodust[indx].zstrong_12oh_pt05 gt 8.2) and $
      (sdssnodust[indx].zstrong_niiha gt -1.0),nup)
;   up = where((sdssnodust[indx].zstrong_12oh_niiha_pettini gt 8.2),nup)
    if (nup ne 0L) then begin
       xup = sdssnodust[indx[up]].zstrong_12oh_m91
       xerrup = sdssnodust[indx[up]].zstrong_12oh_m91_err

       yup = sdssnodust[indx[up]].zstrong_12oh_pt05
       yerrup = sdssnodust[indx[up]].zstrong_12oh_pt05_err

       x = xup & xerr = xerrup
       y = yup & yerr = yerrup
    endif

    lo = where((sdssnodust[indx].zstrong_12oh_pt05 lt 7.95) and $
      (sdssnodust[indx].zstrong_niiha lt -1.0),nlo)
;   lo = where((sdssnodust[indx].zstrong_12oh_niiha_pettini lt 7.95),nlo)
    if (nlo ne 0L) then begin
       xlo = sdssnodust[indx[lo]].zstrong_12oh_m91
       xerrlo = sdssnodust[indx[lo]].zstrong_12oh_m91_err

       ylo = sdssnodust[indx[lo]].zstrong_12oh_pt05
       yerrlo = sdssnodust[indx[lo]].zstrong_12oh_pt05_err

       x = [xup,xlo] & xerr = [xerrup,xerrlo]
       y = [yup,ylo] & yerr = [yerrup,yerrlo]
    endif

    upstats = im_stats(yup-xup)
    upxstr = strtrim(string(upstats.median_rej,format='(F12.2)'),2)+' ('+$
      strtrim(string(upstats.mean_rej,format='(F12.2)'),2)+$
      '\pm'+strtrim(string(upstats.sigma_rej,format='(F12.2)'),2)+')'

    lostats = im_stats(ylo-xlo)
    loxstr = strtrim(string(lostats.median_rej,format='(F12.2)'),2)+' ('+$
      strtrim(string(lostats.mean_rej,format='(F12.2)'),2)+$
      '\pm'+strtrim(string(lostats.sigma_rej,format='(F12.2)'),2)+')'

    sdss_lineplot, x, y, xerr, yerr, plottype=1, postscript=postscript, $
      xtitle='', ytitle='', xrange=xrange, yrange=yrange, legendtype=0, $
      /right, /top, position=pos[*,2], charsize=charsize_2, ytickname=replicate(' ',10), $
      /noerase
    djs_oplot, !x.crange, !y.crange, thick=postthick, color=djs_icolor(talkcolor), line=2

    legend, textoidl(upxstr+' (Upper)'), /left, /top, box=0, $
      charsize=charsize_1, charthick=postthick, textcolor=djs_icolor(talkcolor)
;   legend, textoidl('(c)'), /left, /top, box=0, $
;     charsize=charsize_2, charthick=postthick, textcolor=djs_icolor(talkcolor)

    im_openclose, postscript=postscript, /close    

; ------------------------------------------------------------
; inter-compare empirical abundances [HII,Integrated,SDSS]
; ------------------------------------------------------------

    psname = '12oh_empirical_intercompare_6panel'
    im_openclose, pspath+psname, postscript=postscript, encapsulated=encapsulated, xsize=13.7, ysize=9.4

    pagemaker, nx=3, ny=2, xspace=0.0, yspace=0.0, width=4.0*[1,1,1], $
      height=4.0*[1,1], xmargin=[1.3,0.4], ymargin=[0.4,1.0], $
      xpage=13.7, ypage=9.4, position=pos, /normal

;   if keyword_set(dotalk) then polyfill, [1,1,0,0,1], [1,0,0,1,1], /normal, $
;     color=djs_icolor('black')

    xrange = ohrange4
    yrange = xrange

; --------------------------------------------------    
; Panel 1 - HII Regions - 12+log(O/H) [N II]/Ha vs 12+log(O/H) ([O III]/Hb)/([N II]/Ha)
; --------------------------------------------------    

    indx = where((hii.zstrong_12oh_niiha_pettini gt -900.0) and $
      (hii.zstrong_12oh_oiiinii_pettini gt -900.0),nindx)
 
    x = hii[indx].zstrong_12oh_niiha_pettini
    xerr = hii[indx].zstrong_12oh_niiha_pettini_err
 
    y = hii[indx].zstrong_12oh_oiiinii_pettini
    yerr = hii[indx].zstrong_12oh_oiiinii_pettini_err

    xtitle = '12+log(O/H) [N2]'   ; [N II]/H\alpha'
    ytitle = '12+log(O/H) [O3N2]' ; ([O III]/H\beta)/([N II]/H\alpha)'

    stats = im_stats(y-x,/verbose)
    xstr = strtrim(string(stats.median_rej,format='(F12.2)'),2)+' ('+$
      strtrim(string(stats.mean_rej,format='(F12.2)'),2)+$
      '\pm'+strtrim(string(stats.sigma_rej,format='(F12.2)'),2)+')'

    atlas1d_lineplot, x, y, xerr, yerr, plottype=1, postscript=postscript, $
      xtitle='', ytitle=ytitle, xrange=xrange, yrange=yrange, legendtype=0, charsize=charsize_6, $
      xstyle=3, ystyle=3, /right, /top, position=pos[*,0], xtickname=replicate(' ',10), hiipsize=0.5, $
      hiicolor='purple'
    djs_oplot, !x.crange, !y.crange, thick=postthick, color=djs_icolor(talkcolor), line=2
 
    legend, textoidl(xstr), /right, /bottom, box=0, charsize=charsize_6, $
      charthick=postthick, textcolor=djs_icolor(talkcolor)
    legend, '(a)', /left, /top, box=0, charsize=charsize_6, charthick=postthick, $
      textcolor=djs_icolor(talkcolor)

; --------------------------------------------------    
; Panel 2 - Integrated - 12+log(O/H) [N II]/Ha vs 12+log(O/H) ([O III]/Hb)/([N II]/Ha)
; --------------------------------------------------    

; Atlas     
    
    indx = where((intnodust.zstrong_12oh_niiha_pettini gt -900.0) and $
      (intnodust.zstrong_12oh_oiiinii_pettini gt -900.0),nindx)
 
    x = intnodust[indx].zstrong_12oh_niiha_pettini
    xerr = intnodust[indx].zstrong_12oh_niiha_pettini_err
 
    y = intnodust[indx].zstrong_12oh_oiiinii_pettini
    yerr = intnodust[indx].zstrong_12oh_oiiinii_pettini_err

    xtitle = '12+log(O/H) [N2]'   ; [N II]/H\alpha'
    ytitle = '12+log(O/H) [O3N2]' ; ([O III]/H\beta)/([N II]/H\alpha)'

    stats = im_stats(y-x,/verbose,/no_head)
    xstr = strtrim(string(stats.median_rej,format='(F12.2)'),2)+' ('+$
      strtrim(string(stats.mean_rej,format='(F12.2)'),2)+$
      '\pm'+strtrim(string(stats.sigma_rej,format='(F12.2)'),2)+')'
    
    atlas1d_lineplot, x, y, xerr, yerr, plottype=1, postscript=postscript, $
      xtitle='', ytitle='', xrange=xrange, yrange=yrange, legendtype=0, charsize=charsize_6, $
      xstyle=3, ystyle=3, /right, /top, position=pos[*,1], /noerase, ytickname=replicate(' ',10), $
      xtickname=replicate(' ',10)
    djs_oplot, !x.crange, !y.crange, thick=postthick, color=djs_icolor(talkcolor), line=2

    legend, textoidl(xstr), /right, /bottom, box=0, charsize=charsize_6, $
      charthick=postthick, textcolor=djs_icolor(talkcolor)
    legend, '(b)', /left, /top, box=0, charsize=charsize_6, charthick=postthick, $
      textcolor=djs_icolor(talkcolor)

; --------------------------------------------------    
; Panel 3 - SDSS - 12+log(O/H) [N II]/Ha vs 12+log(O/H) ([O III]/Hb)/([N II]/Ha)
; --------------------------------------------------    

    indx = where((sdssnodust.zstrong_12oh_niiha_pettini gt -900.0) and $
      (sdssnodust.zstrong_12oh_oiiinii_pettini gt -900.0),nindx)
 
    x = sdssnodust[indx].zstrong_12oh_niiha_pettini
    xerr = sdssnodust[indx].zstrong_12oh_niiha_pettini_err
 
    y = sdssnodust[indx].zstrong_12oh_oiiinii_pettini
    yerr = sdssnodust[indx].zstrong_12oh_oiiinii_pettini_err

    xtitle = '12+log(O/H) [N2]'   ; [N II]/H\alpha'
    ytitle = '12+log(O/H) [O3N2]' ; ([O III]/H\beta)/([N II]/H\alpha)'

    stats = im_stats(y-x,/verbose,/no_head)
    xstr = strtrim(string(stats.median_rej,format='(F12.2)'),2)+' ('+$
      strtrim(string(stats.mean_rej,format='(F12.2)'),2)+$
      '\pm'+strtrim(string(stats.sigma_rej,format='(F12.2)'),2)+')'
    
    sdss_lineplot, x, y, xerr, yerr, plottype=1, postscript=postscript, $
      xtitle='', ytitle='', xrange=xrange, yrange=yrange, legendtype=0, $
      xstyle=3, ystyle=3, /right, /top, position=pos[*,2], /noerase, charsize=charsize_6, $
      ytickname=replicate(' ',10), xtickname=replicate(' ',10), hiipsize=0.5
    djs_oplot, !x.crange, !y.crange, thick=postthick, color=djs_icolor(talkcolor), line=2
 
    legend, textoidl(xstr), /right, /bottom, box=0, charsize=charsize_6, $
      charthick=postthick, textcolor=djs_icolor(talkcolor)
    legend, '(c)', /left, /top, box=0, charsize=charsize_6, charthick=postthick, $
      textcolor=djs_icolor(talkcolor)

; --------------------------------------------------    
; Panel 4 - HII Regions - 12+log(O/H) [N II]/Ha vs 12+log(O/H) P-method
; --------------------------------------------------    

    indx = where((hii.zstrong_12oh_niiha_pettini gt -900.0) and $
      (hii.zstrong_12oh_pt05 gt -900.0) and $
      (hii.zstrong_12oh_pt05 gt -900.0),nindx)
 
    up = where((hii[indx].zstrong_12oh_pt05 gt 8.2) and (hii[indx].zstrong_12oh_niiha_pettini gt 8.2),nup)
    if (nup ne 0L) then begin
       xup = hii[indx[up]].zstrong_12oh_niiha_pettini
       xerrup = hii[indx[up]].zstrong_12oh_niiha_pettini_err

       yup = hii[indx[up]].zstrong_12oh_pt05
       yerrup = hii[indx[up]].zstrong_12oh_pt05_err

       x = xup & xerr = xerrup
       y = yup & yerr = yerrup
    endif

    lo = where((hii[indx].zstrong_12oh_pt05 lt 7.95) and (hii[indx].zstrong_12oh_niiha_pettini lt 8.2),nlo)
    if (nlo ne 0L) then begin
       xlo = hii[indx[lo]].zstrong_12oh_niiha_pettini
       xerrlo = hii[indx[lo]].zstrong_12oh_niiha_pettini_err

       ylo = hii[indx[lo]].zstrong_12oh_pt05
       yerrlo = hii[indx[lo]].zstrong_12oh_pt05_err

       x = [xup,xlo] & xerr = [xerrup,xerrlo]
       y = [yup,ylo] & yerr = [yerrup,yerrlo]
    endif

    stats = im_stats(y-x,/verbose,/no_head)
    xstr = strtrim(string(stats.median_rej,format='(F12.2)'),2)+' ('+$
      strtrim(string(stats.mean_rej,format='(F12.2)'),2)+$
      '\pm'+strtrim(string(stats.sigma_rej,format='(F12.2)'),2)+')'
    
    atlas1d_lineplot, x, y, xerr, yerr, plottype=1, postscript=postscript, $
      xtitle=xtitle, ytitle=ytitle, xrange=xrange, yrange=yrange, legendtype=0, $
      xstyle=3, ystyle=3, /right, /top, position=pos[*,3], /noerase, charsize=charsize_6, $
      hiipsize=0.5, hiicolor='purple'
    djs_oplot, !x.crange, !y.crange, thick=postthick, color=djs_icolor(talkcolor), line=2

    legend, textoidl(xstr), /right, /bottom, box=0, charsize=charsize_6, $
      charthick=postthick, textcolor=djs_icolor(talkcolor)
    legend, '(d)', /left, /top, box=0, charsize=charsize_6, charthick=postthick, $
      textcolor=djs_icolor(talkcolor)

; --------------------------------------------------    
; Panel 5 - Integrated - 12+log(O/H) [N II]/Ha vs 12+log(O/H) P-method
; --------------------------------------------------    

; Atlas
    
    indx = where((intnodust.zstrong_12oh_niiha_pettini gt -900.0) and $
      (intnodust.zstrong_12oh_pt05 gt -900.0) and $
      (intnodust.zstrong_12oh_pt05 gt -900.0),nindx)

    xtitle = '12+log(O/H) [N2]'
    ytitle = '12+log(O/H) [PT05]'

    up = where((intnodust[indx].zstrong_12oh_pt05 gt 8.2) and $
      (intnodust[indx].zstrong_12oh_niiha_pettini gt 8.2),nup)
    if (nup ne 0L) then begin
       xup = intnodust[indx[up]].zstrong_12oh_niiha_pettini
       xerrup = intnodust[indx[up]].zstrong_12oh_niiha_pettini_err

       yup = intnodust[indx[up]].zstrong_12oh_pt05
       yerrup = intnodust[indx[up]].zstrong_12oh_pt05_err

       x = xup & xerr = xerrup
       y = yup & yerr = yerrup
    endif

    lo = where((intnodust[indx].zstrong_12oh_pt05 lt 7.95) and $
      (intnodust[indx].zstrong_12oh_niiha_pettini lt 8.2),nlo)
    if (nlo ne 0L) then begin
       xlo = intnodust[indx[lo]].zstrong_12oh_niiha_pettini
       xerrlo = intnodust[indx[lo]].zstrong_12oh_niiha_pettini_err

       ylo = intnodust[indx[lo]].zstrong_12oh_pt05
       yerrlo = intnodust[indx[lo]].zstrong_12oh_pt05_err

       x = [xup,xlo] & xerr = [xerrup,xerrlo]
       y = [yup,ylo] & yerr = [yerrup,yerrlo]
    endif

    stats = im_stats(y-x,/verbose,/no_head)
    xstr = strtrim(string(stats.median_rej,format='(F12.2)'),2)+' ('+$
      strtrim(string(stats.mean_rej,format='(F12.2)'),2)+$
      '\pm'+strtrim(string(stats.sigma_rej,format='(F12.2)'),2)+')'
    
    atlas1d_lineplot, x, y, xerr, yerr, plottype=1, postscript=postscript, $
      xtitle=xtitle, ytitle='', xrange=xrange, yrange=yrange, legendtype=0, charsize=charsize_6, $
      xstyle=3, ystyle=3, /right, /top, position=pos[*,4], /noerase, $
      ytickname=replicate(' ',10)
    djs_oplot, !x.crange, !y.crange, thick=postthick, color=djs_icolor(talkcolor), line=2

    legend, textoidl(xstr), /right, /bottom, box=0, charsize=charsize_6, $
      charthick=postthick, textcolor=djs_icolor(talkcolor)
    legend, '(e)', /left, /top, box=0, charsize=charsize_6, charthick=postthick, $
      textcolor=djs_icolor(talkcolor)

; --------------------------------------------------    
; Panel 6 - SDSS - 12+log(O/H) [N II]/Ha vs 12+log(O/H) P-method
; --------------------------------------------------    

    indx = where((sdssnodust.zstrong_12oh_niiha_pettini gt -900.0) and $
      (sdssnodust.zstrong_12oh_pt05 gt -900.0) and $
      (sdssnodust.zstrong_12oh_pt05 gt -900.0),nindx)
 
    up = where((sdssnodust[indx].zstrong_12oh_pt05 gt 8.2) and $
      (sdssnodust[indx].zstrong_12oh_niiha_pettini gt 8.2),nup)
    if (nup ne 0L) then begin
       xup = sdssnodust[indx[up]].zstrong_12oh_niiha_pettini
       xerrup = sdssnodust[indx[up]].zstrong_12oh_niiha_pettini_err

       yup = sdssnodust[indx[up]].zstrong_12oh_pt05
       yerrup = sdssnodust[indx[up]].zstrong_12oh_pt05_err

       x = xup & xerr = xerrup
       y = yup & yerr = yerrup
    endif

    lo = where((sdssnodust[indx].zstrong_12oh_pt05 lt 7.95) and $
      (sdssnodust[indx].zstrong_12oh_niiha_pettini lt 8.2),nlo)
    if (nlo ne 0L) then begin
       xlo = sdssnodust[indx[lo]].zstrong_12oh_niiha_pettini
       xerrlo = sdssnodust[indx[lo]].zstrong_12oh_niiha_pettini_err

       ylo = sdssnodust[indx[lo]].zstrong_12oh_pt05
       yerrlo = sdssnodust[indx[lo]].zstrong_12oh_pt05_err

       x = [xup,xlo] & xerr = [xerrup,xerrlo]
       y = [yup,ylo] & yerr = [yerrup,yerrlo]
    endif

    stats = im_stats(y-x,/verbose,/no_head)
    xstr = strtrim(string(stats.median_rej,format='(F12.2)'),2)+' ('+$
      strtrim(string(stats.mean_rej,format='(F12.2)'),2)+$
      '\pm'+strtrim(string(stats.sigma_rej,format='(F12.2)'),2)+')'
    
    sdss_lineplot, x, y, xerr, yerr, plottype=1, postscript=postscript, $
      xtitle=xtitle, ytitle='', xrange=xrange, yrange=yrange, legendtype=0, charsize=charsize_6, $
      xstyle=3, ystyle=3, /right, /top, position=pos[*,5], /noerase, $
      ytickname=replicate(' ',10)

    djs_oplot, !x.crange, !y.crange, thick=postthick, color=djs_icolor(talkcolor), line=2

    legend, textoidl(xstr), /right, /bottom, box=0, charsize=charsize_6, $
      charthick=postthick, textcolor=djs_icolor(talkcolor)
    legend, '(f)', /left, /top, box=0, charsize=charsize_6, charthick=postthick, $
      textcolor=djs_icolor(talkcolor)

    im_openclose, postscript=postscript, /close    

; ------------------------------------------------------------
; O32 vs 12+log(O/H) [O3N2]+constant
; ------------------------------------------------------------

    psname = 'o32_vs_12oh_o3n2'
    im_openclose, pspath+psname, postscript=postscript, xsize=8.5, ysize=8.5, encapsulated=encapsulated

    pagemaker, nx=1, ny=1, xspace=0, yspace=0.0, width=7.1, height=7.1, $
      xmargin=[1.1,0.3], ymargin=[0.3,1.1], xpage=8.5, ypage=8.5, $
      position=pos, /normal

    indx = where((sdssdust.zstrong_12oh_oiiinii_niiha gt -900.0) and (sdssnodust.zstrong_o32 gt -900.0),nindx)

    x = sdssnodust[indx].zstrong_o32
    xerr = sdssnodust[indx].zstrong_o32_err

    y = sdssdust[indx].zstrong_12oh_oiiinii_niiha+offset
    yerr = sdssdust[indx].zstrong_12oh_oiiinii_niiha_err

    xtitle = 'log O_{32}'
    ytitle = '12 + log (O/H) [O3N2] + '+string(offset,format='(F4.2)')

    xrange = o32range
    yrange = ohrange11

    if keyword_set(hiiregions) then begin
       good = where((hii.zstrong_o32 gt -900.0) and (hii.zstrong_12oh_oiiinii_niiha gt -900))
       xregion = hii[good].zstrong_o32
       xerrregion = hii[good].zstrong_o32
       yregion = hii[good].zstrong_12oh_oiiinii_niiha+offset
       yerrregion = hii[good].zstrong_12oh_oiiinii_niiha_err
    endif
       
    sdss_lineplot, x, y, xerr, yerr, plottype=1, postscript=postscript, $
      xtitle=xtitle, ytitle=ytitle, xrange=xrange, yrange=yrange, legendtype=0, $
      /right, /top, charsize=charsize_8, xregion=xregion, yregion=yregion, $
      xerrregion=xerrregion, yerrregion=yerrregion, position=pos[*,0]

    legend, textoidl(xstr), /right, /top, box=0, charsize=charsize_6, charthick=postthick

    im_openclose, postscript=postscript, /close    

stop    

return
end    
