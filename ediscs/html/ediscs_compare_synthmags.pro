pro ediscs_compare_synthmags, ediscs, encapsulated=encapsulated, postscript=postscript, cleanpng=cleanpng
; jm07apr16nyu - compare synthesized and photometric magnitudes

; rsync -auv --delete html/ howdy:"public_html/research/ediscs"

    if (n_elements(ediscs) eq 0L) then ediscs = ediscs_read_ancillary()
    
    if keyword_set(postscript) then begin
       postthick = 4.0
       postthick2 = 2.0
    endif else begin
       postthick = 2.0
       postthick2 = 1.0
       im_window, 0, xratio=0.45, /square
    endelse

    htmlbase = 'compare_synthmags'
    html_path = ediscs_path(/html)
    pspath = html_path+htmlbase+'/'

; clean up old PS and PNG files

    if keyword_set(cleanpng) then begin
       splog, 'Deleting all PNG and PS files in '+pspath
       spawn, ['rm -f '+pspath+'*.png*'], /sh
       spawn, ['rm -f '+pspath+'*.*ps'], /sh
    endif
    
; ###########################################################################
; PAPER PLOTS
; ###########################################################################

; ---------------------------------------------------------------------------    

    psname = 'synth_vr_vs_phot_vr_paper'
    im_openclose, pspath+psname, postscript=postscript, xsize=6.9, ysize=7.5, encapsulated=encapsulated
    pagemaker, nx=1, ny=2, xspace=0, yspace=1.1, width=5.0, height=2.5*[1,1], $
      xmargin=[1.5,0.4], ymargin=[0.3,1.1], xpage=6.9, ypage=7.5, $
      position=pos, /normal

    indx = where((ediscs.synth_mobs_ab[2] gt 0.0) and (ediscs.kcorr_mobs_ab_aper[2] gt 0.0) and $
      (ediscs.synth_mobs_ab[1] gt 0.0) and (ediscs.kcorr_mobs_ab_aper[1] gt 0.0) and $
      (ediscs.minwave le 4800.0) and (ediscs.maxwave ge 7500.0),nindx)
    flag = where((ediscs[indx].targetflag eq 3) or (ediscs[indx].photsplit eq 1),comp=good,ncomp=ngood)
    run3 = where(strmatch(ediscs[indx[good]].specfile,'*m0[5-8]*'),nrun3)
    run4 = where(strmatch(ediscs[indx[good]].specfile,'*m09*') or strmatch(ediscs[indx[good]].specfile,'*m10*') or $
      strmatch(ediscs[indx[good]].specfile,'*m11*'),nrun4)

    x = ediscs[indx].kcorr_mobs_ab_aper[1] - ediscs[indx].kcorr_mobs_ab_aper[2]
    y = ediscs[indx].synth_mobs_ab[1] - ediscs[indx].synth_mobs_ab[2]
    xerr = sqrt(ediscs[indx].kcorr_mobs_ab_aper_err[1]^2+ediscs[indx].kcorr_mobs_ab_aper_err[2]^2)

    x2 = x
    y2 = y-x
    x2err = xerr
    y2err = xerr

    x3 = ediscs[indx].kcorr_mobs_ab_aper[2]
    y3 = y2
    x3err = ediscs[indx].kcorr_mobs_ab_aper_err[2]
    y3err = xerr
    y3err = xerr
    
    stats = im_stats(y[good]-x[good],sigrej=3.0)
    stats_run3 = im_stats(y[good[run3]]-x[good[run3]],sigrej=3.0)
    stats_run4 = im_stats(y[good[run4]]-x[good[run4]],sigrej=3.0)
    xstr = strtrim(string(stats.mean_rej,format='(F12.3)'),2)+$
      '\pm'+strtrim(string(stats.sigma_rej,format='(F12.3)'),2)+$
      ' ('+strtrim(string(stats.median_rej,format='(F12.3)'),2)+')'
    xstr_run3 = strtrim(string(stats_run3.mean_rej,format='(F12.3)'),2)+$
      '\pm'+strtrim(string(stats_run3.sigma_rej,format='(F12.3)'),2)+$
      ' ('+strtrim(string(stats_run3.median_rej,format='(F12.3)'),2)+')'
    xstr_run4 = strtrim(string(stats_run4.mean_rej,format='(F12.3)'),2)+$
      '\pm'+strtrim(string(stats_run4.sigma_rej,format='(F12.3)'),2)+$
      ' ('+strtrim(string(stats_run4.median_rej,format='(F12.3)'),2)+')'
    splog, psname+' - All:  '+xstr+     ', # = '+string(stats.npts_rej,format='(I0)')+'/'+string(ngood,format='(I0)')
    splog, psname+' - Run3: '+xstr_run3+', # = '+string(stats_run3.npts_rej,format='(I0)')+'/'+string(nrun3,format='(I0)')
    splog, psname+' - Run4: '+xstr_run4+', # = '+string(stats_run4.npts_rej,format='(I0)')+'/'+string(nrun4,format='(I0)')

    xtitle = '(V-R)_{phot}'
    ytitle = '(V-R)_{synth}'
    xrange = [0.0,1.3]
    yrange = xrange

    xtitle2 = xtitle
    ytitle2 = ytitle+'-'+xtitle
    xrange2 = xrange
    yrange2 = [-0.5,0.5]

    xtitle3 = 'R_{phot}'
    ytitle3 = ytitle2
    xrange3 = [20.5,24.5]
    yrange3 = yrange2

    xpivot = mean(x2[good])
    coeff = linfit(x2[good]-xpivot,y2[good],measure_errors=y2err[good],sigma=coeff_err)
    coeff[0] = coeff[0] - coeff[1]*xpivot
    xfit = findgen(((2.0)-(-1.0))/0.01)*0.01+(-1.0)
    yfit = poly(xfit,coeff)
    fitstr = '\Delta(V-R) = ('+strtrim(string(coeff[0],format='(F12.3)'),2)+'\pm'+$
      strtrim(string(coeff_err[0],format='(F12.3)'),2)+') + ('+$
      strtrim(string(coeff[1],format='(F12.3)'),2)+'\pm'+strtrim(string(coeff_err[1],format='(F12.3)'),2)+$
      ') (V-R)_{phot}'
    splog, psname, fitstr
    
    djs_plot, [0], [0], /nodata, xsty=3, ysty=3, xrange=xrange2, yrange=yrange2, $
      xtitle=xtitle2, ytitle=ytitle2, charsize=2.0, charthick=postthick, $
      xthick=postthick, ythick=postthick, position=pos[*,0]
    plotsym, 0, 1.2, fill=1
    djs_oplot, x2[good], y2[good], ps=8
    if (nrun4 ne 0L) then djs_oplot, x2[good[run4]], y2[good[run4]], ps=7, sym=2, thick=postthick, color='red'
    plotsym, 0, 1.2, fill=0, thick=postthick, color=djs_icolor('red')
    djs_oplot, !x.crange, [0,0], line=0, thick=postthick
    legend, '(a)', /left, /top, box=0, charsize=2.0, charthick=postthick, margin=0

    djs_plot, [0], [0], /nodata, /noerase, xsty=3, ysty=3, xrange=xrange3, yrange=yrange3, $
      xtitle=xtitle3, ytitle=ytitle3, charsize=2.0, charthick=postthick, $
      xthick=postthick, ythick=postthick, position=pos[*,1]
    plotsym, 0, 1.2, fill=1
    djs_oplot, x3[good], y3[good], ps=8
    if (nrun4 ne 0L) then djs_oplot, x3[good[run4]], y3[good[run4]], ps=7, sym=2, thick=postthick, color='red'
    plotsym, 0, 1.2, fill=0, thick=postthick, color=djs_icolor('red')
    djs_oplot, !x.crange, [0,0], line=0, thick=postthick
    legend, '(b)', /left, /top, box=0, charsize=2.0, charthick=postthick, margin=0

    im_openclose, postscript=postscript, /close

; ---------------------------------------------------------------------------    
    
    psname = 'synth_ri_vs_phot_ri_paper'
    im_openclose, pspath+psname, postscript=postscript, xsize=6.9, ysize=7.5, encapsulated=encapsulated
    pagemaker, nx=1, ny=2, xspace=0, yspace=1.1, width=5.0, height=2.5*[1,1], $
      xmargin=[1.5,0.4], ymargin=[0.3,1.1], xpage=6.9, ypage=7.5, $
      position=pos, /normal

    indx = where((ediscs.synth_mobs_ab[2] gt 0.0) and (ediscs.kcorr_mobs_ab_aper[2] gt 0.0) and $
      (ediscs.synth_mobs_ab[3] gt 0.0) and (ediscs.kcorr_mobs_ab_aper[3] gt 0.0) and $
      (ediscs.minwave le 5700.0) and (ediscs.maxwave ge 8700.0),nindx)
    flag = where((ediscs[indx].targetflag eq 3) or (ediscs[indx].photsplit eq 1),comp=good,ncomp=ngood)
    run3 = where(strmatch(ediscs[indx[good]].specfile,'*m0[5-8]*'),nrun3)
    run4 = where(strmatch(ediscs[indx[good]].specfile,'*m09*') or strmatch(ediscs[indx[good]].specfile,'*m10*') or $
      strmatch(ediscs[indx[good]].specfile,'*m11*'),nrun4)

    x = ediscs[indx].kcorr_mobs_ab_aper[2] - ediscs[indx].kcorr_mobs_ab_aper[3]
    y = ediscs[indx].synth_mobs_ab[2] - ediscs[indx].synth_mobs_ab[3]
    xerr = sqrt(ediscs[indx].kcorr_mobs_ab_aper_err[2]^2+ediscs[indx].kcorr_mobs_ab_aper_err[3]^2)

    x2 = x
    y2 = y-x
    x2err = xerr
    y2err = xerr

    x3 = ediscs[indx].kcorr_mobs_ab_aper[3]
    y3 = y2
    x3err = ediscs[indx].kcorr_mobs_ab_aper_err[3]
    y3err = xerr
    y3err = xerr

    stats = im_stats(y[good]-x[good],sigrej=3.0)
    stats_run3 = im_stats(y[good[run3]]-x[good[run3]],sigrej=3.0)
    stats_run4 = im_stats(y[good[run4]]-x[good[run4]],sigrej=3.0)
    xstr = strtrim(string(stats.mean_rej,format='(F12.3)'),2)+$
      '\pm'+strtrim(string(stats.sigma_rej,format='(F12.3)'),2)+$
      ' ('+strtrim(string(stats.median_rej,format='(F12.3)'),2)+')'
    xstr_run3 = strtrim(string(stats_run3.mean_rej,format='(F12.3)'),2)+$
      '\pm'+strtrim(string(stats_run3.sigma_rej,format='(F12.3)'),2)+$
      ' ('+strtrim(string(stats_run3.median_rej,format='(F12.3)'),2)+')'
    xstr_run4 = strtrim(string(stats_run4.mean_rej,format='(F12.3)'),2)+$
      '\pm'+strtrim(string(stats_run4.sigma_rej,format='(F12.3)'),2)+$
      ' ('+strtrim(string(stats_run4.median_rej,format='(F12.3)'),2)+')'
    splog, psname+' - All:  '+xstr+     ', # = '+string(stats.npts_rej,format='(I0)')+'/'+string(ngood,format='(I0)')
    splog, psname+' - Run3: '+xstr_run3+', # = '+string(stats_run3.npts_rej,format='(I0)')+'/'+string(nrun3,format='(I0)')
    splog, psname+' - Run4: '+xstr_run4+', # = '+string(stats_run4.npts_rej,format='(I0)')+'/'+string(nrun4,format='(I0)')

    xtitle = '(R-I)_{phot}'
    ytitle = '(R-I)_{synth}'
    xrange = [0.0,1.3]
    yrange = xrange

    xtitle2 = xtitle
    ytitle2 = ytitle+'-'+xtitle
    xrange2 = xrange
    yrange2 = [-0.5,0.5]

    xtitle3 = 'I_{phot}'
    ytitle3 = ytitle2
    xrange3 = [20.5,24.5] ; [20.5,24.0]
    yrange3 = yrange2

    xpivot = mean(x2[good])
    coeff = linfit(x2[good]-xpivot,y2[good],measure_errors=y2err[good],sigma=coeff_err)
    coeff[0] = coeff[0] - coeff[1]*xpivot
    xfit = findgen(((2.0)-(-1.0))/0.01)*0.01+(-1.0)
    yfit = poly(xfit,coeff)
    fitstr = '\Delta(R-I) = ('+strtrim(string(coeff[0],format='(F12.3)'),2)+'\pm'+$
      strtrim(string(coeff_err[0],format='(F12.3)'),2)+') + ('+$
      strtrim(string(coeff[1],format='(F12.3)'),2)+'\pm'+strtrim(string(coeff_err[1],format='(F12.3)'),2)+$
      ') (R-I)_{phot}'
    splog, psname, fitstr
    
    djs_plot, [0], [0], /nodata, xsty=3, ysty=3, xrange=xrange2, yrange=yrange2, $
      xtitle=xtitle2, ytitle=ytitle2, charsize=2.0, charthick=postthick, $
      xthick=postthick, ythick=postthick, position=pos[*,0]
    plotsym, 0, 1.2, fill=1
    djs_oplot, x2[good], y2[good], ps=8
    if (nrun4 ne 0L) then djs_oplot, x2[good[run4]], y2[good[run4]], ps=7, sym=2, thick=postthick, color='red'
    plotsym, 0, 1.2, fill=0, thick=postthick, color=djs_icolor('red')
    djs_oplot, !x.crange, [0,0], line=0, thick=postthick
    legend, '(c)', /left, /top, box=0, charsize=2.0, charthick=postthick, margin=0

    djs_plot, [0], [0], /nodata, /noerase, xsty=3, ysty=3, xrange=xrange3, yrange=yrange3, $
      xtitle=xtitle3, ytitle=ytitle3, charsize=2.0, charthick=postthick, $
      xthick=postthick, ythick=postthick, position=pos[*,1]
    plotsym, 0, 1.2, fill=1
    djs_oplot, x3[good], y3[good], ps=8
    if (nrun4 ne 0L) then djs_oplot, x3[good[run4]], y3[good[run4]], ps=7, sym=2, thick=postthick, color='red'
    plotsym, 0, 1.2, fill=0, thick=postthick, color=djs_icolor('red')
    djs_oplot, !x.crange, [0,0], line=0, thick=postthick
    legend, '(d)', /left, /top, box=0, charsize=2.0, charthick=postthick, margin=0
    
    im_openclose, postscript=postscript, /close

;   if keyword_set(postscript) and keyword_set(encapsulated) then im_ps2html, htmlbase, $
;     html_path=html_path, cleanpng=0, npscols=3, psfind=0

; ###########################################################################
; END PAPER PLOTS
; ###########################################################################
    
; ---------------------------------------------------------------------------    

    psname = 'synth_ri_vs_phot_ri'
    im_openclose, pspath+psname, postscript=postscript, xsize=8.5, ysize=11.0, encapsulated=encapsulated
    pagemaker, nx=1, ny=3, xspace=0, yspace=[1.0,1.0], width=6.6, height=[2.5,2.5,2.5], $
      xmargin=[1.5,0.4], ymargin=[0.4,1.1], xpage=8.5, ypage=11.0, $
      position=pos, /normal
;   im_openclose, pspath+psname, postscript=postscript, xsize=8.5, ysize=9.0, encapsulated=encapsulated
;   pagemaker, nx=1, ny=2, xspace=0, yspace=1.0, width=6.6, height=[3.6,3.0], $
;     xmargin=[1.5,0.4], ymargin=[0.3,1.1], xpage=8.5, ypage=9.0, $
;     position=pos, /normal

;   indx = where((ediscs.synth_mobs_ab[2] gt 0.0) and (ediscs.kcorr_mobs_ab[2] gt 0.0) and $
;     (ediscs.synth_mobs_ab[3] gt 0.0) and (ediscs.kcorr_mobs_ab[3] gt 0.0) and $
;     (ediscs.minwave gt 5700.0) and (ediscs.maxwave lt 8700.0),nindx)
    indx = where((ediscs.synth_mobs_ab[2] gt 0.0) and (ediscs.kcorr_mobs_ab_aper[2] gt 0.0) and $
      (ediscs.synth_mobs_ab[3] gt 0.0) and (ediscs.kcorr_mobs_ab_aper[3] gt 0.0),nindx)
    flag = where((ediscs[indx].targetflag eq 3) or (ediscs[indx].photsplit eq 1),comp=good,ncomp=ngood)

    x = ediscs[indx].kcorr_mobs_ab_aper[2] - ediscs[indx].kcorr_mobs_ab_aper[3]
    y = ediscs[indx].synth_mobs_ab[2] - ediscs[indx].synth_mobs_ab[3]
    xerr = sqrt(ediscs[indx].kcorr_mobs_ab_aper_err[2]^2+ediscs[indx].kcorr_mobs_ab_aper_err[3]^2)

    x2 = x
    y2 = y-x
    x2err = xerr
    y2err = xerr

    x3 = ediscs[indx].kcorr_mobs_ab_aper[3]
    y3 = y2
    x3err = ediscs[indx].kcorr_mobs_ab_aper_err[3]
    y3err = xerr
    y3err = xerr
    
    stats = im_stats(y[good]-x[good],sigrej=3.0)
    xstr = strtrim(string(stats.mean_rej,format='(F12.3)'),2)+$
      '\pm'+strtrim(string(stats.sigma_rej,format='(F12.3)'),2)+$
      ' ('+strtrim(string(stats.median_rej,format='(F12.3)'),2)+')'
    splog, psname+': '+xstr+', # = '+string(ngood,format='(I0)')

    xtitle = '(R-I)_{phot}'
    ytitle = '(R-I)_{synth}'
    xrange = minmax(x[good])
;   xrange = [min(x[good])<min(y[good]),max(x[good])>max(y[good])]
    yrange = xrange

    xtitle2 = xtitle
    ytitle2 = '\Delta(R-I)'
;   ytitle2 = xtitle+'-'+ytitle
;   ytitle2 = '\Delta(R-I)='+xtitle+'-'+ytitle
    xrange2 = xrange
    yrange2 = (min(abs(y2[good]))>max(abs(y2[good])))*[-1,1]
;   yrange2 = (min(abs(y2))>max(abs(y2)))*[-1,1]

    xtitle3 = 'I_{phot}'
    ytitle3 = ytitle2
    xrange3 = minmax(x3[good])
    yrange3 = yrange2

    xpivot = mean(x2[good])
    coeff = linfit(x2[good]-xpivot,y2[good],measure_errors=y2err[good],sigma=coeff_err)
    coeff[0] = coeff[0] - coeff[1]*xpivot
    xfit = findgen(((2.0)-(-1.0))/0.01)*0.01+(-1.0)
;   xfit = findgen((xrange2[1]-xrange2[0])/0.01)*0.01+xrange2[0]
    yfit = poly(xfit,coeff)
    
    djs_plot, [0], [0], /nodata, xsty=3, ysty=3, xrange=xrange, yrange=yrange, $
      xtitle=xtitle, ytitle=ytitle, charsize=2.0, charthick=postthick, $
      xthick=postthick, ythick=postthick, position=pos[*,0]
    plotsym, 0, 0.5, fill=1
    djs_oplot, x[good], y[good], ps=8, errthick=postthick2
;   oploterror, x[good], y[good], yerr[good], ps=3, errthick=postthick2
    plotsym, 0, 1.2, fill=0, thick=postthick, color=djs_icolor('red')
;   djs_oplot, x[flag], y[flag], ps=8
    djs_oplot, !x.crange, !y.crange, line=0, thick=postthick
    legend, '(a)', /left, /top, box=0, charsize=2.0, charthick=postthick
;   legend, textoidl('<\Delta(R-I)> = '+xstr+' mag'), /left, /top, box=0, $
;     charsize=1.6, charthick=postthick

    djs_plot, [0], [0], /nodata, /noerase, xsty=3, ysty=3, xrange=xrange2, yrange=yrange2, $
      xtitle=xtitle2, ytitle=ytitle2, charsize=2.0, charthick=postthick, $
      xthick=postthick, ythick=postthick, position=pos[*,1]
    plotsym, 0, 0.5, fill=1
    djs_oplot, x2[good], y2[good], ps=8
    plotsym, 0, 1.2, fill=0, thick=postthick, color=djs_icolor('red')
;   djs_oplot, x2[flag], y2[flag], ps=8
    djs_oplot, !x.crange, [0,0], line=0, thick=postthick
;   djs_oplot, xfit, yfit, line=2, thick=postthick
    legend, '(b)', /left, /top, box=0, charsize=2.0, charthick=postthick

;   legend, textoidl(['\Delta(R-I) = ('+strtrim(string(coeff[0],format='(F12.3)'),2)+'\pm'+$
;     strtrim(string(coeff_err[0],format='(F12.3)'),2)+') + ('+$
;     strtrim(string(coeff[1],format='(F12.3)'),2)+'\pm'+strtrim(string(coeff_err[1],format='(F12.3)'),2)+$
;     ') (R-I)_{phot}']), /right, /bottom, $
;     charsize=1.3, charthick=postthick, box=0, line=2, thick=postthick, /clear

    djs_plot, [0], [0], /nodata, /noerase, xsty=3, ysty=3, xrange=xrange3, yrange=yrange3, $
      xtitle=xtitle3, ytitle=ytitle3, charsize=2.0, charthick=postthick, $
      xthick=postthick, ythick=postthick, position=pos[*,2]
    plotsym, 0, 0.5, fill=1
    djs_oplot, x3[good], y3[good], ps=8
    plotsym, 0, 1.2, fill=0, thick=postthick, color=djs_icolor('red')
;   djs_oplot, x3[flag], y3[flag], ps=8
    djs_oplot, !x.crange, [0,0], line=0, thick=postthick
    legend, '(c)', /left, /top, box=0, charsize=2.0, charthick=postthick
    
    im_openclose, postscript=postscript, /close

; ---------------------------------------------------------------------------    
    
    psname = 'synth_ri_vs_phot_ri_subset'
    im_openclose, pspath+psname, postscript=postscript, xsize=8.5, ysize=11.0, encapsulated=encapsulated
    pagemaker, nx=1, ny=3, xspace=0, yspace=[1.0,1.0], width=6.6, height=[2.5,2.5,2.5], $
      xmargin=[1.5,0.4], ymargin=[0.4,1.1], xpage=8.5, ypage=11.0, $
      position=pos, /normal

    indx = where((ediscs.synth_mobs_ab[2] gt 0.0) and (ediscs.kcorr_mobs_ab_aper[2] gt 0.0) and $
      (ediscs.synth_mobs_ab[3] gt 0.0) and (ediscs.kcorr_mobs_ab_aper[3] gt 0.0) and $
      (ediscs.minwave le 5700.0) and (ediscs.maxwave ge 8700.0),nindx)
    flag = where((ediscs[indx].targetflag eq 3) or (ediscs[indx].photsplit eq 1),comp=good,ncomp=ngood)

    x = ediscs[indx].kcorr_mobs_ab_aper[2] - ediscs[indx].kcorr_mobs_ab_aper[3]
    y = ediscs[indx].synth_mobs_ab[2] - ediscs[indx].synth_mobs_ab[3]
    xerr = sqrt(ediscs[indx].kcorr_mobs_ab_aper_err[2]^2+ediscs[indx].kcorr_mobs_ab_aper_err[3]^2)

    x2 = x
    y2 = y-x
    x2err = xerr
    y2err = xerr

    x3 = ediscs[indx].kcorr_mobs_ab_aper[3]
    y3 = y2
    x3err = ediscs[indx].kcorr_mobs_ab_aper_err[3]
    y3err = xerr
    y3err = xerr

;   w = where(x3[good] lt 22.6)
;   stats = im_stats(y[good[w]]-x[good[w]],sigrej=10.0,mask=mask)
    stats = im_stats(y[good]-x[good],sigrej=3.0,mask=mask)
    xstr = strtrim(string(stats.mean_rej,format='(F12.3)'),2)+$
      '\pm'+strtrim(string(stats.sigma_rej,format='(F12.3)'),2)+$
      ' ('+strtrim(string(stats.median_rej,format='(F12.3)'),2)+')'
    splog, psname+': '+xstr+', # = '+string(ngood,format='(I0)')

    xtitle = '(R-I)_{phot}'
    ytitle = '(R-I)_{synth}'
    xrange = [0.0,1.2]
;   xrange = minmax(x[good])
;   xrange = [min(x)<min(y),max(x)>max(y)] ; FROM ABOVE!
    yrange = xrange

    xtitle2 = xtitle
    ytitle2 = '\Delta(R-I)'
;   ytitle2 = xtitle+'-'+ytitle
;   ytitle2 = '\Delta(R-I)='+xtitle+'-'+ytitle
    xrange2 = xrange
    yrange2 = [-0.6,0.6]
;   yrange2 = (min(abs(y2[good]))>max(abs(y2[good])))*[-1,1]
;   xrange2 = xrange ; FROM ABOVE!
;   yrange2 = (min(abs(y2[good]))>max(abs(y2[good])))*[-1,1] ; FROM ABOVE!
;   yrange2 = (min(abs(y2))>max(abs(y2)))*[-1,1]

    xtitle3 = 'I_{phot}'
    ytitle3 = ytitle2
    xrange3 = minmax(x3[good])
    yrange3 = yrange2
;   xrange3 = minmax(x3) ; FROM ABOVE!
;   yrange3 = yrange2    ; FROM ABOVE!

    xpivot = mean(x2[good])
    coeff = linfit(x2[good]-xpivot,y2[good],measure_errors=y2err[good],sigma=coeff_err)
    coeff[0] = coeff[0] - coeff[1]*xpivot
    xfit = findgen(((2.0)-(-1.0))/0.01)*0.01+(-1.0)
;   xfit = findgen((xrange2[1]-xrange2[0])/0.01)*0.01+xrange2[0]
    yfit = poly(xfit,coeff)
    fitstr = '\Delta(R-I) = ('+strtrim(string(coeff[0],format='(F12.3)'),2)+'\pm'+$
      strtrim(string(coeff_err[0],format='(F12.3)'),2)+') + ('+$
      strtrim(string(coeff[1],format='(F12.3)'),2)+'\pm'+strtrim(string(coeff_err[1],format='(F12.3)'),2)+$
      ') (R-I)_{phot}'
    splog, psname, fitstr
    
    djs_plot, [0], [0], /nodata, xsty=3, ysty=3, xrange=xrange, yrange=yrange, $
      xtitle=xtitle, ytitle=ytitle, charsize=2.0, charthick=postthick, $
      xthick=postthick, ythick=postthick, position=pos[*,0]
    plotsym, 0, 1.2, fill=1
    djs_oplot, x[good], y[good], ps=8, errthick=postthick2
;   oploterror, x[good], y[good], yerr[good], ps=3, errthick=postthick2
    plotsym, 0, 1.2, fill=0, thick=postthick, color=djs_icolor('red')
;   djs_oplot, x[flag], y[flag], ps=8
    djs_oplot, !x.crange, !y.crange, line=0, thick=postthick
    legend, '(a)', /left, /top, box=0, charsize=2.0, charthick=postthick
;   legend, textoidl('<\Delta(R-I)> = '+xstr+' mag'), /left, /top, box=0, $
;     charsize=1.6, charthick=postthick
;   legend, textoidl('5700 \AA < \lambda < 8700 \AA'), /right, /bottom, box=0, $
;     charsize=1.6, charthick=postthick

    djs_plot, [0], [0], /nodata, /noerase, xsty=3, ysty=3, xrange=xrange2, yrange=yrange2, $
      xtitle=xtitle2, ytitle=ytitle2, charsize=2.0, charthick=postthick, $
      xthick=postthick, ythick=postthick, position=pos[*,1]
    plotsym, 0, 1.2, fill=1
    djs_oplot, x2[good], y2[good], ps=8
    plotsym, 0, 1.2, fill=0, thick=postthick, color=djs_icolor('red')
;   djs_oplot, x2[flag], y2[flag], ps=8
    djs_oplot, !x.crange, [0,0], line=0, thick=postthick
;   djs_oplot, xfit, yfit, line=2, thick=postthick
    legend, '(b)', /left, /top, box=0, charsize=2.0, charthick=postthick
    legend, textoidl(fitstr), /right, /bottom, $
      charsize=1.3, charthick=postthick, box=0, line=2, thick=postthick;, clear=1

    djs_plot, [0], [0], /nodata, /noerase, xsty=3, ysty=3, xrange=xrange3, yrange=yrange3, $
      xtitle=xtitle3, ytitle=ytitle3, charsize=2.0, charthick=postthick, $
      xthick=postthick, ythick=postthick, position=pos[*,2]
    plotsym, 0, 1.2, fill=1
    djs_oplot, x3[good], y3[good], ps=8
;   djs_oplot, x3[good[where(mask eq 0)]], y3[good[where(mask eq 0)]], ps=7, color='red', thick=3, sym=2.0
    plotsym, 0, 1.2, fill=0, thick=postthick, color=djs_icolor('red')
;   djs_oplot, x3[flag], y3[flag], ps=8
    djs_oplot, !x.crange, [0,0], line=0, thick=postthick
    legend, '(c)', /left, /top, box=0, charsize=2.0, charthick=postthick
    
    im_openclose, postscript=postscript, /close

; ---------------------------------------------------------------------------    
    
    psname = 'synth_vr_vs_phot_vr'
    im_openclose, pspath+psname, postscript=postscript, xsize=8.5, ysize=11.0, encapsulated=encapsulated
    pagemaker, nx=1, ny=3, xspace=0, yspace=[1.0,1.0], width=6.6, height=[2.5,2.5,2.5], $
      xmargin=[1.5,0.4], ymargin=[0.4,1.1], xpage=8.5, ypage=11.0, $
      position=pos, /normal

    indx = where((ediscs.synth_mobs_ab[2] gt 0.0) and (ediscs.kcorr_mobs_ab_aper[2] gt 0.0) and $
      (ediscs.synth_mobs_ab[1] gt 0.0) and (ediscs.kcorr_mobs_ab_aper[1] gt 0.0),nindx)
    flag = where((ediscs[indx].targetflag eq 3) or (ediscs[indx].photsplit eq 1),comp=good,ncomp=ngood)

    x = ediscs[indx].kcorr_mobs_ab_aper[1] - ediscs[indx].kcorr_mobs_ab_aper[2]
    y = ediscs[indx].synth_mobs_ab[1] - ediscs[indx].synth_mobs_ab[2]
    xerr = sqrt(ediscs[indx].kcorr_mobs_ab_aper_err[1]^2+ediscs[indx].kcorr_mobs_ab_aper_err[2]^2)

    x2 = x
    y2 = y-x
    x2err = xerr
    y2err = xerr

    x3 = ediscs[indx].kcorr_mobs_ab_aper[2]
    y3 = y2
    x3err = ediscs[indx].kcorr_mobs_ab_aper_err[2]
    y3err = xerr
    y3err = xerr
    
    stats = im_stats(y[good]-x[good],sigrej=3.0)
    xstr = strtrim(string(stats.mean_rej,format='(F12.3)'),2)+$
      '\pm'+strtrim(string(stats.sigma_rej,format='(F12.3)'),2)+$
      ' ('+strtrim(string(stats.median_rej,format='(F12.3)'),2)+')'
    splog, psname+': '+xstr+', # = '+string(ngood,format='(I0)')

    xtitle = '(V-R)_{phot}'
    ytitle = '(V-R)_{synth}'
    xrange = minmax(x[good])
;   xrange = [min(x[good])<min(y[good]),max(x[good])>max(y[good])]
    yrange = xrange

    xtitle2 = xtitle
    ytitle2 = '\Delta(V-R)'
    xrange2 = xrange
    yrange2 = (min(abs(y2[good]))>max(abs(y2[good])))*[-1,1]

    xtitle3 = 'R_{phot}'
    ytitle3 = ytitle2
    xrange3 = minmax(x3[good])
    yrange3 = yrange2

    djs_plot, [0], [0], /nodata, xsty=3, ysty=3, xrange=xrange, yrange=yrange, $
      xtitle=xtitle, ytitle=ytitle, charsize=2.0, charthick=postthick, $
      xthick=postthick, ythick=postthick, position=pos[*,0]
    plotsym, 0, 0.5, fill=1
    djs_oplot, x[good], y[good], ps=8, errthick=postthick2
;   oploterror, x[good], y[good], yerr[good], ps=3, errthick=postthick2
    plotsym, 0, 1.2, fill=0, thick=postthick, color=djs_icolor('red')
;   djs_oplot, x[flag], y[flag], ps=8
    djs_oplot, !x.crange, !y.crange, line=0, thick=postthick
    legend, '(a)', /left, /top, box=0, charsize=2.0, charthick=postthick
;   legend, textoidl('<\Delta(V-R)> = '+xstr+' mag'), /left, /top, box=0, $
;     charsize=1.6, charthick=postthick

    djs_plot, [0], [0], /nodata, /noerase, xsty=3, ysty=3, xrange=xrange2, yrange=yrange2, $
      xtitle=xtitle2, ytitle=ytitle2, charsize=2.0, charthick=postthick, $
      xthick=postthick, ythick=postthick, position=pos[*,1]
    plotsym, 0, 0.5, fill=1
    djs_oplot, x2[good], y2[good], ps=8
    plotsym, 0, 1.2, fill=0, thick=postthick, color=djs_icolor('red')
;   djs_oplot, x2[flag], y2[flag], ps=8
    djs_oplot, !x.crange, [0,0], line=0, thick=postthick
    legend, '(b)', /left, /top, box=0, charsize=2.0, charthick=postthick

    djs_plot, [0], [0], /nodata, /noerase, xsty=3, ysty=3, xrange=xrange3, yrange=yrange3, $
      xtitle=xtitle3, ytitle=ytitle3, charsize=2.0, charthick=postthick, $
      xthick=postthick, ythick=postthick, position=pos[*,2]
    plotsym, 0, 0.5, fill=1
    djs_oplot, x3[good], y3[good], ps=8
    plotsym, 0, 1.2, fill=0, thick=postthick, color=djs_icolor('red')
;   djs_oplot, x3[flag], y3[flag], ps=8
    djs_oplot, !x.crange, [0,0], line=0, thick=postthick
    legend, '(c)', /left, /top, box=0, charsize=2.0, charthick=postthick

    im_openclose, postscript=postscript, /close

; ---------------------------------------------------------------------------    

    psname = 'synth_vr_vs_phot_vr_subset'
    im_openclose, pspath+psname, postscript=postscript, xsize=8.5, ysize=11.0, encapsulated=encapsulated
    pagemaker, nx=1, ny=3, xspace=0, yspace=[1.0,1.0], width=6.6, height=[2.5,2.5,2.5], $
      xmargin=[1.5,0.4], ymargin=[0.4,1.1], xpage=8.5, ypage=11.0, $
      position=pos, /normal

    indx = where((ediscs.synth_mobs_ab[2] gt 0.0) and (ediscs.kcorr_mobs_ab_aper[2] gt 0.0) and $
      (ediscs.synth_mobs_ab[1] gt 0.0) and (ediscs.kcorr_mobs_ab_aper[1] gt 0.0) and $
      (ediscs.minwave le 4800.0) and (ediscs.maxwave ge 7500.0),nindx)
    flag = where((ediscs[indx].targetflag eq 3) or (ediscs[indx].photsplit eq 1),comp=good,ncomp=ngood)

    x = ediscs[indx].kcorr_mobs_ab_aper[1] - ediscs[indx].kcorr_mobs_ab_aper[2]
    y = ediscs[indx].synth_mobs_ab[1] - ediscs[indx].synth_mobs_ab[2]
    xerr = sqrt(ediscs[indx].kcorr_mobs_ab_aper_err[1]^2+ediscs[indx].kcorr_mobs_ab_aper_err[2]^2)

    x2 = x
    y2 = y-x
    x2err = xerr
    y2err = xerr

    x3 = ediscs[indx].kcorr_mobs_ab_aper[2]
    y3 = y2
    x3err = ediscs[indx].kcorr_mobs_ab_aper_err[2]
    y3err = xerr
    y3err = xerr
    
    stats = im_stats(y[good]-x[good],sigrej=3.0)
    xstr = strtrim(string(stats.mean_rej,format='(F12.3)'),2)+$
      '\pm'+strtrim(string(stats.sigma_rej,format='(F12.3)'),2)+$
      ' ('+strtrim(string(stats.median_rej,format='(F12.3)'),2)+')'
    splog, psname+': '+xstr+', # = '+string(ngood,format='(I0)')

    xtitle = '(V-R)_{phot}'
    ytitle = '(V-R)_{synth}'
    xrange = [0.0,1.3]
;   xrange = minmax(x[good])
;   xrange = [min(x)<min(y),max(x)>max(y)]
    yrange = xrange

    xtitle2 = xtitle
    ytitle2 = '\Delta(V-R)'
    xrange2 = xrange
;   yrange2 = (min(abs(y2[good]))>max(abs(y2[good])))*[-1,1]
    yrange2 = [-0.6,0.6]

    xtitle3 = 'R_{phot}'
    ytitle3 = ytitle2
    xrange3 = minmax(x3[good])
    yrange3 = yrange2

    xpivot = mean(x2[good])
    coeff = linfit(x2[good]-xpivot,y2[good],measure_errors=y2err[good],sigma=coeff_err)
    coeff[0] = coeff[0] - coeff[1]*xpivot
    xfit = findgen(((2.0)-(-1.0))/0.01)*0.01+(-1.0)
;   xfit = findgen((xrange2[1]-xrange2[0])/0.01)*0.01+xrange2[0]
    yfit = poly(xfit,coeff)
    fitstr = '\Delta(V-R) = ('+strtrim(string(coeff[0],format='(F12.3)'),2)+'\pm'+$
      strtrim(string(coeff_err[0],format='(F12.3)'),2)+') + ('+$
      strtrim(string(coeff[1],format='(F12.3)'),2)+'\pm'+strtrim(string(coeff_err[1],format='(F12.3)'),2)+$
      ') (V-R)_{phot}'
    splog, psname, fitstr
    
    djs_plot, [0], [0], /nodata, xsty=3, ysty=3, xrange=xrange, yrange=yrange, $
      xtitle=xtitle, ytitle=ytitle, charsize=2.0, charthick=postthick, $
      xthick=postthick, ythick=postthick, position=pos[*,0]
    plotsym, 0, 1.2, fill=1
    djs_oplot, x[good], y[good], ps=8, errthick=postthick2
;   oploterror, x[good], y[good], yerr[good], ps=3, errthick=postthick2
    plotsym, 0, 1.2, fill=0, thick=postthick, color=djs_icolor('red')
;   djs_oplot, x[flag], y[flag], ps=8
    djs_oplot, !x.crange, !y.crange, line=0, thick=postthick
    legend, '(a)', /left, /top, box=0, charsize=2.0, charthick=postthick
;   legend, textoidl('<\Delta(V-R)> = '+xstr+' mag'), /left, /top, box=0, $
;     charsize=1.6, charthick=postthick

    djs_plot, [0], [0], /nodata, /noerase, xsty=3, ysty=3, xrange=xrange2, yrange=yrange2, $
      xtitle=xtitle2, ytitle=ytitle2, charsize=2.0, charthick=postthick, $
      xthick=postthick, ythick=postthick, position=pos[*,1]
    plotsym, 0, 1.2, fill=1
    djs_oplot, x2[good], y2[good], ps=8
    plotsym, 0, 1.2, fill=0, thick=postthick, color=djs_icolor('red')
;   djs_oplot, x2[flag], y2[flag], ps=8
    djs_oplot, !x.crange, [0,0], line=0, thick=postthick
    legend, '(b)', /left, /top, box=0, charsize=2.0, charthick=postthick

    djs_plot, [0], [0], /nodata, /noerase, xsty=3, ysty=3, xrange=xrange3, yrange=yrange3, $
      xtitle=xtitle3, ytitle=ytitle3, charsize=2.0, charthick=postthick, $
      xthick=postthick, ythick=postthick, position=pos[*,2]
    plotsym, 0, 1.2, fill=1
    djs_oplot, x3[good], y3[good], ps=8
    plotsym, 0, 1.2, fill=0, thick=postthick, color=djs_icolor('red')
;   djs_oplot, x3[flag], y3[flag], ps=8
    djs_oplot, !x.crange, [0,0], line=0, thick=postthick
    legend, '(c)', /left, /top, box=0, charsize=2.0, charthick=postthick

    im_openclose, postscript=postscript, /close

; ---------------------------------------------------------------------------    

    psname = 'synth_bv_vs_phot_bv'
    im_openclose, pspath+psname, postscript=postscript, xsize=8.5, ysize=11.0, encapsulated=encapsulated
    pagemaker, nx=1, ny=3, xspace=0, yspace=[1.0,1.0], width=6.6, height=[2.5,2.5,2.5], $
      xmargin=[1.5,0.4], ymargin=[0.4,1.1], xpage=8.5, ypage=11.0, $
      position=pos, /normal

    indx = where((ediscs.synth_mobs_ab[0] gt 0.0) and (ediscs.kcorr_mobs_ab_aper[0] gt 0.0) and $
      (ediscs.synth_mobs_ab[1] gt 0.0) and (ediscs.kcorr_mobs_ab_aper[1] gt 0.0),nindx)
    flag = where((ediscs[indx].targetflag eq 3) or (ediscs[indx].photsplit eq 1),comp=good,ncomp=ngood)
    
    x = ediscs[indx].kcorr_mobs_ab_aper[0] - ediscs[indx].kcorr_mobs_ab_aper[1]
    y = ediscs[indx].synth_mobs_ab[0] - ediscs[indx].synth_mobs_ab[1]
    xerr = sqrt(ediscs[indx].kcorr_mobs_ab_aper_err[0]^2+ediscs[indx].kcorr_mobs_ab_aper_err[1]^2)

    x2 = x
    y2 = y-x
    x2err = xerr
    y2err = xerr

    x3 = ediscs[indx].kcorr_mobs_ab_aper[1]
    y3 = y2
    x3err = ediscs[indx].kcorr_mobs_ab_aper_err[1]
    y3err = xerr
    y3err = xerr
    
    stats = im_stats(y[good]-x[good],sigrej=3.0)
    xstr = strtrim(string(stats.mean_rej,format='(F12.3)'),2)+$
      '\pm'+strtrim(string(stats.sigma_rej,format='(F12.3)'),2)+$
      ' ('+strtrim(string(stats.median_rej,format='(F12.3)'),2)+')'
    splog, psname+': '+xstr+', # = '+string(ngood,format='(I0)')

    xtitle = '(B-V)_{phot}'
    ytitle = '(B-V)_{synth}'
;   xrange = minmax(x[good])
    xrange = [min(x)<min(y),max(x)>max(y)]
    yrange = xrange

    xtitle2 = xtitle
    ytitle2 = '\Delta(B-V)'
    xrange2 = xrange
    yrange2 = (min(abs(y2[good]))>max(abs(y2[good])))*[-1,1]

    xtitle3 = 'V_{phot}'
    ytitle3 = ytitle2
    xrange3 = minmax(x3[good])
    yrange3 = yrange2

    djs_plot, [0], [0], /nodata, xsty=3, ysty=3, xrange=xrange, yrange=yrange, $
      xtitle=xtitle, ytitle=ytitle, charsize=2.0, charthick=postthick, $
      xthick=postthick, ythick=postthick, position=pos[*,0]
    plotsym, 0, 0.5, fill=1
    djs_oplot, x[good], y[good], ps=8, errthick=postthick2
;   oploterror, x[good], y[good], yerr[good], ps=3, errthick=postthick2
    plotsym, 0, 1.2, fill=0, thick=postthick, color=djs_icolor('red')
;   djs_oplot, x[flag], y[flag], ps=8
    djs_oplot, !x.crange, !y.crange, line=0, thick=postthick
    legend, '(a)', /left, /top, box=0, charsize=2.0, charthick=postthick
;   legend, textoidl('<\Delta(B-V)> = '+xstr+' mag'), /left, /top, box=0, $
;     charsize=1.6, charthick=postthick

    djs_plot, [0], [0], /nodata, /noerase, xsty=3, ysty=3, xrange=xrange2, yrange=yrange2, $
      xtitle=xtitle2, ytitle=ytitle2, charsize=2.0, charthick=postthick, $
      xthick=postthick, ythick=postthick, position=pos[*,1]
    plotsym, 0, 0.5, fill=1
    djs_oplot, x2[good], y2[good], ps=8
    plotsym, 0, 1.2, fill=0, thick=postthick, color=djs_icolor('red')
;   djs_oplot, x2[flag], y2[flag], ps=8
    djs_oplot, !x.crange, [0,0], line=0, thick=postthick
    legend, '(b)', /left, /top, box=0, charsize=2.0, charthick=postthick

    djs_plot, [0], [0], /nodata, /noerase, xsty=3, ysty=3, xrange=xrange3, yrange=yrange3, $
      xtitle=xtitle3, ytitle=ytitle3, charsize=2.0, charthick=postthick, $
      xthick=postthick, ythick=postthick, position=pos[*,2]
    plotsym, 0, 0.5, fill=1
    djs_oplot, x3[good], y3[good], ps=8
    plotsym, 0, 1.2, fill=0, thick=postthick, color=djs_icolor('red')
;   djs_oplot, x3[flag], y3[flag], ps=8
    djs_oplot, !x.crange, [0,0], line=0, thick=postthick
    legend, '(c)', /left, /top, box=0, charsize=2.0, charthick=postthick

    im_openclose, postscript=postscript, /close

;; ---------------------------------------------------------------------------    
;
;    psname = 'synth_b_vs_phot_b'
;    im_openclose, pspath+psname, postscript=postscript, xsize=8.5, ysize=8.0, encapsulated=encapsulated
;
;    pagemaker, nx=1, ny=1, xspace=0, yspace=0.0, width=6.6, height=6.6, $
;      xmargin=[1.5,0.4], ymargin=[0.3,1.1], xpage=8.5, ypage=8.0, $
;      position=pos, /normal
;
;    indx = where((ediscs.synth_mobs_ab[0] gt 0.0) and (ediscs.kcorr_mobs_ab_aper[0] gt 0.0),nindx)
;    x = ediscs[indx].synth_mobs_ab[0]
;    y = ediscs[indx].kcorr_mobs_ab_aper[0]
;    yerr = ediscs[indx].kcorr_mobs_ab_aper_err[0]
;    
;    stats = im_stats(y-x,sigrej=3.0)
;    xstr = strtrim(string(stats.mean_rej,format='(F12.3)'),2)+$
;      '\pm'+strtrim(string(stats.sigma_rej,format='(F12.3)'),2)+$
;      ' ('+strtrim(string(stats.median_rej,format='(F12.3)'),2)+')'
;
;    xtitle = 'B_{synth}'
;    ytitle = 'B_{phot}'
;
;    xrange = [min(x)<min(y),max(x)>max(y)]
;    yrange = xrange
;    
;    djs_plot, [0], [0], /nodata, xsty=3, ysty=3, xrange=xrange, yrange=yrange, $
;      xtitle=xtitle, ytitle=ytitle, charsize=2.0, charthick=postthick, $
;      xthick=postthick, ythick=postthick, position=pos[*,0]
;;   oplot, x, y, ps=8
;    oploterror, x, y, yerr, ps=3, errthick=postthick2
;    djs_oplot, !x.crange, !y.crange, line=0, thick=postthick
;    legend, textoidl('\Delta(B) = '+xstr+' mag'), /left, /top, box=0, $
;      charsize=1.6, charthick=postthick
;
;    im_openclose, postscript=postscript, /close
;
;; ---------------------------------------------------------------------------    
;
;    psname = 'synth_v_vs_phot_v'
;    im_openclose, pspath+psname, postscript=postscript, xsize=8.5, ysize=8.0, encapsulated=encapsulated
;
;    pagemaker, nx=1, ny=1, xspace=0, yspace=0.0, width=6.6, height=6.6, $
;      xmargin=[1.5,0.4], ymargin=[0.3,1.1], xpage=8.5, ypage=8.0, $
;      position=pos, /normal
;
;    indx = where((ediscs.synth_mobs_ab[1] gt 0.0) and (ediscs.kcorr_mobs_ab_aper[1] gt 0.0),nindx)
;    x = ediscs[indx].synth_mobs_ab[1]
;    y = ediscs[indx].kcorr_mobs_ab_aper[1]
;    yerr = ediscs[indx].kcorr_mobs_ab_aper_err[1]
;    
;    stats = im_stats(y-x,sigrej=3.0)
;    xstr = strtrim(string(stats.mean_rej,format='(F12.3)'),2)+$
;      '\pm'+strtrim(string(stats.sigma_rej,format='(F12.3)'),2)+$
;      ' ('+strtrim(string(stats.median_rej,format='(F12.3)'),2)+')'
;
;    xtitle = 'V_{synth}'
;    ytitle = 'V_{phot}'
;
;    xrange = [min(x)<min(y),max(x)>max(y)]
;    yrange = xrange
;    
;    djs_plot, [0], [0], /nodata, xsty=3, ysty=3, xrange=xrange, yrange=yrange, $
;      xtitle=xtitle, ytitle=ytitle, charsize=2.0, charthick=postthick, $
;      xthick=postthick, ythick=postthick, position=pos[*,0]
;;   oplot, x, y, ps=8
;    oploterror, x, y, yerr, ps=3, errthick=postthick2
;    djs_oplot, !x.crange, !y.crange, line=0, thick=postthick
;    legend, textoidl('\Delta(V) = '+xstr+' mag'), /left, /top, box=0, $
;      charsize=1.6, charthick=postthick
;
;    im_openclose, postscript=postscript, /close
;
;; ---------------------------------------------------------------------------    
;
;    psname = 'synth_i_vs_phot_i'
;    im_openclose, pspath+psname, postscript=postscript, xsize=8.5, ysize=8.0, encapsulated=encapsulated
;
;    pagemaker, nx=1, ny=1, xspace=0, yspace=0.0, width=6.6, height=6.6, $
;      xmargin=[1.5,0.4], ymargin=[0.3,1.1], xpage=8.5, ypage=8.0, $
;      position=pos, /normal
;
;    indx = where((ediscs.synth_mobs_ab[3] gt 0.0) and (ediscs.kcorr_mobs_ab_aper[3] gt 0.0),nindx)
;    x = ediscs[indx].synth_mobs_ab[3]
;    y = ediscs[indx].kcorr_mobs_ab_aper[3]
;    yerr = ediscs[indx].kcorr_mobs_ab_aper_err[3]
;    
;    stats = im_stats(y-x,sigrej=3.0)
;    xstr = strtrim(string(stats.mean_rej,format='(F12.3)'),2)+$
;      '\pm'+strtrim(string(stats.sigma_rej,format='(F12.3)'),2)+$
;      ' ('+strtrim(string(stats.median_rej,format='(F12.3)'),2)+')'
;
;    xtitle = 'I_{synth}'
;    ytitle = 'I_{phot}'
;
;    xrange = [min(x)<min(y),max(x)>max(y)]
;    yrange = xrange
;    
;    djs_plot, [0], [0], /nodata, xsty=3, ysty=3, xrange=xrange, yrange=yrange, $
;      xtitle=xtitle, ytitle=ytitle, charsize=2.0, charthick=postthick, $
;      xthick=postthick, ythick=postthick, position=pos[*,0]
;;   oplot, x, y, ps=8
;    oploterror, x, y, yerr, ps=3, errthick=postthick2
;    djs_oplot, !x.crange, !y.crange, line=0, thick=postthick
;    legend, textoidl('\Delta(I) = '+xstr+' mag'), /left, /top, box=0, $
;      charsize=1.6, charthick=postthick
;
;    im_openclose, postscript=postscript, /close
;
;; ---------------------------------------------------------------------------    
;
;    psname = 'synth_r_vs_phot_r'
;    im_openclose, pspath+psname, postscript=postscript, xsize=8.5, ysize=8.0, encapsulated=encapsulated
;
;    pagemaker, nx=1, ny=1, xspace=0, yspace=0.0, width=6.6, height=6.6, $
;      xmargin=[1.5,0.4], ymargin=[0.3,1.1], xpage=8.5, ypage=8.0, $
;      position=pos, /normal
;
;    indx = where((ediscs.synth_mobs_ab[2] gt 0.0) and (ediscs.kcorr_mobs_ab_aper[2] gt 0.0),nindx)
;    x = ediscs[indx].synth_mobs_ab[2]
;    y = ediscs[indx].kcorr_mobs_ab_aper[2]
;    yerr = ediscs[indx].kcorr_mobs_ab_aper_err[2]
;    
;    stats = im_stats(y-x,sigrej=3.0)
;    xstr = strtrim(string(stats.mean_rej,format='(F12.3)'),2)+$
;      '\pm'+strtrim(string(stats.sigma_rej,format='(F12.3)'),2)+$
;      ' ('+strtrim(string(stats.median_rej,format='(F12.3)'),2)+')'
;
;    xtitle = 'R_{synth}'
;    ytitle = 'R_{phot}'
;
;    xrange = [min(x)<min(y),max(x)>max(y)]
;    yrange = xrange
;    
;    djs_plot, [0], [0], /nodata, xsty=3, ysty=3, xrange=xrange, yrange=yrange, $
;      xtitle=xtitle, ytitle=ytitle, charsize=2.0, charthick=postthick, $
;      xthick=postthick, ythick=postthick, position=pos[*,0]
;;   oplot, x, y, ps=8
;    oploterror, x, y, yerr, ps=3, errthick=postthick2
;    djs_oplot, !x.crange, !y.crange, line=0, thick=postthick
;    legend, textoidl('\Delta(R) = '+xstr+' mag'), /left, /top, box=0, $
;      charsize=1.6, charthick=postthick
;
;    im_openclose, postscript=postscript, /close

;; ---------------------------------------------------------------------------    

    if keyword_set(postscript) and keyword_set(encapsulated) then im_ps2html, htmlbase, $
      html_path=html_path, cleanpng=0, npscols=3, psfind=0

return
end
    
