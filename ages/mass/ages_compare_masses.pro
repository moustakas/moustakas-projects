;+
; NAME:
;       AGES_COMPARE_MASSES
;
; PURPOSE:
;       Compare various stellar mass estimates.
;
; CALLING SEQUENCE:
;
;
; INPUTS:
;
;
; OPTIONAL INPUTS:
;
;
; KEYWORD PARAMETERS:
;
;
; OUTPUTS:
;
;
; OPTIONAL OUTPUTS:
;
;
; PROCEDURES USED:
;
;
; COMMENTS:
;
;
; EXAMPLES:
;
;
; MODIFICATION HISTORY:
;       J. Moustakas, 2005 Mar 01, U of A - written
;-

pro ages_compare_masses, ages, BwRIK, BwRI, postscript=postscript, paper=paper, $
  cleanpng=cleanpng, _extra=extra

    htmlbase = 'compare_masses'

    html_path = ages_path(/web)
    pspath = html_path+htmlbase+'/'
    
    datapath = ages_path(/mass)
    paperpath = ages_path(/papers)+'FIG_PAPER1/'
    
    if (n_elements(ages) eq 0L) then ages = read_ages()
    if (n_elements(BwRIK) eq 0L) then BwRIK = mrdfits(datapath+'ages_BwRIK_isedfit_mass.fits.gz',1,/silent)
    if (n_elements(BwRI) eq 0L) then BwRI = mrdfits(datapath+'ages_BwRI_isedfit_mass.fits.gz',1,/silent)

    if keyword_set(paper) then postscript = 1L

    if keyword_set(postscript) then begin
       postthick = 8.0 
    endif else begin
       im_window, 0, /square
       postthick = 2.0
    endelse

; clean up old PS and PNG files

    if keyword_set(cleanpng) then begin
       splog, 'Deleting all PNG and PS files in '+pspath
       spawn, ['rm -f '+pspath+'*.png*'], /sh
       spawn, ['rm -f '+pspath+'*.ps'], /sh
    endif

; --------------------------------------------------    
; M(isedfit - BwRIK) vs M(isedfit - BwRI)
; --------------------------------------------------    

; what effect does it have on the stellar mass estimates to take away
; the K-band photometry?
    
    psname = 'mass_isedfit_BwRIK_vs_mass_isedfit_BwRI'
    im_openclose, pspath+psname, postscript=postscript, /encapsulated
    
    indx = where((BwRIK.isedfit_dof eq 3.0) and (BwRIK.isedfit_chi2min lt 2.0),nindx)
    ngalaxy = n_elements(BwRIK)
    
    x = BwRIK[indx].isedfit_mass
    xerr = BwRIK[indx].isedfit_mass_err
    
    y = BwRI[indx].isedfit_mass
    yerr = BwRI[indx].isedfit_mass_err
    
    residuals = x - y
    residuals_err = sqrt(xerr^2 + yerr^2)

    stats = im_stats(residuals,/verbose)
    xstr = strtrim(string(stats.median,format='(F12.2)'),2)+'\pm'+strtrim(string(stats.sigma_rej,format='(F12.2)'),2)

    xtitle = 'log (M / M'+sunsymbol()+') [isedfit - BwRIK]'
    ytitle = 'log (M / M'+sunsymbol()+') [isedfit - BwRI]'
    ytitle2 = 'Residuals'
    
    xrange = [min(x)<min(y),max(x)>max(y)]*[0.95,1.05]
    yrange = xrange

    xrange2 = xrange
    inrange = where((x gt xrange2[0]) and (x lt xrange[1]))
    yrange2 = max(abs(residuals[inrange]))*[-1.1,1.1]

    pagemaker, nx=1, ny=2, yspace=0, xspace=0, xmargin=[1.4,0.3], $
      ymargin=[0.3,1.3], position=pos, height=[6.0,3.4], /normal

    ages_lineplot, x, y, xerr, yerr, plottype=1, postscript=postscript, $
      xtitle='', ytitle=ytitle, xrange=xrange, yrange=yrange, legendtype=0, $
      /right, /top, position=pos[*,0], yminor=3, charsize=1.8, xtickname=replicate(' ',10)
    djs_oplot, !x.crange, !y.crange, line=0, thick=postthick
    legend, string(nindx,format='(I0)')+'/'+string(ngalaxy,format='(I0)'), /right, $
      /bottom, box=0, charsize=1.8, charthick=postthick
    legend, textoidl('\chi^{2}_{\nu} < 2'), /left, /top, $
      box=0, charsize=1.8, charthick=postthick
    
    ages_lineplot, x, residuals, xerr, residuals_err, plottype=1, $
      postscript=postscript, xtitle=xtitle, ytitle=ytitle2, xrange=xrange2, $
      yrange=yrange2, legendtype=0, /right, /top, position=pos[*,1], /noerase, $
      yminor=3, charsize=1.8
    djs_oplot, !x.crange, [0,0], line=0, thick=postthick

    legend, textoidl(xstr), /right, /bottom, box=0, charsize=1.8, charthick=postthick

    im_openclose, postscript=postscript, /close    

; --------------------------------------------------    
; M(isedfit) vs M(Bell - (g-r)_synth & r)
; --------------------------------------------------    

    psname = 'mass_isedfit_vs_mass_bell_gr_synth_r'
    im_openclose, pspath+psname, postscript=postscript, /encapsulated
    
    indx = where((ages.isedfit_mass gt -900.0) and (ages.mass_synth_gr_r gt -900) and $
      (ages.isedfit_chi2min lt 2.0),nindx)
    ngalaxy = n_elements(ages)
    
    x = ages[indx].isedfit_mass
    xerr = ages[indx].isedfit_mass_err
    
    y = ages[indx].mass_synth_gr_r
    yerr = ages[indx].mass_synth_gr_r_err
    
    residuals = x - y
    residuals_err = sqrt(xerr^2 + yerr^2)

    stats = im_stats(residuals,/verbose)
    xstr = strtrim(string(stats.median,format='(F12.2)'),2)+'\pm'+strtrim(string(stats.sigma_rej,format='(F12.2)'),2)

    xtitle = 'log (M / M'+sunsymbol()+') [isedfit]'
    ytitle = 'log (M / M'+sunsymbol()+') [Bell - (g-r)_{synth} & r]'
    ytitle2 = 'Residuals'
    
    xrange = [min(x)<min(y),max(x)>max(y)]*[0.95,1.05]
    yrange = xrange

    xrange2 = xrange
    inrange = where((x gt xrange2[0]) and (x lt xrange[1]))
    yrange2 = max(abs(residuals[inrange]))*[-1.1,1.1]

    pagemaker, nx=1, ny=2, yspace=0, xspace=0, xmargin=[1.4,0.3], $
      ymargin=[0.3,1.3], position=pos, height=[6.0,3.4], /normal

    ages_lineplot, x, y, xerr, yerr, plottype=1, postscript=postscript, $
      xtitle='', ytitle=ytitle, xrange=xrange, yrange=yrange, legendtype=0, $
      /right, /top, position=pos[*,0], yminor=3, charsize=1.8, xtickname=replicate(' ',10)
    djs_oplot, !x.crange, !y.crange, line=0, thick=postthick
    legend, string(nindx,format='(I0)')+'/'+string(ngalaxy,format='(I0)'), /right, $
      /bottom, box=0, charsize=1.8, charthick=postthick
    legend, textoidl('\chi^{2}_{\nu} < 2'), /left, /top, $
      box=0, charsize=1.8, charthick=postthick
    
    ages_lineplot, x, residuals, xerr, residuals_err, plottype=1, $
      postscript=postscript, xtitle=xtitle, ytitle=ytitle2, xrange=xrange2, $
      yrange=yrange2, legendtype=0, /right, /top, position=pos[*,1], /noerase, $
      yminor=3, charsize=1.8
    djs_oplot, !x.crange, [0,0], line=0, thick=postthick

    legend, textoidl(xstr), /right, /bottom, box=0, charsize=1.8, charthick=postthick

    im_openclose, postscript=postscript, /close    

; --------------------------------------------------    
; M(isedfit) vs M(Bell - g-r & r)
; --------------------------------------------------    

    psname = 'mass_isedfit_vs_mass_bell_gr_r'
    im_openclose, pspath+psname, postscript=postscript, /encapsulated
    
    indx = where((ages.isedfit_mass gt -900.0) and (ages.mass_gr_r gt -900) and $
      (ages.isedfit_chi2min lt 2.0),nindx)
    ngalaxy = n_elements(ages)
    
    x = ages[indx].isedfit_mass
    xerr = ages[indx].isedfit_mass_err
    
    y = ages[indx].mass_gr_r
    yerr = ages[indx].mass_gr_r_err
    
    residuals = x - y
    residuals_err = sqrt(xerr^2 + yerr^2)

    stats = im_stats(residuals,/verbose)
    xstr = strtrim(string(stats.median,format='(F12.2)'),2)+'\pm'+strtrim(string(stats.sigma_rej,format='(F12.2)'),2)

    xtitle = 'log (M / M'+sunsymbol()+') [isedfit]'
    ytitle = 'log (M / M'+sunsymbol()+') [Bell - (g-r) & r]'
    ytitle2 = 'Residuals'
    
    xrange = [min(x)<min(y),max(x)>max(y)]*[0.95,1.05]
    yrange = xrange

    xrange2 = xrange
    inrange = where((x gt xrange2[0]) and (x lt xrange[1]))
    yrange2 = max(abs(residuals[inrange]))*[-1.1,1.1]

    pagemaker, nx=1, ny=2, yspace=0, xspace=0, xmargin=[1.4,0.3], $
      ymargin=[0.3,1.3], position=pos, height=[6.0,3.4], /normal

    ages_lineplot, x, y, xerr, yerr, plottype=1, postscript=postscript, $
      xtitle='', ytitle=ytitle, xrange=xrange, yrange=yrange, legendtype=0, $
      /right, /top, position=pos[*,0], yminor=3, charsize=1.8, xtickname=replicate(' ',10)
    djs_oplot, !x.crange, !y.crange, line=0, thick=postthick
    legend, string(nindx,format='(I0)')+'/'+string(ngalaxy,format='(I0)'), /right, $
      /bottom, box=0, charsize=1.8, charthick=postthick
    legend, textoidl('\chi^{2}_{\nu} < 2'), /left, /top, $
      box=0, charsize=1.8, charthick=postthick
    
    ages_lineplot, x, residuals, xerr, residuals_err, plottype=1, $
      postscript=postscript, xtitle=xtitle, ytitle=ytitle2, xrange=xrange2, $
      yrange=yrange2, legendtype=0, /right, /top, position=pos[*,1], /noerase, $
      yminor=3, charsize=1.8
    djs_oplot, !x.crange, [0,0], line=0, thick=postthick

    legend, textoidl(xstr), /right, /bottom, box=0, charsize=1.8, charthick=postthick

    im_openclose, postscript=postscript, /close    

; --------------------------------------------------    
; M(isedfit) vs M(Bell - V-H & H)
; --------------------------------------------------    

    psname = 'mass_isedfit_vs_mass_bell_VHH'
    im_openclose, pspath+psname, postscript=postscript, /encapsulated
    
    indx = where((ages.isedfit_mass gt -900.0) and (ages.mass_vh_h gt -900) and $
      (ages.isedfit_chi2min lt 2.0),nindx)
    ngalaxy = n_elements(ages)
    
    x = ages[indx].isedfit_mass
    xerr = ages[indx].isedfit_mass_err
    
    y = ages[indx].mass_vh_h
    yerr = ages[indx].mass_vh_h_err
    
    residuals = x - y
    residuals_err = sqrt(xerr^2 + yerr^2)

    stats = im_stats(residuals,/verbose)
    xstr = strtrim(string(stats.median,format='(F12.2)'),2)+'\pm'+strtrim(string(stats.sigma_rej,format='(F12.2)'),2)

    xtitle = 'log (M / M'+sunsymbol()+') [isedfit]'
    ytitle = 'log (M / M'+sunsymbol()+') [Bell - (V-H) & H]'
    ytitle2 = 'Residuals'
    
    xrange = [min(x)<min(y),max(x)>max(y)]*[0.95,1.05]
    yrange = xrange

    xrange2 = xrange
    inrange = where((x gt xrange2[0]) and (x lt xrange[1]))
    yrange2 = max(abs(residuals[inrange]))*[-1.1,1.1]

    pagemaker, nx=1, ny=2, yspace=0, xspace=0, xmargin=[1.4,0.3], $
      ymargin=[0.3,1.3], position=pos, height=[6.0,3.4], /normal

    ages_lineplot, x, y, xerr, yerr, plottype=1, postscript=postscript, $
      xtitle='', ytitle=ytitle, xrange=xrange, yrange=yrange, legendtype=0, $
      /right, /top, position=pos[*,0], yminor=3, charsize=1.8, xtickname=replicate(' ',10)
    djs_oplot, !x.crange, !y.crange, line=0, thick=postthick
    legend, string(nindx,format='(I0)')+'/'+string(ngalaxy,format='(I0)'), /right, $
      /bottom, box=0, charsize=1.8, charthick=postthick
    legend, textoidl('\chi^{2}_{\nu} < 2'), /left, /top, $
      box=0, charsize=1.8, charthick=postthick
    
    ages_lineplot, x, residuals, xerr, residuals_err, plottype=1, $
      postscript=postscript, xtitle=xtitle, ytitle=ytitle2, xrange=xrange2, $
      yrange=yrange2, legendtype=0, /right, /top, position=pos[*,1], /noerase, $
      yminor=3, charsize=1.8
    djs_oplot, !x.crange, [0,0], line=0, thick=postthick

    legend, textoidl(xstr), /right, /bottom, box=0, charsize=1.8, charthick=postthick

    im_openclose, postscript=postscript, /close    

; --------------------------------------------------    
; M_err (Monte Carlo) vs M_err (Chi2)
; --------------------------------------------------    

    psname = 'mass_err_MC_vs_mass_err_dchi2'
    im_openclose, pspath+psname, postscript=postscript, /encapsulated
    
    indx = where((ages.isedfit_mass_err gt -900.0) and (ages.isedfit_dchi2_mass_err gt -900) and $
      (ages.isedfit_mass gt -900) and (ages.isedfit_chi2min lt 2.0),nindx)
    ngalaxy = n_elements(ages)
    
    x = ages[indx].isedfit_mass

    y1 = ages[indx].isedfit_mass_err
    y2 = ages[indx].isedfit_dchi2_mass_err
    y = y1-y2
    
    stats = im_stats(y,/verbose)
    xstr = strtrim(string(stats.median,format='(F12.2)'),2)+'\pm'+strtrim(string(stats.sigma_rej,format='(F12.2)'),2)

    xtitle = 'log (M / M'+sunsymbol()+') [isedfit]'
    ytitle = '\sigma(M) [Monte-Carlo] - \sigma(M) [\Delta\chi^{2}]'
    
    xrange = minmax(x)*[0.95,1.05]
    yrange = max(abs(y))*[-1.05,1.05]

    ages_lineplot, x, y, x*0.0, y*0.0, plottype=1, postscript=postscript, $
      xtitle=xtitle, ytitle=ytitle, xrange=xrange, yrange=yrange, legendtype=0, $
      /right, /top, yminor=3, charsize=1.8
    djs_oplot, !x.crange, [0,0], line=0, thick=postthick

    legend, string(nindx,format='(I0)')+'/'+string(ngalaxy,format='(I0)'), /right, $
      /bottom, box=0, charsize=1.8, charthick=postthick
    legend, textoidl('\chi^{2}_{\nu} < 2'), /left, /top, $
      box=0, charsize=1.8, charthick=postthick
    legend, textoidl(xstr), /right, /top, box=0, charsize=1.8, charthick=postthick

    im_openclose, postscript=postscript, /close    

; ------------------------------------------------------------
; inter-compare various stellar mass estimates
; ------------------------------------------------------------

    psname = 'ages_stellar_mass'
    im_openclose, pspath+psname, postscript=postscript, /encapsulated

    pagemaker, nx=2, ny=2, position=pos, /normal, xspace=0.0, $
      xmargin=[1.4,0.2], ymargin=[0.2,1.1], yspace=0.0

    xrange = [7,14]
    yrange = xrange
    multicharsize = 1.5

; --------------------------------------------------    
; Panel 1 - M[(V-H) & H] vs M[(g-r) & r]
; --------------------------------------------------    

    indx = where((ages.mass_VH_H gt -900.0) and (ages.mass_gr_r gt -900.0),nindx)
 
    x = ages[indx].mass_VH_H
    xerr = ages[indx].mass_VH_H_err

    y = ages[indx].mass_gr_r
    yerr = ages[indx].mass_gr_r_err

    xtitle = 'log (M / M'+sunsymbol()+') [(V-H) & H]'
    ytitle = 'log (M / M'+sunsymbol()+') [(g-r) & r]'

    stats = im_stats(x-y,/verbose)
    xstr = strtrim(string(stats.median,format='(F12.2)'),2)+'\pm'+strtrim(string(stats.sigma_rej,format='(F12.2)'),2)
    
    ages_lineplot, x, y, xerr, yerr, plottype=1, postscript=postscript, $
      xtitle='', ytitle=ytitle, xrange=xrange, yrange=yrange, legendtype=0, charsize=multicharsize, $
      xstyle=3, ystyle=3, /right, /top, position=pos[*,0], xtickname=replicate(' ',10), $
      xminor=3, yminor=3
    djs_oplot, !x.crange, !y.crange, thick=postthick
 
;   legend, '('+string(nindx,format='(I0)')+')', /left, /top, box=0, charsize=multicharsize, charthick=postthick
    legend, textoidl(xstr), /right, /bottom, box=0, charsize=multicharsize, charthick=postthick
    legend, '(a)', /left, /top, box=0, charsize=multicharsize, charthick=postthick

; --------------------------------------------------    
; Panel 2 - M[L(K)] vs BdJ: M[(g-r) & r]
; --------------------------------------------------    

    indx = where((ages.K_lum gt -900.0) and (ages.mass_gr_r gt -900.0),nindx)
 
    x = ages[indx].K_lum; + alog10(0.8)
    xerr = ages[indx].K_lum_err

    y = ages[indx].mass_gr_r
    yerr = ages[indx].mass_gr_r_err

    xtitle = 'log (M / M'+sunsymbol()+') [L(K)]'
    ytitle = 'log (M / M'+sunsymbol()+') [(g-r) & r]'

    stats = im_stats(x-y,/verbose,/no_head)
    xstr = strtrim(string(stats.median,format='(F12.2)'),2)+'\pm'+strtrim(string(stats.sigma_rej,format='(F12.2)'),2)
    
    ages_lineplot, x, y, xerr, yerr, plottype=1, postscript=postscript, $
      xtitle=xtitle, ytitle='', xrange=xrange, yrange=yrange, legendtype=0, charsize=multicharsize, $
      xstyle=3, ystyle=3, /right, /top, position=pos[*,1], /noerase, ytickname=replicate(' ',10), $
      xminor=3, yminor=3
    djs_oplot, !x.crange, !y.crange, thick=postthick

;   legend, '('+string(nindx,format='(I0)')+')', /left, /top, box=0, charsize=multicharsize, charthick=postthick
    legend, textoidl(xstr), /right, /bottom, box=0, charsize=multicharsize, charthick=postthick
    legend, '(b)', /left, /top, box=0, charsize=multicharsize, charthick=postthick

; --------------------------------------------------    
; Panel 3 - M[BdJ: (V-H) & H] vs M[L(K)]
; --------------------------------------------------    

    indx = where((ages.mass_VH_H gt -900.0) and (ages.K_lum gt -900.0),nindx)
 
    x = ages[indx].mass_VH_H
    xerr = ages[indx].mass_VH_H_err

    y = ages[indx].K_lum; + alog10(0.8)
    yerr = ages[indx].K_lum_err

    xtitle = 'log (M / M'+sunsymbol()+') [(V-H) & H]'
    ytitle = 'log (M / M'+sunsymbol()+') [L(K)]'

    stats = im_stats(x-y,/verbose,/no_head)
    xstr = strtrim(string(stats.median,format='(F12.2)'),2)+'\pm'+strtrim(string(stats.sigma_rej,format='(F12.2)'),2)
    
    ages_lineplot, x, y, xerr, yerr, plottype=1, postscript=postscript, $
      xtitle=xtitle, ytitle=ytitle, xrange=xrange, yrange=yrange, legendtype=0, $
      xstyle=3, ystyle=3, /right, /top, position=pos[*,2], /noerase, charsize=multicharsize, $
      xminor=3, yminor=3
    djs_oplot, !x.crange, !y.crange, thick=postthick

;   legend, '('+string(nindx,format='(I0)')+')', /left, /top, box=0, charsize=multicharsize, charthick=postthick
    legend, textoidl(xstr), /right, /bottom, box=0, charsize=multicharsize, charthick=postthick
    legend, '(c)', /left, /top, box=0, charsize=multicharsize, charthick=postthick
    
; --------------------------------------------------    
; Panel 4 - No Data
; --------------------------------------------------    

    im_openclose, postscript=postscript, /close    

; --------------------------------------------------    
; GENERATE THE WEB PAGE
; --------------------------------------------------    
    
    if keyword_set(postscript) then $
      im_ps2html, htmlbase, html_path=html_path, cleanpng=0, _extra=extra

; --------------------------------------------------    
; SELECT PLOTS FOR THE PAPER HERE
; --------------------------------------------------    

    if keyword_set(paper) then begin

       splog, 'Writing paper plots to '+paperpath+'.'
       paperplots = [$
         'mass_isedfit_BwRIK_vs_mass_isedfit_BwRI',$
         'mass_isedfit_vs_mass_bell_gr_synth_r'$         
         ]+'.*ps'

       for k = 0L, n_elements(paperplots)-1L do $
         spawn, ['/bin/cp -f '+html_path+htmlbase+'/'+paperplots[k]+' '+paperpath], /sh

    endif
    
    
stop    
    
return
end
    
