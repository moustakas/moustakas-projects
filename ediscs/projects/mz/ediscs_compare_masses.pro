pro ediscs_compare_masses, ediscs, isedfit, isedfit_maxold, bamford, postscript=postscript, $
  encapsulated=encapsulated, paper=paper, cleanpng=cleanpng, _extra=extra
; jm07sep03yu

    htmlbase = 'compare_masses'
    html_path = ediscs_path(/html)
    pspath = html_path+htmlbase+'/'
    
    analysis_path = ediscs_path(/analysis)
    isedfitpath = ediscs_path(/isedfit)
    
    if (n_elements(ediscs) eq 0L) then ediscs = read_ediscs()
    if (n_elements(isedfit) eq 0L) then isedfit = mrdfits(isedfitpath+'ediscs_isedfit.fits.gz',1,/silent)
    if (n_elements(isedfit_maxold) eq 0L) then isedfit_maxold = mrdfits(isedfitpath+'ediscs_isedfit_maxold.fits.gz',1,/silent)
    if (n_elements(bamford) eq 0L) then bamford = mrdfits(analysis_path+'Mstar_selected.fits.gz',1,/silent)

    if keyword_set(paper) then postscript = 1L

    if keyword_set(postscript) then begin
       postthick1 = 5.0 
    endif else begin
       im_window, 0, /square
       postthick1 = 2.0
    endelse

; clean up old PS and PNG files

    if keyword_set(cleanpng) then begin
       splog, 'Deleting all PNG and PS files in '+pspath
       spawn, ['rm -f '+pspath+'*.png*'], /sh
       spawn, ['rm -f '+pspath+'*.ps'], /sh
    endif

; ------------------------------------------------------------
; inter-compare various stellar mass estimates
; ------------------------------------------------------------

    psname = 'ediscs_compare_masses'
    im_openclose, pspath+psname, postscript=postscript, encapsulated=encapsulated, xsize=11.0, ysize=11.0

    pagemaker, nx=3, ny=3, position=pos, /normal, xspace=0.0, $
      xmargin=[1.4,0.2], ymargin=[0.2,1.1], yspace=0.0, xpage=11.0, ypage=11.0

    xrange = [6.4,14.0]
    yrange = xrange
    charsize1 = 1.3
    plotsym, 0, 0.7, /fill

    indx = where((ediscs.kcorr_mass gt -900.0) and (bamford.mstar_k_br gt -900.0) and $
      (isedfit.mass gt -900.0) and ((isedfit_maxold.mass gt 0.0) or $
      (isedfit_maxold.maxold_mass gt 0.0)),nindx)

    mass_kcorrect = alog10(ediscs[indx].kcorr_mass)
    mass_bamford = bamford[indx].mstar_k_br[0]
    mass_isedfit = alog10(isedfit[indx].mass)

    mass_isedfit_maxold = fltarr(nindx)
    indx1 = where((isedfit_maxold[indx].mass gt 0.0) and (isedfit_maxold[indx].maxold_mass le 0.0),nindx1)
    if (nindx1 ne 0L) then mass_isedfit_maxold[indx1] = alog10(isedfit_maxold[indx[indx1]].mass)
    indx2 = where((isedfit_maxold[indx].mass le 0.0) and (isedfit_maxold[indx].maxold_mass gt 0.0),nindx2)
    if (nindx2 ne 0L) then mass_isedfit_maxold[indx2] = alog10(isedfit_maxold[indx[indx2]].maxold_mass)
    indx3 = where((isedfit_maxold[indx].mass gt 0.0) and (isedfit_maxold[indx].maxold_mass gt 0.0),nindx3)
    if (nindx3 ne 0L) then mass_isedfit_maxold[indx3] = alog10(isedfit_maxold[indx[indx3]].maxold_mass+isedfit_maxold[indx[indx3]].mass)
    
; --------------------------------------------------    
; Panel 0
; --------------------------------------------------    

    xtitle = 'log (M/M'+sunsymbol()+') [K-correct]'
    ytitle = 'log (M/M'+sunsymbol()+') [(B-R) & K]'

    stats = im_stats(mass_bamford-mass_kcorrect,verbose=0,sigrej=3.0)
    xstr = '\Delta = '+strtrim(string(stats.mean_rej,format='(F12.3)'),2)+$
      '\pm'+strtrim(string(stats.sigma_rej,format='(F12.3)'),2)+' dex'

    djs_plot, [0], [0], /nodata, xsty=3, ysty=3, xthick=postthick1, ythick=postthick1, $
      charthick=postthick1, charsize=charsize1, xrange=xrange, yrange=yrange, $
      xtitle='', ytitle=ytitle, position=pos[*,0], xtickname=replicate(' ',10)
    djs_oplot, !x.crange, !y.crange, thick=postthick1
    djs_oplot, mass_kcorrect, mass_bamford, ps=8
 
    legend, textoidl(xstr), /right, /bottom, box=0, charsize=charsize1, charthick=postthick1
    legend, '(a)', /left, /top, box=0, charsize=charsize1, charthick=postthick1
    
; --------------------------------------------------    
; Panel 1
; --------------------------------------------------    

    xtitle = 'log (M/M'+sunsymbol()+') [isedfit-maxold]'
    ytitle = 'log (M/M'+sunsymbol()+') [(B-R) & K]'

    stats = im_stats(mass_bamford-mass_isedfit_maxold,verbose=0,sigrej=3.0)
    xstr = '\Delta = '+strtrim(string(stats.mean_rej,format='(F12.3)'),2)+$
      '\pm'+strtrim(string(stats.sigma_rej,format='(F12.3)'),2)+' dex'

    djs_plot, [0], [0], /nodata, /noerase, xsty=3, ysty=3, xthick=postthick1, ythick=postthick1, $
      charthick=postthick1, charsize=charsize1, xrange=xrange, yrange=yrange, $
      xtitle='', ytitle='', position=pos[*,1], xtickname=replicate(' ',10), ytickname=replicate(' ',10)
    djs_oplot, !x.crange, !y.crange, thick=postthick1
    djs_oplot, mass_isedfit_maxold, mass_bamford, ps=8
 
    legend, textoidl(xstr), /right, /bottom, box=0, charsize=charsize1, charthick=postthick1
    legend, '(b)', /left, /top, box=0, charsize=charsize1, charthick=postthick1
    
; --------------------------------------------------    
; Panel 2
; --------------------------------------------------    

    xtitle = 'log (M/M'+sunsymbol()+') [isedfit]'
    ytitle = 'log (M/M'+sunsymbol()+') [(B-R) & K]'

    stats = im_stats(mass_bamford-mass_isedfit,verbose=0,sigrej=3.0)
    xstr = '\Delta = '+strtrim(string(stats.mean_rej,format='(F12.3)'),2)+$
      '\pm'+strtrim(string(stats.sigma_rej,format='(F12.3)'),2)+' dex'

    djs_plot, [0], [0], /nodata, /noerase, xsty=3, ysty=3, xthick=postthick1, ythick=postthick1, $
      charthick=postthick1, charsize=charsize1, xrange=xrange, yrange=yrange, $
      xtitle=xtitle, ytitle='', position=pos[*,2], ytickname=replicate(' ',10)
    djs_oplot, !x.crange, !y.crange, thick=postthick1
    djs_oplot, mass_isedfit, mass_bamford, ps=8
 
    legend, textoidl(xstr), /right, /bottom, box=0, charsize=charsize1, charthick=postthick1
    legend, '(c)', /left, /top, box=0, charsize=charsize1, charthick=postthick1

; --------------------------------------------------    
; Panel 3
; --------------------------------------------------    

    xtitle = 'log (M/M'+sunsymbol()+') [K-correct]'
    ytitle = 'log (M/M'+sunsymbol()+') [isedfit]'

    stats = im_stats(mass_isedfit-mass_kcorrect,verbose=0,sigrej=3.0)
    xstr = '\Delta = '+strtrim(string(stats.mean_rej,format='(F12.3)'),2)+$
      '\pm'+strtrim(string(stats.sigma_rej,format='(F12.3)'),2)+' dex'

    djs_plot, [0], [0], /nodata, /noerase, xsty=3, ysty=3, xthick=postthick1, ythick=postthick1, $
      charthick=postthick1, charsize=charsize1, xrange=xrange, yrange=yrange, $
      xtitle='', ytitle=ytitle, position=pos[*,3], xtickname=replicate(' ',10)
    djs_oplot, !x.crange, !y.crange, thick=postthick1
    djs_oplot, mass_kcorrect, mass_isedfit, ps=8
 
    legend, textoidl(xstr), /right, /bottom, box=0, charsize=charsize1, charthick=postthick1
    legend, '(d)', /left, /top, box=0, charsize=charsize1, charthick=postthick1

; --------------------------------------------------    
; Panel 4
; --------------------------------------------------    
    
    xtitle = 'log (M/M'+sunsymbol()+') [isedfit-maxold]'
    ytitle = 'log (M/M'+sunsymbol()+') [isedfit]'

    stats = im_stats(mass_isedfit-mass_isedfit_maxold,verbose=0,sigrej=3.0)
    xstr = '\Delta = '+strtrim(string(stats.mean_rej,format='(F12.3)'),2)+$
      '\pm'+strtrim(string(stats.sigma_rej,format='(F12.3)'),2)+' dex'

    djs_plot, [0], [0], /nodata, /noerase, xsty=3, ysty=3, xthick=postthick1, ythick=postthick1, $
      charthick=postthick1, charsize=charsize1, xrange=xrange, yrange=yrange, $
      xtitle=xtitle, ytitle='', position=pos[*,4], ytickname=replicate(' ',10)
    djs_oplot, !x.crange, !y.crange, thick=postthick1
    djs_oplot, mass_isedfit_maxold, mass_isedfit, ps=8
 
    legend, textoidl(xstr), /right, /bottom, box=0, charsize=charsize1, charthick=postthick1
    legend, '(e)', /left, /top, box=0, charsize=charsize1, charthick=postthick1

; --------------------------------------------------    
; Panel 6
; --------------------------------------------------    

    xtitle = 'log (M/M'+sunsymbol()+') [K-correct]'
    ytitle = 'log (M/M'+sunsymbol()+') [isedfit-maxold]'

    stats = im_stats(mass_isedfit_maxold-mass_kcorrect,verbose=0,sigrej=3.0)
    xstr = '\Delta = '+strtrim(string(stats.mean_rej,format='(F12.3)'),2)+$
      '\pm'+strtrim(string(stats.sigma_rej,format='(F12.3)'),2)+' dex'

    djs_plot, [0], [0], /nodata, /noerase, xsty=3, ysty=3, xthick=postthick1, ythick=postthick1, $
      charthick=postthick1, charsize=charsize1, xrange=xrange, yrange=yrange, $
      xtitle=xtitle, ytitle=ytitle, position=pos[*,6]
    djs_oplot, !x.crange, !y.crange, thick=postthick1
    djs_oplot, mass_kcorrect, mass_isedfit_maxold, ps=8
 
    legend, textoidl(xstr), /right, /bottom, box=0, charsize=charsize1, charthick=postthick1
    legend, '(f)', /left, /top, box=0, charsize=charsize1, charthick=postthick1

    im_openclose, postscript=postscript, /close    

; --------------------------------------------------    
; GENERATE THE WEB PAGE
; --------------------------------------------------    
    
    if keyword_set(postscript) and keyword_set(encapsulated) then $
      im_ps2html, htmlbase, html_path=html_path, cleanpng=0, _extra=extra

; --------------------------------------------------    
; SELECT PLOTS FOR THE PAPER HERE
; --------------------------------------------------    

    if keyword_set(paper) then begin
;      paperpath = ediscs_path(/papers)+'FIG_PAPER1/'
;      splog, 'Writing paper plots to '+paperpath+'.'
;      paperplots = [$
;        ]+'.*ps'
;      for k = 0L, n_elements(paperplots)-1L do $
;        spawn, ['/bin/cp -f '+html_path+htmlbase+'/'+paperplots[k]+' '+paperpath], /sh
    endif
    
    
stop    
    
return
end
    
