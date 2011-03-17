; ------------------------------------------------------------
; MB versus abundance gradient (dex/R25) from Pilyugin 2004 
; ------------------------------------------------------------

    psname = 'MB_vs_gradient_r25'
    im_openclose, pspath+psname, postscript=postscript, xsize=8.5, ysize=7.5

    pagemaker, nx=1, ny=1, xspace=0, yspace=0.0, width=6.6, height=6.6, $
      xmargin=[1.5,0.4], ymargin=[0.3,1.1], xpage=8.5, ypage=8.0, $
      position=pos, /normal

    if keyword_set(dotalk) then polyfill, [1,1,0,0,1], [1,0,0,1,1], /normal, $
      color=djs_icolor('black')

    p04 = read_04pilyugin()
    indx = where((p04.m_b gt -900) and (p04.ohgradient_r25 gt -900) and $
      (p04.ohgradient_good),nindx)

    xrange = [-16,-24] ; mbrange
    yrange = gradientrange_r25

    xtitle = 'M_{B}'
    ytitle = 'Oxygen Gradient [dex/R_{25}]'
    
    plotsym, 8, 1.0, /fill
    djs_plot, [0], [0], /nodata, xtitle=xtitle, ytitle=ytitle, position=pos[*,0], $
      xrange=xrange, yrange=yrange, xstyle=3, ystyle=3, color=djs_icolor(talkcolor), $
      xthick=postthick, ythick=postthick, charsize=2.0, charthick=postthick, $
      noerase=keyword_set(dotalk)
;   djs_oplot, p04[indx].m_b, p04[indx].ohgradient_r25, ps=8, color='dark green'
    oploterror, p04[indx].m_b, p04[indx].ohgradient_r25, p04[indx].m_b_err, $
      p04[indx].logoh_rms, ps=8, color=djs_icolor('dark green'), $
      errcolor=djs_icolor('dark green'), thick=postthick, errthick=postthick

    im_openclose, postscript=postscript, /close    
    
; ------------------------------------------------------------
; mass versus abundance gradient (dex/R25) from Pilyugin 2004 
; ------------------------------------------------------------

    psname = 'mass_vs_gradient_r25'
    im_openclose, pspath+psname, postscript=postscript, xsize=8.5, ysize=7.5

    pagemaker, nx=1, ny=1, xspace=0, yspace=0.0, width=6.6, height=6.6, $
      xmargin=[1.5,0.4], ymargin=[0.3,1.1], xpage=8.5, ypage=8.0, $
      position=pos, /normal

    if keyword_set(dotalk) then polyfill, [1,1,0,0,1], [1,0,0,1,1], /normal, $
      color=djs_icolor('black')

    p04 = read_04pilyugin()
    indx = where((p04.mass_bv_b gt -900) and (p04.ohgradient_r25 gt -900) and $
      (p04.ohgradient_good),nindx)

    xrange = massrange2
    yrange = gradientrange_r25

    xtitle = 'log (M / M'+sunsymbol()+')'
    ytitle = 'Oxygen Gradient [dex/R_{25}]'
    
    plotsym, 8, 1.0, /fill
    djs_plot, [0], [0], /nodata, xtitle=xtitle, ytitle=ytitle, position=pos[*,0], $
      xrange=xrange, yrange=yrange, xstyle=3, ystyle=3, color=djs_icolor(talkcolor), $
      xthick=postthick, ythick=postthick, charsize=2.0, charthick=postthick, $
      noerase=keyword_set(dotalk)
;   djs_oplot, p04[indx].m_b, p04[indx].ohgradient_r25, ps=8, color='dark green'
    oploterror, p04[indx].mass_bv_b, p04[indx].ohgradient_r25, p04[indx].mass_bv_b_err, $
      p04[indx].logoh_rms, ps=8, color=djs_icolor('dark green'), $
      errcolor=djs_icolor('dark green'), thick=postthick, errthick=postthick

    im_openclose, postscript=postscript, /close    
    
; ------------------------------------------------------------
; MB versus abundance gradient (dex/kpc) from Pilyugin 2004 
; ------------------------------------------------------------

    psname = 'MB_vs_gradient_kpc'
    im_openclose, pspath+psname, postscript=postscript, xsize=8.5, ysize=7.5

    pagemaker, nx=1, ny=1, xspace=0, yspace=0.0, width=6.6, height=6.6, $
      xmargin=[1.5,0.4], ymargin=[0.3,1.1], xpage=8.5, ypage=8.0, $
      position=pos, /normal

    if keyword_set(dotalk) then polyfill, [1,1,0,0,1], [1,0,0,1,1], /normal, $
      color=djs_icolor('black')

    p04 = read_04pilyugin()
    indx = where((p04.m_b gt -900) and (p04.ohgradient_kpc gt -900) and $
      (p04.ohgradient_good),nindx)

    xrange = [-16,-24] ; mbrange
    yrange = gradientrange_kpc

    xtitle = 'M_{B}'
    ytitle = 'Oxygen Gradient [dex/kpc]'
    
    plotsym, 8, 1.0, /fill
    djs_plot, [0], [0], /nodata, xtitle=xtitle, ytitle=ytitle, position=pos[*,0], $
      xrange=xrange, yrange=yrange, xstyle=3, ystyle=3, color=djs_icolor(talkcolor), $
      xthick=postthick, ythick=postthick, charsize=2.0, charthick=postthick, $
      noerase=keyword_set(dotalk)
    djs_oplot, p04[indx].m_b, p04[indx].ohgradient_kpc, ps=8, color='dark green'
;   oploterror, p04[indx].m_b, p04[indx].ohgradient_kpc, p04[indx].m_b_err, p04[indx].logoh_rms, ps=8, $
;     color=djs_icolor('dark green'), errcolor=djs_icolor('dark green'), $
;     thick=postthick, errthick=postthick

stop
    
    im_openclose, postscript=postscript, /close    
    
; ------------------------------------------------------------
; mass versus abundance gradient (dex/kpc) from Pilyugin 2004 
; ------------------------------------------------------------

    psname = 'mass_vs_gradient_kpc'
    im_openclose, pspath+psname, postscript=postscript, xsize=8.5, ysize=7.5

    pagemaker, nx=1, ny=1, xspace=0, yspace=0.0, width=6.6, height=6.6, $
      xmargin=[1.5,0.4], ymargin=[0.3,1.1], xpage=8.5, ypage=8.0, $
      position=pos, /normal

    if keyword_set(dotalk) then polyfill, [1,1,0,0,1], [1,0,0,1,1], /normal, $
      color=djs_icolor('black')

    p04 = read_04pilyugin()
    indx = where((p04.mass_bv_b gt -900) and (p04.ohgradient_kpc gt -900) and $
      (p04.ohgradient_good),nindx)

    xrange = massrange2
    yrange = gradientrange_kpc

    xtitle = 'log (M / M'+sunsymbol()+')'
    ytitle = 'Oxygen Gradient [dex/kpc]'
    
    plotsym, 8, 1.0, /fill
    djs_plot, [0], [0], /nodata, xtitle=xtitle, ytitle=ytitle, position=pos[*,0], $
      xrange=xrange, yrange=yrange, xstyle=3, ystyle=3, color=djs_icolor(talkcolor), $
      xthick=postthick, ythick=postthick, charsize=2.0, charthick=postthick, $
      noerase=keyword_set(dotalk)
    djs_oplot, p04[indx].mass_bv_b, p04[indx].ohgradient_kpc, ps=8, color='dark green'
;   oploterror, p04[indx].mass_bv_b, p04[indx].ohgradient_kpc, p04[indx].mass_bv_b_err, p04[indx].logoh_rms, ps=8, $
;     color=djs_icolor('dark green'), errcolor=djs_icolor('dark green'), $
;     thick=postthick, errthick=postthick

    im_openclose, postscript=postscript, /close    
    
; ------------------------------------------------------------
; Mass vs 12 + log (O/H)
; ------------------------------------------------------------

    psname = 'mass_vs_12oh_oiiinii'
    im_openclose, pspath+psname, postscript=postscript, xsize=8.5, ysize=7.5

    pagemaker, nx=1, ny=1, xspace=0, yspace=0.0, width=5.8, height=5.5, $
      xmargin=[2.2,0.5], ymargin=[0.9,1.1], xpage=8.5, ypage=7.5, $
      position=pos, /normal

    indx = where((sdssnodust.zstrong_12oh_oiiinii_moustakas gt -900) and (sdssnodust.mass_bv_b gt -900),nindx)

    x = sdssnodust[indx].mass_bv_b
    xerr = sdssnodust[indx].mass_bv_b_err

    y = sdssnodust[indx].zstrong_12oh_oiiinii_moustakas
    yerr = sdssnodust[indx].zstrong_12oh_oiiinii_moustakas_err

; Atlas

    indxatlas = where((atlasnodust.zstrong_12oh_oiiinii_moustakas gt -900) and (atlasnodust.mass_bv_b gt -900),nindxatlas)

    xatlas = atlasnodust[indxatlas].mass_bv_b
    xerratlas = atlasnodust[indxatlas].mass_bv_b_err

    yatlas = atlasnodust[indxatlas].zstrong_12oh_oiiinii_moustakas
    yerratlas = atlasnodust[indxatlas].zstrong_12oh_oiiinii_moustakas_err

; NFGS

    indxnfgs = where((nfgsnodust.zstrong_12oh_oiiinii_moustakas gt -900) and (nfgsnodust.mass_bv_b gt -900),nindxnfgs)

    xnfgs = nfgsnodust[indxnfgs].mass_bv_b
    xerrnfgs = nfgsnodust[indxnfgs].mass_bv_b_err

    ynfgs = nfgsnodust[indxnfgs].zstrong_12oh_oiiinii_moustakas
    yerrnfgs = nfgsnodust[indxnfgs].zstrong_12oh_oiiinii_moustakas_err

    xtitle = 'log (M / M'+sunsymbol()+')'
    ytitle = ohtitle+' ([O III]/H\beta)/([N II]/H\alpha)'

    xrange = massrange
    yrange = ohrange

    sdss_lineplot, x, y, xerr, yerr, plottype=3, postscript=postscript, $
      xtitle=xtitle, ytitle=ytitle, xrange=xrange, yrange=yrange, legendtype=0, $
      /right, /top, position=pos[*,0], charsize=charsize_8, $
      xatlas=xatlas, xerratlas=xerratlas, yatlas=yatlas, yerratlas=yerratlas, $
      xnfgs=xnfgs, xerrnfgs=xerrnfgs, ynfgs=ynfgs, yerrnfgs=yerrnfgs

; overplot the SDSS Tremonti result

    massaxis = findgen((11.5-8.5)/0.01+1)*0.01+8.5
    djs_oplot, massaxis, -1.492 + 1.847*massaxis - 0.08026*massaxis^2, $
      line=0, thick=postthick, color='dark green'
    
    im_openclose, postscript=postscript, /close    
    
; ------------------------------------------------------------
; O32 versus 12+log(O/H) [N II]/Ha - Integrated + SDSS
; ------------------------------------------------------------

    psname = 'o32_vs_12oh_niiha'
    im_openclose, pspath+psname, postscript=postscript, xsize=8.5, ysize=4.5, /encapsulated

    pagemaker, nx=2, ny=1, height=3.5, width=[3.4,3.4], xmargin=[1.1,0.3], $
      ymargin=[0.2,0.9], xspace=0.0, yspace=0, xpage=8.5, ypage=4.5, $
      position=pos, /normal

    OHup = 8.3
    OHlo = 8.0

    OHup_pilyugin = 8.2
    OHlo_pilyugin = 7.95

; Integrated + HII Regions    
    
    indx = where((atlasnodust.zstrong_12oh_niiha_pettini gt -900) and (atlasnodust.zstrong_o32 gt -900.0),nindx)

    x = atlasnodust[indx].zstrong_o32
    xerr = atlasnodust[indx].zstrong_o32_err

    y = atlasnodust[indx].zstrong_12oh_niiha_pettini
    yerr = atlasnodust[indx].zstrong_12oh_niiha_pettini_err

    indxnfgs = where((nfgsnodust.zstrong_12oh_niiha_pettini gt -900) and (nfgsnodust.zstrong_o32 gt -900.0),nindxnfgs)

    xnfgs = nfgsnodust[indxnfgs].zstrong_o32
    xerrnfgs = nfgsnodust[indxnfgs].zstrong_o32_err

    ynfgs = nfgsnodust[indxnfgs].zstrong_12oh_niiha_pettini
    yerrnfgs = nfgsnodust[indxnfgs].zstrong_12oh_niiha_pettini_err

    xtitle = 'log O_{32}'
    ytitle = ohtitle+' [N II]/H\alpha'

    xrange = o32range2
    yrange = ohrange

    good = where((hii.ZSTRONG_o32 gt -900.0) and (hii.zstrong_12oh_niiha_pettini gt -900))
    xregion = hii[good].ZSTRONG_o32 & xerrregion = hii[good].ZSTRONG_o32_err
    yregion = hii[good].zstrong_12oh_niiha_pettini & yerrregion = hii[good].zstrong_12oh_niiha_pettini_err

;   xbig = x & ybig = y
;   xbig = [x,xnfgs] & ybig = [y,ynfgs]
;   xbig = [x,xregion] & ybig = [y,yregion]
    xbig = [x,xnfgs,xregion] & ybig = [y,ynfgs,yregion]

    atlas1d_lineplot, x, y, xerr, yerr, plottype=3, postscript=postscript, $
      xtitle=xtitle, ytitle=ytitle, xrange=xrange, yrange=yrange, legendtype=0, $
      /right, /top, position=pos[*,0], charsize=charsize_6, $
      xregion=xregion, yregion=yregion, xerrregion=xerrregion, yerrregion=yerrregion, $
      xnfgs=xnfgs, ynfgs=ynfgs, xerrnfgs=xerrnfgs, yerrnfgs=yerrnfgs

; fit a line to the turn-around region

    minoh = 7.7
    maxoh = 8.5
    
    inrange = where((ybig gt minoh) and (ybig lt maxoh),ninrange)

;   coeff = linfit(xbig[inrange],ybig[inrange],sigma=coeff_err)
    coeff = robust_linefit(xbig[inrange],ybig[inrange],junk,junk2,coeff_err,/bisect)
    sixlin, xbig[inrange], ybig[inrange], a, siga, b, sigb
    coeff = [a[2],b[2]] ; Ordinary Least Squares Bisector
    coeff_err = [siga[2],sigb[2]]

    yfitdata = poly(xbig,coeff)

    bigxfit = findgen((!x.crange[1]-!x.crange[0])/0.01+1L)*0.01+!x.crange[0]
    bigyfit = poly(bigxfit,coeff)    

    fitrange = where((bigyfit gt minoh) and (bigyfit lt maxoh))
    xfit = bigxfit[fitrange]
    yfit = bigyfit[fitrange]
    
    djs_oplot, xfit, yfit, thick=postthick, line=0

    splog, 'Coefficients: '
    niceprint, coeff, coeff_err
 
    turnaround = where((yfitdata gt OHlo) and (yfitdata lt OHup),nturn)
    splog, 'Residuals in the turn-around region:'
    stats = im_stats(yfitdata[turnaround]-ybig[turnaround],/verbose)
    
; compute some numbers to calibrate the M91 calibration    

    o32up = interpol(xbig,ybig,OHup)
    o32lo = interpol(xbig,ybig,OHlo)
 
    oplot, !x.crange, OHup*[1,1], line=1, thick=postthick
    oplot, !x.crange, OHlo*[1,1], line=1, thick=postthick
 
    xyouts, -1.2, 7.6, 'Lower Branch', align=0.0, /data, charsize=charsize_6, $
      charthick=postthick
    xyouts, 1.3, 9.1, 'Upper Branch', align=1.0, /data, charsize=charsize_6, $
      charthick=postthick
 
    splog, '12+log(O/H) = '+string(OHlo,format='(F5.3)')+' --> '+string(o32lo,format='(F5.3)')
    splog, '12+log(O/H) = '+string(OHup,format='(F5.3)')+' --> '+string(o32up,format='(F5.3)')

; compute some numbers to calibrate the Pilyugin calibration    
    
    o32up_pilyugin = interpol(xbig,ybig,OHup_pilyugin)
    o32lo_pilyugin = interpol(xbig,ybig,OHlo_pilyugin)
    splog, '12+log(O/H) = '+string(OHlo_pilyugin,format='(F5.3)')+' --> '+string(o32lo_pilyugin,format='(F5.3)')
    splog, '12+log(O/H) = '+string(OHup_pilyugin,format='(F5.3)')+' --> '+string(o32up_pilyugin,format='(F5.3)')
    
; compute running statistics in the range of interest
    
;   result = im_medxbin(xbig,ybig,0.15);,minx=o32up,maxx=o32lo)
;   oploterror, result.meanx, result.meany, replicate(result.binsz/2.0,result.nbins), result.sigy, $
;     ps=4, thick=postthick, errthick=postthick

; SDSS
    
    indxsdss = where((sdssnodust.zstrong_12oh_niiha_pettini gt -900) and (sdssnodust.zstrong_o32 gt -900.0),nindxsdss)

    xsdss = sdssnodust[indxsdss].zstrong_o32
    xerrsdss = sdssnodust[indxsdss].zstrong_o32_err

    ysdss = sdssnodust[indxsdss].zstrong_12oh_niiha_pettini
    yerrsdss = sdssnodust[indxsdss].zstrong_12oh_niiha_pettini_err

    sdss_lineplot, xsdss, ysdss, xerrsdss, yerrsdss, plottype=3, postscript=postscript, $
      xtitle=xtitle, ytitle='', xrange=xrange, yrange=yrange, legendtype=0, $
      /right, /top, position=pos[*,1], ytickname=replicate(' ',10), /noerase, $
      charsize=charsize_6

; overplot the bisector    
    
    djs_oplot, xfit, yfit, thick=postthick, line=0

    yfitsdssdata = poly(xsdss,coeff)
    turnaround = where((yfitsdssdata gt OHlo) and (yfitsdssdata lt OHup),nturn)
    splog, 'SDSS residuals in the turn-around region:'
    stats = im_stats(yfitsdssdata[turnaround]-ysdss[turnaround],/verbose)
    
    oplot, !x.crange, OHup*[1,1], line=1, thick=postthick
    oplot, !x.crange, OHlo*[1,1], line=1, thick=postthick
 
    xyouts, -1.2, 7.6, 'Lower Branch', align=0.0, /data, charsize=charsize_6, $
      charthick=postthick
    xyouts, 1.3, 9.1, 'Upper Branch', align=1.0, /data, charsize=charsize_6, $
      charthick=postthick

    im_openclose, postscript=postscript, /close    

; ------------------------------------------------------------
; 12+log(O/H) P01_O32 versus 12+log(O/H) ([O III]/Hb)/([N II]/Ha) - Integrated + SDSS
; ------------------------------------------------------------

    psname = '12oh_p01o32_vs_12oh_oiiinii_niiha'
    im_openclose, pspath+psname, postscript=postscript, xsize=8.5, ysize=4.5, /encapsulated

    pagemaker, nx=2, ny=1, height=3.5, width=[3.4,3.4], xmargin=[1.1,0.3], $
      ymargin=[0.2,0.9], xspace=0.0, yspace=0, xpage=8.5, ypage=4.5, $
      position=pos, /normal

; Integrated + HII Regions    

    indx = where((atlasnodust.zstrong_12oh_p01_o32 gt -900) and (atlasnodust.zstrong_12oh_oiiinii_niiha gt -900),nindx)

    x = atlasnodust[indx].zstrong_12oh_p01_o32
    xerr = atlasnodust[indx].zstrong_12oh_p01_o32_err

    y = atlasnodust[indx].zstrong_12oh_oiiinii_niiha
    yerr = atlasnodust[indx].zstrong_12oh_oiiinii_niiha_err

    indxnfgs = where((nfgsnodust.zstrong_12oh_p01_o32 gt -900) and (nfgsnodust.zstrong_12oh_oiiinii_niiha gt -900),nindxnfgs)

    xnfgs = nfgsnodust[indxnfgs].zstrong_12oh_p01_o32
    xerrnfgs = nfgsnodust[indxnfgs].zstrong_12oh_p01_o32_err

    ynfgs = nfgsnodust[indxnfgs].zstrong_12oh_oiiinii_niiha
    yerrnfgs = nfgsnodust[indxnfgs].zstrong_12oh_oiiinii_niiha_err

    residuals = [y,ynfgs]-[x,xnfgs]
    stats = im_stats(residuals,/verbose)
    xstr = strtrim(string(stats.median,format='(F12.2)'),2)+'\pm'+strtrim(string(stats.sig68mean,format='(F12.2)'),2)

    xtitle = ohtitle+' R_{23}/P01 O_{32}'
    ytitle = ohtitle+' ([O III]/H\beta)/([N II]/H\alpha)'

    xrange = ohrange
    yrange = xrange

;   good = where((hii.zstrong_12oh_oiiinii_niiha gt -900) and (hii.zstrong_12oh_p01_o32 gt -900))
;   xregion = hii[good].zstrong_12oh_p01_o32 & xerrregion = hii[good].zstrong_12oh_p01_o32_err
;   yregion = hii[good].zstrong_12oh_oiiinii_niiha & yerrregion = hii[good].zstrong_12oh_oiiinii_niiha_err

    atlas1d_lineplot, x, y, xerr, yerr, plottype=3, postscript=postscript, $
      xtitle='', ytitle=ytitle, xrange=xrange, yrange=yrange, legendtype=0, $
      /right, /top, position=pos[*,0], charsize=charsize_2, $
;     xregion=xregion, yregion=yregion, xerrregion=xerrregion, yerrregion=yerrregion, $
      xnfgs=xnfgs, ynfgs=ynfgs, xerrnfgs=xerrnfgs, yerrnfgs=yerrnfgs
    djs_oplot, !x.crange, !y.crange, line=0, thick=postthick
    legend, textoidl(xstr), /left, /top, box=0, charsize=charsize_2, charthick=postthick

; SDSS
    
;   indxsdss = where((sdssnodust.zstrong_12oh_p01_o32 gt -900) and (sdssnodust.zstrong_12oh_oiiinii_niiha gt -900),nindxsdss)
;
;   xsdss = sdssnodust[indxsdss].zstrong_12oh_p01_o32
;   xerrsdss = sdssnodust[indxsdss].zstrong_12oh_p01_o32_err
;
;   ysdss = sdssnodust[indxsdss].zstrong_12oh_oiiinii_niiha
;   yerrsdss = sdssnodust[indxsdss].zstrong_12oh_oiiinii_niiha_err
;
;   residuals = ysdss-xsdss
;   stats = im_stats(residuals,/verbose)
;   xstr = strtrim(string(stats.median,format='(F12.2)'),2)+'\pm'+strtrim(string(stats.sig68mean,format='(F12.2)'),2)
;
;   sdss_lineplot, xsdss, ysdss, xerrsdss, yerrsdss, plottype=3, postscript=postscript, $
;     xtitle='', ytitle='', xrange=xrange, yrange=yrange, legendtype=0, $
;     /right, /top, position=pos[*,1], ytickname=replicate(' ',10), /noerase, $
;     charsize=charsize_2
;   djs_oplot, !x.crange, !y.crange, line=0, thick=postthick
;   legend, textoidl(xstr), /left, /top, box=0, charsize=charsize_2, charthick=postthick

; the title
    
    xyouts, 0.5*(pos[2,1]-pos[0,0])+pos[0,0], pos[1,0]*0.3, textoidl(xtitle), align=0.5, $
      charsize=charsize_2, charthick=postthick, /normal

    im_openclose, postscript=postscript, /close    

; ------------------------------------------------------------
; 12+log(O/H) P01 Upper versus 12+log(O/H) ([O III]/Hb)/([N II]/Ha) - Integrated + SDSS

; Don't bother comparing P01 Lower because there are too few objects
; in our galaxy sample with sufficiently low abundances
; ------------------------------------------------------------

    psname = '12oh_p01_upper_vs_12oh_oiiinii_niiha'
    im_openclose, pspath+psname, postscript=postscript, xsize=8.5, ysize=4.5, /encapsulated

    pagemaker, nx=2, ny=1, height=3.5, width=[3.4,3.4], xmargin=[1.1,0.3], $
      ymargin=[0.2,0.9], xspace=0.0, yspace=0, xpage=8.5, ypage=4.5, $
      position=pos, /normal

; Integrated + HII Regions    

    ohcut = -900.0 ; 8.2
    
    indx = where((atlasnodust.zstrong_12oh_p01_upper gt -900) and (atlasnodust.zstrong_12oh_oiiinii_niiha gt ohcut),nindx)

    x = atlasnodust[indx].zstrong_12oh_p01_upper
    xerr = atlasnodust[indx].zstrong_12oh_p01_upper_err

    y = atlasnodust[indx].zstrong_12oh_oiiinii_niiha
    yerr = atlasnodust[indx].zstrong_12oh_oiiinii_niiha_err

    indxnfgs = where((nfgsnodust.zstrong_12oh_p01_upper gt -900) and (nfgsnodust.zstrong_12oh_oiiinii_niiha gt ohcut),nindxnfgs)

    xnfgs = nfgsnodust[indxnfgs].zstrong_12oh_p01_upper
    xerrnfgs = nfgsnodust[indxnfgs].zstrong_12oh_p01_upper_err

    ynfgs = nfgsnodust[indxnfgs].zstrong_12oh_oiiinii_niiha
    yerrnfgs = nfgsnodust[indxnfgs].zstrong_12oh_oiiinii_niiha_err

    residuals = [y,ynfgs]-[x,xnfgs]
    stats = im_stats(residuals,/verbose)
    xstr = strtrim(string(stats.median,format='(F12.2)'),2)+'\pm'+strtrim(string(stats.sig68mean,format='(F12.2)'),2)

    xtitle = ohtitle+' Pilyugin Upper Branch'
    ytitle = ohtitle+' ([O III]/H\beta)/([N II]/H\alpha)'

    xrange = ohrange
    yrange = xrange

    good = where((hii.zstrong_12oh_oiiinii_niiha gt ohcut) and (hii.zstrong_12oh_p01_upper gt -900))
    xregion = hii[good].zstrong_12oh_p01_upper & xerrregion = hii[good].zstrong_12oh_p01_upper_err
    yregion = hii[good].zstrong_12oh_oiiinii_niiha & yerrregion = hii[good].zstrong_12oh_oiiinii_niiha_err

    xbig = [x,xnfgs,xregion]
    ybig = [y,ynfgs,yregion]

    atlas1d_lineplot, x, y, xerr, yerr, plottype=3, postscript=postscript, $
      xtitle='', ytitle=ytitle, xrange=xrange, yrange=yrange, legendtype=0, $
      /right, /top, position=pos[*,0], charsize=charsize_2, $
      xregion=xregion, yregion=yregion, xerrregion=xerrregion, yerrregion=yerrregion, $
      xnfgs=xnfgs, ynfgs=ynfgs, xerrnfgs=xerrnfgs, yerrnfgs=yerrnfgs
    djs_oplot, !x.crange, !y.crange, line=0, thick=postthick
    legend, textoidl(xstr), /left, /top, box=0, charsize=charsize_2, charthick=postthick

; SDSS
    
    indxsdss = where((sdssnodust.zstrong_12oh_p01_upper gt -900) and (sdssnodust.zstrong_12oh_oiiinii_niiha gt ohcut),nindxsdss)

    xsdss = sdssnodust[indxsdss].zstrong_12oh_p01_upper
    xerrsdss = sdssnodust[indxsdss].zstrong_12oh_p01_upper_err

    ysdss = sdssnodust[indxsdss].zstrong_12oh_oiiinii_niiha
    yerrsdss = sdssnodust[indxsdss].zstrong_12oh_oiiinii_niiha_err

    residuals = ysdss-xsdss
    stats = im_stats(residuals,/verbose)
    xstr = strtrim(string(stats.median,format='(F12.2)'),2)+'\pm'+strtrim(string(stats.sig68mean,format='(F12.2)'),2)

    sdss_lineplot, xsdss, ysdss, xerrsdss, yerrsdss, plottype=3, postscript=postscript, $
      xtitle='', ytitle='', xrange=xrange, yrange=yrange, legendtype=0, $
      /right, /top, position=pos[*,1], ytickname=replicate(' ',10), /noerase, $
      charsize=charsize_2
    djs_oplot, !x.crange, !y.crange, line=0, thick=postthick
    legend, textoidl(xstr), /left, /top, box=0, charsize=charsize_2, charthick=postthick

; the title
    
    xyouts, 0.5*(pos[2,1]-pos[0,0])+pos[0,0], pos[1,0]*0.3, textoidl(xtitle), align=0.5, $
      charsize=charsize_2, charthick=postthick, /normal

    im_openclose, postscript=postscript, /close    

; ------------------------------------------------------------
; 12+log(O/H) P01 Upper versus 12+log(O/H) [N II]/Ha - [Integrated + SDSS]
; ------------------------------------------------------------

    psname = '12oh_p01_upper_vs_12oh_niiha'
    im_openclose, pspath+psname, postscript=postscript, xsize=8.5, ysize=4.5, /encapsulated

    pagemaker, nx=2, ny=1, height=3.5, width=[3.4,3.4], xmargin=[1.1,0.3], $
      ymargin=[0.2,0.9], xspace=0.0, yspace=0, xpage=8.5, ypage=4.5, $
      position=pos, /normal

; Integrated + HII Regions    

    ohcut = -900.0 ; 8.2
    
    indx = where((atlasnodust.zstrong_12oh_p01_upper gt -900) and (atlasnodust.zstrong_12oh_niiha_pettini gt ohcut),nindx)

    x = atlasnodust[indx].zstrong_12oh_p01_upper
    xerr = atlasnodust[indx].zstrong_12oh_p01_upper_err

    y = atlasnodust[indx].zstrong_12oh_niiha_pettini
    yerr = atlasnodust[indx].zstrong_12oh_niiha_pettini_err

    indxnfgs = where((nfgsnodust.zstrong_12oh_p01_upper gt -900) and (nfgsnodust.zstrong_12oh_niiha_pettini gt ohcut),nindxnfgs)

    xnfgs = nfgsnodust[indxnfgs].zstrong_12oh_p01_upper
    xerrnfgs = nfgsnodust[indxnfgs].zstrong_12oh_p01_upper_err

    ynfgs = nfgsnodust[indxnfgs].zstrong_12oh_niiha_pettini
    yerrnfgs = nfgsnodust[indxnfgs].zstrong_12oh_niiha_pettini_err

    residuals = [y,ynfgs]-[x,xnfgs]
    stats = im_stats(residuals,/verbose)
    xstr = strtrim(string(stats.median,format='(F12.2)'),2)+'\pm'+strtrim(string(stats.sig68mean,format='(F12.2)'),2)

    xtitle = ohtitle+' Pilyugin Upper Branch'
    ytitle = ohtitle+' [N II]/H\alpha'

    xrange = ohrange
    yrange = xrange

    good = where((hii.zstrong_12oh_niiha_pettini gt ohcut) and (hii.zstrong_12oh_p01_upper gt -900))
    xregion = hii[good].zstrong_12oh_p01_upper & xerrregion = hii[good].zstrong_12oh_p01_upper_err
    yregion = hii[good].zstrong_12oh_niiha_pettini & yerrregion = hii[good].zstrong_12oh_niiha_pettini_err

    xbig = [x,xnfgs,xregion]
    ybig = [y,ynfgs,yregion]

    atlas1d_lineplot, x, y, xerr, yerr, plottype=3, postscript=postscript, $
      xtitle='', ytitle=ytitle, xrange=xrange, yrange=yrange, legendtype=0, $
      /right, /top, position=pos[*,0], charsize=charsize_2, $
      xregion=xregion, yregion=yregion, xerrregion=xerrregion, yerrregion=yerrregion, $
      xnfgs=xnfgs, ynfgs=ynfgs, xerrnfgs=xerrnfgs, yerrnfgs=yerrnfgs
    djs_oplot, !x.crange, !y.crange, line=0, thick=postthick
    legend, textoidl(xstr), /left, /top, box=0, charsize=charsize_2, charthick=postthick

; SDSS
    
    indxsdss = where((sdssnodust.zstrong_12oh_p01_upper gt -900) and (sdssnodust.zstrong_12oh_niiha_pettini gt ohcut),nindxsdss)

    xsdss = sdssnodust[indxsdss].zstrong_12oh_p01_upper
    xerrsdss = sdssnodust[indxsdss].zstrong_12oh_p01_upper_err

    ysdss = sdssnodust[indxsdss].zstrong_12oh_niiha_pettini
    yerrsdss = sdssnodust[indxsdss].zstrong_12oh_niiha_pettini_err

    residuals = ysdss-xsdss
    stats = im_stats(residuals,/verbose)
    xstr = strtrim(string(stats.median,format='(F12.2)'),2)+'\pm'+strtrim(string(stats.sig68mean,format='(F12.2)'),2)

    sdss_lineplot, xsdss, ysdss, xerrsdss, yerrsdss, plottype=3, postscript=postscript, $
      xtitle='', ytitle='', xrange=xrange, yrange=yrange, legendtype=0, $
      /right, /top, position=pos[*,1], ytickname=replicate(' ',10), /noerase, $
      charsize=charsize_2
    djs_oplot, !x.crange, !y.crange, line=0, thick=postthick
    legend, textoidl(xstr), /left, /top, box=0, charsize=charsize_2, charthick=postthick

; the title
    
    xyouts, 0.5*(pos[2,1]-pos[0,0])+pos[0,0], pos[1,0]*0.3, textoidl(xtitle), align=0.5, $
      charsize=charsize_2, charthick=postthick, /normal

    im_openclose, postscript=postscript, /close    

; ------------------------------------------------------------
; 12+log(O/H) R23/M91_O32 versus 12+log(O/H) ([O III]/Hb)/([N II]/Ha) - Integrated + SDSS
; ------------------------------------------------------------

    psname = '12oh_m91o32_vs_12oh_oiiinii_niiha'
    im_openclose, pspath+psname, postscript=postscript, xsize=8.5, ysize=4.5, /encapsulated

    pagemaker, nx=2, ny=1, height=3.5, width=[3.4,3.4], xmargin=[1.1,0.3], $
      ymargin=[0.2,0.9], xspace=0.0, yspace=0, xpage=8.5, ypage=4.5, $
      position=pos, /normal

; Integrated + HII Regions    

    indx = where((atlasnodust.zstrong_12oh_m91_o32 gt -900) and (atlasnodust.zstrong_12oh_oiiinii_niiha gt -900),nindx)

    x = atlasnodust[indx].zstrong_12oh_m91_o32
    xerr = atlasnodust[indx].zstrong_12oh_m91_o32_err

    y = atlasnodust[indx].zstrong_12oh_oiiinii_niiha
    yerr = atlasnodust[indx].zstrong_12oh_oiiinii_niiha_err

    indxnfgs = where((nfgsnodust.zstrong_12oh_m91_o32 gt -900) and (nfgsnodust.zstrong_12oh_oiiinii_niiha gt -900),nindxnfgs)

    xnfgs = nfgsnodust[indxnfgs].zstrong_12oh_m91_o32
    xerrnfgs = nfgsnodust[indxnfgs].zstrong_12oh_m91_o32_err

    ynfgs = nfgsnodust[indxnfgs].zstrong_12oh_oiiinii_niiha
    yerrnfgs = nfgsnodust[indxnfgs].zstrong_12oh_oiiinii_niiha_err

    residuals = [y,ynfgs]-[x,xnfgs]
    stats = im_stats(residuals,/verbose)
    xstr = strtrim(string(stats.median,format='(F12.2)'),2)+'\pm'+strtrim(string(stats.sig68mean,format='(F12.2)'),2)

    xtitle = ohtitle+' R_{23}/M91 O_{32}'
    ytitle = ohtitle+' ([O III]/H\beta)/([N II]/H\alpha)'

    xrange = ohrange
    yrange = xrange

    good = where((hii.zstrong_12oh_oiiinii_niiha gt -900) and (hii.zstrong_12oh_m91_o32 gt -900))
    xregion = hii[good].zstrong_12oh_m91_o32 & xerrregion = hii[good].zstrong_12oh_m91_o32_err
    yregion = hii[good].zstrong_12oh_oiiinii_niiha & yerrregion = hii[good].zstrong_12oh_oiiinii_niiha_err

    xbig = [x,xnfgs,xregion]
    ybig = [y,ynfgs,yregion]

    atlas1d_lineplot, x, y, xerr, yerr, plottype=3, postscript=postscript, $
      xtitle='', ytitle=ytitle, xrange=xrange, yrange=yrange, legendtype=0, $
      /right, /top, position=pos[*,0], charsize=charsize_2, $
      xregion=xregion, yregion=yregion, xerrregion=xerrregion, yerrregion=yerrregion, $
      xnfgs=xnfgs, ynfgs=ynfgs, xerrnfgs=xerrnfgs, yerrnfgs=yerrnfgs
    djs_oplot, !x.crange, !y.crange, line=0, thick=postthick
    legend, textoidl(xstr), /left, /top, box=0, charsize=charsize_2, charthick=postthick

; SDSS
    
    indxsdss = where((sdssnodust.zstrong_12oh_m91_o32 gt -900) and (sdssnodust.zstrong_12oh_oiiinii_niiha gt -900),nindxsdss)

    xsdss = sdssnodust[indxsdss].zstrong_12oh_m91_o32
    xerrsdss = sdssnodust[indxsdss].zstrong_12oh_m91_o32_err

    ysdss = sdssnodust[indxsdss].zstrong_12oh_oiiinii_niiha
    yerrsdss = sdssnodust[indxsdss].zstrong_12oh_oiiinii_niiha_err

    residuals = ysdss-xsdss
    stats = im_stats(residuals,/verbose)
    xstr = strtrim(string(stats.median,format='(F12.2)'),2)+'\pm'+strtrim(string(stats.sig68mean,format='(F12.2)'),2)

    sdss_lineplot, xsdss, ysdss, xerrsdss, yerrsdss, plottype=3, postscript=postscript, $
      xtitle='', ytitle='', xrange=xrange, yrange=yrange, legendtype=0, $
      /right, /top, position=pos[*,1], ytickname=replicate(' ',10), /noerase, $
      charsize=charsize_2
    djs_oplot, !x.crange, !y.crange, line=0, thick=postthick
    legend, textoidl(xstr), /left, /top, box=0, charsize=charsize_2, charthick=postthick

; the title
    
    xyouts, 0.5*(pos[2,1]-pos[0,0])+pos[0,0], pos[1,0]*0.3, textoidl(xtitle), align=0.5, $
      charsize=charsize_2, charthick=postthick, /normal

    im_openclose, postscript=postscript, /close    

; ------------------------------------------------------------
; 12+log(O/H) R23/M91 Upper versus 12+log(O/H) ([O III]/Hb)/([N II]/Ha) - Integrated + SDSS
; ------------------------------------------------------------

    psname = '12oh_m91_upper_vs_12oh_oiiinii_niiha'
    im_openclose, pspath+psname, postscript=postscript, xsize=8.5, ysize=4.5, /encapsulated

    pagemaker, nx=2, ny=1, height=3.5, width=[3.4,3.4], xmargin=[1.1,0.3], $
      ymargin=[0.2,0.9], xspace=0.0, yspace=0, xpage=8.5, ypage=4.5, $
      position=pos, /normal

; Integrated + HII Regions    

    ohcut = 8.3
    
    indx = where((atlasnodust.zstrong_12oh_m91_upper gt -900) and (atlasnodust.zstrong_12oh_oiiinii_niiha gt ohcut),nindx)
;   indx = where((atlasnodust.zstrong_12oh_m91_upper gt -900) and (atlasnodust.zstrong_12oh_oiiinii_niiha gt -900.0),nindx)

    x = atlasnodust[indx].zstrong_12oh_m91_upper
    xerr = atlasnodust[indx].zstrong_12oh_m91_upper_err

    y = atlasnodust[indx].zstrong_12oh_oiiinii_niiha
    yerr = atlasnodust[indx].zstrong_12oh_oiiinii_niiha_err

    indxnfgs = where((nfgsnodust.zstrong_12oh_m91_upper gt -900) and (nfgsnodust.zstrong_12oh_oiiinii_niiha gt ohcut),nindxnfgs)
;   indxnfgs = where((nfgsnodust.zstrong_12oh_m91_upper gt -900) and (nfgsnodust.zstrong_12oh_oiiinii_niiha gt -900.0),nindxnfgs)

    xnfgs = nfgsnodust[indxnfgs].zstrong_12oh_m91_upper
    xerrnfgs = nfgsnodust[indxnfgs].zstrong_12oh_m91_upper_err

    ynfgs = nfgsnodust[indxnfgs].zstrong_12oh_oiiinii_niiha
    yerrnfgs = nfgsnodust[indxnfgs].zstrong_12oh_oiiinii_niiha_err

    residuals = [y,ynfgs]-[x,xnfgs]
    stats = im_stats(residuals,/verbose)
    xstr = strtrim(string(stats.median,format='(F12.2)'),2)+'\pm'+strtrim(string(stats.sig68mean,format='(F12.2)'),2)

    xtitle = ohtitle+' R_{23}/M91 Upper Branch'
    ytitle = ohtitle+' ([O III]/H\beta)/([N II]/H\alpha)'

    xrange = ohrange
    yrange = xrange

    good = where((hii.zstrong_12oh_oiiinii_niiha gt ohcut) and (hii.zstrong_12oh_m91_upper gt -900))
;   good = where((hii.zstrong_12oh_oiiinii_niiha gt -900.0) and (hii.zstrong_12oh_m91_upper gt -900))
    xregion = hii[good].zstrong_12oh_m91_upper & xerrregion = hii[good].zstrong_12oh_m91_upper_err
    yregion = hii[good].zstrong_12oh_oiiinii_niiha & yerrregion = hii[good].zstrong_12oh_oiiinii_niiha_err

    xbig = [x,xnfgs,xregion]
    ybig = [y,ynfgs,yregion]

    atlas1d_lineplot, x, y, xerr, yerr, plottype=3, postscript=postscript, $
      xtitle='', ytitle=ytitle, xrange=xrange, yrange=yrange, legendtype=0, $
      /right, /top, position=pos[*,0], charsize=charsize_2, $
      xregion=xregion, yregion=yregion, xerrregion=xerrregion, yerrregion=yerrregion, $
      xnfgs=xnfgs, ynfgs=ynfgs, xerrnfgs=xerrnfgs, yerrnfgs=yerrnfgs
    djs_oplot, !x.crange, !y.crange, line=0, thick=postthick
    legend, textoidl(xstr), /left, /top, box=0, charsize=charsize_2, charthick=postthick

; SDSS
    
    indxsdss = where((sdssnodust.zstrong_12oh_m91_upper gt -900) and (sdssnodust.zstrong_12oh_oiiinii_niiha gt ohcut),nindxsdss)
;   indxsdss = where((sdssnodust.zstrong_12oh_m91_upper gt -900) and (sdssnodust.zstrong_12oh_oiiinii_niiha gt -900.0),nindxsdss)

    xsdss = sdssnodust[indxsdss].zstrong_12oh_m91_upper
    xerrsdss = sdssnodust[indxsdss].zstrong_12oh_m91_upper_err

    ysdss = sdssnodust[indxsdss].zstrong_12oh_oiiinii_niiha
    yerrsdss = sdssnodust[indxsdss].zstrong_12oh_oiiinii_niiha_err

    residuals = ysdss-xsdss
    stats = im_stats(residuals,/verbose)
    xstr = strtrim(string(stats.median,format='(F12.2)'),2)+'\pm'+strtrim(string(stats.sig68mean,format='(F12.2)'),2)

    sdss_lineplot, xsdss, ysdss, xerrsdss, yerrsdss, plottype=3, postscript=postscript, $
      xtitle='', ytitle='', xrange=xrange, yrange=yrange, legendtype=0, $
      /right, /top, position=pos[*,1], ytickname=replicate(' ',10), /noerase, $
      charsize=charsize_2
    djs_oplot, !x.crange, !y.crange, line=0, thick=postthick
    legend, textoidl(xstr), /left, /top, box=0, charsize=charsize_2, charthick=postthick

; the title
    
    xyouts, 0.5*(pos[2,1]-pos[0,0])+pos[0,0], pos[1,0]*0.3, textoidl(xtitle), align=0.5, $
      charsize=charsize_2, charthick=postthick, /normal
    
    im_openclose, postscript=postscript, /close    

; ------------------------------------------------------------
; 12+log(O/H) R23/M91 Lower versus 12+log(O/H) ([O III]/Hb)/([N II]/Ha) - Integrated + SDSS
; ------------------------------------------------------------

    psname = '12oh_m91_lower_vs_12oh_oiiinii_niiha'
    im_openclose, pspath+psname, postscript=postscript, xsize=8.5, ysize=4.5, /encapsulated

    pagemaker, nx=2, ny=1, height=3.5, width=[3.4,3.4], xmargin=[1.1,0.3], $
      ymargin=[0.2,0.9], xspace=0.0, yspace=0, xpage=8.5, ypage=4.5, $
      position=pos, /normal

; Integrated + HII Regions    

    ohcut = 8.2
    
    indx = where((atlasnodust.zstrong_12oh_m91_lower gt -900) and (atlasnodust.zstrong_12oh_oiiinii_niiha lt ohcut) and $
      (atlasnodust.zstrong_12oh_oiiinii_niiha gt -900),nindx)

    x = atlasnodust[indx].zstrong_12oh_m91_lower
    xerr = atlasnodust[indx].zstrong_12oh_m91_lower_err

    y = atlasnodust[indx].zstrong_12oh_oiiinii_niiha
    yerr = atlasnodust[indx].zstrong_12oh_oiiinii_niiha_err

    indxnfgs = where((nfgsnodust.zstrong_12oh_m91_lower gt -900) and (nfgsnodust.zstrong_12oh_oiiinii_niiha lt ohcut) and $
      (nfgsnodust.zstrong_12oh_oiiinii_niiha gt -900),nindxnfgs)

    xnfgs = nfgsnodust[indxnfgs].zstrong_12oh_m91_lower
    xerrnfgs = nfgsnodust[indxnfgs].zstrong_12oh_m91_lower_err

    ynfgs = nfgsnodust[indxnfgs].zstrong_12oh_oiiinii_niiha
    yerrnfgs = nfgsnodust[indxnfgs].zstrong_12oh_oiiinii_niiha_err

    residuals = [y,ynfgs]-[x,xnfgs]
    stats = im_stats(residuals,/verbose)
    xstr = strtrim(string(stats.median,format='(F12.2)'),2)+'\pm'+strtrim(string(stats.sig68mean,format='(F12.2)'),2)

    xtitle = ohtitle+' R_{23}/M91 Lower Branch'
    ytitle = ohtitle+' ([O III]/H\beta)/([N II]/H\alpha)'

    xrange = ohrange
    yrange = xrange

    good = where((hii.zstrong_12oh_oiiinii_niiha gt -900) and (hii.zstrong_12oh_oiiinii_niiha lt ohcut) and (hii.zstrong_12oh_m91_lower gt -900))
    xregion = hii[good].zstrong_12oh_m91_lower & xerrregion = hii[good].zstrong_12oh_m91_lower_err
    yregion = hii[good].zstrong_12oh_oiiinii_niiha & yerrregion = hii[good].zstrong_12oh_oiiinii_niiha_err

    xbig = [x,xnfgs,xregion]
    ybig = [y,ynfgs,yregion]

    atlas1d_lineplot, x, y, xerr, yerr, plottype=3, postscript=postscript, $
      xtitle='', ytitle=ytitle, xrange=xrange, yrange=yrange, legendtype=0, $
      /right, /top, position=pos[*,0], charsize=charsize_2, $
      xregion=xregion, yregion=yregion, xerrregion=xerrregion, yerrregion=yerrregion, $
      xnfgs=xnfgs, ynfgs=ynfgs, xerrnfgs=xerrnfgs, yerrnfgs=yerrnfgs
    djs_oplot, !x.crange, !y.crange, line=0, thick=postthick
    legend, textoidl(xstr), /left, /top, box=0, charsize=charsize_2, charthick=postthick

; SDSS
    
    indxsdss = where((sdssnodust.zstrong_12oh_m91_lower gt -900) and (sdssnodust.zstrong_12oh_oiiinii_niiha lt ohcut) and $
      (sdssnodust.zstrong_12oh_oiiinii_niiha gt -900),nindxsdss)

    xsdss = sdssnodust[indxsdss].zstrong_12oh_m91_lower
    xerrsdss = sdssnodust[indxsdss].zstrong_12oh_m91_lower_err

    ysdss = sdssnodust[indxsdss].zstrong_12oh_oiiinii_niiha
    yerrsdss = sdssnodust[indxsdss].zstrong_12oh_oiiinii_niiha_err

    residuals = ysdss-xsdss
    stats = im_stats(residuals,/verbose)
    xstr = strtrim(string(stats.median,format='(F12.2)'),2)+'\pm'+strtrim(string(stats.sig68mean,format='(F12.2)'),2)

    sdss_lineplot, xsdss, ysdss, xerrsdss, yerrsdss, plottype=3, postscript=postscript, $
      xtitle='', ytitle='', xrange=xrange, yrange=yrange, legendtype=0, $
      /right, /top, position=pos[*,1], ytickname=replicate(' ',10), /noerase, $
      charsize=charsize_2
    djs_oplot, !x.crange, !y.crange, line=0, thick=postthick
    legend, textoidl(xstr), /left, /top, box=0, charsize=charsize_2, charthick=postthick

; the title
    
    xyouts, 0.5*(pos[2,1]-pos[0,0])+pos[0,0], pos[1,0]*0.3, textoidl(xtitle), align=0.5, $
      charsize=charsize_2, charthick=postthick, /normal
    
    im_openclose, postscript=postscript, /close    

; ------------------------------------------------------------
; 12+log(O/H) ([O III]/Hb)/([N II]/Ha) versus 12+log(O/H) [N II]/Ha - Integrated + SDSS
; ------------------------------------------------------------

    psname = '12oh_oiiinii_niiha_vs_12oh_niiha'
    im_openclose, pspath+psname, postscript=postscript, xsize=8.5, ysize=4.5, /encapsulated

    pagemaker, nx=2, ny=1, height=3.5, width=[3.4,3.4], xmargin=[1.1,0.3], $
      ymargin=[0.2,0.9], xspace=0.0, yspace=0, xpage=8.5, ypage=4.5, $
      position=pos, /normal

; Integrated + HII Regions    
    
    indx = where((atlasnodust.zstrong_12oh_niiha_pettini gt -900) and (atlasnodust.zstrong_12oh_oiiinii_pettini gt -900.0),nindx)

    x = atlasnodust[indx].zstrong_12oh_oiiinii_pettini
    xerr = atlasnodust[indx].zstrong_12oh_oiiinii_pettini_err

    y = atlasnodust[indx].zstrong_12oh_niiha_pettini
    yerr = atlasnodust[indx].zstrong_12oh_niiha_pettini_err

    indxnfgs = where((nfgsnodust.zstrong_12oh_niiha_pettini gt -900) and (nfgsnodust.zstrong_12oh_oiiinii_pettini gt -900.0),nindxnfgs)

    xnfgs = nfgsnodust[indxnfgs].zstrong_12oh_oiiinii_pettini
    xerrnfgs = nfgsnodust[indxnfgs].zstrong_12oh_oiiinii_pettini_err

    ynfgs = nfgsnodust[indxnfgs].zstrong_12oh_niiha_pettini
    yerrnfgs = nfgsnodust[indxnfgs].zstrong_12oh_niiha_pettini_err

    residuals = [y,ynfgs]-[x,xnfgs]
    stats = im_stats(residuals,/verbose)
    xstr = strtrim(string(stats.median,format='(F12.2)'),2)+'\pm'+strtrim(string(stats.sig68mean,format='(F12.2)'),2)

    xtitle = ohtitle+' ([O III]/H\beta)/([N II]/H\alpha)'
    ytitle = ohtitle+' [N II]/H\alpha'

    xrange = ohrange
    yrange = xrange

    good = where((hii.zstrong_12oh_oiiinii_pettini gt -900.0) and (hii.zstrong_12oh_niiha_pettini gt -900))
    xregion = hii[good].zstrong_12oh_oiiinii_pettini & xerrregion = hii[good].zstrong_12oh_oiiinii_pettini_err
    yregion = hii[good].zstrong_12oh_niiha_pettini & yerrregion = hii[good].zstrong_12oh_niiha_pettini_err

    xbig = [x,xnfgs,xregion]
    ybig = [y,ynfgs,yregion]

    atlas1d_lineplot, x, y, xerr, yerr, plottype=3, postscript=postscript, $
      xtitle='', ytitle=ytitle, xrange=xrange, yrange=yrange, legendtype=0, $
      /right, /top, position=pos[*,0], charsize=charsize_6, $
      xregion=xregion, yregion=yregion, xerrregion=xerrregion, yerrregion=yerrregion, $
      xnfgs=xnfgs, ynfgs=ynfgs, xerrnfgs=xerrnfgs, yerrnfgs=yerrnfgs
    djs_oplot, !x.crange, !y.crange, line=0, thick=postthick
    legend, textoidl(xstr), /left, /top, box=0, charsize=charsize_6, charthick=postthick

; SDSS
    
    indxsdss = where((sdssnodust.zstrong_12oh_niiha_pettini gt -900) and (sdssnodust.zstrong_12oh_oiiinii_pettini gt -900.0),nindxsdss)

    xsdss = sdssnodust[indxsdss].zstrong_12oh_oiiinii_pettini
    xerrsdss = sdssnodust[indxsdss].zstrong_12oh_oiiinii_pettini_err

    ysdss = sdssnodust[indxsdss].zstrong_12oh_niiha_pettini
    yerrsdss = sdssnodust[indxsdss].zstrong_12oh_niiha_pettini_err

    residuals = ysdss-xsdss
    stats = im_stats(residuals,/verbose)
    xstr = strtrim(string(stats.median,format='(F12.2)'),2)+'\pm'+strtrim(string(stats.sig68mean,format='(F12.2)'),2)

    sdss_lineplot, xsdss, ysdss, xerrsdss, yerrsdss, plottype=3, postscript=postscript, $
      xtitle='', ytitle='', xrange=xrange, yrange=yrange, legendtype=0, $
      /right, /top, position=pos[*,1], ytickname=replicate(' ',10), /noerase, $
      charsize=charsize_6
    djs_oplot, !x.crange, !y.crange, line=0, thick=postthick
    legend, textoidl(xstr), /left, /top, box=0, charsize=charsize_6, charthick=postthick

; the title
    
    xyouts, 0.5*(pos[2,1]-pos[0,0])+pos[0,0], pos[1,0]*0.2, textoidl(xtitle), align=0.5, $
      charsize=charsize_6, charthick=postthick, /normal
    
    im_openclose, postscript=postscript, /close    

; ------------------------------------------------------------
; 12+log(O/H) ([O III]/Hb)/([N II]/Ha) and [N II]/Ha versus 12+log(O/H) [Strong]
; ------------------------------------------------------------

    psname = '12oh_oiiinii_niiha_vs_12oh_strong'
    im_openclose, pspath+psname, postscript=postscript, /encapsulated, xsize=8.5, ysize=10.5

    pagemaker, nx=1, ny=3, height=[3.0,3.0,3.0], width=6.5, xmargin=[1.5,0.5], $
      ymargin=[0.4,1.1], xspace=0, yspace=0, xpage=8.5, ypage=10.5, $
      position=pos, /normal

    xtitle = ohtitle+' ([O III]/H\beta)/([N II]/H\alpha)'
;   xtitle = '12 + log (O/H)'

    xrange = ohrange3
    yrange = xrange

; --------------------------------------------------    
; KD02-Combined
; --------------------------------------------------    

    indx = where((atlasnodust.zstrong_12oh_oiiinii_niiha gt -900) and (atlasnodust.ZSTRONG_12oh_kd02_niioii gt -900.0),nindx)

    x = atlasnodust[indx].ZSTRONG_12oh_oiiinii_niiha
    xerr = atlasnodust[indx].ZSTRONG_12oh_oiiinii_niiha_err

    y = atlasnodust[indx].zstrong_12oh_kd02_niioii
    yerr = atlasnodust[indx].zstrong_12oh_kd02_niioii_err

    indxnfgs = where((nfgsnodust.zstrong_12oh_oiiinii_niiha gt -900) and (nfgsnodust.ZSTRONG_12oh_kd02_niioii gt -900.0),nindxnfgs)

    xnfgs = nfgsnodust[indxnfgs].ZSTRONG_12oh_oiiinii_niiha
    xerrnfgs = nfgsnodust[indxnfgs].ZSTRONG_12oh_oiiinii_niiha_err

    ynfgs = nfgsnodust[indxnfgs].zstrong_12oh_kd02_niioii
    yerrnfgs = nfgsnodust[indxnfgs].zstrong_12oh_kd02_niioii_err

    xbig = [x,xnfgs]
    ybig = [y,ynfgs]
    
    ytitle = ohtitle+' [N II]/[O II]'
;   ytitle = '12 + log (O/H)'

    stats = im_stats(xbig-ybig,/verbose)
    xstr = strtrim(string(stats.median,format='(F12.2)'),2)+'\pm'+strtrim(string(stats.sigma_rej,format='(F12.2)'),2)
    
    atlas1d_lineplot, x, y, xerr, yerr, plottype=1, postscript=postscript, $
      xtitle='', ytitle=ytitle, xrange=xrange, yrange=yrange, legendtype=0, $
      xstyle=3, ystyle=3, /right, /top, position=pos[*,0], charsize=charsize_4, $
      xtickname=replicate(' ',10), xnfgs=xnfgs, ynfgs=ynfgs, $
      xerrnfgs=xerrnfgs, yerrnfgs=yerrnfgs
    djs_oplot, !x.crange, !y.crange, thick=postthick

    legend, textoidl(xstr), /right, /bottom, box=0, charsize=charsize_4, charthick=postthick
;   label = ['(a) '+textoidl('12+log(O/H) from [N II]/H\alpha')]
;   legend, label, charthick=postthick, charsize=charsize_4, $
;     box=0, /left, /top
    legend, '(a)', /left, /top, box=0, charsize=charsize_4, charthick=postthick
    
; --------------------------------------------------    
; R23/M91
; --------------------------------------------------    

    indx = where((atlasnodust.zstrong_12oh_oiiinii_niiha gt -900) and (atlasnodust.ZSTRONG_12oh_m91_upper gt -900.0) and $
      (atlasnodust.ZSTRONG_12oh_m91_lower gt -900.0),nindx)

    x = atlasnodust[indx].zstrong_12oh_oiiinii_niiha
    xerr = atlasnodust[indx].zstrong_12oh_oiiinii_niiha_err
    
    y = x*0.0
    yerr = xerr*0.0
    
    for idiv = 0L, ndiv-1L do begin

       up = where((x gt div[idiv]),nup,comp=lo,ncomp=nlo)
       if (nup ne 0L) then begin
          y[up] = atlasnodust[indx[up]].zstrong_12oh_m91_upper
          yerr[up] = atlasnodust[indx[up]].zstrong_12oh_m91_upper_err
       endif

       if (nlo ne 0L) then begin
          y[lo] = atlasnodust[indx[lo]].zstrong_12oh_m91_lower
          yerr[lo] = atlasnodust[indx[lo]].zstrong_12oh_m91_lower_err
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
    splog, 'o3n2+n2 vs m91: '+string(mindiv,format='(F4.2)')

    up = where((atlasnodust[indx].zstrong_12oh_oiiinii_niiha gt mindiv),nup)
    if (nup ne 0L) then begin
       x = atlasnodust[indx[up]].zstrong_12oh_oiiinii_niiha
       xerr = atlasnodust[indx[up]].zstrong_12oh_oiiinii_niiha_err

       y = atlasnodust[indx[up]].zstrong_12oh_m91_upper
       yerr = atlasnodust[indx[up]].zstrong_12oh_m91_upper_err
    endif

    lo = where((atlasnodust[indx].zstrong_12oh_oiiinii_niiha lt mindiv),nlo)
    if (nlo ne 0L) then begin
       x = [x,atlasnodust[indx[lo]].zstrong_12oh_oiiinii_niiha]
       xerr = [xerr,atlasnodust[indx[lo]].zstrong_12oh_oiiinii_niiha_err]

       y = [y,atlasnodust[indx[lo]].zstrong_12oh_m91_lower]
       yerr = [yerr,atlasnodust[indx[lo]].zstrong_12oh_m91_lower_err]
    endif

    indxnfgs = where((nfgsnodust.zstrong_12oh_oiiinii_niiha gt -900) and (nfgsnodust.ZSTRONG_12oh_m91_upper gt -900.0) and $
      (nfgsnodust.ZSTRONG_12oh_m91_lower gt -900.0),nindxnfgs)

    xnfgs = nfgsnodust[indxnfgs].zstrong_12oh_oiiinii_niiha
    xnfgserr = nfgsnodust[indxnfgs].zstrong_12oh_oiiinii_niiha_err
    
    ynfgs = xnfgs*0.0
    ynfgserr = xnfgserr*0.0
       
    for idiv = 0L, ndiv-1L do begin

       up = where((xnfgs gt div[idiv]),nup,comp=lo,ncomp=nlo)
       if (nup ne 0L) then begin
          ynfgs[up] = nfgsnodust[indxnfgs[up]].zstrong_12oh_m91_upper
          ynfgserr[up] = nfgsnodust[indxnfgs[up]].zstrong_12oh_m91_upper_err
       endif

       if (nlo ne 0L) then begin
          ynfgs[lo] = nfgsnodust[indxnfgs[lo]].zstrong_12oh_m91_lower
          ynfgserr[lo] = nfgsnodust[indxnfgs[lo]].zstrong_12oh_m91_lower_err
       endif

       resid[idiv] = stddev(xnfgs-ynfgs)

    endfor

    minresid = min(resid,minindxnfgs)
    mindiv = div[minindxnfgs]
    splog, 'o3n2+n2 vs m91 [NFGS]: '+string(mindiv,format='(F4.2)')

    up = where((nfgsnodust[indxnfgs].zstrong_12oh_oiiinii_niiha gt mindiv),nup)
    if (nup ne 0L) then begin
       xnfgs = nfgsnodust[indxnfgs[up]].zstrong_12oh_oiiinii_niiha
       xnfgserr = nfgsnodust[indxnfgs[up]].zstrong_12oh_oiiinii_niiha_err

       ynfgs = nfgsnodust[indxnfgs[up]].zstrong_12oh_m91_upper
       ynfgserr = nfgsnodust[indxnfgs[up]].zstrong_12oh_m91_upper_err
    endif

    lo = where((nfgsnodust[indxnfgs].zstrong_12oh_oiiinii_niiha lt mindiv),nlo)
    if (nlo ne 0L) then begin
       xnfgs = [xnfgs,nfgsnodust[indxnfgs[lo]].zstrong_12oh_oiiinii_niiha]
       xnfgserr = [xnfgserr,nfgsnodust[indxnfgs[lo]].zstrong_12oh_oiiinii_niiha_err]

       ynfgs = [ynfgs,nfgsnodust[indxnfgs[lo]].zstrong_12oh_m91_lower]
       ynfgserr = [ynfgserr,nfgsnodust[indxnfgs[lo]].zstrong_12oh_m91_lower_err]
    endif

    xbig = [x,xnfgs]
    ybig = [y,ynfgs]

    ytitle = ohtitle+' R_{23}/M91'
;   ytitle = '12 + log (O/H)'

    stats = im_stats(x-y,/verbose,/no_head)
    xstr = strtrim(string(stats.median,format='(F12.2)'),2)+'\pm'+strtrim(string(stats.sigma_rej,format='(F12.2)'),2)
    
    atlas1d_lineplot, x, y, xerr, yerr, plottype=1, postscript=postscript, $
      xtitle='', ytitle=ytitle, xrange=xrange, yrange=yrange, legendtype=0, $
      xstyle=3, ystyle=3, /right, /top, position=pos[*,1], /noerase, charsize=charsize_4, $
      xtickname=replicate(' ',10), xnfgs=xnfgs, ynfgs=ynfgs, $
      xerrnfgs=xerrnfgs, yerrnfgs=yerrnfgs
    djs_oplot, !x.crange, !y.crange, thick=postthick

    legend, textoidl(xstr), /right, /bottom, box=0, charsize=charsize_4, charthick=postthick
;   label = ['(b) '+textoidl('O/H R_{23}/M91')]
;   legend, label, charthick=postthick, charsize=charsize_4, $
;     box=0, /left, /top
    legend, '(b)', /left, /top, box=0, charsize=charsize_4, charthick=postthick
    
; --------------------------------------------------    
; KD02-Combined
; --------------------------------------------------    

    indx = where((atlasnodust.ZSTRONG_12oh_oiiinii_niiha gt -900.0) and (atlasnodust.zstrong_12oh_kd02_combined gt -900.0),nindx)
 
    x = atlasnodust[indx].zstrong_12oh_oiiinii_niiha
    xerr = atlasnodust[indx].zstrong_12oh_oiiinii_niiha_err

    y = atlasnodust[indx].zstrong_12oh_kd02_combined
    yerr = atlasnodust[indx].zstrong_12oh_kd02_combined_err

    indxnfgs = where((nfgsnodust.ZSTRONG_12oh_oiiinii_niiha gt -900.0) and $
      (nfgsnodust.zstrong_12oh_kd02_combined gt -900.0),nindxnfgs)
    
    xnfgs = nfgsnodust[indxnfgs].zstrong_12oh_oiiinii_niiha
    xerrnfgs = nfgsnodust[indxnfgs].zstrong_12oh_oiiinii_niiha_err

    ynfgs = nfgsnodust[indxnfgs].zstrong_12oh_kd02_combined
    yerrnfgs = nfgsnodust[indxnfgs].zstrong_12oh_kd02_combined_err

    xbig = [x,xnfgs]
    ybig = [y,ynfgs]
       
    ytitle = ohtitle+' KD02'
;   ytitle = '12 + log (O/H)'

    stats = im_stats(x-y,/verbose,/no_head)
    xstr = strtrim(string(stats.median,format='(F12.2)'),2)+'\pm'+strtrim(string(stats.sigma_rej,format='(F12.2)'),2)
    
    atlas1d_lineplot, x, y, xerr, yerr, plottype=1, postscript=postscript, $
      xtitle=xtitle, ytitle=ytitle, xrange=xrange, yrange=yrange, legendtype=0, charsize=charsize_4, $
      xstyle=3, ystyle=3, /right, /top, position=pos[*,2], /noerase, xnfgs=xnfgs, ynfgs=ynfgs, $
      xerrnfgs=xerrnfgs, yerrnfgs=yerrnfgs
    djs_oplot, !x.crange, !y.crange, thick=postthick
 
    legend, textoidl(xstr), /right, /bottom, box=0, charsize=charsize_4, charthick=postthick
;   label = ['(c) '+textoidl('O/H KD02')]
;   legend, label, charthick=postthick, charsize=charsize_4, $
;     box=0, /left, /top
    legend, '(c)', /left, /top, box=0, charsize=charsize_4, charthick=postthick
    
; the title

;   xyouts, pos[0,0]*0.25, 0.5*(pos[3,0]-pos[1,2])+pos[1,2], textoidl(ytitle), $
;     orientation=90, align=0.5, charsize=charsize_4, charthick=postthick, /normal
;   xyouts, 0.5*(pos[2,1]-pos[0,0])+pos[0,0], pos[1,2]*0.15, textoidl(xtitle), align=0.5, $
;     charsize=charsize_4, charthick=postthick, /normal

    im_openclose, postscript=postscript, /close    

; ------------------------------------------------------------
; inter-compare empirical and strong-line calibrations
; ------------------------------------------------------------

    psname = '12oh_strong_vs_12oh_empirical'
    im_openclose, pspath+psname, postscript=postscript, /encapsulated

    pagemaker, nx=3, ny=3, position=pos, /normal, xspace=0.0, $
      xmargin=[1.0,0.1], ymargin=[0.2,1.1], yspace=0.0

    xrange = ohrange
    yrange = xrange

; --------------------------------------------------    
; Panel 1 - o3n2+n2 versus KD02-Combined
; --------------------------------------------------    

    indx = where((atlasnodust.ZSTRONG_12oh_oiiinii_niiha gt -900.0) and (atlasnodust.zstrong_12oh_kd02_combined gt -900.0),nindx)
 
    x = atlasnodust[indx].zstrong_12oh_oiiinii_niiha
    xerr = atlasnodust[indx].zstrong_12oh_oiiinii_niiha_err

    y = atlasnodust[indx].zstrong_12oh_kd02_combined
    yerr = atlasnodust[indx].zstrong_12oh_kd02_combined_err

    xtitle = ohtitle+' ([O III]/H\beta)/([N II]/H\alpha)'
    ytitle = ohtitle+' KD02'

    stats = im_stats(x-y,/verbose)
    xstr = strtrim(string(stats.median,format='(F12.2)'),2)+'\pm'+strtrim(string(stats.sigma_rej,format='(F12.2)'),2)
    
    atlas1d_lineplot, x, y, xerr, yerr, plottype=1, postscript=postscript, $
      xtitle='', ytitle=ytitle, xrange=xrange, yrange=yrange, legendtype=0, charsize=charsize_0, $
      xstyle=3, ystyle=3, /right, /top, position=pos[*,0], xtickname=replicate(' ',10)
    djs_oplot, !x.crange, !y.crange, thick=postthick
 
    legend, textoidl(xstr), /right, /bottom, box=0, charsize=charsize_0, charthick=postthick
    legend, '(a)', /left, /top, box=0, charsize=charsize_0, charthick=postthick
    
; --------------------------------------------------    
; Panel 2 - n2O2 versus KD02-Combined
; --------------------------------------------------    

    indx = where((atlasnodust.ZSTRONG_12oh_kd02_niioii gt -900.0) and (atlasnodust.zstrong_12oh_kd02_combined gt -900.0),nindx)

    x = atlasnodust[indx].zstrong_12oh_kd02_niioii
    xerr = atlasnodust[indx].zstrong_12oh_kd02_niioii_err

    y = atlasnodust[indx].zstrong_12oh_kd02_combined
    yerr = atlasnodust[indx].zstrong_12oh_kd02_combined_err

    xtitle = ohtitle+' [N II]/[O II]'
    ytitle = ohtitle+' KD02'

    stats = im_stats(x-y,/verbose,/no_head)
    xstr = strtrim(string(stats.median,format='(F12.2)'),2)+'\pm'+strtrim(string(stats.sigma_rej,format='(F12.2)'),2)
    
    atlas1d_lineplot, x, y, xerr, yerr, plottype=1, postscript=postscript, $
      xtitle='', ytitle='', xrange=xrange, yrange=yrange, legendtype=0, charsize=charsize_0, $
      xstyle=3, ystyle=3, /right, /top, position=pos[*,1], /noerase, ytickname=replicate(' ',10), $
      xtickname=replicate(' ',10)
    djs_oplot, !x.crange, !y.crange, thick=postthick

    legend, textoidl(xstr), /right, /bottom, box=0, charsize=charsize_0, charthick=postthick
    legend, '(b)', /left, /top, box=0, charsize=charsize_0, charthick=postthick

; --------------------------------------------------    
; Panel 3 - m91 versus KD02-Combined
; --------------------------------------------------    
 
    indx = where((atlasnodust.ZSTRONG_12oh_m91_upper gt -900.0) and (atlasnodust.ZSTRONG_12oh_m91_lower gt -900.0) and $
      (atlasnodust.zstrong_12oh_kd02_combined gt -900.0),nindx)

    y = atlasnodust[indx].zstrong_12oh_kd02_combined
    yerr = atlasnodust[indx].zstrong_12oh_kd02_combined_err

    x = y*0.0
    xerr = yerr*0.0
    
    for idiv = 0L, ndiv-1L do begin

       up = where((atlasnodust[indx].zstrong_12oh_kd02_combined gt div[idiv]),nup,comp=lo,ncomp=nlo)
       if (nup ne 0L) then begin
          x[up] = atlasnodust[indx[up]].zstrong_12oh_m91_upper
          xerr[up] = atlasnodust[indx[up]].zstrong_12oh_m91_upper_err
       endif

       if (nlo ne 0L) then begin
          x[lo] = atlasnodust[indx[lo]].zstrong_12oh_m91_lower
          xerr[lo] = atlasnodust[indx[lo]].zstrong_12oh_m91_lower_err
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
    splog, 'm91 vs KD02-Combined: '+string(mindiv,format='(F4.2)')

    up = where((atlasnodust[indx].zstrong_12oh_kd02_combined gt mindiv),nup)
    if (nup ne 0L) then begin
       x = atlasnodust[indx[up]].zstrong_12oh_m91_upper
       xerr = atlasnodust[indx[up]].zstrong_12oh_m91_upper_err

       y = atlasnodust[indx[up]].zstrong_12oh_kd02_combined
       yerr = atlasnodust[indx[up]].zstrong_12oh_kd02_combined_err
    endif

    lo = where((atlasnodust[indx].zstrong_12oh_kd02_combined lt mindiv),nlo)
    if (nlo ne 0L) then begin
       x = [x,atlasnodust[indx[lo]].zstrong_12oh_m91_lower]
       xerr = [xerr,atlasnodust[indx[lo]].zstrong_12oh_m91_lower_err]

       y = [y,atlasnodust[indx[lo]].zstrong_12oh_kd02_combined]
       yerr = [yerr,atlasnodust[indx[lo]].zstrong_12oh_kd02_combined_err]
    endif

    xtitle = ohtitle+' R_{23}/M91'
    ytitle = ohtitle+' KD02'

    stats = im_stats(x-y,/verbose,/no_head)
    xstr = strtrim(string(stats.median,format='(F12.2)'),2)+'\pm'+strtrim(string(stats.sigma_rej,format='(F12.2)'),2)
    
    atlas1d_lineplot, x, y, xerr, yerr, plottype=1, postscript=postscript, $
      xtitle=xtitle, ytitle='', xrange=xrange, yrange=yrange, legendtype=0, $
      xstyle=3, ystyle=3, /right, /top, position=pos[*,2], /noerase, charsize=charsize_0, $
      ytickname=replicate(' ',10)
    djs_oplot, !x.crange, !y.crange, thick=postthick
 
    legend, textoidl(xstr), /right, /bottom, box=0, charsize=charsize_0, charthick=postthick
    legend, '(c)', /left, /top, box=0, charsize=charsize_0, charthick=postthick

; --------------------------------------------------    
; Panel 4 - o3n2+n2 versus m91
; --------------------------------------------------    

    indx = where((atlasnodust.zstrong_12oh_oiiinii_niiha gt -900) and (atlasnodust.ZSTRONG_12oh_m91_upper gt -900.0) and $
      (atlasnodust.ZSTRONG_12oh_m91_lower gt -900.0),nindx)

    x = atlasnodust[indx].zstrong_12oh_oiiinii_niiha
    xerr = atlasnodust[indx].zstrong_12oh_oiiinii_niiha_err
    
    y = x*0.0
    yerr = xerr*0.0
    
    for idiv = 0L, ndiv-1L do begin

       up = where((x gt div[idiv]),nup,comp=lo,ncomp=nlo)
       if (nup ne 0L) then begin
          y[up] = atlasnodust[indx[up]].zstrong_12oh_m91_upper
          yerr[up] = atlasnodust[indx[up]].zstrong_12oh_m91_upper_err
       endif

       if (nlo ne 0L) then begin
          y[lo] = atlasnodust[indx[lo]].zstrong_12oh_m91_lower
          yerr[lo] = atlasnodust[indx[lo]].zstrong_12oh_m91_lower_err
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
    splog, 'o3n2+n2 vs m91: '+string(mindiv,format='(F4.2)')

    up = where((atlasnodust[indx].zstrong_12oh_oiiinii_niiha gt mindiv),nup)
    if (nup ne 0L) then begin
       x = atlasnodust[indx[up]].zstrong_12oh_oiiinii_niiha
       xerr = atlasnodust[indx[up]].zstrong_12oh_oiiinii_niiha_err

       y = atlasnodust[indx[up]].zstrong_12oh_m91_upper
       yerr = atlasnodust[indx[up]].zstrong_12oh_m91_upper_err
    endif

    lo = where((atlasnodust[indx].zstrong_12oh_oiiinii_niiha lt mindiv),nlo)
    if (nlo ne 0L) then begin
       x = [x,atlasnodust[indx[lo]].zstrong_12oh_oiiinii_niiha]
       xerr = [xerr,atlasnodust[indx[lo]].zstrong_12oh_oiiinii_niiha_err]

       y = [y,atlasnodust[indx[lo]].zstrong_12oh_m91_lower]
       yerr = [yerr,atlasnodust[indx[lo]].zstrong_12oh_m91_lower_err]
    endif

    xtitle = ohtitle+' ([O III]/H\beta)/([N II]/H\alpha)'
    ytitle = ohtitle+' R_{23}/M91'

    stats = im_stats(x-y,/verbose,/no_head)
    xstr = strtrim(string(stats.median,format='(F12.2)'),2)+'\pm'+strtrim(string(stats.sigma_rej,format='(F12.2)'),2)
    
    atlas1d_lineplot, x, y, xerr, yerr, plottype=1, postscript=postscript, $
      xtitle='', ytitle=ytitle, xrange=xrange, yrange=yrange, legendtype=0, $
      xstyle=3, ystyle=3, /right, /top, position=pos[*,3], /noerase, charsize=charsize_0, $
      xtickname=replicate(' ',10)
    djs_oplot, !x.crange, !y.crange, thick=postthick

    legend, textoidl(xstr), /right, /bottom, box=0, charsize=charsize_0, charthick=postthick
    legend, '(d)', /left, /top, box=0, charsize=charsize_0, charthick=postthick
    
; --------------------------------------------------    
; Panel 5 - n2O2 versus m91
; --------------------------------------------------    

    indx = where((atlasnodust.ZSTRONG_12oh_kd02_niioii gt -900.0) and (atlasnodust.zstrong_12oh_m91_upper gt -900.0) and $
      (atlasnodust.zstrong_12oh_m91_lower gt -900.0),nindx)
    
    x = atlasnodust[indx].zstrong_12oh_kd02_niioii
    xerr = atlasnodust[indx].zstrong_12oh_kd02_niioii_err
    
    y = x*0.0
    yerr = xerr*0.0
    
    for idiv = 0L, ndiv-1L do begin

       up = where((x gt div[idiv]),nup,comp=lo,ncomp=nlo)
       if (nup ne 0L) then begin
          y[up] = atlasnodust[indx[up]].zstrong_12oh_m91_upper
          yerr[up] = atlasnodust[indx[up]].zstrong_12oh_m91_upper_err
       endif

       if (nlo ne 0L) then begin
          y[lo] = atlasnodust[indx[lo]].zstrong_12oh_m91_lower
          yerr[lo] = atlasnodust[indx[lo]].zstrong_12oh_m91_lower_err
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
    splog, 'n2O2 vs m91: '+string(mindiv,format='(F4.2)')

    up = where((atlasnodust[indx].zstrong_12oh_kd02_niioii gt mindiv),nup)
    if (nup ne 0L) then begin
       x = atlasnodust[indx[up]].zstrong_12oh_kd02_niioii
       xerr = atlasnodust[indx[up]].zstrong_12oh_kd02_niioii_err

       y = atlasnodust[indx[up]].zstrong_12oh_m91_upper
       yerr = atlasnodust[indx[up]].zstrong_12oh_m91_upper_err
    endif

    lo = where((atlasnodust[indx].zstrong_12oh_kd02_niioii lt mindiv),nlo)
    if (nlo ne 0L) then begin
       x = [x,atlasnodust[indx[lo]].zstrong_12oh_kd02_niioii]
       xerr = [xerr,atlasnodust[indx[lo]].zstrong_12oh_kd02_niioii_err]

       y = [y,atlasnodust[indx[lo]].zstrong_12oh_m91_lower]
       yerr = [yerr,atlasnodust[indx[lo]].zstrong_12oh_m91_lower_err]
    endif

    xtitle = ohtitle+' [N II]/[O II]'
    ytitle = ohtitle+' R_{23}/M91'

    stats = im_stats(x-y,/verbose,/no_head)
    xstr = strtrim(string(stats.median,format='(F12.2)'),2)+'\pm'+strtrim(string(stats.sigma_rej,format='(F12.2)'),2)
    
    atlas1d_lineplot, x, y, xerr, yerr, plottype=1, postscript=postscript, $
      xtitle=xtitle, ytitle='', xrange=xrange, yrange=yrange, legendtype=0, charsize=charsize_0, $
      xstyle=3, ystyle=3, /right, /top, position=pos[*,4], /noerase, $
      ytickname=replicate(' ',10)
    djs_oplot, !x.crange, !y.crange, thick=postthick

    legend, textoidl(xstr), /right, /bottom, box=0, charsize=charsize_0, charthick=postthick
    legend, '(e)', /left, /top, box=0, charsize=charsize_0, charthick=postthick
    
; --------------------------------------------------    
; Panel 6 - No Data
; --------------------------------------------------    

; --------------------------------------------------    
; Panel 7 - o3n2+n2 versus n2O2
; --------------------------------------------------    

    indx = where((atlasnodust.zstrong_12oh_oiiinii_niiha gt -900) and (atlasnodust.ZSTRONG_12oh_kd02_niioii gt -900.0),nindx)

    x = atlasnodust[indx].ZSTRONG_12oh_oiiinii_niiha
    xerr = atlasnodust[indx].ZSTRONG_12oh_oiiinii_niiha_err

    y = atlasnodust[indx].zstrong_12oh_kd02_niioii
    yerr = atlasnodust[indx].zstrong_12oh_kd02_niioii_err

    xtitle = ohtitle+' ([O III]/H\beta)/([N II]/H\alpha)'
    ytitle = ohtitle+' [N II]/[O II]'

    stats = im_stats(x-y,/verbose,/no_head)
    xstr = strtrim(string(stats.median,format='(F12.2)'),2)+'\pm'+strtrim(string(stats.sigma_rej,format='(F12.2)'),2)
    
    atlas1d_lineplot, x, y, xerr, yerr, plottype=1, postscript=postscript, $
      xtitle=xtitle, ytitle=ytitle, xrange=xrange, yrange=yrange, legendtype=0, $
      xstyle=3, ystyle=3, /right, /top, position=pos[*,6], /noerase, charsize=charsize_0
    djs_oplot, !x.crange, !y.crange, thick=postthick

    legend, textoidl(xstr), /right, /bottom, box=0, charsize=charsize_0, charthick=postthick
    legend, '(f)', /left, /top, box=0, charsize=charsize_0, charthick=postthick
    
; --------------------------------------------------    
; Panel 8 - No Data
; --------------------------------------------------    

; --------------------------------------------------    
; Panel 9 - No Data
; --------------------------------------------------    

    im_openclose, postscript=postscript, /close    

; ------------------------------------------------------------
; R23 versus 12+log(O/H) [N II]/Ha - Integrated + SDSS
; ------------------------------------------------------------

    psname = 'r23_vs_12oh_niiha'
    im_openclose, pspath+psname, postscript=postscript, xsize=8.5, ysize=4.5, /encapsulated

    pagemaker, nx=2, ny=1, height=3.5, width=[3.4,3.4], xmargin=[1.1,0.3], $
      ymargin=[0.2,0.9], xspace=0.0, yspace=0, xpage=8.5, ypage=4.5, $
      position=pos, /normal

; Integrated + HII Regions    
    
    indx = where((atlasnodust.zstrong_12oh_niiha_pettini gt -900) and (atlasnodust.zstrong_r23 gt -900.0),nindx)

    x = atlasnodust[indx].zstrong_r23
    xerr = atlasnodust[indx].zstrong_r23_err

    y = atlasnodust[indx].zstrong_12oh_niiha_pettini
    yerr = atlasnodust[indx].zstrong_12oh_niiha_pettini_err

    indxnfgs = where((nfgsnodust.zstrong_12oh_niiha_pettini gt -900) and (nfgsnodust.zstrong_r23 gt -900.0),nindxnfgs)

    xnfgs = nfgsnodust[indxnfgs].zstrong_r23
    xerrnfgs = nfgsnodust[indxnfgs].zstrong_r23_err

    ynfgs = nfgsnodust[indxnfgs].zstrong_12oh_niiha_pettini
    yerrnfgs = nfgsnodust[indxnfgs].zstrong_12oh_niiha_pettini_err

    xtitle = 'log R_{23}'
    ytitle = ohtitle+' [N II]/H\alpha'

    xrange = r23range
    yrange = ohrange

    good = where((hii.ZSTRONG_r23 gt -900.0) and (hii.zstrong_12oh_niiha_pettini gt -900))
    xregion = hii[good].ZSTRONG_r23 & xerrregion = hii[good].ZSTRONG_r23
    yregion = hii[good].zstrong_12oh_niiha_pettini & yerrregion = hii[good].zstrong_12oh_niiha_pettini_err

    xbig = [x,xnfgs,xregion]
    ybig = [y,ynfgs,yregion]

    atlas1d_lineplot, x, y, xerr, yerr, plottype=3, postscript=postscript, $
      xtitle=xtitle, ytitle=ytitle, xrange=xrange, yrange=yrange, legendtype=0, $
      /right, /top, position=pos[*,0], charsize=charsize_6, $
      xregion=xregion, yregion=yregion, xerrregion=xerrregion, yerrregion=yerrregion, $
      xnfgs=xnfgs, ynfgs=ynfgs, xerrnfgs=xerrnfgs, yerrnfgs=yerrnfgs

; SDSS
    
    indxsdss = where((sdssnodust.zstrong_12oh_niiha_pettini gt -900) and (sdssnodust.zstrong_r23 gt -900.0),nindxsdss)

    xsdss = sdssnodust[indxsdss].zstrong_r23
    xerrsdss = sdssnodust[indxsdss].zstrong_r23_err

    ysdss = sdssnodust[indxsdss].zstrong_12oh_niiha_pettini
    yerrsdss = sdssnodust[indxsdss].zstrong_12oh_niiha_pettini_err

    sdss_lineplot, xsdss, ysdss, xerrsdss, yerrsdss, plottype=3, postscript=postscript, $
      xtitle=xtitle, ytitle='', xrange=xrange, yrange=yrange, legendtype=0, $
      /right, /top, position=pos[*,1], ytickname=replicate(' ',10), /noerase, $
      charsize=charsize_6

    im_openclose, postscript=postscript, /close    

; ------------------------------------------------------------
; P versus 12+log(O/H) [N II]/Ha - Integrated + SDSS
; ------------------------------------------------------------

    psname = 'p_vs_12oh_niiha'
    im_openclose, pspath+psname, postscript=postscript, xsize=8.5, ysize=4.5, /encapsulated

    pagemaker, nx=2, ny=1, height=3.5, width=[3.4,3.4], xmargin=[1.1,0.3], $
      ymargin=[0.2,0.9], xspace=0.0, yspace=0, xpage=8.5, ypage=4.5, $
      position=pos, /normal

; Integrated + HII Regions    
    
    indx = where((atlasnodust.zstrong_12oh_niiha_pettini gt -900) and (atlasnodust.zstrong_p gt -900.0),nindx)

    x = atlasnodust[indx].zstrong_p
    xerr = atlasnodust[indx].zstrong_p_err

    y = atlasnodust[indx].zstrong_12oh_niiha_pettini
    yerr = atlasnodust[indx].zstrong_12oh_niiha_pettini_err

    indxnfgs = where((nfgsnodust.zstrong_12oh_niiha_pettini gt -900) and (nfgsnodust.zstrong_p gt -900.0),nindxnfgs)

    xnfgs = nfgsnodust[indxnfgs].zstrong_p
    xerrnfgs = nfgsnodust[indxnfgs].zstrong_p_err

    ynfgs = nfgsnodust[indxnfgs].zstrong_12oh_niiha_pettini
    yerrnfgs = nfgsnodust[indxnfgs].zstrong_12oh_niiha_pettini_err

    xtitle = 'log P'
    ytitle = ohtitle+' [N II]/H\alpha'

    xrange = prange
    yrange = ohrange

    good = where((hii.ZSTRONG_p gt -900.0) and (hii.zstrong_12oh_niiha_pettini gt -900))
    xregion = hii[good].ZSTRONG_p & xerrregion = hii[good].ZSTRONG_p
    yregion = hii[good].zstrong_12oh_niiha_pettini & yerrregion = hii[good].zstrong_12oh_niiha_pettini_err

    xbig = [x,xnfgs,xregion]
    ybig = [y,ynfgs,yregion]

    atlas1d_lineplot, x, y, xerr, yerr, plottype=3, postscript=postscript, $
      xtitle=xtitle, ytitle=ytitle, xrange=xrange, yrange=yrange, legendtype=0, $
      /right, /top, position=pos[*,0], charsize=charsize_6, $
      xregion=xregion, yregion=yregion, xerrregion=xerrregion, yerrregion=yerrregion, $
      xnfgs=xnfgs, ynfgs=ynfgs, xerrnfgs=xerrnfgs, yerrnfgs=yerrnfgs

; fit a bisector    
    
;   inrange = where((ybig gt 7.3) and (ybig lt 8.9)) ; range of the Denicolo calibration
;   
;   sixlin, xbig[inrange], ybig[inrange], a, siga, b, sigb
;   ii = 2                      ; Ordinary Least Squares Bisector
;   coeff = [a[ii],b[ii]]
;   coeff_err = [siga[2],sigb[2]]
;
;   axis = findgen(((1.0)-(0.0))/0.01)*0.01+(0.0)
;   yfit = poly(axis,coeff)    
;   keep = where((yfit ge 7.3) and (yfit le 8.9))
;   axis = axis[keep] & yfit = yfit[keep]
;
;   djs_oplot, axis, yfit, thick=postthick, line=0
;   
;   splog, 'log P vs 12+log(O/H) [n2] coefficients: '
;   niceprint, coeff, coeff_err, minmax(yfit)
;
;   OHup = 8.5
;   OHlo = 8.2
;   
;   oplot, !x.crange, OHup*[1,1], line=1, thick=postthick
;   oplot, !x.crange, OHlo*[1,1], line=1, thick=postthick
;
;   xyouts, 0.0, 7.7, 'Lower Branch', align=0.0, /data, charsize=2.0, $
;     charthick=postthick
;   xyouts, 1.0, 9.5, 'Upper Branch', align=1.0, /data, charsize=2.0, $
;     charthick=postthick

; SDSS
    
    indxsdss = where((sdssnodust.zstrong_12oh_niiha_pettini gt -900) and (sdssnodust.zstrong_p gt -900.0),nindxsdss)

    xsdss = sdssnodust[indxsdss].zstrong_p
    xerrsdss = sdssnodust[indxsdss].zstrong_p_err

    ysdss = sdssnodust[indxsdss].zstrong_12oh_niiha_pettini
    yerrsdss = sdssnodust[indxsdss].zstrong_12oh_niiha_pettini_err

    sdss_lineplot, xsdss, ysdss, xerrsdss, yerrsdss, plottype=3, postscript=postscript, $
      xtitle=xtitle, ytitle='', xrange=xrange, yrange=yrange, legendtype=0, $
      /right, /top, position=pos[*,1], ytickname=replicate(' ',10), /noerase, $
      charsize=charsize_6

    im_openclose, postscript=postscript, /close    

; ------------------------------------------------------------
; 12+log(O/H) R23/M91_O32 vs EW{12+log(O/H) R23/M91_O32}
; ------------------------------------------------------------

    psname = '12oh_m91o32_vs_12oh_ewm91o32'
    im_openclose, pspath+psname, postscript=postscript, xsize=8.5, ysize=4.5, /encapsulated

    pagemaker, nx=2, ny=1, height=3.5, width=[3.4,3.4], xmargin=[1.1,0.3], $
      ymargin=[0.2,0.9], xspace=0.0, yspace=0, xpage=8.5, ypage=4.5, $
      position=pos, /normal

; raw line fluxes    
    
    indx = where((atlasdust.zstrong_12oh_m91_o32 gt -900) and (atlasdust.zstrong_ew_12oh_m91_o32 gt -900.0),nindx)

    x = atlasdust[indx].zstrong_12oh_m91_o32
    xerr = atlasdust[indx].zstrong_12oh_m91_o32_err

    y = atlasdust[indx].zstrong_ew_12oh_m91_o32
    yerr = atlasdust[indx].zstrong_ew_12oh_m91_o32_err

    indxnfgs = where((nfgsdust.zstrong_12oh_m91_o32 gt -900) and (nfgsdust.zstrong_ew_12oh_m91_o32 gt -900.0),nindxnfgs)

    xnfgs = nfgsdust[indxnfgs].zstrong_12oh_m91_o32
    xerrnfgs = nfgsdust[indxnfgs].zstrong_12oh_m91_o32_err

    ynfgs = nfgsdust[indxnfgs].zstrong_ew_12oh_m91_o32
    yerrnfgs = nfgsdust[indxnfgs].zstrong_ew_12oh_m91_o32_err

    residuals = [y,ynfgs]-[x,xnfgs]
    stats = im_stats(residuals,/verbose)
    xstr = strtrim(string(stats.median,format='(F12.2)'),2)+'\pm'+strtrim(string(stats.sig68mean,format='(F12.2)'),2)

    xtitle = ohtitle+' [Observed]'
    ytitle = 'EW['+ohtitle+']'

    xrange = ohrange2
    yrange = xrange

    atlas1d_lineplot, x, y, xerr, yerr, plottype=3, postscript=postscript, $
      xtitle=xtitle, ytitle=ytitle, xrange=xrange, yrange=yrange, legendtype=0, $
      /right, /top, position=pos[*,0], charsize=charsize_6, $
      xnfgs=xnfgs, ynfgs=ynfgs, xerrnfgs=xerrnfgs, yerrnfgs=yerrnfgs
    djs_oplot, !x.crange, !y.crange, line=0, thick=postthick
    legend, textoidl(xstr), /left, /top, box=0, charsize=charsize_6, charthick=postthick

; reddening-corrected fluxes
    
    indx = where((atlasnodust.zstrong_12oh_m91_o32 gt -900) and (atlasdust.zstrong_ew_12oh_m91_o32 gt -900.0),nindx)

    x = atlasnodust[indx].zstrong_12oh_m91_o32
    xerr = atlasnodust[indx].zstrong_12oh_m91_o32_err

    y = atlasdust[indx].zstrong_ew_12oh_m91_o32
    yerr = atlasdust[indx].zstrong_ew_12oh_m91_o32_err

    indxnfgs = where((nfgsnodust.zstrong_12oh_m91_o32 gt -900) and (nfgsdust.zstrong_ew_12oh_m91_o32 gt -900.0),nindxnfgs)

    xnfgs = nfgsnodust[indxnfgs].zstrong_12oh_m91_o32
    xerrnfgs = nfgsnodust[indxnfgs].zstrong_12oh_m91_o32_err

    ynfgs = nfgsdust[indxnfgs].zstrong_ew_12oh_m91_o32
    yerrnfgs = nfgsdust[indxnfgs].zstrong_ew_12oh_m91_o32_err

    residuals = [y,ynfgs]-[x,xnfgs]
    stats = im_stats(residuals,/verbose)
    xstr = strtrim(string(stats.median,format='(F12.2)'),2)+'\pm'+strtrim(string(stats.sig68mean,format='(F12.2)'),2)

    xtitle = ohtitle+' [Corrected]'

    atlas1d_lineplot, x, y, xerr, yerr, plottype=3, postscript=postscript, $
      xtitle=xtitle, ytitle='', xrange=xrange, yrange=yrange, legendtype=0, $
      /right, /top, position=pos[*,1], /noerase, ytickname=replicate(' ',10), charsize=charsize_6, $
      xnfgs=xnfgs, ynfgs=ynfgs, xerrnfgs=xerrnfgs, yerrnfgs=yerrnfgs
    djs_oplot, !x.crange, !y.crange, line=0, thick=postthick
    legend, textoidl(xstr), /left, /top, box=0, charsize=charsize_6, charthick=postthick

    im_openclose, postscript=postscript, /close    

; ------------------------------------------------------------
; R23 vs EWR23 - Integrated
; ------------------------------------------------------------

    psname = 'r23_vs_ewr23'
    im_openclose, pspath+psname, postscript=postscript, xsize=8.5, ysize=4.5, /encapsulated

    pagemaker, nx=2, ny=1, height=3.5, width=[3.4,3.4], xmargin=[1.1,0.3], $
      ymargin=[0.2,0.9], xspace=0.0, yspace=0, xpage=8.5, ypage=4.5, $
      position=pos, /normal

; raw line fluxes    
    
    indx = where((atlasdust.zstrong_r23 gt -900) and (atlasdust.zstrong_ew_r23 gt -900.0),nindx)

    x = atlasdust[indx].zstrong_r23
    xerr = atlasdust[indx].zstrong_r23_err

    y = atlasdust[indx].zstrong_ew_r23
    yerr = atlasdust[indx].zstrong_ew_r23_err

    indxnfgs = where((nfgsdust.zstrong_r23 gt -900) and (nfgsdust.zstrong_ew_r23 gt -900.0),nindxnfgs)

    xnfgs = nfgsdust[indxnfgs].zstrong_r23
    xerrnfgs = nfgsdust[indxnfgs].zstrong_r23_err

    ynfgs = nfgsdust[indxnfgs].zstrong_ew_r23
    yerrnfgs = nfgsdust[indxnfgs].zstrong_ew_r23_err

    residuals = [y,ynfgs]-[x,xnfgs]
    stats = im_stats(residuals,/verbose)
    xstr = strtrim(string(stats.median,format='(F12.2)'),2)+'\pm'+strtrim(string(stats.sig68mean,format='(F12.2)'),2)

;   xtitle = 'log (R_{23})_{obs}'
    xtitle = 'log (R_{23}) [Observed]'
    ytitle = 'log (EWR_{23})'

    xrange = R23range
    yrange = R23range

    atlas1d_lineplot, x, y, xerr, yerr, plottype=3, postscript=postscript, $
      xtitle=xtitle, ytitle=ytitle, xrange=xrange, yrange=yrange, legendtype=0, $
      /right, /top, position=pos[*,0], charsize=charsize_6, $
      xnfgs=xnfgs, ynfgs=ynfgs, xerrnfgs=xerrnfgs, yerrnfgs=yerrnfgs
    djs_oplot, !x.crange, !y.crange, line=0, thick=postthick
    legend, textoidl(xstr), /left, /top, box=0, charsize=charsize_6, charthick=postthick

; reddening-corrected fluxes
    
    indx = where((atlasnodust.zstrong_r23 gt -900) and (atlasdust.zstrong_ew_r23 gt -900.0),nindx)

    x = atlasnodust[indx].zstrong_r23
    xerr = atlasnodust[indx].zstrong_r23_err

    y = atlasdust[indx].zstrong_ew_r23
    yerr = atlasdust[indx].zstrong_ew_r23_err

    indxnfgs = where((nfgsnodust.zstrong_r23 gt -900) and (nfgsdust.zstrong_ew_r23 gt -900.0),nindxnfgs)

    xnfgs = nfgsnodust[indxnfgs].zstrong_r23
    xerrnfgs = nfgsnodust[indxnfgs].zstrong_r23_err

    ynfgs = nfgsdust[indxnfgs].zstrong_ew_r23
    yerrnfgs = nfgsdust[indxnfgs].zstrong_ew_r23_err

    residuals = [y,ynfgs]-[x,xnfgs]
    stats = im_stats(residuals,/verbose)
    xstr = strtrim(string(stats.median,format='(F12.2)'),2)+'\pm'+strtrim(string(stats.sig68mean,format='(F12.2)'),2)

;   xtitle = 'log (R_{23})_{cor}'
    xtitle = 'log (R_{23}) [Corrected]'

    atlas1d_lineplot, x, y, xerr, yerr, plottype=3, postscript=postscript, $
      xtitle=xtitle, ytitle='', xrange=xrange, yrange=yrange, legendtype=0, $
      /right, /top, position=pos[*,1], /noerase, ytickname=replicate(' ',10), $
      charsize=charsize_6, xnfgs=xnfgs, ynfgs=ynfgs, xerrnfgs=xerrnfgs, yerrnfgs=yerrnfgs
    djs_oplot, !x.crange, !y.crange, line=0, thick=postthick
    legend, textoidl(xstr), /left, /top, box=0, charsize=charsize_6, charthick=postthick

    im_openclose, postscript=postscript, /close    

; ------------------------------------------------------------
; R23 vs EWR23 - SDSS
; ------------------------------------------------------------

    psname = 'sdss_r23_vs_ewr23'
    im_openclose, pspath+psname, postscript=postscript, xsize=8.5, ysize=4.5, /encapsulated

    pagemaker, nx=2, ny=1, height=3.5, width=[3.4,3.4], xmargin=[1.1,0.3], $
      ymargin=[0.2,0.9], xspace=0.0, yspace=0, xpage=8.5, ypage=4.5, $
      position=pos, /normal

; raw line fluxes    
    
    indx = where((sdssdust.zstrong_r23 gt -900) and (sdssdust.zstrong_ew_r23 gt -900.0),nindx)

    x = sdssdust[indx].zstrong_r23
    xerr = sdssdust[indx].zstrong_r23_err

    y = sdssdust[indx].zstrong_ew_r23
    yerr = sdssdust[indx].zstrong_ew_r23_err

    residuals = x-y
    stats = im_stats(residuals,/verbose)
    xstr = strtrim(string(stats.median,format='(F12.2)'),2)+'\pm'+strtrim(string(stats.sig68mean,format='(F12.2)'),2)

;   xtitle = 'log (R_{23})_{obs}'
    xtitle = 'log (R_{23}) [Observed]'
    ytitle = 'log (EWR_{23})'

    xrange = R23range
    yrange = R23range

    sdss_lineplot, x, y, xerr, yerr, plottype=3, postscript=postscript, $
      xtitle=xtitle, ytitle=ytitle, xrange=xrange, yrange=yrange, legendtype=0, $
      /right, /top, position=pos[*,0], charsize=charsize_6
    djs_oplot, !x.crange, !y.crange, line=0, thick=postthick
    legend, textoidl(xstr), /left, /top, box=0, charsize=charsize_6, charthick=postthick

; reddening-corrected fluxes
    
    indx = where((sdssnodust.zstrong_r23 gt -900) and (sdssdust.zstrong_ew_r23 gt -900.0),nindx)

    x = sdssnodust[indx].zstrong_r23
    xerr = sdssnodust[indx].zstrong_r23_err

    y = sdssdust[indx].zstrong_ew_r23
    yerr = sdssdust[indx].zstrong_ew_r23_err

    residuals = x-y
    stats = im_stats(residuals,/verbose)
    xstr = strtrim(string(stats.median,format='(F12.2)'),2)+'\pm'+strtrim(string(stats.sig68mean,format='(F12.2)'),2)

;   xtitle = 'log (R_{23})_{cor}'
    xtitle = 'log (R_{23}) [Corrected]'

    sdss_lineplot, x, y, xerr, yerr, plottype=3, postscript=postscript, $
      xtitle=xtitle, ytitle='', xrange=xrange, yrange=yrange, legendtype=0, $
      /right, /top, position=pos[*,1], /noerase, ytickname=replicate(' ',10), $
      charsize=charsize_6
    djs_oplot, !x.crange, !y.crange, line=0, thick=postthick
    legend, textoidl(xstr), /left, /top, box=0, charsize=charsize_6, charthick=postthick

    im_openclose, postscript=postscript, /close    

; ------------------------------------------------------------
; O32 vs EWO32 - Integrated
; ------------------------------------------------------------

    psname = 'o32_vs_ewo32'
    im_openclose, pspath+psname, postscript=postscript, xsize=8.5, ysize=4.5, /encapsulated

    pagemaker, nx=2, ny=1, height=3.5, width=[3.4,3.4], xmargin=[1.1,0.3], $
      ymargin=[0.2,0.9], xspace=0.0, yspace=0, xpage=8.5, ypage=4.5, $
      position=pos, /normal

; raw line fluxes    
    
    indx = where((atlasdust.zstrong_o32 gt -900) and (atlasdust.zstrong_ew_o32 gt -900.0),nindx)

    x = atlasdust[indx].zstrong_o32
    xerr = atlasdust[indx].zstrong_o32_err

    y = atlasdust[indx].zstrong_ew_o32
    yerr = atlasdust[indx].zstrong_ew_o32_err

    indxnfgs = where((nfgsdust.zstrong_o32 gt -900) and (nfgsdust.zstrong_ew_o32 gt -900.0),nindxnfgs)

    xnfgs = nfgsdust[indxnfgs].zstrong_o32
    xerrnfgs = nfgsdust[indxnfgs].zstrong_o32_err

    ynfgs = nfgsdust[indxnfgs].zstrong_ew_o32
    yerrnfgs = nfgsdust[indxnfgs].zstrong_ew_o32_err

    residuals = [y,ynfgs]-[x,xnfgs]
    stats = im_stats(residuals,/verbose)
    xstr = strtrim(string(stats.median,format='(F12.2)'),2)+'\pm'+strtrim(string(stats.sig68mean,format='(F12.2)'),2)

;   xtitle = 'log (O_{32})_{obs}'
    xtitle = 'log (O_{32}) [Observed]'
    ytitle = 'log (EWO_{32})'

    xrange = o32range2
    yrange = o32range2

    atlas1d_lineplot, x, y, xerr, yerr, plottype=3, postscript=postscript, $
      xtitle=xtitle, ytitle=ytitle, xrange=xrange, yrange=yrange, legendtype=0, $
      /right, /top, position=pos[*,0], charsize=charsize_6, $
      xnfgs=xnfgs, ynfgs=ynfgs, xerrnfgs=xerrnfgs, yerrnfgs=yerrnfgs
    djs_oplot, !x.crange, !y.crange, line=0, thick=postthick
    legend, textoidl(xstr), /left, /top, box=0, charsize=charsize_6, charthick=postthick

; reddening-corrected fluxes
    
    indx = where((atlasnodust.zstrong_o32 gt -900) and (atlasdust.zstrong_ew_o32 gt -900.0),nindx)

    x = atlasnodust[indx].zstrong_o32
    xerr = atlasnodust[indx].zstrong_o32_err

    y = atlasdust[indx].zstrong_ew_o32
    yerr = atlasdust[indx].zstrong_ew_o32_err

    indxnfgs = where((nfgsnodust.zstrong_o32 gt -900) and (nfgsdust.zstrong_ew_o32 gt -900.0),nindxnfgs)

    xnfgs = nfgsnodust[indxnfgs].zstrong_o32
    xerrnfgs = nfgsnodust[indxnfgs].zstrong_o32_err

    ynfgs = nfgsdust[indxnfgs].zstrong_ew_o32
    yerrnfgs = nfgsdust[indxnfgs].zstrong_ew_o32_err
    
    residuals = [y,ynfgs]-[x,xnfgs]
    stats = im_stats(residuals,/verbose)
    xstr = strtrim(string(stats.median,format='(F12.2)'),2)+'\pm'+strtrim(string(stats.sig68mean,format='(F12.2)'),2)

;   xtitle = 'log (O_{32})_{cor}'
    xtitle = 'log (O_{32}) [Corrected]'

    atlas1d_lineplot, x, y, xerr, yerr, plottype=3, postscript=postscript, $
      xtitle=xtitle, ytitle='', xrange=xrange, yrange=yrange, legendtype=0, $
      /right, /top, position=pos[*,1], /noerase, ytickname=replicate(' ',10), $
      charsize=charsize_6, xnfgs=xnfgs, ynfgs=ynfgs, xerrnfgs=xerrnfgs, yerrnfgs=yerrnfgs
    djs_oplot, !x.crange, !y.crange, line=0, thick=postthick
    legend, textoidl(xstr), /left, /top, box=0, charsize=charsize_6, charthick=postthick

    im_openclose, postscript=postscript, /close    

; ------------------------------------------------------------
; O32 vs EWO32 - SDSS
; ------------------------------------------------------------

    psname = 'sdss_o32_vs_ewo32'
    im_openclose, pspath+psname, postscript=postscript, xsize=8.5, ysize=4.5, /encapsulated

    pagemaker, nx=2, ny=1, height=3.5, width=[3.4,3.4], xmargin=[1.1,0.3], $
      ymargin=[0.2,0.9], xspace=0.0, yspace=0, xpage=8.5, ypage=4.5, $
      position=pos, /normal

; raw line fluxes    
    
    indx = where((sdssdust.zstrong_o32 gt -900) and (sdssdust.zstrong_ew_o32 gt -900.0),nindx)

    x = sdssdust[indx].zstrong_o32
    xerr = sdssdust[indx].zstrong_o32_err

    y = sdssdust[indx].zstrong_ew_o32
    yerr = sdssdust[indx].zstrong_ew_o32_err

    residuals = x-y
    stats = im_stats(residuals,/verbose)
    xstr = strtrim(string(stats.median,format='(F12.2)'),2)+'\pm'+strtrim(string(stats.sig68mean,format='(F12.2)'),2)

;   xtitle = 'log (O_{32})_{obs}'
    xtitle = 'log (O_{32}) [Observed]'
    ytitle = 'log (EWO_{32})'

    xrange = o32range
    yrange = o32range

    sdss_lineplot, x, y, xerr, yerr, plottype=3, postscript=postscript, $
      xtitle=xtitle, ytitle=ytitle, xrange=xrange, yrange=yrange, legendtype=0, $
      /right, /top, position=pos[*,0], charsize=charsize_6
    djs_oplot, !x.crange, !y.crange, line=0, thick=postthick
    legend, textoidl(xstr), /left, /top, box=0, charsize=charsize_6, charthick=postthick

; reddening-corrected fluxes
    
    indx = where((sdssnodust.zstrong_o32 gt -900) and (sdssdust.zstrong_ew_o32 gt -900.0),nindx)

    x = sdssnodust[indx].zstrong_o32
    xerr = sdssnodust[indx].zstrong_o32_err

    y = sdssdust[indx].zstrong_ew_o32
    yerr = sdssdust[indx].zstrong_ew_o32_err

    residuals = x-y
    stats = im_stats(residuals,/verbose)
    xstr = strtrim(string(stats.median,format='(F12.2)'),2)+'\pm'+strtrim(string(stats.sig68mean,format='(F12.2)'),2)

;   xtitle = 'log (O_{32})_{cor}'
    xtitle = 'log (O_{32}) [Corrected]'

    sdss_lineplot, x, y, xerr, yerr, plottype=3, postscript=postscript, $
      xtitle=xtitle, ytitle='', xrange=xrange, yrange=yrange, legendtype=0, $
      /right, /top, position=pos[*,1], /noerase, ytickname=replicate(' ',10), $
      charsize=charsize_6
    djs_oplot, !x.crange, !y.crange, line=0, thick=postthick
    legend, textoidl(xstr), /left, /top, box=0, charsize=charsize_6, charthick=postthick

    im_openclose, postscript=postscript, /close    

; --------------------------------------------------    
; GENERATE THE WEB PAGE
; --------------------------------------------------    
    
    if keyword_set(postscript) then im_ps2html, htmlbase, html_path=html_path, $
      cleanpng=0, npscols=3, _extra=extra

stop    
    
; --------------------------------------------------    
; SELECT PLOTS FOR THE PAPER HERE
; --------------------------------------------------    

    if keyword_set(paper) then begin

       splog, 'Writing paper plots to '+paperpath+'.'
       paperplots = [$
         'integrated_12oh_oiiinii_niiha_vs_12oh_strong',$
         'o32_vs_12oh_niiha'$
         ]
       for k = 0L, n_elements(paperplots)-1L do $
         spawn, ['/bin/cp -f '+html_path+htmlbase+'/'+paperplots[k]+' '+paperpath], /sh

    endif

; ###########################################################################    
; End Paper Plots
; ###########################################################################    

stop    
    
; ###########################################################################    
; Begin Older Plots
; ###########################################################################    

; ------------------------------------------------------------
; 3-panel 12+log(O/H) [T04] vs 12+log(O/H) Empirical, Moustakas et
; al. (2005) calibrations - residuals
; ------------------------------------------------------------

    psname = 'delta_12oh_T04_empirical_3panel'
    im_openclose, pspath+psname, postscript=postscript, encapsulated=encapsulated, xsize=8.5, ysize=10.5

    pagemaker, nx=1, ny=3, height=3.0*[1,1,1], width=6.8, xmargin=[1.3,0.4], $
      ymargin=[0.4,1.1], xspace=0, yspace=0, xpage=8.5, ypage=10.5, $
      position=pos, /normal

;   if keyword_set(dotalk) then polyfill, [1,1,0,0,1], [1,0,0,1,1], /normal, $
;     color=djs_icolor('black')

    xtitle = ohtitle+' [T04]'
    ytitle = '\Delta[log(O/H)]'

    xrange = ohrange2
    yrange = [-1.2,1.2]

; --------------------------------------------------    
; ([O III]/Hb)/([N II]/Ha)
; --------------------------------------------------    

    indx = where((sdssnodust.zstrong_12oh_oiiinii_moustakas gt -900) and (sdssancillary.tremonti_oh gt 0.0),nindx)

    x = sdssancillary[indx].tremonti_oh
    xerr = sdssancillary[indx].tremonti_oh_err

    y = sdssnodust[indx].zstrong_12oh_oiiinii_moustakas
    yerr = sdssnodust[indx].zstrong_12oh_oiiinii_moustakas_err

    residuals = x-y
    residuals_err = sqrt(xerr^2+yerr^2)

    stats = im_stats(residuals)
    xstr = strtrim(string(stats.median,format='(F12.2)'),2)+'\pm'+strtrim(string(stats.sig68mean,format='(F12.2)'),2)

    sdss_lineplot, x, residuals, xerr, residuals_err, plottype=1, postscript=postscript, $
      xtitle='', ytitle=ytitle, xrange=xrange, yrange=yrange, legendtype=0, $
      /right, /top, position=pos[*,0], charsize=charsize_4, xtickname=replicate(' ',10)
    djs_oplot, !x.crange, [0,0], thick=postthick, color=djs_icolor(talkcolor), line=2

    legend, textoidl(xstr), /right, /bottom, box=0, charsize=charsize_4, $
      charthick=postthick, textcolor=djs_icolor(talkcolor)
    legend, textoidl('(a) ([O III]/H\beta)/([N II]/H\alpha)'), /left, /top, box=0, $
      charsize=charsize_4, charthick=postthick, textcolor=djs_icolor(talkcolor)

; --------------------------------------------------    
; [N II]/Ha
; --------------------------------------------------    

    indx = where((sdssnodust.zstrong_12oh_niiha_moustakas gt -900) and (sdssancillary.tremonti_oh gt 0.0),nindx)

    x = sdssancillary[indx].tremonti_oh
    xerr = sdssancillary[indx].tremonti_oh_err

    y = sdssnodust[indx].zstrong_12oh_niiha_moustakas
    yerr = sdssnodust[indx].zstrong_12oh_niiha_moustakas_err

    residuals = x-y
    residuals_err = sqrt(xerr^2+yerr^2)

    stats = im_stats(residuals)
    xstr = strtrim(string(stats.median,format='(F12.2)'),2)+'\pm'+strtrim(string(stats.sig68mean,format='(F12.2)'),2)

    sdss_lineplot, x, residuals, xerr, residuals_err, plottype=1, postscript=postscript, $
      xtitle='', ytitle=ytitle, xrange=xrange, yrange=yrange, legendtype=0, /noerase, $
      /right, /top, position=pos[*,1], charsize=charsize_4, xtickname=replicate(' ',10)
    djs_oplot, !x.crange, [0,0], thick=postthick, color=djs_icolor(talkcolor), line=2

    legend, textoidl(xstr), /right, /bottom, box=0, charsize=charsize_4, $
      charthick=postthick, textcolor=djs_icolor(talkcolor)
    legend, textoidl('(b) [N II]/H\alpha'), /left, /top, box=0, $
      charsize=charsize_4, charthick=postthick, textcolor=djs_icolor(talkcolor)

; --------------------------------------------------    
; P-method
; --------------------------------------------------    

    indx = where((sdssancillary.tremonti_oh gt 0.0) and (sdssnodust.zstrong_12oh_p01_upper gt -900.0) and $
      (sdssnodust.zstrong_12oh_p01_lower gt -900.0),nindx)

    up = where((sdssancillary[indx].tremonti_oh gt 8.2),nup)
    if (nup ne 0L) then begin
       x = sdssancillary[indx[up]].tremonti_oh
       xerr = sdssancillary[indx[up]].tremonti_oh_err

       y = sdssnodust[indx[up]].zstrong_12oh_p01_upper
       yerr = sdssnodust[indx[up]].zstrong_12oh_p01_upper_err
    endif

    lo = where((sdssancillary[indx].tremonti_oh lt 7.95),nlo)
    if (nlo ne 0L) then begin
       x = [x,sdssancillary[indx[lo]].tremonti_oh]
       xerr = [xerr,[sdssnodust[indx[lo]].tremonti_p84-sdssnodust[indx[lo]].tremonti_p16]/2.0]

       y = [y,sdssnodust[indx[lo]].zstrong_12oh_p01_lower]
       yerr = [yerr,sdssnodust[indx[lo]].zstrong_12oh_p01_lower_err]
    endif

    residuals = x-y
    residuals_err = sqrt(xerr^2+yerr^2)

    stats = im_stats(residuals)
    xstr = strtrim(string(stats.median,format='(F12.2)'),2)+'\pm'+strtrim(string(stats.sig68mean,format='(F12.2)'),2)

    sdss_lineplot, x, residuals, xerr, residuals_err, plottype=1, postscript=postscript, $
      xtitle=xtitle, ytitle=ytitle, xrange=xrange, yrange=yrange, legendtype=0, $
      /right, /top, position=pos[*,2], charsize=charsize_4, /noerase
    djs_oplot, !x.crange, [0,0], thick=postthick, color=djs_icolor(talkcolor), line=2

    legend, textoidl(xstr), /right, /bottom, box=0, charsize=charsize_4, $
      charthick=postthick, textcolor=djs_icolor(talkcolor)
    legend, textoidl('(c) P-method'), /left, /top, box=0, $
      charsize=charsize_4, charthick=postthick, textcolor=djs_icolor(talkcolor)

    im_openclose, postscript=postscript, /close    
    
; ------------------------------------------------------------
; 3-panel 12+log(O/H) [T04] vs 12+log(O/H) Empirical, Pettini & Pagel
; (2004) calibrations 
; ------------------------------------------------------------

    psname = '12oh_T04_vs_12oh_empirical_pettini'
    im_openclose, pspath+psname, postscript=postscript, encapsulated=encapsulated, xsize=8.5, ysize=10.5

    pagemaker, nx=1, ny=3, height=3.0*[1,1,1], width=6.8, xmargin=[1.3,0.4], $
      ymargin=[0.4,1.1], xspace=0, yspace=0, xpage=8.5, ypage=10.5, $
      position=pos, /normal

    if keyword_set(dotalk) then polyfill, [1,1,0,0,1], [1,0,0,1,1], /normal, $
      color=djs_icolor('black')

    xtitle = ohtitle+' [T04]'

    xrange = ohrange2
    yrange = xrange

; --------------------------------------------------    
; ([O III]/Hb)/([N II]/Ha)
; --------------------------------------------------    

    indx = where((sdssnodust.zstrong_12oh_oiiinii_pettini gt -900) and (sdssancillary.tremonti_oh gt 0.0),nindx)

    x = sdssancillary[indx].tremonti_oh
    xerr = sdssancillary[indx].tremonti_oh_err

    y = sdssnodust[indx].zstrong_12oh_oiiinii_pettini
    yerr = sdssnodust[indx].zstrong_12oh_oiiinii_pettini_err

    residuals = x-y
    stats = im_stats(residuals)
    xstr = strtrim(string(stats.median,format='(F12.2)'),2)+'\pm'+strtrim(string(stats.sig68mean,format='(F12.2)'),2)

    ytitle = ohtitle

    sdss_lineplot, x, y, xerr, yerr, plottype=1, postscript=postscript, $
      xtitle='', ytitle=ytitle, xrange=xrange, yrange=yrange, legendtype=0, $
      /right, /top, position=pos[*,0], charsize=charsize_4, xtickname=replicate(' ',10)
    djs_oplot, !x.crange, !y.crange, thick=postthick, color=djs_icolor(talkcolor), line=2

    legend, textoidl(xstr), /right, /bottom, box=0, charsize=charsize_4, $
      charthick=postthick, textcolor=djs_icolor(talkcolor)
    legend, textoidl('(a) ([O III]/H\beta)/([N II]/H\alpha) (Pettini & Pagel 2004)'), /left, /top, box=0, $
      charsize=charsize_4, charthick=postthick, textcolor=djs_icolor(talkcolor)

; --------------------------------------------------    
; [N II]/Ha
; --------------------------------------------------    

    indx = where((sdssnodust.zstrong_12oh_niiha_pettini gt -900) and (sdssancillary.tremonti_oh gt 0.0),nindx)

    x = sdssancillary[indx].tremonti_oh
    xerr = sdssancillary[indx].tremonti_oh_err

    y = sdssnodust[indx].zstrong_12oh_niiha_pettini
    yerr = sdssnodust[indx].zstrong_12oh_niiha_pettini_err

    residuals = x-y
    stats = im_stats(residuals)
    xstr = strtrim(string(stats.median,format='(F12.2)'),2)+'\pm'+strtrim(string(stats.sig68mean,format='(F12.2)'),2)

    ytitle = ohtitle

    sdss_lineplot, x, y, xerr, yerr, plottype=1, postscript=postscript, $
      xtitle='', ytitle=ytitle, xrange=xrange, yrange=yrange, legendtype=0, /noerase, $
      /right, /top, position=pos[*,1], charsize=charsize_4, xtickname=replicate(' ',10)
    djs_oplot, !x.crange, !y.crange, thick=postthick, color=djs_icolor(talkcolor), line=2

    legend, textoidl(xstr), /right, /bottom, box=0, charsize=charsize_4, $
      charthick=postthick, textcolor=djs_icolor(talkcolor)
    legend, textoidl('(b) [N II]/H\alpha (Pettini & Pagel 2004)'), /left, /top, box=0, $
      charsize=charsize_4, charthick=postthick, textcolor=djs_icolor(talkcolor)

; --------------------------------------------------    
; P-method
; --------------------------------------------------    

    indx = where((sdssancillary.tremonti_oh gt 0.0) and (sdssnodust.zstrong_12oh_p01_upper gt -900.0) and $
      (sdssnodust.zstrong_12oh_p01_lower gt -900.0),nindx)

    up = where((sdssancillary[indx].tremonti_oh gt 8.2),nup)
    if (nup ne 0L) then begin
       x = sdssancillary[indx[up]].tremonti_oh
       xerr = sdssancillary[indx[up]].tremonti_oh_err

       y = sdssnodust[indx[up]].zstrong_12oh_p01_upper
       yerr = sdssnodust[indx[up]].zstrong_12oh_p01_upper_err
    endif

    lo = where((sdssancillary[indx].tremonti_oh lt 7.95),nlo)
    if (nlo ne 0L) then begin
       x = [x,sdssancillary[indx[lo]].tremonti_oh]
       xerr = [xerr,[sdssnodust[indx[lo]].tremonti_p84-sdssnodust[indx[lo]].tremonti_p16]/2.0]

       y = [y,sdssnodust[indx[lo]].zstrong_12oh_p01_lower]
       yerr = [yerr,sdssnodust[indx[lo]].zstrong_12oh_p01_lower_err]
    endif

    residuals = x-y
    stats = im_stats(residuals)
    xstr = strtrim(string(stats.median,format='(F12.2)'),2)+'\pm'+strtrim(string(stats.sig68mean,format='(F12.2)'),2)

    ytitle = ohtitle

    sdss_lineplot, x, y, xerr, yerr, plottype=1, postscript=postscript, $
      xtitle=xtitle, ytitle=ytitle, xrange=xrange, yrange=yrange, legendtype=0, $
      /right, /top, position=pos[*,2], charsize=charsize_4, /noerase
    djs_oplot, !x.crange, !y.crange, thick=postthick, color=djs_icolor(talkcolor), line=2

    legend, textoidl(xstr), /right, /bottom, box=0, charsize=charsize_4, $
      charthick=postthick, textcolor=djs_icolor(talkcolor)
    legend, textoidl('(c) P-method'), /left, /top, box=0, $
      charsize=charsize_4, charthick=postthick, textcolor=djs_icolor(talkcolor)

    im_openclose, postscript=postscript, /close    
    
; ------------------------------------------------------------
; 3-panel 12+log(O/H) [T04] vs 12+log(O/H) Empirical, Pettini & Pagel
; (2004) calibrations - residuals
; ------------------------------------------------------------

    psname = '12oh_T04_vs_12oh_empirical_pettini_residuals'
    im_openclose, pspath+psname, postscript=postscript, encapsulated=encapsulated, xsize=8.5, ysize=10.5

    pagemaker, nx=1, ny=3, height=3.0*[1,1,1], width=6.8, xmargin=[1.3,0.4], $
      ymargin=[0.4,1.1], xspace=0, yspace=0, xpage=8.5, ypage=10.5, $
      position=pos, /normal

    if keyword_set(dotalk) then polyfill, [1,1,0,0,1], [1,0,0,1,1], /normal, $
      color=djs_icolor('black')

    xtitle = ohtitle+' [T04]'
    ytitle = '\Delta[log(O/H)]'

    xrange = ohrange2
    yrange = [-1.2,1.2]

; --------------------------------------------------    
; ([O III]/Hb)/([N II]/Ha)
; --------------------------------------------------    

    indx = where((sdssnodust.zstrong_12oh_oiiinii_moustakas gt -900) and (sdssancillary.tremonti_oh gt 0.0),nindx)

    x = sdssancillary[indx].tremonti_oh
    xerr = sdssancillary[indx].tremonti_oh_err

    y = sdssnodust[indx].zstrong_12oh_oiiinii_moustakas
    yerr = sdssnodust[indx].zstrong_12oh_oiiinii_moustakas_err

    residuals = x-y
    residuals_err = sqrt(xerr^2+yerr^2)

    stats = im_stats(residuals)
    xstr = strtrim(string(stats.median,format='(F12.2)'),2)+'\pm'+strtrim(string(stats.sig68mean,format='(F12.2)'),2)

    sdss_lineplot, x, residuals, xerr, residuals_err, plottype=1, postscript=postscript, $
      xtitle='', ytitle=ytitle, xrange=xrange, yrange=yrange, legendtype=0, $
      /right, /top, position=pos[*,0], charsize=charsize_4, xtickname=replicate(' ',10)
    djs_oplot, !x.crange, [0,0], thick=postthick, color=djs_icolor(talkcolor), line=2

    legend, textoidl(xstr), /right, /bottom, box=0, charsize=charsize_4, $
      charthick=postthick, textcolor=djs_icolor(talkcolor)
    legend, textoidl('(a) ([O III]/H\beta)/([N II]/H\alpha) (Pettini & Pagel 2004)'), /left, /top, box=0, $
      charsize=charsize_4, charthick=postthick, textcolor=djs_icolor(talkcolor)

; --------------------------------------------------    
; [N II]/Ha
; --------------------------------------------------    

    indx = where((sdssnodust.zstrong_12oh_niiha_moustakas gt -900) and (sdssancillary.tremonti_oh gt 0.0),nindx)

    x = sdssancillary[indx].tremonti_oh
    xerr = sdssancillary[indx].tremonti_oh_err

    y = sdssnodust[indx].zstrong_12oh_niiha_moustakas
    yerr = sdssnodust[indx].zstrong_12oh_niiha_moustakas_err

    residuals = x-y
    residuals_err = sqrt(xerr^2+yerr^2)

    stats = im_stats(residuals)
    xstr = strtrim(string(stats.median,format='(F12.2)'),2)+'\pm'+strtrim(string(stats.sig68mean,format='(F12.2)'),2)

    sdss_lineplot, x, residuals, xerr, residuals_err, plottype=1, postscript=postscript, $
      xtitle='', ytitle=ytitle, xrange=xrange, yrange=yrange, legendtype=0, /noerase, $
      /right, /top, position=pos[*,1], charsize=charsize_4, xtickname=replicate(' ',10)
    djs_oplot, !x.crange, [0,0], thick=postthick, color=djs_icolor(talkcolor), line=2

    legend, textoidl(xstr), /right, /bottom, box=0, charsize=charsize_4, $
      charthick=postthick, textcolor=djs_icolor(talkcolor)
    legend, textoidl('(b) [N II]/H\alpha (Pettini & Pagel 2004)'), /left, /top, box=0, $
      charsize=charsize_4, charthick=postthick, textcolor=djs_icolor(talkcolor)

; --------------------------------------------------    
; P-method
; --------------------------------------------------    

    indx = where((sdssancillary.tremonti_oh gt 0.0) and (sdssnodust.zstrong_12oh_p01_upper gt -900.0) and $
      (sdssnodust.zstrong_12oh_p01_lower gt -900.0),nindx)

    up = where((sdssancillary[indx].tremonti_oh gt 8.2),nup)
    if (nup ne 0L) then begin
       x = sdssancillary[indx[up]].tremonti_oh
       xerr = sdssancillary[indx[up]].tremonti_oh_err

       y = sdssnodust[indx[up]].zstrong_12oh_p01_upper
       yerr = sdssnodust[indx[up]].zstrong_12oh_p01_upper_err
    endif

    lo = where((sdssancillary[indx].tremonti_oh lt 7.95),nlo)
    if (nlo ne 0L) then begin
       x = [x,sdssancillary[indx[lo]].tremonti_oh]
       xerr = [xerr,sdssancillary[indx[lo]].tremonti_oh_err]

       y = [y,sdssnodust[indx[lo]].zstrong_12oh_p01_lower]
       yerr = [yerr,sdssnodust[indx[lo]].zstrong_12oh_p01_lower_err]
    endif

    residuals = x-y
    residuals_err = sqrt(xerr^2+yerr^2)

    stats = im_stats(residuals)
    xstr = strtrim(string(stats.median,format='(F12.2)'),2)+'\pm'+strtrim(string(stats.sig68mean,format='(F12.2)'),2)

    sdss_lineplot, x, residuals, xerr, residuals_err, plottype=1, postscript=postscript, $
      xtitle=xtitle, ytitle=ytitle, xrange=xrange, yrange=yrange, legendtype=0, $
      /right, /top, position=pos[*,2], charsize=charsize_4, /noerase
    djs_oplot, !x.crange, [0,0], thick=postthick, color=djs_icolor(talkcolor), line=2

    legend, textoidl(xstr), /right, /bottom, box=0, charsize=charsize_4, $
      charthick=postthick, textcolor=djs_icolor(talkcolor)
    legend, textoidl('(c) P-method'), /left, /top, box=0, $
      charsize=charsize_4, charthick=postthick, textcolor=djs_icolor(talkcolor)

    im_openclose, postscript=postscript, /close    

; ------------------------------------------------------------
; D(4000) vs alpha_beta_3
; ------------------------------------------------------------

    psname = 'D4000_vs_ab3_models'
    im_openclose, pspath+psname, postscript=postscript

    x = sfhgrid.d4000_narrow
    xerr = x*0.0

    y = sfhgrid.oiii_5007_continuum/sfhgrid.h_beta_continuum ; alpha_beta_3
    yerr = y*0.0
    
    xtitle = 'D_{n}(4000)'
    ytitle = '\alpha_{\beta3}'

    xrange = [1.0,2.3]
    yrange = [0.85,1.15]
    
    atlas1d_lineplot, x, y, xerr, yerr, plottype=3, postscript=postscript, $
      xtitle=xtitle, ytitle=ytitle, xrange=xrange, yrange=yrange, legendtype=0, $
      /left, /top

    im_openclose, postscript=postscript, /close
    
; ------------------------------------------------------------
; D(4000) vs alpha_2_beta
; ------------------------------------------------------------

    psname = 'D4000_vs_a2b'
    im_openclose, pspath+psname, postscript=postscript

    lineratio, atlasnodust, 'D4000_NARROW', '', 'OII_3727_CONTINUUM', 'H_BETA_CONTINUUM', $
      x, xerr, y, yerr, index=indx, nindex=nindx, snrcut=snrcut, /nolog
    
    xmodel = sfhgrid.d4000_narrow
    ymodel = sfhgrid.oii_3727_continuum/sfhgrid.h_beta_continuum ; alpha_2_beta
    srt = sort(xmodel)
    
    xtitle = 'D_{n}(4000)'
    ytitle = '\alpha_{2\beta}'

    xrange = [0.65,2.3]
    yrange = [0.2,1.6]
    
    atlas1d_lineplot, x, y, xerr, yerr, plottype=3, postscript=postscript, $
      xtitle=xtitle, ytitle=ytitle, xrange=xrange, yrange=yrange, legendtype=0, $
      /left, /top

    djs_oplot, xmodel[srt], ymodel[srt], ps=4
;   djs_oplot, xmodel[srt], ymodel[srt]
;   polyfill, xmodel[srt], ymodel[srt], /fill, color=djs_icolor('grey')
    
    im_openclose, postscript=postscript, /close
    
; ------------------------------------------------------------
; B-V vs alpha_2_beta
; ------------------------------------------------------------

    psname = 'BV_vs_a2b_models'
    im_openclose, pspath+psname, postscript=postscript

    x = sfhgrid.B-sfhgrid.V
    xerr = x*0.0

    y = sfhgrid.oii_3727_continuum/sfhgrid.h_beta_continuum ; alpha_2_beta
    yerr = y*0.0
    
    xtitle = 'B-V'
    ytitle = '\alpha_{2\beta}'

    xrange = [-0.1,1.1]
    yrange = [0.2,1.6]
    
    atlas1d_lineplot, x, y, xerr, yerr, plottype=3, postscript=postscript, $
      xtitle=xtitle, ytitle=ytitle, xrange=xrange, yrange=yrange, legendtype=0, $
      /left, /top

    im_openclose, postscript=postscript, /close
    
; ------------------------------------------------------------
; 12+log (O/H) [N II]/Ha versus 12+log (O/H) R23/M91
; ------------------------------------------------------------

    psname = 'sdss_12oh_niiha_vs_12oh_m91_upper'
    im_openclose, pspath+psname, postscript=postscript, xsize=8.5, ysize=7.5, /encapsulated

    pagemaker, nx=1, ny=1, xspace=0, yspace=0.0, width=5.8, height=5.5, $
      xmargin=[2.2,0.5], ymargin=[0.9,1.1], xpage=8.5, ypage=7.5, $
      position=pos, /normal

    indx = where((sdssnodust.zstrong_12oh_niiha_pettini gt -900) and (sdssnodust.zstrong_12oh_m91_upper gt 0.0),nindx)

    x = sdssnodust[indx].zstrong_12oh_niiha_pettini
    xerr = sdssnodust[indx].zstrong_12oh_niiha_pettini_err

    y = sdssnodust[indx].zstrong_12oh_m91_upper
    yerr = sdssnodust[indx].zstrong_12oh_m91_upper_err

    residuals = x-y
    stats = im_stats(residuals,/verbose)
    xstr = strtrim(string(stats.median,format='(F12.2)'),2)+'\pm'+strtrim(string(stats.sig68mean,format='(F12.2)'),2)

    xtitle = ohtitle+' [N II]/H\alpha'
    ytitle = ohtitle+' R_{23}/M91 Upper'

    xrange = ohrange2
    yrange = xrange

    sdss_lineplot, x, y, xerr, yerr, plottype=3, postscript=postscript, $
      xtitle=xtitle, ytitle=ytitle, xrange=xrange, yrange=yrange, legendtype=0, $
      /right, /top, position=pos[*,0], charsize=charsize_8
    djs_oplot, !x.crange, !y.crange, line=0, thick=postthick
    legend, textoidl(xstr), /left, /top, box=0, charsize=charsize_8, charthick=postthick

    im_openclose, postscript=postscript, /close    

;; --------------------------------------------------    
;; GENERATE TALK PLOTS
;; --------------------------------------------------    
;
;    if keyword_set(talk) then begin
;       if (n_elements(dotalk) eq 0L) then abundances, atlasdust, atlasnodust, $
;         nfgsdust=nfgsdust, nfgsnodust=nfgsnodust, hii=hii, sdssdust=sdssdust, $
;         sdssnodust=sdssnodust, suffix='.talk', encapsulated=0, /talk, /postscript, $
;         /dotalk, _extra=extra else return
;    endif

; --------------------------------------------------    
; GENERATE ENCAPSULATED POSTSCRIPT
; --------------------------------------------------    

;   if keyword_set(postscript) then if (n_elements(doencapsulated) eq 0L) then $
;     abundances, atlasdust, atlasnodust, nfgsdust=nfgsdust, nfgsnodust=nfgsnodust, $
;       hii=hii, sdssdust=sdssdust, sdssnodust=sdssnodust, postscript=postscript, $
;       /encapsulated, /doencapsulated, talk=0, _extra=extra else return
    
