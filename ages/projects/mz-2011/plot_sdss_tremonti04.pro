pro plot_sdss_tremonti04, tdust, tnodust, tancillary, postscript=postscript
; jm06arp30uofa - make some plots

    if (n_elements(tdust) eq 0L) then tdust = read_sdss_tremonti04_sample(sdssnodust=tnodust,sdssancillary=tancillary)

    if keyword_set(postscript) then begin
       postthick = 5.0 
       postthick2 = 8.0 
       dfpsplot, ages_path(/projects)+'mz/plot_sdss_tremonti04.ps', /color, /square
    endif else begin
       postthick = 2.0
       postthick2 = 2.0
       im_window, 0, xratio=0.6, /square
    endelse

    imf_kroupa_to_salpeter = alog10(1.5) ; Brinchmann et al. (2004) - differs from Bell et al. 2003 (0.3 dex)!

    tremonti_mass = findgen((11.5-8.3)/0.01+1)*0.01+8.3
    tremonti_mass_coeff = [-1.492,1.847,0.08026]
    tremonti_mass_oh = tremonti_mass_coeff[0] + tremonti_mass_coeff[1]*(tremonti_mass) - $ ; Kroupa
      tremonti_mass_coeff[2]*(tremonti_mass)^2
    tremonti_mass_oh_salp = tremonti_mass_coeff[0] + tremonti_mass_coeff[1]*(tremonti_mass+imf_kroupa_to_salpeter) - $ ; Salpeter
      tremonti_mass_coeff[2]*(tremonti_mass+imf_kroupa_to_salpeter)^2

    tremonti_mb = findgen(((-10.0)-(-30.0))/0.01)*0.01+(-30.0)
    tremonti_mb_coeff = [5.129,-0.202]
;   tremonti_mb_coeff = [5.276,-0.186]
    tremonti_mb_coeff_err = [0.018,0.001]
    tremonti_mb_oh = poly(tremonti_mb,tremonti_mb_coeff)

    tremonti_mb_coeffstr = '('+strtrim(string(tremonti_mb_coeff[0],format='(F12.3)'),2)+$
      '+/-'+strtrim(string(tremonti_mb_coeff_err[0],format='(F12.3)'),2)+') + ('+$
      strtrim(string(tremonti_mb_coeff[1],format='(F12.3)'),2)+'+/-'+$
      strtrim(string(tremonti_mb_coeff_err[1],format='(F12.3)'),2)+') M_{B}'

    mbaxis = findgen(((-10.0)-(-30.0))/0.01)*0.01+(-30.0)

; ---------------------------------------------------------------------------    
; LZ relation split by inclination [KK04/fluxcor]
; ---------------------------------------------------------------------------    

;   psname = 'ages_sdss_oiihb_vs_oiiihb'
;   im_openclose, pspath+psname, postscript=postscript, xsize=8.5, ysize=4.8, encapsulated=encapsulated

    pagemaker, nx=2, ny=1, height=3.5, width=[3.5,3.5], xmargin=[1.2,0.3], $
      ymargin=[0.3,1.0], xspace=0.0, yspace=0, position=pos, /normal, $
      xpage=8.5, ypage=4.8

    xtitle = 'M_{B, cor} [mag]'
    ytitle = '12 + log (O/H)_{cor} [KK04]'

    xrange = [-16.5,-23.1]
    yrange = [8.5,9.3]

    indx = where((tancillary.m_b gt -900.0) and (tnodust.zstrong_12oh_kk04 gt -900),nindx)

    mb = tancillary[indx].m_b - tancillary[indx].a_b_incl
    mberr = tancillary[indx].m_b_err
    oh = tnodust[indx].zstrong_12oh_kk04
    oherr = tnodust[indx].zstrong_12oh_kk04_err

; -------------------------    
; fit to everything
; -------------------------    

    sixlin, mb, oh, a, siga, b, sigb
    allcoeff = [a[2],b[2]] & allcoeff_err = [siga[2],sigb[2]]
    allcoeffstr = 'All    : 12+log (O/H)_{cor} = ('+strtrim(string(allcoeff[0],format='(F12.3)'),2)+$
      '+/-'+strtrim(string(allcoeff_err[0],format='(F12.3)'),2)+') + ('+$
      strtrim(string(allcoeff[1],format='(F12.3)'),2)+'+/-'+strtrim(string(allcoeff_err[1],format='(F12.3)'),2)+') M_{B}'
    allcoeffstr2 = 'All ('+strtrim(string(allcoeff[1],format='(F12.3)'),2)+')'
    splog, allcoeffstr

; -------------------------    
; edge-on
; -------------------------    

    edgeon = where(tancillary[indx].inclination gt 45.0)
    
    sixlin, mb[edgeon], oh[edgeon], a, siga, b, sigb
    edgecoeff = [a[2],b[2]] & edgecoeff_err = [siga[2],sigb[2]]

    edgecoeffstr = 'Edge-On: 12+log (O/H)_{cor} = ('+strtrim(string(edgecoeff[0],format='(F12.3)'),2)+$
      '+/-'+strtrim(string(edgecoeff_err[0],format='(F12.3)'),2)+') + ('+$
      strtrim(string(edgecoeff[1],format='(F12.3)'),2)+'+/-'+strtrim(string(edgecoeff_err[1],format='(F12.3)'),2)+') M_{B}'
    edgecoeffstr2 = 'Edge-On ('+strtrim(string(edgecoeff[1],format='(F12.3)'),2)+')'
    splog, edgecoeffstr

; -------------------------    
; face-on
; -------------------------    

    faceon = where(tancillary[indx].inclination lt 45.0)
    
    sixlin, mb[faceon], oh[faceon], a, siga, b, sigb
    facecoeff = [a[2],b[2]] & facecoeff_err = [siga[2],sigb[2]]

    facecoeffstr = 'Face-On: 12+log (O/H)_{cor} = ('+strtrim(string(facecoeff[0],format='(F12.3)'),2)+$
      '+/-'+strtrim(string(facecoeff_err[0],format='(F12.3)'),2)+') + ('+$
      strtrim(string(facecoeff[1],format='(F12.3)'),2)+'+/-'+strtrim(string(facecoeff_err[1],format='(F12.3)'),2)+') M_{B}'
    facecoeffstr2 = 'Face-On ('+strtrim(string(facecoeff[1],format='(F12.3)'),2)+')'
    splog, facecoeffstr
    print
    
; -------------------------    
; make the plot!
; -------------------------    

    sdss_lineplot, mb[edgeon], oh[edgeon], mberr[edgeon], oherr[edgeon], plottype=1, postscript=postscript, $
      xtitle=xtitle, ytitle=ytitle, xrange=xrange, yrange=yrange, position=pos[*,0], $
      charsize=1.5, /xreverse, xminor=5
    legend, '(a) Edge-On (i>45)', /left, /top, box=0, charsize=1.5, charthick=postthick

    djs_oplot, mbaxis, poly(mbaxis,allcoeff), line=0, thick=postthick2;, color='blue'
    djs_oplot, mbaxis, poly(mbaxis,edgecoeff), line=2, thick=postthick2, color='red'
    djs_oplot, mbaxis, poly(mbaxis,facecoeff), line=3, thick=postthick2, color='blue'

; -------------------------    
; face-on
; -------------------------    

    bin = where(tancillary[indx].inclination lt 45.0)
    
    sdss_lineplot, mb[faceon], oh[faceon], mberr[faceon], oherr[faceon], plottype=1, postscript=postscript, $
      xtitle=xtitle, ytitle='', xrange=xrange, yrange=yrange, position=pos[*,1], $
      charsize=1.5, /xreverse, /noerase, ytickname=replicate(' ',10), xminor=5
    legend, '(b) Face-On (i<45)', /left, /top, box=0, charsize=1.5, charthick=postthick

    djs_oplot, mbaxis, poly(mbaxis,allcoeff), line=0, thick=postthick2;, color='blue'
    djs_oplot, mbaxis, poly(mbaxis,edgecoeff), line=2, thick=postthick2, color='red'
    djs_oplot, mbaxis, poly(mbaxis,facecoeff), line=3, thick=postthick2, color='blue'

    legend, textoidl([allcoeffstr2,edgecoeffstr2,facecoeffstr2]), /right, /bottom, box=0, charsize=1.1, $
      charthick=postthick, thick=postthick2, line=[0,2,3], color=djs_icolor(['','red','blue']), $
      textcolor=djs_icolor(['','red','blue'])
;   legend, ['All','Edge-On','Face-On'], /right, /bottom, box=0, charsize=1.1, $
;     charthick=postthick, thick=postthick2, line=[0,2,3], color=djs_icolor(['','red','blue']), $
;     textcolor=djs_icolor(['','red','blue'])
    
;   im_openclose, postscript=postscript, /close
    if (not keyword_set(postscript)) then cc = get_kbrd(1)

; ---------------------------------------------------------------------------    
; LZ relation split by inclination [KK04/EW]
; ---------------------------------------------------------------------------    

;   psname = 'ages_sdss_oiihb_vs_oiiihb'
;   im_openclose, pspath+psname, postscript=postscript, xsize=8.5, ysize=4.8, encapsulated=encapsulated

    pagemaker, nx=2, ny=1, height=3.5, width=[3.5,3.5], xmargin=[1.2,0.3], $
      ymargin=[0.3,1.0], xspace=0.0, yspace=0, position=pos, /normal, $
      xpage=8.5, ypage=4.8

    xtitle = 'M_{B, cor} [mag]'
    ytitle = '12 + log (O/H)_{EW} [KK04]'

    xrange = [-16.5,-23.1]
    yrange = [8.5,9.3]

    indx = where((tancillary.m_b gt -900.0) and (tdust.zstrong_ew_12oh_kk04 gt -900),nindx)

    mb = tancillary[indx].m_b - tancillary[indx].a_b_incl
    mberr = tancillary[indx].m_b_err
    oh = tdust[indx].zstrong_ew_12oh_kk04
    oherr = tdust[indx].zstrong_ew_12oh_kk04_err

; -------------------------    
; fit to everything
; -------------------------    

    sixlin, mb, oh, a, siga, b, sigb
    allcoeff = [a[2],b[2]] & allcoeff_err = [siga[2],sigb[2]]
    allcoeffstr = 'All    : 12+log (O/H)_{EW} = ('+strtrim(string(allcoeff[0],format='(F12.3)'),2)+$
      '+/-'+strtrim(string(allcoeff_err[0],format='(F12.3)'),2)+') + ('+$
      strtrim(string(allcoeff[1],format='(F12.3)'),2)+'+/-'+strtrim(string(allcoeff_err[1],format='(F12.3)'),2)+') M_{B}'
    allcoeffstr2 = 'All ('+strtrim(string(allcoeff[1],format='(F12.3)'),2)+')'
    splog, allcoeffstr

; -------------------------    
; edge-on
; -------------------------    

    edgeon = where(tancillary[indx].inclination gt 45.0)
    
    sixlin, mb[edgeon], oh[edgeon], a, siga, b, sigb
    edgecoeff = [a[2],b[2]] & edgecoeff_err = [siga[2],sigb[2]]

    edgecoeffstr = 'Edge-On: 12+log (O/H)_{EW} = ('+strtrim(string(edgecoeff[0],format='(F12.3)'),2)+$
      '+/-'+strtrim(string(edgecoeff_err[0],format='(F12.3)'),2)+') + ('+$
      strtrim(string(edgecoeff[1],format='(F12.3)'),2)+'+/-'+strtrim(string(edgecoeff_err[1],format='(F12.3)'),2)+') M_{B}'
    edgecoeffstr2 = 'Edge-On ('+strtrim(string(edgecoeff[1],format='(F12.3)'),2)+')'
    splog, edgecoeffstr

; -------------------------    
; face-on
; -------------------------    

    faceon = where(tancillary[indx].inclination lt 45.0)
    
    sixlin, mb[faceon], oh[faceon], a, siga, b, sigb
    facecoeff = [a[2],b[2]] & facecoeff_err = [siga[2],sigb[2]]

    facecoeffstr = 'Face-On: 12+log (O/H)_{EW} = ('+strtrim(string(facecoeff[0],format='(F12.3)'),2)+$
      '+/-'+strtrim(string(facecoeff_err[0],format='(F12.3)'),2)+') + ('+$
      strtrim(string(facecoeff[1],format='(F12.3)'),2)+'+/-'+strtrim(string(facecoeff_err[1],format='(F12.3)'),2)+') M_{B}'
    facecoeffstr2 = 'Face-On ('+strtrim(string(facecoeff[1],format='(F12.3)'),2)+')'
    splog, facecoeffstr
    print
    
; -------------------------    
; make the plot!
; -------------------------    

    sdss_lineplot, mb[edgeon], oh[edgeon], mberr[edgeon], oherr[edgeon], plottype=1, postscript=postscript, $
      xtitle=xtitle, ytitle=ytitle, xrange=xrange, yrange=yrange, position=pos[*,0], $
      charsize=1.5, /xreverse, xminor=5
    legend, '(a) Edge-On (i>45)', /left, /top, box=0, charsize=1.5, charthick=postthick

    djs_oplot, mbaxis, poly(mbaxis,allcoeff), line=0, thick=postthick2;, color='blue'
    djs_oplot, mbaxis, poly(mbaxis,edgecoeff), line=2, thick=postthick2, color='red'
    djs_oplot, mbaxis, poly(mbaxis,facecoeff), line=3, thick=postthick2, color='blue'

; -------------------------    
; face-on
; -------------------------    

    bin = where(tancillary[indx].inclination lt 45.0)
    
    sdss_lineplot, mb[faceon], oh[faceon], mberr[faceon], oherr[faceon], plottype=1, postscript=postscript, $
      xtitle=xtitle, ytitle='', xrange=xrange, yrange=yrange, position=pos[*,1], $
      charsize=1.5, /xreverse, /noerase, ytickname=replicate(' ',10), xminor=5
    legend, '(b) Face-On (i<45)', /left, /top, box=0, charsize=1.5, charthick=postthick

    djs_oplot, mbaxis, poly(mbaxis,allcoeff), line=0, thick=postthick2;, color='blue'
    djs_oplot, mbaxis, poly(mbaxis,edgecoeff), line=2, thick=postthick2, color='red'
    djs_oplot, mbaxis, poly(mbaxis,facecoeff), line=3, thick=postthick2, color='blue'

    legend, textoidl([allcoeffstr2,edgecoeffstr2,facecoeffstr2]), /right, /bottom, box=0, charsize=1.1, $
      charthick=postthick, thick=postthick2, line=[0,2,3], color=djs_icolor(['','red','blue']), $
      textcolor=djs_icolor(['','red','blue'])
;   legend, ['All','Edge-On','Face-On'], /right, /bottom, box=0, charsize=1.1, $
;     charthick=postthick, thick=postthick2, line=[0,2,3], color=djs_icolor(['','red','blue']), $
;     textcolor=djs_icolor(['','red','blue'])
    
;   im_openclose, postscript=postscript, /close
    if (not keyword_set(postscript)) then cc = get_kbrd(1)

; ---------------------------------------------------------------------------    
; LZ (Figure 4)
; ---------------------------------------------------------------------------    
    
;   psname = 't04_lz_sdss'
;   im_openclose, psname, postscript=postscript, xsize=8.5, ysize=8.5, encapsulated=encapsulated

    pagemaker, nx=1, ny=1, xspace=0, yspace=0.0, width=6.5, height=6.5, $
      xmargin=[1.5,0.5], ymargin=[0.9,1.1], xpage=8.5, ypage=8.5, $
      position=pos, /normal

    indx = where((tancillary.m_b gt -900.0) and (tancillary.tremonti_oh gt -900),nindx)

    mb = tancillary[indx].m_b
    mberr = tancillary[indx].m_b_err
    oh = tancillary[indx].tremonti_oh
    oherr = tancillary[indx].tremonti_oh_err
    
    xtitle = 'M_{B, obs} [mag] [K-correct]'
    ytitle = '12 + log (O/H)_{cor} [Tremonti]'

    xrange = [-14,-24]
    yrange = [8.0,9.5]

    sdss_lineplot, mb, oh, mberr, oherr, plottype=1, postscript=postscript, $
      xtitle=xtitle, ytitle=ytitle, xrange=xrange, yrange=yrange, position=pos, $
      charsize=1.8, /xreverse

    sixlin, mb, oh, a, siga, b, sigb
    lzcoeff = [a[2],b[2]] & lzcoeff_err = [siga[2],sigb[2]]
    lzcoeffstr = '('+strtrim(string(lzcoeff[0],format='(F12.3)'),2)+$
      '+/-'+strtrim(string(lzcoeff_err[0],format='(F12.3)'),2)+') + ('+$
      strtrim(string(lzcoeff[1],format='(F12.3)'),2)+'+/-'+strtrim(string(lzcoeff_err[1],format='(F12.3)'),2)+') M_{B}'

    lzfit = poly(mbaxis,lzcoeff)
    djs_oplot, mbaxis, lzfit, line=0, thick=postthick, color='red'

    djs_oplot, tremonti_mb, tremonti_mb_oh, line=2, thick=postthick, color='blue'

    legend, textoidl(['LZ : '+lzcoeffstr,'T04: '+tremonti_mb_coeffstr]), /left, /top, box=0, charsize=1.1, $
      charthick=postthick, thick=postthick2, line=[0,2], color=djs_icolor(['red','blue']), $
      textcolor=djs_icolor(['red','blue'])

    splog, 'LZ : '+lzcoeffstr
    splog, 'T04: '+tremonti_mb_coeffstr
    print

;   im_openclose, postscript=postscript, /close
    if (not keyword_set(postscript)) then cc = get_kbrd(1)

; ---------------------------------------------------------------------------    
; LZ: KK04/EWs
; ---------------------------------------------------------------------------    
    
    pagemaker, nx=1, ny=1, xspace=0, yspace=0.0, width=6.5, height=6.5, $
      xmargin=[1.5,0.5], ymargin=[0.9,1.1], xpage=8.5, ypage=8.5, $
      position=pos, /normal

    indx = where((tancillary.m_b gt -900.0) and (tdust.zstrong_ew_12oh_kk04 gt -900),nindx)

    mb = tancillary[indx].m_b
    mberr = tancillary[indx].m_b_err
    oh = tdust[indx].zstrong_ew_12oh_kk04
    oherr = tdust[indx].zstrong_ew_12oh_kk04_err
    
    xtitle = 'M_{B, obs} [mag] [K-correct]'
    ytitle = '12 + log (O/H)_{EW} [KK04]'

    xrange = [-14,-24]
    yrange = [8.0,9.5]

    sdss_lineplot, mb, oh, mberr, oherr, plottype=1, postscript=postscript, $
      xtitle=xtitle, ytitle=ytitle, xrange=xrange, yrange=yrange, position=pos, $
      charsize=1.8, /xreverse

    sixlin, mb, oh, a, siga, b, sigb
    lzcoeff = [a[2],b[2]] & lzcoeff_err = [siga[2],sigb[2]]
    lzcoeffstr = '('+strtrim(string(lzcoeff[0],format='(F12.3)'),2)+$
      '+/-'+strtrim(string(lzcoeff_err[0],format='(F12.3)'),2)+') + ('+$
      strtrim(string(lzcoeff[1],format='(F12.3)'),2)+'+/-'+strtrim(string(lzcoeff_err[1],format='(F12.3)'),2)+') M_{B}'

    lzfit = poly(mbaxis,lzcoeff)
    djs_oplot, mbaxis, lzfit, line=0, thick=postthick, color='red'

    djs_oplot, tremonti_mb, tremonti_mb_oh, line=2, thick=postthick, color='blue'

    legend, textoidl(['LZ : '+lzcoeffstr,'T04: '+tremonti_mb_coeffstr]), /left, /top, box=0, charsize=1.1, $
      charthick=postthick, thick=postthick2, line=[0,2], color=djs_icolor(['red','blue']), $
      textcolor=djs_icolor(['red','blue'])

    splog, 'LZ : '+lzcoeffstr
    splog, 'T04: '+tremonti_mb_coeffstr
    print

;   im_openclose, postscript=postscript, /close
    if (not keyword_set(postscript)) then cc = get_kbrd(1)

; ---------------------------------------------------------------------------    
; MZ - Kauffmann masses
; ---------------------------------------------------------------------------    
    
    pagemaker, nx=1, ny=1, xspace=0, yspace=0.0, width=6.5, height=6.5, $
      xmargin=[1.5,0.5], ymargin=[0.9,1.1], xpage=8.5, ypage=8.5, $
      position=pos, /normal

    indx = where((tancillary.kauffmann_mass gt -900.0) and (tdust.zstrong_ew_12oh_kk04 gt -900),nindx)

    mass = tancillary[indx].kauffmann_mass ; Salpeter
    masserr = tancillary[indx].kauffmann_mass_err
    oh = tdust[indx].zstrong_ew_12oh_kk04
    oherr = tdust[indx].zstrong_ew_12oh_kk04_err
    
    xtitle = 'log M [M'+sunsymbol()+']  [Kauffmann, Salpeter]'
    ytitle = '12 + log (O/H)_{EW} [KK04]'

    xrange = [7.8,11.8]
    yrange = [8.0,9.5]

    sdss_lineplot, mass, oh, masserr, oherr, plottype=1, postscript=postscript, $
      xtitle=xtitle, ytitle=ytitle, xrange=xrange, yrange=yrange, position=pos, $
      charsize=1.8

    result = im_medxbin(mass,oh,0.1,minx=8.57,minpts=100L);,/verbose)
    coeff = poly_fit(result.binctr,result.medy,2)
    massaxis = findgen((max(result.binctr)-min(result.binctr))/0.01+1)*0.01+min(result.binctr)
    t04fit = poly(massaxis,coeff)
    
    im_symbols, 106, psize=1.1, thick=postthick, fill=1
    djs_oplot, result.binctr, result.medy, ps=8
    djs_oplot, massaxis, t04fit, line=0, thick=postthick, color='red'

    djs_oplot, tremonti_mass, tremonti_mass_oh_salp, line=2, thick=postthick, color='blue'
    niceprint, coeff, tremonti_mass_coeff

    if (not keyword_set(postscript)) then cc = get_kbrd(1)

; ---------------------------------------------------------------------------    
; MZ - Kcorrect masses
; ---------------------------------------------------------------------------    
    
    pagemaker, nx=1, ny=1, xspace=0, yspace=0.0, width=6.5, height=6.5, $
      xmargin=[1.5,0.5], ymargin=[0.9,1.1], xpage=8.5, ypage=8.5, $
      position=pos, /normal

    indx = where((tancillary.kcorr_mass gt -900.0) and (tdust.zstrong_ew_12oh_kk04 gt -900),nindx)

    mass = tancillary[indx].kcorr_mass ; Salpeter
    masserr = mass*0.0
    oh = tdust[indx].zstrong_ew_12oh_kk04
    oherr = tdust[indx].zstrong_ew_12oh_kk04_err
    
    xtitle = 'log M [M'+sunsymbol()+']  [K-correct, Salpeter]'
    ytitle = '12 + log (O/H)_{EW} [KK04]'

    xrange = [7.8,11.8]
    yrange = [8.0,9.5]

    sdss_lineplot, mass, oh, masserr, oherr, plottype=1, postscript=postscript, $
      xtitle=xtitle, ytitle=ytitle, xrange=xrange, yrange=yrange, position=pos, $
      charsize=1.8

    result = im_medxbin(mass,oh,0.1,minx=8.57,minpts=100L);,/verbose)
    coeff = poly_fit(result.binctr,result.medy,2)
    massaxis = findgen((max(result.binctr)-min(result.binctr))/0.01+1)*0.01+min(result.binctr)
    t04fit = poly(massaxis,coeff)
    
    im_symbols, 106, psize=1.1, thick=postthick, fill=1
    djs_oplot, result.binctr, result.medy, ps=8
    djs_oplot, massaxis, t04fit, line=0, thick=postthick, color='red'

    djs_oplot, tremonti_mass, tremonti_mass_oh_salp, line=2, thick=postthick, color='blue'
    niceprint, coeff, tremonti_mass_coeff

    if (not keyword_set(postscript)) then cc = get_kbrd(1)

; ---------------------------------------------------------------------------    
; MZ - Kauffmann masses (Figure 6)
; ---------------------------------------------------------------------------    
    
    pagemaker, nx=1, ny=1, xspace=0, yspace=0.0, width=6.5, height=6.5, $
      xmargin=[1.5,0.5], ymargin=[0.9,1.1], xpage=8.5, ypage=8.5, $
      position=pos, /normal

    indx = where((tancillary.kauffmann_mass gt -900.0) and (tancillary.tremonti_oh gt -900),nindx)

    mass = tancillary[indx].kauffmann_mass - imf_kroupa_to_salpeter ; Salpeter-->Kroupa
    masserr = tancillary[indx].kauffmann_mass_err
    oh = tancillary[indx].tremonti_oh
    oherr = tancillary[indx].tremonti_oh_err
    
    xtitle = 'log M [M'+sunsymbol()+']  [Kauffmann, Kroupa]'
    ytitle = '12 + log (O/H) [Tremonti]'

    xrange = [7.8,11.8]
    yrange = [8.0,9.5]

    sdss_lineplot, mass, oh, masserr, oherr, plottype=1, postscript=postscript, $
      xtitle=xtitle, ytitle=ytitle, xrange=xrange, yrange=yrange, position=pos, $
      charsize=1.8

    result = im_medxbin(mass,oh,0.1,minx=8.57,minpts=100L);,/verbose)
    coeff = poly_fit(result.binctr,result.medy,2)
    massaxis = findgen((max(result.binctr)-min(result.binctr))/0.01+1)*0.01+min(result.binctr)
    t04fit = poly(massaxis,coeff)
    
    im_symbols, 106, psize=1.1, thick=postthick, fill=1
    djs_oplot, result.binctr, result.medy, ps=8
;   djs_oplot, result.binctr, result.sigy95, line=0, thick=postthick
;   djs_oplot, result.binctr, result.sigy84, line=0, thick=postthick
;   djs_oplot, result.binctr, result.sigy16, line=0, thick=postthick
;   djs_oplot, result.binctr, result.sigy05, line=0, thick=postthick
    djs_oplot, massaxis, t04fit, line=0, thick=postthick, color='red'

    xchristy = [8.57,8.67,8.76,8.86,8.96,9.06,9.16,9.26,9.36,9.46,9.57,9.66,9.76,9.86,9.96,10.06,$
                10.16,10.26,10.36,10.46,10.56,10.66,10.76,10.86,10.95,11.05,11.15,11.25]
    ychristy = [8.44,8.48,8.57,8.61,8.63,8.66,8.68,8.71,8.74,8.78,8.82,8.84,8.87,8.90,8.94,$
                8.97,8.99,9.01,9.03,9.05,9.07,9.08,9.09,9.10,9.11,9.11,9.12,9.12]
    im_symbols, 108, psize=1.0, thick=postthick2, fill=1, color=djs_icolor('grey')
    djs_oplot, xchristy, ychristy, ps=8

    djs_oplot, tremonti_mass, tremonti_mass_oh, line=2, thick=postthick, color='blue'
    niceprint, coeff, tremonti_mass_coeff

;   im_openclose, postscript=postscript, /close    
    if (not keyword_set(postscript)) then cc = get_kbrd(1)

; ---------------------------------------------------------------------------    
; MZ - Kcorrect masses (Figure 6)
; ---------------------------------------------------------------------------    
    
;   psname = 't04_mz_kcorrect'
;   im_openclose, psname, postscript=postscript, xsize=8.5, ysize=8.5, encapsulated=encapsulated

    pagemaker, nx=1, ny=1, xspace=0, yspace=0.0, width=6.5, height=6.5, $
      xmargin=[1.5,0.5], ymargin=[0.9,1.1], xpage=8.5, ypage=8.5, $
      position=pos, /normal

    indx = where((tancillary.kcorr_mass gt -900.0) and (tancillary.tremonti_oh gt -900),nindx)

    mass = tancillary[indx].kcorr_mass - imf_kroupa_to_salpeter ; Salpeter-->Kroupa
    masserr = mass*0.0
    oh = tancillary[indx].tremonti_oh
    oherr = tancillary[indx].tremonti_oh_err
    
    xtitle = 'log M [M'+sunsymbol()+']  [K-correct, Kroupa]'
    ytitle = '12 + log (O/H) [Tremonti]'

    xrange = [7.8,11.8]
    yrange = [8.0,9.5]

    sdss_lineplot, mass, oh, masserr, oherr, plottype=1, postscript=postscript, $
      xtitle=xtitle, ytitle=ytitle, xrange=xrange, yrange=yrange, position=pos, $
      charsize=1.8

    result = im_medxbin(mass,oh,0.1,minx=8.57,minpts=100L);,/verbose)
    coeff = poly_fit(result.binctr,result.medy,2)
    massaxis = findgen((max(result.binctr)-min(result.binctr))/0.01+1)*0.01+min(result.binctr)
    t04fit = poly(massaxis,coeff)
    
    im_symbols, 106, psize=1.1, thick=postthick, fill=1
    djs_oplot, result.binctr, result.medy, ps=8
    djs_oplot, massaxis, t04fit, line=0, thick=postthick, color='red'

    xchristy = [8.57,8.67,8.76,8.86,8.96,9.06,9.16,9.26,9.36,9.46,9.57,9.66,9.76,9.86,9.96,10.06,$
                10.16,10.26,10.36,10.46,10.56,10.66,10.76,10.86,10.95,11.05,11.15,11.25]
    ychristy = [8.44,8.48,8.57,8.61,8.63,8.66,8.68,8.71,8.74,8.78,8.82,8.84,8.87,8.90,8.94,$
                8.97,8.99,9.01,9.03,9.05,9.07,9.08,9.09,9.10,9.11,9.11,9.12,9.12]
    im_symbols, 108, psize=1.0, thick=postthick2, fill=1, color=djs_icolor('grey')
    djs_oplot, xchristy, ychristy, ps=8

    djs_oplot, tremonti_mass, tremonti_mass_oh, line=2, thick=postthick, color='blue'
    niceprint, coeff, tremonti_mass_coeff

;   im_openclose, postscript=postscript, /close    
    if (not keyword_set(postscript)) then cc = get_kbrd(1)

; ---------------------------------------------------------------------------    
; Kcorrect vs Bell & de Jong
; ---------------------------------------------------------------------------    
    
;   psname = 't04_kcorrect_vs_BdJ'
;   im_openclose, psname, postscript=postscript, xsize=8.5, ysize=8.5, encapsulated=encapsulated

    pagemaker, nx=1, ny=1, xspace=0, yspace=0.0, width=6.5, height=6.5, $
      xmargin=[1.5,0.5], ymargin=[0.9,1.1], xpage=8.5, ypage=8.5, $
      position=pos, /normal

    indx = where((tancillary.kcorr_mass gt -900.0) and (tancillary.mass_sdss_gr_r gt -900),nindx)

    x = tancillary[indx].kcorr_mass ; Salpeter
    xerr = mass*0.0
    y = tancillary[indx].mass_sdss_gr_r ; Salpeter
    yerr = tancillary[indx].mass_sdss_gr_r_err

    stats = im_stats(y-x)
    xstr = strtrim(string(stats.median_rej,format='(F12.2)'),2)+' ('+$
      strtrim(string(stats.mean_rej,format='(F12.2)'),2)+$
      '\pm'+strtrim(string(stats.sigma_rej,format='(F12.2)'),2)+')'

    xtitle = 'log M [M'+sunsymbol()+']  [K-correct, Salpeter]'
    ytitle = 'log M [M'+sunsymbol()+']  [Bell & de Jong, Salpeter]'

    xrange = [7.8,11.8]
    yrange = xrange

    sdss_lineplot, x, y, xerr, yerr, plottype=1, postscript=postscript, $
      xtitle=xtitle, ytitle=ytitle, xrange=xrange, yrange=yrange, position=pos, $
      charsize=1.8
    djs_oplot, !x.crange, !y.crange, line=0, thick=postthick
    legend, textoidl('\Delta(log M) = '+xstr), /left, /top, box=0, charsize=1.8, $
      charthick=postthick

    sixlin, x, y, a, siga, b, sigb
    coeff = [a[2],b[2]] & coeff_err = [siga[2],sigb[2]]
    yfit = poly(massaxis,coeff)
    djs_oplot, massaxis, yfit, line=0, thick=postthick, color='red'

    niceprint, coeff
    
;   im_openclose, postscript=postscript, /close    
    if (not keyword_set(postscript)) then cc = get_kbrd(1)

; ---------------------------------------------------------------------------    
; Kcorrect vs Kauffmann
; ---------------------------------------------------------------------------    
    
;   psname = 't04_kcorrect_vs_kauffmann'
;   im_openclose, psname, postscript=postscript, xsize=8.5, ysize=8.5, encapsulated=encapsulated

    pagemaker, nx=1, ny=1, xspace=0, yspace=0.0, width=6.5, height=6.5, $
      xmargin=[1.5,0.5], ymargin=[0.9,1.1], xpage=8.5, ypage=8.5, $
      position=pos, /normal

    indx = where((tancillary.kcorr_mass gt -900.0) and (tancillary.kauffmann_mass gt -900),nindx)

    x = tancillary[indx].kcorr_mass - imf_kroupa_to_salpeter ; Salpeter-->Kroupa
    xerr = mass*0.0
    y = tancillary[indx].kauffmann_mass - imf_kroupa_to_salpeter
    yerr = tancillary[indx].kauffmann_mass_err

    stats = im_stats(y-x)
    xstr = strtrim(string(stats.median_rej,format='(F12.2)'),2)+' ('+$
      strtrim(string(stats.mean_rej,format='(F12.2)'),2)+$
      '\pm'+strtrim(string(stats.sigma_rej,format='(F12.2)'),2)+')'

    xtitle = 'log M [M'+sunsymbol()+']  [K-correct, Kroupa]'
    ytitle = 'log M [M'+sunsymbol()+']  [Kauffmann, Kroupa]'

    xrange = [7.8,11.8]
    yrange = xrange

    sdss_lineplot, x, y, xerr, yerr, plottype=1, postscript=postscript, $
      xtitle=xtitle, ytitle=ytitle, xrange=xrange, yrange=yrange, position=pos, $
      charsize=1.8
    djs_oplot, !x.crange, !y.crange, line=0, thick=postthick
    legend, textoidl('\Delta(log M) = '+xstr), /left, /top, box=0, charsize=1.8, $
      charthick=postthick

    sixlin, x, y, a, siga, b, sigb
    coeff = [a[2],b[2]] & coeff_err = [siga[2],sigb[2]]
    yfit = poly(massaxis,coeff)
    djs_oplot, massaxis, yfit, line=0, thick=postthick, color='red'

    niceprint, coeff
    
;   im_openclose, postscript=postscript, /close
    if (not keyword_set(postscript)) then cc = get_kbrd(1)
    
    if keyword_set(postscript) then dfpsclose

stop    
    
return
end
    
