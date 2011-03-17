pro zintegrated, gradients, samplehii, atlasdust, atlasnodust, postscript=postscript, $
  blackwhite=blackwhite, encapsulated=encapsulated, paper=paper, cleanpng=cleanpng, $
  _extra=extra
; jm05mar11uofa - originally written 
; jm05sep09uofa - revamped
; jm06mar26uofa - major re-write    

    if keyword_set(paper) then begin
       postscript = 1L
       encapsulated = 1L
       cmyk = 1L
    endif

    htmlbase = 'zintegrated'

    datapath = atlas_path(/projects)+'zintegrated/'
    html_path = atlas_path(/web)+'analysis/'
    pspath = atlas_path(/web)+'analysis/'+htmlbase+'/'
    paperpath = atlas_path(/papers)+'zintegrated/FIG_ZINTEGRATED'

    rr25_int = 0.4
    
; setup some plotting variables

    @'xyrange_zintegrated'

    if keyword_set(postscript) then begin
       postthick = 5.0
       postthick2 = 8.0
       postthick3 = 4.0
       defaultcolor = 'black'
    endif else begin
       blackwhite = 0L
       postthick = 2.0
       postthick2 = postthick
       postthick3 = postthick
       im_window, 0, xratio=0.5, /square
       defaultcolor = 'white'
    endelse

    if keyword_set(blackwhite) then begin
       pspath = pspath+'blackwhite/'
       suffix = 'blackwhite'
       color = 0L
    endif else begin
       suffix = ''
       color = 0L
    endelse
    
    psizeall = 1.1
    psize = 1.4
    psize2 = 1.8
    psize3 = 2.8
    psize4 = 2.0
    psize5 = 1.1
    hiipsize = 0.7
    digpsize = 1.5

    if keyword_set(blackwhite) then obscolor = 'black' else obscolor = 'red'
    if keyword_set(blackwhite) then corcolor = 'black' else corcolor = 'blue'
    if keyword_set(blackwhite) then ewcolor = 'black' else ewcolor = 'dark green'
    if keyword_set(blackwhite) then kkpcolor = 'dark gray' else kkpcolor = 'dark gray'

    if keyword_set(blackwhite) then intcolor = 'dark gray' else intcolor = 'dark green'
    if keyword_set(blackwhite) then digcolor = 'dark gray' else digcolor = 'blue' ; 'dark orchid'
    if keyword_set(blackwhite) then hiicolor = 'gray' else hiicolor = 'gray'
    if keyword_set(blackwhite) then atlascolor = 'black' else atlascolor = 'red'

; clean up old PS and PNG files

    if keyword_set(cleanpng) then begin
       splog, 'Deleting all PNG and PS files in '+pspath
       spawn, ['rm -f '+pspath+'*.png*'], /sh
       spawn, ['rm -f '+pspath+'*ps'], /sh
    endif

; restore the abundance gradient fitting results, read a limited
; sample of HII regions, the complete integrated spectral atlas, and
; some DIG measurements 

    if (n_elements(gradients) eq 0L) then gradients = mrdfits(datapath+'zintegrated_gradients.fits.gz',1,/silent)
    if (n_elements(samplehii) eq 0L) then samplehii = read_hii_regions(/limitedrefs)
    if (n_elements(atlasdust) eq 0L) then atlasdust = read_integrated(atlasnodust=atlasnodust)

    kkp = rsex(datapath+'99kkp.dat') ; read the data from KKP99

    h03_m33 = rsex('/home/ioannis/catalogs/03hoopes/03hoopes_m33.dat')
    h03_m51 = rsex('/home/ioannis/catalogs/03hoopes/03hoopes_m51.dat')
    h03 = struct_append(h03_m33,h03_m51)
    w97 = rsex('/home/ioannis/catalogs/97wang/97wang.dat')
    g99_niisii = rsex('/home/ioannis/catalogs/99galarza/99galarza_dig_niisii.dat')
    g99_oiioiii = rsex('/home/ioannis/catalogs/99galarza/99galarza_dig_oiioiii.dat')

    ngalaxy = n_elements(gradients)
    galaxy = strtrim(gradients.galaxy,2)

; ---------------------------------------------------------------------------
; [S II]/Ha vs [N II]/Ha - Full Integrated Sample
; ---------------------------------------------------------------------------

    psname = 'all_siiha_vs_niiha'
    im_openclose, pspath+psname, postscript=postscript, xsize=8.5, ysize=8.0, encapsulated=encapsulated, cmyk=cmyk

    pagemaker, nx=1, ny=1, xspace=0, yspace=0.0, width=6.6, height=6.6, $
      xmargin=[1.5,0.4], ymargin=[0.3,1.1], xpage=8.5, ypage=8.0, $
      position=pos, /normal

; integrated spectra, HII galaxies    

    sf = where(strtrim(atlasdust.bpt_class,2) eq 'HII')
    sfatlasdust = atlasdust[sf]
    sfatlasnodust = atlasnodust[sf]
    
    cutindx = cmset_op(strtrim(sfatlasdust.galaxy,2),'and',/not2,strtrim(gradients.galaxy,2),/index)
    lineratio, sfatlasnodust[cutindx], 'SII_6716', 'H_ALPHA', 'NII_6584', 'H_ALPHA', $
      x, xerr, y, yerr, snrcut=3.0, index=indx, nindex=nindx

    cutindxint = cmset_op(strtrim(sfatlasdust.galaxy,2),'and',strtrim(gradients.galaxy,2),/index)
    lineratio, sfatlasnodust[cutindxint], 'SII_6716', 'H_ALPHA', 'NII_6584', 'H_ALPHA', $
      xint, xerrint, yint, yerrint, snrcut=3.0, index=indxint, nindex=nindxint
    
    xtitle = 'log ([S II] \lambda6716/H\alpha)_{cor}'
    ytitle = 'log ([N II] \lambda6584/H\alpha)_{cor}'

    xrange = siiharange2
    yrange = niiharange3

    good = where((samplehii.sii_6716_h_alpha gt -900.0) and (samplehii.nii_6584_h_alpha gt -900.0))
    xregion = samplehii[good].sii_6716_h_alpha & xerrregion = samplehii[good].sii_6716_h_alpha_err
    yregion = samplehii[good].nii_6584_h_alpha & yerrregion = samplehii[good].nii_6584_h_alpha_err

    atlas1d_lineplot, x, y, xerr, yerr, plottype=1, postscript=postscript, $
      xtitle=xtitle, ytitle=ytitle, xrange=xrange, yrange=yrange, legendtype=0, $
      /right, /top, position=pos[*,0], atlaspsize=psizeall, atlassym=104, atlascolor=atlascolor, $
      atlasfill=1, symthick=postthick, postthick=postthick, charsize=charsize_9, blackwhite=blackwhite
    legend, '(a)', /left, /top, box=0, charsize=singlecharsize, charthick=postthick

; overplot the HII regions

    im_symbols, 108, fill=1L, psize=hiipsize, color=fsc_color(hiicolor,25), thick=postthick2
    djs_oplot, xregion, yregion, psym=8

; overplot the small sample    
    
    im_symbols, 106, fill=1L, psize=1.7, color=fsc_color(intcolor,73), thick=postthick2
    djs_oplot, xint, yint, psym=8

; overplot DIG points
 
    djs_oplot, alog10(h03.siiha), alog10(h03.niiha), ps=7, syms=digpsize, $
      thick=postthick, color=fsc_color(digcolor,72)

    djs_oplot, g99_niisii.siiha, g99_niisii.niiha, ps=7, syms=digpsize, $
      thick=postthick, color=fsc_color(digcolor,72)

;   gw97 = where((w97.siiha gt -900) and (w97.niiha gt -900),ngw97)
;   djs_oplot, alog10(w97[gw97].siiha), alog10(w97[gw97].niiha), ps=7, syms=digpsize, $
;     thick=postthick, color=fsc_color('dark orchid',75)
    
;   plot_kewley_grids, model=8L, plotnumber=27L, /overplot, labeltype=1L, /siioffset

    im_openclose, postscript=postscript, /close

; ---------------------------------------------------------------------------
; [O III]/Hb vs [O II]/Hb - Full Integrated Sample
; ---------------------------------------------------------------------------

    psname = 'all_oiiihb_vs_oiihb'
    im_openclose, pspath+psname, postscript=postscript, xsize=8.5, ysize=8.0, encapsulated=encapsulated, cmyk=cmyk

    pagemaker, nx=1, ny=1, xspace=0, yspace=0.0, width=6.6, height=6.6, $
      xmargin=[1.5,0.4], ymargin=[0.3,1.1], xpage=8.5, ypage=8.0, $
      position=pos, /normal

; integrated spectra, HII galaxies    

    sf = where(strtrim(atlasdust.bpt_class,2) eq 'HII')
    sfatlasdust = atlasdust[sf]
    sfatlasnodust = atlasnodust[sf]
    
    cutindx = cmset_op(strtrim(sfatlasdust.galaxy,2),'and',/not2,strtrim(gradients.galaxy,2),/index)
    lineratio, sfatlasnodust[cutindx], 'OIII_5007', 'H_BETA', 'OII_3727', 'H_BETA', $
      x, xerr, y, yerr, snrcut=3.0, index=indx, nindex=nindx

    cutindxint = cmset_op(strtrim(sfatlasdust.galaxy,2),'and',strtrim(gradients.galaxy,2),/index)
    lineratio, sfatlasnodust[cutindxint], 'OIII_5007', 'H_BETA', 'OII_3727', 'H_BETA', $
      xint, xerrint, yint, yerrint, snrcut=3.0, index=indxint, nindex=nindxint
    
    xtitle = 'log ([O III] \lambda5007/H\beta)_{cor}'
    ytitle = 'log ([O II] \lambda3727/H\beta)_{cor}'

    xrange = oiiihbrange2
    yrange = oiihbrange2

    good = where((samplehii.oii_h_beta gt -900.0) and (samplehii.oiii_5007_h_beta gt -900.0))
    xregion = samplehii[good].oiii_5007_h_beta & xerrregion = samplehii[good].oiii_5007_h_beta_err
    yregion = samplehii[good].oii_h_beta & yerrregion = samplehii[good].oii_h_beta_err

    atlas1d_lineplot, x, y, xerr, yerr, plottype=1, postscript=postscript, $
      xtitle=xtitle, ytitle=ytitle, xrange=xrange, yrange=yrange, legendtype=0, $
      /right, /top, position=pos[*,0], atlaspsize=psizeall, atlassym=104, atlascolor=atlascolor, $
      atlasfill=1, symthick=postthick, postthick=postthick, charsize=charsize_9, blackwhite=blackwhite
    legend, '(b)', /left, /top, box=0, charsize=singlecharsize, charthick=postthick

; overplot the HII regions

    im_symbols, 108, fill=1L, psize=hiipsize, color=fsc_color(hiicolor,25), thick=postthick2
    djs_oplot, xregion, yregion, psym=8

; overplot the small sample    
    
    im_symbols, 106, fill=1L, psize=1.7, color=fsc_color(intcolor,73), thick=postthick2
    djs_oplot, xint, yint, psym=8

; overplot DIG points
 
    djs_oplot, g99_oiioiii.oiiihb, g99_oiioiii.oiiihb+g99_oiioiii.oiioiii, $
      ps=7, syms=digpsize, thick=postthick, color=fsc_color(digcolor,72)
    
;   gg99 = where((g99.oiihb gt -900) and (g99.oiiihb gt -900),ngg99)
;   djs_oplot, alog10(g99[gg99].oiiihb), alog10(g99[gg99].oiihb), ps=7, syms=digpsize, $
;     thick=postthick, color=fsc_color('dark orchid',75)

;   plot_kewley_grids, model=8L, plotnumber=29L, labeltype=1L, /siioffset

    im_openclose, postscript=postscript, /close

; ---------------------------------------------------------------------------
; Gradient Slope vs Characteristic Radius
; ---------------------------------------------------------------------------

    psname = 'gradient_vs_rr25_char'
    im_openclose, pspath+psname, postscript=postscript, xsize=8.5, ysize=8.0, encapsulated=encapsulated

    pagemaker, nx=1, ny=1, xspace=0, yspace=0.0, width=6.6, height=6.6, $
      xmargin=[1.5,0.4], ymargin=[0.3,1.1], xpage=8.5, ypage=8.0, $
      position=pos, /normal

    xtitle = textoidl('Gradient Slope [dex \rho_{25}^{-1}]')
    ytitle = textoidl('Characteristic \rho / \rho_{25}')

    xrange = sloperange1
    yrange = [0.0,1.1]

    plot, [0], [0], /nodata, xthick=postthick, ythick=postthick, $
      charsize=charsize_9, charthick=postthick, thick=postthick, $
      xtitle=xtitle, ytitle=ytitle, xstyle=3, color=djs_icolor(defaultcolor), $
      ystyle=1, xrange=xrange, yrange=yrange, position=pos
    oplot, !x.crange, [0.4,0.4], line=2, thick=postthick

    im_symbols, 108, psize=psize2, fill=1, thick=postthick2, color=fsc_color(corcolor,74)
    m91good = where((gradients.hii_m91_slope[0] gt -900.0) and (gradients.int_cor_rr25_m91 gt -900.0),nm91good)
;   niceprint, gradients[m91good].galaxy, gradients[m91good].hii_m91_slope[0], gradients[m91good].int_cor_rr25_m91
    stats = im_stats(gradients[m91good].int_cor_rr25_m91,/verbose)
    oplot, gradients[m91good].hii_m91_slope[0], gradients[m91good].int_cor_rr25_m91, ps=8

    im_symbols, 106, psize=psize2, fill=0, thick=postthick2, color=fsc_color(obscolor,75)
    pt05good = where((gradients.hii_pt05_slope[0] gt -900.0) and (gradients.int_cor_rr25_pt05 gt -900.0),npt05good)
;   niceprint, gradients[pt05good].galaxy, gradients[pt05good].hii_pt05_slope[0], gradients[pt05good].int_cor_rr25_pt05
    stats = im_stats(gradients[pt05good].int_cor_rr25_pt05,/verbose)
    oplot, gradients[pt05good].hii_pt05_slope[0], gradients[pt05good].int_cor_rr25_pt05, ps=8

;   niceprint, gradients[m91good].galaxy, gradients[m91good].hii_m91_slope[0], gradients[m91good].int_cor_rr25_m91, $
;     gradients[pt05good].galaxy, gradients[pt05good].hii_pt05_slope[0], gradients[pt05good].int_cor_rr25_pt05

    im_openclose, postscript=postscript, /close

; ---------------------------------------------------------------------------
; Slope [M91] vs Slope [PT05]
; ---------------------------------------------------------------------------

    psname = 'slope_m91_vs_slope_pt05'
    im_openclose, pspath+psname, postscript=postscript, xsize=8.5, ysize=8.0, encapsulated=encapsulated, cmyk=cmyk

    pagemaker, nx=1, ny=1, xspace=0, yspace=0.0, width=6.6, height=6.6, $
      xmargin=[1.5,0.4], ymargin=[0.3,1.1], xpage=8.5, ypage=8.0, $
      position=pos, /normal

    xtitle = 'M91 Gradient Slope [dex \rho_{25}^{-1}]'
    ytitle = 'PT05 Gradient Slope [dex \rho_{25}^{-1}]'

    xrange = sloperange2
    yrange = xrange

    x = gradients.hii_m91_slope[0]
    xerr = gradients.hii_m91_slope[1]
    y = gradients.hii_pt05_slope[0]
    yerr = gradients.hii_pt05_slope[1]

    resid = abs(y-x) & srt = sort(resid)
;   niceprint, gradients[srt].galaxy, resid[srt]
    
    atlas1d_lineplot, x, y, xerr, yerr, plottype=2, postscript=postscript, $
      xtitle=xtitle, ytitle=ytitle, xrange=xrange, yrange=yrange, legendtype=0, $
      /right, /top, position=pos[*,0], atlaspsize=psize2, atlascolor=defaultcolor, $
      atlasfill=1, symthick=postthick2, postthick=postthick
    oplot, !x.crange, !y.crange, line=0, thick=postthick

; plot an inset comparing the characteristic abundances

    xtitle = 'Characteristic 12+log(O/H)_{M91}'
    ytitle = 'Characteristic 12+log(O/H)_{PT05}'

    xrange = ohrange3
    yrange = xrange

    x = gradients.hii_m91_log12oh_char[0]
    xerr = gradients.hii_m91_log12oh_char[1]
    y = gradients.hii_pt05_log12oh_char[0]
    yerr = gradients.hii_pt05_log12oh_char[1]

    resid = abs(y-x) & srt = sort(resid)
;   niceprint, gradients[srt].galaxy, resid[srt]
    
    xhoff = 0.45 & yhoff = 0.08 & xsize = 0.28 & ysize = 0.23
    insetpos = [pos[0]+xhoff,pos[1]+yhoff,pos[0]+xhoff+xsize,pos[1]+yhoff+ysize]

    atlas1d_lineplot, x, y, xerr, yerr, plottype=1, postscript=postscript, $
      xtitle=xtitle, ytitle=ytitle, xrange=xrange, yrange=yrange, legendtype=0, $
      /right, /top, position=insetpos[*,0], atlaspsize=psize5, atlascolor=defaultcolor, $
      atlasfill=1, symthick=postthick2, postthick=postthick, /noerase, $
      charsize=charsize_1, charthick=postthick3
    oplot, !x.crange, !y.crange, line=0, thick=postthick
    oplot, !x.crange, !y.crange-median(resid), line=2, thick=postthick

    im_openclose, postscript=postscript, /close

; ---------------------------------------------------------------------------
; 12+log(O/H) Residuals [Int. minus Char.] vs various quantities - M91
; ---------------------------------------------------------------------------

    psname = '12oh_m91_residuals'
    im_openclose, pspath+psname, postscript=postscript, xsize=8.5, ysize=8.5, encapsulated=encapsulated

    pagemaker, nx=2, ny=2, yspace=1.0, xspace=0.0, width=3.2*[1,1], height=3.0*[1,1], $
      xmargin=[1.1,0.5], ymargin=[0.4,1.1], xpage=8.5, ypage=8.5, position=pos, /normal

    ytitle = textoidl('\Deltalog(O/H)_{M91}')
    yrange = ohresidrange

; residuals - observed
    
    good = where((gradients.int_obs_log12oh_m91[0] gt -900.0) and (gradients.hii_m91_log12oh_char[0] gt -900.0) and $
      (gradients.int_log12oh_lower_limit eq 0),ngood)
    oh_char_obs = gradients[good].hii_m91_log12oh_char[0]
    oh_int_obs = gradients[good].int_obs_log12oh_m91[0]
    delta_oh_obs = oh_int_obs - oh_char_obs
    slope_obs = gradients[good].hii_m91_slope[0]
    ebv_obs = gradients[good].int_ebv[0]
    incl_obs = gradients[good].incl
    hasb_obs = gradients[good].hasb

    limit = where((gradients.int_obs_log12oh_m91[0] gt -900.0) and (gradients.hii_m91_log12oh_char[0] gt -900.0) and $
      (gradients.int_log12oh_lower_limit eq 1),nlimit)
    oh_char_obs_limit = gradients[limit].hii_m91_log12oh_char[0]
    oh_int_obs_limit = gradients[limit].int_obs_log12oh_m91[0]
    delta_oh_obs_limit = oh_int_obs_limit - oh_char_obs_limit
    slope_obs_limit = gradients[limit].hii_m91_slope[0]
    ebv_obs_limit = gradients[limit].int_ebv[0]
    incl_obs_limit = gradients[limit].incl
    hasb_obs_limit = gradients[limit].hasb
    
; residuals - corrected
    
    good = where((gradients.int_cor_log12oh_m91[0] gt -900.0) and (gradients.hii_m91_log12oh_char[0] gt -900.0) and $
      (gradients.int_log12oh_lower_limit eq 0),ngood)
    oh_char_cor = gradients[good].hii_m91_log12oh_char[0]
    oh_int_cor = gradients[good].int_cor_log12oh_m91[0]
    delta_oh_cor = oh_int_cor - oh_char_cor
    slope_cor = gradients[good].hii_m91_slope[0]
    ebv_cor = gradients[good].int_ebv[0]
    incl_cor = gradients[good].incl
    hasb_cor = gradients[good].hasb

    limit = where((gradients.int_cor_log12oh_m91[0] gt -900.0) and (gradients.hii_m91_log12oh_char[0] gt -900.0) and $
      (gradients.int_log12oh_lower_limit eq 1),nlimit)
    oh_char_cor_limit = gradients[limit].hii_m91_log12oh_char[0]
    oh_int_cor_limit = gradients[limit].int_cor_log12oh_m91[0]
    delta_oh_cor_limit = oh_int_cor_limit - oh_char_cor_limit
    slope_cor_limit = gradients[limit].hii_m91_slope[0]
    ebv_cor_limit = gradients[limit].int_ebv[0]
    incl_cor_limit = gradients[limit].incl
    hasb_cor_limit = gradients[limit].hasb
    
; residuals - EW
    
    good = where((gradients.int_ew_log12oh_m91[0] gt -900.0) and (gradients.hii_m91_log12oh_char[0] gt -900.0) and $
      (gradients.int_log12oh_lower_limit eq 0),ngood)
    oh_char_ew = gradients[good].hii_m91_log12oh_char[0]
    oh_int_ew = gradients[good].int_ew_log12oh_m91[0]
    delta_oh_ew = oh_int_ew - oh_char_ew
    slope_ew = gradients[good].hii_m91_slope[0]
    ebv_ew = gradients[good].int_ebv[0]
    incl_ew = gradients[good].incl
    hasb_ew = gradients[good].hasb

    limit = where((gradients.int_ew_log12oh_m91[0] gt -900.0) and (gradients.hii_m91_log12oh_char[0] gt -900.0) and $
      (gradients.int_log12oh_lower_limit eq 1),nlimit)
    oh_char_ew_limit = gradients[limit].hii_m91_log12oh_char[0]
    oh_int_ew_limit = gradients[limit].int_ew_log12oh_m91[0]
    delta_oh_ew_limit = oh_int_ew_limit - oh_char_ew_limit
    slope_ew_limit = gradients[limit].hii_m91_slope[0]
    ebv_ew_limit = gradients[limit].int_ebv[0]
    incl_ew_limit = gradients[limit].incl
    hasb_ew_limit = gradients[limit].hasb
    
; slope
    
    xrange = sloperange1
    xtitle = textoidl('Gradient Slope [dex \rho_{25}^{-1}]')

    plot, [0], [0], /nodata, xthick=postthick, ythick=postthick, $
      charsize=charsize_6, charthick=postthick, thick=postthick, $
      xtitle=xtitle, ytitle=ytitle, xstyle=3, ystyle=3, xrange=xrange, $
      yrange=yrange, position=pos[*,0], xtickinterval=0.5, color=djs_icolor(defaultcolor)
    oplot, !x.crange, [0,0], line=0, thick=postthick
    legend, '(a)', /left, /top, box=0, charsize=charsize_6, charthick=postthick
 
    im_symbols, 108, psize=psize2, fill=0, thick=postthick2, color=fsc_color(obscolor,75)
    oplot, slope_obs, delta_oh_obs, ps=8
    im_symbols, 112, psize=psize3, fill=0, thick=postthick2, color=fsc_color(obscolor,75)
    oplot, slope_obs_limit, delta_oh_obs_limit, ps=8

    im_symbols, 108, psize=psize2, fill=1, thick=postthick2, color=fsc_color(corcolor,74)
    oplot, slope_cor, delta_oh_cor, ps=8
    im_symbols, 112, psize=psize3, fill=0, thick=postthick2, color=fsc_color(corcolor,74)
    oplot, slope_cor_limit, delta_oh_cor_limit, ps=8

    im_symbols, 106, psize=psize2, fill=0, thick=postthick2, color=fsc_color(ewcolor,73)
    oplot, slope_ew, delta_oh_ew, ps=8
    im_symbols, 112, psize=psize3, fill=0, thick=postthick2, color=fsc_color(ewcolor,73)
    oplot, slope_ew_limit, delta_oh_ew_limit, ps=8

; overplot the KKP99 data

    im_symbols, 105, psize=psize4, fill=1, thick=postthick2, color=fsc_color(kkpcolor,71)
    djs_oplot, kkp.gradient, kkp.oh12_global-kkp.oh12_char, ps=8

; ebv
    
    xrange = ebvrange
    xtitle = 'E(B-V) [mag]'

    plot, [0], [0], /nodata, /noerase, xthick=postthick, ythick=postthick, $
      charsize=charsize_6, charthick=postthick, thick=postthick, $
      xtitle=xtitle, ytitle='', xstyle=3, ystyle=11, xrange=xrange, $
      yrange=yrange, position=pos[*,1], ytickname=replicate(' ',10), color=djs_icolor(defaultcolor)
    axis, /yaxis, yrange=yrange, ysty=3, charsize=charsize_6, charthick=postthick, $
      ythick=postthick, ytitle=textoidl(ytitle)
    oplot, !x.crange, [0,0], line=0, thick=postthick
    legend, '(b)', /left, /top, box=0, charsize=charsize_6, charthick=postthick
 
    im_symbols, 108, psize=psize2, fill=0, thick=postthick2, color=fsc_color(obscolor,75)
    oplot, ebv_obs, delta_oh_obs, ps=8
    im_symbols, 112, psize=psize3, fill=0, thick=postthick2, color=fsc_color(obscolor,75)
    oplot, ebv_obs_limit, delta_oh_obs_limit, ps=8

    im_symbols, 108, psize=psize2, fill=1, thick=postthick2, color=fsc_color(corcolor,74)
    oplot, ebv_cor, delta_oh_cor, ps=8
    im_symbols, 112, psize=psize3, fill=0, thick=postthick2, color=fsc_color(corcolor,74)
    oplot, ebv_cor_limit, delta_oh_cor_limit, ps=8

    im_symbols, 106, psize=psize2, fill=0, thick=postthick2, color=fsc_color(ewcolor,73)
    oplot, ebv_ew, delta_oh_ew, ps=8
    im_symbols, 112, psize=psize3, fill=0, thick=postthick2, color=fsc_color(ewcolor,73)
    oplot, ebv_ew_limit, delta_oh_ew_limit, ps=8

; inclination
    
    xrange = inclrange
    xtitle = 'Inclination [degrees]'

    plot, [0], [0], /nodata, /noerase, xthick=postthick, ythick=postthick, $
      charsize=charsize_6, charthick=postthick, thick=postthick, $
      xtitle=xtitle, ytitle=ytitle, xstyle=3, ystyle=3, xrange=xrange, $
      yrange=yrange, position=pos[*,2], color=djs_icolor(defaultcolor)
    oplot, !x.crange, [0,0], line=0, thick=postthick
    legend, '(c)', /left, /top, box=0, charsize=charsize_6, charthick=postthick
 
    im_symbols, 108, psize=psize2, fill=0, thick=postthick2, color=fsc_color(obscolor,75)
    oplot, incl_obs, delta_oh_obs, ps=8
    im_symbols, 112, psize=psize3, fill=0, thick=postthick2, color=fsc_color(obscolor,75)
    oplot, incl_obs_limit, delta_oh_obs_limit, ps=8

    im_symbols, 108, psize=psize2, fill=1, thick=postthick2, color=fsc_color(corcolor,74)
    oplot, incl_cor, delta_oh_cor, ps=8
    im_symbols, 112, psize=psize3, fill=0, thick=postthick2, color=fsc_color(corcolor,74)
    oplot, incl_cor_limit, delta_oh_cor_limit, ps=8

    im_symbols, 106, psize=psize2, fill=0, thick=postthick2, color=fsc_color(ewcolor,73)
    oplot, incl_ew, delta_oh_ew, ps=8
    im_symbols, 112, psize=psize3, fill=0, thick=postthick2, color=fsc_color(ewcolor,73)
    oplot, incl_ew_limit, delta_oh_ew_limit, ps=8

; SB(H-alpha)
    
    xrange = hasbrange
    xtitle = textoidl('log [\Sigma(H\alpha)] [erg s^{-1} pc^{-2}]')

    plot, [0], [0], /nodata, /noerase, xthick=postthick, ythick=postthick, $
      charsize=charsize_6, charthick=postthick, thick=postthick, $
      xtitle=xtitle, ytitle='', xstyle=3, ystyle=11, xrange=xrange, $
      yrange=yrange, position=pos[*,3], ytickname=replicate(' ',10), color=djs_icolor(defaultcolor)
    axis, /yaxis, yrange=yrange, ysty=3, charsize=charsize_6, charthick=postthick, $
      ythick=postthick, ytitle=textoidl(ytitle)
    oplot, !x.crange, [0,0], line=0, thick=postthick
    legend, '(d)', /left, /top, box=0, charsize=charsize_6, charthick=postthick
 
    im_symbols, 108, psize=psize2, fill=0, thick=postthick2, color=fsc_color(obscolor,75)
    oplot, hasb_obs, delta_oh_obs, ps=8
    im_symbols, 112, psize=psize3, fill=0, thick=postthick2, color=fsc_color(obscolor,75)
    oplot, hasb_obs_limit, delta_oh_obs_limit, ps=8

    im_symbols, 108, psize=psize2, fill=1, thick=postthick2, color=fsc_color(corcolor,74)
    oplot, hasb_cor, delta_oh_cor, ps=8
    im_symbols, 112, psize=psize3, fill=0, thick=postthick2, color=fsc_color(corcolor,74)
    oplot, hasb_cor_limit, delta_oh_cor_limit, ps=8

    im_symbols, 106, psize=psize2, fill=0, thick=postthick2, color=fsc_color(ewcolor,73)
    oplot, hasb_ew, delta_oh_ew, ps=8
    im_symbols, 112, psize=psize3, fill=0, thick=postthick2, color=fsc_color(ewcolor,73)
    oplot, hasb_ew_limit, delta_oh_ew_limit, ps=8

    im_openclose, postscript=postscript, /close

; ---------------------------------------------------------------------------
; 12+log(O/H) Residuals [Int. minus Char.] vs various quantities - PT05
; ---------------------------------------------------------------------------

    psname = '12oh_pt05_residuals'
    im_openclose, pspath+psname, postscript=postscript, xsize=8.5, ysize=8.5, encapsulated=encapsulated

    pagemaker, nx=2, ny=2, yspace=1.0, xspace=0.0, width=3.2*[1,1], height=3.0*[1,1], $
      xmargin=[1.1,0.5], ymargin=[0.4,1.1], xpage=8.5, ypage=8.5, position=pos, /normal

    ytitle = '\Deltalog(O/H)_{PT05}'
    yrange = ohresidrange

; residuals - observed
    
    good = where((gradients.int_obs_log12oh_pt05[0] gt -900.0) and (gradients.hii_pt05_log12oh_char[0] gt -900.0) and $
      (gradients.int_log12oh_lower_limit eq 0),ngood)
    oh_char_obs = gradients[good].hii_pt05_log12oh_char[0]
    oh_int_obs = gradients[good].int_obs_log12oh_pt05[0]
    delta_oh_obs = oh_int_obs - oh_char_obs
    slope_obs = gradients[good].hii_pt05_slope[0]
    ebv_obs = gradients[good].int_ebv[0]
    incl_obs = gradients[good].incl
    hasb_obs = gradients[good].hasb

    limit = where((gradients.int_obs_log12oh_pt05[0] gt -900.0) and (gradients.hii_pt05_log12oh_char[0] gt -900.0) and $
      (gradients.int_log12oh_lower_limit eq 1),nlimit)
    oh_char_obs_limit = gradients[limit].hii_pt05_log12oh_char[0]
    oh_int_obs_limit = gradients[limit].int_obs_log12oh_pt05[0]
    delta_oh_obs_limit = oh_int_obs_limit - oh_char_obs_limit
    slope_obs_limit = gradients[limit].hii_pt05_slope[0]
    ebv_obs_limit = gradients[limit].int_ebv[0]
    incl_obs_limit = gradients[limit].incl
    hasb_obs_limit = gradients[limit].hasb
    
; residuals - corrected
    
    good = where((gradients.int_cor_log12oh_pt05[0] gt -900.0) and (gradients.hii_pt05_log12oh_char[0] gt -900.0) and $
      (gradients.int_log12oh_lower_limit eq 0),ngood)
    oh_char_cor = gradients[good].hii_pt05_log12oh_char[0]
    oh_int_cor = gradients[good].int_cor_log12oh_pt05[0]
    delta_oh_cor = oh_int_cor - oh_char_cor
    slope_cor = gradients[good].hii_pt05_slope[0]
    ebv_cor = gradients[good].int_ebv[0]
    incl_cor = gradients[good].incl
    hasb_cor = gradients[good].hasb

    limit = where((gradients.int_cor_log12oh_pt05[0] gt -900.0) and (gradients.hii_pt05_log12oh_char[0] gt -900.0) and $
      (gradients.int_log12oh_lower_limit eq 1),nlimit)
    oh_char_cor_limit = gradients[limit].hii_pt05_log12oh_char[0]
    oh_int_cor_limit = gradients[limit].int_cor_log12oh_pt05[0]
    delta_oh_cor_limit = oh_int_cor_limit - oh_char_cor_limit
    slope_cor_limit = gradients[limit].hii_pt05_slope[0]
    ebv_cor_limit = gradients[limit].int_ebv[0]
    incl_cor_limit = gradients[limit].incl
    hasb_cor_limit = gradients[limit].hasb
    
; residuals - EW
    
    good = where((gradients.int_ew_log12oh_pt05[0] gt -900.0) and (gradients.hii_pt05_log12oh_char[0] gt -900.0) and $
      (gradients.int_log12oh_lower_limit eq 0),ngood)
    oh_char_ew = gradients[good].hii_pt05_log12oh_char[0]
    oh_int_ew = gradients[good].int_ew_log12oh_pt05[0]
    delta_oh_ew = oh_int_ew - oh_char_ew
    slope_ew = gradients[good].hii_pt05_slope[0]
    ebv_ew = gradients[good].int_ebv[0]
    incl_ew = gradients[good].incl
    hasb_ew = gradients[good].hasb

    limit = where((gradients.int_ew_log12oh_pt05[0] gt -900.0) and (gradients.hii_pt05_log12oh_char[0] gt -900.0) and $
      (gradients.int_log12oh_lower_limit eq 1),nlimit)
    oh_char_ew_limit = gradients[limit].hii_pt05_log12oh_char[0]
    oh_int_ew_limit = gradients[limit].int_ew_log12oh_pt05[0]
    delta_oh_ew_limit = oh_int_ew_limit - oh_char_ew_limit
    slope_ew_limit = gradients[limit].hii_pt05_slope[0]
    ebv_ew_limit = gradients[limit].int_ebv[0]
    incl_ew_limit = gradients[limit].incl
    hasb_ew_limit = gradients[limit].hasb
    
; slope
    
    xrange = sloperange1
    xtitle = textoidl('Gradient Slope [dex \rho_{25}^{-1}]')

    plot, [0], [0], /nodata, xthick=postthick, ythick=postthick, $
      charsize=charsize_6, charthick=postthick, thick=postthick, $
      xtitle=xtitle, ytitle=textoidl(ytitle), xstyle=3, ystyle=3, xrange=xrange, $
      yrange=yrange, position=pos[*,0], xtickinterval=0.5, color=djs_icolor(defaultcolor)
    oplot, !x.crange, [0,0], line=0, thick=postthick
    legend, '(a)', /left, /top, box=0, charsize=charsize_6, charthick=postthick
 
    im_symbols, 108, psize=psize2, fill=0, thick=postthick2, color=fsc_color(obscolor,75)
    oplot, slope_obs, delta_oh_obs, ps=8
    im_symbols, 112, psize=psize3, fill=0, thick=postthick2, color=fsc_color(obscolor,75)
    oplot, slope_obs_limit, delta_oh_obs_limit, ps=8

    im_symbols, 108, psize=psize2, fill=1, thick=postthick2, color=fsc_color(corcolor,74)
    oplot, slope_cor, delta_oh_cor, ps=8
    im_symbols, 112, psize=psize3, fill=0, thick=postthick2, color=fsc_color(corcolor,74)
    oplot, slope_cor_limit, delta_oh_cor_limit, ps=8

    im_symbols, 106, psize=psize2, fill=0, thick=postthick2, color=fsc_color(ewcolor,73)
    oplot, slope_ew, delta_oh_ew, ps=8
    im_symbols, 112, psize=psize3, fill=0, thick=postthick2, color=fsc_color(ewcolor,73)
    oplot, slope_ew_limit, delta_oh_ew_limit, ps=8

; ebv
    
    xrange = ebvrange
    xtitle = 'E(B-V) [mag]'

    plot, [0], [0], /nodata, /noerase, xthick=postthick, ythick=postthick, $
      charsize=charsize_6, charthick=postthick, thick=postthick, $
      xtitle=xtitle, ytitle='', xstyle=3, ystyle=11, xrange=xrange, $
      yrange=yrange, position=pos[*,1], ytickname=replicate(' ',10), color=djs_icolor(defaultcolor)
    axis, /yaxis, yrange=yrange, ysty=3, charsize=charsize_6, charthick=postthick, $
      ythick=postthick, ytitle=textoidl(ytitle)
    oplot, !x.crange, [0,0], line=0, thick=postthick
    legend, '(b)', /left, /top, box=0, charsize=charsize_6, charthick=postthick
 
    im_symbols, 108, psize=psize2, fill=0, thick=postthick2, color=fsc_color(obscolor,75)
    oplot, ebv_obs, delta_oh_obs, ps=8
    im_symbols, 112, psize=psize3, fill=0, thick=postthick2, color=fsc_color(obscolor,75)
    oplot, ebv_obs_limit, delta_oh_obs_limit, ps=8

    im_symbols, 108, psize=psize2, fill=1, thick=postthick2, color=fsc_color(corcolor,74)
    oplot, ebv_cor, delta_oh_cor, ps=8
    im_symbols, 112, psize=psize3, fill=0, thick=postthick2, color=fsc_color(corcolor,74)
    oplot, ebv_cor_limit, delta_oh_cor_limit, ps=8

    im_symbols, 106, psize=psize2, fill=0, thick=postthick2, color=fsc_color(ewcolor,73)
    oplot, ebv_ew, delta_oh_ew, ps=8
    im_symbols, 112, psize=psize3, fill=0, thick=postthick2, color=fsc_color(ewcolor,73)
    oplot, ebv_ew_limit, delta_oh_ew_limit, ps=8

; inclination
    
    xrange = inclrange
    xtitle = 'Inclination [degrees]'

    plot, [0], [0], /nodata, /noerase, xthick=postthick, ythick=postthick, $
      charsize=charsize_6, charthick=postthick, thick=postthick, $
      xtitle=xtitle, ytitle=textoidl(ytitle), xstyle=3, ystyle=3, xrange=xrange, $
      yrange=yrange, position=pos[*,2], color=djs_icolor(defaultcolor)
    oplot, !x.crange, [0,0], line=0, thick=postthick
    legend, '(c)', /left, /top, box=0, charsize=charsize_6, charthick=postthick
 
    im_symbols, 108, psize=psize2, fill=0, thick=postthick2, color=fsc_color(obscolor,75)
    oplot, incl_obs, delta_oh_obs, ps=8
    im_symbols, 112, psize=psize3, fill=0, thick=postthick2, color=fsc_color(obscolor,75)
    oplot, incl_obs_limit, delta_oh_obs_limit, ps=8

    im_symbols, 108, psize=psize2, fill=1, thick=postthick2, color=fsc_color(corcolor,74)
    oplot, incl_cor, delta_oh_cor, ps=8
    im_symbols, 112, psize=psize3, fill=0, thick=postthick2, color=fsc_color(corcolor,74)
    oplot, incl_cor_limit, delta_oh_cor_limit, ps=8

    im_symbols, 106, psize=psize2, fill=0, thick=postthick2, color=fsc_color(ewcolor,73)
    oplot, incl_ew, delta_oh_ew, ps=8
    im_symbols, 112, psize=psize3, fill=0, thick=postthick2, color=fsc_color(ewcolor,73)
    oplot, incl_ew_limit, delta_oh_ew_limit, ps=8

; SB(H-alpha)
    
    xrange = hasbrange
    xtitle = textoidl('log [\Sigma(H\alpha)] [erg s^{-1} pc^{-2}]')

    plot, [0], [0], /nodata, /noerase, xthick=postthick, ythick=postthick, $
      charsize=charsize_6, charthick=postthick, thick=postthick, $
      xtitle=xtitle, ytitle='', xstyle=3, ystyle=11, xrange=xrange, $
      yrange=yrange, position=pos[*,3], ytickname=replicate(' ',10), color=djs_icolor(defaultcolor)
    axis, /yaxis, yrange=yrange, ysty=3, charsize=charsize_6, charthick=postthick, $
      ythick=postthick, ytitle=textoidl(ytitle)
    oplot, !x.crange, [0,0], line=0, thick=postthick
    legend, '(d)', /left, /top, box=0, charsize=charsize_6, charthick=postthick
 
    im_symbols, 108, psize=psize2, fill=0, thick=postthick2, color=fsc_color(obscolor,75)
    oplot, hasb_obs, delta_oh_obs, ps=8
    im_symbols, 112, psize=psize3, fill=0, thick=postthick2, color=fsc_color(obscolor,75)
    oplot, hasb_obs_limit, delta_oh_obs_limit, ps=8

    im_symbols, 108, psize=psize2, fill=1, thick=postthick2, color=fsc_color(corcolor,74)
    oplot, hasb_cor, delta_oh_cor, ps=8
    im_symbols, 112, psize=psize3, fill=0, thick=postthick2, color=fsc_color(corcolor,74)
    oplot, hasb_cor_limit, delta_oh_cor_limit, ps=8

    im_symbols, 106, psize=psize2, fill=0, thick=postthick2, color=fsc_color(ewcolor,73)
    oplot, hasb_ew, delta_oh_ew, ps=8
    im_symbols, 112, psize=psize3, fill=0, thick=postthick2, color=fsc_color(ewcolor,73)
    oplot, hasb_ew_limit, delta_oh_ew_limit, ps=8

    im_openclose, postscript=postscript, /close

; ---------------------------------------------------------------------------
; 12+log(O/H) [Characteristic] versus 12+log(O/H) [Integrated] - M91
; ---------------------------------------------------------------------------

    psname = 'int_12oh_vs_char_12oh_m91'
    im_openclose, pspath+psname, postscript=postscript, xsize=8.5, ysize=8.5, encapsulated=encapsulated

    pagemaker, nx=1, ny=1, xspace=0, yspace=0.0, width=7.0, height=7.0, $
      xmargin=[1.1,0.4], ymargin=[0.4,1.1], xpage=8.5, ypage=8.5, $
      position=pos, /normal

    xtitle = 'Characteristic 12 + log (O/H)'
    ytitle = 'Integrated 12 + log (O/H)'

    xrange = ohrange4
    yrange = xrange

; M91 - distinguish between observed, corrected, and EW abundances
    
    plot, [0], [0], /nodata, xthick=postthick, ythick=postthick, $
      charsize=charsize_9, charthick=postthick, thick=postthick, $
      xtitle=xtitle, ytitle=ytitle, xstyle=3, ystyle=3, xrange=xrange, $
      yrange=yrange, position=pos, color=djs_icolor(defaultcolor)
    oplot, !x.crange, !y.crange, line=0, thick=postthick
    legend, 'McGaugh (1991) O/H Calibration', /left, /top, box=0, charsize=charsize_7, charthick=postthick

    medcharerr = median(gradients.hii_m91_log12oh_char[1])
    medinterr = [gradients.int_obs_log12oh_m91[1],gradients.int_cor_log12oh_m91[1],gradients.int_ew_log12oh_m91[1]]
;   xohoff = 0.15 & yohoff = 0.22
;   oploterror, !x.crange[1]-xohoff, !y.crange[0]+yohoff, sqrt(medcharerr^2+0.1^2), $
;     sqrt(medinterr^2+0.1^2), ps=3, /data, /nohat, errstyle=0, thick=postthick, errthick=postthick
    xohoff = 0.1 & yohoff = 0.1
    oploterror, !x.crange[1]-xohoff, !y.crange[0]+yohoff, medcharerr, medinterr, $
      ps=3, /data, /nohat, errstyle=0, thick=postthick, errthick=postthick
    
    good = where((gradients.int_obs_log12oh_m91[0] gt -900.0) and (gradients.hii_m91_log12oh_char[0] gt -900.0) and $
      (gradients.int_log12oh_lower_limit eq 0),ngood)
    im_symbols, 108, psize=psize2, fill=0, thick=postthick2, color=fsc_color(obscolor,75)
    oplot, gradients[good].hii_m91_log12oh_char[0], gradients[good].int_obs_log12oh_m91[0], ps=8
    limit = where((gradients.int_obs_log12oh_m91[0] gt -900.0) and (gradients.hii_m91_log12oh_char[0] gt -900.0) and $
      (gradients.int_log12oh_lower_limit eq 1),nlimit)
    im_symbols, 112, psize=psize3, fill=0, thick=postthick2, color=fsc_color(obscolor,75)
    oplot, gradients[limit].hii_m91_log12oh_char[0], gradients[limit].int_obs_log12oh_m91[0], ps=8

    good = where((gradients.int_cor_log12oh_m91[0] gt -900.0) and (gradients.hii_m91_log12oh_char[0] gt -900.0) and $
      (gradients.int_log12oh_lower_limit eq 0),ngood)
    im_symbols, 108, psize=psize2, fill=1, thick=postthick2, color=fsc_color(corcolor,74)
    oplot, gradients[good].hii_m91_log12oh_char[0], gradients[good].int_cor_log12oh_m91[0], ps=8
    limit = where((gradients.int_cor_log12oh_m91[0] gt -900.0) and (gradients.hii_m91_log12oh_char[0] gt -900.0) and $
      (gradients.int_log12oh_lower_limit eq 1),nlimit)
    im_symbols, 112, psize=psize3, fill=0, thick=postthick2, color=fsc_color(corcolor,74)
    oplot, gradients[limit].hii_m91_log12oh_char[0], gradients[limit].int_cor_log12oh_m91[0], ps=8

    good = where((gradients.int_ew_log12oh_m91[0] gt -900.0) and (gradients.hii_m91_log12oh_char[0] gt -900.0) and $
      (gradients.int_log12oh_lower_limit eq 0),ngood)
    im_symbols, 106, psize=psize2, fill=0, thick=postthick2, color=fsc_color(ewcolor,73)
    oplot, gradients[good].hii_m91_log12oh_char[0], gradients[good].int_ew_log12oh_m91[0], ps=8
    limit = where((gradients.int_ew_log12oh_m91[0] gt -900.0) and (gradients.hii_m91_log12oh_char[0] gt -900.0) and $
      (gradients.int_log12oh_lower_limit eq 1),nlimit)
    im_symbols, 112, psize=psize3, fill=0, thick=postthick2, color=fsc_color(ewcolor,73)
    oplot, gradients[limit].hii_m91_log12oh_char[0], gradients[limit].int_ew_log12oh_m91[0], ps=8

; overplot the KKP99 data

    im_symbols, 105, psize=psize4, fill=1, thick=postthick2, color=fsc_color(kkpcolor,71)
    djs_oplot, kkp.oh12_char-0.17, kkp.oh12_global-0.17, ps=8

    im_openclose, postscript=postscript, /close

; ---------------------------------------------------------------------------
; 12+log(O/H) [Characteristic] versus 12+log(O/H) [Integrated] - PT05
; ---------------------------------------------------------------------------

    psname = 'int_12oh_vs_char_12oh_pt05'
    im_openclose, pspath+psname, postscript=postscript, xsize=8.5, ysize=8.5, encapsulated=encapsulated

    pagemaker, nx=1, ny=1, xspace=0, yspace=0.0, width=7.0, height=7.0, $
      xmargin=[1.1,0.4], ymargin=[0.4,1.1], xpage=8.5, ypage=8.5, $
      position=pos, /normal

    xtitle = 'Characteristic 12 + log (O/H)'
    ytitle = 'Integrated 12 + log (O/H)'

    xrange = ohrange4
    yrange = xrange

; PT05 - distinguish between observed, corrected, and EW abundances
    
    plot, [0], [0], /nodata, xthick=postthick, ythick=postthick, $
      charsize=charsize_9, charthick=postthick, thick=postthick, $
      xtitle=xtitle, ytitle=ytitle, xstyle=3, ystyle=3, xrange=xrange, $
      yrange=yrange, position=pos, color=djs_icolor(defaultcolor)
    oplot, !x.crange, !y.crange, line=0, thick=postthick
    legend, 'Pilyugin & Thuan (2005) O/H Calibration', /left, /top, box=0, charsize=charsize_7, charthick=postthick

;   medcharerr = median(gradients.hii_pt05_log12oh_char[1])
;   medinterr = [gradients.int_obs_log12oh_pt05[1],gradients.int_cor_log12oh_pt05[1],gradients.int_ew_log12oh_pt05[1]]
;   xohoff = 0.15 & yohoff = 0.22
;   oploterror, !x.crange[1]-xohoff, !y.crange[0]+yohoff, sqrt(medcharerr^2+0.1^2), $
;     sqrt(medinterr^2+0.1^2), ps=3, /data, /nohat, errstyle=0, thick=postthick, errthick=postthick
    xohoff = 0.1 & yohoff = 0.1
    oploterror, !x.crange[1]-xohoff, !y.crange[0]+yohoff, medcharerr, medinterr, $
      ps=3, /data, /nohat, errstyle=0, thick=postthick, errthick=postthick
    
    im_symbols, 108, psize=psize2, fill=0, thick=postthick2, color=fsc_color(obscolor,75)
    good = where((gradients.int_obs_log12oh_pt05[0] gt -900.0) and (gradients.hii_pt05_log12oh_char[0] gt -900.0) and $
      (gradients.int_log12oh_lower_limit eq 0),ngood)
    oplot, gradients[good].hii_pt05_log12oh_char[0], gradients[good].int_obs_log12oh_pt05[0], ps=8
    limit = where((gradients.int_obs_log12oh_pt05[0] gt -900.0) and (gradients.hii_pt05_log12oh_char[0] gt -900.0) and $
      (gradients.int_log12oh_lower_limit eq 1),nlimit)
    im_symbols, 112, psize=psize3, fill=0, thick=postthick2, color=fsc_color(obscolor,75)
    oplot, gradients[limit].hii_pt05_log12oh_char[0], gradients[limit].int_obs_log12oh_pt05[0], ps=8

    im_symbols, 108, psize=psize2, fill=1, thick=postthick2, color=fsc_color(corcolor,74)
    good = where((gradients.int_cor_log12oh_pt05[0] gt -900.0) and (gradients.hii_pt05_log12oh_char[0] gt -900.0) and $
      (gradients.int_log12oh_lower_limit eq 0),ngood)
    oplot, gradients[good].hii_pt05_log12oh_char[0], gradients[good].int_cor_log12oh_pt05[0], ps=8
    limit = where((gradients.int_cor_log12oh_pt05[0] gt -900.0) and (gradients.hii_pt05_log12oh_char[0] gt -900.0) and $
      (gradients.int_log12oh_lower_limit eq 1),nlimit)
    im_symbols, 112, psize=psize3, fill=0, thick=postthick2, color=fsc_color(corcolor,74)
    oplot, gradients[limit].hii_pt05_log12oh_char[0], gradients[limit].int_cor_log12oh_pt05[0], ps=8

    im_symbols, 106, psize=psize2, fill=0, thick=postthick2, color=fsc_color(ewcolor,73)
    good = where((gradients.int_ew_log12oh_pt05[0] gt -900.0) and (gradients.hii_pt05_log12oh_char[0] gt -900.0) and $
      (gradients.int_log12oh_lower_limit eq 0),ngood)
    oplot, gradients[good].hii_pt05_log12oh_char[0], gradients[good].int_ew_log12oh_pt05[0], ps=8
    limit = where((gradients.int_ew_log12oh_pt05[0] gt -900.0) and (gradients.hii_pt05_log12oh_char[0] gt -900.0) and $
      (gradients.int_log12oh_lower_limit eq 1),nlimit)
    im_symbols, 112, psize=psize3, fill=0, thick=postthick2, color=fsc_color(ewcolor,73)
    oplot, gradients[limit].hii_pt05_log12oh_char[0], gradients[limit].int_ew_log12oh_pt05[0], ps=8

    im_openclose, postscript=postscript, /close

; ---------------------------------------------------------------------------
; 12+log(O/H) [Characteristic] versus 12+log(O/H) [Integrated] - M91+PT05
; ---------------------------------------------------------------------------

    psname = 'int_12oh_vs_char_12oh'
    im_openclose, pspath+psname, postscript=postscript, xsize=8.5, ysize=8.5, encapsulated=encapsulated

    pagemaker, nx=1, ny=2, height=[3.5,3.5], width=7.0, xmargin=[1.2,0.3], $
      ymargin=[0.4,1.1], xspace=0.0, yspace=0.0, xpage=8.5, ypage=8.5, $
      position=pos, /normal

    xtitle = 'Characteristic 12 + log (O/H)'
    ytitle = 'Integrated 12 + log (O/H)'

    xrange = ohrange4
    yrange = xrange

; M91 - distinguish between observed, corrected, and EW abundances
    
    plot, [0], [0], /nodata, xthick=postthick, ythick=postthick, $
      charsize=charsize_7, charthick=postthick, thick=postthick, $
      xtitle='', ytitle='', xstyle=3, ystyle=3, xrange=xrange, $
      yrange=yrange, position=pos[*,0], xtickname=replicate(' ',10), color=djs_icolor(defaultcolor)
    oplot, !x.crange, !y.crange, line=0, thick=postthick
    legend, '(a) McGaugh (1991) O/H Calibration', /left, /top, box=0, charsize=charsize_5, charthick=postthick

    medcharerr = median(gradients.hii_m91_log12oh_char[1])
    medinterr = [gradients.int_obs_log12oh_m91[1],gradients.int_cor_log12oh_m91[1],gradients.int_ew_log12oh_m91[1]]
;   xohoff = 0.15 & yohoff = 0.22
;   oploterror, !x.crange[1]-xohoff, !y.crange[0]+yohoff, sqrt(medcharerr^2+0.1^2), $
;     sqrt(medinterr^2+0.1^2), ps=3, /data, /nohat, errstyle=0, thick=postthick, errthick=postthick
    xohoff = 0.1 & yohoff = 0.1
    oploterror, !x.crange[1]-xohoff, !y.crange[0]+yohoff, medcharerr, medinterr, $
      ps=3, /data, /nohat, errstyle=0, thick=postthick, errthick=postthick
    
    good = where((gradients.int_obs_log12oh_m91[0] gt -900.0) and (gradients.hii_m91_log12oh_char[0] gt -900.0) and $
      (gradients.int_log12oh_lower_limit eq 0),ngood)
    im_symbols, 108, psize=psize2, fill=0, thick=postthick2, color=fsc_color(obscolor,75)
    oplot, gradients[good].hii_m91_log12oh_char[0], gradients[good].int_obs_log12oh_m91[0], ps=8
    limit = where((gradients.int_obs_log12oh_m91[0] gt -900.0) and (gradients.hii_m91_log12oh_char[0] gt -900.0) and $
      (gradients.int_log12oh_lower_limit eq 1),nlimit)
    im_symbols, 112, psize=psize3, fill=0, thick=postthick2, color=fsc_color(obscolor,75)
    oplot, gradients[limit].hii_m91_log12oh_char[0], gradients[limit].int_obs_log12oh_m91[0], ps=8

    good = where((gradients.int_cor_log12oh_m91[0] gt -900.0) and (gradients.hii_m91_log12oh_char[0] gt -900.0) and $
      (gradients.int_log12oh_lower_limit eq 0),ngood)
    im_symbols, 108, psize=psize2, fill=1, thick=postthick2, color=fsc_color(corcolor,74)
    oplot, gradients[good].hii_m91_log12oh_char[0], gradients[good].int_cor_log12oh_m91[0], ps=8
    limit = where((gradients.int_cor_log12oh_m91[0] gt -900.0) and (gradients.hii_m91_log12oh_char[0] gt -900.0) and $
      (gradients.int_log12oh_lower_limit eq 1),nlimit)
    im_symbols, 112, psize=psize3, fill=0, thick=postthick2, color=fsc_color(corcolor,74)
    oplot, gradients[limit].hii_m91_log12oh_char[0], gradients[limit].int_cor_log12oh_m91[0], ps=8

    good = where((gradients.int_ew_log12oh_m91[0] gt -900.0) and (gradients.hii_m91_log12oh_char[0] gt -900.0) and $
      (gradients.int_log12oh_lower_limit eq 0),ngood)
    im_symbols, 106, psize=psize2, fill=0, thick=postthick2, color=fsc_color(ewcolor,73)
    oplot, gradients[good].hii_m91_log12oh_char[0], gradients[good].int_ew_log12oh_m91[0], ps=8
    limit = where((gradients.int_ew_log12oh_m91[0] gt -900.0) and (gradients.hii_m91_log12oh_char[0] gt -900.0) and $
      (gradients.int_log12oh_lower_limit eq 1),nlimit)
    im_symbols, 112, psize=psize3, fill=0, thick=postthick2, color=fsc_color(ewcolor,73)
    oplot, gradients[limit].hii_m91_log12oh_char[0], gradients[limit].int_ew_log12oh_m91[0], ps=8

; overplot the KKP99 data

    im_symbols, 105, psize=psize4, fill=1, thick=postthick2, color=fsc_color(kkpcolor,71)
    djs_oplot, kkp.oh12_char-0.17, kkp.oh12_global-0.17, ps=8

; y-title

    xyouts, pos[0,0]*0.4, (pos[2,0]-pos[1,1])/2.0+pos[1,1], ytitle, /normal, $
      orientation=90.0, align=0.5, charsize=charsize_7, charthick=postthick
    
; PT05 - distinguish between observed, corrected, and EW abundances
    
    plot, [0], [0], /nodata, /noerase, xthick=postthick, ythick=postthick, $
      charsize=charsize_7, charthick=postthick, thick=postthick, $
      xtitle=xtitle, ytitle='', xstyle=3, ystyle=3, xrange=xrange, $
      yrange=yrange, position=pos[*,1], color=djs_icolor(defaultcolor)
    oplot, !x.crange, !y.crange, line=0, thick=postthick
    legend, '(b) Pilyugin & Thuan (2005) O/H Calibration', /left, /top, box=0, charsize=charsize_5, charthick=postthick

;   medcharerr = median(gradients.hii_pt05_log12oh_char[1])
;   medinterr = [gradients.int_obs_log12oh_pt05[1],gradients.int_cor_log12oh_pt05[1],gradients.int_ew_log12oh_pt05[1]]
;   xohoff = 0.2 & yohoff = 0.25
;   oploterror, !x.crange[1]-xohoff, !y.crange[0]+yohoff, sqrt(medcharerr^2+0.1^2), $
;     sqrt(medinterr^2+0.1^2), ps=3, /data, /nohat, errstyle=0, thick=postthick, errthick=postthick
    xohoff = 0.1 & yohoff = 0.1
    oploterror, !x.crange[1]-xohoff, !y.crange[0]+yohoff, medcharerr, medinterr, $
      ps=3, /data, /nohat, errstyle=0, thick=postthick, errthick=postthick
    
    im_symbols, 108, psize=psize2, fill=0, thick=postthick2, color=fsc_color(obscolor,75)
    good = where((gradients.int_obs_log12oh_pt05[0] gt -900.0) and (gradients.hii_pt05_log12oh_char[0] gt -900.0) and $
      (gradients.int_log12oh_lower_limit eq 0),ngood)
    oplot, gradients[good].hii_pt05_log12oh_char[0], gradients[good].int_obs_log12oh_pt05[0], ps=8
    limit = where((gradients.int_obs_log12oh_pt05[0] gt -900.0) and (gradients.hii_pt05_log12oh_char[0] gt -900.0) and $
      (gradients.int_log12oh_lower_limit eq 1),nlimit)
    im_symbols, 112, psize=psize3, fill=0, thick=postthick2, color=fsc_color(obscolor,75)
    oplot, gradients[limit].hii_pt05_log12oh_char[0], gradients[limit].int_obs_log12oh_pt05[0], ps=8

    im_symbols, 108, psize=psize2, fill=1, thick=postthick2, color=fsc_color(corcolor,74)
    good = where((gradients.int_cor_log12oh_pt05[0] gt -900.0) and (gradients.hii_pt05_log12oh_char[0] gt -900.0) and $
      (gradients.int_log12oh_lower_limit eq 0),ngood)
    oplot, gradients[good].hii_pt05_log12oh_char[0], gradients[good].int_cor_log12oh_pt05[0], ps=8
    limit = where((gradients.int_cor_log12oh_pt05[0] gt -900.0) and (gradients.hii_pt05_log12oh_char[0] gt -900.0) and $
      (gradients.int_log12oh_lower_limit eq 1),nlimit)
    im_symbols, 112, psize=psize3, fill=0, thick=postthick2, color=fsc_color(corcolor,74)
    oplot, gradients[limit].hii_pt05_log12oh_char[0], gradients[limit].int_cor_log12oh_pt05[0], ps=8

    im_symbols, 106, psize=psize2, fill=0, thick=postthick2, color=fsc_color(ewcolor,73)
    good = where((gradients.int_ew_log12oh_pt05[0] gt -900.0) and (gradients.hii_pt05_log12oh_char[0] gt -900.0) and $
      (gradients.int_log12oh_lower_limit eq 0),ngood)
    oplot, gradients[good].hii_pt05_log12oh_char[0], gradients[good].int_ew_log12oh_pt05[0], ps=8
    limit = where((gradients.int_ew_log12oh_pt05[0] gt -900.0) and (gradients.hii_pt05_log12oh_char[0] gt -900.0) and $
      (gradients.int_log12oh_lower_limit eq 1),nlimit)
    im_symbols, 112, psize=psize3, fill=0, thick=postthick2, color=fsc_color(ewcolor,73)
    oplot, gradients[limit].hii_pt05_log12oh_char[0], gradients[limit].int_ew_log12oh_pt05[0], ps=8

    im_openclose, postscript=postscript, /close

; --------------------------------------------------    
; SELECT PLOTS FOR THE PAPER HERE
; --------------------------------------------------    

    if keyword_set(paper) then begin

       splog, 'Writing paper plots to '+paperpath+'.'
       paperplots = [$
         'int_12oh_vs_char_12oh_m91',$
         'int_12oh_vs_char_12oh_pt05',$
         '12oh_m91_residuals', $
         '12oh_pt05_residuals', $
         'slope_m91_vs_slope_pt05',$
         'all_siiha_vs_niiha',$
         'all_oiiihb_vs_oiihb'$
         ]+'.eps'

       for k = 0L, n_elements(paperplots)-1L do $
         spawn, ['/bin/cp -f '+html_path+htmlbase+'/'+paperplots[k]+' '+paperpath], /sh

    endif
    
; --------------------------------------------------    
; GENERATE THE WEB PAGE
; --------------------------------------------------    
    
    if keyword_set(postscript) and keyword_set(encapsulated) then begin
       im_ps2html, htmlbase, html_path=html_path, npscols=3, $
         cleanpng=0, _extra=extra
    endif

stop

return
end

;; ---------------------------------------------------------------------------
;; R23 vs O32 - Full Integrated Sample
;; ---------------------------------------------------------------------------
;
;    psname = 'all_oiiioiiihb_vs_oiiioii'
;    im_openclose, pspath+psname, postscript=postscript, xsize=8.5, ysize=8.0, encapsulated=encapsulated, cmyk=cmyk
;
;    pagemaker, nx=1, ny=1, xspace=0, yspace=0.0, width=6.6, height=6.6, $
;      xmargin=[1.5,0.4], ymargin=[0.3,1.1], xpage=8.5, ypage=8.0, $
;      position=pos, /normal
;
;; integrated spectra, HII galaxies    
;
;    sf = where(strtrim(atlasdust.bpt_class,2) eq 'HII')
;    sfatlasdust = atlasdust[sf]
;    sfatlasnodust = atlasnodust[sf]
;    
;    cutindx = cmset_op(strtrim(sfatlasdust.galaxy,2),'and',/not2,strtrim(gradients.galaxy,2),/index)
;    lineratio, sfatlasnodust[cutindx], ['OII_3727','OIII_4959','OIII_5007'], 'H_BETA', $
;      ['OIII_4959','OIII_5007'], 'OII_3727', x, xerr, y, yerr, snrcut=3.0, index=indx, nindex=nindx
;
;    cutindxint = cmset_op(strtrim(sfatlasdust.galaxy,2),'and',strtrim(gradients.galaxy,2),/index)
;    lineratio, sfatlasnodust[cutindxint], ['OII_3727','OIII_4959','OIII_5007'], 'H_BETA', $
;      ['OIII_4959','OIII_5007'], 'OII_3727', xint, xerrint, yint, yerrint, snrcut=3.0, index=indxint, nindex=nindxint
;    
;    xtitle = 'log (R_{23})_{cor}'
;    ytitle = 'log (O_{32})_{cor}'
;;   xtitle = 'log {[O II] \lambda3727+[O III] \lambda5007)/H\beta}_{cor}'
;;   ytitle = 'log ([O III] \lambda5007/[O II] \lambda3727)_{cor}'
;
;    xrange = r23range
;    yrange = o32range
;
;    good = where((samplehii.zstrong_r23 gt -900.0) and (samplehii.zstrong_o32 gt -900.0))
;    xregion = samplehii[good].zstrong_r23 & xerrregion = samplehii[good].zstrong_r23_err
;    yregion = samplehii[good].zstrong_o32 & yerrregion = samplehii[good].zstrong_o32_err
;
;;   good = where((samplehii.oii_h_beta gt -900.0) and (samplehii.oiii_5007_h_beta gt -900.0))
;;   xregion = samplehii[good].oii_h_beta + samplehii[good].oiii_5007_h_beta
;;   xerrregion = sqrt(samplehii[good].oii_h_beta_err^2 + samplehii[good].oiii_5007_h_beta_err^2)
;;   yregion = samplehii[good].oiii_5007_h_beta - samplehii[good].oii_h_beta
;;   yerrregion = sqrt(samplehii[good].oii_h_beta_err^2 + samplehii[good].oiii_5007_h_beta_err^2)
;
;    atlas1d_lineplot, x, y, xerr, yerr, plottype=1, postscript=postscript, $
;      xtitle=xtitle, ytitle=ytitle, xrange=xrange, yrange=yrange, legendtype=0, $
;      /right, /top, position=pos[*,0], atlaspsize=psizeall, atlassym=104, atlascolor=atlascolor, $
;      atlasfill=1, symthick=postthick, postthick=postthick, charsize=charsize_8, blackwhite=blackwhite
;    legend, '(b)', /left, /top, box=0, charsize=2.0, charthick=postthick
;
;; overplot the HII regions
;
;    im_symbols, 108, fill=1L, psize=hiipsize, color=fsc_color(hiicolor,25), thick=postthick2
;    djs_oplot, xregion, yregion, psym=8
;
;; overplot the small sample    
;    
;    im_symbols, 106, fill=1L, psize=1.7, color=fsc_color(intcolor,73), thick=postthick2
;    djs_oplot, xint, yint, psym=8
;
;; overplot DIG points
;
;    digr23 = alog10(1.33*10^g99_oiioiii.oiiihb*(1.0 + 10^g99_oiioiii.oiioiii))
;    digo32 = alog10(1.33) - g99_oiioiii.oiioiii
;    
;    djs_oplot, digr23, digo32, ps=7, syms=digpsize, thick=postthick, color=fsc_color(digcolor,72)
;    
;;   gg99 = where((g99.oiihb gt -900) and (g99.oiiihb gt -900),ngg99)
;;   djs_oplot, alog10(g99[gg99].oiiihb), alog10(g99[gg99].oiihb), ps=7, syms=digpsize, $
;;     thick=postthick, color=fsc_color('dark orchid',75)
;
;;   plot_kewley_grids, model=8L, plotnumber=29L, labeltype=1L, /siioffset
;
;    im_openclose, postscript=postscript, /close
;
;; ---------------------------------------------------------------------------
;; [O II]/[O III] vs [N II]/Ha - Full Integrated Sample
;; ---------------------------------------------------------------------------
;
;    psname = 'all_oiioiii_vs_niiha'
;    im_openclose, pspath+psname, postscript=postscript, xsize=8.5, ysize=8.0, encapsulated=encapsulated
;
;    pagemaker, nx=1, ny=1, xspace=0, yspace=0.0, width=6.6, height=6.6, $
;      xmargin=[1.5,0.4], ymargin=[0.3,1.1], xpage=8.5, ypage=8.0, $
;      position=pos, /normal
;
;; integrated spectra, HII galaxies    
;    
;    indx = where(strtrim(atlasdust.bpt_class,2) eq 'HII',nindx)
;    lineratio, atlasnodust[indx], 'OII_3727', 'OIII_5007', 'NII_6584', 'H_ALPHA', $
;      x, xerr, y, yerr, snrcut=3.0, index=cutindx, nindex=ncutindx
;
;    indx2 = lindgen(ncutindx)
;;   indx1 = cmset_op(atlasdust[indx[cutindx]].atlas_id,'and',intdust.atlas_id,/index)
;;   indx2 = cmset_op(atlasdust[indx[cutindx]].atlas_id,'and',/not2,intdust.atlas_id,/index)
;
;    xtitle = 'log ([O II] \lambda3727/[O III] \lambda5007)'
;    ytitle = 'log ([N II] \lambda6584/H\alpha)'
;
;    xrange = oiioiiirange
;    yrange = niiharange3
;
;    good = where((samplehii.oii_oiii_5007 gt -900.0) and (samplehii.nii_6584_h_alpha gt -900.0))
;    xregion = samplehii[good].oii_oiii_5007 & xerrregion = samplehii[good].oii_oiii_5007_err
;    yregion = samplehii[good].nii_6584_h_alpha & yerrregion = samplehii[good].nii_6584_h_alpha_err
;
;    atlas1d_lineplot, x[indx2], y[indx2], xerr[indx2], yerr[indx2], plottype=1, postscript=postscript, $
;      xtitle=xtitle, ytitle=ytitle, xrange=xrange, yrange=yrange, legendtype=0, $
;      /right, /top, position=pos[*,0], atlaspsize=psizeall, atlascolor=atlascolor, $
;      atlasfill=1, symthick=postthick2, postthick=postthick;, hiipsize=hiipsize, hiicolor=hiicolor, $
;;     xregion=xregion, yregion=yregion, xerrregion=xerrregion, yerrregion=yerrregion
;
;;   im_symbols, 106, fill=1L, psize=psizeall, color=fsc_color('dark green',50), thick=postthick2
;;   djs_oplot, x[indx1], y[indx1], psym=8
;    
;; overplot DIG points
; 
;;   gg99 = where((g99.oiihb gt -900) and (g99.oiiihb gt -900) and (g99.niiha gt -900),ngg99)
;;   djs_oplot, alog10(g99[gg99].oiihb/g99[gg99].oiiihb), alog10(g99[gg99].niiha), ps=7, syms=digpsize, $
;;     thick=postthick, color=fsc_color('dark orchid',75)
;
;; overplot the HII regions
;
;    im_symbols, 108, fill=1L, psize=hiipsize, color=fsc_color(hiicolor,25), thick=postthick2
;    djs_oplot, xregion, yregion, psym=8
; 
;;   plot_kewley_grids, model=8L, plotnumber=30L, /overplot, labeltype=1L, /siioffset
;    
;; legend
;
;    legend, '(a)', /left, /top, box=0, charsize=2.0, charthick=postthick
;    
;    im_openclose, postscript=postscript, /close
;
;; ---------------------------------------------------------------------------
;; [N II]/[S II] vs [N II]/Ha - Full Integrated Sample
;; ---------------------------------------------------------------------------
;
;    psname = 'all_niisii_vs_niiha'
;    im_openclose, pspath+psname, postscript=postscript, xsize=8.5, ysize=8.0, encapsulated=encapsulated
;
;    pagemaker, nx=1, ny=1, xspace=0, yspace=0.0, width=6.6, height=6.6, $
;      xmargin=[1.5,0.4], ymargin=[0.3,1.1], xpage=8.5, ypage=8.0, $
;      position=pos, /normal
;
;; integrated spectra, HII galaxies    
;    
;    indx = where(strtrim(atlasdust.bpt_class,2) eq 'HII')
;    lineratio, atlasnodust[indx], 'NII_6584', 'SII_6716', 'NII_6584', 'H_ALPHA', $
;      x, xerr, y, yerr, snrcut=1.0
;    
;    xtitle = 'log ([N II] \lambda6584/[S II] \lambda6716)'
;    ytitle = 'log ([N II] \lambda6584/H\alpha)'
;
;    xrange = niisiirange
;    yrange = niiharange3
;
;    good = where((samplehii.sii_6716_h_alpha gt -900.0) and (samplehii.nii_6584_h_alpha gt -900.0)) 
;    xregion = samplehii[good].nii_6584_h_alpha - samplehii[good].sii_6716_h_alpha
;    xerrregion = sqrt(samplehii[good].nii_6584_h_alpha_err^2 + samplehii[good].sii_6716_h_alpha_err^2)
;    yregion = samplehii[good].nii_6584_h_alpha
;    yerrregion = samplehii[good].nii_6584_h_alpha_err
;
;    atlas1d_lineplot, x, y, xerr, yerr, plottype=1, postscript=postscript, $
;      xtitle=xtitle, ytitle=ytitle, xrange=xrange, yrange=yrange, legendtype=0, $
;      /right, /top, position=pos[*,0], atlaspsize=psizeall, atlascolor=atlascolor, atlasfill=1, $
;      postthick=postthick     ;, hiipsize=hiipsize, hiicolor=hiicolor, $
;;     xregion=xregion, yregion=yregion, xerrregion=xerrregion, yerrregion=yerrregion
;
;; overplot the HII regions
;
;    im_symbols, 108, fill=1L, psize=hiipsize, color=fsc_color(hiicolor,25), thick=postthick2
;    djs_oplot, xregion, yregion, psym=8
; 
;; overplot DIG points
; 
;    djs_oplot, alog10(h03.niiha/h03.siiha), alog10(h03.niiha), ps=7, syms=digpsize, $
;      thick=postthick, color=fsc_color('dark orchid',75)
;
;;   gg99 = where((g99.siiha gt -900) and (g99.niiha gt -900),ngg99)
;;   djs_oplot, alog10(g99[gg99].siiha), alog10(g99[gg99].niiha), ps=7, syms=digpsize, $
;;     thick=postthick, color=fsc_color('dark orchid',75)
;
;; legend
;
;    legend, '(a)', /left, /top, box=0, charsize=2.0, charthick=postthick
;    
;    im_openclose, postscript=postscript, /close

