pro zintegrated_talkplots, gradients, samplehii, atlasdust, atlasnodust, postscript=postscript
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
    pspath = atlas_path(/web)+'analysis/'+htmlbase+'/talk/'

    rr25_int = 0.4
    
; setup some plotting variables

    @'xyrange_zintegrated'

    if keyword_set(postscript) then begin
       postthick = 8.0
       postthick2 = 8.0
       defaultcolor = 'black'
       talkcolor = 'white'
    endif else begin
       postthick = 2.0
       postthick2 = postthick
       im_window, 0, xratio=0.5, /square
       defaultcolor = 'white'
       talkcolor = ''
    endelse

    psize2 = 1.8
    psize3 = 2.9
    psize4 = 2.0
    hiipsize = 0.5
    digpsize = 1.5
    hiicolor = 'dark grey' ; 'blue'
;   atlascolor = 'grey'
    
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
; 12+log(O/H) Residuals [Int. minus Char.] vs various quantities - M91
; ---------------------------------------------------------------------------

    psname = '12oh_m91_residuals'
    im_openclose, pspath+psname, postscript=postscript, xsize=8.5, ysize=6.0, encapsulated=encapsulated

    pagemaker, nx=2, ny=2, yspace=1.0, xspace=0.0, width=3.2*[1,1], height=1.75*[1,1], $
      xmargin=[1.1,0.5], ymargin=[0.4,1.1], xpage=8.5, ypage=6.0, position=pos, /normal

    polyfill, [1,1,0,0,1], [1,0,0,1,1], /normal, color=djs_icolor('black')

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
      charsize=charsize_5, charthick=postthick, thick=postthick, $
      xtitle=xtitle, ytitle=ytitle, xstyle=3, ystyle=3, xrange=xrange, $
      yrange=yrange, position=pos[*,0], xtickinterval=0.5, $
      color=djs_icolor(talkcolor), /noerase
    oplot, !x.crange, [0,0], line=0, thick=postthick, color=djs_icolor(talkcolor)
    legend, '(a)', /left, /top, box=0, charsize=charsize_5, charthick=postthick, $
      color=djs_icolor(talkcolor), textcolor=djs_icolor(talkcolor)
 
; overplot the KKP99 data

    im_symbols, 105, psize=psize4, fill=1, thick=postthick2, color=djs_icolor('magenta blue')
    djs_oplot, kkp.gradient, kkp.oh12_global-kkp.oh12_char, ps=8

    im_symbols, 108, psize=psize2, fill=0, thick=postthick2, color=djs_icolor('red')
    oplot, slope_obs, delta_oh_obs, ps=8
    im_symbols, 112, psize=psize3, fill=0, thick=postthick2, color=djs_icolor('red')
    oplot, slope_obs_limit, delta_oh_obs_limit, ps=8

    im_symbols, 108, psize=psize2, fill=1, color=djs_icolor('blue')
    oplot, slope_cor, delta_oh_cor, ps=8
    im_symbols, 112, psize=psize3, fill=0, thick=postthick2, color=djs_icolor('blue')
    oplot, slope_cor_limit, delta_oh_cor_limit, ps=8

    im_symbols, 106, psize=psize2, fill=0, thick=postthick2, color=djs_icolor('green')
    oplot, slope_ew, delta_oh_ew, ps=8
    im_symbols, 112, psize=psize3, fill=0, thick=postthick2, color=djs_icolor('green')
    oplot, slope_ew_limit, delta_oh_ew_limit, ps=8

; ebv
    
    xrange = ebvrange
    xtitle = 'E(B-V) [mag]'

    plot, [0], [0], /nodata, /noerase, xthick=postthick, ythick=postthick, $
      charsize=charsize_5, charthick=postthick, thick=postthick, $
      xtitle=xtitle, ytitle='', xstyle=3, ystyle=11, xrange=xrange, $
      yrange=yrange, position=pos[*,1], ytickname=replicate(' ',10), $
      color=djs_icolor(talkcolor)
    axis, /yaxis, yrange=yrange, ysty=3, charsize=charsize_5, charthick=postthick, $
      ythick=postthick, ytitle=textoidl(ytitle), color=djs_icolor(talkcolor)
    oplot, !x.crange, [0,0], line=0, thick=postthick, color=djs_icolor(talkcolor)
    legend, '(b)', /left, /top, box=0, charsize=charsize_5, charthick=postthick, $
      color=djs_icolor(talkcolor), textcolor=djs_icolor(talkcolor)
 
    im_symbols, 108, psize=psize2, fill=0, thick=postthick2, color=djs_icolor('red')
    oplot, ebv_obs, delta_oh_obs, ps=8
    im_symbols, 112, psize=psize3, fill=0, thick=postthick2, color=djs_icolor('red')
    oplot, ebv_obs_limit, delta_oh_obs_limit, ps=8

    im_symbols, 108, psize=psize2, fill=1, color=djs_icolor('blue')
    oplot, ebv_cor, delta_oh_cor, ps=8
    im_symbols, 112, psize=psize3, fill=0, thick=postthick2, color=djs_icolor('blue')
    oplot, ebv_cor_limit, delta_oh_cor_limit, ps=8

    im_symbols, 106, psize=psize2, fill=0, thick=postthick2, color=djs_icolor('green')
    oplot, ebv_ew, delta_oh_ew, ps=8
    im_symbols, 112, psize=psize3, fill=0, thick=postthick2, color=djs_icolor('green')
    oplot, ebv_ew_limit, delta_oh_ew_limit, ps=8

; inclination
    
    xrange = inclrange
    xtitle = 'Inclination [degrees]'

    plot, [0], [0], /nodata, /noerase, xthick=postthick, ythick=postthick, $
      charsize=charsize_5, charthick=postthick, thick=postthick, $
      xtitle=xtitle, ytitle=ytitle, xstyle=3, ystyle=3, xrange=xrange, $
      yrange=yrange, position=pos[*,2], color=djs_icolor(talkcolor)
    oplot, !x.crange, [0,0], line=0, thick=postthick, color=djs_icolor(talkcolor)
    legend, '(c)', /left, /top, box=0, charsize=charsize_5, charthick=postthick, $
      color=djs_icolor(talkcolor), textcolor=djs_icolor(talkcolor)
 
    im_symbols, 108, psize=psize2, fill=0, thick=postthick2, color=djs_icolor('red')
    oplot, incl_obs, delta_oh_obs, ps=8
    im_symbols, 112, psize=psize3, fill=0, thick=postthick2, color=djs_icolor('red')
    oplot, incl_obs_limit, delta_oh_obs_limit, ps=8

    im_symbols, 108, psize=psize2, fill=1, color=djs_icolor('blue')
    oplot, incl_cor, delta_oh_cor, ps=8
    im_symbols, 112, psize=psize3, fill=0, thick=postthick2, color=djs_icolor('blue')
    oplot, incl_cor_limit, delta_oh_cor_limit, ps=8

    im_symbols, 106, psize=psize2, fill=0, thick=postthick2, color=djs_icolor('green')
    oplot, incl_ew, delta_oh_ew, ps=8
    im_symbols, 112, psize=psize3, fill=0, thick=postthick2, color=djs_icolor('green')
    oplot, incl_ew_limit, delta_oh_ew_limit, ps=8

; SB(H-alpha)
    
    xrange = hasbrange
    xtitle = textoidl('log [\Sigma(H\alpha)] [erg s^{-1} pc^{-2}]')

    plot, [0], [0], /nodata, /noerase, xthick=postthick, ythick=postthick, $
      charsize=charsize_5, charthick=postthick, thick=postthick, $
      xtitle=xtitle, ytitle='', xstyle=3, ystyle=11, xrange=xrange, $
      yrange=yrange, position=pos[*,3], ytickname=replicate(' ',10), $
      color=djs_icolor(talkcolor)
    axis, /yaxis, yrange=yrange, ysty=3, charsize=charsize_5, charthick=postthick, $
      ythick=postthick, ytitle=textoidl(ytitle), color=djs_icolor(talkcolor)
    oplot, !x.crange, [0,0], line=0, thick=postthick, color=djs_icolor(talkcolor)
    legend, '(d)', /left, /top, box=0, charsize=charsize_5, charthick=postthick, $
      color=djs_icolor(talkcolor), textcolor=djs_icolor(talkcolor)
 
    im_symbols, 108, psize=psize2, fill=0, thick=postthick2, color=djs_icolor('red')
    oplot, hasb_obs, delta_oh_obs, ps=8
    im_symbols, 112, psize=psize3, fill=0, thick=postthick2, color=djs_icolor('red')
    oplot, hasb_obs_limit, delta_oh_obs_limit, ps=8

    im_symbols, 108, psize=psize2, fill=1, color=djs_icolor('blue')
    oplot, hasb_cor, delta_oh_cor, ps=8
    im_symbols, 112, psize=psize3, fill=0, thick=postthick2, color=djs_icolor('blue')
    oplot, hasb_cor_limit, delta_oh_cor_limit, ps=8

    im_symbols, 106, psize=psize2, fill=0, thick=postthick2, color=djs_icolor('green')
    oplot, hasb_ew, delta_oh_ew, ps=8
    im_symbols, 112, psize=psize3, fill=0, thick=postthick2, color=djs_icolor('green')
    oplot, hasb_ew_limit, delta_oh_ew_limit, ps=8

    im_openclose, postscript=postscript, /close
    if keyword_set(postscript) then spawn, ['convert '+pspath+psname+'.ps '+pspath+psname+'.png'], /sh

; ---------------------------------------------------------------------------
; 12+log(O/H) [Characteristic] versus 12+log(O/H) [Integrated] - M91+PT05
; ---------------------------------------------------------------------------

    psname = 'int_12oh_vs_char_12oh'
    im_openclose, pspath+psname, postscript=postscript, xsize=8.5, ysize=7.5, encapsulated=encapsulated

    pagemaker, nx=1, ny=2, height=[3.0,3.0], width=7.0, xmargin=[1.2,0.3], $
      ymargin=[0.4,1.1], xspace=0.0, yspace=0.0, xpage=8.5, ypage=7.5, $
      position=pos, /normal

    polyfill, [1,1,0,0,1], [1,0,0,1,1], /normal, color=djs_icolor('black')

    xtitle = 'Characteristic 12 + log (O/H)'
    ytitle = 'Integrated 12 + log (O/H)'

    xrange = ohrange4
    yrange = xrange

; M91 - distinguish between observed, corrected, and EW abundances
    
    plot, [0], [0], /nodata, xthick=postthick, ythick=postthick, $
      charsize=charsize_7, charthick=postthick, thick=postthick, $
      xtitle='', ytitle='', xstyle=3, ystyle=3, xrange=xrange, $
      yrange=yrange, position=pos[*,0], xtickname=replicate(' ',10), $
      color=djs_icolor(talkcolor), /noerase
    oplot, !x.crange, !y.crange, line=0, thick=postthick, color=djs_icolor(talkcolor)
    legend, '(a) McGaugh (1991) O/H Calibration', /left, /top, box=0, charsize=charsize_5, $
      charthick=postthick, color=djs_icolor(talkcolor), textcolor=djs_icolor(talkcolor)

; overplot the KKP99 data

    im_symbols, 105, psize=psize4, fill=1, thick=postthick2, color=djs_icolor('magenta blue')
    djs_oplot, kkp.oh12_char-0.17, kkp.oh12_global-0.17, ps=8, color='magenta blue'

; error bar    
    
    medcharerr = median(gradients.hii_m91_log12oh_char[1])
    medinterr = [gradients.int_obs_log12oh_m91[1],gradients.int_cor_log12oh_m91[1],gradients.int_ew_log12oh_m91[1]]
    xohoff = 0.15 & yohoff = 0.22
    oploterror, !x.crange[1]-xohoff, !y.crange[0]+yohoff, sqrt(medcharerr^2+0.1^2), $
      sqrt(medinterr^2+0.1^2), ps=3, /data, /nohat, errstyle=0, thick=postthick, $
      errthick=postthick, color=djs_icolor(talkcolor), errcolor=djs_icolor(talkcolor)
;   oploterror, !x.crange[1]-xohoff, !y.crange[0]+yohoff, medcharerr, medinterr, $
;     ps=3, /data, /nohat, errstyle=0, thick=postthick, errthick=postthick
    
    good = where((gradients.int_obs_log12oh_m91[0] gt -900.0) and (gradients.hii_m91_log12oh_char[0] gt -900.0) and $
      (gradients.int_log12oh_lower_limit eq 0),ngood)
    im_symbols, 108, psize=psize2, fill=0, thick=postthick2, color=djs_icolor('red')
    oplot, gradients[good].hii_m91_log12oh_char[0], gradients[good].int_obs_log12oh_m91[0], ps=8
    limit = where((gradients.int_obs_log12oh_m91[0] gt -900.0) and (gradients.hii_m91_log12oh_char[0] gt -900.0) and $
      (gradients.int_log12oh_lower_limit eq 1),nlimit)
    im_symbols, 112, psize=psize3, fill=0, thick=postthick2, color=djs_icolor('red')
    oplot, gradients[limit].hii_m91_log12oh_char[0], gradients[limit].int_obs_log12oh_m91[0], ps=8

    good = where((gradients.int_cor_log12oh_m91[0] gt -900.0) and (gradients.hii_m91_log12oh_char[0] gt -900.0) and $
      (gradients.int_log12oh_lower_limit eq 0),ngood)
    im_symbols, 108, psize=psize2, fill=1, color=djs_icolor('blue')
    oplot, gradients[good].hii_m91_log12oh_char[0], gradients[good].int_cor_log12oh_m91[0], ps=8
    limit = where((gradients.int_cor_log12oh_m91[0] gt -900.0) and (gradients.hii_m91_log12oh_char[0] gt -900.0) and $
      (gradients.int_log12oh_lower_limit eq 1),nlimit)
    im_symbols, 112, psize=psize3, fill=0, thick=postthick2, color=djs_icolor('blue')
    oplot, gradients[limit].hii_m91_log12oh_char[0], gradients[limit].int_cor_log12oh_m91[0], ps=8

    good = where((gradients.int_ew_log12oh_m91[0] gt -900.0) and (gradients.hii_m91_log12oh_char[0] gt -900.0) and $
      (gradients.int_log12oh_lower_limit eq 0),ngood)
    im_symbols, 106, psize=psize2, fill=0, thick=postthick2, color=djs_icolor('green')
    oplot, gradients[good].hii_m91_log12oh_char[0], gradients[good].int_ew_log12oh_m91[0], ps=8
    limit = where((gradients.int_ew_log12oh_m91[0] gt -900.0) and (gradients.hii_m91_log12oh_char[0] gt -900.0) and $
      (gradients.int_log12oh_lower_limit eq 1),nlimit)
    im_symbols, 112, psize=psize3, fill=0, thick=postthick2, color=djs_icolor('green')
    oplot, gradients[limit].hii_m91_log12oh_char[0], gradients[limit].int_ew_log12oh_m91[0], ps=8

; y-title

    xyouts, pos[0,0]*0.4, (pos[2,0]-pos[1,1])/2.0+pos[1,1], ytitle, /normal, $
      orientation=90.0, align=0.5, charsize=charsize_7, charthick=postthick, $
      color=djs_icolor(talkcolor)
    
; PT05 - distinguish between observed, corrected, and EW abundances
    
    plot, [0], [0], /nodata, /noerase, xthick=postthick, ythick=postthick, $
      charsize=charsize_7, charthick=postthick, thick=postthick, $
      xtitle=xtitle, ytitle='', xstyle=3, ystyle=3, xrange=xrange, $
      yrange=yrange, position=pos[*,1], color=djs_icolor(talkcolor)
    oplot, !x.crange, !y.crange, line=0, thick=postthick, color=djs_icolor(talkcolor)
    legend, '(b) Pilyugin & Thuan (2005) O/H Calibration', /left, /top, box=0, charsize=charsize_5, $
      charthick=postthick, color=djs_icolor(talkcolor), textcolor=djs_icolor(talkcolor)

    medcharerr = median(gradients.hii_pt05_log12oh_char[1])
    medinterr = [gradients.int_obs_log12oh_pt05[1],gradients.int_cor_log12oh_pt05[1],gradients.int_ew_log12oh_pt05[1]]
    xohoff = 0.15 & yohoff = 0.2
    oploterror, !x.crange[1]-xohoff, !y.crange[0]+yohoff, sqrt(medcharerr^2+0.1^2), $
      sqrt(medinterr^2+0.1^2), ps=3, /data, /nohat, errstyle=0, thick=postthick, $
      errthick=postthick, color=djs_icolor(talkcolor), errcolor=djs_icolor(talkcolor)
    
    im_symbols, 108, psize=psize2, fill=0, thick=postthick2, color=djs_icolor('red')
    good = where((gradients.int_obs_log12oh_pt05[0] gt -900.0) and (gradients.hii_pt05_log12oh_char[0] gt -900.0) and $
      (gradients.int_log12oh_lower_limit eq 0),ngood)
    oplot, gradients[good].hii_pt05_log12oh_char[0], gradients[good].int_obs_log12oh_pt05[0], ps=8
    limit = where((gradients.int_obs_log12oh_pt05[0] gt -900.0) and (gradients.hii_pt05_log12oh_char[0] gt -900.0) and $
      (gradients.int_log12oh_lower_limit eq 1),nlimit)
    im_symbols, 112, psize=psize3, fill=0, thick=postthick2, color=djs_icolor('red')
    oplot, gradients[limit].hii_pt05_log12oh_char[0], gradients[limit].int_obs_log12oh_pt05[0], ps=8

    im_symbols, 108, psize=psize2, fill=1, color=djs_icolor('blue')
    good = where((gradients.int_cor_log12oh_pt05[0] gt -900.0) and (gradients.hii_pt05_log12oh_char[0] gt -900.0) and $
      (gradients.int_log12oh_lower_limit eq 0),ngood)
    oplot, gradients[good].hii_pt05_log12oh_char[0], gradients[good].int_cor_log12oh_pt05[0], ps=8
    limit = where((gradients.int_cor_log12oh_pt05[0] gt -900.0) and (gradients.hii_pt05_log12oh_char[0] gt -900.0) and $
      (gradients.int_log12oh_lower_limit eq 1),nlimit)
    im_symbols, 112, psize=psize3, fill=0, thick=postthick2, color=djs_icolor('blue')
    oplot, gradients[limit].hii_pt05_log12oh_char[0], gradients[limit].int_cor_log12oh_pt05[0], ps=8

    im_symbols, 106, psize=psize2, fill=0, thick=postthick2, color=djs_icolor('green')
    good = where((gradients.int_ew_log12oh_pt05[0] gt -900.0) and (gradients.hii_pt05_log12oh_char[0] gt -900.0) and $
      (gradients.int_log12oh_lower_limit eq 0),ngood)
    oplot, gradients[good].hii_pt05_log12oh_char[0], gradients[good].int_ew_log12oh_pt05[0], ps=8
    limit = where((gradients.int_ew_log12oh_pt05[0] gt -900.0) and (gradients.hii_pt05_log12oh_char[0] gt -900.0) and $
      (gradients.int_log12oh_lower_limit eq 1),nlimit)
    im_symbols, 112, psize=psize3, fill=0, thick=postthick2, color=djs_icolor('green')
    oplot, gradients[limit].hii_pt05_log12oh_char[0], gradients[limit].int_ew_log12oh_pt05[0], ps=8

    im_openclose, postscript=postscript, /close
    if keyword_set(postscript) then spawn, ['convert '+pspath+psname+'.ps '+pspath+psname+'.png'], /sh

return
end
