; ------------------------------------------------------------
; MB vs 12+log(O/H)
; ------------------------------------------------------------

    psname = 'MB_vs_12oh_kk04_r23'
    im_openclose, pspath+psname, postscript=postscript, xsize=8.5, ysize=8.0

    pagemaker, nx=1, ny=1, xspace=0, yspace=0.0, width=6.6, height=6.6, $
      xmargin=[1.5,0.4], ymargin=[0.3,1.1], xpage=8.5, ypage=8.0, $
      position=pos, /normal

    mbcut = -18.0
    
    indx = where((sdssnodust.zstrong_12oh_kk04_r23_upper gt -900) and $
      (sdssancillary.M_B gt -900),nindx)

    x = sdssancillary[indx].M_B
    xerr = sdssancillary[indx].M_B_err

    y = sdssnodust[indx].zstrong_12oh_kk04_r23_upper
    yerr = sdssnodust[indx].zstrong_12oh_kk04_r23_upper_err

; Atlas

    indxatlas = where((atlasnodust.zstrong_12oh_kk04_r23_upper gt -900) and $
      (atlasnodust.M_B gt -900),nindxatlas)
;     (atlasnodust.M_B gt -900) and (atlasnodust.m_b lt mbcut),nindxatlas)

    atlast = atlasnodust[indxatlas].lit_t
    atlastype = atlasnodust[indxatlas].lit_type

    xatlas = atlasnodust[indxatlas].M_B
    xerratlas = atlasnodust[indxatlas].M_B_err

    yatlas = xatlas*0.0
    yerratlas = xerratlas*0.0
    
    lower = where((strmatch(atlastype,'*Irr*',/fold) eq 1B) and (xatlas gt mbcut),nlower,comp=upper,ncomp=nupper)

    if (nlower ne 0L) then begin
       yatlas[lower] = atlasnodust[indxatlas[lower]].zstrong_12oh_kk04_r23_lower
       yerratlas[lower] = atlasnodust[indxatlas[lower]].zstrong_12oh_kk04_r23_lower_err
    endif
    if (nupper ne 0L) then begin
       yatlas[upper] = atlasnodust[indxatlas[upper]].zstrong_12oh_kk04_r23_upper
       yerratlas[upper] = atlasnodust[indxatlas[upper]].zstrong_12oh_kk04_r23_upper_err
    endif

;   yatlas = atlasnodust[indxatlas].zstrong_12oh_kk04_r23_upper
;   yerratlas = atlasnodust[indxatlas].zstrong_12oh_kk04_r23_upper_err

    ohindx = where(atlasnodust[indxatlas].lit_log12oh gt -900.0,nohindx,comp=strong)
    yatlas[ohindx] = atlasnodust[indxatlas[ohindx]].lit_log12oh
    yerratlas[ohindx] = atlasnodust[indxatlas[ohindx]].lit_log12oh_err

; NFGS

    indxnfgs = where((nfgsnodust.zstrong_12oh_kk04_r23_upper gt -900) and $
      (nfgsnodust.M_B gt -900),nindxnfgs)
;     (nfgsnodust.M_B gt -900) and (nfgsnodust.m_b lt mbcut),nindxnfgs)

    nfgst = nfgsnodust[indxnfgs].lit_t
    nfgstype = nfgsnodust[indxnfgs].lit_type

    xnfgs = nfgsnodust[indxnfgs].M_B
    xerrnfgs = nfgsnodust[indxnfgs].M_B_err

    ynfgs = xnfgs*0.0
    yerrnfgs = xerrnfgs*0.0
    
    lower = where((strmatch(nfgstype,'*Irr*',/fold) eq 1B) and (xnfgs gt mbcut),nlower,comp=upper,ncomp=nupper)

    if (nlower ne 0L) then begin
       ynfgs[lower] = nfgsnodust[indxnfgs[lower]].zstrong_12oh_kk04_r23_lower
       yerrnfgs[lower] = nfgsnodust[indxnfgs[lower]].zstrong_12oh_kk04_r23_lower_err
    endif
    if (nupper ne 0L) then begin
       ynfgs[upper] = nfgsnodust[indxnfgs[upper]].zstrong_12oh_kk04_r23_upper
       yerrnfgs[upper] = nfgsnodust[indxnfgs[upper]].zstrong_12oh_kk04_r23_upper_err
    endif

;   ynfgs = nfgsnodust[indxnfgs].zstrong_12oh_kk04_r23_upper
;   yerrnfgs = nfgsnodust[indxnfgs].zstrong_12oh_kk04_r23_upper_err

    xatlas = [xatlas,xnfgs]
    yatlas = [yatlas,ynfgs]
    type = [atlast,nfgst]

    pec = where(type eq 10.0,comp=normal)
    
    xtitle = 'M_{B}'
    ytitle = ohtitle+' from KK04'
;   ytitle = ohtitle+' ([O III]/H\beta)/([N II]/H\alpha)'

    xrange = MBrange
;   yrange = [8.2,9.5] ; ohrange
    yrange = ohrange7

    sdss_lineplot, x, y, xerr, yerr, plottype=1, postscript=postscript, $
      xtitle=xtitle, ytitle=ytitle, xrange=xrange, yrange=yrange, legendtype=0, $
      /right, /top, position=pos[*,0], charsize=charsize_8, /xreverse, $
;     xatlas=xatlas[normal], xerratlas=xerratlas[normal], yatlas=yatlas[normal], yerratlas=yerratlas[normal] ;, $
      xatlas=xatlas[strong], xerratlas=xerratlas[strong], yatlas=yatlas[strong], yerratlas=yerratlas[strong] ;, $
;     xatlas=xatlas, xerratlas=xerratlas, yatlas=yatlas, yerratlas=yerratlas;, $
;     xnfgs=xnfgs, xerrnfgs=xerrnfgs, ynfgs=ynfgs, yerrnfgs=yerrnfgs

    plotsym, 0, 1, /fill
;   djs_oplot, xatlas[pec], yatlas[pec], color='red', ps=8
    djs_oplot, xatlas[ohindx], yatlas[ohindx], color='red', ps=8
    
; overplot the SDSS Tremonti result

;   massaxis = findgen((11.5-8.5)/0.01+1)*0.01+8.5
;   djs_oplot, massaxis, -1.492 + 1.847*massaxis - 0.08026*massaxis^2, $
;     line=0, thick=postthick, color='dark green'

    MBaxis = findgen(((-16)-(-21.5))/0.01+1)*0.01-21.5
;   djs_oplot, MBaxis+0.168, 5.238 - 0.185*MBaxis, line=0, thick=postthick, color='dark green'

; KZ99    

;   djs_oplot, MBaxis, 6.637 - 0.108*MBaxis, line=0, thick=postthick, color='red'

    sixlin, x, y, a, siga, b, sigb
    cc = [a[2],b[2]] 
    djs_oplot, mbaxis, poly(mbaxis,cc), thick=postthick, color='dark green'

    sixlin, xatlas, yatlas, a, siga, b, sigb
    cc = [a[2],b[2]]
;   djs_oplot, mbaxis, poly(mbaxis,cc), line=0, color='green', thick=postthick

    im_openclose, postscript=postscript, /close    
    
; ------------------------------------------------------------
; Mass vs 12+log (O/H)
; ------------------------------------------------------------

    psname = 'mass_vs_12oh_kk04_r23'
    im_openclose, pspath+psname, postscript=postscript, xsize=8.5, ysize=8.0

    pagemaker, nx=1, ny=1, xspace=0, yspace=0.0, width=6.6, height=6.6, $
      xmargin=[1.5,0.4], ymargin=[0.3,1.1], xpage=8.5, ypage=8.0, $
      position=pos, /normal

    masscut = 8.5

    indx = where((sdssnodust.zstrong_12oh_kk04_r23_upper gt -900) and $
      (sdssancillary.kauffmann_mass gt -900),nindx)
;     (sdssancillary.mass_bv_b gt -900),nindx)

    x = sdssancillary[indx].kauffmann_mass
    xerr = x*0.0
;   x = sdssancillary[indx].mass_bv_b
;   xerr = sdssancillary[indx].mass_bv_b_err

    y = sdssnodust[indx].zstrong_12oh_kk04_r23_upper
    yerr = sdssnodust[indx].zstrong_12oh_kk04_r23_upper_err

; Atlas

    indxatlas = where((atlasnodust.zstrong_12oh_kk04_r23_upper gt -900) and $
      (atlasnodust.continuum_total_mass gt 10^masscut) and $
;     (atlasnodust.mass_bv_b gt -900) and (atlasnodust.mass_bv_b gt masscut) and $
      (strmatch(atlasnodust.lit_type,'*Irr*',/fold) eq 0B),nindxatlas)
;   indxatlas = where((atlasnodust.zstrong_12oh_kk04_r23_upper gt 8.4) and $
;     (atlasnodust.mass_bv_b gt -900) and (atlasnodust.mass_bv_b gt masscut),nindxatlas)
;     (atlasnodust.mass_bv_b gt -900) and (atlasnodust.mass_bv_b lt mbcut),nindxatlas)

    atlast = atlasnodust[indxatlas].lit_t
    atlastype = atlasnodust[indxatlas].lit_type

    xatlas = alog10(atlasnodust[indxatlas].continuum_total_mass)
    xerratlas = xatlas*0.0
;   xatlas = atlasnodust[indxatlas].mass_bv_b
;   xerratlas = atlasnodust[indxatlas].mass_bv_b_err

;   yatlas = xatlas*0.0
;   yerratlas = xerratlas*0.0
;   
;   lower = where((strmatch(atlastype,'*Irr*',/fold) eq 1B) and (xatlas gt mbcut),nlower,comp=upper,ncomp=nupper)
;   if (nlower ne 0L) then begin
;      yatlas[lower] = atlasnodust[indxatlas[lower]].zstrong_12oh_kk04_r23_lower
;      yerratlas[lower] = atlasnodust[indxatlas[lower]].zstrong_12oh_kk04_r23_lower_err
;   endif
;   if (nupper ne 0L) then begin
;      yatlas[upper] = atlasnodust[indxatlas[upper]].zstrong_12oh_kk04_r23_upper
;      yerratlas[upper] = atlasnodust[indxatlas[upper]].zstrong_12oh_kk04_r23_upper_err
;   endif

    yatlas = atlasnodust[indxatlas].zstrong_12oh_kk04_r23_upper
    yerratlas = atlasnodust[indxatlas].zstrong_12oh_kk04_r23_upper_err

;   ohindx = where(atlasnodust[indxatlas].lit_log12oh gt -900.0,nohindx,comp=strong)
;   yatlas[ohindx] = atlasnodust[indxatlas[ohindx]].lit_log12oh
;   yerratlas[ohindx] = atlasnodust[indxatlas[ohindx]].lit_log12oh_err

; NFGS

    indxnfgs = where((nfgsnodust.zstrong_12oh_kk04_r23_upper gt -900.0) and $
      (nfgsnodust.continuum_total_mass gt 10^masscut) and $
;     (nfgsnodust.mass_bv_b gt -900) and (nfgsnodust.mass_bv_b gt masscut) and $
      (strmatch(nfgsnodust.lit_type,'*Irr*',/fold) eq 0B),nindxnfgs)
;   indxnfgs = where((nfgsnodust.zstrong_12oh_kk04_r23_upper gt 8.4) and $
;     (nfgsnodust.mass_bv_b gt -900) and (nfgsnodust.mass_bv_b gt masscut),nindxnfgs)
;     (nfgsnodust.mass_bv_b gt -900) and (nfgsnodust.mass_bv_b lt mbcut),nindxnfgs)

    nfgst = nfgsnodust[indxnfgs].lit_t
    nfgstype = nfgsnodust[indxnfgs].lit_type

    xnfgs = nfgsnodust[indxnfgs].continuum_total_mass
    xerrnfgs = xnfgs*0.0
;   xnfgs = nfgsnodust[indxnfgs].mass_bv_b
;   xerrnfgs = nfgsnodust[indxnfgs].mass_bv_b_err

;   ynfgs = xnfgs*0.0
;   yerrnfgs = xerrnfgs*0.0
;   
;   lower = where((strmatch(nfgstype,'*Irr*',/fold) eq 1B) and (xnfgs gt mbcut),nlower,comp=upper,ncomp=nupper)
;   if (nlower ne 0L) then begin
;      ynfgs[lower] = nfgsnodust[indxnfgs[lower]].zstrong_12oh_kk04_r23_lower
;      yerrnfgs[lower] = nfgsnodust[indxnfgs[lower]].zstrong_12oh_kk04_r23_lower_err
;   endif
;   if (nupper ne 0L) then begin
;      ynfgs[upper] = nfgsnodust[indxnfgs[upper]].zstrong_12oh_kk04_r23_upper
;      yerrnfgs[upper] = nfgsnodust[indxnfgs[upper]].zstrong_12oh_kk04_r23_upper_err
;   endif

    ynfgs = nfgsnodust[indxnfgs].zstrong_12oh_kk04_r23_upper
    yerrnfgs = nfgsnodust[indxnfgs].zstrong_12oh_kk04_r23_upper_err

    xatlas = [xatlas,xnfgs]
    yatlas = [yatlas,ynfgs]
    type = [atlast,nfgst]

    pec = where(type eq 10.0,comp=normal)
    
    xtitle = 'log (M / M'+sunsymbol()+')'
    ytitle = ohtitle+' from KK04'
;   ytitle = ohtitle+' ([O III]/H\beta)/([N II]/H\alpha)'

    xrange = [7.5,12.5] ; massrange
;   yrange = [8.2,9.5] ; ohrange
    yrange = [8.1,9.6] ; ohrange
;   yrange = ohrange7

    sdss_lineplot, x, y, xerr, yerr, plottype=1, postscript=postscript, $
      xtitle=xtitle, ytitle=ytitle, xrange=xrange, yrange=yrange, legendtype=0, $
      /right, /top, position=pos[*,0], charsize=charsize_8, $
;     xatlas=xatlas[normal], xerratlas=xerratlas[normal], yatlas=yatlas[normal], yerratlas=yerratlas[normal] ;, $
;     xatlas=xatlas[strong], xerratlas=xerratlas[strong], yatlas=yatlas[strong], yerratlas=yerratlas[strong] ;, $
      xatlas=xatlas, xerratlas=xerratlas, yatlas=yatlas, yerratlas=yerratlas;, $
;     xnfgs=xnfgs, xerrnfgs=xerrnfgs, ynfgs=ynfgs, yerrnfgs=yerrnfgs

    plotsym, 0, 1, /fill
;   djs_oplot, xatlas[pec], yatlas[pec], color='red', ps=8
;   djs_oplot, xatlas[ohindx], yatlas[ohindx], color='red', ps=8
    
; overplot the SDSS Tremonti result

    offset1 = 0.02 ; <-- NOTE!
    offset2 = 0.0 ; <-- NOTE!
;   offset1 = -0.09 ; <-- NOTE!
;   offset2 = 0.23 ; <-- NOTE!
    
    massaxis = findgen((11.5-8.5)/0.01+1)*0.01+8.5
    djs_oplot, massaxis, -1.492 + 1.847*massaxis - 0.08026*massaxis^2 + offset1, $
      line=0, thick=postthick, color='dark green'

    im_openclose, postscript=postscript, /close    
    
