;+
; NAME:
;       IR SFRS
;
; PURPOSE:
;       Investigate IR star formation rate indicators. 
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
;       J. Moustakas, 2004 Jun 10, U of A
;-

pro ir_sfrs, atlasdust, atlasnodust, nfgsdust, nfgsnodust, atlasdust_agn, $
  atlasnodust_agn, nfgsdust_agn, nfgsnodust_agn, postscript=postscript, $
  paper=paper, cleanpng=cleanpng, _extra=extra

    lsun = 3.826D33              ; [erg/s]
    dsun = 1.496D13              ; Earth-Sun distance [cm]
    light = 2.99792458D18        ; speed of light [A/s]
    mpc2cm = 3.086D24    ; [cm/Mpc]
    haconst = 7.9D-42            ; K98 conversion L(Ha) --> SFR(Ha)
    irconst = 4.5D-44            ; K98 conversion L(IR) --> SFR(IR)
    loghaconst = alog10(haconst)
    logUconst = alog10(7.9D-30)  ; Cram et al. (1998) conversion L(U) --> SFR(U)
    hasfrconstoffset = 41.0
    Uhasfrconstoffset = 43.0
    
    LBnorm = 1E10   ; [L_B_sun]
    sfrnorm = 1.0   ; [M_sun/yr]
    elumnorm = 1D41 ; [erg/s]
    Ulumnorm = 1D43 ; [erg/s]

    syserr = 0.0                ; systematic flux error [%]
;   syserr = 4.0                ; systematic flux error [%]
    nfgssyserr = 6.0
    
    HaHb = 2.86
;   HaHb = return_tbalmer(/HaHb)
    
    Uinfo = im_filterspecs(filterlist='bessell_U.par')
    Uconstant = Uinfo.weff*Uinfo.vega_flam
    
    if keyword_set(paper) then postscript = 1L
    if keyword_set(postscript) then begin
       postthick = 8.0 
       postthick2 = 12.0
    endif else begin
       postthick = 2.0
       postthick2 = 2.0
    endelse

    htmlbase = 'ir_sfrs'

    html_path = atlas_path(/web)+'analysis/'
    pspath = atlas_path(/web)+'analysis/'+htmlbase+'/'
    paperpath = atlas_path(/papers)+'sfrs/FIG_SFRS/'
    latexpath = atlas_path(/papers)+'sfrs/'
    sfrspath = atlas_path(/projects)+'sfrs/'

    if (n_elements(snrcut) eq 0L) then snrcut = 3.0
    snrcut_highz = 0.0

; read the data and the models    
    
    if (n_elements(atlasdust) eq 0L) then atlasdust = read_atlas_sfrs_sample(atlasnodust=atlasnodust)
    if (n_elements(atlasdust_agn) eq 0L) then atlasdust_agn = read_atlas_sfrs_sample(atlasnodust=atlasnodust_agn,/agn)
    if (n_elements(nfgsdust) eq 0L) then nfgsdust = read_nfgs_sfrs_sample(nfgsnodust=nfgsnodust)
    if (n_elements(nfgsdust_agn) eq 0L) then nfgsdust_agn = read_nfgs_sfrs_sample(nfgsnodust=nfgsnodust_agn,/agn)

; initialize plotting variables

    @'xyrange_sfrs'

    Zsun_old = 8.9
    Zsun_new = 8.7

; clean up old PS and PNG files

    if keyword_set(cleanpng) then begin
       splog, 'Deleting all PNG and PS files in '+pspath
       spawn, ['rm -f '+pspath+'*.png*'], /sh
       spawn, ['rm -f '+pspath+'*ps'], /sh
    endif

    if (not keyword_set(postscript)) then im_window, 0, /square

; ###########################################################################    
; Paper Plots
; ###########################################################################    

; ------------------------------------------------------------
; L(B) vs SFR(Ha)_cor/SFR(IR) and SFR(Ha)_cor/[SFR(Ha)_obs+SFR(IR)]
; ------------------------------------------------------------
    
    psname = 'lb_vs_sfr_ha_sfr_hair_2panel'
    im_openclose, pspath+psname, postscript=postscript, /encapsulated, xsize=7.0, ysize=7.5

    pagemaker, nx=1, ny=2, height=[3.0,3.0], width=5.5, xmargin=[1.1,0.4], $
      ymargin=[0.4,1.1], xspace=0, yspace=0, xpage=7.0, ypage=7.5, $
      position=pos, /normal

; SF galaxies    
    
    indx = where((atlasnodust.b_lum_obs gt -900) and (atlasnodust.ir_flux gt -900) and $
      (atlasnodust.ebv_hahb_err gt 0.0),nindx)

    x = atlasdust[indx].b_lum_obs
    xerr = atlasdust[indx].b_lum_obs_err
    
    lir = atlasdust[indx].ir_flux ; [erg/s/cm2]
    lir_err = sqrt(atlasdust[indx].ir_flux_err^2 + (atlasdust[indx].ir_flux*0.15)^2)

    indxnfgs = where((nfgsdust.b_lum_obs gt -900) and (nfgsdust.ir_flux gt -900) and $
      (nfgsnodust.ebv_hahb_err gt 0.0),nindxnfgs)

    xnfgs = nfgsnodust[indxnfgs].b_lum_obs
    xerrnfgs = nfgsnodust[indxnfgs].b_lum_obs_err
    
    lirnfgs = nfgsdust[indxnfgs].ir_flux
    lirnfgs_err = sqrt(nfgsdust[indxnfgs].ir_flux_err^2 + (nfgsdust[indxnfgs].ir_flux*0.15)^2)
    
; AGN
    
    indx_agn = where((atlasnodust_agn.b_lum_obs gt -900) and (atlasnodust_agn.ir_flux gt -900) and $
      (atlasnodust_agn.ebv_hahb_err gt 0.0),nindx_agn)

    x_agn = atlasdust_agn[indx_agn].b_lum_obs
    xerr_agn = atlasdust_agn[indx_agn].b_lum_obs_err
    
    lir_agn = atlasdust_agn[indx_agn].ir_flux ; [erg/s/cm2]
    lir_err_agn = sqrt(atlasdust_agn[indx_agn].ir_flux_err^2 + (atlasdust_agn[indx_agn].ir_flux*0.15)^2)

    indxnfgs_agn = where((nfgsdust_agn.b_lum_obs gt -900) and (nfgsdust_agn.ir_flux gt -900) and $
      (nfgsnodust_agn.ebv_hahb_err gt 0.0),nindxnfgs_agn)

    xnfgs_agn = nfgsnodust_agn[indxnfgs_agn].b_lum_obs
    xerrnfgs_agn = nfgsnodust_agn[indxnfgs_agn].b_lum_obs_err
    
    lirnfgs_agn = nfgsdust_agn[indxnfgs_agn].ir_flux
    lirnfgs_err_agn = sqrt(nfgsdust_agn[indxnfgs_agn].ir_flux_err^2 + (nfgsdust_agn[indxnfgs_agn].ir_flux*0.15)^2)
    
; plotting variables    
    
    xrange = lbrange2
    yrange = sfrHaIRrange

    xtitle = 'log L(B) [L(B)'+sunsymbol()+']'

    lhalir = alog10(4.5D-44/7.9D-42)

; ##########################################################
; Panel 1
; ##########################################################

; SF galaxies    
    
    ha = atlasnodust[indx].h_alpha[0]
    ha_err = atlasnodust[indx].h_alpha[1]
    
    y = alog10(ha/lir) - lhalir
    yerr = im_compute_error(ha,ha_err,lir,lir_err,/log)

    hanfgs = nfgsnodust[indxnfgs].h_alpha[0]
    hanfgs_err = nfgsnodust[indxnfgs].h_alpha[1]
    
    ynfgs = alog10(hanfgs/lirnfgs) - lhalir
    yerrnfgs = im_compute_error(hanfgs,hanfgs_err,lirnfgs,lirnfgs_err,/log)

; AGN    
    
    ha_agn = atlasnodust_agn[indx_agn].h_alpha[0]
    ha_err_agn = atlasnodust_agn[indx_agn].h_alpha[1]
    
    y_agn = alog10(ha_agn/lir_agn) - lhalir
    yerr_agn = im_compute_error(ha_agn,ha_err_agn,lir_agn,lir_err_agn,/log)

    hanfgs_agn = nfgsnodust_agn[indxnfgs_agn].h_alpha[0]
    hanfgs_err_agn = nfgsnodust_agn[indxnfgs_agn].h_alpha[1]
    
    ynfgs_agn = alog10(hanfgs_agn/lirnfgs_agn) - lhalir
    yerrnfgs_agn = im_compute_error(hanfgs_agn,hanfgs_err_agn,lirnfgs_agn,lirnfgs_err_agn,/log)

; combine the two samples - SF

    xbig = [x,xnfgs] & xerrbig = [xerr,xerrnfgs]
    ybig = [y,ynfgs] & yerrbig = [yerr,yerrnfgs]
    
; combine the two samples - AGN

    xbig_agn = [x_agn,xnfgs_agn] & xerrbig_agn = [xerr_agn,xerrnfgs_agn]
    ybig_agn = [y_agn,ynfgs_agn] & yerrbig_agn = [yerr_agn,yerrnfgs_agn]

    stats = im_stats(ybig,/verbose)
    xstr = strtrim(string(stats.median_rej,format='(F12.2)'),2)+' ('+$
      strtrim(string(stats.mean_rej,format='(F12.2)'),2)+$
      '\pm'+strtrim(string(stats.sigma_rej,format='(F12.2)'),2)+')'

    ytitle = 'log \psi(H\alpha)_{cor}/\psi(IR)'

    atlas1d_lineplot, xbig, ybig, xerrbig, yerrbig, plottype=3, postscript=postscript, $
      xtitle='', ytitle=ytitle, xrange=xrange, yrange=yrange, legendtype=0, $
      charsize=charsize_5, position=pos[*,0], yminor=3, atlaspsize=1.0, xtickname=replicate(' ',10)
    atlas1d_lineplot, xbig_agn, ybig_agn, xerrbig_agn, yerrbig_agn, postscript=postscript, $
      atlaspsize=0.8, atlasfill=0, /overplot, atlascolor='red'
    djs_oplot, !x.crange, [0,0], line=0, thick=postthick

    legend, '(a)', /left, /top, box=0, charsize=charsize_3, charthick=postthick
    legend, textoidl(xstr), /right, /top, box=0, charsize=charsize_3, charthick=postthick

; ############################################################
; Panel 2
; ############################################################

; SF galaxies    
    
    ha = atlasnodust[indx].h_alpha[0]
    ha_err = atlasnodust[indx].h_alpha[1]

    haobs = atlasdust[indx].h_alpha[0]
    haobs_err = atlasdust[indx].h_alpha[1]

    y1 = haconst*ha
    y1_err = haconst*ha_err

    y2 = haconst*haobs+irconst*lir
    y2_err = sqrt((haconst*haobs_err)^2+(irconst*lir_err)^2)
    
    y = alog10(y1/y2)
    yerr = im_compute_error(y1,y1_err,y2,y2_err,/log)

    hanfgs = nfgsnodust[indxnfgs].h_alpha[0]
    hanfgs_err = nfgsnodust[indxnfgs].h_alpha[1]
    
    haobsnfgs = nfgsdust[indxnfgs].h_alpha[0]
    haobsnfgs_err = nfgsdust[indxnfgs].h_alpha[1]
    
    y1nfgs = haconst*hanfgs
    y1nfgs_err = haconst*hanfgs_err

    y2nfgs = haconst*haobsnfgs+irconst*lirnfgs
    y2nfgs_err = sqrt((haconst*haobsnfgs_err)^2+(irconst*lirnfgs_err)^2)
    
    ynfgs = alog10(y1nfgs/y2nfgs)
    yerrnfgs = im_compute_error(y1nfgs,y1nfgs_err,y2nfgs,y2nfgs_err,/log)

; AGN    

    ha_agn = atlasnodust_agn[indx_agn].h_alpha[0]
    ha_err_agn = atlasnodust_agn[indx_agn].h_alpha[1]

    haobs_agn = atlasdust_agn[indx_agn].h_alpha[0]
    haobs_err_agn = atlasdust_agn[indx_agn].h_alpha[1]

    y1_agn = haconst*ha_agn
    y1_err_agn = haconst*ha_err_agn

    y2_agn = haconst*haobs_agn+irconst*lir_agn
    y2_err_agn = sqrt((haconst*haobs_err_agn)^2+(irconst*lir_err_agn)^2)
    
    y_agn = alog10(y1_agn/y2_agn)
    yerr_agn = im_compute_error(y1_agn,y1_err_agn,y2_agn,y2_err_agn,/log)

    hanfgs_agn = nfgsnodust_agn[indxnfgs_agn].h_alpha[0]
    hanfgs_err_agn = nfgsnodust_agn[indxnfgs_agn].h_alpha[1]
    
    haobsnfgs_agn = nfgsdust_agn[indxnfgs_agn].h_alpha[0]
    haobsnfgs_err_agn = nfgsdust_agn[indxnfgs_agn].h_alpha[1]
    
    y1nfgs_agn = haconst*hanfgs_agn
    y1nfgs_err_agn = haconst*hanfgs_err_agn

    y2nfgs_agn = haconst*haobsnfgs_agn+irconst*lirnfgs_agn
    y2nfgs_err_agn = sqrt((haconst*haobsnfgs_err_agn)^2+(irconst*lirnfgs_err_agn)^2)
    
    ynfgs_agn = alog10(y1nfgs_agn/y2nfgs_agn)
    yerrnfgs_agn = im_compute_error(y1nfgs_agn,y1nfgs_err_agn,y2nfgs_agn,y2nfgs_err_agn,/log)

; combine the two samples - SF

    xbig = [x,xnfgs] & xerrbig = [xerr,xerrnfgs]
    ybig = [y,ynfgs] & yerrbig = [yerr,yerrnfgs]
    
; combine the two samples - AGN

    xbig_agn = [x_agn,xnfgs_agn] & xerrbig_agn = [xerr_agn,xerrnfgs_agn]
    ybig_agn = [y_agn,ynfgs_agn] & yerrbig_agn = [yerr_agn,yerrnfgs_agn]

    stats = im_stats(ybig,/verbose)
    xstr = strtrim(string(stats.median_rej,format='(F12.2)'),2)+' ('+$
      strtrim(string(stats.mean_rej,format='(F12.2)'),2)+$
      '\pm'+strtrim(string(stats.sigma_rej,format='(F12.2)'),2)+')'

    ytitle = 'log {\psi(H\alpha)_{cor}/[\psi(H\alpha)_{obs}+\psi(IR)]}'

    atlas1d_lineplot, xbig, ybig, xerrbig, yerrbig, plottype=3, postscript=postscript, $
      xtitle=xtitle, ytitle=ytitle, xrange=xrange, yrange=yrange, legendtype=0, $
      charsize=charsize_5, position=pos[*,1], yminor=3, atlaspsize=1.0, /noerase
    atlas1d_lineplot, xbig_agn, ybig_agn, xerrbig_agn, yerrbig_agn, postscript=postscript, $
      atlaspsize=0.8, atlasfill=0, /overplot, atlascolor='red'
;   djs_oplot, !x.crange, [0,0], line=0, thick=postthick

    legend, '(b)', /left, /top, box=0, charsize=charsize_3, charthick=postthick
    legend, textoidl(xstr), /right, /top, box=0, charsize=charsize_3, charthick=postthick

;   xyouts, pos[0]*0.5, (pos[2,0]-pos[1,1])/2.0+pos[1,1], textoidl(ytitle), align=0.5, $
;     orientation=90, charsize=charsize_7, charthick=postthick, /normal
    
    im_openclose, postscript=postscript, /close

; ------------------------------------------------------------
; L(B) vs L(12)/[SFR(Ha)_obs+SFR(IR)]
; ------------------------------------------------------------
    
    psname = 'lb_vs_l12_sfr_hair'
    im_openclose, pspath+psname, postscript=postscript, xsize=8.5, ysize=8.0, /encapsulated

    pagemaker, nx=1, ny=1, xspace=0, yspace=0.0, width=6.6, height=6.6, $
      xmargin=[1.5,0.4], ymargin=[0.3,1.1], xpage=8.5, ypage=8.0, $
      position=pos, /normal

; SF galaxies
    
    indx = where((atlasdust.b_lum_obs gt -900) and (atlasdust.ir_flux gt -900) and $
      (atlasdust.iras_12 gt -900) and (atlasnodust.ebv_hahb_err gt 0.0),nindx)

    x = atlasdust[indx].b_lum_obs
    xerr = atlasdust[indx].b_lum_obs_err
    
    l12 = atlasdust[indx].iras_12*1D-23*light/12E4^2 ; [erg/s/cm2]
    l12_err = atlasdust[indx].iras_12_err*1D-23*light/12E4^2

    lir = atlasdust[indx].ir_flux ; [erg/s/cm2]
    lir_err = sqrt(atlasdust[indx].ir_flux_err^2 + (atlasdust[indx].ir_flux*0.15)^2)

    indxnfgs = where((nfgsdust.b_lum_obs gt -900) and (nfgsdust.ir_flux gt -900) and $
      (nfgsdust.iras_12 gt -900) and (nfgsnodust.ebv_hahb_err gt 0.0),nindxnfgs)

    xnfgs = nfgsnodust[indxnfgs].b_lum_obs
    xerrnfgs = nfgsnodust[indxnfgs].b_lum_obs_err
    
    l12nfgs = nfgsdust[indxnfgs].iras_12*1D-23*light/12E4^2 ; [erg/s/cm2]
    l12nfgs_err = nfgsdust[indxnfgs].iras_12_err*1D-23*light/12E4^2

    lirnfgs = nfgsdust[indxnfgs].ir_flux
    lirnfgs_err = sqrt(nfgsdust[indxnfgs].ir_flux_err^2 + (nfgsdust[indxnfgs].ir_flux*0.15)^2)
    
; AGN
    
    indx_agn = where((atlasdust_agn.b_lum_obs gt -900) and (atlasdust_agn.ir_flux gt -900) and $
      (atlasdust_agn.iras_12 gt -900) and (atlasnodust_agn.ebv_hahb_err gt 0.0),nindx_agn)

    x_agn = atlasdust_agn[indx_agn].b_lum_obs
    xerr_agn = atlasdust_agn[indx_agn].b_lum_obs_err
    
    l12_agn = atlasdust_agn[indx_agn].iras_12*1D-23*light/12E4^2 ; [erg/s/cm2]
    l12_err_agn = atlasdust_agn[indx_agn].iras_12_err*1D-23*light/12E4^2

    lir_agn = atlasdust_agn[indx_agn].ir_flux
    lir_err_agn = atlasdust_agn[indx_agn].ir_flux_err

    indxnfgs_agn = where((nfgsdust_agn.b_lum_obs gt -900) and (nfgsdust_agn.ir_flux gt -900) and $
      (nfgsdust_agn.iras_12 gt -900) and (nfgsnodust_agn.ebv_hahb_err gt 0.0),nindxnfgs_agn)

    xnfgs_agn = nfgsdust_agn[indxnfgs_agn].b_lum_obs
    xerrnfgs_agn = nfgsdust_agn[indxnfgs_agn].b_lum_obs_err
    
    l12nfgs_agn = nfgsdust_agn[indxnfgs_agn].iras_12*1D-23*light/12E4^2 ; [erg/s/cm2]
    l12nfgs_err_agn = nfgsdust_agn[indxnfgs_agn].iras_12_err*1D-23*light/12E4^2

    lirnfgs_agn = nfgsdust_agn[indxnfgs_agn].ir_flux
    lirnfgs_err_agn = nfgsdust_agn[indxnfgs_agn].ir_flux_err

; plotting variables    
    
    xrange = lbrange2
    yrange = [-1.1,1.3]

    xtitle = 'log L(B) [L(B)'+sunsymbol()+']'
    ytitle = 'log {10^{-37} L(12\mu'+'m'+') / [\psi(H\alpha)_{obs}+\psi(IR)]}'

; SF galaxies    

    y1 = haconst*atlasdust[indx].h_alpha[0]+irconst*lir
    y1_err = sqrt((haconst*atlasdust[indx].h_alpha[1])^2+(irconst*lir_err)^2)
    
    y = alog10(l12/y1) - 37.0
    yerr = im_compute_error(l12,l12_err,y1,y1_err,/log)

    y1nfgs = haconst*nfgsdust[indxnfgs].h_alpha[0]+irconst*lirnfgs
    y1nfgs_err = sqrt((haconst*nfgsdust[indxnfgs].h_alpha[1])^2+(irconst*lirnfgs_err)^2)
    
    ynfgs = alog10(l12nfgs/y1nfgs) - 37.0
    yerrnfgs = im_compute_error(l12nfgs,l12nfgs_err,y1nfgs,y1nfgs_err,/log)

; AGN
    
    y1_agn = haconst*atlasdust_agn[indx_agn].h_alpha[0]+irconst*lir_agn
    y1_err_agn = sqrt((haconst*atlasdust_agn[indx_agn].h_alpha[1])^2+(irconst*lir_err_agn)^2)
    
    y_agn = alog10(l12_agn/y1_agn) - 37.0
    yerr_agn = im_compute_error(l12_agn,l12_err_agn,y1_agn,y1_err_agn,/log)

    y1nfgs_agn = haconst*nfgsdust_agn[indxnfgs_agn].h_alpha[0]+irconst*lirnfgs_agn - 37.0
    y1nfgs_err_agn = sqrt((haconst*nfgsdust_agn[indxnfgs_agn].h_alpha[1])^2+(irconst*lirnfgs_err_agn)^2)
    
    ynfgs_agn = alog10(l12nfgs_agn/y1nfgs_agn)
    yerrnfgs_agn = im_compute_error(l12nfgs_agn,l12nfgs_err_agn,y1nfgs_agn,y1nfgs_err_agn,/log)

; combine the two samples - SF

    xbig = [x,xnfgs] & xerrbig = [xerr,xerrnfgs]
    ybig = [y,ynfgs] & yerrbig = [yerr,yerrnfgs]
    
; combine the two samples - AGN

    xbig_agn = [x_agn,xnfgs_agn] & xerrbig_agn = [xerr_agn,xerrnfgs_agn]
    ybig_agn = [y_agn,ynfgs_agn] & yerrbig_agn = [yerr_agn,yerrnfgs_agn]

    stats = im_stats(ybig,/verbose)
    xstr = strtrim(string(stats.median_rej,format='(F12.2)'),2)+' ('+$
      strtrim(string(stats.mean_rej,format='(F12.2)'),2)+$
      '\pm'+strtrim(string(stats.sigma_rej,format='(F12.2)'),2)+')'

    atlas1d_lineplot, xbig, ybig, xerrbig, yerrbig, plottype=3, postscript=postscript, $
      xtitle=xtitle, ytitle=ytitle, xrange=xrange, yrange=yrange, legendtype=0, $
      charsize=charsize_7, position=pos[*,0], yminor=3, atlaspsize=1.0
    atlas1d_lineplot, xbig_agn, ybig_agn, xerrbig_agn, yerrbig_agn, postscript=postscript, $
      atlaspsize=0.8, atlasfill=0, /overplot, atlascolor='red'

    legend, textoidl(xstr), /right, /top, box=0, charsize=charsize_4, charthick=postthick
    
    im_openclose, postscript=postscript, /close

; ------------------------------------------------------------
; L(B) vs L(25)/[SFR(Ha)_obs+SFR(IR)]
; ------------------------------------------------------------
    
    psname = 'lb_vs_l25_sfr_hair'
    im_openclose, pspath+psname, postscript=postscript, xsize=8.5, ysize=8.0, /encapsulated

    pagemaker, nx=1, ny=1, xspace=0, yspace=0.0, width=6.6, height=6.6, $
      xmargin=[1.5,0.4], ymargin=[0.3,1.1], xpage=8.5, ypage=8.0, $
      position=pos, /normal

; SF galaxies
    
    indx = where((atlasdust.b_lum_obs gt -900) and (atlasdust.ir_flux gt -900) and $
      (atlasdust.iras_25 gt -900) and (atlasnodust.ebv_hahb_err gt 0.0),nindx)

    x = atlasdust[indx].b_lum_obs
    xerr = atlasdust[indx].b_lum_obs_err
    
    l25 = atlasdust[indx].iras_25*1D-23*light/25E4^2 ; [erg/s/cm2]
    l25_err = atlasdust[indx].iras_25_err*1D-23*light/25E4^2

    lir = atlasdust[indx].ir_flux ; [erg/s/cm2]
    lir_err = sqrt(atlasdust[indx].ir_flux_err^2 + (atlasdust[indx].ir_flux*0.15)^2)

    indxnfgs = where((nfgsdust.b_lum_obs gt -900) and (nfgsdust.ir_flux gt -900) and $
      (nfgsdust.iras_25 gt -900) and (nfgsnodust.ebv_hahb_err gt 0.0),nindxnfgs)

    xnfgs = nfgsnodust[indxnfgs].b_lum_obs
    xerrnfgs = nfgsnodust[indxnfgs].b_lum_obs_err
    
    l25nfgs = nfgsdust[indxnfgs].iras_25*1D-23*light/25E4^2 ; [erg/s/cm2]
    l25nfgs_err = nfgsdust[indxnfgs].iras_25_err*1D-23*light/25E4^2

    lirnfgs = nfgsdust[indxnfgs].ir_flux
    lirnfgs_err = sqrt(nfgsdust[indxnfgs].ir_flux_err^2 + (nfgsdust[indxnfgs].ir_flux*0.15)^2)
    
; AGN
    
    indx_agn = where((atlasdust_agn.b_lum_obs gt -900) and (atlasdust_agn.ir_flux gt -900) and $
      (atlasdust_agn.iras_25 gt -900) and (atlasnodust_agn.ebv_hahb_err gt 0.0),nindx_agn)

    x_agn = atlasdust_agn[indx_agn].b_lum_obs
    xerr_agn = atlasdust_agn[indx_agn].b_lum_obs_err
    
    l25_agn = atlasdust_agn[indx_agn].iras_25*1D-23*light/25E4^2 ; [erg/s/cm2]
    l25_err_agn = atlasdust_agn[indx_agn].iras_25_err*1D-23*light/25E4^2

    lir_agn = atlasdust_agn[indx_agn].ir_flux
    lir_err_agn = atlasdust_agn[indx_agn].ir_flux_err

    indxnfgs_agn = where((nfgsdust_agn.b_lum_obs gt -900) and (nfgsdust_agn.ir_flux gt -900) and $
      (nfgsdust_agn.iras_25 gt -900) and (nfgsnodust_agn.ebv_hahb_err gt 0.0),nindxnfgs_agn)

    xnfgs_agn = nfgsdust_agn[indxnfgs_agn].b_lum_obs
    xerrnfgs_agn = nfgsdust_agn[indxnfgs_agn].b_lum_obs_err
    
    l25nfgs_agn = nfgsdust_agn[indxnfgs_agn].iras_25*1D-23*light/25E4^2 ; [erg/s/cm2]
    l25nfgs_err_agn = nfgsdust_agn[indxnfgs_agn].iras_25_err*1D-23*light/25E4^2

    lirnfgs_agn = nfgsdust_agn[indxnfgs_agn].ir_flux
    lirnfgs_err_agn = nfgsdust_agn[indxnfgs_agn].ir_flux_err

; plotting variables    
    
    xrange = lbrange2
    yrange = [-1.1,1.3]

    xtitle = 'log L(B) [L(B)'+sunsymbol()+']'
    ytitle = 'log {10^{-37} L(25\mu'+'m'+') / [\psi(H\alpha)_{obs}+\psi(IR)]}'

; SF galaxies    

    y1 = haconst*atlasdust[indx].h_alpha[0]+irconst*lir
    y1_err = sqrt((haconst*atlasdust[indx].h_alpha[1])^2+(irconst*lir_err)^2)
    
    y = alog10(l25/y1) - 37.0
    yerr = im_compute_error(l25,l25_err,y1,y1_err,/log)

    y1nfgs = haconst*nfgsdust[indxnfgs].h_alpha[0]+irconst*lirnfgs
    y1nfgs_err = sqrt((haconst*nfgsdust[indxnfgs].h_alpha[1])^2+(irconst*lirnfgs_err)^2)
    
    ynfgs = alog10(l25nfgs/y1nfgs) - 37.0
    yerrnfgs = im_compute_error(l25nfgs,l25nfgs_err,y1nfgs,y1nfgs_err,/log)

; AGN
    
    y1_agn = haconst*atlasdust_agn[indx_agn].h_alpha[0]+irconst*lir_agn
    y1_err_agn = sqrt((haconst*atlasdust_agn[indx_agn].h_alpha[1])^2+(irconst*lir_err_agn)^2)
    
    y_agn = alog10(l25_agn/y1_agn) - 37.0
    yerr_agn = im_compute_error(l25_agn,l25_err_agn,y1_agn,y1_err_agn,/log)

    y1nfgs_agn = haconst*nfgsdust_agn[indxnfgs_agn].h_alpha[0]+irconst*lirnfgs_agn - 37.0
    y1nfgs_err_agn = sqrt((haconst*nfgsdust_agn[indxnfgs_agn].h_alpha[1])^2+(irconst*lirnfgs_err_agn)^2)
    
    ynfgs_agn = alog10(l25nfgs_agn/y1nfgs_agn)
    yerrnfgs_agn = im_compute_error(l25nfgs_agn,l25nfgs_err_agn,y1nfgs_agn,y1nfgs_err_agn,/log)

; combine the two samples - SF

    xbig = [x,xnfgs] & xerrbig = [xerr,xerrnfgs]
    ybig = [y,ynfgs] & yerrbig = [yerr,yerrnfgs]
    
; combine the two samples - AGN

    xbig_agn = [x_agn,xnfgs_agn] & xerrbig_agn = [xerr_agn,xerrnfgs_agn]
    ybig_agn = [y_agn,ynfgs_agn] & yerrbig_agn = [yerr_agn,yerrnfgs_agn]

    stats = im_stats(ybig,/verbose)
    xstr = strtrim(string(stats.median_rej,format='(F12.2)'),2)+' ('+$
      strtrim(string(stats.mean_rej,format='(F12.2)'),2)+$
      '\pm'+strtrim(string(stats.sigma_rej,format='(F12.2)'),2)+')'

    atlas1d_lineplot, xbig, ybig, xerrbig, yerrbig, plottype=3, postscript=postscript, $
      xtitle=xtitle, ytitle=ytitle, xrange=xrange, yrange=yrange, legendtype=0, $
      charsize=charsize_7, position=pos[*,0], yminor=3, atlaspsize=1.0
    atlas1d_lineplot, xbig_agn, ybig_agn, xerrbig_agn, yerrbig_agn, postscript=postscript, $
      atlaspsize=0.8, atlasfill=0, /overplot, atlascolor='red'

    legend, textoidl(xstr), /right, /top, box=0, charsize=charsize_4, charthick=postthick
    
    im_openclose, postscript=postscript, /close

; ------------------------------------------------------------
; 2-panel - L(B) vs SFR(Ha)/SFR(IR)
; ------------------------------------------------------------
    
    psname = 'lb_vs_sfr_ha_sfr_ir_2panel'
    im_openclose, pspath+psname, postscript=postscript, /encapsulated, xsize=7.0, ysize=7.5

    pagemaker, nx=1, ny=2, height=[3.0,3.0], width=5.5, xmargin=[1.1,0.4], $
      ymargin=[0.4,1.1], xspace=0, yspace=0, xpage=7.0, ypage=7.5, $
      position=pos, /normal

; SF galaxies    
    
    indx = where((atlasnodust.b_lum_obs gt -900) and (atlasnodust.ir_flux gt -900) and $
      (atlasnodust.ebv_hahb_err gt 0.0),nindx)

    x = atlasdust[indx].b_lum_obs
    xerr = atlasdust[indx].b_lum_obs_err
    
    lir = atlasdust[indx].ir_flux ; [erg/s/cm2]
    lir_err = sqrt(atlasdust[indx].ir_flux_err^2 + (atlasdust[indx].ir_flux*0.15)^2)

    indxnfgs = where((nfgsdust.b_lum_obs gt -900) and (nfgsdust.ir_flux gt -900) and (nfgsnodust.ebv_hahb_err gt 0.0),nindxnfgs)

    xnfgs = nfgsnodust[indxnfgs].b_lum_obs
    xerrnfgs = nfgsnodust[indxnfgs].b_lum_obs_err
    
    lirnfgs = nfgsdust[indxnfgs].ir_flux
    lirnfgs_err = sqrt(nfgsdust[indxnfgs].ir_flux_err^2 + (nfgsdust[indxnfgs].ir_flux*0.15)^2)
    
; AGN
    
    indx_agn = where((atlasnodust_agn.b_lum_obs gt -900) and (atlasnodust_agn.ir_flux gt -900) and $
      (atlasnodust_agn.ebv_hahb_err gt 0.0),nindx_agn)

    x_agn = atlasdust_agn[indx_agn].b_lum_obs
    xerr_agn = atlasdust_agn[indx_agn].b_lum_obs_err
    
    lir_agn = atlasdust_agn[indx_agn].ir_flux ; [erg/s/cm2]
    lir_err_agn = sqrt(atlasdust_agn[indx_agn].ir_flux_err^2 + (atlasdust_agn[indx_agn].ir_flux*0.15)^2)

    indxnfgs_agn = where((nfgsdust_agn.b_lum_obs gt -900) and (nfgsdust_agn.ir_flux gt -900) and $
      (nfgsnodust_agn.ebv_hahb_err gt 0.0),nindxnfgs_agn)

    xnfgs_agn = nfgsnodust_agn[indxnfgs_agn].b_lum_obs
    xerrnfgs_agn = nfgsnodust_agn[indxnfgs_agn].b_lum_obs_err
    
    lirnfgs_agn = nfgsdust_agn[indxnfgs_agn].ir_flux
    lirnfgs_err_agn = sqrt(nfgsdust_agn[indxnfgs_agn].ir_flux_err^2 + (nfgsdust_agn[indxnfgs_agn].ir_flux*0.15)^2)
    
; plotting variables    
    
    xrange = lbrange2
    yrange = sfrHaIRrange

    xtitle = 'log L(B) [L(B)'+sunsymbol()+']'

    lhalir = alog10(4.5D-44/7.9D-42)

; ##########################################################
; Panel 1: No reddening correction
; ##########################################################

; SF galaxies    
    
    ha = atlasdust[indx].h_alpha[0]
    ha_err = atlasdust[indx].h_alpha[1]
    
    y = alog10(ha/lir) - lhalir
    yerr = im_compute_error(ha,ha_err,lir,lir_err,/log)

    hanfgs = nfgsdust[indxnfgs].h_alpha[0]
    hanfgs_err = nfgsdust[indxnfgs].h_alpha[1]
    
    ynfgs = alog10(hanfgs/lirnfgs) - lhalir
    yerrnfgs = im_compute_error(hanfgs,hanfgs_err,lirnfgs,lirnfgs_err,/log)

; AGN

    ha_agn = atlasdust_agn[indx_agn].h_alpha[0]
    ha_err_agn = atlasdust_agn[indx_agn].h_alpha[1]
    
    y_agn = alog10(ha_agn/lir_agn) - lhalir
    yerr_agn = im_compute_error(ha_agn,ha_err_agn,lir_agn,lir_err_agn,/log)

    hanfgs_agn = nfgsdust_agn[indxnfgs_agn].h_alpha[0]
    hanfgs_err_agn = nfgsdust_agn[indxnfgs_agn].h_alpha[1]
    
    ynfgs_agn = alog10(hanfgs_agn/lirnfgs_agn) - lhalir
    yerrnfgs_agn = im_compute_error(hanfgs_agn,hanfgs_err_agn,lirnfgs_agn,lirnfgs_err_agn,/log)

; combine the two samples - SF

    xbig = [x,xnfgs] & xerrbig = [xerr,xerrnfgs]
    ybig = [y,ynfgs] & yerrbig = [yerr,yerrnfgs]
    
; combine the two samples - AGN

    xbig_agn = [x_agn,xnfgs_agn] & xerrbig_agn = [xerr_agn,xerrnfgs_agn]
    ybig_agn = [y_agn,ynfgs_agn] & yerrbig_agn = [yerr_agn,yerrnfgs_agn]

    stats = im_stats(ybig,/verbose)
    xstr = strtrim(string(stats.median_rej,format='(F12.2)'),2)+' ('+$
      strtrim(string(stats.mean_rej,format='(F12.2)'),2)+$
      '\pm'+strtrim(string(stats.sigma_rej,format='(F12.2)'),2)+')'

    ytitle = 'log \psi(H\alpha)_{obs}/\psi(IR)'

    atlas1d_lineplot, xbig, ybig, xerrbig, yerrbig, plottype=3, postscript=postscript, $
      xtitle='', ytitle=ytitle, xrange=xrange, yrange=yrange, legendtype=0, $
      charsize=charsize_5, position=pos[*,0], yminor=3, atlaspsize=1.0, xtickname=replicate(' ',10)
    atlas1d_lineplot, xbig_agn, ybig_agn, xerrbig_agn, yerrbig_agn, postscript=postscript, $
      atlaspsize=0.8, atlasfill=0, /overplot, atlascolor='red'
    djs_oplot, !x.crange, [0,0], line=0, thick=postthick

    legend, '(a)', /left, /top, box=0, charsize=charsize_3, charthick=postthick
    legend, textoidl(xstr), /right, /top, box=0, charsize=charsize_3, charthick=postthick

; ############################################################
; Panel 2: Good Reddening Correction
; ############################################################

; SF galaxies    
    
    ha = atlasnodust[indx].h_alpha[0]
    ha_err = atlasnodust[indx].h_alpha[1]

    y = alog10(ha/lir) - lhalir
    yerr = im_compute_error(ha,ha_err,lir,lir_err,/log)

    hanfgs = nfgsnodust[indxnfgs].h_alpha[0]
    hanfgs_err = nfgsnodust[indxnfgs].h_alpha[1]

    ynfgs = alog10(hanfgs/lirnfgs) - lhalir
    yerrnfgs = im_compute_error(hanfgs,hanfgs_err,lirnfgs,lirnfgs_err,/log)

; AGN

    ha_agn = atlasnodust_agn[indx_agn].h_alpha[0]
    ha_err_agn = atlasnodust_agn[indx_agn].h_alpha[1]

    y_agn = alog10(ha_agn/lir_agn) - lhalir
    yerr_agn = im_compute_error(ha_agn,ha_err_agn,lir_agn,lir_err_agn,/log)

    hanfgs_agn = nfgsnodust_agn[indxnfgs_agn].h_alpha[0]
    hanfgs_err_agn = nfgsnodust_agn[indxnfgs_agn].h_alpha[1]

    ynfgs_agn = alog10(hanfgs_agn/lirnfgs_agn) - lhalir
    yerrnfgs_agn = im_compute_error(hanfgs_agn,hanfgs_err_agn,lirnfgs_agn,lirnfgs_err_agn,/log)

; combine the two samples - SF

    xbig = [x,xnfgs] & xerrbig = [xerr,xerrnfgs]
    ybig = [y,ynfgs] & yerrbig = [yerr,yerrnfgs]
    
; combine the two samples - AGN

    xbig_agn = [x_agn,xnfgs_agn] & xerrbig_agn = [xerr_agn,xerrnfgs_agn]
    ybig_agn = [y_agn,ynfgs_agn] & yerrbig_agn = [yerr_agn,yerrnfgs_agn]

    stats = im_stats(ybig,/verbose)
    xstr = strtrim(string(stats.median_rej,format='(F12.2)'),2)+' ('+$
      strtrim(string(stats.mean_rej,format='(F12.2)'),2)+$
      '\pm'+strtrim(string(stats.sigma_rej,format='(F12.2)'),2)+')'

    ytitle = 'log \psi(H\alpha)_{cor}/\psi(IR)'

    atlas1d_lineplot, xbig, ybig, xerrbig, yerrbig, plottype=3, postscript=postscript, $
      xtitle=xtitle, ytitle=ytitle, xrange=xrange, yrange=yrange, legendtype=0, $
      charsize=charsize_5, position=pos[*,1], yminor=3, atlaspsize=1.0, /noerase
    atlas1d_lineplot, xbig_agn, ybig_agn, xerrbig_agn, yerrbig_agn, postscript=postscript, $
      atlaspsize=0.8, atlasfill=0, /overplot, atlascolor='red'
    djs_oplot, !x.crange, [0,0], line=0, thick=postthick

    legend, '(b)', /left, /top, box=0, charsize=charsize_3, charthick=postthick
    legend, textoidl(xstr), /right, /top, box=0, charsize=charsize_3, charthick=postthick
    
    im_openclose, postscript=postscript, /close

; ------------------------------------------------------------
; EW(Ha) vs L(Ha)/L(12) - 2-panel
; ------------------------------------------------------------
    
    psname = 'ewha_vs_lha_l12'
    im_openclose, pspath+psname, postscript=postscript, xsize=8.5, ysize=8.0, /encapsulated

    pagemaker, nx=1, ny=1, xspace=0, yspace=0.0, width=6.6, height=6.6, $
      xmargin=[1.5,0.4], ymargin=[0.3,1.1], xpage=8.5, ypage=8.0, $
      position=pos, /normal

; SF galaxies
    
    indx = where((atlasdust.h_alpha_ew[0]/atlasdust.h_alpha_ew[1] gt snrcut) and (atlasdust.iras_12 gt -900) and $
      (atlasnodust.ebv_hahb_err gt 0.0),nindx)

    x = atlasdust[indx].h_alpha_ew[0]
    xerr = atlasdust[indx].h_alpha_ew[1]
    
    lir = atlasdust[indx].iras_12*1D-23*light/12E4^2 ; [erg/s/cm2]
    lir_err = atlasdust[indx].iras_12_err*1D-23*light/12E4^2

    indxnfgs = where((nfgsdust.h_alpha_ew[0]/nfgsdust.h_alpha_ew[1] gt snrcut) and (nfgsdust.iras_12 gt -900) and $
      (nfgsnodust.ebv_hahb_err gt 0.0),nindxnfgs)

    xnfgs = nfgsdust[indxnfgs].h_alpha_ew[0]
    xerrnfgs = nfgsdust[indxnfgs].h_alpha_ew[1]
    
    lirnfgs = nfgsdust[indxnfgs].iras_12*1D-23*light/12E4^2 ; [erg/s/cm2]
    lirnfgs_err = nfgsdust[indxnfgs].iras_12_err*1D-23*light/12E4^2

; AGN
    
    indx_agn = where((atlasdust_agn.h_alpha_ew[0]/atlasdust_agn.h_alpha_ew[1] gt snrcut) and $
      (atlasdust_agn.iras_12 gt -900) and (atlasnodust_agn.ebv_hahb_err gt 0.0),nindx_agn)

    x_agn = atlasdust_agn[indx_agn].h_alpha_ew[0]
    xerr_agn = atlasdust_agn[indx_agn].h_alpha_ew[1]
    
    lir_agn = atlasdust_agn[indx_agn].iras_12*1D-23*light/12E4^2 ; [erg/s/cm2]
    lir_err_agn = atlasdust_agn[indx_agn].iras_12_err*1D-23*light/12E4^2

    indxnfgs_agn = where((nfgsdust_agn.h_alpha_ew[0]/nfgsdust_agn.h_alpha_ew[1] gt -900) and $
      (nfgsdust_agn.iras_12 gt -900) and (nfgsnodust_agn.ebv_hahb_err gt 0.0),nindxnfgs_agn)

    xnfgs_agn = nfgsdust_agn[indxnfgs_agn].h_alpha_ew[0]
    xerrnfgs_agn = nfgsdust_agn[indxnfgs_agn].h_alpha_ew[1]
    
    lirnfgs_agn = nfgsdust_agn[indxnfgs_agn].iras_12*1D-23*light/12E4^2 ; [erg/s/cm2]
    lirnfgs_err_agn = nfgsdust_agn[indxnfgs_agn].iras_12_err*1D-23*light/12E4^2

; plotting variables    
    
    xrange = ewharangenolog
    yrange = [1.7,5.3]

    xtitle = 'EW(H\alpha) ['+angstrom()+']'
    ytitle = 'log [L(H\alpha)_{cor}/L(12\mu'+'m'+')]'

; SF galaxies    
    
    ha = atlasnodust[indx].h_alpha[0]
    ha_err = atlasnodust[indx].h_alpha[1]

    y = alog10(ha/lir)
    yerr = im_compute_error(ha,ha_err,lir,lir_err,/log)

    hanfgs = nfgsnodust[indxnfgs].h_alpha[0]
    hanfgs_err = nfgsnodust[indxnfgs].h_alpha[1]

    ynfgs = alog10(hanfgs/lirnfgs)
    yerrnfgs = im_compute_error(hanfgs,hanfgs_err,lirnfgs,lirnfgs_err,/log)

; AGN
    
    ha_agn = atlasnodust_agn[indx_agn].h_alpha[0]
    ha_err_agn = atlasnodust_agn[indx_agn].h_alpha[1]
    
    y_agn = alog10(ha_agn/lir_agn)
    yerr_agn = im_compute_error(ha_agn,ha_err_agn,lir_agn,lir_err_agn,/log)

    hanfgs_agn = nfgsnodust_agn[indxnfgs_agn].h_alpha[0]
    hanfgs_err_agn = nfgsnodust_agn[indxnfgs_agn].h_alpha[1]
    
    ynfgs_agn = alog10(hanfgs_agn/lirnfgs_agn)
    yerrnfgs_agn = im_compute_error(hanfgs_agn,hanfgs_err_agn,lirnfgs_agn,lirnfgs_err_agn,/log)

; combine the two samples - SF

    xbig = [x,xnfgs] & xerrbig = [xerr,xerrnfgs]
    ybig = [y,ynfgs] & yerrbig = [yerr,yerrnfgs]
    
; combine the two samples - AGN

    xbig_agn = [x_agn,xnfgs_agn] & xerrbig_agn = [xerr_agn,xerrnfgs_agn]
    ybig_agn = [y_agn,ynfgs_agn] & yerrbig_agn = [yerr_agn,yerrnfgs_agn]

    stats = im_stats(ybig,/verbose)
    xstr = strtrim(string(stats.median_rej,format='(F12.2)'),2)+' ('+$
      strtrim(string(stats.mean_rej,format='(F12.2)'),2)+$
      '\pm'+strtrim(string(stats.sigma_rej,format='(F12.2)'),2)+')'

    atlas1d_lineplot, xbig, ybig, xerrbig, yerrbig, plottype=1, postscript=postscript, $
      xtitle=xtitle, ytitle=ytitle, xrange=xrange, yrange=yrange, legendtype=0, $
      charsize=charsize_7, position=pos[*,0], yminor=3, atlaspsize=1.0, /xlog
    atlas1d_lineplot, xbig_agn, ybig_agn, xerrbig_agn, yerrbig_agn, postscript=postscript, $
      atlaspsize=0.8, atlasfill=0, /overplot, atlascolor='red'

    im_openclose, postscript=postscript, /close

; ------------------------------------------------------------
; 12 + log (O/H) vs L(Ha)/L(12) - 2-panel
; ------------------------------------------------------------
    
    psname = '12oh_oiiinii_niiha_vs_lha_l12'
    im_openclose, pspath+psname, postscript=postscript, xsize=8.5, ysize=8.0, /encapsulated

    pagemaker, nx=1, ny=1, xspace=0, yspace=0.0, width=6.6, height=6.6, $
      xmargin=[1.5,0.4], ymargin=[0.3,1.1], xpage=8.5, ypage=8.0, $
      position=pos, /normal

; SF galaxies
    
    indx = where((atlasnodust.zstrong_12oh_oiiinii_niiha gt -900) and (atlasnodust.iras_12 gt -900) and $
      (atlasnodust.ebv_hahb_err gt 0.0),nindx)

    x = atlasdust[indx].zstrong_12oh_oiiinii_niiha
    xerr = atlasdust[indx].zstrong_12oh_oiiinii_niiha_err
    
    lir = atlasdust[indx].iras_12*1D-23*light/12E4^2 ; [erg/s/cm2]
    lir_err = atlasdust[indx].iras_12_err*1D-23*light/12E4^2

    indxnfgs = where((nfgsdust.zstrong_12oh_oiiinii_niiha gt -900) and (nfgsnodust.iras_12 gt -900) and $
      (nfgsnodust.ebv_hahb_err gt 0.0),nindxnfgs)

    xnfgs = nfgsnodust[indxnfgs].zstrong_12oh_oiiinii_niiha
    xerrnfgs = nfgsnodust[indxnfgs].zstrong_12oh_oiiinii_niiha_err
    
    lirnfgs = nfgsdust[indxnfgs].iras_12*1D-23*light/12E4^2 ; [erg/s/cm2]
    lirnfgs_err = nfgsdust[indxnfgs].iras_12_err*1D-23*light/12E4^2

; AGN
    
    indx_agn = where((atlasnodust_agn.zstrong_12oh_oiiinii_niiha gt -900) and (atlasnodust_agn.iras_12 gt -900) and $
      (atlasnodust_agn.ebv_hahb_err gt 0.0),nindx_agn)

    x_agn = atlasdust_agn[indx_agn].zstrong_12oh_oiiinii_niiha
    xerr_agn = atlasdust_agn[indx_agn].zstrong_12oh_oiiinii_niiha_err
    
    lir_agn = atlasdust_agn[indx_agn].iras_12*1D-23*light/12E4^2 ; [erg/s/cm2]
    lir_err_agn = atlasdust_agn[indx_agn].iras_12_err*1D-23*light/12E4^2

    indxnfgs_agn = where((nfgsdust_agn.zstrong_12oh_oiiinii_niiha gt -900) and (nfgsnodust_agn.iras_12 gt -900) and $
      (nfgsnodust_agn.ebv_hahb_err gt 0.0),nindxnfgs_agn)

    xnfgs_agn = nfgsnodust_agn[indxnfgs_agn].zstrong_12oh_oiiinii_niiha
    xerrnfgs_agn = nfgsnodust_agn[indxnfgs_agn].zstrong_12oh_oiiinii_niiha_err
    
    lirnfgs_agn = nfgsdust_agn[indxnfgs_agn].iras_12*1D-23*light/12E4^2 ; [erg/s/cm2]
    lirnfgs_err_agn = nfgsdust_agn[indxnfgs_agn].iras_12_err*1D-23*light/12E4^2

; plotting variables    
    
    xrange = ohrange3
    yrange = [1.7,5.3]

    xtitle = '12 + log (O/H)'
    ytitle = 'log [L(H\alpha)_{cor}/L(12\mu'+'m'+')]'

; SF galaxies    
    
    ha = atlasnodust[indx].h_alpha[0]
    ha_err = atlasnodust[indx].h_alpha[1]

    y = alog10(ha/lir)
    yerr = im_compute_error(ha,ha_err,lir,lir_err,/log)

    hanfgs = nfgsnodust[indxnfgs].h_alpha[0]
    hanfgs_err = nfgsnodust[indxnfgs].h_alpha[1]

    ynfgs = alog10(hanfgs/lirnfgs)
    yerrnfgs = im_compute_error(hanfgs,hanfgs_err,lirnfgs,lirnfgs_err,/log)

; AGN
    
    ha_agn = atlasnodust_agn[indx_agn].h_alpha[0]
    ha_err_agn = atlasnodust_agn[indx_agn].h_alpha[1]
    
    y_agn = alog10(ha_agn/lir_agn)
    yerr_agn = im_compute_error(ha_agn,ha_err_agn,lir_agn,lir_err_agn,/log)

    hanfgs_agn = nfgsnodust_agn[indxnfgs_agn].h_alpha[0]
    hanfgs_err_agn = nfgsnodust_agn[indxnfgs_agn].h_alpha[1]
    
    ynfgs_agn = alog10(hanfgs_agn/lirnfgs_agn)
    yerrnfgs_agn = im_compute_error(hanfgs_agn,hanfgs_err_agn,lirnfgs_agn,lirnfgs_err_agn,/log)

; combine the two samples - SF

    xbig = [x,xnfgs] & xerrbig = [xerr,xerrnfgs]
    ybig = [y,ynfgs] & yerrbig = [yerr,yerrnfgs]
    
; combine the two samples - AGN

    xbig_agn = [x_agn,xnfgs_agn] & xerrbig_agn = [xerr_agn,xerrnfgs_agn]
    ybig_agn = [y_agn,ynfgs_agn] & yerrbig_agn = [yerr_agn,yerrnfgs_agn]

    stats = im_stats(ybig,/verbose)
    xstr = strtrim(string(stats.median_rej,format='(F12.2)'),2)+' ('+$
      strtrim(string(stats.mean_rej,format='(F12.2)'),2)+$
      '\pm'+strtrim(string(stats.sigma_rej,format='(F12.2)'),2)+')'

    atlas1d_lineplot, xbig, ybig, xerrbig, yerrbig, plottype=3, postscript=postscript, $
      xtitle=xtitle, ytitle=ytitle, xrange=xrange, yrange=yrange, legendtype=0, $
      charsize=charsize_7, position=pos[*,0], yminor=3, atlaspsize=1.0
    atlas1d_lineplot, xbig_agn, ybig_agn, xerrbig_agn, yerrbig_agn, postscript=postscript, $
      atlaspsize=0.8, atlasfill=0, /overplot, atlascolor='red'

    im_openclose, postscript=postscript, /close

; ------------------------------------------------------------
; EW(OII) vs L(Ha)/L(12) - 2-panel
; ------------------------------------------------------------
    
    psname = 'ewoii_vs_lha_l12'
    im_openclose, pspath+psname, postscript=postscript, xsize=8.5, ysize=8.0, /encapsulated

    pagemaker, nx=1, ny=1, xspace=0, yspace=0.0, width=6.6, height=6.6, $
      xmargin=[1.5,0.4], ymargin=[0.3,1.1], xpage=8.5, ypage=8.0, $
      position=pos, /normal

; SF galaxies
    
    indx = where((atlasdust.oii_3727_ew[0]/atlasdust.oii_3727_ew[1] gt snrcut) and (atlasdust.iras_12 gt -900) and $
      (atlasnodust.ebv_hahb_err gt 0.0),nindx)

    x = atlasdust[indx].oii_3727_ew[0]
    xerr = atlasdust[indx].oii_3727_ew[1]
    
    lir = atlasdust[indx].iras_12*1D-23*light/12E4^2 ; [erg/s/cm2]
    lir_err = atlasdust[indx].iras_12_err*1D-23*light/12E4^2

    indxnfgs = where((nfgsdust.oii_3727_ew[0]/nfgsdust.oii_3727_ew[1] gt snrcut) and (nfgsdust.iras_12 gt -900) and $
      (nfgsnodust.ebv_hahb_err gt 0.0),nindxnfgs)

    xnfgs = nfgsdust[indxnfgs].oii_3727_ew[0]
    xerrnfgs = nfgsdust[indxnfgs].oii_3727_ew[1]
    
    lirnfgs = nfgsdust[indxnfgs].iras_12*1D-23*light/12E4^2 ; [erg/s/cm2]
    lirnfgs_err = nfgsdust[indxnfgs].iras_12_err*1D-23*light/12E4^2

; AGN
    
    indx_agn = where((atlasdust_agn.oii_3727_ew[0]/atlasdust_agn.oii_3727_ew[1] gt snrcut) and $
      (atlasdust_agn.iras_12 gt -900) and (atlasnodust_agn.ebv_hahb_err gt 0.0),nindx_agn)

    x_agn = atlasdust_agn[indx_agn].oii_3727_ew[0]
    xerr_agn = atlasdust_agn[indx_agn].oii_3727_ew[1]
    
    lir_agn = atlasdust_agn[indx_agn].iras_12*1D-23*light/12E4^2 ; [erg/s/cm2]
    lir_err_agn = atlasdust_agn[indx_agn].iras_12_err*1D-23*light/12E4^2

    indxnfgs_agn = where((nfgsdust_agn.oii_3727_ew[0]/nfgsdust_agn.oii_3727_ew[1] gt -900) and $
      (nfgsdust_agn.iras_12 gt -900) and (nfgsnodust_agn.ebv_hahb_err gt 0.0),nindxnfgs_agn)

    xnfgs_agn = nfgsdust_agn[indxnfgs_agn].oii_3727_ew[0]
    xerrnfgs_agn = nfgsdust_agn[indxnfgs_agn].oii_3727_ew[1]
    
    lirnfgs_agn = nfgsdust_agn[indxnfgs_agn].iras_12*1D-23*light/12E4^2 ; [erg/s/cm2]
    lirnfgs_err_agn = nfgsdust_agn[indxnfgs_agn].iras_12_err*1D-23*light/12E4^2

; plotting variables    
    
    xrange = ewoiirangenolog
    yrange = [1.7,5.3]

    xtitle = 'EW([O II] \lambda3727) ['+angstrom()+']'
    ytitle = 'log [L(H\alpha)_{cor}/L(12\mu'+'m'+')]'

; SF galaxies    
    
    ha = atlasnodust[indx].h_alpha[0]
    ha_err = atlasnodust[indx].h_alpha[1]

    y = alog10(ha/lir)
    yerr = im_compute_error(ha,ha_err,lir,lir_err,/log)

    hanfgs = nfgsnodust[indxnfgs].h_alpha[0]
    hanfgs_err = nfgsnodust[indxnfgs].h_alpha[1]

    ynfgs = alog10(hanfgs/lirnfgs)
    yerrnfgs = im_compute_error(hanfgs,hanfgs_err,lirnfgs,lirnfgs_err,/log)

; AGN
    
    ha_agn = atlasnodust_agn[indx_agn].h_alpha[0]
    ha_err_agn = atlasnodust_agn[indx_agn].h_alpha[1]
    
    y_agn = alog10(ha_agn/lir_agn)
    yerr_agn = im_compute_error(ha_agn,ha_err_agn,lir_agn,lir_err_agn,/log)

    hanfgs_agn = nfgsnodust_agn[indxnfgs_agn].h_alpha[0]
    hanfgs_err_agn = nfgsnodust_agn[indxnfgs_agn].h_alpha[1]
    
    ynfgs_agn = alog10(hanfgs_agn/lirnfgs_agn)
    yerrnfgs_agn = im_compute_error(hanfgs_agn,hanfgs_err_agn,lirnfgs_agn,lirnfgs_err_agn,/log)

; combine the two samples - SF

    xbig = [x,xnfgs] & xerrbig = [xerr,xerrnfgs]
    ybig = [y,ynfgs] & yerrbig = [yerr,yerrnfgs]
    
; combine the two samples - AGN

    xbig_agn = [x_agn,xnfgs_agn] & xerrbig_agn = [xerr_agn,xerrnfgs_agn]
    ybig_agn = [y_agn,ynfgs_agn] & yerrbig_agn = [yerr_agn,yerrnfgs_agn]

    stats = im_stats(ybig,/verbose)
    xstr = strtrim(string(stats.median_rej,format='(F12.2)'),2)+' ('+$
      strtrim(string(stats.mean_rej,format='(F12.2)'),2)+$
      '\pm'+strtrim(string(stats.sigma_rej,format='(F12.2)'),2)+')'

    atlas1d_lineplot, xbig, ybig, xerrbig, yerrbig, plottype=3, postscript=postscript, $
      xtitle=xtitle, ytitle=ytitle, xrange=xrange, yrange=yrange, legendtype=0, $
      charsize=charsize_7, position=pos[*,0], yminor=3, atlaspsize=1.0, /xlog
    atlas1d_lineplot, xbig_agn, ybig_agn, xerrbig_agn, yerrbig_agn, postscript=postscript, $
      atlaspsize=0.8, atlasfill=0, /overplot, atlascolor='red'

;   legend, textoidl(xstr), /right, /top, box=0, charsize=charsize_3, charthick=postthick
    
    im_openclose, postscript=postscript, /close

; ------------------------------------------------------------
; L(IR) vs L(Ha)/L(12) - 2-panel
; ------------------------------------------------------------
    
    psname = 'lir_vs_lha_l12_2panel'
    im_openclose, pspath+psname, postscript=postscript, /encapsulated, xsize=7.0, ysize=7.5

    pagemaker, nx=1, ny=2, height=[3.0,3.0], width=5.5, xmargin=[1.1,0.4], $
      ymargin=[0.4,1.1], xspace=0, yspace=0, xpage=7.0, ypage=7.5, $
      position=pos, /normal

; SF galaxies
    
    indx = where((atlasnodust.ir_flux gt -900) and (atlasnodust.iras_12 gt -900) and $
      (atlasnodust.ebv_hahb_err gt 0.0),nindx)

    x = atlasdust[indx].ir_lum
    xerr = atlasdust[indx].ir_lum_err
    
    lir = atlasdust[indx].iras_12*1D-23*light/12E4^2 ; [erg/s/cm2]
    lir_err = atlasdust[indx].iras_12_err*1D-23*light/12E4^2

    indxnfgs = where((nfgsdust.ir_flux gt -900) and (nfgsnodust.iras_12 gt -900) and $
      (nfgsnodust.ebv_hahb_err gt 0.0),nindxnfgs)

    xnfgs = nfgsnodust[indxnfgs].ir_lum
    xerrnfgs = nfgsnodust[indxnfgs].ir_lum_err
    
    lirnfgs = nfgsdust[indxnfgs].iras_12*1D-23*light/12E4^2 ; [erg/s/cm2]
    lirnfgs_err = nfgsdust[indxnfgs].iras_12_err*1D-23*light/12E4^2

; AGN
    
    indx_agn = where((atlasnodust_agn.ir_flux gt -900) and (atlasnodust_agn.iras_12 gt -900) and $
      (atlasnodust_agn.ebv_hahb_err gt 0.0),nindx_agn)

    x_agn = atlasdust_agn[indx_agn].ir_lum
    xerr_agn = atlasdust_agn[indx_agn].ir_lum_err
    
    lir_agn = atlasdust_agn[indx_agn].iras_12*1D-23*light/12E4^2 ; [erg/s/cm2]
    lir_err_agn = atlasdust_agn[indx_agn].iras_12_err*1D-23*light/12E4^2

    indxnfgs_agn = where((nfgsdust_agn.ir_flux gt -900) and (nfgsnodust_agn.iras_12 gt -900) and $
      (nfgsnodust_agn.ebv_hahb_err gt 0.0),nindxnfgs_agn)

    xnfgs_agn = nfgsnodust_agn[indxnfgs_agn].ir_lum
    xerrnfgs_agn = nfgsnodust_agn[indxnfgs_agn].ir_lum_err
    
    lirnfgs_agn = nfgsdust_agn[indxnfgs_agn].iras_12*1D-23*light/12E4^2 ; [erg/s/cm2]
    lirnfgs_err_agn = nfgsdust_agn[indxnfgs_agn].iras_12_err*1D-23*light/12E4^2

; plotting variables    
    
    xrange = lirrange2
    yrange = [1.7,5.3]

    xtitle = 'log L(IR) [L'+sunsymbol()+']'

; ##########################################################
; Panel 1: No reddening correction
; ##########################################################

; SF galaxies    
    
    ha = atlasdust[indx].h_alpha[0]
    ha_err = atlasdust[indx].h_alpha[1]
    
    y = alog10(ha/lir)
    yerr = im_compute_error(ha,ha_err,lir,lir_err,/log)

    hanfgs = nfgsdust[indxnfgs].h_alpha[0]
    hanfgs_err = nfgsdust[indxnfgs].h_alpha[1]
    
    ynfgs = alog10(hanfgs/lirnfgs)
    yerrnfgs = im_compute_error(hanfgs,hanfgs_err,lirnfgs,lirnfgs_err,/log)

; AGN
    
    ha_agn = atlasdust_agn[indx_agn].h_alpha[0]
    ha_err_agn = atlasdust_agn[indx_agn].h_alpha[1]
    
    y_agn = alog10(ha_agn/lir_agn)
    yerr_agn = im_compute_error(ha_agn,ha_err_agn,lir_agn,lir_err_agn,/log)

    hanfgs_agn = nfgsdust_agn[indxnfgs_agn].h_alpha[0]
    hanfgs_err_agn = nfgsdust_agn[indxnfgs_agn].h_alpha[1]
    
    ynfgs_agn = alog10(hanfgs_agn/lirnfgs_agn)
    yerrnfgs_agn = im_compute_error(hanfgs_agn,hanfgs_err_agn,lirnfgs_agn,lirnfgs_err_agn,/log)

; combine the two samples - SF

    xbig = [x,xnfgs] & xerrbig = [xerr,xerrnfgs]
    ybig = [y,ynfgs] & yerrbig = [yerr,yerrnfgs]
    
; combine the two samples - AGN

    xbig_agn = [x_agn,xnfgs_agn] & xerrbig_agn = [xerr_agn,xerrnfgs_agn]
    ybig_agn = [y_agn,ynfgs_agn] & yerrbig_agn = [yerr_agn,yerrnfgs_agn]

    stats = im_stats(ybig,/verbose)
    xstr = strtrim(string(stats.median_rej,format='(F12.2)'),2)+' ('+$
      strtrim(string(stats.mean_rej,format='(F12.2)'),2)+$
      '\pm'+strtrim(string(stats.sigma_rej,format='(F12.2)'),2)+')'

    ytitle = 'log [L(H\alpha)_{obs}/L(12\mu'+'m'+')]'

    atlas1d_lineplot, xbig, ybig, xerrbig, yerrbig, plottype=3, postscript=postscript, $
      xtitle='', ytitle=ytitle, xrange=xrange, yrange=yrange, legendtype=0, $
      charsize=charsize_7, position=pos[*,0], xtickname=replicate(' ',10), yminor=3, $
      atlaspsize=1.0
    atlas1d_lineplot, xbig_agn, ybig_agn, xerrbig_agn, yerrbig_agn, postscript=postscript, $
      atlaspsize=0.8, atlasfill=0, /overplot, atlascolor='red'

    legend, '(a)', /left, /top, box=0, charsize=charsize_3, charthick=postthick
    legend, textoidl(xstr), /right, /top, box=0, charsize=charsize_3, charthick=postthick

; ############################################################
; Panel 2: Good Reddening Correction
; ############################################################

; SF galaxies    
    
    ha = atlasnodust[indx].h_alpha[0]
    ha_err = atlasnodust[indx].h_alpha[1]

    y = alog10(ha/lir)
    yerr = im_compute_error(ha,ha_err,lir,lir_err,/log)

    hanfgs = nfgsnodust[indxnfgs].h_alpha[0]
    hanfgs_err = nfgsnodust[indxnfgs].h_alpha[1]

    ynfgs = alog10(hanfgs/lirnfgs)
    yerrnfgs = im_compute_error(hanfgs,hanfgs_err,lirnfgs,lirnfgs_err,/log)

; AGN
    
    ha_agn = atlasnodust_agn[indx_agn].h_alpha[0]
    ha_err_agn = atlasnodust_agn[indx_agn].h_alpha[1]
    
    y_agn = alog10(ha_agn/lir_agn)
    yerr_agn = im_compute_error(ha_agn,ha_err_agn,lir_agn,lir_err_agn,/log)

    hanfgs_agn = nfgsnodust_agn[indxnfgs_agn].h_alpha[0]
    hanfgs_err_agn = nfgsnodust_agn[indxnfgs_agn].h_alpha[1]
    
    ynfgs_agn = alog10(hanfgs_agn/lirnfgs_agn)
    yerrnfgs_agn = im_compute_error(hanfgs_agn,hanfgs_err_agn,lirnfgs_agn,lirnfgs_err_agn,/log)

; combine the two samples - SF

    xbig = [x,xnfgs] & xerrbig = [xerr,xerrnfgs]
    ybig = [y,ynfgs] & yerrbig = [yerr,yerrnfgs]
    
; combine the two samples - AGN

    xbig_agn = [x_agn,xnfgs_agn] & xerrbig_agn = [xerr_agn,xerrnfgs_agn]
    ybig_agn = [y_agn,ynfgs_agn] & yerrbig_agn = [yerr_agn,yerrnfgs_agn]

    stats = im_stats(ybig,/verbose)
    xstr = strtrim(string(stats.median_rej,format='(F12.2)'),2)+' ('+$
      strtrim(string(stats.mean_rej,format='(F12.2)'),2)+$
      '\pm'+strtrim(string(stats.sigma_rej,format='(F12.2)'),2)+')'

    ytitle = 'log [L(H\alpha)_{cor}/L(12\mu'+'m'+')]'

    atlas1d_lineplot, xbig, ybig, xerrbig, yerrbig, plottype=3, postscript=postscript, $
      xtitle=xtitle, ytitle=ytitle, xrange=xrange, yrange=yrange, legendtype=0, $
      charsize=charsize_7, position=pos[*,1], yminor=3, atlaspsize=1.0, /noerase
    atlas1d_lineplot, xbig_agn, ybig_agn, xerrbig_agn, yerrbig_agn, postscript=postscript, $
      atlaspsize=0.8, atlasfill=0, /overplot, atlascolor='red'

    legend, '(b)', /left, /top, box=0, charsize=charsize_3, charthick=postthick
;   legend, textoidl('Extinction-Corrected H\alpha'), /left, /bottom, box=0, charsize=charsize_3, charthick=postthick
    legend, textoidl(xstr), /right, /top, box=0, charsize=charsize_3, charthick=postthick
    
    im_openclose, postscript=postscript, /close

; ------------------------------------------------------------
; L(IR) vs L(Ha)/L(25) - 2-panel
; ------------------------------------------------------------
    
    psname = 'lir_vs_lha_l25_2panel'
    im_openclose, pspath+psname, postscript=postscript, /encapsulated, xsize=7.0, ysize=7.5

    pagemaker, nx=1, ny=2, height=[3.0,3.0], width=5.5, xmargin=[1.1,0.4], $
      ymargin=[0.4,1.1], xspace=0, yspace=0, xpage=7.0, ypage=7.5, $
      position=pos, /normal

; SF galaxies
    
    indx = where((atlasnodust.ir_flux gt -900) and (atlasnodust.iras_25 gt -900) and $
      (atlasnodust.ebv_hahb_err gt 0.0),nindx)

    x = atlasdust[indx].ir_lum
    xerr = atlasdust[indx].ir_lum_err
    
    lir = atlasdust[indx].iras_25*1D-23*light/25E4^2 ; [erg/s/cm2]
    lir_err = atlasdust[indx].iras_25_err*1D-23*light/25E4^2

    indxnfgs = where((nfgsdust.ir_flux gt -900) and (nfgsnodust.iras_25 gt -900) and $
      (nfgsnodust.ebv_hahb_err gt 0.0),nindxnfgs)

    xnfgs = nfgsnodust[indxnfgs].ir_lum
    xerrnfgs = nfgsnodust[indxnfgs].ir_lum_err
    
    lirnfgs = nfgsdust[indxnfgs].iras_25*1D-23*light/25E4^2 ; [erg/s/cm2]
    lirnfgs_err = nfgsdust[indxnfgs].iras_25_err*1D-23*light/25E4^2

; AGN
    
    indx_agn = where((atlasnodust_agn.ir_flux gt -900) and (atlasnodust_agn.iras_25 gt -900) and $
      (atlasnodust_agn.ebv_hahb_err gt 0.0),nindx_agn)

    x_agn = atlasdust_agn[indx_agn].ir_lum
    xerr_agn = atlasdust_agn[indx_agn].ir_lum_err
    
    lir_agn = atlasdust_agn[indx_agn].iras_25*1D-23*light/25E4^2 ; [erg/s/cm2]
    lir_err_agn = atlasdust_agn[indx_agn].iras_25_err*1D-23*light/25E4^2

    indxnfgs_agn = where((nfgsdust_agn.ir_flux gt -900) and (nfgsnodust_agn.iras_25 gt -900) and $
      (nfgsnodust_agn.ebv_hahb_err gt 0.0),nindxnfgs_agn)

    xnfgs_agn = nfgsnodust_agn[indxnfgs_agn].ir_lum
    xerrnfgs_agn = nfgsnodust_agn[indxnfgs_agn].ir_lum_err
    
    lirnfgs_agn = nfgsdust_agn[indxnfgs_agn].iras_25*1D-23*light/25E4^2 ; [erg/s/cm2]
    lirnfgs_err_agn = nfgsdust_agn[indxnfgs_agn].iras_25_err*1D-23*light/25E4^2

; plotting variables    
    
    xrange = lirrange2
    yrange = [0.5,6.7]

    xtitle = 'log L(IR) [L'+sunsymbol()+']'

; ##########################################################
; Panel 1: No reddening correction
; ##########################################################

; SF galaxies    
    
    ha = atlasdust[indx].h_alpha[0]
    ha_err = atlasdust[indx].h_alpha[1]
    
    y = alog10(ha/lir)
    yerr = im_compute_error(ha,ha_err,lir,lir_err,/log)

    hanfgs = nfgsdust[indxnfgs].h_alpha[0]
    hanfgs_err = nfgsdust[indxnfgs].h_alpha[1]
    
    ynfgs = alog10(hanfgs/lirnfgs)
    yerrnfgs = im_compute_error(hanfgs,hanfgs_err,lirnfgs,lirnfgs_err,/log)

; AGN
    
    ha_agn = atlasdust_agn[indx_agn].h_alpha[0]
    ha_err_agn = atlasdust_agn[indx_agn].h_alpha[1]
    
    y_agn = alog10(ha_agn/lir_agn)
    yerr_agn = im_compute_error(ha_agn,ha_err_agn,lir_agn,lir_err_agn,/log)

    hanfgs_agn = nfgsdust_agn[indxnfgs_agn].h_alpha[0]
    hanfgs_err_agn = nfgsdust_agn[indxnfgs_agn].h_alpha[1]
    
    ynfgs_agn = alog10(hanfgs_agn/lirnfgs_agn)
    yerrnfgs_agn = im_compute_error(hanfgs_agn,hanfgs_err_agn,lirnfgs_agn,lirnfgs_err_agn,/log)

; combine the two samples - SF

    xbig = [x,xnfgs] & xerrbig = [xerr,xerrnfgs]
    ybig = [y,ynfgs] & yerrbig = [yerr,yerrnfgs]
    
; combine the two samples - AGN

    xbig_agn = [x_agn,xnfgs_agn] & xerrbig_agn = [xerr_agn,xerrnfgs_agn]
    ybig_agn = [y_agn,ynfgs_agn] & yerrbig_agn = [yerr_agn,yerrnfgs_agn]

    stats = im_stats(ybig,/verbose)
    xstr = strtrim(string(stats.median_rej,format='(F12.2)'),2)+' ('+$
      strtrim(string(stats.mean_rej,format='(F12.2)'),2)+$
      '\pm'+strtrim(string(stats.sigma_rej,format='(F12.2)'),2)+')'

    ytitle = 'log [L(H\alpha)_{obs}/L(25\mu'+'m'+')]'

    atlas1d_lineplot, xbig, ybig, xerrbig, yerrbig, plottype=3, postscript=postscript, $
      xtitle='', ytitle=ytitle, xrange=xrange, yrange=yrange, legendtype=0, $
      charsize=charsize_7, position=pos[*,0], xtickname=replicate(' ',10), yminor=3, $
      atlaspsize=1.0
    atlas1d_lineplot, xbig_agn, ybig_agn, xerrbig_agn, yerrbig_agn, postscript=postscript, $
      atlaspsize=0.8, atlasfill=0, /overplot, atlascolor='red'

    legend, '(a)', /left, /top, box=0, charsize=charsize_3, charthick=postthick
    legend, textoidl('Observed H\alpha'), /left, /bottom, box=0, charsize=charsize_3, charthick=postthick
    legend, textoidl(xstr), /right, /top, box=0, charsize=charsize_3, charthick=postthick

; ############################################################
; Panel 2: Good Reddening Correction
; ############################################################

; SF galaxies    
    
    ha = atlasnodust[indx].h_alpha[0]
    ha_err = atlasnodust[indx].h_alpha[1]

    y = alog10(ha/lir)
    yerr = im_compute_error(ha,ha_err,lir,lir_err,/log)

    hanfgs = nfgsnodust[indxnfgs].h_alpha[0]
    hanfgs_err = nfgsnodust[indxnfgs].h_alpha[1]

    ynfgs = alog10(hanfgs/lirnfgs)
    yerrnfgs = im_compute_error(hanfgs,hanfgs_err,lirnfgs,lirnfgs_err,/log)

; AGN
    
    ha_agn = atlasnodust_agn[indx_agn].h_alpha[0]
    ha_err_agn = atlasnodust_agn[indx_agn].h_alpha[1]
    
    y_agn = alog10(ha_agn/lir_agn)
    yerr_agn = im_compute_error(ha_agn,ha_err_agn,lir_agn,lir_err_agn,/log)

    hanfgs_agn = nfgsnodust_agn[indxnfgs_agn].h_alpha[0]
    hanfgs_err_agn = nfgsnodust_agn[indxnfgs_agn].h_alpha[1]
    
    ynfgs_agn = alog10(hanfgs_agn/lirnfgs_agn)
    yerrnfgs_agn = im_compute_error(hanfgs_agn,hanfgs_err_agn,lirnfgs_agn,lirnfgs_err_agn,/log)

; combine the two samples - SF

    xbig = [x,xnfgs] & xerrbig = [xerr,xerrnfgs]
    ybig = [y,ynfgs] & yerrbig = [yerr,yerrnfgs]
    
; combine the two samples - AGN

    xbig_agn = [x_agn,xnfgs_agn] & xerrbig_agn = [xerr_agn,xerrnfgs_agn]
    ybig_agn = [y_agn,ynfgs_agn] & yerrbig_agn = [yerr_agn,yerrnfgs_agn]

    stats = im_stats(ybig,/verbose)
    xstr = strtrim(string(stats.median_rej,format='(F12.2)'),2)+' ('+$
      strtrim(string(stats.mean_rej,format='(F12.2)'),2)+$
      '\pm'+strtrim(string(stats.sigma_rej,format='(F12.2)'),2)+')'

    ytitle = 'log [L(H\alpha)_{cor}/L(25\mu'+'m'+')]'

    atlas1d_lineplot, xbig, ybig, xerrbig, yerrbig, plottype=3, postscript=postscript, $
      xtitle=xtitle, ytitle=ytitle, xrange=xrange, yrange=yrange, legendtype=0, $
      charsize=charsize_7, position=pos[*,1], yminor=3, atlaspsize=1.0, /noerase
    atlas1d_lineplot, xbig_agn, ybig_agn, xerrbig_agn, yerrbig_agn, postscript=postscript, $
      atlaspsize=0.8, atlasfill=0, /overplot, atlascolor='red'

    legend, '(b)', /left, /top, box=0, charsize=charsize_3, charthick=postthick
;   legend, textoidl('Extinction-Corrected H\alpha'), /left, /bottom, box=0, charsize=charsize_3, charthick=postthick
    legend, textoidl(xstr), /right, /top, box=0, charsize=charsize_3, charthick=postthick
    
    im_openclose, postscript=postscript, /close

; ------------------------------------------------------------
; L(IR) vs L(Ha)/L(IR) - 2-panel
; ------------------------------------------------------------
    
    psname = 'lir_vs_lhalir_2panel'
    im_openclose, pspath+psname, postscript=postscript, /encapsulated, xsize=7.0, ysize=7.5

    pagemaker, nx=1, ny=2, height=[3.0,3.0], width=5.5, xmargin=[1.1,0.4], $
      ymargin=[0.4,1.1], xspace=0, yspace=0, xpage=7.0, ypage=7.5, $
      position=pos, /normal

; SF galaxies
    
    indx = where((atlasnodust.ir_flux gt -900) and (atlasnodust.ebv_hahb_err gt 0.0),nindx)

    x = atlasdust[indx].ir_lum
    xerr = atlasdust[indx].ir_lum_err
    
    lir = atlasdust[indx].ir_flux ; [erg/s/cm2]
    lir_err = atlasdust[indx].ir_flux_err

    indxnfgs = where((nfgsdust.ir_flux gt -900) and (nfgsnodust.ebv_hahb_err gt 0.0),nindxnfgs)

    xnfgs = nfgsnodust[indxnfgs].ir_lum
    xerrnfgs = nfgsnodust[indxnfgs].ir_lum_err
    
    lirnfgs = nfgsdust[indxnfgs].ir_flux ; [erg/s/cm2]
    lirnfgs_err = nfgsdust[indxnfgs].ir_flux_err

; AGN
    
    indx_agn = where((atlasnodust_agn.ir_flux gt -900) and (atlasnodust_agn.ebv_hahb_err gt 0.0),nindx_agn)

    x_agn = atlasdust_agn[indx_agn].ir_lum
    xerr_agn = atlasdust_agn[indx_agn].ir_lum_err
    
    lir_agn = atlasdust_agn[indx_agn].ir_flux ; [erg/s/cm2]
    lir_err_agn = atlasdust_agn[indx_agn].ir_flux_err

    indxnfgs_agn = where((nfgsdust_agn.ir_flux gt -900) and (nfgsnodust_agn.ebv_hahb_err gt 0.0),nindxnfgs_agn)

    xnfgs_agn = nfgsnodust_agn[indxnfgs_agn].ir_lum
    xerrnfgs_agn = nfgsnodust_agn[indxnfgs_agn].ir_lum_err
    
    lirnfgs_agn = nfgsdust_agn[indxnfgs_agn].ir_flux ; [erg/s/cm2]
    lirnfgs_err_agn = nfgsdust_agn[indxnfgs_agn].ir_flux_err

; plotting variables    
    
    xrange = lirrange2
    yrange = [-5.5,0.5]

    xtitle = 'log L(IR) [L'+sunsymbol()+']'
    ytitle = 'log [L(H\alpha)/L(IR)]'

; ##########################################################
; Panel 1: No reddening correction
; ##########################################################

; SF galaxies    
    
    ha = atlasdust[indx].h_alpha[0]
    ha_err = atlasdust[indx].h_alpha[1]
    
    y = alog10(ha/lir)
    yerr = im_compute_error(ha,ha_err,lir,lir_err,/log)

    hanfgs = nfgsdust[indxnfgs].h_alpha[0]
    hanfgs_err = nfgsdust[indxnfgs].h_alpha[1]
    
    ynfgs = alog10(hanfgs/lirnfgs)
    yerrnfgs = im_compute_error(hanfgs,hanfgs_err,lirnfgs,lirnfgs_err,/log)

; AGN
    
    ha_agn = atlasdust_agn[indx_agn].h_alpha[0]
    ha_err_agn = atlasdust_agn[indx_agn].h_alpha[1]
    
    y_agn = alog10(ha_agn/lir_agn)
    yerr_agn = im_compute_error(ha_agn,ha_err_agn,lir_agn,lir_err_agn,/log)

    hanfgs_agn = nfgsdust_agn[indxnfgs_agn].h_alpha[0]
    hanfgs_err_agn = nfgsdust_agn[indxnfgs_agn].h_alpha[1]
    
    ynfgs_agn = alog10(hanfgs_agn/lirnfgs_agn)
    yerrnfgs_agn = im_compute_error(hanfgs_agn,hanfgs_err_agn,lirnfgs_agn,lirnfgs_err_agn,/log)

; combine the two samples - SF

    xbig = [x,xnfgs] & xerrbig = [xerr,xerrnfgs]
    ybig = [y,ynfgs] & yerrbig = [yerr,yerrnfgs]
    
; combine the two samples - AGN

    xbig_agn = [x_agn,xnfgs_agn] & xerrbig_agn = [xerr_agn,xerrnfgs_agn]
    ybig_agn = [y_agn,ynfgs_agn] & yerrbig_agn = [yerr_agn,yerrnfgs_agn]

    stats = im_stats(ybig,/verbose)
    xstr = strtrim(string(stats.median_rej,format='(F12.2)'),2)+' ('+$
      strtrim(string(stats.mean_rej,format='(F12.2)'),2)+$
      '\pm'+strtrim(string(stats.sigma_rej,format='(F12.2)'),2)+')'

    ytitle = 'log [L(H\alpha)_{obs}/L(IR)]'

    atlas1d_lineplot, xbig, ybig, xerrbig, yerrbig, plottype=3, postscript=postscript, $
      xtitle='', ytitle=ytitle, xrange=xrange, yrange=yrange, legendtype=0, $
      charsize=charsize_7, position=pos[*,0], xtickname=replicate(' ',10), yminor=3, $
      atlaspsize=1.0
    atlas1d_lineplot, xbig_agn, ybig_agn, xerrbig_agn, yerrbig_agn, postscript=postscript, $
      atlaspsize=0.8, atlasfill=0, /overplot, atlascolor='red'

    legend, '(a)', /left, /top, box=0, charsize=charsize_3, charthick=postthick
    legend, textoidl(xstr), /right, /top, box=0, charsize=charsize_3, charthick=postthick

; ############################################################
; Panel 2: Good Reddening Correction
; ############################################################

; SF galaxies    
    
    ha = atlasnodust[indx].h_alpha[0]
    ha_err = atlasnodust[indx].h_alpha[1]

    y = alog10(ha/lir)
    yerr = im_compute_error(ha,ha_err,lir,lir_err,/log)

    hanfgs = nfgsnodust[indxnfgs].h_alpha[0]
    hanfgs_err = nfgsnodust[indxnfgs].h_alpha[1]

    ynfgs = alog10(hanfgs/lirnfgs)
    yerrnfgs = im_compute_error(hanfgs,hanfgs_err,lirnfgs,lirnfgs_err,/log)

; AGN
    
    ha_agn = atlasnodust_agn[indx_agn].h_alpha[0]
    ha_err_agn = atlasnodust_agn[indx_agn].h_alpha[1]
    
    y_agn = alog10(ha_agn/lir_agn)
    yerr_agn = im_compute_error(ha_agn,ha_err_agn,lir_agn,lir_err_agn,/log)

    hanfgs_agn = nfgsnodust_agn[indxnfgs_agn].h_alpha[0]
    hanfgs_err_agn = nfgsnodust_agn[indxnfgs_agn].h_alpha[1]
    
    ynfgs_agn = alog10(hanfgs_agn/lirnfgs_agn)
    yerrnfgs_agn = im_compute_error(hanfgs_agn,hanfgs_err_agn,lirnfgs_agn,lirnfgs_err_agn,/log)

; combine the two samples - SF

    xbig = [x,xnfgs] & xerrbig = [xerr,xerrnfgs]
    ybig = [y,ynfgs] & yerrbig = [yerr,yerrnfgs]
    
; combine the two samples - AGN

    xbig_agn = [x_agn,xnfgs_agn] & xerrbig_agn = [xerr_agn,xerrnfgs_agn]
    ybig_agn = [y_agn,ynfgs_agn] & yerrbig_agn = [yerr_agn,yerrnfgs_agn]

    stats = im_stats(ybig,/verbose)
    xstr = strtrim(string(stats.median_rej,format='(F12.2)'),2)+' ('+$
      strtrim(string(stats.mean_rej,format='(F12.2)'),2)+$
      '\pm'+strtrim(string(stats.sigma_rej,format='(F12.2)'),2)+')'

    ytitle = 'log [L(H\alpha)_{cor}/L(IR)]'

    atlas1d_lineplot, xbig, ybig, xerrbig, yerrbig, plottype=3, postscript=postscript, $
      xtitle=xtitle, ytitle=ytitle, xrange=xrange, yrange=yrange, legendtype=0, $
      charsize=charsize_7, position=pos[*,1], yminor=3, atlaspsize=1.0, /noerase
    atlas1d_lineplot, xbig_agn, ybig_agn, xerrbig_agn, yerrbig_agn, postscript=postscript, $
      atlaspsize=0.8, atlasfill=0, /overplot, atlascolor='red'

    legend, '(b)', /left, /top, box=0, charsize=charsize_3, charthick=postthick
;   legend, textoidl('Extinction-Corrected H\alpha'), /left, /bottom, box=0, charsize=charsize_3, charthick=postthick
    legend, textoidl(xstr), /right, /top, box=0, charsize=charsize_3, charthick=postthick
    
    im_openclose, postscript=postscript, /close

; --------------------------------------------------    
; GENERATE THE WEB PAGE
; --------------------------------------------------    
    
    if keyword_set(postscript) then begin

       im_ps2html, htmlbase, html_path=html_path, cleanpng=0, npscols=3, _extra=extra
       
    endif

stop

; ------------------------------------------------------------
; L(B) vs L(IR)/[SFR(Ha)_obs+SFR(IR)]
; ------------------------------------------------------------
    
    psname = 'lb_vs_lir_sfr_hair'
    im_openclose, pspath+psname, postscript=postscript, xsize=8.5, ysize=8.0, /encapsulated

    pagemaker, nx=1, ny=1, xspace=0, yspace=0.0, width=6.6, height=6.6, $
      xmargin=[1.5,0.4], ymargin=[0.3,1.1], xpage=8.5, ypage=8.0, $
      position=pos, /normal

; SF galaxies
    
    indx = where((atlasdust.b_lum_obs gt -900) and (atlasdust.ir_flux gt -900) and $
      (atlasnodust.ebv_hahb_err gt 0.0),nindx)

    x = atlasdust[indx].b_lum_obs
    xerr = atlasdust[indx].b_lum_obs_err
    
    lir = atlasdust[indx].ir_flux ; [erg/s/cm2]
    lir_err = sqrt(atlasdust[indx].ir_flux_err^2 + (atlasdust[indx].ir_flux*0.15)^2)

    indxnfgs = where((nfgsdust.b_lum_obs gt -900) and (nfgsdust.ir_flux gt -900) and $
      (nfgsnodust.ebv_hahb_err gt 0.0),nindxnfgs)

    xnfgs = nfgsnodust[indxnfgs].b_lum_obs
    xerrnfgs = nfgsnodust[indxnfgs].b_lum_obs_err
    
    lirnfgs = nfgsdust[indxnfgs].ir_flux
    lirnfgs_err = sqrt(nfgsdust[indxnfgs].ir_flux_err^2 + (nfgsdust[indxnfgs].ir_flux*0.15)^2)
    
; AGN
    
    indx_agn = where((atlasdust_agn.b_lum_obs gt -900) and (atlasdust_agn.ir_flux gt -900) and $
      (atlasnodust_agn.ebv_hahb_err gt 0.0),nindx_agn)

    x_agn = atlasdust_agn[indx_agn].b_lum_obs
    xerr_agn = atlasdust_agn[indx_agn].b_lum_obs_err
    
    lir_agn = atlasdust_agn[indx_agn].ir_flux
    lir_err_agn = atlasdust_agn[indx_agn].ir_flux_err

    indxnfgs_agn = where((nfgsdust_agn.b_lum_obs gt -900) and (nfgsdust_agn.ir_flux gt -900) and $
      (nfgsnodust_agn.ebv_hahb_err gt 0.0),nindxnfgs_agn)

    xnfgs_agn = nfgsdust_agn[indxnfgs_agn].b_lum_obs
    xerrnfgs_agn = nfgsdust_agn[indxnfgs_agn].b_lum_obs_err
    
    lirnfgs_agn = nfgsdust_agn[indxnfgs_agn].ir_flux
    lirnfgs_err_agn = nfgsdust_agn[indxnfgs_agn].ir_flux_err

; plotting variables    
    
    xrange = lbrange2
    yrange = [-1.1,1.3]

    xtitle = 'log L(B) [L(B)'+sunsymbol()+']'
    ytitle = 'log {10^{-43} L(IR) / [\psi(H\alpha)_{obs}+\psi(IR)]}'

; SF galaxies    

    y1 = haconst*atlasdust[indx].h_alpha[0]+irconst*lir
    y1_err = sqrt((haconst*atlasdust[indx].h_alpha[1])^2+(irconst*lir_err)^2)
    
    y = alog10(lir/y1) - 43.0
    yerr = im_compute_error(lir,lir_err,y1,y1_err,/log)

    y1nfgs = haconst*nfgsdust[indxnfgs].h_alpha[0]+irconst*lirnfgs
    y1nfgs_err = sqrt((haconst*nfgsdust[indxnfgs].h_alpha[1])^2+(irconst*lirnfgs_err)^2)
    
    ynfgs = alog10(lirnfgs/y1nfgs) - 43.0
    yerrnfgs = im_compute_error(lirnfgs,lirnfgs_err,y1nfgs,y1nfgs_err,/log)

; AGN
    
    y1_agn = haconst*atlasdust_agn[indx_agn].h_alpha[0]+irconst*lir_agn
    y1_err_agn = sqrt((haconst*atlasdust_agn[indx_agn].h_alpha[1])^2+(irconst*lir_err_agn)^2)
    
    y_agn = alog10(lir_agn/y1_agn) - 43.0
    yerr_agn = im_compute_error(lir_agn,lir_err_agn,y1_agn,y1_err_agn,/log)

    y1nfgs_agn = haconst*nfgsdust_agn[indxnfgs_agn].h_alpha[0]+irconst*lirnfgs_agn - 43.0
    y1nfgs_err_agn = sqrt((haconst*nfgsdust_agn[indxnfgs_agn].h_alpha[1])^2+(irconst*lirnfgs_err_agn)^2)
    
    ynfgs_agn = alog10(lirnfgs_agn/y1nfgs_agn)
    yerrnfgs_agn = im_compute_error(lirnfgs_agn,lirnfgs_err_agn,y1nfgs_agn,y1nfgs_err_agn,/log)

; combine the two samples - SF

    xbig = [x,xnfgs] & xerrbig = [xerr,xerrnfgs]
    ybig = [y,ynfgs] & yerrbig = [yerr,yerrnfgs]
    
; combine the two samples - AGN

    xbig_agn = [x_agn,xnfgs_agn] & xerrbig_agn = [xerr_agn,xerrnfgs_agn]
    ybig_agn = [y_agn,ynfgs_agn] & yerrbig_agn = [yerr_agn,yerrnfgs_agn]

    stats = im_stats(ybig,/verbose)
    xstr = strtrim(string(stats.median_rej,format='(F12.2)'),2)+' ('+$
      strtrim(string(stats.mean_rej,format='(F12.2)'),2)+$
      '\pm'+strtrim(string(stats.sigma_rej,format='(F12.2)'),2)+')'

    atlas1d_lineplot, xbig, ybig, xerrbig, yerrbig, plottype=3, postscript=postscript, $
      xtitle=xtitle, ytitle=ytitle, xrange=xrange, yrange=yrange, legendtype=0, $
      charsize=charsize_7, position=pos[*,0], yminor=3, atlaspsize=1.0
    atlas1d_lineplot, xbig_agn, ybig_agn, xerrbig_agn, yerrbig_agn, postscript=postscript, $
      atlaspsize=0.8, atlasfill=0, /overplot, atlascolor='red'

    legend, textoidl(xstr), /right, /top, box=0, charsize=charsize_4, charthick=postthick
    
    im_openclose, postscript=postscript, /close

return
end    
