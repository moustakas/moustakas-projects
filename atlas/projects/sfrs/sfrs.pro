;+
; NAME:
;       SFRS
;
; PURPOSE:
;       Investigate optical star-formation rate indicators. 
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
;       J. Moustakas, 2003-2005, U of A
;-

pro sfrs, atlasdust, atlasnodust, nfgsdust, nfgsnodust, atlasdust_agn, $
  atlasnodust_agn, nfgsdust_agn, nfgsnodust_agn, sdssdust, sdssnodust, $
  sdssancillary, hii, postscript=postscript, paper=paper, cleanpng=cleanpng, $
  blackwhite=blackwhite, encapsulated=encapsulated, _extra=extra

; sfrs,atlasdust,atlasnodust,nfgsdust,nfgsnodust,atlasdust_agn,atlasnodust_agn,nfgsdust_agn,nfgsnodust_agn,sdssdust,sdssnodust,sdssancillary,hii

; define the cosmology and some constants

    loadct, 0, /silent
    
    red, h100=0.7, omega_0=0.3, omega_lambda=0.7
    
    lsun = 3.826D33              ; [erg/s]
    light = 2.99792458D18        ; speed of light [A/s]
    haconst = 7.9D-42            ; K98 conversion L(Ha) --> SFR(Ha)
    irconst = 4.5D-44            ; K98 conversion L(IR) --> SFR(IR)
    loghaconst = alog10(haconst)
    logUconst = alog10(7.9D-30)  ; Cram et al. (1998) conversion L(U) --> SFR(U)
    hasfrconstoffset = 41.0D
    Uhasfrconstoffset = 42.0 ; 43.0D
    mbolsun = 4.74         ; bolometric absolute solar magnitude [mag]
  
    LBnorm = 1E10   ; [L_B_sun]
    sfrnorm = 1.0   ; [M_sun/yr]
    elumnorm = 1D41 ; [erg/s]
    Ulumnorm = 1D42 ; [erg/s]

    syserr = 0.0                ; systematic flux error [%]
    nfgssyserr = 0.0
    ewcut = 1.0
    ewcut2 = 1.0

    Uminerr = alog10(1+0.15) ; [dex] minimum synthesized U-band uncertainty [15%]
    Bminerr = 0.05           ; minimum SDSS B-band uncertainty 
    
    Uinfo = im_filterspecs(filterlist='bessell_U.par')
    Uconstant = Uinfo.weff*Uinfo.vega_flam

    lusun = lsun*10^(-0.4*(uinfo.solarmags-mbolsun))
    
    HaHb = 2.86
;   HaHb = return_tbalmer(/HaHb)

    kl_ha = k_lambda(6563.0,/odonnell)
    kl_hb = k_lambda(4861.0,/odonnell)
    kl_oii = k_lambda(3727.0,/odonnell)
    kl_oiii = k_lambda(5007.0,/odonnell)
    kl_u = k_lambda(Uinfo.weff,/odonnell)
    
    bigbprime_ha = 0.4*(-2.5*kl_ha/(kl_hb-kl_ha))
    bigbprime_oii = 0.4*(-2.5*kl_oii/(kl_hb-kl_ha))
    bigbprime_oiii = 0.4*(-2.5*kl_oiii/(kl_hb-kl_ha))
    bigbprime_u = 0.4*(-2.5*kl_u/(kl_hb-kl_ha))
    
; define path names    
    
    htmlbase = 'sfrs'

    html_path = atlas_path(/web)+'analysis/'
    pspath = atlas_path(/web)+'analysis/'+htmlbase+'/'
    paperpath = atlas_path(/papers)+'sfrs/FIG_SFRS/'
    latexpath = atlas_path(/papers)+'sfrs/'
    sfrspath = atlas_path(/projects)+'sfrs/'

    if keyword_set(paper) then postscript = 1L
    if keyword_set(postscript) then begin
       postthick = 8.0 
       postthick2 = 12.0
       postthick3 = 5.0
       symthick = 5.0
       highz_errthick = 4.0
    endif else begin
       blackwhite = 0L
       postthick = 2.0
       postthick2 = 2.0
       postthick3 = 2.0
       symthick = 2.0
       highz_errthick = 1.0
       im_window, 0, /square
    endelse

    if keyword_set(blackwhite) then begin
       pspath = pspath+'blackwhite/'
       suffix = 'blackwhite'
       color = 0L
    endif else begin
       suffix = ''
       color = 0L
    endelse
    
    if (n_elements(snrcut) eq 0L) then snrcut = 3.0
    snrcut_highz = 0.0

; read the data and the models    
    
    if (n_elements(atlasdust) eq 0L) then atlasdust = read_atlas_sfrs_sample(atlasnodust=atlasnodust)
    if (n_elements(atlasdust_agn) eq 0L) then atlasdust_agn = read_atlas_sfrs_sample(atlasnodust=atlasnodust_agn,/agn)
    if (n_elements(nfgsdust) eq 0L) then nfgsdust = read_nfgs_sfrs_sample(nfgsnodust=nfgsnodust)
    if (n_elements(nfgsdust_agn) eq 0L) then nfgsdust_agn = read_nfgs_sfrs_sample(nfgsnodust=nfgsnodust_agn,/agn)
    if (n_elements(sdssdust) eq 0L) then sdssdust = read_sdss_sfrs_sample(sdssnodust=sdssnodust,sdssancillary=sdssancillary)
    if (n_elements(hii) eq 0L) then hii = read_hii_regions(/samplerefs)

    sbgrids = read_kewley_grids(model=3,Z=Z,U=U)  ; model=1

; read the SFR calibrations    
    
    atlas_sfrs = mrdfits(atlas_path(/projects)+'sfrs/atlas_sfrs.fits',1,/silent)
    atlas_sfrs_agn = mrdfits(atlas_path(/projects)+'sfrs/atlas_sfrs_agn.fits',1,/silent)
    nfgs_sfrs = mrdfits(atlas_path(/projects)+'sfrs/nfgs_sfrs.fits',1,/silent)
    nfgs_sfrs_agn = mrdfits(atlas_path(/projects)+'sfrs/nfgs_sfrs_agn.fits',1,/silent)

;   sfrs = nfgs_sfrs
;   sfrs_agn = nfgs_sfrs_agn
    sfrs = struct_append(atlas_sfrs,nfgs_sfrs)
    sfrs_agn = struct_append(atlas_sfrs_agn,nfgs_sfrs_agn)

    hgsnrcut = 7.0
    hbewcut = 10.0
    
    hgsnrcut_sfrs = atlas_sfrs[0].hgsnrcut
    hbewcut_sfrs = atlas_sfrs[0].hbewcut
    stats = mrdfits(atlas_path(/projects)+'sfrs/stats_sfrs.fits',1,/silent)

; read the high-redshift data

    g99 = read_99glazebrook()
    h02 = read_02hicks()
    t02 = read_02tresse()
    lilly03 = read_03lilly()
    liang04 = read_04liang(datanodust=liang04nodust)
    m05 = read_05maier(datanodust=m05nodust)
    shap05 = read_05shapley(datanodust=shap05nodust)
    sava05 = read_05savaglio(datanodust=sava05nodust)

    k03 = read_03kobulnicky()
    k04 = read_04kobulnicky()
    k04 = struct_append(k03,k04)

    if keyword_set(blackwhite) then localcolor = 'medium gray' else localcolor = 'medium gray'
    localsym = 108
    localpsize = 1.0
    
    if keyword_set(blackwhite) then g99color = 'dark gray' else g99color = 'navy'
    g99sym = 122
    g99psize = 1.5

    if keyword_set(blackwhite) then t02color = 'dark gray' else t02color = 'dark orchid'
    t02sym = 104
    t02psize = 1.3

    if keyword_set(blackwhite) then h02color = 'dark gray' else h02color = 'dark green'
    h02sym = 106
    h02psize = 1.1

    if keyword_set(blackwhite) then liang04color = 'dark gray' else liang04color = 'navy'
    liang04sym = 106
    liang04psize = 1.0

    if keyword_set(blackwhite) then lilly03color = 'dark gray' else lilly03color = 'dodger blue'
    lilly03sym = 122
    lilly03psize = 1.3

    if keyword_set(blackwhite) then k04color = 'dark gray' else k04color = 'black'
    k04sym = '104'
    k04psize = 0.8

    if keyword_set(blackwhite) then m05color = 'dark gray' else m05color = 'red'
    m05sym = 105
    m05psize = 1.3
    
    if keyword_set(blackwhite) then shap05color = 'dark gray' else shap05color = 'dark green'
    shap05sym = 122
    shap05psize = 1.8
    
    if keyword_set(blackwhite) then sava05color = 'dark gray' else sava05color = 'dark orchid'
    sava05sym = 115
    sava05psize = 1.8
    
; initialize plotting variables

    @'xyrange_sfrs'

    Zsun_old = 8.9
    Zsun_new = 8.7

; output sfr calibration structure for Hb and [O II]

    template_sfr = {$
      loglb:  0.0, $
      mb:     0.0, $
      p25:    0.0, $
      p50:    0.0, $
      p75:    0.0, $
      mean:   0.0, $
      stddev: 0.0}
    
; clean up old PS and PNG files

    if keyword_set(cleanpng) then begin
       splog, 'Deleting all PNG and PS files in '+pspath
       spawn, ['rm -f '+pspath+'*.png*'], /sh
       spawn, ['rm -f '+pspath+'*ps'], /sh
    endif

; ##################################################
; BEGIN PAPER PLOTS
; ##################################################

; ------------------------------------------------------------
; E(B-V) [Ha/Hb] vs E(B-V) [Hb/Hg] - SDSS + Integrated
; ------------------------------------------------------------

    psname = 'sdss_ebv_hahb_vs_ebv_hbhg'
    im_openclose, pspath+psname, postscript=postscript, xsize=8.5, ysize=8.0, encapsulated=encapsulated
    
    pagemaker, nx=1, ny=1, xspace=0, yspace=0.0, width=6.6, height=6.6, $
      xmargin=[1.5,0.4], ymargin=[0.3,1.1], xpage=8.5, ypage=8.0, $
      position=pos, /normal

    if keyword_set(blackwhite) then begin
       atlascolor = 'dark gray'
       nfgscolor = 'black'
    endif else begin
       atlascolor = 'blue'
       nfgscolor = 'red'
    endelse
       
    indx = where((sdssnodust.ebv_hahb_err gt 0.0) and (sdssnodust.ebv_hbhg_err gt 0.0) and $
      (sdssdust.h_gamma[0]/sdssdust.h_gamma[1] gt hgsnrcut) and (sdssdust.h_beta_ew[0] gt hbewcut),nindx)
;   indx = where((sdssnodust.ebv_hahb_err gt 0.0) and (sdssnodust.ebv_hbhg_err gt 0.0) and $
;     (sdssdust.h_gamma[0]/sdssdust.h_gamma[1] gt hgsnrcut),nindx)
    
    x = sdssnodust[indx].ebv_hahb
    xerr = sdssnodust[indx].ebv_hahb_err

    y = sdssnodust[indx].ebv_hbhg
    yerr = sdssnodust[indx].ebv_hbhg_err

    indxatlas = where((atlasnodust.ebv_hahb_err gt 0.0) and (atlasnodust.ebv_hbhg_err gt 0.0) and $
      (atlasdust.h_gamma[0]/atlasdust.h_gamma[1] gt hgsnrcut) and (atlasdust.h_beta_ew[0] gt hbewcut),nindxatlas)
    
    xatlas = atlasnodust[indxatlas].ebv_hahb
    xerratlas = atlasnodust[indxatlas].ebv_hahb_err

    yatlas = atlasnodust[indxatlas].ebv_hbhg
    yerratlas = atlasnodust[indxatlas].ebv_hbhg_err
    
    indxnfgs = where((nfgsnodust.ebv_hahb_err gt 0.0) and (nfgsnodust.ebv_hbhg_err gt 0.0) and $
      (nfgsdust.h_gamma[0]/nfgsdust.h_gamma[1] gt hgsnrcut) and (nfgsdust.h_beta_ew[0] gt hbewcut),nindxnfgs)
    
    xnfgs = nfgsnodust[indxnfgs].ebv_hahb
    xerrnfgs = nfgsnodust[indxnfgs].ebv_hahb_err

    ynfgs = nfgsnodust[indxnfgs].ebv_hbhg
    yerrnfgs = nfgsnodust[indxnfgs].ebv_hbhg_err

    splog, 'SDSS, ATLAS, NFGS:'
    junk = im_stats(x-y,sigrej=3.0,/verbose)
    junk = im_stats(yatlas-xatlas,sigrej=3.0,/verbose)
    junk = im_stats(ynfgs-xnfgs,sigrej=3.0,/verbose)
    
    xtitle = 'E(B-V) using H\alpha/H\beta [mag]'
    ytitle = 'E(B-V) using H\beta/H\gamma [mag]'
    
    xrange = [-0.01,1.0]
    yrange = xrange

    plotsym, 0, 0.2, /fill
    sdss_lineplot, x, y, xerr, yerr, plottype=3, postscript=postscript, $
      xtitle=xtitle, ytitle=ytitle, xrange=xrange, yrange=yrange, legendtype=0, $
      yfrac=1.3, position=pos[*,0], blackwhite=blackwhite;, psym=8;, $
;     xatlas=xatlas, yatlas=yatlas, xerratlas=xerratlas, yerratlas=yerratlas, $
;     xnfgs=xnfgs, ynfgs=ynfgs, xerrnfgs=xerrnfgs, yerrnfgs=yerrnfgs
;   djs_oplot, [!x.crange[0],0.9], [!y.crange[0],0.9], line=0, thick=postthick3
    djs_oplot, !x.crange, !y.crange, line=0, thick=postthick

    col = fsc_color(atlascolor,!d.table_size-3)
    im_symbols, 106, fill=1, psize=0.9, color=col, thick=symthick
    oploterror, xatlas, yatlas, xerratlas, yerratlas, psym=8, errcolor=col
    
    col = fsc_color(nfgscolor,!d.table_size-4)
    im_symbols, 105, fill=1, psize=1.2, color=col, thick=symthick
    oploterror, xnfgs, ynfgs, xerrnfgs, yerrnfgs, psym=8, errcolor=col

;   legend, [textoidl('S/N(H\gamma) > '+string(hgsnrcut,format='(I0)')),$
;     textoidl('EW(H\beta) > '+string(hbewcut,format='(I0)'))+' '+angstrom()], $
;     /right, /top, box=0, charsize=charsize_5, charthick=postthick, clear=keyword_set(postscript), $
;     spacing=2.0
    
    xyouts, 0.4, 0.95, textoidl('S/N(H\gamma) > '+string(hgsnrcut,format='(I0)')), $
      /data, charsize=charsize_5, charthick=postthick
    xyouts, 0.4, 0.90, textoidl('EW(H\beta) > '+string(hbewcut,format='(I0)'))+' '+angstrom(), $
      /data, charsize=charsize_5, charthick=postthick
    
    im_openclose, postscript=postscript, /close

; ------------------------------------------------------------
; 3-panel - L(B) vs [O II]/Ha - Integrated
; ------------------------------------------------------------

    psname = 'lb_vs_oiiha_3panel'
    im_openclose, pspath+psname, postscript=postscript, xsize=8.5, ysize=9.5, encapsulated=encapsulated

    pagemaker, nx=1, ny=3, height=[2.5,2.5,2.5], width=7.0, xmargin=[1.1,0.4], $
      ymargin=[0.9,1.1], xspace=0, yspace=0, xpage=8.5, ypage=9.5, $
      position=pos, /normal

    medbin = 0.5
    minpts = 10
    minx = 7.75

    indx = where((atlasdust.oii_3727[0]/atlasdust.oii_3727[1] gt snrcut) and $
      (atlasnodust.b_lum_obs gt -900.0) and (atlasnodust.ehbha_err gt 0.0),nindx)

    x = atlasnodust[indx].b_lum_obs
    xerr = atlasnodust[indx].b_lum_obs_err
    xabs = atlasnodust[indx].m_b_obs

    indxnfgs = where((nfgsdust.oii_3727[0]/nfgsdust.oii_3727[1] gt snrcut) and $
      (nfgsnodust.b_lum_obs gt -900.0) and (nfgsnodust.ehbha_err gt 0.0),nindxnfgs)

    xnfgs = nfgsnodust[indxnfgs].b_lum_obs
    xerrnfgs = nfgsnodust[indxnfgs].b_lum_obs_err

    xtitle = 'log L(B) [L(B)'+sunsymbol()+']'
    ytitle = 'log ([O II]/H\alpha)'

    xrange = LBrange
    yrange = oiihacorrange

; ##########################################################
; Panel 1: [O II]_obs, Ha_obs
; ##########################################################

    y1 = atlasdust[indx].oii_3727[0]
    y1err = atlasdust[indx].oii_3727[1]

    y2 = atlasdust[indx].h_alpha[0]
    y2err = atlasdust[indx].h_alpha[1]

    y = alog10(y1/y2)
    yerr = im_compute_error(y1,y1err,y2,y2err,/log)

; NFGS    
    
    y1nfgs = nfgsdust[indxnfgs].oii_3727[0]
    y1errnfgs = nfgsdust[indxnfgs].oii_3727[1]

    y2nfgs = nfgsdust[indxnfgs].h_alpha[0]
    y2errnfgs = nfgsdust[indxnfgs].h_alpha[1]

    ynfgs = alog10(y1nfgs/y2nfgs)
    yerrnfgs = im_compute_error(y1nfgs,y1errnfgs,y2nfgs,y2errnfgs,/log)

    stats = im_stats([y,ynfgs],/verbose,/baremin)
    xstr = strtrim(string(stats.median_rej,format='(F12.2)'),2)+' ('+$
      strtrim(string(stats.mean_rej,format='(F12.2)'),2)+$
      '\pm'+strtrim(string(stats.sigma_rej,format='(F12.2)'),2)+')'
    rcor = r_correlate([x,xnfgs],[y,ynfgs],zd=zd,probd=probd)
    splog, 'Spearman rank test: ', rcor, probd, zd

    ytitle = 'log ([O II]/H\alpha)_{obs}'

    atlas1d_lineplot, x, y, xerr, yerr, plottype=3, postscript=postscript, $
      ytitle=ytitle, xrange=xrange, yrange=yrange, legendtype=0, blackwhite=blackwhite, $
      xstyle=11, /right, /top, xnfgs=xnfgs, ynfgs=ynfgs, xerrnfgs=xerrnfgs, $
      yerrnfgs=yerrnfgs, position=pos[*,0], charsize=charsize_4, xtickname=replicate(' ',10)
    axis, /xaxis, xrange=interpol(xabs,x,!x.crange), xthick=postthick, xstyle=1, $
      charsize=charsize_4, charthick=postthick
    im_xyouts_title, xtitle='M_{B} [mag]', charsize=charsize_4, charthick=postthick

    legend, '(a)', /left, /top, box=0, charsize=charsize_3, charthick=postthick
    legend, textoidl(xstr), /right, /top, box=0, charsize=charsize_3, charthick=postthick

; ##########################################################
; Panel 2: [O II]_cor, Ha_cor
; ##########################################################

;   y1 = atlasnodust[indx].oii_3727[0]
;   y1err = atlasnodust[indx].oii_3727[1]
;   y2 = atlasnodust[indx].h_alpha[0]
;   y2err = atlasnodust[indx].h_alpha[1]
;   y = alog10(y1/y2)
;   yerr = im_compute_error(y1,y1err,y2,y2err,/log)

    oii_dust = atlasdust[indx].oii_3727[0]
    oii_dust_err = atlasdust[indx].oii_3727[1]
    oii_nodust = atlasnodust[indx].oii_3727[0]

    hb_dust = atlasdust[indx].h_beta[0]
    hb_dust_err = atlasdust[indx].h_beta[1]

    ha_dust = atlasdust[indx].h_alpha[0]
    ha_dust_err = atlasdust[indx].h_alpha[1]
    ha_nodust = atlasnodust[indx].h_alpha[0]

    y1 = oii_nodust/ha_nodust
    y1err = sqrt((y1*(1+bigbprime_oii-bigbprime_ha)*(ha_dust_err/ha_dust))^2 + $
      (y1*(bigbprime_oii-bigbprime_ha)*(hb_dust_err/hb_dust))^2 + (y1*oii_dust_err/oii_dust)^2)

    y = alog10(y1)
    yerr = y1err/y1/alog(10.0)

;   y1nfgs = nfgsnodust[indxnfgs].oii_3727[0]
;   y1errnfgs = nfgsnodust[indxnfgs].oii_3727[1]
;   y2nfgs = nfgsnodust[indxnfgs].h_alpha[0]
;   y2errnfgs = nfgsnodust[indxnfgs].h_alpha[1]
;   ynfgs = alog10(y1nfgs/y2nfgs)
;   yerrnfgs = im_compute_error(y1nfgs,y1errnfgs,y2nfgs,y2errnfgs,/log)

    oiinfgs_dust = nfgsdust[indxnfgs].oii_3727[0]
    oiinfgs_dust_err = nfgsdust[indxnfgs].oii_3727[1]
    oiinfgs_nodust = nfgsnodust[indxnfgs].oii_3727[0]

    hbnfgs_dust = nfgsdust[indxnfgs].h_beta[0]
    hbnfgs_dust_err = nfgsdust[indxnfgs].h_beta[1]

    hanfgs_dust = nfgsdust[indxnfgs].h_alpha[0]
    hanfgs_dust_err = nfgsdust[indxnfgs].h_alpha[1]
    hanfgs_nodust = nfgsnodust[indxnfgs].h_alpha[0]

    y1nfgs = oiinfgs_nodust/hanfgs_nodust
    y1errnfgs = sqrt((y1nfgs*(1+bigbprime_oii-bigbprime_ha)*(hanfgs_dust_err/hanfgs_dust))^2 + $
      (y1nfgs*(bigbprime_oii-bigbprime_ha)*(hbnfgs_dust_err/hbnfgs_dust))^2 + (y1nfgs*oiinfgs_dust_err/oiinfgs_dust)^2)
    
    ynfgs = alog10(y1nfgs)
    yerrnfgs = y1errnfgs/y1nfgs/alog(10.0)
    
    stats = im_stats([y,ynfgs],/verbose,/baremin,/no_head)
    xstr = strtrim(string(stats.median_rej,format='(F12.2)'),2)+' ('+$
      strtrim(string(stats.mean_rej,format='(F12.2)'),2)+$
      '\pm'+strtrim(string(stats.sigma_rej,format='(F12.2)'),2)+')'
    rcor = r_correlate([x,xnfgs],[y,ynfgs],zd=zd,probd=probd)
    splog, 'Spearman rank test: ', rcor, probd, zd

    ytitle = 'log ([O II]/H\alpha)_{cor}'

    atlas1d_lineplot, x, y, xerr, yerr, plottype=3, postscript=postscript, $
      xtitle='', ytitle=ytitle, xrange=xrange, yrange=yrange, legendtype=0, /noerase, $
      xstyle=3, /right, /top, xtickname=replicate(' ',10), blackwhite=blackwhite, $
      xnfgs=xnfgs, ynfgs=ynfgs, xerrnfgs=xerrnfgs, $
      yerrnfgs=yerrnfgs, position=pos[*,1], charsize=charsize_4
    
    print & running = im_medxbin([x,xnfgs],[y,ynfgs],medbin,minx=minx,minpts=minpts,/verbose) & print

    legend, '(b)', /left, /top, box=0, charsize=charsize_3, charthick=postthick
    legend, textoidl(xstr), /right, /top, box=0, charsize=charsize_3, charthick=postthick

;   w = where(10^[x,xnfgs] ge 2E10 and 10^[x,xnfgs] le 5E10)
;;  w = where([x,xnfgs] ge 10.25 and [x,xnfgs] le 10.75)
;   junk = im_stats(([y,ynfgs])[w],/ver)

; ##########################################################
; Panel 3: [O II]_obs, SFR(Ha)
; ##########################################################

;   y1 = atlasdust[indx].oii_3727_lum[0]
;   y1err = atlasdust[indx].oii_3727_lum[1]
;   y2 = atlasnodust[indx].h_alpha_lum[0]
;   y2err = atlasnodust[indx].h_alpha_lum[1]
;   y = y1 - y2 - loghaconst - hasfrconstoffset
;   yerr = sqrt(y1err^2 + y2err^2)

; account for the covariance in the Ha reddening correction

    oii_dust = atlasdust[indx].oii_3727[0]
    oii_dust_err = atlasdust[indx].oii_3727[1]

    hb_dust = atlasdust[indx].h_beta[0]
    hb_dust_err = atlasdust[indx].h_beta[1]

    ha_dust = atlasdust[indx].h_alpha[0]
    ha_dust_err = atlasdust[indx].h_alpha[1]
    ha_nodust = atlasnodust[indx].h_alpha[0]
    ha_nodust_err = sqrt((ha_nodust*(1-bigbprime_ha)*ha_dust_err/ha_dust)^2 + (bigbprime_ha*ha_nodust*hb_dust_err/hb_dust)^2)

    y = alog10(oii_dust/ha_nodust) - loghaconst - hasfrconstoffset
    yerr = im_compute_error(oii_dust,oii_dust_err,ha_nodust,ha_nodust_err,/log)

; NFGS    
    
    oiinfgs_dust = nfgsdust[indxnfgs].oii_3727[0]
    oiinfgs_dust_err = nfgsdust[indxnfgs].oii_3727[1]

    hbnfgs_dust = nfgsdust[indxnfgs].h_beta[0]
    hbnfgs_dust_err = nfgsdust[indxnfgs].h_beta[1]

    hanfgs_dust = nfgsdust[indxnfgs].h_alpha[0]
    hanfgs_dust_err = nfgsdust[indxnfgs].h_alpha[1]
    hanfgs_nodust = nfgsnodust[indxnfgs].h_alpha[0]
    hanfgs_nodust_err = sqrt((hanfgs_nodust*(1-bigbprime_ha)*hanfgs_dust_err/hanfgs_dust)^2 + $
      (bigbprime_ha*hanfgs_nodust*hbnfgs_dust_err/hbnfgs_dust)^2)

    ynfgs = alog10(oiinfgs_dust/hanfgs_nodust) - loghaconst - hasfrconstoffset
    yerrnfgs = im_compute_error(oiinfgs_dust,oiinfgs_dust_err,hanfgs_nodust,hanfgs_nodust_err,/log)

;   y1nfgs = nfgsdust[indxnfgs].oii_3727_lum[0]
;   y1errnfgs = nfgsdust[indxnfgs].oii_3727_lum[1]
;   y2nfgs = nfgsnodust[indxnfgs].h_alpha_lum[0]
;   y2errnfgs = nfgsnodust[indxnfgs].h_alpha_lum[1]
;   ynfgs = y1nfgs - y2nfgs - loghaconst - hasfrconstoffset
;   yerrnfgs = sqrt(y1errnfgs^2 + y2errnfgs^2)
    
    stats = im_stats([y,ynfgs],/verbose,/baremin,/no_head)
    xstr = strtrim(string(stats.median_rej,format='(F12.2)'),2)+' ('+$
      strtrim(string(stats.mean_rej,format='(F12.2)'),2)+$
      '\pm'+strtrim(string(stats.sigma_rej,format='(F12.2)'),2)+')'
    rcor = r_correlate([x,xnfgs],[y,ynfgs],zd=zd,probd=probd)
    splog, 'Spearman rank test: ', rcor, probd, zd

    ytitle = 'log [10^{-41} L([O II])_{obs}/\psi(H\alpha)]'
    yrange2 = sfroiiharange - hasfrconstoffset

    atlas1d_lineplot, x, y, xerr, yerr, plottype=3, postscript=postscript, $
      xtitle=xtitle, ytitle=ytitle, xrange=xrange, yrange=yrange2, legendtype=0, /noerase, $
      xstyle=3, /right, /top, blackwhite=blackwhite, $
      xnfgs=xnfgs, ynfgs=ynfgs, xerrnfgs=xerrnfgs, $
      yerrnfgs=yerrnfgs, position=pos[*,2], charsize=charsize_4

    print & running = im_medxbin([x,xnfgs],[y,ynfgs],medbin,minx=minx,minpts=minpts,/verbose) & print

    legend, '(c)', /left, /top, box=0, charsize=charsize_3, charthick=postthick
    legend, textoidl(xstr), /right, /top, box=0, charsize=charsize_3, charthick=postthick

;   w = where([x,xnfgs] lt 9)
;   w2 = where([x,xnfgs] gt 10 and [x,xnfgs] lt 11)
;   w3 = where([x,xnfgs] gt 9.5)
;   junk = im_stats(([y,ynfgs])[w],/ver)
;   junk = im_stats(([y,ynfgs])[w2],/ver)
;   junk = im_stats(([y,ynfgs])[w3],/ver)

;   w = where([x,xnfgs] ge 10.0 and [x,xnfgs] le 10.5)
;   junk = im_stats(([y,ynfgs])[w],/verbose)

    im_openclose, postscript=postscript, /close

; ------------------------------------------------------------
; 2-panel - D4000 vs L(U)/L(Ha) - Integrated
; ------------------------------------------------------------

    psname = 'd4000_vs_lu_lha_2panel'
    im_openclose, pspath+psname, postscript=postscript, xsize=8.5, ysize=7.5, encapsulated=encapsulated

    pagemaker, nx=1, ny=2, height=[3.0,3.0], width=7.0, xmargin=[1.2,0.3], $
      ymargin=[0.4,1.1], xspace=0, yspace=0, xpage=8.5, ypage=7.5, $
      position=pos, /normal

    medbin = 0.1
    minpts = 10
    minx = 1.0

    indx = where((atlasdust.d4000_narrow_model[1] gt 0.0) and (atlasdust.synth_u_obs gt -900) and $
      (atlasnodust.ebv_hahb_err gt 0),nindx)

    x = atlasdust[indx].d4000_narrow_model[0]
    xerr = atlasdust[indx].d4000_narrow_model[1]
    xlum = atlasdust[indx].b_lum_obs

    indxnfgs = where((nfgsdust.d4000_narrow_model[1] gt 0.0) and (nfgsdust.synth_u_obs gt -900) and $
      (nfgsnodust.ebv_hahb_err gt 0),nindxnfgs)

    xnfgs = nfgsdust[indxnfgs].d4000_narrow_model[0]
    xerrnfgs = nfgsdust[indxnfgs].d4000_narrow_model[1]
    xlumnfgs = nfgsdust[indxnfgs].b_lum_obs

    xtitle = 'D_{n}(4000)'
    ytitle = 'log [L(U)/L(H\alpha)] '

    xrange = D4000range2
    yrange = Uhacorrange
    
; ##########################################################
; Panel 1: U observed, Ha observed
; ##########################################################

    y1 = Uconstant*10^(-0.4*atlasdust[indx].synth_u_obs)
    y1err = 0.4*sqrt(atlasdust[indx].synth_u_obs_err^2+Uminerr^2)*alog(10.0)*y1
    
    y2 = atlasdust[indx].h_alpha[0]
    y2err = atlasdust[indx].h_alpha[1]

    y = alog10(y1/y2)
    yerr = im_compute_error(y1,y1err,y2,y2err,/log)

    y1nfgs = Uconstant*10^(-0.4*nfgsdust[indxnfgs].synth_u_obs)
    y1errnfgs = 0.4*sqrt(nfgsdust[indxnfgs].synth_u_err^2+Uminerr^2)*alog(10.0)*y1
    
    y2nfgs = nfgsdust[indxnfgs].h_alpha[0]
    y2errnfgs = nfgsdust[indxnfgs].h_alpha[1]

    ynfgs = alog10(y1nfgs/y2nfgs)
    yerrnfgs = im_compute_error(y1nfgs,y1errnfgs,y2nfgs,y2errnfgs,/log)
    
    stats = im_stats([y,ynfgs],/verbose,/baremin)
    xstr = strtrim(string(stats.median_rej,format='(F12.2)'),2)+' ('+$
      strtrim(string(stats.mean_rej,format='(F12.2)'),2)+$
      '\pm'+strtrim(string(stats.sigma_rej,format='(F12.2)'),2)+')'
    rcor = r_correlate([x,xnfgs],[y,ynfgs],zd=zd,probd=probd)
    splog, 'Spearman rank test: ', rcor, probd, zd

    ytitle = 'log [L(U)/L(H\alpha)]_{obs} '

    atlas1d_lineplot, x, y, xerr, yerr, plottype=3, postscript=postscript, $
      ytitle=ytitle, xrange=xrange, yrange=yrange, legendtype=0, $
      /right, /top, xnfgs=xnfgs, ynfgs=ynfgs, xerrnfgs=xerrnfgs, blackwhite=blackwhite, $
      yerrnfgs=yerrnfgs, position=pos[*,0], charsize=charsize_5, xtickname=replicate(' ',10)

    print & running = im_medxbin([x,xnfgs],[y,ynfgs],medbin,minx=minx,minpts=minpts,/verbose) & print
;   oploterror, running.binctr, running.medy, replicate(running.binsz/2.0,n_elements(running.medy)), $
;     running.sigy, ps=-4, xsty=3, ysty=3, yrange=yrange

    legend, '(a)', /left, /top, box=0, charsize=charsize_3, charthick=postthick
    legend, textoidl(xstr), /right, /top, box=0, charsize=charsize_3, charthick=postthick

; ##########################################################
; Panel 2: U observed, Ha corrected
; ##########################################################

    u_dust = Uconstant*10^(-0.4*atlasdust[indx].synth_u_obs)
    u_dust_err = 0.4*sqrt(atlasdust[indx].synth_u_obs_err^2+Uminerr^2)*alog(10.0)*u_dust

    hb_dust = atlasdust[indx].h_beta[0]
    hb_dust_err = atlasdust[indx].h_beta[1]

    ha_dust = atlasdust[indx].h_alpha[0]
    ha_dust_err = atlasdust[indx].h_alpha[1]
    ha_nodust = atlasnodust[indx].h_alpha[0]
    ha_nodust_err = sqrt((ha_nodust*(1-bigbprime_ha)*ha_dust_err/ha_dust)^2 + (bigbprime_ha*ha_nodust*hb_dust_err/hb_dust)^2)

    y = alog10(u_dust/ha_nodust) - loghaconst - Uhasfrconstoffset
    yerr = im_compute_error(u_dust,u_dust_err,ha_nodust,ha_nodust_err,/log)
    
; the code below is wrong!!!  the U-band luminosity was normalized to
; the U-band solar luminosity, while the H-alpha luminosity was
; normalized to the bolometric luminosity!

;   cor = (atlasdust[indx].synth_u-atlasdust[indx].synth_u_obs) > 0.0 ; <-- NOTE!
;   y1 = atlasdust[indx].synth_u_lum + 0.4*(mbolsun-uinfo.solarmags) - 0.4*cor
;;  y1 = atlasdust[indx].synth_u_lum - 0.4*cor + 0.4*(mbolsun-uinfo.solarmags)
;   y1err = atlasdust[indx].synth_u_lum_err
;   y2 = atlasnodust[indx].h_alpha_lum[0]
;   y2err = atlasnodust[indx].h_alpha_lum[1]
;   y = y1 - y2 - loghaconst - Uhasfrconstoffset
;   yerr = sqrt(y1err^2 + y2err^2)

;;  y2 = atlasdust[indx].synth_u_lum
;;  y3 = atlasdust[indx].synth_u_lum + 0.4*(mbolsun-uinfo.solarmags)
;;  y4 = 0.4*(mbolsun-atlasdust[indx].synth_M_u)   
    
; NFGS    
    
    unfgs_dust = Uconstant*10^(-0.4*nfgsdust[indxnfgs].synth_u_obs)
    unfgs_dust_err = 0.4*sqrt(nfgsdust[indxnfgs].synth_u_obs_err^2+Uminerr^2)*alog(10.0)*unfgs_dust

    hbnfgs_dust = nfgsdust[indxnfgs].h_beta[0]
    hbnfgs_dust_err = nfgsdust[indxnfgs].h_beta[1]

    hanfgs_dust = nfgsdust[indxnfgs].h_alpha[0]
    hanfgs_dust_err = nfgsdust[indxnfgs].h_alpha[1]
    hanfgs_nodust = nfgsnodust[indxnfgs].h_alpha[0]
    hanfgs_nodust_err = sqrt((hanfgs_nodust*(1-bigbprime_ha)*hanfgs_dust_err/hanfgs_dust)^2 + $
      (bigbprime_ha*hanfgs_nodust*hbnfgs_dust_err/hbnfgs_dust)^2)

    ynfgs = alog10(unfgs_dust/hanfgs_nodust) - loghaconst - Uhasfrconstoffset
    yerrnfgs = im_compute_error(unfgs_dust,unfgs_dust_err,hanfgs_nodust,hanfgs_nodust_err,/log)

;   cor = (nfgsdust[indxnfgs].synth_u-nfgsdust[indxnfgs].synth_u_obs) > 0.0 ; <-- NOTE!
;   y1nfgs = nfgsdust[indxnfgs].synth_u_lum + 0.4*(mbolsun-uinfo.solarmags) - 0.4*cor
;   y1errnfgs = nfgsdust[indxnfgs].synth_u_lum_err
;   y2nfgs = nfgsnodust[indxnfgs].h_alpha_lum[0]
;   y2errnfgs = nfgsnodust[indxnfgs].h_alpha_lum[1]
;   ynfgs = y1nfgs - y2nfgs - loghaconst - Uhasfrconstoffset
;   yerrnfgs = sqrt(y1errnfgs^2 + y2errnfgs^2)
    
; stats
    
    stats = im_stats([y,ynfgs],/verbose)
    xstr = strtrim(string(stats.median_rej,format='(F12.2)'),2)+' ('+$
      strtrim(string(stats.mean_rej,format='(F12.2)'),2)+$
      '\pm'+strtrim(string(stats.sigma_rej,format='(F12.2)'),2)+')'
    rcor = r_correlate([x,xnfgs],[y,ynfgs],zd=zd,probd=probd)
    splog, 'Spearman rank test: ', rcor, probd, zd

; derive the SFR conversion factor for the full sample

    constant = 10.0^(-(stats.median_rej+Uhasfrconstoffset))
    constant_err = constant*alog(10.0)*stats.sigma_rej
    
    splog, 'SFR(U) coefficients: '
    print, constant, constant_err

; derive the SFR conversion factor for just luminous galaxies 

    w = where([xlum,xlumnfgs] gt 9.5)
    lumstats = im_stats(([y,ynfgs])[w])
    
    constant2 = 10.0^(-(lumstats.median_rej+Uhasfrconstoffset))
    constant2_err = constant*alog(10.0)*lumstats.sigma_rej
    
    splog, 'SFR(U) coefficients for log LB > 9.5: '
    print, constant2, constant2_err

; now move on with your life    

    ytitle = 'log [10^{-42} L(U)_{obs}/\psi(H\alpha)]'
    yrange2 = sfrUharange - Uhasfrconstoffset

    atlas1d_lineplot, x, y, xerr, yerr, plottype=3, postscript=postscript, $
      xrange=xrange, yrange=yrange2, legendtype=0, /noerase, xtitle=xtitle, $
      /right, /top, xnfgs=xnfgs, ynfgs=ynfgs, xerrnfgs=xerrnfgs, $
      yerrnfgs=yerrnfgs, position=pos[*,1], charsize=charsize_5, ytitle=ytitle, $
      blackwhite=blackwhite

    print & running = im_medxbin([x,xnfgs],[y,ynfgs],medbin,minx=minx,minpts=minpts,/verbose) & print
;   oploterror, running.binctr, running.medy, replicate(running.binsz/2.0,n_elements(running.medy)), $
;     running.sigy, ps=-4, xsty=3, ysty=3, yrange=yrange
;   plot, [xlum,xlumnfgs], [y,ynfgs]-interpol(running.medy,running.binctr,[x,xnfgs]), ps=4, xr=[8.5,11.5]
;   plot, [atlasnodust[indx].ebv_hahb,nfgsnodust[indxnfgs].ebv_hahb], $
;     [y,ynfgs]-interpol(running.medy,running.binctr,[x,xnfgs]), ps=4, xr=[0,1.2]

;   plot, [xlum,xlumnfgs], [y,ynfgs], ps=4, xr=[8.5,11.5], yrange=yrange2, xsty=3, ysty=3
;   print & running = im_medxbin([xlum,xlumnfgs],[y,ynfgs],1.5,minx=7.75,minpts=minpts,/verbose) & print
;   oploterror, running.binctr, running.medy, replicate(running.binsz/2.0,n_elements(running.medy)), $
;     running.sigy, ps=-4, xsty=3, ysty=3, yrange=yrange

; overplot the Cram et al. SFR calibration

    legend, '(b)', /left, /top, box=0, charsize=charsize_3, charthick=postthick
    legend, textoidl(xstr), /right, /top, box=0, charsize=charsize_3, charthick=postthick

    im_openclose, postscript=postscript, /close

; ------------------------------------------------------------
; 2-panel - D4000 vs L(U)/L(Ha) - SDSS
; ------------------------------------------------------------

    psname = 'sdss_d4000_vs_lu_lha_2panel'
    im_openclose, pspath+psname, postscript=postscript, xsize=8.5, ysize=7.5, encapsulated=encapsulated

    pagemaker, nx=1, ny=2, height=[3.0,3.0], width=7.0, xmargin=[1.2,0.3], $
      ymargin=[0.4,1.1], xspace=0, yspace=0, xpage=8.5, ypage=7.5, $
      position=pos, /normal

    medbin = 0.1
    minpts = 10
    minx = 1.0

    indx = where((sdssdust.model_d4000_narrow[1] gt 0.0) and (sdssancillary.fiber_u gt -900) and $
      finite(sdssancillary.fiber_u) and (sdssnodust.ebv_hahb_err gt 0),nindx)

    x = sdssdust[indx].model_d4000_narrow[0]
    xerr = sdssdust[indx].model_d4000_narrow[1]
    xlum = sdssancillary[indx].b_lum

    xtitle = 'D_{n}(4000)'
    ytitle = 'log [L(U)/L(H\alpha)] '

    xrange = D4000range2
    yrange = Uhacorrange
    
; ##########################################################
; Panel 1: U observed, Ha observed
; ##########################################################

    y1 = Uconstant*10^(-0.4*sdssancillary[indx].fiber_U)
    y1err = 0.4*sdssancillary[indx].fiber_U_err*alog(10.0)*y1
    
    y2 = sdssdust[indx].h_alpha[0]
    y2err = sdssdust[indx].h_alpha[1]

    y = alog10(y1/y2)
    yerr = im_compute_error(y1,y1err,y2,y2err,/log)

    stats = im_stats(y,/verbose,/baremin)
    xstr = strtrim(string(stats.median_rej,format='(F12.2)'),2)+' ('+$
      strtrim(string(stats.mean_rej,format='(F12.2)'),2)+$
      '\pm'+strtrim(string(stats.sigma_rej,format='(F12.2)'),2)+')'
    rcor = r_correlate(x,y,zd=zd,probd=probd)
    splog, 'Spearman rank test: ', rcor, probd, zd

    ytitle = 'log [L(U)/L(H\alpha)]_{obs} '

    sdss_lineplot, x, y, xerr, yerr, plottype=3, postscript=postscript, $
      ytitle=ytitle, xrange=xrange, yrange=yrange, legendtype=0, $
      /right, /top, position=pos[*,0], charsize=charsize_5, xtickname=replicate(' ',10)

    print & running = im_medxbin(x,y,medbin,minx=minx,minpts=minpts,/verbose) & print
;   oploterror, running.binctr, running.medy, replicate(running.binsz/2.0,n_elements(running.medy)), $
;     running.sigy, ps=-4, xsty=3, ysty=3, yrange=yrange

    legend, '(a)', /left, /top, box=0, charsize=charsize_3, charthick=postthick
    legend, textoidl(xstr), /right, /top, box=0, charsize=charsize_3, charthick=postthick

; ##########################################################
; Panel 2: U observed, Ha corrected
; ##########################################################

    u_dust = Uconstant*10^(-0.4*sdssancillary[indx].fiber_u)
    u_dust_err = 0.4*sdssancillary[indx].fiber_u_err*alog(10.0)*u_dust

    hb_dust = sdssdust[indx].h_beta[0]
    hb_dust_err = sdssdust[indx].h_beta[1]

    ha_dust = sdssdust[indx].h_alpha[0]
    ha_dust_err = sdssdust[indx].h_alpha[1]
    ha_nodust = sdssnodust[indx].h_alpha[0]
    ha_nodust_err = sqrt((ha_nodust*(1-bigbprime_ha)*ha_dust_err/ha_dust)^2 + (bigbprime_ha*ha_nodust*hb_dust_err/hb_dust)^2)

    y = alog10(u_dust/ha_nodust) - loghaconst - Uhasfrconstoffset
    yerr = im_compute_error(u_dust,u_dust_err,ha_nodust,ha_nodust_err,/log)

; wrong!    
    
;   y1 = sdssancillary[indx].fiber_U_lum
;   y1err = sdssancillary[indx].fiber_U_lum_err
;   y2 = sdssnodust[indx].h_alpha_lum[0]
;   y2err = sdssnodust[indx].h_alpha_lum[1]
;   y = y1 - y2 - loghaconst - Uhasfrconstoffset
;   yerr = sqrt(y1err^2 + y2err^2)

    stats = im_stats(y,/verbose)
    xstr = strtrim(string(stats.median_rej,format='(F12.2)'),2)+' ('+$
      strtrim(string(stats.mean_rej,format='(F12.2)'),2)+$
      '\pm'+strtrim(string(stats.sigma_rej,format='(F12.2)'),2)+')'
    rcor = r_correlate(x,y,zd=zd,probd=probd)
    splog, 'Spearman rank test: ', rcor, probd, zd

; derive the SFR conversion factor for the full sample

    constant = 10.0^(-(stats.median_rej+Uhasfrconstoffset))
    constant_err = constant*alog(10.0)*stats.sigma_rej
    
    splog, 'SFR(U) coefficients: '
    print, constant, constant_err

; derive the SFR conversion factor for just luminous galaxies 

    w = where(xlum gt 9.5)
    lumstats = im_stats(y[w])
    
    constant2 = 10.0^(-(lumstats.median_rej+Uhasfrconstoffset))
    constant2_err = constant*alog(10.0)*lumstats.sigma_rej
    
    splog, 'SFR(U) coefficients for log LB > 9.5: '
    print, constant2, constant2_err

; make the plot    
    
    ytitle = 'log [10^{-42} L(U)_{obs}/\psi(H\alpha)]'
    yrange2 = sfrUharange - Uhasfrconstoffset

    sdss_lineplot, x, y, xerr, yerr, plottype=3, postscript=postscript, xtitle=xtitle, $
      xrange=xrange, yrange=yrange2, legendtype=0, /noerase, $
      /right, /top, position=pos[*,1], charsize=charsize_5, ytitle=ytitle
    
    print & running = im_medxbin(x,y,medbin,minx=minx,minpts=minpts,/verbose) & print
;   oploterror, running.binctr, running.medy, replicate(running.binsz/2.0,n_elements(running.medy)), $
;     running.sigy, ps=-4, xsty=3, ysty=3, yrange=yrange
;   plot, xlum, y-interpol(running.medy,running.binctr,x), ps=3, xr=[8.5,11.5]
;   plot, sdssnodust[indx].ebv_hahb, y-interpol(running.medy,running.binctr,x), ps=3, xr=[0,1.1]

    legend, '(b)', /left, /top, box=0, charsize=charsize_3, charthick=postthick
    legend, textoidl(xstr), /right, /top, box=0, charsize=charsize_3, charthick=postthick

    im_openclose, postscript=postscript, /close

; ------------------------------------------------------------
; SFR(Ha) vs various SFR([O II]) calibrations - residual plots
; ------------------------------------------------------------

    psname = 'sfr_ha_vs_sfr_oii_4panel'
    im_openclose, pspath+psname, postscript=postscript, xsize=8.5, ysize=6.5, encapsulated=encapsulated

    pagemaker, nx=2, ny=2, yspace=0, xspace=0, width=3.5*[1,1], height=2.5*[1,1], $
      xmargin=[1.1,0.4], ymargin=[0.4,1.1], xpage=8.5, ypage=6.5, position=pos, /normal

    if keyword_set(blackwhite) then begin
       sfcolor = 'gray'
       agncolor = 'black'
    endif else begin
       sfcolor = 'sky blue'
       agncolor = 'navy'
    endelse
       
    xtitle = 'log \psi(H\alpha) ['+sfr_units()+']'
    ytitle = 'log \psi([O II])/\psi(H\alpha)'
    
    xrange = sfrharange
    yrange = residrange_sfrs
    
; ##########################################################
; Panel 1
; ##########################################################

; SF
    
    indx = where((sfrs.sfr_ha gt -900) and (sfrs.sfr_oii_k98 gt -900),nindx)

    x = sfrs[indx].sfr_ha
    xerr = sfrs[indx].sfr_ha_err
    
    y1 = sfrs[indx].sfr_ha
    y1err = sfrs[indx].sfr_ha_err
    
    y2 = sfrs[indx].sfr_oii_k98
    y2err = sfrs[indx].sfr_oii_k98_err

    resid = y2-y1
    residerr = sqrt(y1err^2+y2err^2)

; AGN
    
    indx_agn = where((sfrs_agn.sfr_ha gt -900) and (sfrs_agn.sfr_oii_k98 gt -900),nindx_agn)

    x_agn = sfrs_agn[indx_agn].sfr_ha
    xerr_agn = sfrs_agn[indx_agn].sfr_ha_err
    
    y1_agn = sfrs_agn[indx_agn].sfr_ha
    y1err_agn = sfrs_agn[indx_agn].sfr_ha_err
    
    y2_agn = sfrs_agn[indx_agn].sfr_oii_k98
    y2err_agn = sfrs_agn[indx_agn].sfr_oii_k98_err

    resid_agn = y2_agn-y1_agn
    residerr_agn = sqrt(y1err_agn^2+y2err_agn^2)

; make the plot    

    stats = im_stats(resid,/verbose,/baremin)
    xstr = strtrim(string(stats.median_rej,format='(F12.2)'),2)+' ('+$
      strtrim(string(stats.mean_rej,format='(F12.2)'),2)+$
      '\pm'+strtrim(string(stats.sigma_rej,format='(F12.2)'),2)+')'
    
; SF    
    
    atlas1d_lineplot, x, resid, xerr, residerr, plottype=1, postscript=postscript, $
      xtitle='', ytitle=ytitle, xrange=xrange, yrange=yrange, legendtype=0, blackwhite=blackwhite, $
      charsize=charsize_5, xtickname=replicate(' ',10), position=pos[*,0], $
      atlaspsize=0.6, ytickinterval=1.0, atlascolor=sfcolor

; AGN    
    
    atlas1d_lineplot, x_agn, resid_agn, xerr_agn, residerr_agn, blackwhite=blackwhite, $
      postscript=postscript, atlaspsize=0.6, atlasfill=0, /overplot, atlascolor=agncolor, $
      coloroffset=10
    djs_oplot, !x.crange, [0,0], line=0, thick=postthick

    legend, '(a)', /left, /top, box=0, charsize=charsize_3, $
      charthick=postthick, clear=keyword_set(postscript)
    legend, 'Kennicutt (1998)', /left, /bottom, box=0, charsize=charsize_2, $
      charthick=postthick, clear=keyword_set(postscript)
    legend, textoidl(xstr), /right, /top, box=0, charsize=charsize_2, $
      charthick=postthick, clear=keyword_set(postscript)

; ##########################################################
; Panel 2
; ##########################################################

; SF
    
    indx = where((sfrs.sfr_ha gt -900) and (sfrs.sfr_oii_hbhg gt -900),nindx)

    x = sfrs[indx].sfr_ha
    xerr = sfrs[indx].sfr_ha_err
    
    y1 = sfrs[indx].sfr_ha
    y1err = sfrs[indx].sfr_ha_err
    
    y2 = sfrs[indx].sfr_oii_hbhg
    y2err = sfrs[indx].sfr_oii_hbhg_err

    resid = y2-y1
    residerr = sqrt(y1err^2+y2err^2)

; AGN
    
    indx_agn = where((sfrs_agn.sfr_ha gt -900) and (sfrs_agn.sfr_oii_hbhg gt -900),nindx_agn)

    x_agn = sfrs_agn[indx_agn].sfr_ha
    xerr_agn = sfrs_agn[indx_agn].sfr_ha_err
    
    y1_agn = sfrs_agn[indx_agn].sfr_ha
    y1err_agn = sfrs_agn[indx_agn].sfr_ha_err
    
    y2_agn = sfrs_agn[indx_agn].sfr_oii_hbhg
    y2err_agn = sfrs_agn[indx_agn].sfr_oii_hbhg_err

    resid_agn = y2_agn-y1_agn
    residerr_agn = sqrt(y1err_agn^2+y2err_agn^2)

; make the plot    

    stats = im_stats(resid,/verbose,/baremin)
    xstr = strtrim(string(stats.median_rej,format='(F12.2)'),2)+' ('+$
      strtrim(string(stats.mean_rej,format='(F12.2)'),2)+$
      '\pm'+strtrim(string(stats.sigma_rej,format='(F12.2)'),2)+')'
    
; SF    
    
    atlas1d_lineplot, x, resid, xerr, residerr, plottype=1, postscript=postscript, $
      xtitle='', ytitle='', xrange=xrange, yrange=yrange, legendtype=0, blackwhite=blackwhite, $
      charsize=charsize_5, xtickname=replicate(' ',10), position=pos[*,1], $
      ytickname=replicate(' ',10), atlaspsize=0.6, $
      ytickinterval=1.0, /noerase, atlascolor=sfcolor

; AGN    
    
    atlas1d_lineplot, x_agn, resid_agn, xerr_agn, residerr_agn, blackwhite=blackwhite, $
      postscript=postscript, atlaspsize=0.6, atlasfill=0, /overplot, atlascolor=agncolor, $
      coloroffset=10
    djs_oplot, !x.crange, [0,0], line=0, thick=postthick
    
    legend, '(b)', /left, /top, box=0, charsize=charsize_3, $
      charthick=postthick, clear=keyword_set(postscript)
    legend, textoidl('Extinction-Corrected using H\beta/H\gamma'), /left, /bottom, box=0, charsize=charsize_2, $
      charthick=postthick, clear=keyword_set(postscript)
    legend, textoidl(xstr), /right, /top, box=0, charsize=charsize_2, $
      charthick=postthick, clear=keyword_set(postscript)

; ##########################################################
; Panel 3
; ##########################################################

; SF
    
    indx = where((sfrs.sfr_ha gt -900) and (sfrs.sfr_oii_k04_theory gt -900),nindx)

    x = sfrs[indx].sfr_ha
    xerr = sfrs[indx].sfr_ha_err
    
    y1 = sfrs[indx].sfr_ha
    y1err = sfrs[indx].sfr_ha_err
    
    y2 = sfrs[indx].sfr_oii_k04_theory
    y2err = sfrs[indx].sfr_oii_k04_theory_err

    resid = y2-y1
    residerr = sqrt(y1err^2+y2err^2)

; AGN
    
    indx_agn = where((sfrs_agn.sfr_ha gt -900) and (sfrs_agn.sfr_oii_k04_theory gt -900),nindx_agn)

    x_agn = sfrs_agn[indx_agn].sfr_ha
    xerr_agn = sfrs_agn[indx_agn].sfr_ha_err
    
    y1_agn = sfrs_agn[indx_agn].sfr_ha
    y1err_agn = sfrs_agn[indx_agn].sfr_ha_err
    
    y2_agn = sfrs_agn[indx_agn].sfr_oii_k04_theory
    y2err_agn = sfrs_agn[indx_agn].sfr_oii_k04_theory_err

    resid_agn = y2_agn-y1_agn
    residerr_agn = sqrt(y1err_agn^2+y2err_agn^2)

; make the plot    

    stats = im_stats(resid,/verbose,/baremin)
    xstr = strtrim(string(stats.median_rej,format='(F12.2)'),2)+' ('+$
      strtrim(string(stats.mean_rej,format='(F12.2)'),2)+$
      '\pm'+strtrim(string(stats.sigma_rej,format='(F12.2)'),2)+')'
    
; SF    
    
    atlas1d_lineplot, x, resid, xerr, residerr, plottype=1, postscript=postscript, $
      xtitle=xtitle, ytitle=ytitle, xrange=xrange, yrange=yrange, legendtype=0, blackwhite=blackwhite, $
      charsize=charsize_5, position=pos[*,2], atlaspsize=0.6, $
      ytickinterval=1.0, /noerase, atlascolor=sfcolor

; AGN    
    
    atlas1d_lineplot, x_agn, resid_agn, xerr_agn, residerr_agn, blackwhite=blackwhite, $
      postscript=postscript, atlaspsize=0.6, atlasfill=0, /overplot, atlascolor=agncolor, $
      coloroffset=10
    djs_oplot, !x.crange, [0,0], line=0, thick=postthick

    legend, '(c)', /left, /top, box=0, charsize=charsize_3, $
      charthick=postthick, clear=keyword_set(postscript)
    legend, 'Kewley et al. (2004)', /left, /bottom, box=0, charsize=charsize_2, $
      charthick=postthick, clear=keyword_set(postscript)
    legend, textoidl(xstr), /right, /top, box=0, charsize=charsize_2, $
      charthick=postthick, clear=keyword_set(postscript)

; ##########################################################
; Panel 4
; ##########################################################

; SF
    
    indx = where((sfrs.sfr_ha gt -900) and (sfrs.sfr_oii_best gt -900),nindx)

    x = sfrs[indx].sfr_ha
    xerr = sfrs[indx].sfr_ha_err
    
    y1 = sfrs[indx].sfr_ha
    y1err = sfrs[indx].sfr_ha_err
    
    y2 = sfrs[indx].sfr_oii_best
    y2err = sfrs[indx].sfr_oii_best_err

    resid = y2-y1
    residerr = sqrt(y1err^2+y2err^2)

; AGN
    
    indx_agn = where((sfrs_agn.sfr_ha gt -900) and (sfrs_agn.sfr_oii_best gt -900),nindx_agn)

    x_agn = sfrs_agn[indx_agn].sfr_ha
    xerr_agn = sfrs_agn[indx_agn].sfr_ha_err
    
    y1_agn = sfrs_agn[indx_agn].sfr_ha
    y1err_agn = sfrs_agn[indx_agn].sfr_ha_err
    
    y2_agn = sfrs_agn[indx_agn].sfr_oii_best
    y2err_agn = sfrs_agn[indx_agn].sfr_oii_best_err

    resid_agn = y2_agn-y1_agn
    residerr_agn = sqrt(y1err_agn^2+y2err_agn^2)

; make the plot    

    stats = im_stats(resid,/verbose,/baremin)
    xstr = strtrim(string(stats.median_rej,format='(F12.2)'),2)+' ('+$
      strtrim(string(stats.mean_rej,format='(F12.2)'),2)+$
      '\pm'+strtrim(string(stats.sigma_rej,format='(F12.2)'),2)+')'
    
; SF    
    
    atlas1d_lineplot, x, resid, xerr, residerr, plottype=1, postscript=postscript, $
      xtitle=xtitle, ytitle='', xrange=xrange, yrange=yrange, legendtype=0, blackwhite=blackwhite, $
      charsize=charsize_5, position=pos[*,3], ytickname=replicate(' ',10), $
      atlaspsize=0.6, ytickinterval=1.0, /noerase, atlascolor=sfcolor

; AGN    
    
    atlas1d_lineplot, x_agn, resid_agn, xerr_agn, residerr_agn, blackwhite=blackwhite, $
      postscript=postscript, atlaspsize=0.6, atlasfill=0, /overplot, atlascolor=agncolor, $
      coloroffset=10
    djs_oplot, !x.crange, [0,0], line=0, thick=postthick

    legend, '(d)', /left, /top, box=0, charsize=charsize_3, $
      charthick=postthick, clear=keyword_set(postscript)
    legend, textoidl('L(B) vs \psi/L_{obs} Median Relation'), /left, /bottom, box=0, charsize=charsize_2, $
      charthick=postthick, clear=keyword_set(postscript)
    legend, textoidl(xstr), /right, /top, box=0, charsize=charsize_2, $
      charthick=postthick, clear=keyword_set(postscript)

    im_openclose, postscript=postscript, /close    

; ------------------------------------------------------------
; High-Redshift L(B) vs Hb/Ha
; ------------------------------------------------------------

    psname = 'highz_lb_vs_hahb'
    im_openclose, pspath+psname, postscript=postscript, xsize=8.5, ysize=8.3, encapsulated=encapsulated

    pagemaker, nx=1, ny=1, xspace=0, yspace=0.0, width=6.3, height=6.3, $
      xmargin=[1.1,1.1], ymargin=[0.9,1.1], xpage=8.5, ypage=8.3, $
      position=pos, /normal

; SF galaxies    
    
    cut = where((atlasdust.b_lum_obs gt -900.0) and (atlasdust.h_beta_ew_uncor[0] ge ewcut2))
    lineratio, atlasdust[cut], '', '', 'H_ALPHA', 'H_BETA', $
      dum1, dum2, y, yerr, index=indx, nindex=nindx, snrcut=snrcut, $
      xsyserr=syserr, ysyserr=syserr

    x = atlasdust[cut[indx]].b_lum_obs
    xerr = atlasdust[cut[indx]].b_lum_obs_err
    xabs = atlasdust[cut[indx]].m_b_obs
    yebv = atlasnodust[cut[indx]].ebv_hahb

    cutnfgs = where((nfgsdust.b_lum_obs gt -900.0) and (nfgsdust.h_beta_ew_uncor[0] ge ewcut2))
    lineratio, nfgsdust[cutnfgs], '', '', 'H_ALPHA', 'H_BETA', $
      dum1, dum2, ynfgs, yerrnfgs, index=indxnfgs, snrcut=snrcut, $
      xsyserr=nfgssyserr, ysyserr=nfgssyserr
    xnfgs = nfgsdust[cutnfgs[indxnfgs]].b_lum_obs
    xerrnfgs = nfgsdust[cutnfgs[indxnfgs]].b_lum_obs_err

    xbig = [x,xnfgs]
    ybig = [y,ynfgs]
    xerrbig = [xerr,xerrnfgs]
    yerrbig = [yerr,yerrnfgs]

; AGN
    
    cut_agn = where((atlasdust_agn.b_lum_obs gt -900.0) and (atlasdust_agn.h_beta_ew_uncor[0] ge ewcut2))
    lineratio, atlasdust_agn[cut_agn], '', '', 'H_ALPHA', 'H_BETA', $
      dum1, dum2, y_agn, yerr_agn, index=indx_agn, nindex=nindx_agn, snrcut=snrcut, $
      xsyserr=syserr, ysyserr=syserr

    x_agn = atlasdust_agn[cut_agn[indx_agn]].b_lum_obs
    xerr_agn = atlasdust_agn[cut_agn[indx_agn]].b_lum_obs_err

    cutnfgs_agn = where((nfgsdust_agn.b_lum_obs gt -900.0) and (nfgsdust_agn.h_beta_ew_uncor[0] ge ewcut2))
    lineratio, nfgsdust_agn[cutnfgs_agn], '', '', 'H_ALPHA', 'H_BETA', $
      dum1, dum2, ynfgs_agn, yerrnfgs_agn, index=indxnfgs_agn, snrcut=snrcut, $
      xsyserr=nfgssyserr, ysyserr=nfgssyserr
    xnfgs_agn = nfgsdust_agn[cutnfgs_agn[indxnfgs_agn]].b_lum_obs
    xerrnfgs_agn = nfgsdust_agn[cutnfgs_agn[indxnfgs_agn]].b_lum_obs_err

    xbig_agn = [x_agn,xnfgs_agn]
    ybig_agn = [y_agn,ynfgs_agn]
    xerrbig_agn = [xerr_agn,xerrnfgs_agn]
    yerrbig_agn = [yerr_agn,yerrnfgs_agn]

; average error bar
    
    avgerrlocal = djs_median([yerrbig,yerrbig_agn])

; ####################
; Maier (2005)
; ####################
    
    cutm05 = where(m05.b_lum gt -900.0)
    lineratio, m05[cutm05], '', '', 'H_ALPHA', 'H_BETA', dum1, dum2, $
      ym05, yerrm05, index=indxm05, nindex=nindxm05, snrcut=snrcut_highz

    xm05 = m05[cutm05[indxm05]].b_lum
    xerrm05 = m05[cutm05[indxm05]].b_lum_err

    lhbm05 = m05[cutm05[indxm05]].h_beta_lum[0] + alog10(lsun)
    hasfrm05 = m05nodust[cutm05[indxm05]].sfr_h_alpha

    m05agn = where(strmatch(m05[cutm05[indxm05]].bpt_class,'*agn*',/fold) eq 1B)

    avgerrm05 = djs_median(yerrm05)
    splog, 'Maier: redshift = ', minmax(m05[cutm05[indxm05]].z_obj), median(m05[cutm05[indxm05]].z_obj)

    niceprint, xm05, ym05, m05[cutm05[indxm05]].bpt_class

; ####################
; Shapley et al. (2005)
; ####################
    
    cutshap05 = where(shap05.b_lum gt -900.0)
    lineratio, shap05[cutshap05], '', '', 'H_ALPHA', 'H_BETA', dum1, dum2, $
      yshap05, yerrshap05, index=indxshap05, nindex=nindxshap05, snrcut=snrcut_highz

    xshap05 = shap05[cutshap05[indxshap05]].b_lum
    xerrshap05 = shap05[cutshap05[indxshap05]].b_lum_err

    lhbshap05 = shap05[cutshap05[indxshap05]].h_beta_lum[0] + alog10(lsun)
    hasfrshap05 = shap05nodust[cutshap05[indxshap05]].sfr_h_alpha

    shap05agn = where(strmatch(shap05[cutshap05[indxshap05]].bpt_class,'*agn*',/fold) eq 1B)

    avgerrshap05 = djs_median(yerrshap05)
    splog, 'Shapley: redshift = ', minmax(shap05[cutshap05[indxshap05]].z_obj), median(shap05[cutshap05[indxshap05]].z_obj)

    niceprint, xshap05, yshap05, shap05[cutshap05[indxshap05]].bpt_class

; ####################
; Savaglio et al. (2005)
; ####################

    indxsava05 = where((sava05.b_lum gt -900.0) and (sava05nodust.ebv_hbhg_err gt 0.0) and (sava05nodust.ebv_hbhg gt 0.0))
    factor = -0.4*(kl_ha-kl_hb)*sava05nodust[indxsava05].ebv_hbhg
    factor_err = -0.4*(kl_ha-kl_hb)*sava05nodust[indxsava05].ebv_hbhg_err

    ysava05 = alog10(hahb) + factor
    yerrsava05 = factor_err/factor/alog(10.0)
    
    xsava05 = sava05[indxsava05].b_lum
    xerrsava05 = sava05[indxsava05].b_lum_err

    avgerrsava05 = djs_median(yerrsava05)
    splog, 'Savaglio: redshift = ', minmax(sava05[indxsava05].z_obj), median(sava05[indxsava05].z_obj)

; ####################
; Liang et al. (2004)
; ####################

    indxliang04 = where((liang04.b_lum gt -900.0) and (liang04nodust.ebv_hbhg_err gt 0.0))
    factor = -0.4*(kl_ha-kl_hb)*liang04nodust[indxliang04].ebv_hbhg
    factor_err = -0.4*(kl_ha-kl_hb)*liang04nodust[indxliang04].ebv_hbhg_err

    yliang04 = alog10(hahb) + factor
    yerrliang04 = factor_err/factor/alog(10.0)
    
    xliang04 = liang04[indxliang04].b_lum
    xerrliang04 = liang04[indxliang04].b_lum_err

    lirliang04 = liang04[indxliang04].lir
    good = where(liang04[indxliang04].lir gt -900.0)
    splog, 'Liang: L(IR) = ', minmax(liang04[indxliang04[good]].lir), median(liang04[indxliang04[good]].lir)

;   g = where(atlasdust.l_ir_l_b gt -900 and atlasnodust.ebv_hahb_err gt 0.0)
;   plot, alog10(atlasdust[g].l_ir_l_b), atlasnodust[g].ebv_hahb, ps=4, xsty=3, ysty=3
;   djs_oplot, lirliang04[good]-(xliang04[good]+0.4*(4.74-5.42)), yliang04[good], ps=7, color='red'

    avgerrliang04 = djs_median(yerrliang04)
    splog, 'Liang: redshift = ', minmax(liang04[indxliang04].z_obj), median(liang04[indxliang04].z_obj)

; ####################
    
    xtitle = 'log L(B) [L(B)'+sunsymbol()+']'
    ytitle = 'log (H\alpha/H\beta)_{obs}'

    xrange = LBrange
    yrange = hahbrange

    djs_plot, [0], [0], /nodata, xthick=postthick, ythick=postthick, charsize=charsize_8, $
      charthick=postthick, thick=postthick, xtitle=xtitle, ytitle=ytitle, $
      ymargin=[4,3], xstyle=11, ystyle=11, xrange=xrange, yrange=yrange, position=pos[*,0]
    axis, /xaxis, xrange=interpol(xabs,x,!x.crange), xthick=postthick, xstyle=1, $
      charsize=charsize_8, charthick=postthick
    axis, /yaxis, yrange=interpol(yebv,y,!y.crange), ythick=postthick, ystyle=1, $
      charsize=charsize_8, charthick=postthick
    im_xyouts_title, xtitle='M_{B} [mag]', charsize=charsize_8, charthick=postthick
    im_xyouts_title, ytitle='E(B-V) [mag]', charsize=charsize_8, charthick=postthick, xspacing=9.0

    im_symbols, localsym, psize=localpsize, /fill, color=fsc_color(localcolor,!d.table_size-75)
    djs_oplot, xbig, ybig, ps=8, color=fsc_color(localcolor,!d.table_size-75)

    im_symbols, localsym, psize=localpsize, fill=0, color=fsc_color(localcolor,!d.table_size-75), thick=symthick
    djs_oplot, xbig_agn, ybig_agn, ps=8, color=fsc_color(localcolor,!d.table_size-75)

; put liang on bottom because darkest

    im_symbols, liang04sym, psize=liang04psize, /fill, color=fsc_color(liang04color,!d.table_size-79), thick=postthick
    djs_oplot, xliang04, yliang04, ps=8, color=fsc_color(liang04color,!d.table_size-79), thick=postthick

    im_symbols, m05sym, psize=m05psize, /fill, color=fsc_color(m05color,!d.table_size-76), thick=postthick
    djs_oplot, xm05, ym05, ps=8, color=fsc_color(m05color,!d.table_size-76), thick=postthick

    im_symbols, shap05sym, psize=shap05psize, /fill, color=fsc_color(shap05color,!d.table_size-77), thick=postthick
    djs_oplot, xshap05, yshap05, ps=8, color=fsc_color(shap05color,!d.table_size-77), thick=postthick

    im_symbols, sava05sym, psize=sava05psize, /fill, color=fsc_color(sava05color,!d.table_size-78), thick=postthick
    djs_oplot, xsava05, ysava05, ps=8, color=fsc_color(sava05color,!d.table_size-78), thick=postthick

    djs_oplot, !x.crange, alog10(HaHb)*[1,1], line=0, thick=postthick

; overplot the average error bars for the high-z points and the local 
; sample

;   ypos = 0.25
    ypos = 0.9
    
    im_symbols, localsym, psize=localpsize, /fill, color=fsc_color(localcolor,!d.table_size-75)
    oploterror, 7.2, ypos, 0.0, avgerrlocal, /data, thick=postthick, $
      errthick=postthick, psym=8, errcolor=fsc_color(localcolor,!d.table_size-75)
       
    im_symbols, m05sym, psize=m05psize, /fill, color=fsc_color(m05color,!d.table_size-76), thick=postthick
    oploterror, 7.5, ypos, 0.0, avgerrm05, /data, thick=postthick, $
      errthick=postthick, psym=8, errcolor=fsc_color(m05color,!d.table_size-76)
    
    im_symbols, shap05sym, psize=shap05psize, /fill, color=fsc_color(shap05color,!d.table_size-77), thick=postthick
    oploterror, 7.8, ypos, 0.0, avgerrshap05, /data, thick=postthick, $
      errthick=postthick, psym=8, errcolor=fsc_color(shap05color,!d.table_size-77)
    
    im_symbols, sava05sym, psize=sava05psize, /fill, color=fsc_color(sava05color,!d.table_size-78), thick=postthick
    oploterror, 8.1, ypos, 0.0, avgerrsava05, /data, thick=postthick, $
      errthick=postthick, psym=8, errcolor=fsc_color(sava05color,!d.table_size-78)
    
    im_symbols, liang04sym, psize=liang04psize, /fill, color=fsc_color(liang04color,!d.table_size-79), thick=postthick
    oploterror, 8.4, ypos, 0.0, avgerrliang04, /data, thick=postthick, $
      errthick=postthick, psym=8, errcolor=fsc_color(liang04color,!d.table_size-79)

; give the scatter in the SFRs using our empirical Hb conversion
; versus the extinction-corrected Ha SFRs for the Maier and Shapley
; samples

    hbsfr = hb_sfr([xm05,xshap05],[lhbm05,lhbshap05],/log)
    hasfr = [hasfrm05,hasfrshap05]
    
    good = where(hasfr gt -900.0)
    hbsfr = hbsfr[good]
    hasfr = hasfr[good]
    stats = im_stats(hbsfr-hasfr,/verbose)
;   plot, hasfr, hbsfr-hasfr, ps=4, xsty=3, ysty=3

; legend

    label = ['Local Sample - Star-Forming','Local Sample - AGN',$
      'Maier et al. (2005)','Shapley et al. (2005)',$
      'Savaglio et al. (2005)','Liang et al. (2004)']
    psym = [localsym,localsym,m05sym,shap05sym,sava05sym,liang04sym]
    fill = [1,0,1,1,1,1]
    color = fsc_color([localcolor,localcolor,m05color,shap05color,sava05color,liang04color],!d.table_size-[75,75,76,77,78,79])
    postthick1 = postthick
    im_legend, textoidl(label), /left, /bottom, box=0, charsize=charsize_2, $
      charthick=postthick1, psym=psym, fill=fill, symsize=1.3, $
      spacing=1.8, thick=postthick1;, textcolor=djs_icolor(replicate('',6))

    im_openclose, postscript=postscript, /close    

; ------------------------------------------------------------
; L(B) versus L([O III])_obs/SFR(Ha)
; ------------------------------------------------------------
    
    psname = 'lb_vs_loiii_sfr_ha_2panel'
    im_openclose, pspath+psname, postscript=postscript, xsize=8.5, ysize=5.5, encapsulated=encapsulated

    pagemaker, nx=2, ny=1, height=3.5, width=3.5, xmargin=[1.2,0.3], $
      ymargin=[0.9,1.1], xspace=0, yspace=0, xpage=8.5, ypage=5.5, $
      position=pos, /normal

    medbin = 0.5
    minpts = 10
    minx = 7.75

    indx = where((atlasdust.oiii_5007[0]/atlasdust.oiii_5007[1] gt snrcut) and $
      (atlasnodust.b_lum_obs gt -900.0) and (atlasnodust.ehbha_err gt 0.0),nindx)

    x = atlasnodust[indx].b_lum_obs
    xerr = atlasnodust[indx].b_lum_obs_err
    xabs = atlasnodust[indx].m_b_obs

    indxnfgs = where((nfgsdust.oiii_5007[0]/nfgsdust.oiii_5007[1] gt snrcut) and $
      (nfgsnodust.b_lum_obs gt -900.0) and (nfgsnodust.ehbha_err gt 0.0),nindxnfgs)

    xnfgs = nfgsnodust[indxnfgs].b_lum_obs
    xerrnfgs = nfgsnodust[indxnfgs].b_lum_obs_err

    xtitle = 'log L(B) [L(B)'+sunsymbol()+']'
    ytitle = 'log [10^{-41} L([O III] \lambda5007)_{obs}/\psi(H\alpha)]'

    xrange = LBrange
    yrange = sfroiiiharange - hasfrconstoffset

; Atlas    

;   y1 = atlasdust[indx].oiii_5007_lum[0]
;   y1err = atlasdust[indx].oiii_5007_lum[1]
;   y2 = atlasnodust[indx].h_alpha_lum[0]
;   y2err = atlasnodust[indx].h_alpha_lum[1]
;   y = y1 - y2 - loghaconst - hasfrconstoffset
;   yerr = sqrt(y1err^2 + y2err^2)

; account for the covariance in the Ha reddening correction

    oiii_dust = atlasdust[indx].oiii_5007[0]
    oiii_dust_err = atlasdust[indx].oiii_5007[1]

    hb_dust = atlasdust[indx].h_beta[0]
    hb_dust_err = atlasdust[indx].h_beta[1]

    ha_dust = atlasdust[indx].h_alpha[0]
    ha_dust_err = atlasdust[indx].h_alpha[1]
    ha_nodust = atlasnodust[indx].h_alpha[0]
    ha_nodust_err = sqrt((ha_nodust*(1-bigbprime_ha)*ha_dust_err/ha_dust)^2 + (bigbprime_ha*ha_nodust*hb_dust_err/hb_dust)^2)

    y = alog10(oiii_dust/ha_nodust) - loghaconst - hasfrconstoffset
    yerr = im_compute_error(oiii_dust,oiii_dust_err,ha_nodust,ha_nodust_err,/log)

; NFGS    
    
;   y1nfgs = nfgsdust[indxnfgs].oiii_5007_lum[0]
;   y1errnfgs = nfgsdust[indxnfgs].oiii_5007_lum[1]
;   y2nfgs = nfgsnodust[indxnfgs].h_alpha_lum[0]
;   y2errnfgs = nfgsnodust[indxnfgs].h_alpha_lum[1]
;   ynfgs = y1nfgs - y2nfgs - loghaconst - hasfrconstoffset
;   yerrnfgs = sqrt(y1errnfgs^2 + y2errnfgs^2)
    
    oiiinfgs_dust = nfgsdust[indxnfgs].oiii_5007[0]
    oiiinfgs_dust_err = nfgsdust[indxnfgs].oiii_5007[1]

    hbnfgs_dust = nfgsdust[indxnfgs].h_beta[0]
    hbnfgs_dust_err = nfgsdust[indxnfgs].h_beta[1]

    hanfgs_dust = nfgsdust[indxnfgs].h_alpha[0]
    hanfgs_dust_err = nfgsdust[indxnfgs].h_alpha[1]
    hanfgs_nodust = nfgsnodust[indxnfgs].h_alpha[0]
    hanfgs_nodust_err = sqrt((hanfgs_nodust*(1-bigbprime_ha)*hanfgs_dust_err/hanfgs_dust)^2 + $
      (bigbprime_ha*hanfgs_nodust*hbnfgs_dust_err/hbnfgs_dust)^2)

    ynfgs = alog10(oiiinfgs_dust/hanfgs_nodust) - loghaconst - hasfrconstoffset
    yerrnfgs = im_compute_error(oiiinfgs_dust,oiiinfgs_dust_err,hanfgs_nodust,hanfgs_nodust_err,/log)

    stats = im_stats([y,ynfgs],/verbose,/baremin,/no_head)
    xstr = strtrim(string(stats.median_rej,format='(F12.2)'),2)+' ('+$
      strtrim(string(stats.mean_rej,format='(F12.2)'),2)+$
      '\pm'+strtrim(string(stats.sigma_rej,format='(F12.2)'),2)+')'
    rcor = r_correlate([x,xnfgs],[y,ynfgs],zd=zd,probd=probd)
    splog, 'Spearman rank test: ', rcor, probd, zd

    atlas1d_lineplot, x, y, xerr, yerr, plottype=3, postscript=postscript, $
      xtitle=xtitle, ytitle=ytitle, xrange=xrange, yrange=yrange, legendtype=0, $
      xstyle=11, /right, /top, /errorleft, blackwhite=blackwhite, $
      xnfgs=xnfgs, ynfgs=ynfgs, xerrnfgs=xerrnfgs, $
      yerrnfgs=yerrnfgs, position=pos[*,0], charsize=charsize_5
    axis, /xaxis, xrange=interpol(xabs,x,!x.crange), xthick=postthick, xstyle=1, $
      charsize=charsize_5, charthick=postthick
    im_xyouts_title, xtitle='M_{B} [mag]', charsize=charsize_5, charthick=postthick

    print & running = im_medxbin([x,xnfgs],[y,ynfgs],medbin,minx=minx,minpts=minpts,/verbose) & print

    legend, textoidl(xstr), /right, /top, box=0, charsize=charsize_3, charthick=postthick

;   w = where([x,xnfgs] lt 9)
;   w2 = where([x,xnfgs] le 9.5)
;   w3 = where([x,xnfgs] gt 9.5)
;   junk = im_stats(([y,ynfgs])[w],/ver)
;   junk = im_stats(([y,ynfgs])[w2],/ver)
;   junk = im_stats(([y,ynfgs])[w3],/ver)

;   w = where([x,xnfgs] ge 10.25 and [x,xnfgs] le 10.75)
;   junk = im_stats(([y,ynfgs])[w],/verbose)

; SDSS

    indx = where((sdssdust.oiii_5007[0]/sdssdust.oiii_5007[1] gt snrcut) and $
      (sdssancillary.b_lum gt -900.0) and (sdssnodust.ehbha_err gt 0.0),nindx)

    x = sdssancillary[indx].b_lum
    xerr = sdssancillary[indx].b_lum_err
    xabs = sdssancillary[indx].m_b

;   y1 = sdssdust[indx].oiii_5007_lum[0]
;   y1err = sdssdust[indx].oiii_5007_lum[1]
;   y2 = sdssnodust[indx].h_alpha_lum[0]
;   y2err = sdssnodust[indx].h_alpha_lum[1]
;   y = y1 - y2 - loghaconst - hasfrconstoffset
;   yerr = sqrt(y1err^2 + y2err^2)

    oiii_dust = sdssdust[indx].oiii_5007[0]
    oiii_dust_err = sdssdust[indx].oiii_5007[1]

    hb_dust = sdssdust[indx].h_beta[0]
    hb_dust_err = sdssdust[indx].h_beta[1]

    ha_dust = sdssdust[indx].h_alpha[0]
    ha_dust_err = sdssdust[indx].h_alpha[1]
    ha_nodust = sdssnodust[indx].h_alpha[0]
    ha_nodust_err = sqrt((ha_nodust*(1-bigbprime_ha)*ha_dust_err/ha_dust)^2 + (bigbprime_ha*ha_nodust*hb_dust_err/hb_dust)^2)

    y = alog10(oiii_dust/ha_nodust) - loghaconst - hasfrconstoffset
    yerr = im_compute_error(oiii_dust,oiii_dust_err,ha_nodust,ha_nodust_err,/log)

    stats = im_stats(y,/verbose,/baremin,/no_head)
    xstr = strtrim(string(stats.median_rej,format='(F12.2)'),2)+' ('+$
      strtrim(string(stats.mean_rej,format='(F12.2)'),2)+$
      '\pm'+strtrim(string(stats.sigma_rej,format='(F12.2)'),2)+')'

    xtitle = 'log L(B) [L(B)'+sunsymbol()+']'
    ytitle = 'log [10^{-41} L([O III] \lambda5007)_{obs}/\psi(H\alpha)]'

    xrange = LBrange
    yrange = sfroiiiharange - hasfrconstoffset

    sdss_lineplot, x, y, xerr, yerr, plottype=3, postscript=postscript, $
      xtitle=xtitle, ytitle='', xrange=xrange, yrange=yrange, legendtype=0, $
      charsize=charsize_5, xstyle=11, position=pos[*,1], /noerase, $
      ytickname=replicate(' ',10), /errorleft
    axis, /xaxis, xrange=interpol(xabs,x,!x.crange), xthick=postthick, xstyle=1, $
      charsize=charsize_5, charthick=postthick
    im_xyouts_title, xtitle='M_{B} [mag]', charsize=charsize_5, charthick=postthick
    
    legend, textoidl(xstr), /right, /top, box=0, charsize=charsize_3, charthick=postthick

    im_openclose, postscript=postscript, /close

; ------------------------------------------------------------
; 12+log(O/H) vs ([O II]/Ha)_cor [Integrated + SDSS]
; ------------------------------------------------------------

    psname = '12oh_oiiinii_niiha_vs_oiiha_2panel'
    im_openclose, pspath+psname, postscript=postscript, xsize=8.5, ysize=5.0, encapsulated=encapsulated

    pagemaker, nx=2, ny=1, height=3.5, width=3.5, xmargin=[1.2,0.3], $
      ymargin=[0.5,1.0], xspace=0, yspace=0, xpage=8.5, ypage=5.0, $
      position=pos, /normal

    if keyword_set(blackwhite) then ugridcolor = '' else ugridcolor = 'dark green'

    xtitle = '12 + log (O/H)'
    ytitle = 'log ([O II]/H\alpha)_{cor}'

    xrange = hiiohrange
    yrange = oiiharange

; Atlas

    indx = where((atlasdust.zstrong_12oh_oiiinii_niiha gt -900) and $
      (atlasnodust.ehbha_err gt 0.0) and (atlasdust.oii_3727[0]/atlasdust.oii_3727[1] gt snrcut),nindx)

    x = atlasdust[indx].zstrong_12oh_oiiinii_niiha
    xerr = atlasdust[indx].zstrong_12oh_oiiinii_niiha_err

    oii_dust = atlasdust[indx].oii_3727[0]
    oii_dust_err = atlasdust[indx].oii_3727[1]
    oii_nodust = atlasnodust[indx].oii_3727[0]

    hb_dust = atlasdust[indx].h_beta[0]
    hb_dust_err = atlasdust[indx].h_beta[1]

    ha_dust = atlasdust[indx].h_alpha[0]
    ha_dust_err = atlasdust[indx].h_alpha[1]
    ha_nodust = atlasnodust[indx].h_alpha[0]

    y1 = oii_nodust/ha_nodust
    y1err = sqrt((y1*(1+bigbprime_oii-bigbprime_ha)*(ha_dust_err/ha_dust))^2 + $
      (y1*(bigbprime_oii-bigbprime_ha)*(hb_dust_err/hb_dust))^2 + (y1*oii_dust_err/oii_dust)^2)

    y = alog10(y1)
    yerr = y1err/y1/alog(10.0)

; NFGS    
    
    indxnfgs = where((nfgsdust.zstrong_12oh_oiiinii_niiha gt -900) and $
      (nfgsnodust.ehbha_err gt 0.0) and (nfgsdust.oii_3727[0]/nfgsdust.oii_3727[1] gt snrcut),nindxnfgs)

    xnfgs = nfgsdust[indxnfgs].zstrong_12oh_oiiinii_niiha
    xerrnfgs = nfgsdust[indxnfgs].zstrong_12oh_oiiinii_niiha_err

    oiinfgs_dust = nfgsdust[indxnfgs].oii_3727[0]
    oiinfgs_dust_err = nfgsdust[indxnfgs].oii_3727[1]
    oiinfgs_nodust = nfgsnodust[indxnfgs].oii_3727[0]

    hbnfgs_dust = nfgsdust[indxnfgs].h_beta[0]
    hbnfgs_dust_err = nfgsdust[indxnfgs].h_beta[1]

    hanfgs_dust = nfgsdust[indxnfgs].h_alpha[0]
    hanfgs_dust_err = nfgsdust[indxnfgs].h_alpha[1]
    hanfgs_nodust = nfgsnodust[indxnfgs].h_alpha[0]

    y1nfgs = oiinfgs_nodust/hanfgs_nodust
    y1errnfgs = sqrt((y1nfgs*(1+bigbprime_oii-bigbprime_ha)*(hanfgs_dust_err/hanfgs_dust))^2 + $
      (y1nfgs*(bigbprime_oii-bigbprime_ha)*(hbnfgs_dust_err/hbnfgs_dust))^2 + (y1nfgs*oiinfgs_dust_err/oiinfgs_dust)^2)
    
    ynfgs = alog10(y1nfgs)
    yerrnfgs = y1errnfgs/y1nfgs/alog(10.0)
    
;   indx1 = where((atlasnodust.zstrong_12oh_oiiinii_niiha gt -900),nindx1)
;   lineratio, atlasnodust[indx1], '', '', 'OII_3727', 'H_ALPHA', $
;     dum1, dum2, y, yerr, index=indx2, nindex=nindx2, snrcut=snrcut, $
;     xsyserr=syserr, ysyserr=syserr
;   x = atlasdust[indx1[indx2]].zstrong_12oh_oiiinii_niiha
;   xerr = atlasdust[indx1[indx2]].zstrong_12oh_oiiinii_niiha_err
;   indx1nfgs = where(nfgsnodust.zstrong_12oh_oiiinii_niiha gt -900,nindx1nfgs)
;   lineratio, nfgsnodust[indx1nfgs], '', '', 'OII_3727', 'H_ALPHA', $
;     dum1, dum2, ynfgs, yerrnfgs, index=indx2nfgs, nindex=nindx2nfgs, snrcut=snrcut, $
;     xsyserr=nfgssyserr, ysyserr=nfgssyserr
;   xnfgs = nfgsdust[indx1nfgs[indx2nfgs]].zstrong_12oh_oiiinii_niiha
;   xerrnfgs = nfgsdust[indx1nfgs[indx2nfgs]].zstrong_12oh_oiiinii_niiha_err

    xbig = [x,xnfgs]
    ybig = [y,ynfgs]

;   w = where((xbig gt 8.3) and (xbig lt 8.5))
;   junk = im_stats(ybig[w],/verbose)

; HII regions    
    
    indx1 = where((hii.zstrong_12oh_oiiinii_niiha gt -900) and (hii.oii_h_alpha gt -900),nindx1)

    xregion = hii[indx1].zstrong_12oh_oiiinii_niiha
    xerrregion = hii[indx1].zstrong_12oh_oiiinii_niiha_err
    yregion = hii[indx1].oii_h_alpha
    yerrregion = hii[indx1].oii_h_alpha_err

    atlas1d_lineplot, x, y, xerr, yerr, plottype=3, postscript=postscript, $n
      xtitle=xtitle, ytitle=ytitle, xrange=xrange, yrange=yrange, legendtype=0, $
      /left, /top, position=pos[*,0], charsize=charsize_5, yfrac=1.5, /errorleft, $
      xregion=xregion, yregion=yregion, xerrregion=xerrregion, yerrregion=yerrregion, $
      xnfgs=xnfgs, ynfgs=ynfgs, xerrnfgs=xerrnfgs, yerrnfgs=yerrnfgs, blackwhite=blackwhite

    if keyword_set(postscript) then thisthick = 6.0 else thickthick = 2.0

    plot_kewley_grids, plotnumber=17, model=3, labeltype=4, /log12oh, /overplot, $
      /noZgrid, thick=3.0, Ugridcolor=ugridcolor, Ulinestyle=0, $
      postscript=postscript, Ulabelcolor='', Ucharsize=0.8, charthick=thisthick
    xyouts, 9.1, -0.18, 'log U', align=0.5, /data, charsize=0.9, charthick=thisthick

; overlay the three metallicity regions

    r12 = 8.15
    r23 = 8.7

    oplot, r12*[1,1], [-1.25,0.75], line=1, thick=postthick
    oplot, r23*[1,1], [-1.25,0.75], line=1, thick=postthick
 
;   xyouts, 7.7, 0.65, 'R1', align=0.5, /data, charsize=charsize_5, charthick=postthick
;   xyouts, 8.4, 0.65, 'R2', align=0.5, /data, charsize=charsize_5, charthick=postthick
;   xyouts, 9.0, 0.65, 'R3', align=0.5, /data, charsize=charsize_5, charthick=postthick

    region1 = where(xbig lt r12,nregion1)
    splog, 'Integrated Sample: Region 1, 2, 3 separately: '
    stats = im_stats(ybig[region1],sigrej=3.0,/verbose,/baremin)

    regioniiha = where((xbig gt r12) and (xbig lt r23),nregioniiha)
    stats = im_stats(ybig[regioniiha],sigrej=3.0,/no_head,/verbose,/baremin)

    region3 = where(xbig gt r23,nregion3)
    stats = im_stats(ybig[region3],sigrej=3.0,/no_head,/verbose,/baremin)

    splog, 'Integrated Sample: Regions 1&2&3 combined: '
    stats = im_stats(ybig,/verbose,/baremin,/no_head)

    splog, 'Integrated Sample: Regions 2&3 combined: '
    if (nregioniiha ne 0L) and (nregion3 ne 0L) then $
      stats = im_stats([ybig[regioniiha],ybig[region3]],/verbose,/baremin,/no_head)

; SDSS    

    indx = where((sdssdust.zstrong_12oh_oiiinii_niiha gt -900) and $
      (sdssnodust.ehbha_err gt 0.0) and (sdssdust.oii_3727[0]/sdssdust.oii_3727[1] gt snrcut),nindx)

    x = sdssdust[indx].zstrong_12oh_oiiinii_niiha
    xerr = sdssdust[indx].zstrong_12oh_oiiinii_niiha_err

    oii_dust = sdssdust[indx].oii_3727[0]
    oii_dust_err = sdssdust[indx].oii_3727[1]
    oii_nodust = sdssnodust[indx].oii_3727[0]

    hb_dust = sdssdust[indx].h_beta[0]
    hb_dust_err = sdssdust[indx].h_beta[1]

    ha_dust = sdssdust[indx].h_alpha[0]
    ha_dust_err = sdssdust[indx].h_alpha[1]
    ha_nodust = sdssnodust[indx].h_alpha[0]

    y1 = oii_nodust/ha_nodust
    y1err = sqrt((y1*(1+bigbprime_oii-bigbprime_ha)*(ha_dust_err/ha_dust))^2 + $
      (y1*(bigbprime_oii-bigbprime_ha)*(hb_dust_err/hb_dust))^2 + (y1*oii_dust_err/oii_dust)^2)

    y = alog10(y1)
    yerr = y1err/y1/alog(10.0)
    
;   indx1 = where(sdssdust.zstrong_12oh_oiiinii_niiha gt -900,nindx1)
;   lineratio, sdssnodust[indx1], '', '', 'OII_3727', 'H_ALPHA', $
;     dum1, dum2, y, yerr, index=indx2, nindex=nindx2, snrcut=snrcut
;   x = sdssdust[indx1[indx2]].zstrong_12oh_oiiinii_niiha
;   xerr = sdssdust[indx1[indx2]].zstrong_12oh_oiiinii_niiha_err

    xbig = x
    ybig = y
    
    sdss_lineplot, x, y, xerr, yerr, plottype=3, postscript=postscript, $
      xtitle=xtitle, ytitle='', xrange=xrange, yrange=yrange, legendtype=0, $
      /left, /top, position=pos[*,1], yfrac=1.5, /errorleft, $
      ytickname=replicate(' ',10), /noerase, charsize=charsize_5;, $
;     xregion=xregion, yregion=yregion, xerrregion=xerrregion, yerrregion=yerrregion

    if keyword_set(postscript) then thisthick = 6.0 else thickthick = 2.0

    plot_kewley_grids, plotnumber=17, model=3, labeltype=4, /log12oh, /overplot, $
      /noZgrid, thick=3.0, Ugridcolor=ugridcolor, Ulinestyle=0, $
      postscript=postscript, Ulabelcolor='', Ucharsize=0.8, charthick=thisthick
    xyouts, 9.1, -0.18, 'log U', align=0.5, /data, charsize=0.9, charthick=thisthick

; overlay the three metallicity regions

    r12 = 8.15
    r23 = 8.7

    oplot, r12*[1,1], [-1.25,0.75], line=1, thick=postthick
    oplot, r23*[1,1], [-1.25,0.75], line=1, thick=postthick
 
;   xyouts, 7.7, 0.65, 'R1', align=0.5, /data, charsize=charsize_5, charthick=postthick
;   xyouts, 8.4, 0.65, 'R2', align=0.5, /data, charsize=charsize_5, charthick=postthick
;   xyouts, 9.0, 0.65, 'R3', align=0.5, /data, charsize=charsize_5, charthick=postthick

    splog, 'SDSS: Region 1, 2, 3 separately: '
    region1 = where(xbig lt r12,nregion1)
    if (nregion1 ne 0L) then stats = im_stats(ybig[region1],/verbose,/baremin)

    regioniiha = where((xbig gt r12) and (xbig lt r23),nregioniiha)
    if (nregioniiha ne 0L) then stats = im_stats(ybig[regioniiha],/verbose,/baremin,/no_head)

    region3 = where(xbig gt r23,nregion3)
    if (nregion3 ne 0L) then stats = im_stats(ybig[region3],/verbose,/baremin,/no_head)

    splog, 'SDSS Sample: Regions 1&2&3 combined: '
    stats = im_stats(ybig,/verbose,/baremin,/no_head)

    splog, 'SDSS: Regions 2&3 combined: '
    if (nregioniiha ne 0L) and (nregion3 ne 0L) then $
      stats = im_stats([ybig[regioniiha],ybig[region3]],/verbose,/baremin,/no_head)
    
    im_openclose, postscript=postscript, /close    

; ------------------------------------------------------------
; 3-panel - L(B) vs [O II]/Ha - SDSS
; ------------------------------------------------------------

    psname = 'sdss_lb_vs_oiiha_3panel'
    im_openclose, pspath+psname, postscript=postscript, xsize=8.5, ysize=9.5, encapsulated=encapsulated

    pagemaker, nx=1, ny=3, height=[2.5,2.5,2.5], width=7.0, xmargin=[1.1,0.4], $
      ymargin=[0.9,1.1], xspace=0, yspace=0, xpage=8.5, ypage=9.5, $
      position=pos, /normal

    indx = where((sdssdust.oii_3727[0]/sdssdust.oii_3727[1] gt snrcut) and $
      (sdssancillary.b_lum gt -900.0) and (sdssnodust.ehbha_err gt 0),nindx)

    x = sdssancillary[indx].b_lum
    xerr = sdssancillary[indx].b_lum_err
    xabs = sdssancillary[indx].m_b

    xtitle = 'log L(B) [L(B)'+sunsymbol()+']'
    ytitle = 'log ([O II]/H\alpha)'

    xrange = LBrange
    yrange = oiihacorrange

; ##########################################################
; Panel 1: [O II]_obs, Ha_obs
; ##########################################################

    y1 = sdssdust[indx].oii_3727[0]
    y1err = sdssdust[indx].oii_3727[1]

    y2 = sdssdust[indx].h_alpha[0]
    y2err = sdssdust[indx].h_alpha[1]

    y = alog10(y1/y2)
    yerr = im_compute_error(y1,y1err,y2,y2err,/log)

    stats = im_stats(y,/verbose,/baremin)
    xstr = strtrim(string(stats.median_rej,format='(F12.2)'),2)+' ('+$
      strtrim(string(stats.mean_rej,format='(F12.2)'),2)+$
      '\pm'+strtrim(string(stats.sigma_rej,format='(F12.2)'),2)+')'

    ytitle = 'log ([O II]/H\alpha)_{obs}'

    sdss_lineplot, x, y, xerr, yerr, plottype=3, postscript=postscript, $
      ytitle=ytitle, xrange=xrange, yrange=yrange, legendtype=0, $
      xstyle=11, /right, /top, position=pos[*,0], charsize=charsize_4, xtickname=replicate(' ',10)
    axis, /xaxis, xrange=interpol(xabs,x,!x.crange), xthick=postthick, xstyle=1, $
      charsize=charsize_4, charthick=postthick
    im_xyouts_title, xtitle='M_{B} [mag]', charsize=charsize_4, charthick=postthick

    legend, '(a)', /left, /top, box=0, charsize=charsize_3, charthick=postthick
    legend, textoidl(xstr), /right, /top, box=0, charsize=charsize_3, charthick=postthick

; ##########################################################
; Panel 2: [O II]_cor, Ha_cor
; ##########################################################

;   y1 = sdssnodust[indx].oii_3727[0]
;   y1err = sdssnodust[indx].oii_3727[1]
;   y2 = sdssnodust[indx].h_alpha[0]
;   y2err = sdssnodust[indx].h_alpha[1]
;   y = alog10(y1/y2)
;   yerr = im_compute_error(y1,y1err,y2,y2err,/log)

    oii_dust = sdssdust[indx].oii_3727[0]
    oii_dust_err = sdssdust[indx].oii_3727[1]
    oii_nodust = sdssnodust[indx].oii_3727[0]

    hb_dust = sdssdust[indx].h_beta[0]
    hb_dust_err = sdssdust[indx].h_beta[1]

    ha_dust = sdssdust[indx].h_alpha[0]
    ha_dust_err = sdssdust[indx].h_alpha[1]
    ha_nodust = sdssnodust[indx].h_alpha[0]

    y1 = oii_nodust/ha_nodust
    y1err = sqrt((y1*(1+bigbprime_oii-bigbprime_ha)*(ha_dust_err/ha_dust))^2 + $
      (y1*(bigbprime_oii-bigbprime_ha)*(hb_dust_err/hb_dust))^2 + (y1*oii_dust_err/oii_dust)^2)

    y = alog10(y1)
    yerr = y1err/y1/alog(10.0)

    stats = im_stats(y,/verbose,/baremin,/no_head)
    xstr = strtrim(string(stats.median_rej,format='(F12.2)'),2)+' ('+$
      strtrim(string(stats.mean_rej,format='(F12.2)'),2)+$
      '\pm'+strtrim(string(stats.sigma_rej,format='(F12.2)'),2)+')'

    ytitle = 'log ([O II]/H\alpha)_{cor}'

    sdss_lineplot, x, y, xerr, yerr, plottype=3, postscript=postscript, $
      xtitle='', ytitle=ytitle, xrange=xrange, yrange=yrange, $
      legendtype=0, /noerase, xtickname=replicate(' ',10), $
      xstyle=3, /right, /top, position=pos[*,1], charsize=charsize_4
    
    legend, '(b)', /left, /top, box=0, charsize=charsize_3, charthick=postthick
    legend, textoidl(xstr), /right, /top, box=0, charsize=charsize_3, charthick=postthick

;   w = where(x ge 10.25 and x le 10.75)
    w = where(10^x ge 2E10 and 10^x le 5E10)
    junk = im_stats(y[w],/verbose)

; ##########################################################
; Panel 3: [O II]_obs, SFR(Ha)
; ##########################################################

;   y1 = sdssdust[indx].oii_3727_lum[0]
;   y1err = sdssdust[indx].oii_3727_lum[1]
;   y2 = sdssnodust[indx].h_alpha_lum[0]
;   y2err = sdssnodust[indx].h_alpha_lum[1]
;   y = y1 - y2 - loghaconst - hasfrconstoffset
;   yerr = sqrt(y1err^2 + y2err^2)

; account for the covariance in the Ha reddening correction

    oii_dust = sdssdust[indx].oii_3727[0]
    oii_dust_err = sdssdust[indx].oii_3727[1]

    hb_dust = sdssdust[indx].h_beta[0]
    hb_dust_err = sdssdust[indx].h_beta[1]

    ha_dust = sdssdust[indx].h_alpha[0]
    ha_dust_err = sdssdust[indx].h_alpha[1]
    ha_nodust = sdssnodust[indx].h_alpha[0]
    ha_nodust_err = sqrt((ha_nodust*(1-bigbprime_ha)*ha_dust_err/ha_dust)^2 + (bigbprime_ha*ha_nodust*hb_dust_err/hb_dust)^2)

    y = alog10(oii_dust/ha_nodust) - loghaconst - hasfrconstoffset
    yerr = im_compute_error(oii_dust,oii_dust_err,ha_nodust,ha_nodust_err,/log)

    stats = im_stats(y,/verbose,/baremin,/no_head)
    xstr = strtrim(string(stats.median_rej,format='(F12.2)'),2)+' ('+$
      strtrim(string(stats.mean_rej,format='(F12.2)'),2)+$
      '\pm'+strtrim(string(stats.sigma_rej,format='(F12.2)'),2)+')'

    ytitle = 'log [10^{-41} L([O II])_{obs}/\psi(H\alpha)]'
    yrange2 = sfroiiharange - hasfrconstoffset

    sdss_lineplot, x, y, xerr, yerr, plottype=3, postscript=postscript, $
      xtitle=xtitle, ytitle=ytitle, xrange=xrange, yrange=yrange2, legendtype=0, /noerase, $
      xstyle=3, /right, /top, position=pos[*,2], charsize=charsize_4
    
    legend, '(c)', /left, /top, box=0, charsize=charsize_3, charthick=postthick
    legend, textoidl(xstr), /right, /top, box=0, charsize=charsize_3, charthick=postthick

;   w = where(x ge 10.0 and x le 10.5)
;   junk = im_stats(y[w],/verbose)

    im_openclose, postscript=postscript, /close

; ------------------------------------------------------------
; 3-panel - L(B) vs Hb/Ha - SDSS
; ------------------------------------------------------------
    
    psname = 'sdss_lb_vs_hbha_3panel'
    im_openclose, pspath+psname, postscript=postscript, xsize=8.5, ysize=9.5, encapsulated=encapsulated

    pagemaker, nx=1, ny=3, height=[2.5,2.5,2.5], width=7.0, xmargin=[1.1,0.4], $
      ymargin=[0.9,1.1], xspace=0, yspace=0, xpage=8.5, ypage=9.5, $
      position=pos, /normal

    medbin = 0.5
    minpts = 100
    minx = 7.75
    
    indx = where((sdssancillary.b_lum gt -900) and (sdssnodust.ebv_hahb_err gt 0.0) and $
      (sdssdust.h_beta_ew_uncor[0] ge ewcut),nindx)

    x = sdssancillary[indx].b_lum
    xabs = sdssancillary[indx].m_b
    xerr = sdssancillary[indx].b_lum_err
    
    xrange = LBrange
    yrange = hbharange; + alog10(HaHb)

    xtitle = 'log L(B) [L(B)'+sunsymbol()+']'
    ytitle = 'log (H\beta/H\alpha)'

; ##########################################################
; Panel 1: No absorption correction, A(Hb)=0, A(Ha)=0
; ##########################################################

    hb = sdssdust[indx].h_beta_uncor[0]
    hb_err = sdssdust[indx].h_beta_uncor[1]
    
    ha = sdssdust[indx].h_alpha[0]
    ha_err = sdssdust[indx].h_alpha[1]
    
    y = alog10(hb/ha); + alog10(HaHb)
    yerr = im_compute_error(hb,hb_err,ha,ha_err,/log)

    stats = im_stats(y,/verbose,/baremin)
    xstr = strtrim(string(stats.median_rej,format='(F12.2)'),2)+' ('+$
      strtrim(string(stats.mean_rej,format='(F12.2)'),2)+$
      '\pm'+strtrim(string(stats.sigma_rej,format='(F12.2)'),2)+')'

    ytitle = 'log (H\beta/H\alpha)_{obs}'

    sdss_lineplot, x, y, xerr, yerr, plottype=3, postscript=postscript, $
      xtitle='', ytitle=ytitle, xrange=xrange, yrange=yrange, legendtype=0, $
      charsize=charsize_5, xstyle=11, xtickname=replicate(' ',10), position=pos[*,0], /silent
    axis, /xaxis, xrange=interpol(xabs,x,!x.crange), xthick=postthick, xstyle=1, $
      charsize=charsize_5, charthick=postthick
    im_xyouts_title, xtitle='M_{B} [mag]', charsize=charsize_5, charthick=postthick
    djs_oplot, !x.crange, alog10(1.0/HaHb)*[1,1], line=2, thick=postthick

    legend, '(a)', /left, /top, box=0, charsize=charsize_3, charthick=postthick
    legend, textoidl(xstr), /right, /top, box=0, charsize=charsize_3, charthick=postthick

    label = ['H\beta Uncorrected for Stellar Absorption']
    legend, textoidl(label), /left, /bottom, box=0, charsize=charsize_3, $
      charthick=postthick

; ############################################################
; Panel 2: Absorption-Corrected, A(Hb)=0, A(Ha)=0
; ############################################################

    hb = sdssdust[indx].h_beta[0]
    hb_err = sdssdust[indx].h_beta[1]
    
    ha = sdssdust[indx].h_alpha[0]
    ha_err = sdssdust[indx].h_alpha[1]
    
    y = alog10(hb/ha); + alog10(HaHb)
    yerr = im_compute_error(hb,hb_err,ha,ha_err,/log)

    stats = im_stats(y,/verbose,/baremin,/no_head)
    xstr = strtrim(string(stats.median_rej,format='(F12.2)'),2)+' ('+$
      strtrim(string(stats.mean_rej,format='(F12.2)'),2)+$
      '\pm'+strtrim(string(stats.sigma_rej,format='(F12.2)'),2)+')'

    ytitle = 'log (H\beta/H\alpha)_{obs}'

    sdss_lineplot, x, y, xerr, yerr, plottype=3, postscript=postscript, $
      xtitle='', ytitle=ytitle, xrange=xrange, yrange=yrange, legendtype=0, $
      charsize=charsize_5, xtickname=replicate(' ',10), position=pos[*,1], /noerase, /silent
    djs_oplot, !x.crange, alog10(1.0/HaHb)*[1,1], line=2, thick=postthick

    legend, '(b)', /left, /top, box=0, charsize=charsize_3, charthick=postthick
    legend, textoidl(xstr), /right, /top, box=0, charsize=charsize_3, charthick=postthick

    label = ['H\beta Corrected for Stellar Absorption']
    legend, textoidl(label), /left, /bottom, box=0, charsize=charsize_3, $
      charthick=postthick
    
; ###############################################################
; Panel 3: Absorption-Corrected, A(Hb)=0, Individual A(Ha)
; ###############################################################

;   hb = sdssdust[indx].h_beta[0]
;   hb_err = sdssdust[indx].h_beta[1]
;   ha = sdssnodust[indx].h_alpha[0]
;   ha_err = sdssnodust[indx].h_alpha[1]
;   y = alog10(hb/ha) - loghaconst - hasfrconstoffset
;   yerr = im_compute_error(hb,hb_err,ha,ha_err,/log)
;   y = hb - ha - loghaconst - hasfrconstoffset
;   yerr = sqrt(ha_err^2 + hb_err^2)

; account for the covariance in the Ha reddening correction

    hb_dust = sdssdust[indx].h_beta[0]
    hb_dust_err = sdssdust[indx].h_beta[1]
    ha_dust = sdssdust[indx].h_alpha[0]
    ha_dust_err = sdssdust[indx].h_alpha[1]
    ha_nodust = sdssnodust[indx].h_alpha[0]

    yratio = hb_dust/ha_nodust
    yratio_err = sqrt(((1-bigbprime_ha)*yratio)^2*((ha_dust_err/ha_dust)^2+(hb_dust_err/hb_dust)^2))

    y = alog10(yratio) - loghaconst - hasfrconstoffset
    yerr = yratio_err/yratio/alog(10.0)
    
    stats = im_stats(y,/verbose,/baremin,/no_head)
    xstr = strtrim(string(stats.median_rej,format='(F12.2)'),2)+' ('+$
      strtrim(string(stats.mean_rej,format='(F12.2)'),2)+$
      '\pm'+strtrim(string(stats.sigma_rej,format='(F12.2)'),2)+')'

    ytitle = 'log [10^{-41} L(H\beta)_{obs}/\psi(H\alpha)]'
    yrange = sfrhbharange - hasfrconstoffset

    sdss_lineplot, x, y, xerr, yerr, plottype=3, postscript=postscript, $
      xtitle=xtitle, ytitle=ytitle, xrange=xrange, yrange=yrange, legendtype=0, $
      charsize=charsize_5, position=pos[*,2], /noerase, /silent
    
    print & running = im_medxbin(x,y,medbin,minx=minx,minpts=minpts,/verbose) & print

    legend, '(c)', /left, /top, box=0, charsize=charsize_3, charthick=postthick
    legend, textoidl(xstr), /right, /top, box=0, charsize=charsize_3, charthick=postthick

    label = ['H\beta Corrected for Stellar Absorption']
    legend, textoidl(label), /left, /bottom, box=0, charsize=charsize_3, $
      charthick=postthick

    im_openclose, postscript=postscript, /close

; ------------------------------------------------------------
; 3-panel - L(B) vs Hb/Ha - Integrated
; ------------------------------------------------------------
    
    psname = 'lb_vs_hbha_3panel'
    im_openclose, pspath+psname, postscript=postscript, xsize=8.5, ysize=9.5, encapsulated=encapsulated

    pagemaker, nx=1, ny=3, height=[2.5,2.5,2.5], width=7.0, xmargin=[1.1,0.4], $
      ymargin=[0.9,1.1], xspace=0, yspace=0, xpage=8.5, ypage=9.5, $
      position=pos, /normal

    medbin = 0.5
    minpts = 10
    minx = 7.75
    
    indx = where((atlasdust.b_lum_obs gt -900) and (atlasnodust.ebv_hahb_err gt 0.0) and $
      (atlasdust.h_beta_ew_uncor[0] ge ewcut),nindx)

    ewstatcor = 4.0*atlasdust[indx].babs_h_beta_continuum[0]
    
    x = atlasdust[indx].b_lum_obs
    xerr = atlasdust[indx].b_lum_obs_err
    xabs = atlasdust[indx].m_b_obs
    
    indxnfgs = where((nfgsnodust.b_lum_obs gt -900.0) and (nfgsnodust.ebv_hahb_err gt 0.0) and $
      (nfgsdust.h_beta_ew_uncor[0] ge ewcut),nindxnfgs)

    ewstatcornfgs = 4.0*nfgsdust[indxnfgs].babs_h_beta_continuum[0]

    xnfgs = nfgsdust[indxnfgs].b_lum_obs
    xerrnfgs = nfgsdust[indxnfgs].b_lum_obs_err
    
    xrange = LBrange
    yrange = hbharange; + alog10(HaHb)

    xtitle = 'log L(B) [L(B)'+sunsymbol()+']'
    ytitle = 'log (H\beta/H\alpha)'

; ##########################################################
; Panel 1: No absorption correction, A(Hb)=0, A(Ha)=0
; ##########################################################

    hb = atlasdust[indx].h_beta_uncor[0]
    hb_err = atlasdust[indx].h_beta_uncor[1]
    
    ha = atlasdust[indx].h_alpha[0]
    ha_err = atlasdust[indx].h_alpha[1]
    
    y = alog10(hb/ha); + alog10(HaHb)
    yerr = im_compute_error(hb,hb_err,ha,ha_err,/log)

    hbnfgs = nfgsdust[indxnfgs].h_beta_uncor[0]
    hbnfgs_err = nfgsdust[indxnfgs].h_beta_uncor[1]
    
    hanfgs = nfgsdust[indxnfgs].h_alpha[0]
    hanfgs_err = nfgsdust[indxnfgs].h_alpha[1]
    
    ynfgs = alog10(hbnfgs/hanfgs); + alog10(HaHb)
    yerrnfgs = im_compute_error(hbnfgs,hbnfgs_err,hanfgs,hanfgs_err,/log)

    stats = im_stats([y,ynfgs],/verbose,/baremin)
    xstr = strtrim(string(stats.median_rej,format='(F12.2)'),2)+' ('+$
      strtrim(string(stats.mean_rej,format='(F12.2)'),2)+$
      '\pm'+strtrim(string(stats.sigma_rej,format='(F12.2)'),2)+')'
    
    ytitle = 'log (H\beta/H\alpha)_{obs}'

    atlas1d_lineplot, x, y, xerr, yerr, plottype=3, postscript=postscript, $
      xtitle='', ytitle=ytitle, xrange=xrange, yrange=yrange, legendtype=0, $
      charsize=charsize_5, xstyle=11, xtickname=replicate(' ',10), position=pos[*,0], $
      xnfgs=xnfgs, ynfgs=ynfgs, xerrnfgs=xerrnfgs, yerrnfgs=yerrnfgs, blackwhite=blackwhite
    axis, /xaxis, xrange=interpol(xabs,x,!x.crange), xthick=postthick, xstyle=1, $
      charsize=charsize_5, charthick=postthick
    im_xyouts_title, xtitle='M_{B} [mag]', charsize=charsize_5, charthick=postthick
    djs_oplot, !x.crange, alog10(1.0/HaHb)*[1,1], line=2, thick=postthick

    legend, '(a)', /left, /top, box=0, charsize=charsize_3, charthick=postthick
    legend, textoidl(xstr), /right, /top, box=0, charsize=charsize_3, charthick=postthick

    label = ['H\beta Uncorrected for Stellar Absorption']
    legend, textoidl(label), /left, /bottom, box=0, charsize=charsize_3, $
      charthick=postthick

; ############################################################
; Panel 2: Absorption-Corrected, A(Hb)=0, A(Ha)=0
; ############################################################

    hb = atlasdust[indx].h_beta[0]
    hb_err = atlasdust[indx].h_beta[1]
    
    ha = atlasdust[indx].h_alpha[0]
    ha_err = atlasdust[indx].h_alpha[1]
    
    y = alog10(hb/ha); + alog10(HaHb)
    yerr = im_compute_error(hb,hb_err,ha,ha_err,/log)

    hbnfgs = nfgsdust[indxnfgs].h_beta[0]
    hbnfgs_err = nfgsdust[indxnfgs].h_beta[1]
    
    hanfgs = nfgsdust[indxnfgs].h_alpha[0]
    hanfgs_err = nfgsdust[indxnfgs].h_alpha[1]
    
    ynfgs = alog10(hbnfgs/hanfgs); + alog10(HaHb)
    yerrnfgs = im_compute_error(hbnfgs,hbnfgs_err,hanfgs,hanfgs_err,/log)

    stats = im_stats([y,ynfgs],/verbose,/baremin,/no_head)
    xstr = strtrim(string(stats.median_rej,format='(F12.2)'),2)+' ('+$
      strtrim(string(stats.mean_rej,format='(F12.2)'),2)+$
      '\pm'+strtrim(string(stats.sigma_rej,format='(F12.2)'),2)+')'

    ytitle = 'log (H\beta/H\alpha)_{obs}'

    atlas1d_lineplot, x, y, xerr, yerr, plottype=3, postscript=postscript, $
      xtitle='', ytitle=ytitle, xrange=xrange, yrange=yrange, legendtype=0, $
      charsize=charsize_5, xtickname=replicate(' ',10), blackwhite=blackwhite, $
      position=pos[*,1], /noerase, xnfgs=xnfgs, ynfgs=ynfgs, xerrnfgs=xerrnfgs, yerrnfgs=yerrnfgs
    djs_oplot, !x.crange, alog10(1.0/HaHb)*[1,1], line=2, thick=postthick

    legend, '(b)', /left, /top, box=0, charsize=charsize_3, charthick=postthick
    legend, textoidl(xstr), /right, /top, box=0, charsize=charsize_3, charthick=postthick

    label = ['H\beta Corrected for Stellar Absorption']
    legend, textoidl(label), /left, /bottom, box=0, charsize=charsize_3, $
      charthick=postthick

; ############################################################
; Panel 3: Absorption-Corrected, A(Hb)=0, Individual A(Ha) 
; ############################################################

; account for the covariance in the Ha reddening correction
    
;   hb = atlasdust[indx].h_beta_uncor[0]+ewstatcor
;   hb = atlasdust[indx].h_beta_uncor[0]

    hb_dust = atlasdust[indx].h_beta[0]
    hb_dust_err = atlasdust[indx].h_beta[1]
    ha_dust = atlasdust[indx].h_alpha[0]
    ha_dust_err = atlasdust[indx].h_alpha[1]
    ha_nodust = atlasnodust[indx].h_alpha[0]

    yratio = hb_dust/ha_nodust
    yratio_err = sqrt(((1-bigbprime_ha)*yratio)^2*((ha_dust_err/ha_dust)^2+(hb_dust_err/hb_dust)^2))

    y = alog10(yratio) - loghaconst - hasfrconstoffset
    yerr = yratio_err/yratio/alog(10.0)
    
;   ha = atlasnodust[indx].h_alpha[0]
;   ha_err = atlasdust[indx].h_alpha[1]
;   y = alog10(hb/ha) - loghaconst - hasfrconstoffset
;   yerr = im_compute_error(hb,hb_err,ha,ha_err,/log)
    
;   hb = atlasdust[indx].h_beta_lum[0]
;   hb_err = atlasdust[indx].h_beta_lum[1]
;   ha = atlasnodust[indx].h_alpha_lum[0]
;   ha_err = atlasnodust[indx].h_alpha_lum[1]
;   ebv_err = 0.4*alog(10.0)*kl_ha*median(atlasnodust[indx].ebv_hahb_err) ; median reddening error
;   ha_err = sqrt(ha_err^2 + ebv_err^2)
;   y = hb - ha - loghaconst - hasfrconstoffset
;   yerr = sqrt(ha_err^2 + hb_err^2)

;;  hbnfgs = alog10((nfgsdust[indxnfgs].h_beta_uncor[0]+ewstatcornfgs)*4*!pi*((nfgsdust[indxnfgs].distance*3.086D24)^2)/lsun)
;;  hbnfgs = alog10(nfgsdust[indxnfgs].h_beta_uncor[0]*4*!pi*((nfgsdust[indxnfgs].distance*3.086D24)^2)/lsun)
;   hbnfgs = nfgsdust[indxnfgs].h_beta_lum[0]
;   hbnfgs_err = nfgsdust[indxnfgs].h_beta_lum[1]
;   hanfgs = nfgsnodust[indxnfgs].h_alpha_lum[0]
;   hanfgs_err = nfgsnodust[indxnfgs].h_alpha_lum[1]
;   ynfgs = hbnfgs - hanfgs - loghaconst - hasfrconstoffset
;   yerrnfgs = sqrt(hanfgs_err^2 + hbnfgs_err^2)
       
;   hb = nfgsdust[indx].h_beta_uncor[0]+ewstatcor
;   hb = nfgsdust[indx].h_beta_uncor[0]
;   hbnfgs = nfgsdust[indxnfgs].h_beta[0]
;   hbnfgs_err = nfgsdust[indxnfgs].h_beta[1]

;   hanfgs = nfgsnodust[indxnfgs].h_alpha[0]
;   hanfgs_err = nfgsnodust[indxnfgs].h_alpha[1]
;   ebv_err = 0.4*alog(10.0)*kl_ha*median(nfgsnodust[indxnfgs].ebv_hahb_err)*ha ; median reddening error
;   hanfgs_err = sqrt(hanfgs_err^2 + ebv_err^2)
;   ynfgs = alog10(hbnfgs/hanfgs) - loghaconst - hasfrconstoffset
;   yerrnfgs = im_compute_error(hbnfgs,hbnfgs_err,hanfgs,hanfgs_err,/log)

    hbnfgs_dust = nfgsdust[indxnfgs].h_beta[0]
    hbnfgs_dust_err = nfgsdust[indxnfgs].h_beta[1]
    hanfgs_dust = nfgsdust[indxnfgs].h_alpha[0]
    hanfgs_dust_err = nfgsdust[indxnfgs].h_alpha[1]
    hanfgs_nodust = nfgsnodust[indxnfgs].h_alpha[0]

    yrationfgs = hbnfgs_dust/hanfgs_nodust
    yrationfgs_err = sqrt(((1-bigbprime_ha)*yrationfgs)^2*((hanfgs_dust_err/hanfgs_dust)^2+(hbnfgs_dust_err/hbnfgs_dust)^2))

    ynfgs = alog10(yrationfgs) - loghaconst - hasfrconstoffset
    yerrnfgs = yrationfgs_err/yrationfgs/alog(10.0)
    
    stats = im_stats([y,ynfgs],/verbose,/baremin,/no_head)
    xstr = strtrim(string(stats.median_rej,format='(F12.2)'),2)+' ('+$
      strtrim(string(stats.mean_rej,format='(F12.2)'),2)+$
      '\pm'+strtrim(string(stats.sigma_rej,format='(F12.2)'),2)+')'
        
    ytitle = 'log [10^{-41} L(H\beta)_{obs}/\psi(H\alpha)]'
    yrange = sfrhbharange - hasfrconstoffset

    atlas1d_lineplot, x, y, xerr, yerr, plottype=3, postscript=postscript, $
      xtitle=xtitle, ytitle=ytitle, xrange=xrange, yrange=yrange, legendtype=0, $
      charsize=charsize_5, position=pos[*,2], /noerase, blackwhite=blackwhite, $
      xnfgs=xnfgs, ynfgs=ynfgs, xerrnfgs=xerrnfgs, yerrnfgs=yerrnfgs

    print & running = im_medxbin([x,xnfgs],[y,ynfgs],medbin,minx=minx,minpts=minpts,/verbose) & print
;   oploterror, running.binctr, running.medy, running.sigy, ps=-4, xsty=3, ysty=3, yrange=yrange

    legend, '(c)', /left, /top, box=0, charsize=charsize_3, charthick=postthick
    legend, textoidl(xstr), /right, /top, box=0, charsize=charsize_3, charthick=postthick

    label = ['H\beta Corrected for Stellar Absorption']
    legend, textoidl(label), /left, /bottom, box=0, charsize=charsize_3, $
      charthick=postthick

;   w = where([x,xnfgs] gt 9.5)
;   junk = im_stats(([y,ynfgs])[w],/verbose)

    im_openclose, postscript=postscript, /close

; ------------------------------------------------------------
; 2-panel - L(B) vs SFR(Ha)/SFR(IR)
; ------------------------------------------------------------
    
    psname = 'lb_vs_sfr_hair_2panel'
    im_openclose, pspath+psname, postscript=postscript, xsize=8.5, ysize=8.0, encapsulated=encapsulated

    pagemaker, nx=1, ny=2, height=[3.0,3.0], width=7.0, xmargin=[1.1,0.4], $
      ymargin=[0.9,1.1], xspace=0, yspace=0, xpage=8.5, ypage=8.0, $
      position=pos, /normal

;   im_openclose, pspath+psname, postscript=postscript, xsize=7.0, ysize=8.0, encapsulated=encapsulated
;   pagemaker, nx=1, ny=2, height=[3.0,3.0], width=5.5, xmargin=[1.1,0.4], $
;     ymargin=[0.9,1.1], xspace=0, yspace=0, xpage=7.0, ypage=8.0, $
;     position=pos, /normal

    if keyword_set(blackwhite) then begin
       atlascolor = 'dark gray'
       nfgscolor = 'black'
    endif else begin
       atlascolor = 'blue'
       nfgscolor = 'red'
    endelse
       
    indx = where((atlasnodust.b_lum_obs gt -900) and (atlasnodust.ir_flux gt -900) and (atlasnodust.ebv_hahb_err gt 0.0),nindx)

    x = atlasdust[indx].b_lum_obs
    xerr = atlasdust[indx].b_lum_obs_err
    xabs = atlasdust[indx].m_b_obs
    
    lir = atlasdust[indx].ir_flux ; [erg/s/cm2]
    lir_err = sqrt(atlasdust[indx].ir_flux_err^2 + (atlasdust[indx].ir_flux*0.15)^2)

    indxnfgs = where((nfgsdust.b_lum_obs gt -900) and (nfgsdust.ir_flux gt -900) and (nfgsnodust.ebv_hahb_err gt 0.0),nindxnfgs)

    xnfgs = nfgsnodust[indxnfgs].b_lum_obs
    xerrnfgs = nfgsnodust[indxnfgs].b_lum_obs_err
    
    lirnfgs = nfgsdust[indxnfgs].ir_flux
    lirnfgs_err = sqrt(nfgsdust[indxnfgs].ir_flux_err^2 + (nfgsdust[indxnfgs].ir_flux*0.15)^2)
    
; flag the LIRGS and ULIRGS

    lirgs = where(atlasdust[indx].ir_lum gt 11.0,comp=notlirgs)
    lirgsnfgs = where(nfgsdust[indxnfgs].ir_lum gt 11.0,notlirgsnfgs)

; plotting variables    
    
    xrange = lbrange
    yrange = sfrHaIRrange

    xtitle = 'log L(B) [L(B)'+sunsymbol()+']'
    ytitle = 'log \psi(H\alpha)/\psi(IR)'

    lhalir = alog10(4.5D-44/7.9D-42)

; ##########################################################
; Panel 1: No reddening correction
; ##########################################################

    ha = atlasdust[indx].h_alpha[0]
    ha_err = atlasdust[indx].h_alpha[1]
    
    y = alog10(ha/lir) - lhalir
    yerr = im_compute_error(ha,ha_err,lir,lir_err,/log)

    hanfgs = nfgsdust[indxnfgs].h_alpha[0]
    hanfgs_err = nfgsdust[indxnfgs].h_alpha[1]
    
    ynfgs = alog10(hanfgs/lirnfgs) - lhalir
    yerrnfgs = im_compute_error(hanfgs,hanfgs_err,lirnfgs,lirnfgs_err,/log)

    stats = im_stats([y,ynfgs],/verbose)
    xstr = strtrim(string(stats.median_rej,format='(F12.2)'),2)+' ('+$
      strtrim(string(stats.mean_rej,format='(F12.2)'),2)+$
      '\pm'+strtrim(string(stats.sigma_rej,format='(F12.2)'),2)+')'

    ytitle = 'log \psi(H\alpha)/\psi(IR)'

    atlas1d_lineplot, x, y, xerr, yerr, plottype=3, postscript=postscript, $
      xtitle='', ytitle=ytitle, xrange=xrange, yrange=yrange, legendtype=0, $
      charsize=charsize_7, position=pos[*,0], xnfgs=xnfgs, ynfgs=ynfgs, $
      xerrnfgs=xerrnfgs, yerrnfgs=yerrnfgs, xtickname=replicate(' ',10), $
      yfrac=8, xstyle=11, blackwhite=blackwhite
    axis, /xaxis, xrange=interpol(xabs,x,!x.crange), xthick=postthick, xstyle=1, $
      charsize=charsize_7, charthick=postthick
    im_xyouts_title, xtitle='M_{B} [mag]', charsize=charsize_7, charthick=postthick

    djs_oplot, !x.crange, [0,0], line=0, thick=postthick

    legend, '(a) '+textoidl('Observed H\alpha'), /left, /top, box=0, charsize=charsize_4, charthick=postthick
;   legend, '(a)', /left, /top, box=0, charsize=charsize_4, charthick=postthick
;   legend, textoidl('Observed H\alpha'), /left, /bottom, box=0, charsize=charsize_4, charthick=postthick
    legend, textoidl(xstr), /right, /top, box=0, charsize=charsize_4, charthick=postthick

    w = where([x,xnfgs] ge 9.0 and [x,xnfgs] le 11.0)
;   junk = im_stats(([y,ynfgs])[w],/verbose)

; ############################################################
; Panel 2: Good Reddening Correction
; ############################################################

; account for the covariance in the Ha reddening correction
    
;   ha = atlasnodust[indx].h_alpha[0]
;   ha_err = atlasnodust[indx].h_alpha[1]

    hb_dust = atlasdust[indx].h_beta[0]
    hb_dust_err = atlasdust[indx].h_beta[1]
    ha_dust = atlasdust[indx].h_alpha[0]
    ha_dust_err = atlasdust[indx].h_alpha[1]
    ha_nodust = atlasnodust[indx].h_alpha[0] ; ha_dust*(2.86*hb_dust/ha_dust)^bigbprime_ha
    ha_nodust_err = sqrt((ha_nodust*(1-bigbprime_ha)*ha_dust_err/ha_dust)^2 + (bigbprime_ha*ha_nodust*hb_dust_err/hb_dust)^2)
    
    y = alog10(ha_nodust/lir) - lhalir
    yerr = im_compute_error(ha_nodust,ha_nodust_err,lir,lir_err,/log)
    
;   hanfgs = nfgsnodust[indxnfgs].h_alpha[0]
;   hanfgs_err = nfgsnodust[indxnfgs].h_alpha[1]

    hbnfgs_dust = nfgsdust[indxnfgs].h_beta[0]
    hbnfgs_dust_err = nfgsdust[indxnfgs].h_beta[1]
    hanfgs_dust = nfgsdust[indxnfgs].h_alpha[0]
    hanfgs_dust_err = nfgsdust[indxnfgs].h_alpha[1]
    hanfgs_nodust = hanfgs_dust*(2.86*hbnfgs_dust/hanfgs_dust)^bigbprime_ha
    hanfgs_nodust_err = sqrt((hanfgs_nodust*(1-bigbprime_ha)*hanfgs_dust_err/hanfgs_dust)^2 + $
      (bigbprime_ha*hanfgs_nodust*hbnfgs_dust_err/hbnfgs_dust)^2)
    
    ynfgs = alog10(hanfgs_nodust/lirnfgs) - lhalir
    yerrnfgs = im_compute_error(hanfgs_nodust,hanfgs_nodust_err,lirnfgs,lirnfgs_err,/log)

    stats = im_stats([y,ynfgs],/verbose,/no_head)
    xstr = strtrim(string(stats.median_rej,format='(F12.2)'),2)+' ('+$
      strtrim(string(stats.mean_rej,format='(F12.2)'),2)+$
      '\pm'+strtrim(string(stats.sigma_rej,format='(F12.2)'),2)+')'

    ytitle = 'log \psi(H\alpha)/\psi(IR)'

    atlas1d_lineplot, x, y, xerr, yerr, plottype=3, postscript=postscript, $
      xtitle=xtitle, ytitle=ytitle, xrange=xrange, yrange=yrange, legendtype=0, $
      charsize=charsize_7, position=pos[*,1], /noerase, xnfgs=xnfgs, ynfgs=ynfgs, $
      xerrnfgs=xerrnfgs, yerrnfgs=yerrnfgs, blackwhite=blackwhite
    djs_oplot, !x.crange, [0,0], line=0, thick=postthick

    legend, '(b) '+textoidl('Extinction-Corrected H\alpha'), /left, /top, box=0, charsize=charsize_4, charthick=postthick
;   legend, '(b)', /left, /top, box=0, charsize=charsize_4, charthick=postthick
;   legend, textoidl('Extinction-Corrected H\alpha'), /left, /bottom, box=0, charsize=charsize_4, charthick=postthick
    legend, textoidl(xstr), /right, /top, box=0, charsize=charsize_4, charthick=postthick

    w = where([x,xnfgs] ge 9.0 and [x,xnfgs] le 11.0)
;   junk = im_stats(([y,ynfgs])[w],/verbose)

; overplot the prediction by Hirashita et al. (2003):

    eta = 0.4
    epsilon = 0.5
    frac = 0.6

    hbi_lir = 10^atlasdust[indx].ir_lum
    srt = sort(hbi_lir)
    hbi_lir = hbi_lir[srt]
    hbi_lha = lsun*10^atlasdust[indx].h_alpha_lum[0]
    hbi_lha = hbi_lha[srt]
    
    ir_sfr = 1.79E-10*(1-eta) / (0.13 - 0.085*frac + 0.87*epsilon) * hbi_lir
    ha_sfr = 7.89D-42 / frac * hbi_lha

;   oplot, alog10(ir_sfr), alog10(ha_sfr/ir_sfr), line=2, thick=postthick

; overplot the 60/100 micron ratios

    xratio = alog10(atlasdust[indx].iras_60/atlasdust[indx].iras_100)
    xnfgsratio = alog10(nfgsdust[indxnfgs].iras_60/nfgsdust[indxnfgs].iras_100)
    
    jj = im_stats(xratio,/verbose)
    jj = im_stats(xnfgsratio,/verbose)

    binsize = 0.05

    xbig = [xratio,xnfgsratio]
    plothist, xbig, bin=binsize, xbigbin, ybigbin, /noplot, /halfbin
    yrange3 = minmax(ybigbin)*[1.0,1.1]

; upper panel    

;   histpos = [0.22,0.58,0.42,0.7]
;   xpolypos = [0.22,0.42,0.42,0.22]
;   ypolypos = [0.58,0.58,0.7,0.7]
    
; lower panel    

    histpos = [0.22,0.21,0.42,0.33]
    xpolypos = [0.22,0.42,0.42,0.22]
    ypolypos = [0.21,0.21,0.33,0.33]

    if keyword_set(postscript) then polyfill, xpolypos, ypolypos, /normal, color=djs_icolor('white')
    djs_plot, [0], [0], /nodata, xthick=postthick, ythick=postthick, charsize=1.0, $
      charthick=postthick3, xtitle=textoidl('log [S_{\nu}(60)/S_{\nu}(100)]'), $
      ytitle='Number', ystyle=1, xrange=[-1,0.3], yrange=yrange3, xstyle=1, $
      position=histpos, /noerase, xtickinterval=0.5, ytickinterval=10.0
    plothist, xratio, bin=binsize, /overplot, color=fsc_color(atlascolor,!d.table_size-3), $
      thick=postthick3, /halfbin, line=0
    plothist, xnfgsratio, bin=binsize, thick=1.0, line=0, /halfbin, $
      /overplot, /fill, /fline, forientation=45, fcolor=fsc_color(nfgscolor,!d.table_size-4), $
      color=fsc_color(nfgscolor,!d.table_size-4), fspacing=0.05

    im_openclose, postscript=postscript, /close

; ------------------------------------------------------------
; SFR(Ha) vs various SFR(Hb) calibrations
; ------------------------------------------------------------
    
    psname = 'sfr_ha_vs_sfr_hb_2panel'
    im_openclose, pspath+psname, postscript=postscript, xsize=8.5, ysize=7.5, encapsulated=encapsulated

    pagemaker, nx=1, ny=2, height=[3.0,3.0], width=7.0, xmargin=[1.2,0.3], $
      ymargin=[0.4,1.1], xspace=0, yspace=0, xpage=8.5, ypage=7.5, $
      position=pos, /normal

    if keyword_set(blackwhite) then begin
       sfcolor = 'gray'
       agncolor = 'black'
    endif else begin
       sfcolor = 'sky blue'
       agncolor = 'navy'
    endelse
       
    xrange = sfrharange
    yrange = residrange_hbsfrs

    xtitle = 'log \psi(H\alpha) ['+sfr_units()+']'
    ytitle = 'log \psi(H\beta)/\psi(H\alpha)'
    
; ##########################################################
; Panel 1
; ##########################################################

; SF
    
    indx = where((sfrs.sfr_ha gt -900) and (sfrs.sfr_hb_hbhg gt -900),nindx)

    x = sfrs[indx].sfr_ha
    xerr = sfrs[indx].sfr_ha_err
    
    y1 = sfrs[indx].sfr_ha
    y1err = sfrs[indx].sfr_ha_err
    
    y2 = sfrs[indx].sfr_hb_hbhg
    y2err = sfrs[indx].sfr_hb_hbhg_err

    resid = y2-y1
    residerr = sqrt(y1err^2+y2err^2)

; AGN
    
    indx_agn = where((sfrs_agn.sfr_ha gt -900) and (sfrs_agn.sfr_hb_hbhg gt -900),nindx_agn)

    x_agn = sfrs_agn[indx_agn].sfr_ha
    xerr_agn = sfrs_agn[indx_agn].sfr_ha_err
    
    y1_agn = sfrs_agn[indx_agn].sfr_ha
    y1err_agn = sfrs_agn[indx_agn].sfr_ha_err
    
    y2_agn = sfrs_agn[indx_agn].sfr_hb_hbhg
    y2err_agn = sfrs_agn[indx_agn].sfr_hb_hbhg_err

    resid_agn = y2_agn-y1_agn
    residerr_agn = sqrt(y1err_agn^2+y2err_agn^2)

; make the plot    

    stats = im_stats(resid,/verbose,/baremin)
;   stats = im_stats([resid,resid_agn],/verbose,/baremin)
    xstr = strtrim(string(stats.median_rej,format='(F12.2)'),2)+' ('+$
      strtrim(string(stats.mean_rej,format='(F12.2)'),2)+$
      '\pm'+strtrim(string(stats.sigma_rej,format='(F12.2)'),2)+')'

; SF    
    
    atlas1d_lineplot, x, resid, xerr, residerr, plottype=1, postscript=postscript, $
      xtitle='', ytitle=ytitle, xrange=xrange, yrange=yrange, legendtype=0, blackwhite=blackwhite, $
      charsize=charsize_8, xtickname=replicate(' ',10), position=pos[*,0], $
      atlaspsize=0.8, ytickinterval=1.0, atlascolor=sfcolor
    djs_oplot, !x.crange, [0,0], line=0, thick=postthick

; AGN    
    
    atlas1d_lineplot, x_agn, resid_agn, xerr_agn, residerr_agn, blackwhite=blackwhite, $
      postscript=postscript, atlaspsize=0.8, atlasfill=0, /overplot, atlascolor=agncolor, $
      coloroffset=10

    legend, '(a)', /left, /top, box=0, charsize=charsize_4, $
      charthick=postthick, clear=keyword_set(postscript)
    legend, textoidl('Extinction-Corrected using H\beta/H\gamma'), /left, /bottom, box=0, charsize=charsize_4, $
      charthick=postthick, clear=keyword_set(postscript)
    legend, textoidl(xstr), /right, /top, box=0, charsize=charsize_4, charthick=postthick

; ##########################################################
; Panel 2
; ##########################################################

; SF

    indx = where((sfrs.sfr_ha gt -900) and (sfrs.sfr_hb_best gt -900),nindx)

    x = sfrs[indx].sfr_ha
    xerr = sfrs[indx].sfr_ha_err
    
    y1 = sfrs[indx].sfr_ha
    y1err = sfrs[indx].sfr_ha_err
    
    y2 = sfrs[indx].sfr_hb_best
    y2err = sfrs[indx].sfr_hb_best_err

    resid = y2-y1
    resid_err = sqrt(y1err^2+y2err^2)

; AGN

    indx_agn = where((sfrs_agn.sfr_ha gt -900) and (sfrs_agn.sfr_hb_best gt -900),nindx_agn)

    x_agn = sfrs_agn[indx_agn].sfr_ha
    xerr_agn = sfrs_agn[indx_agn].sfr_ha_err
    
    y1_agn = sfrs_agn[indx_agn].sfr_ha
    y1err_agn = sfrs_agn[indx_agn].sfr_ha_err
    
    y2_agn = sfrs_agn[indx_agn].sfr_hb_best
    y2err_agn = sfrs_agn[indx_agn].sfr_hb_best_err

    resid_agn = y2_agn-y1_agn
    residerr_agn = sqrt(y1err_agn^2+y2err_agn^2)

; make the plot    

    stats = im_stats(resid,/verbose,/baremin)
;   stats = im_stats([resid,resid_agn],/verbose,/baremin)
    xstr = strtrim(string(stats.median_rej,format='(F12.2)'),2)+' ('+$
      strtrim(string(stats.mean_rej,format='(F12.2)'),2)+$
      '\pm'+strtrim(string(stats.sigma_rej,format='(F12.2)'),2)+')'

; SF
    
    atlas1d_lineplot, x, resid, xerr, residerr, plottype=1, postscript=postscript, $
      xtitle=xtitle, ytitle=ytitle, xrange=xrange, yrange=yrange, legendtype=0, $
      charsize=charsize_8, position=pos[*,1], blackwhite=blackwhite, $
      /noerase, atlaspsize=0.8, ytickinterval=1.0, atlascolor=sfcolor
    djs_oplot, !x.crange, [0,0], line=0, thick=postthick

; AGN    
    
    atlas1d_lineplot, x_agn, resid_agn, xerr_agn, residerr_agn, blackwhite=blackwhite, $
      postscript=postscript, atlaspsize=0.8, atlasfill=0, /overplot, atlascolor=agncolor, $
      coloroffset=10

    legend, '(b)', /left, /top, box=0, charsize=charsize_4, $
      charthick=postthick, clear=keyword_set(postscript)
    legend, textoidl('L(B) vs \psi/L_{obs} Median Relation'), /left, /bottom, box=0, charsize=charsize_4, $
      charthick=postthick, clear=keyword_set(postscript)
    legend, textoidl(xstr), /right, /top, box=0, charsize=charsize_4, charthick=postthick

    im_openclose, postscript=postscript, /close    

; ------------------------------------------------------------
; L(B) vs E(Hb-Ha) [Integrated + SDSS]
; ------------------------------------------------------------

    psname = 'lb_vs_ehbha_2panel'
    im_openclose, pspath+psname, postscript=postscript, xsize=8.5, ysize=5.15, encapsulated=encapsulated

    pagemaker, nx=2, ny=1, height=3.15, width=3.15, xmargin=[1.1,1.1], $
      ymargin=[0.9,1.1], xspace=0, yspace=0, xpage=8.5, ypage=5.15, $
      position=pos, /normal

; Integrated    
    
    indx = where((atlasdust.b_lum_obs gt -900.0) and (atlasnodust.ehbha_err gt 0.0),nindx)

    x = atlasdust[indx].b_lum_obs
    xerr = atlasdust[indx].b_lum_obs_err
    xabs = atlasdust[indx].m_b_obs
    
    y = atlasnodust[indx].ehbha
    yerr = atlasnodust[indx].ehbha_err
    yebv = y/(k_lambda(4861,/odonnell)-k_lambda(6563,/odonnell))

    indxnfgs = where((nfgsdust.b_lum_obs gt -900) and (nfgsnodust.ehbha_err gt 0.0),nindxnfgs)

    xnfgs = nfgsdust[indxnfgs].b_lum_obs
    xerrnfgs = nfgsdust[indxnfgs].b_lum_obs_err
    
    ynfgs = nfgsnodust[indxnfgs].ehbha
    yerrnfgs = nfgsnodust[indxnfgs].ehbha_err

    xtitle = 'log L(B) [L(B)'+sunsymbol()+']'
    ytitle = 'E(H\beta-H\alpha) [mag]'

    xrange = LBrange
    yrange = ehbharange
    
    atlas1d_lineplot, x, y, xerr, yerr, plottype=3, postscript=postscript, ysty=3, $
      ytitle=ytitle, xtitle=xtitle, xrange=xrange, yrange=yrange, legendtype=0, $
      xstyle=11, /right, /top, xnfgs=xnfgs, ynfgs=ynfgs, xerrnfgs=xerrnfgs, blackwhite=blackwhite, $
      yerrnfgs=yerrnfgs, position=pos[*,0], charsize=charsize_6
    axis, /xaxis, xrange=interpol(xabs,x,!x.crange), xthick=postthick, xstyle=1, $
      charsize=charsize_6, charthick=postthick
    im_xyouts_title, xtitle='M_{B} [mag]', charsize=charsize_6, charthick=postthick

; SDSS    
    
    indx = where((sdssancillary.b_lum gt -900) and (sdssnodust.ehbha_err gt 0.0),nindx)

    x = sdssancillary[indx].b_lum
    xabs = sdssancillary[indx].m_b
    xerr = sdssancillary[indx].b_lum_err
    
    y = sdssnodust[indx].ehbha
    yerr = sdssnodust[indx].ehbha_err
    yebv = y/(k_lambda(4861,/odonnell)-k_lambda(6563,/odonnell))

    sdss_lineplot, x, y, xerr, yerr, plottype=3, postscript=postscript, $
      ytitle='', xtitle=xtitle, xrange=xrange, yrange=yrange, legendtype=0, $
      xstyle=11, ystyle=11, /right, /top, position=pos[*,1], charsize=charsize_6, /noerase, $
      ytickname=replicate(' ',10)
    axis, /yaxis, yrange=interpol(yebv,y,!y.crange), ythick=postthick, ystyle=1, $
      charsize=charsize_6, charthick=postthick
    axis, /xaxis, xrange=interpol(xabs,x,!x.crange), xthick=postthick, xstyle=1, $
      charsize=charsize_6, charthick=postthick
    im_xyouts_title, ytitle='E(B-V) [mag]', charsize=charsize_6, charthick=postthick
    im_xyouts_title, xtitle='M_{B} [mag]', charsize=charsize_6, charthick=postthick

    im_openclose, postscript=postscript, /close

; ------------------------------------------------------------
; E(Hb-Ha) vs [O II]/Ha [Integrated + SDSS]
; ------------------------------------------------------------

    psname = 'ehbha_vs_oiiha_2panel'
    im_openclose, pspath+psname, postscript=postscript, xsize=8.5, ysize=5.5, encapsulated=encapsulated

    pagemaker, nx=2, ny=1, height=3.5, width=3.5, xmargin=[1.2,0.3], $
      ymargin=[0.9,1.1], xspace=0, yspace=0, xpage=8.5, ypage=5.5, $
      position=pos, /normal

    medbin = 0.05
    minpts = 10
    minx = 0.025

; Integrated    
    
    indx = where((atlasdust.oii_3727[0]/atlasdust.oii_3727[1] gt snrcut) and $
      (atlasnodust.ehbha_err gt 0.0),nindx)

    x = atlasnodust[indx].ehbha
    xerr = atlasnodust[indx].ehbha_err
    xebv = x/(k_lambda(4861,/odonnell)-k_lambda(6563,/odonnell))

    indxnfgs = where((nfgsdust.oii_3727[0]/nfgsdust.oii_3727[1] gt snrcut) and $
      (nfgsnodust.ehbha_err gt 0.0),nindxnfgs)

    xnfgs = nfgsnodust[indxnfgs].ehbha
    xerrnfgs = nfgsnodust[indxnfgs].ehbha_err

    xtitle = 'E(H\beta-H\alpha) [mag]'
    ytitle = 'log ([O II]/H\alpha)_{obs}'

    xrange = ehbharange
    yrange = oiihacorrange2
    
    y1 = atlasdust[indx].oii_3727[0]
    y1err = atlasdust[indx].oii_3727[1]

    y2 = atlasdust[indx].h_alpha[0]
    y2err = atlasdust[indx].h_alpha[1]

    y = alog10(y1/y2)
    yerr = im_compute_error(y1,y1err,y2,y2err,/log)

    y1nfgs = nfgsdust[indxnfgs].oii_3727[0]
    y1errnfgs = nfgsdust[indxnfgs].oii_3727[1]

    y2nfgs = nfgsdust[indxnfgs].h_alpha[0]
    y2errnfgs = nfgsdust[indxnfgs].h_alpha[1]

    ynfgs = alog10(y1nfgs/y2nfgs)
    yerrnfgs = im_compute_error(y1nfgs,y1errnfgs,y2nfgs,y2errnfgs,/log)

    stats = im_stats([y,ynfgs],/verbose,/baremin)
    rcor = r_correlate([x,xnfgs],[y,ynfgs],zd=zd,probd=probd)
;   splog, 'Spearman rank test: ', rcor, probd, zd

    atlas1d_lineplot, x, y, xerr, yerr, plottype=3, postscript=postscript, $
      ytitle=ytitle, xtitle=xtitle, xrange=xrange, yrange=yrange, legendtype=0, $
      xstyle=11, /right, /top, xnfgs=xnfgs, ynfgs=ynfgs, xerrnfgs=xerrnfgs, blackwhite=blackwhite, $
      yerrnfgs=yerrnfgs, position=pos[*,0], charsize=charsize_6
    axis, /xaxis, xrange=interpol(xebv,x,!x.crange), xthick=postthick, xstyle=1, $
      charsize=charsize_6, charthick=postthick
    im_xyouts_title, xtitle='E(B-V) [mag]', charsize=charsize_6, charthick=postthick

    w = where(([x,xnfgs] gt 0.35) and (x lt 0.45))
    junk = im_stats(([y,ynfgs])[w],sigrej=3.0)
    splog, 'Total range: ', 10^(junk.maxrej-junk.minrej)

; mark I Zw 18

;   lowz = speclinefit_locate(atlasnodust[indx],'UGCA166')
;   xyouts, x[lowz]+0.18, y[lowz]-0.03, 'I Zw 18', align=0.5, charsize=charsize_2, charthick=postthick
;   plotsym, 6, 2.5, thick=postthick
;   plots, x[lowz]+0.15, y[lowz], psym=8
    
;   legend, '(a)', /left, /top, box=0, charsize=charsize_2, charthick=postthick
;   print & running = im_medxbin([x,xnfgs],[y,ynfgs],medbin,minx=minx,minpts=minpts,/verbose) & print

; overplot reddening relations

    ratio = 0.0

    x_ebv = findgen((1.3-(-0.1))/0.1)*0.1+(-0.1)
    yodonnell =  -0.4*x_ebv*(k_lambda(3727.0,/odonnell)-k_lambda(6563.0,/odonnell)) + ratio
    ycharlot = -0.4*x_ebv*(k_lambda(3727,/charlot)-k_lambda(6563,/charlot)) + ratio
;   ysmc = -0.4*x_ebv*(k_lambda(3727,/smc)-k_lambda(6563,/smc)) + ratio

    djs_oplot, x_ebv, yodonnell, line=0, thick=postthick
    djs_oplot, x_ebv, ycharlot, line=2, thick=postthick
;   djs_oplot, x_ebv, ysmc, line=2, thick=postthick

    label = ["O'Donnell (1994)",'Charlot & Fall (2000)']
    linestyle = [0,2]
    legend, label, /right, /top, box=0, linestyle=linestyle, $
      charsize=charsize_3, charthick=postthick, thick=postthick

; SDSS    
    
    indx = where((sdssdust.oii_3727[0]/sdssdust.oii_3727[1] gt snrcut) and $
      (sdssnodust.ehbha_err gt 0.0),nindx)

    x = sdssnodust[indx].ehbha
    xerr = sdssnodust[indx].ehbha_err
    xebv = x/(k_lambda(4861,/odonnell)-k_lambda(6563,/odonnell))

    y1 = sdssdust[indx].oii_3727[0]
    y1err = sdssdust[indx].oii_3727[1]

    y2 = sdssdust[indx].h_alpha[0]
    y2err = sdssdust[indx].h_alpha[1]

    y = alog10(y1/y2)
    yerr = im_compute_error(y1,y1err,y2,y2err,/log)

    stats = im_stats(y,/verbose,/baremin,/no_head)
    rcor = r_correlate(x,y,zd=zd,probd=probd)
;   splog, 'Spearman rank test: ', rcor, probd, zd

    sdss_lineplot, x, y, xerr, yerr, plottype=3, postscript=postscript, $
      ytitle='', xtitle=xtitle, xrange=xrange, yrange=yrange, legendtype=0, $
      xstyle=11, /right, /top, position=pos[*,1], charsize=charsize_6, /noerase, $
      ytickname=replicate(' ',10)
    axis, /xaxis, xrange=interpol(xebv,x,!x.crange), xthick=postthick, xstyle=1, $
      charsize=charsize_6, charthick=postthick
    im_xyouts_title, xtitle='E(B-V) [mag]', charsize=charsize_6, charthick=postthick

;   legend, '(b)', /left, /top, box=0, charsize=charsize_2, charthick=postthick
;   print & running = im_medxbin(x,y,medbin,minx=minx,minpts=minpts,/verbose) & print

    w = where((x gt 0.35) and (x lt 0.45))
    junk = im_stats(y[w],sigrej=3.0)
    splog, 'Total range: ', 10^(junk.maxrej-junk.minrej)

; overplot reddening relations

    ratio = -0.16

    x_ebv = findgen((1.3-(-0.1))/0.1)*0.1+(-0.1)
    yodonnell =  -0.4*x_ebv*(k_lambda(3727.0,/odonnell)-k_lambda(6563.0,/odonnell)) + ratio
    ycharlot = -0.4*x_ebv*(k_lambda(3727,/charlot)-k_lambda(6563,/charlot)) + ratio
;   ysmc = -0.4*x_ebv*(k_lambda(3727,/smc)-k_lambda(6563,/smc)) + ratio

    djs_oplot, x_ebv, yodonnell, line=0, thick=postthick
    djs_oplot, x_ebv, ycharlot, line=2, thick=postthick
;   djs_oplot, x_ebv, ysmc, line=2, thick=postthick

    label = ["O'Donnell (1994)",'Charlot & Fall (2000)']
    linestyle = [0,2]
;   legend, label, /right, /top, box=0, linestyle=linestyle, $
;     charsize=charsize_2, charthick=postthick, thick=postthick

    im_openclose, postscript=postscript, /close

; ------------------------------------------------------------
; 4-panel histogram plots - ATLAS/NFGS
; ------------------------------------------------------------

    psname = 'histogram_properties'
    im_openclose, pspath+psname, postscript=postscript, xsize=8.5, ysize=9.25, encapsulated=encapsulated

    pagemaker, nx=2, ny=2, xspace=1.0, yspace=1.0, width=3.1*[1,1], height=3.1*[1,1], $
      xmargin=[1.0,0.3], ymargin=[0.8,1.1], xpage=8.5, ypage=9.25, position=pos, /normal

    if keyword_set(blackwhite) then begin
       atlascolor = ''
       nfgscolor = 'grey'
    endif else begin
       atlascolor = 'blue'
       nfgscolor = 'red'
    endelse

; #######
; Panel 1    
; #######
    
    indx = where(atlasdust.b_lum_obs gt -900.0)
    x = atlasdust[indx].b_lum_obs
    xabs = atlasdust[indx].m_b_obs

    indxnfgs = where(nfgsdust.b_lum_obs gt -900.0)
    xnfgs = nfgsdust[indxnfgs].b_lum_obs
    xabsnfgs = nfgsdust[indxnfgs].m_b_obs

    xbig = [x,xnfgs]
    xabsbig = [xabs,xabsnfgs]

    stats = im_stats(xbig,/verbose)
    stats = im_stats(xabsbig,/verbose)

    xtitle = 'log L(B) [L(B)'+sunsymbol()+']'
    binsize = 0.25

    plothist, xbig, bin=binsize, xbigbin, ybigbin, /noplot, /halfbin
    plothist, x, bin=binsize, xbin, ybin, /noplot, /halfbin
    plothist, xnfgs, bin=binsize, xbinnfgs, ybinnfgs, /noplot, /halfbin

    xrange = LBrange
    yrange = minmax(ybigbin)*[1.0,1.1]

    djs_plot, [0], [0], /nodata, xthick=postthick, ythick=postthick, charsize=charsize_5, $
      charthick=postthick, thick=postthick, xtitle=textoidl(xtitle), ytitle='Number', $
      xstyle=11, ystyle=1, xrange=xrange, yrange=yrange, position=pos[*,0]
    axis, /xaxis, xrange=interpol(xabs,x,!x.crange), xthick=postthick, xstyle=1, $
      charsize=charsize_5, charthick=postthick
    im_xyouts_title, xtitle='M_{B} [mag]', charsize=charsize_5, charthick=postthick

    plothist, x, bin=binsize, /overplot, color=djs_icolor(atlascolor), thick=postthick, /halfbin
    plothist, xnfgs, bin=binsize, thick=postthick, line=0, /halfbin, $
      /overplot, /fill, /fline, forientation=45, fcolor=djs_icolor(nfgscolor), $
      color=djs_icolor(nfgscolor), fspacing=0.05
    plothist, xbig, bin=binsize, thick=postthick, line=2, /overplot, /halfbin

    legend, '(a)', /left, /top, charthick=postthick, charsize=charsize_3, box=0, clear=postscript

; overlay the local galaxy luminosity function

    h = 0.7                     ; Hubble constant
    Msun = 5.42                 ; Solar absolute B magnitude
    
    phistar = 1.6E-2*h^3        ; [/Mpc3]
    Mstar = -19.7 + 5*alog10(h) ; [mag]
    alpha = -1.07
    Lstar = 10^(0.4*Msun)*10^(-0.4*Mstar)

    M = reverse(findgen(((-16.5)-(-22.5))/0.1+1)*0.1+(-22.5))
;   Lum = 10^(findgen((!x.crange[1]-!x.crange[0])/0.01+1)*0.01+!x.crange[0])
    Lum = 10^(0.4*Msun)*10^(-0.4*M)
    
    phi = (phistar/Lstar)*(Lum/Lstar)^alpha*exp(-Lum/Lstar)
;   norm = interpol(ybin,xbin,10.5)/interpol(phi,alog10(Lum),10.5)

    norm = int_tabulated(10^xbin,float(ybin))/int_tabulated(Lum,phi)

;   oplot, alog10(lum), phi*norm, line=0

; #######
; Panel 2
; #######

    kha = k_lambda(6563,/odonnel)    
    
    indx = where(atlasnodust.ebv_hahb_err gt 0)
    xebv = atlasnodust[indx].ebv_hahb
    x = atlasnodust[indx].ebv_hahb*kha

    indxnfgs = where(nfgsnodust.ebv_hahb_err gt 0.0)
    xnfgs = nfgsnodust[indxnfgs].ebv_hahb*kha
    xbig = [x,xnfgs]

    stats = im_stats(xbig,/verbose,/no_head)

    xtitle = 'A(H\alpha) [mag]'

    binsize = 0.15
    
    plothist, xbig, bin=binsize, xbigbin, ybigbin, /noplot, /halfbin
    plothist, x, bin=binsize, xbin, ybin, /noplot, /halfbin
    plothist, xnfgs, bin=binsize, xbinnfgs, ybinnfgs, /noplot, /halfbin

    xrange = AHarange
    yrange = minmax(ybigbin)*[1.0,1.25]

    djs_plot, [0], [0], /nodata, xthick=postthick, ythick=postthick, charsize=charsize_5, $
      charthick=postthick, thick=postthick, xtitle=textoidl(xtitle), ytitle='Number', $
      xstyle=11, ystyle=1, xrange=xrange, yrange=yrange, position=pos[*,1], /noerase
    axis, /xaxis, xrange=interpol(xebv,x,!x.crange), xthick=postthick, xstyle=1, $
      charsize=charsize_5, charthick=postthick
    im_xyouts_title, xtitle='E(B-V) [mag]', charsize=charsize_5, charthick=postthick

    plothist, x, bin=binsize, /overplot, color=djs_icolor(atlascolor), thick=postthick, /halfbin
    plothist, xnfgs, bin=binsize, thick=postthick, line=0, /halfbin, $
      /overplot, /fill, /fline, forientation=45, fcolor=djs_icolor(nfgscolor), $
      color=djs_icolor(nfgscolor), fspacing=0.05
    plothist, xbig, bin=binsize, thick=postthick, line=2, /overplot, /halfbin

    legend, '(b)', /left, /top, charthick=postthick, charsize=charsize_3, box=0, clear=postscript, spacing=0

; #######
; Panel 3
; #######

; in this panel, use the OIIINII metallicity calibration, and the NIIHA
; calibration for objects beyond the range of the OIIINII calibration
    
    indx = where(atlasdust.zstrong_12oh_oiiinii_niiha gt -900,nindx)
    x = atlasdust[indx].zstrong_12oh_oiiinii_niiha

    xsun = x-Zsun_new

    indxnfgs = where(nfgsdust.zstrong_12oh_oiiinii_niiha gt -900)
    xnfgs = nfgsdust[indxnfgs].zstrong_12oh_oiiinii_niiha

    xbig = [x,xnfgs]
    stats = im_stats(xbig,/verbose,/no_head)

    xtitle = '12 + log (O/H)'

    binsize = 0.1
    
    plothist, xbig, bin=binsize, xbigbin, ybigbin, /noplot, /halfbin
    plothist, x, bin=binsize, xbin, ybin, /noplot, /halfbin
    plothist, xnfgs, bin=binsize, xbinnfgs, ybinnfgs, /noplot, /halfbin

    xrange = ohrange
    yrange = minmax(ybigbin)*[1.0,1.1]

    djs_plot, [0], [0], /nodata, xthick=postthick, ythick=postthick, charsize=charsize_5, $
      charthick=postthick, thick=postthick, xtitle=textoidl(xtitle), ytitle='Number', $
      xstyle=1, ystyle=1, xrange=xrange, yrange=yrange, position=pos[*,2], /noerase

    plothist, x, bin=binsize, /overplot, color=djs_icolor(atlascolor), thick=postthick, /halfbin
    plothist, xnfgs, bin=binsize, thick=postthick, line=0, /halfbin, $
      /overplot, /fill, /fline, forientation=45, fcolor=djs_icolor(nfgscolor), $
      color=djs_icolor(nfgscolor), fspacing=0.05
    plothist, xbig, bin=binsize, thick=postthick, line=2, /overplot, /halfbin

    legend, '(c)', /left, /top, charthick=postthick, charsize=charsize_3, box=0, clear=postscript, spacing=0

; #######
; Panel 4
; #######

    lineratio, atlasdust, 'OIII_5007','OII_3727', '', '', $
      x, xerr, dum1, dum2, index=indx, nindex=nindx, snrcut=snrcut, $
      xsyserr=syserr, ysyserr=syserr
    
    lineratio, nfgsdust, 'OIII_5007','OII_3727', '', '', $
      xnfgs, xerrnfgs, dum1, dum2, index=indxnfgs, nindex=nindxnfgs, snrcut=snrcut, $
      xsyserr=nfgssyserr, ysyserr=nfgssyserr

    xbig = [x,xnfgs]
    splog, 'Atlas, NFGS, Combined:'
    stats = im_stats(x,/verbose,/no_head)
    stats = im_stats(xnfgs,/verbose,/no_head)
    stats = im_stats(xbig,/verbose,/no_head)

    xtitle = 'log ([O III]/[O II])_{obs}'

    binsize = 0.1
    
    plothist, xbig, bin=binsize, xbigbin, ybigbin, /noplot, /halfbin
    plothist, x, bin=binsize, xbin, ybin, /noplot, /halfbin
    plothist, xnfgs, bin=binsize, xbinnfgs, ybinnfgs, /noplot, /halfbin

    xrange = oiiioiirange
    yrange = minmax(ybigbin)*[1.0,1.1]

    djs_plot, [0], [0], /nodata, xthick=postthick, ythick=postthick, charsize=charsize_5, $
      charthick=postthick, thick=postthick, xtitle=textoidl(xtitle), ytitle='Number', $
      xstyle=1, ystyle=1, xrange=xrange, yrange=yrange, position=pos[*,3], /noerase
    plothist, x, bin=binsize, /overplot, color=djs_icolor(atlascolor), thick=postthick, /halfbin

    plothist, xnfgs, bin=binsize, thick=postthick, line=0, /halfbin, $
      /overplot, /fill, /fline, forientation=45, fcolor=djs_icolor(nfgscolor), $
      color=djs_icolor(nfgscolor), fspacing=0.05
    plothist, xbig, bin=binsize, thick=postthick, line=2, /overplot, /halfbin

    legend, '(d)', /left, /top, charthick=postthick, charsize=charsize_3, box=0, clear=postscript, spacing=0

    im_openclose, postscript=postscript, /close    

; ------------------------------------------------------------
; 4-panel histogram plots - SDSS
; ------------------------------------------------------------

    psname = 'sdss_histogram_properties'
    im_openclose, pspath+psname, postscript=postscript, xsize=8.5, ysize=9.25, encapsulated=encapsulated

    pagemaker, nx=2, ny=2, xspace=1.0, yspace=1.0, width=3.1*[1,1], height=3.1*[1,1], $
      xmargin=[1.0,0.3], ymargin=[0.8,1.1], xpage=8.5, ypage=9.25, position=pos, /normal

; #######
; Panel 1    
; #######
    
    indx = where(sdssancillary.b_lum gt -900.0)
    x = sdssancillary[indx].b_lum
    xabs = sdssancillary[indx].m_b

    indxatlas = where(atlasdust.b_lum_obs gt -900.0)
    xatlas = atlasdust[indxatlas].b_lum_obs

    indxnfgs = where(nfgsdust.b_lum_obs gt -900.0)
    xnfgs = nfgsdust[indxnfgs].b_lum_obs
    xbig = [xatlas,xnfgs]

    stats = im_stats(x,/verbose)
    stats = im_stats(xabs,/verbose)
    
    xtitle = 'log L(B) [L(B)'+sunsymbol()+']'
    xrange = LBrange
    binsize = 0.25

    im_plothist, x, bin=binsize, xbin, ybin, /noplot, /fraction, /halfbin
    yrange = minmax(ybin)*[1.0,1.1]

    im_plothist, xbig, bin=binsize, xbigbin, ybigbin, /noplot, /fraction, /halfbin
    yrange = [yrange[0],(yrange[1]>max(ybigbin))*1.1]

    djs_plot, [0], [0], /nodata, xthick=postthick, ythick=postthick, charsize=charsize_5, $
      charthick=postthick, thick=postthick, xtitle=textoidl(xtitle), ytitle='Fraction', $
      xstyle=11, ystyle=1, xrange=xrange, yrange=yrange, position=pos[*,0]
    axis, /xaxis, xrange=interpol(xabs,x,!x.crange), xthick=postthick, xstyle=1, $
      charsize=charsize_5, charthick=postthick
    im_xyouts_title, xtitle='M_{B} [mag]', charsize=charsize_5, charthick=postthick

    im_plothist, x, bin=binsize, /overplot, thick=postthick, /fraction, /halfbin; color=djs_icolor('dark green'), 
    im_plothist, xbig, bin=binsize, thick=postthick, line=2, /overplot, /fraction, /halfbin

    legend, '(a)', /left, /top, charthick=postthick, charsize=charsize_3, box=0, clear=postscript, spacing=0

; #######
; Panel 2
; ####### 
    
    kha = k_lambda(6563,/odonnel)    

    indx = where(sdssnodust.ebv_hahb_err gt 0.0)
    xebv = sdssnodust[indx].ebv_hahb
    x = sdssnodust[indx].ebv_hahb*kha

    indxatlas = where(atlasnodust.ebv_hahb_err gt 0.0)
    xatlas = atlasnodust[indxatlas].ebv_hahb*kha

    indxnfgs = where(nfgsnodust.ebv_hahb_err gt 0.0)
    xnfgs = nfgsnodust[indxnfgs].ebv_hahb*kha
    xbig = [xatlas,xnfgs]

    stats = im_stats(x,/verbose,/no_head)

    xtitle = 'A(H\alpha) [mag]'
    xrange = Aharange ; [-0.2,2.5]
    binsize = 0.15
    
    im_plothist, x, bin=binsize, xbin, ybin, /noplot, /fraction, /halfbin
    yrange = minmax(ybin)*[1.0,1.15]

    im_plothist, xbig, bin=binsize, xbigbin, ybigbin, /noplot, /fraction, /halfbin
    yrange = [yrange[0],(yrange[1]>max(ybigbin))*1.1]

    djs_plot, [0], [0], /nodata, xthick=postthick, ythick=postthick, charsize=charsize_5, $
      charthick=postthick, thick=postthick, xtitle=textoidl(xtitle), ytitle='Fraction', $
      xstyle=11, ystyle=1, xrange=xrange, yrange=yrange, position=pos[*,1], /noerase
    axis, /xaxis, xrange=interpol(xebv,x,!x.crange), xthick=postthick, xstyle=1, $
      charsize=charsize_5, charthick=postthick
    im_xyouts_title, xtitle='E(B-V) [mag]', charsize=charsize_5, charthick=postthick

    im_plothist, x, bin=binsize, /overplot, thick=postthick, /fraction, /halfbin; color=djs_icolor('dark green'), 
    im_plothist, xbig, bin=binsize, thick=postthick, line=2, /overplot, /fraction, /halfbin

    legend, '(b)', /left, /top, charthick=postthick, charsize=charsize_3, box=0, clear=postscript, spacing=0

; #######
; Panel 3
; #######

    indx = where(sdssdust.zstrong_12oh_oiiinii_niiha gt -900,nindx)
    x = sdssdust[indx].zstrong_12oh_oiiinii_niiha
    
    xsun = x-Zsun_new
    
    indxatlas = where(atlasdust.zstrong_12oh_oiiinii_niiha gt -900,nindxatlas)
    xatlas = atlasdust[indxatlas].zstrong_12oh_oiiinii_niiha

    indxnfgs = where(nfgsdust.zstrong_12oh_oiiinii_niiha gt -900)
    xnfgs = nfgsdust[indxnfgs].zstrong_12oh_oiiinii_niiha
    xbig = [xatlas,xnfgs]

    stats = im_stats(x,/verbose,/no_head)

    xtitle = '12 + log (O/H)'

    xrange = ohrange
    binsize = 0.1
    
    im_plothist, x, bin=binsize, xbin, ybin, /noplot, /fraction, /halfbin
    yrange = minmax(ybin)*[1.0,1.1]

    im_plothist, xbig, bin=binsize, xbigbin, ybigbin, /noplot, /fraction, /halfbin
    yrange = [yrange[0],(yrange[1]>max(ybigbin))*1.1]
    
    djs_plot, [0], [0], /nodata, xthick=postthick, ythick=postthick, charsize=charsize_5, $
      charthick=postthick, thick=postthick, xtitle=textoidl(xtitle), ytitle='Fraction', $
      xstyle=1, ystyle=1, xrange=xrange, yrange=yrange, position=pos[*,2], /noerase

    im_plothist, x, bin=binsize, /overplot, thick=postthick, /fraction, /halfbin; color=djs_icolor('dark green'), 
    im_plothist, xbig, bin=binsize, thick=postthick, line=2, /overplot, /fraction, /halfbin

    legend, '(c)', /left, /top, charthick=postthick, charsize=charsize_3, box=0, clear=postscript, spacing=0

; #######
; Panel 4
; #######

    lineratio, sdssdust, 'OIII_5007','OII_3727', '', '', $
      x, xerr, dum1, dum2, index=indx, nindex=nindx, snrcut=snrcut
    
    lineratio, atlasdust, 'OIII_5007','OII_3727', '', '', $
      xatlas, xerratlas, dum1, dum2, index=indxatlas, nindex=nindxatlas, snrcut=snrcut, $
      xsyserr=syserr, ysyserr=syserr

    lineratio, nfgsdust, 'OIII_5007','OII_3727', '', '', $
      xnfgs, xerrnfgs, dum1, dum2, index=indxnfgs, nindex=nindxnfgs, snrcut=snrcut, $
      xsyserr=nfgssyserr, ysyserr=nfgssyserr
    xbig = [xatlas,xnfgs]

    stats = im_stats(x,/verbose,/no_head)

    xtitle = 'log ([O III]/[O II])_{obs}'
    xrange = oiiioiirange
    binsize = 0.1
    
    im_plothist, x, bin=binsize, xbin, ybin, /noplot, /fraction, /halfbin
    yrange = minmax(ybin)*[1.0,1.1]

    im_plothist, xbig, bin=binsize, xbigbin, ybigbin, /noplot, /fraction, /halfbin
    yrange = [yrange[0],(yrange[1]>max(ybigbin))*1.1]

    djs_plot, [0], [0], /nodata, xthick=postthick, ythick=postthick, charsize=charsize_5, $
      charthick=postthick, thick=postthick, xtitle=textoidl(xtitle), ytitle='Fraction', $
      xstyle=1, ystyle=1, xrange=xrange, yrange=yrange, position=pos[*,3], /noerase

    im_plothist, x, bin=binsize, /overplot, thick=postthick, /fraction, /halfbin; color=djs_icolor('dark green'), 
    im_plothist, xbig, bin=binsize, thick=postthick, line=2, /overplot, /fraction, /halfbin

    legend, '(d)', /left, /top, charthick=postthick, charsize=charsize_3, box=0, clear=postscript, spacing=0

    im_openclose, postscript=postscript, /close    

; ------------------------------------------------------------
; [N II]/Ha versus [O III]/Hb
; ------------------------------------------------------------
    
    psname = 'niiha_vs_oiiihb'
    im_openclose, pspath+psname, postscript=postscript, xsize=8.5, ysize=8.0, encapsulated=encapsulated

    pagemaker, nx=1, ny=1, xspace=0, yspace=0.0, width=6.6, height=6.6, $
      xmargin=[1.5,0.4], ymargin=[0.3,1.1], xpage=8.5, ypage=8.0, $
      position=pos, /normal

; SF galaxies    
    
    lineratio, atlasdust, 'NII_6584', 'H_ALPHA', 'OIII_5007', 'H_BETA', $
      x, xerr, y, yerr, index=indx, nindex=nindx, snrcut=1.0, $ ; NOTE S/N cut!
      xsyserr=syserr, ysyserr=syserr

    lineratio, nfgsdust, 'NII_6584', 'H_ALPHA', 'OIII_5007', 'H_BETA', $
      xnfgs, xerrnfgs, ynfgs, yerrnfgs, index=indxnfgs, snrcut=1.0, $ ; NOTE S/N cut!
      xsyserr=nfgssyserr, ysyserr=nfgssyserr

; AGN
    
    lineratio, atlasdust_agn, 'NII_6584', 'H_ALPHA', 'OIII_5007', 'H_BETA', $
      x_agn, xerr_agn, y_agn, yerr_agn, index=indx, nindex=nindx_agn, $
      snrcut=1.0, xsyserr=syserr, ysyserr=syserr

    lineratio, nfgsdust_agn, 'NII_6584', 'H_ALPHA', 'OIII_5007', 'H_BETA', $
      xnfgs_agn, xerrnfgs_agn, ynfgs_agn, yerrnfgs_agn, index=indxnfgs_agn, $
      snrcut=1.0, xsyserr=nfgssyserr, ysyserr=nfgssyserr

    xtitle = 'log ([N II] \lambda6584/H\alpha)_{obs}'
    ytitle = 'log ([O III] \lambda5007/H\beta)_{obs}'

    xrange = niiharange
    yrange = oiiihbrange

    good = where((hii.nii_6584_h_alpha gt -900.0) and (hii.oiii_5007_h_beta gt -900.0))
    xregion = hii[good].nii_6584_h_alpha & xerrregion = hii[good].nii_6584_h_alpha_err
    yregion = hii[good].oiii_5007_h_beta & yerrregion = hii[good].oiii_5007_h_beta_err

    atlas1d_lineplot, x, y, xerr, yerr, plottype=3, postscript=postscript, $
      xtitle=xtitle, ytitle=ytitle, xrange=xrange, yrange=yrange, legendtype=0, $
      /right, /top, position=pos[*,0], atlasfill=1, nfgsfill=1, /errorleft, $ 
      xregion=xregion, yregion=yregion, xerrregion=xerrregion, yerrregion=yerrregion, $
      xnfgs=xnfgs, ynfgs=ynfgs, xerrnfgs=xerrnfgs, yerrnfgs=yerrnfgs, blackwhite=blackwhite

; overplot the AGN

    atlas1d_lineplot, x_agn, y_agn, xerr_agn, yerr_agn, $
      postscript=postscript, /overplot, atlasfill=0, $
      nfgsfill=0, xnfgs=xnfgs_agn, ynfgs=ynfgs_agn, $ 
      xerrnfgs=xerrnfgs_agn, yerrnfgs=yerrnfgs_agn, blackwhite=blackwhite

; overplot L. Kewley's starburst mixing line

    models = kewley_bpt_lines(/kauffmann,_extra=extra)
    oplot, models.x_nii, models.y_nii, line=0, thick=postthick

    models = kewley_bpt_lines(_extra=extra)
    oplot, models.x_nii, models.y_nii, line=2, thick=postthick

    im_openclose, postscript=postscript, /close

; ------------------------------------------------------------
; High-Redshift L(B) vs [O III]/[O II]
; ------------------------------------------------------------

    psname = 'highz_lb_vs_oiiioii'
    im_openclose, pspath+psname, postscript=postscript, xsize=8.5, ysize=8.0, encapsulated=encapsulated

    pagemaker, nx=1, ny=1, xspace=0, yspace=0.0, width=6.0, height=6.0, $
      xmargin=[1.3,1.2], ymargin=[0.9,1.1], xpage=8.5, ypage=8.0, $
      position=pos, /normal

; SF galaxies    

    cut = where(atlasnodust.b_lum_obs gt -900.0)
    lineratio, atlasdust[cut], '', '', 'OIII_5007', 'OII_3727', $
      dum1, dum2, y, yerr, index=indx, nindex=nindx, snrcut=snrcut, $
      xsyserr=syserr, ysyserr=syserr

    x = atlasnodust[cut[indx]].b_lum_obs
    xerr = atlasnodust[cut[indx]].b_lum_obs_err
    xabs = atlasnodust[cut[indx]].m_b_obs

    cutnfgs = where(nfgsnodust.b_lum_obs gt -900.0)
    lineratio, nfgsdust[cutnfgs], '', '', 'OIII_5007', 'OII_3727', $
      dum1, dum2, ynfgs, yerrnfgs, index=indxnfgs, snrcut=snrcut, $
      xsyserr=nfgssyserr, ysyserr=nfgssyserr
    xnfgs = nfgsnodust[cutnfgs[indxnfgs]].b_lum_obs
    xerrnfgs = nfgsnodust[cutnfgs[indxnfgs]].b_lum_obs_err

    xbig = [x,xnfgs]
    ybig = [y,ynfgs]
    xerrbig = [xerr,xerrnfgs]
    yerrbig = [yerr,yerrnfgs]
       
; AGN

    cut_agn = where(atlasnodust_agn.b_lum_obs gt -900.0)
    lineratio, atlasdust_agn[cut_agn], '', '', 'OIII_5007', 'OII_3727', $
      dum1, dum2, y_agn, yerr_agn, index=indx_agn, nindex=nindx_agn, snrcut=snrcut, $
      xsyserr=syserr, ysyserr=syserr

    x_agn = atlasnodust_agn[cut_agn[indx_agn]].b_lum_obs
    xerr_agn = atlasnodust_agn[cut_agn[indx_agn]].b_lum_obs_err
    xabs_agn = atlasnodust_agn[cut_agn[indx_agn]].m_b_obs

    cutnfgs_agn = where(nfgsnodust_agn.b_lum_obs gt -900.0)
    lineratio, nfgsdust_agn[cutnfgs_agn], '', '', 'OIII_5007', 'OII_3727', $
      dum1, dum2, ynfgs_agn, yerrnfgs_agn, index=indxnfgs_agn, snrcut=snrcut, $
      xsyserr=nfgssyserr, ysyserr=nfgssyserr
    xnfgs_agn = nfgsnodust_agn[cutnfgs_agn[indxnfgs_agn]].b_lum_obs
    xerrnfgs_agn = nfgsnodust_agn[cutnfgs_agn[indxnfgs_agn]].b_lum_obs_err

    xbig_agn = [x_agn,xnfgs_agn]
    ybig_agn = [y_agn,ynfgs_agn]
    xerrbig_agn = [xerr_agn,xerrnfgs_agn]
    yerrbig_agn = [yerr_agn,yerrnfgs_agn]

; average error bar
    
    avgerrlocal = djs_median([yerrbig,yerrbig_agn])

; ####################
; Maier et al. (2005)
; ####################
    
    cutm05 = where(m05.b_lum gt -900.0)
    lineratio, m05[cutm05], '', '', 'OIII_5007', 'OII_3727', dum1, dum2, $
      ym05, yerrm05, index=indxm05, nindex=nindxm05, snrcut=snrcut_highz

    xm05 = m05[cutm05[indxm05]].b_lum
    xerrm05 = m05[cutm05[indxm05]].b_lum_err

    avgerrm05 = djs_median(yerrm05)
    splog, 'Maier: redshift = ', minmax(m05[cutm05[indxm05]].z_obj), median(m05[cutm05[indxm05]].z_obj)

; ####################
; Savaglio et al. (2005)
; ####################
    
    cutsava05 = where(sava05.b_lum gt -900.0)
    lineratio, sava05[cutsava05], '', '', 'OIII_5007', 'OII_3727', dum1, dum2, $
      ysava05, yerrsava05, index=indxsava05, nindex=nindxsava05, snrcut=snrcut_highz

    xsava05 = sava05[cutsava05[indxsava05]].b_lum
    xerrsava05 = sava05[cutsava05[indxsava05]].b_lum_err

    avgerrsava05 = djs_median(yerrsava05)
    splog, 'Savaglio: redshift = ', minmax(sava05[cutsava05[indxsava05]].z_obj), median(sava05[cutsava05[indxsava05]].z_obj)

; ####################
; Liang et al. (2004)
; ####################
    
    cutliang04 = where(liang04.b_lum gt -900.0)
    lineratio, liang04[cutliang04], '', '', 'OIII_5007', 'OII_3727', dum1, dum2, $
      yliang04, yerrliang04, index=indxliang04, nindex=nindxliang04, snrcut=snrcut_highz

    xliang04 = liang04[cutliang04[indxliang04]].b_lum
    xerrliang04 = liang04[cutliang04[indxliang04]].b_lum_err

    avgerrliang04 = djs_median(yerrliang04)
    splog, 'Liang: redshift = ', minmax(liang04[cutliang04[indxliang04]].z_obj), median(liang04[cutliang04[indxliang04]].z_obj)

; ####################
; Kobulnicky et al. (2003, 2004)
; ####################
    
    cutk04 = where(k04.b_lum gt -900.0)
    lineratio, k04[cutk04], '', '', 'OIII_5007_EW', 'OII_3727_EW', dum1, dum2, $
      yk04, yerrk04, index=indxk04, nindex=nindxk04, snrcut=snrcut_highz

    xk04 = k04[cutk04[indxk04]].b_lum
    xerrk04 = k04[cutk04[indxk04]].b_lum_err

    avgerrk04 = djs_median(yerrk04)
    splog, 'Kobulnicky: redshift = ', minmax(k04[cutk04[indxk04]].z_obj), median(k04[cutk04[indxk04]].z_obj)

; ####################
; Lilly et al. (2003)
; ####################
    
    cutlilly03 = where(lilly03.b_lum gt -900.0)
    lineratio, lilly03[cutlilly03], '', '', 'OIII_5007', 'OII_3727', dum1, dum2, $
      ylilly03, yerrlilly03, index=indxlilly03, nindex=nindxlilly03, snrcut=snrcut_highz

    xlilly03 = lilly03[cutlilly03[indxlilly03]].b_lum
    xerrlilly03 = lilly03[cutlilly03[indxlilly03]].b_lum_err

    avgerrlilly03 = djs_median(yerrlilly03)
    splog, 'Lilly: redshift = ', minmax(lilly03[cutlilly03[indxlilly03]].z_obj), median(lilly03[cutlilly03[indxlilly03]].z_obj)

; ####################

    xtitle = 'log L(B) [L(B)'+sunsymbol()+']'
    ytitle = 'log ([O III] \lambda5007/[O II])_{obs}'

    xrange = LBrange
    yrange = oiiioiirange2

    djs_plot, [0], [0], /nodata, xthick=postthick, ythick=postthick, charsize=charsize_8, $
      charthick=postthick, thick=postthick, xtitle=xtitle, ytitle=ytitle, $
      ymargin=[4,3], xstyle=11, ystyle=11, xrange=xrange, yrange=yrange, position=pos[*,0]
    axis, /yaxis, yrange=interpol(U,sbgrids[3,*].oiii_5007_oii,!y.crange), ythick=postthick, ystyle=1, $
      charsize=charsize_8, charthick=postthick
    axis, /xaxis, xrange=interpol(xabs,x,!x.crange), xthick=postthick, xstyle=1, $
      charsize=charsize_8, charthick=postthick
    im_xyouts_title, xtitle='M_{B} [mag]', charsize=charsize_8, charthick=postthick
    im_xyouts_title, ytitle='log U', charsize=charsize_8, charthick=postthick, xspacing=11.0

    im_symbols, localsym, psize=localpsize, /fill, color=fsc_color(localcolor,!d.table_size-151)
    djs_oplot, xbig, ybig, ps=8, color=fsc_color(localcolor,!d.table_size-151)

    im_symbols, localsym, psize=localpsize, fill=0, color=fsc_color(localcolor,!d.table_size-151), thick=symthick
    djs_oplot, xbig_agn, ybig_agn, ps=8, color=fsc_color(localcolor,!d.table_size-151)

    im_symbols, sava05sym, psize=sava05psize, /fill, color=fsc_color(sava05color,!d.table_size-153), thick=postthick
    djs_oplot, xsava05, ysava05, ps=8, color=fsc_color(sava05color,!d.table_size-153), thick=postthick

    im_symbols, liang04sym, psize=liang04psize, /fill, color=fsc_color(liang04color,!d.table_size-154), thick=postthick
    djs_oplot, xliang04, yliang04, ps=8, color=fsc_color(liang04color,!d.table_size-154), thick=postthick

    im_symbols, lilly03sym, psize=lilly03psize, /fill, color=fsc_color(lilly03color,!d.table_size-155), thick=postthick
    djs_oplot, xlilly03, ylilly03, ps=8, color=fsc_color(lilly03color,!d.table_size-155), thick=postthick

; keep maier on top    
    
    im_symbols, m05sym, psize=m05psize, /fill, color=fsc_color(m05color,!d.table_size-152), thick=postthick
    djs_oplot, xm05, ym05, ps=8, color=fsc_color(m05color,!d.table_size-152), thick=postthick

; overplot the average error bars for the high-z points and the local 
; sample

    ypos = 1.3
    
    im_symbols, localsym, psize=localpsize, /fill, color=fsc_color(localcolor,!d.table_size-151)
    oploterror, 7.2, ypos, 0.0, avgerrlocal, /data, thick=postthick, $
      errthick=postthick, psym=8, errcolor=fsc_color(localcolor,!d.table_size-151)
       
    im_symbols, m05sym, psize=m05psize, /fill, color=fsc_color(m05color,!d.table_size-152), thick=postthick
    oploterror, 7.5, ypos, 0.0, avgerrm05, /data, thick=postthick, $
      errthick=postthick, psym=8, errcolor=fsc_color(m05color,!d.table_size-152)
    
    im_symbols, sava05sym, psize=sava05psize, /fill, color=fsc_color(sava05color,!d.table_size-153), thick=postthick
    oploterror, 7.8, ypos, 0.0, avgerrsava05, /data, thick=postthick, $
      errthick=postthick, psym=8, errcolor=fsc_color(sava05color,!d.table_size-153)
    
    im_symbols, liang04sym, psize=liang04psize, /fill, color=fsc_color(liang04color,!d.table_size-154), thick=postthick
    oploterror, 8.1, ypos, 0.0, avgerrliang04, /data, thick=postthick, $
      errthick=postthick, psym=8, errcolor=fsc_color(liang04color,!d.table_size-154)
    
    im_symbols, lilly03sym, psize=lilly03psize, /fill, color=fsc_color(lilly03color,!d.table_size-155), thick=postthick
    oploterror, 8.4, ypos, 0.0, avgerrlilly03, /data, thick=postthick, $
      errthick=postthick, psym=8, errcolor=fsc_color(lilly03color,!d.table_size-155)

; overplot a typical reddening vector

    ebv = 0.3
    lpos = 11.0
    kl_oiii_oii = k_lambda(5007,/odonnell)-k_lambda(3727,/odonnell)

    oiii_oii_true = -1.7
    oiii_oii_red = oiii_oii_true - 0.4*kl_oiii_oii*ebv
    
    arrow, lpos, oiii_oii_true, lpos, oiii_oii_red, /data, $
      thick=postthick, hsize=-0.5, hthick=postthick
    xyouts, lpos, oiii_oii_true-0.1, 'E(B-V) = '+string(ebv,format='(F3.1)'), $
      charsize=charsize_1, align=0.5, charthick=postthick

; legend    
    
    label = ['Local Sample - Star-Forming','Local Sample - AGN',$
      'Maier et al. (2005)','Savaglio et al. (2005)',$
      'Liang et al. (2004)','Lilly et al. (2003)']
    psym = [localsym,localsym,m05sym,sava05sym,liang04sym,lilly03sym]
    fill = [1,0,1,1,1,1]
    color = fsc_color([localcolor,localcolor,m05color,sava05color,liang04color,lilly03color],!d.table_size-[151,151,152,153,154,155])

    postthick1 = postthick
    im_legend, textoidl(label), /left, /bottom, box=0, charsize=charsize_2, $
      charthick=postthick1, psym=psym, fill=fill, symsize=1.3, $
      spacing=1.8, thick=postthick1, textcolor=djs_icolor(replicate('',6))
    
    im_openclose, postscript=postscript, /close    

; ------------------------------------------------------------
; High-Redshift L(B) vs [O II]/Ha
; ------------------------------------------------------------

    psname = 'highz_lb_vs_oiiha'
    im_openclose, pspath+psname, postscript=postscript, xsize=8.5, ysize=8.8, encapsulated=encapsulated

    pagemaker, nx=1, ny=1, xspace=0, yspace=0.0, width=6.8, height=6.8, $
      xmargin=[1.4,0.3], ymargin=[0.9,1.1], xpage=8.5, ypage=8.8, $
      position=pos, /normal

; SF galaxies    

    cut = where(atlasnodust.b_lum_obs gt -900.0)
    lineratio, atlasdust[cut], '', '', 'OII_3727', 'H_ALPHA', $
      dum1, dum2, y, yerr, index=indx, nindex=nindx, snrcut=snrcut, $
      xsyserr=syserr, ysyserr=syserr

    x = atlasnodust[cut[indx]].b_lum_obs
    xerr = atlasnodust[cut[indx]].b_lum_obs_err
    xabs = atlasnodust[cut[indx]].m_b_obs

    cutnfgs = where(nfgsnodust.b_lum_obs gt -900.0)
    lineratio, nfgsdust[cutnfgs], '', '', 'OII_3727', 'H_ALPHA', $
      dum1, dum2, ynfgs, yerrnfgs, index=indxnfgs, snrcut=snrcut, $
      xsyserr=nfgssyserr, ysyserr=nfgssyserr
    xnfgs = nfgsnodust[cutnfgs[indxnfgs]].b_lum_obs
    xerrnfgs = nfgsnodust[cutnfgs[indxnfgs]].b_lum_obs_err

    xbig = [x,xnfgs]
    ybig = [y,ynfgs]
    xerrbig = [xerr,xerrnfgs]
    yerrbig = [yerr,yerrnfgs]

; AGN

    cut_agn = where(atlasnodust_agn.b_lum_obs gt -900.0)
    lineratio, atlasdust_agn[cut_agn], '', '', 'OII_3727', 'H_ALPHA', $
      dum1, dum2, y_agn, yerr_agn, index=indx_agn, nindex=nindx_agn, snrcut=snrcut, $
      xsyserr=syserr, ysyserr=syserr

    x_agn = atlasnodust_agn[cut_agn[indx_agn]].b_lum_obs
    xerr_agn = atlasnodust_agn[cut_agn[indx_agn]].b_lum_obs_err
    xabs_agn = atlasnodust_agn[cut_agn[indx_agn]].m_b_obs

    cutnfgs_agn = where(nfgsnodust_agn.b_lum_obs gt -900.0)
    lineratio, nfgsdust_agn[cutnfgs_agn], '', '', 'OII_3727', 'H_ALPHA', $
      dum1, dum2, ynfgs_agn, yerrnfgs_agn, index=indxnfgs_agn, snrcut=snrcut, $
      xsyserr=nfgssyserr, ysyserr=nfgssyserr
    xnfgs_agn = nfgsnodust_agn[cutnfgs_agn[indxnfgs_agn]].b_lum_obs
    xerrnfgs_agn = nfgsnodust_agn[cutnfgs_agn[indxnfgs_agn]].b_lum_obs_err

    xbig_agn = [x_agn,xnfgs_agn]
    ybig_agn = [y_agn,ynfgs_agn]
    xerrbig_agn = [xerr_agn,xerrnfgs_agn]
    yerrbig_agn = [yerr_agn,yerrnfgs_agn]

; average error bar
    
    avgerrlocal = djs_median([yerrbig,yerrbig_agn])

; ####################
; Maier et al. (2005)
; ####################
    
    cutm05 = where(m05.b_lum gt -900.0)
    lineratio, m05[cutm05], '', '', 'OII_3727', 'H_ALPHA', dum1, dum2, $
      ym05, yerrm05, index=indxm05, nindex=nindxm05, snrcut=snrcut_highz

    xm05 = m05[cutm05[indxm05]].b_lum
    xerrm05 = m05[cutm05[indxm05]].b_lum_err

    loiim05 = m05[cutm05[indxm05]].oii_3727_lum[0] + alog10(lsun)
    hasfrm05 = m05nodust[cutm05[indxm05]].sfr_h_alpha

    avgerrm05 = djs_median(yerrm05)
    splog, 'Maier: redshift = ', nindxm05, minmax(m05[cutm05[indxm05]].z_obj), $
      median(m05[cutm05[indxm05]].z_obj)

; ####################
; Tresse et al. (2002)    
; ####################
    
    cutt02 = where(t02.b_lum gt -900.0)
    lineratio, t02[cutt02], '', '', 'OII_3727', 'H_ALPHA', dum1, dum2, $
      yt02, yerrt02, index=indxt02, nindex=nindxt02, snrcut=snrcut_highz

    xt02 = t02[cutt02[indxt02]].b_lum
    xerrt02 = t02[cutt02[indxt02]].b_lum_err

    avgerrt02 = djs_median(yerrt02)
    splog, 'Tresse: redshift = ', nindxt02, minmax(t02[cutt02[indxt02]].z_obj), $
      median(t02[cutt02[indxt02]].z_obj)

; ####################
; Hicks et al. (2002)    
; ####################
    
    cuth02 = where(h02.b_lum gt -900.0)
    lineratio, h02[cuth02], '', '', 'OII_3727', 'H_ALPHA', dum1, dum2, $
      yh02, yerrh02, index=indxh02, nindex=nindxh02, snrcut=snrcut_highz

    xh02 = h02[cuth02[indxh02]].b_lum
    xerrh02 = h02[cuth02[indxh02]].b_lum_err

    avgerrh02 = djs_median(yerrh02)
    splog, 'Hicks: redshift = ', nindxh02, minmax(h02[cuth02[indxh02]].z_obj), $
      median(h02[cuth02[indxh02]].z_obj)

; ####################
; Glazebrook et al (1999)
; ####################
    
    cutg99 = where(g99.b_lum gt -900.0)
    lineratio, g99[cutg99], '', '', 'OII_3727', 'H_ALPHA', dum1, dum2, $
      yg99, yerrg99, index=indxg99, nindex=nindxg99, snrcut=snrcut_highz

    xg99 = g99[cutg99[indxg99]].b_lum
    xerrg99 = g99[cutg99[indxg99]].b_lum_err

    avgerrg99 = djs_median(yerrg99)
    splog, 'Glazebrook: redshift = ', nindxg99, minmax(g99[cutg99[indxg99]].z_obj), $
      median(g99[cutg99[indxg99]].z_obj)

; ####################
    
    xtitle = 'log L(B) [L(B)'+sunsymbol()+']'
    ytitle = 'log ([O II]/H\alpha)_{obs}'

    xrange = LBrange
    yrange = [-1.6,1.0] ; oiihacorrange

    djs_plot, [0], [0], /nodata, xthick=postthick, ythick=postthick, charsize=charsize_8, $
      charthick=postthick, thick=postthick, xtitle=xtitle, ytitle=ytitle, $
      ymargin=[4,3], xstyle=11, ystyle=3, xrange=xrange, yrange=yrange, position=pos[*,0]
    axis, /xaxis, xrange=interpol(xabs,x,!x.crange), xthick=postthick, xstyle=1, $
      charsize=charsize_8, charthick=postthick
    im_xyouts_title, xtitle='M_{B} [mag]', charsize=charsize_8, charthick=postthick

    im_symbols, localsym, psize=localpsize, /fill, color=fsc_color(localcolor,!d.table_size-101)
    djs_oplot, xbig, ybig, ps=8, color=fsc_color(localcolor,!d.table_size-101)

    im_symbols, localsym, psize=localpsize, fill=0, color=fsc_color(localcolor,!d.table_size-101), thick=symthick
    djs_oplot, xbig_agn, ybig_agn, ps=8, color=fsc_color(localcolor,!d.table_size-101)

; keep tresse on bottom    
    
    im_symbols, t02sym, psize=t02psize, /fill, color=fsc_color(t02color,!d.table_size-103), thick=postthick
    djs_oplot, xt02, yt02, ps=8, color=fsc_color(t02color,!d.table_size-103), thick=postthick

    im_symbols, m05sym, psize=m05psize, /fill, color=fsc_color(m05color,!d.table_size-102), thick=postthick
    djs_oplot, xm05, ym05, ps=8, color=fsc_color(m05color,!d.table_size-102), thick=postthick

    im_symbols, h02sym, psize=h02psize, /fill, color=fsc_color(h02color,!d.table_size-104), thick=postthick
    djs_oplot, xh02, yh02, ps=8, color=fsc_color(h02color,!d.table_size-104), thick=postthick

    im_symbols, g99sym, psize=g99psize, /fill, color=fsc_color(g99color,!d.table_size-105), thick=postthick
    djs_oplot, xg99, yg99, ps=8, color=fsc_color(g99color,!d.table_size-105), thick=postthick

; overplot the average error bars for the high-z points and the local 
; sample

    ypos = 0.7
    
    im_symbols, localsym, psize=localpsize, /fill, color=fsc_color(localcolor,!d.table_size-101)
    oploterror, 7.2, ypos, 0.0, avgerrlocal, /data, thick=postthick, $
      errthick=postthick, psym=8, errcolor=fsc_color(localcolor,!d.table_size-101)

    im_symbols, m05sym, psize=m05psize, /fill, color=fsc_color(m05color,!d.table_size-102), thick=postthick
    oploterror, 7.5, ypos, 0.0, avgerrm05, /data, thick=postthick, $
      errthick=postthick, psym=8, errcolor=fsc_color(m05color,!d.table_size-102)
    
    im_symbols, t02sym, psize=t02psize, /fill, color=fsc_color(t02color,!d.table_size-103), thick=postthick
    oploterror, 7.8, ypos, 0.0, avgerrt02, /data, thick=postthick, $
      errthick=postthick, psym=8, errcolor=fsc_color(t02color,!d.table_size-103)
    
    im_symbols, h02sym, psize=h02psize, /fill, color=fsc_color(h02color,!d.table_size-104), thick=postthick
    oploterror, 8.1, ypos, 0.0, avgerrh02, /data, thick=postthick, $
      errthick=postthick, psym=8, errcolor=fsc_color(h02color,!d.table_size-104)
    
    im_symbols, g99sym, psize=g99psize, /fill, color=fsc_color(g99color,!d.table_size-105), thick=postthick
    oploterror, 8.4, ypos, 0.0, avgerrg99, /data, thick=postthick, $
      errthick=postthick, psym=8, errcolor=fsc_color(g99color,!d.table_size-105)
    
; legend    
    
    label = ['Local Sample - Star-Forming','Local Sample - AGN',$
      'Maier et al. (2005)','Tresse et al. (2002)',$
      'Hicks et al. (2002)','Glazebrook et al. (1999)']
    psym = [localsym,localsym,m05sym,t02sym,h02sym,g99sym]
    fill = [1,0,1,1,1,1]
    color = fsc_color([localcolor,localcolor,m05color,t02color,h02color,g99color],!d.table_size-[101,101,102,103,104,105])
    postthick1 = postthick
    im_legend, textoidl(label), /left, /bottom, box=0, charsize=charsize_2, $
      charthick=postthick1, psym=psym, fill=fill, $
      spacing=1.8, thick=postthick1, textcolor=djs_icolor(replicate('',6))

    im_openclose, postscript=postscript, /close    

; ------------------------------------------------------------
; L(B) versus SFR(Ha)/L([O II])_obs
; ------------------------------------------------------------
    
    psname = 'lb_vs_sfr_ha_loii_obs'
    im_openclose, pspath+psname, postscript=postscript, xsize=8.5, ysize=8.8, encapsulated=encapsulated

    pagemaker, nx=1, ny=1, xspace=0, yspace=0.0, width=6.8, height=6.8, $
      xmargin=[1.4,0.3], ymargin=[0.9,1.1], xpage=8.5, ypage=8.8, $
      position=pos, /normal

    if keyword_set(blackwhite) then begin
       notlirgscolor = 'gray'
       lirgscolor = 'dark gray'
       nolirdatacolor = 'dark gray'
       medcolor = 'black'
    endif else begin
       notlirgscolor = 'gray'
       lirgscolor = 'medium orchid'
       nolirdatacolor = 'dodger blue'
       medcolor = 'black'
    endelse

    medbin = 0.5
    minpts = 3
    minx = 7.5 ; -3.0
    
; Atlas - SF
    
    indx = where((atlasdust.oii_3727[0]/atlasdust.oii_3727[1] gt snrcut) and $
      (atlasnodust.b_lum_obs gt -900.0) and (atlasnodust.ehbha_err gt 0.0),nindx)

    x = atlasdust[indx].b_lum_obs; - alog10(LBnorm)
    xerr = atlasdust[indx].b_lum_obs_err
    xabs = atlasdust[indx].m_b_obs

    oii = atlasdust[indx].oii_3727_lum[0] + alog10(lsun) - alog10(elumnorm)
    oii_err = atlasdust[indx].oii_3727_lum[1]

    sfr_ha = atlasnodust[indx].sfr_h_alpha
    sfr_ha_err = atlasnodust[indx].sfr_h_alpha_err

    y = sfr_ha - oii
    yerr = sqrt(sfr_ha_err^2 + oii_err^2)

; Atlas - AGN
    
    indx_agn = where((atlasdust_agn.oii_3727[0]/atlasdust_agn.oii_3727[1] gt snrcut) and $
      (atlasnodust_agn.b_lum_obs gt -900.0) and (atlasnodust_agn.ehbha_err gt 0.0),nindx_agn)

    x_agn = atlasdust_agn[indx_agn].b_lum_obs; - alog10(LBnorm)
    xerr_agn = atlasdust_agn[indx_agn].b_lum_obs_err

    oii_agn = atlasdust_agn[indx_agn].oii_3727_lum[0] + alog10(lsun) - alog10(elumnorm)
    oii_err_agn = atlasdust_agn[indx_agn].oii_3727_lum[1]

    sfr_ha_agn = atlasnodust_agn[indx_agn].sfr_h_alpha
    sfr_ha_err_agn = atlasnodust_agn[indx_agn].sfr_h_alpha_err

    y_agn = sfr_ha_agn - oii_agn
    yerr_agn = sqrt(sfr_ha_err_agn^2 + oii_err_agn^2)

; NFGS - SF
    
    indxnfgs = where((nfgsdust.oii_3727[0]/nfgsdust.oii_3727[1] gt snrcut) and $
      (nfgsnodust.b_lum_obs gt -900.0) and (nfgsnodust.ehbha_err gt 0.0),nindxnfgs)

    xnfgs = nfgsdust[indxnfgs].b_lum_obs; - alog10(LBnorm)
    xerrnfgs = nfgsdust[indxnfgs].b_lum_obs_err
    
    oiinfgs = nfgsdust[indxnfgs].oii_3727_lum[0] + alog10(lsun) - alog10(elumnorm)
    oiinfgs_err = nfgsdust[indxnfgs].oii_3727_lum[1]

    sfrnfgs_ha = nfgsnodust[indxnfgs].sfr_h_alpha
    sfrnfgs_ha_err = nfgsnodust[indxnfgs].sfr_h_alpha_err

    ynfgs = sfrnfgs_ha - oiinfgs
    yerrnfgs = sqrt(sfrnfgs_ha_err^2 + oiinfgs_err^2)

; NFGS - AGN
    
    indxnfgs_agn = where((nfgsdust_agn.oii_3727[0]/nfgsdust_agn.oii_3727[1] gt snrcut) and $
      (nfgsnodust_agn.b_lum_obs gt -900.0) and (nfgsnodust_agn.ehbha_err gt 0.0),nindxnfgs_agn)

    xnfgs_agn = nfgsdust_agn[indxnfgs_agn].b_lum_obs; - alog10(LBnorm)
    xerrnfgs_agn = nfgsdust_agn[indxnfgs_agn].b_lum_obs_err
    
    oiinfgs_agn = nfgsdust_agn[indxnfgs_agn].oii_3727_lum[0] + alog10(lsun) - alog10(elumnorm)
    oiinfgs_err_agn = nfgsdust_agn[indxnfgs_agn].oii_3727_lum[1]

    sfrnfgs_ha_agn = nfgsnodust_agn[indxnfgs_agn].sfr_h_alpha
    sfrnfgs_ha_err_agn = nfgsnodust_agn[indxnfgs_agn].sfr_h_alpha_err

    ynfgs_agn = sfrnfgs_ha_agn - oiinfgs_agn
    yerrnfgs_agn = sqrt(sfrnfgs_ha_err_agn^2 + oiinfgs_err_agn^2)

; combine the two samples - SF

    xbig = [x,xnfgs] & xerrbig = [xerr,xerrnfgs]
    ybig = [y,ynfgs] & yerrbig = [yerr,yerrnfgs]
    lir = [atlasdust[indx].ir_lum,nfgsdust[indxnfgs].ir_lum]

    lirgs = where((lir gt -900.0) and (lir gt 11.0))
    notlirgs = where((lir gt -900.0) and (lir lt 11.0))
    nolirdata = where((lir lt -900.0))
    
; combine the two samples - AGN

    xbig_agn = [x_agn,xnfgs_agn] & xerrbig_agn = [xerr_agn,xerrnfgs_agn]
    ybig_agn = [y_agn,ynfgs_agn] & yerrbig_agn = [yerr_agn,yerrnfgs_agn]
    lir_agn = [atlasdust_agn[indx_agn].ir_lum,nfgsdust_agn[indxnfgs_agn].ir_lum]

    lirgs_agn = where((lir_agn gt -900.0) and (lir_agn gt 11.0))
    notlirgs_agn = where((lir_agn gt -900.0) and (lir_agn lt 11.0))
    nolirdata_agn = where((lir_agn lt -900.0))
    
; make the plot

    xtitle = 'log L(B) [L(B)'+sunsymbol()+']'
    ytitle = 'log [10^{41} \psi/L([O II])_{obs}] [erg s^{-1}/'+sfr_units()+']'

    xrange = LBrange; - alog10(LBnorm)
    yrange = alog10([3.0D-42,15D-40]) - alog10(sfrnorm) + alog10(elumnorm)
    
; SF galaxies    
    
    atlas1d_lineplot, xbig[notlirgs], ybig[notlirgs], xerrbig[notlirgs], yerrbig[notlirgs], $
      plottype=1, postscript=postscript, xtitle=xtitle, ytitle=ytitle, xrange=xrange, $
      xstyle=11, yrange=yrange, legendtype=0, /right, /top, yfrac=1.5, xfrac=1.5, charsize=charsize_8, $
      atlaspsize=1.1, position=pos[*,0], atlascolor=notlirgscolor, blackwhite=blackwhite
    axis, /xaxis, xrange=interpol(xabs,x,!x.crange), xthick=postthick, xstyle=1, $
      charsize=charsize_8, charthick=postthick
    im_xyouts_title, xtitle='M_{B} [mag]', charsize=charsize_8, charthick=postthick

    im_symbols, 108, psize=1.4, /fill, color=fsc_color(lirgscolor,!d.table_size-5)
    djs_oplot, xbig[lirgs], ybig[lirgs], ps=8

    im_symbols, 105, psize=1.4, /fill, color=fsc_color(nolirdatacolor,!d.table_size-6)
    djs_oplot, xbig[nolirdata], ybig[nolirdata], ps=8

; AGN galaxies    
    
    atlas1d_lineplot, xbig_agn[notlirgs_agn], ybig_agn[notlirgs_agn], $
      xerrbig_agn[notlirgs_agn], yerrbig_agn[notlirgs_agn], blackwhite=blackwhite, $
      plottype=1, postscript=postscript, xtitle=xtitle, ytitle=ytitle, xrange=xrange, $
      yrange=yrange, legendtype=0, /right, /top, yfrac=1.5, xfrac=1.5, charsize=charsize_8, $
      atlaspsize=1.0, position=pos[*,0], atlasfill=0, /overplot, thick=postthick, atlascolor=notlirgscolor

    im_symbols, 108, psize=1.3, fill=0, color=fsc_color(lirgscolor,!d.table_size-5), thick=postthick
    djs_oplot, xbig_agn[lirgs_agn], ybig_agn[lirgs_agn], ps=8

    im_symbols, 105, psize=1.2, fill=0, color=fsc_color(nolirdatacolor,!d.table_size-6), thick=postthick
    djs_oplot, xbig_agn[nolirdata_agn], ybig_agn[nolirdata_agn], ps=8

    label = ['L(IR) < 10^{11} L'+sunsymbol(),'L(IR) > 10^{11} L'+sunsymbol(),$
      'Undetected with IRAS']
    
;   im_legend, textoidl(label), /left, /top, box=0, charsize=charsize_8, $
;     charthick=postthick, psym=[106,108,105], /fill, symsize=1.3, $
;     spacing=1.7, color=djs_icolor(['blue','red','dark green'])

; fit only to certain data

;   xfit = [xbig,xbig_agn]
;   yfit = [ybig,ybig_agn]
    xfit = xbig
    yfit = ybig
;   xfit = [xbig[notlirgs],xbig[nolirdata]]
;   yfit = [ybig[notlirgs],ybig[nolirdata]]
    
; 9.0 is M_B = -17.07 and 8.7 is M_B = -16.3; 8.57 is -16

    LBlocut = 10.0^6.0  ; 10^7.5
    LBhicut = 10.0^11.0 ; 10^10.8

    LBlo = alog10(LBlocut)>min(xfit)
    LBhi = alog10(LBhicut)<max(xfit)
;   LBlo = alog10(LBlocut/LBnorm)>min(xfit)
;   LBhi = alog10(LBhicut/LBnorm)<max(xfit)

;   djs_oplot, LBlo*[1,1], !y.crange, line=2
;   djs_oplot, LBhi*[1,1], !y.crange, line=2
    
; overlay the running median    

    running = im_medxbin(xfit,yfit,medbin,minx=minx,minpts=minpts,/verbose)

    doit = where((running.binctr gt LBlo) and (running.binctr le LBhi),ndoit)
    lolum = where((running.binctr lt LBlo))

    im_symbols, 108, psize=3.0, thick=postthick, color=medcolor;, /fill
;   oploterror, running.binctr[doit], running.medy[doit], $
;     running.sigy[doit], ps=8, errthick=(postthick-2L)>2L
    oploterror, running.binctr[doit], running.medy[doit], running.binsz/2.0+fltarr(ndoit), $
      running.sigy75[doit]-running.medy[doit], ps=8, errthick=(postthick)>2L, /hibar, $
      color=djs_icolor(medcolor), errcolor=djs_icolor(medcolor)
    oploterror, running.binctr[doit], running.medy[doit], running.binsz/2.0+fltarr(ndoit), $
      running.medy[doit]-running.sigy25[doit], ps=8, errthick=(postthick)>2L, /lobar, $
      color=djs_icolor(medcolor), errcolor=djs_icolor(medcolor)
;   oploterror, running.binctr[doit], running.medy[doit], running.binsz/2.0+fltarr(ndoit), $
;     running.sigyup[doit]-running.medy[doit], ps=8, errthick=(postthick-1L)>2L, /hibar, $
;     color=djs_icolor(medcolor), errcolor=djs_icolor(medcolor)
;   oploterror, running.binctr[doit], running.medy[doit], running.binsz/2.0+fltarr(ndoit), $
;     running.medy[doit]-running.sigylo[doit], ps=8, errthick=(postthick-1L)>2L, /lobar, $
;     color=djs_icolor(medcolor), errcolor=djs_icolor(medcolor)
;   oploterror, running.binctr[doit], running.medy[doit], $
;     running.sigyup[doit]-running.medy[doit], ps=8, errthick=(postthick-1L)>2L, /hibar
;   oploterror, running.binctr[doit], running.medy[doit], $
;     running.medy[doit]-running.sigylo[doit], ps=8, errthick=(postthick-1L)>2L, /lobar

    if keyword_set(paper) then begin

; generate a binary FITS table with the statistics

       oii_sfrfile = 'oii_sfr.fits'

       oii_sfr = replicate(template_sfr,ndoit)
       oii_sfr.loglb  = running.binctr[doit]
       oii_sfr.mb     = interpol(xabs,x,running.binctr[doit])
       oii_sfr.p25    = running.sigy25[doit]
       oii_sfr.p50    = running.medy[doit]
       oii_sfr.p75    = running.sigy75[doit]
       oii_sfr.mean   = running.meany[doit]
       oii_sfr.stddev = running.stddev[doit]
       
       mwrfits, oii_sfr, sfrspath+oii_sfrfile, /create

; generate a latex table with the statistics

       openw, lun, latexpath+'oii_sfr.tex', /get_lun
       printf, lun, '\begin{deluxetable}{ccccccc}'
       printf, lun, '\tabletypesize{\small}'
       printf, lun, '\tablecolumns{7}'
       printf, lun, '\tablecaption{\oiilam{} Star-Formation Rate Conversion Factors \label{table:oii_sfr}}'
       printf, lun, '\tablewidth{0in}'
       printf, lun, '\tablehead{'
       printf, lun, '\colhead{$\log\,\lb$} & '
       printf, lun, '\colhead{\mb} & '
       printf, lun, '\multicolumn{5}{c}{$\log\,[\sfr/\loiiobs]$} \\'

       printf, lun, '\colhead{[\lbsun]} & '
       printf, lun, '\colhead{[mag]} & '
       printf, lun, '\multicolumn{5}{c}{$[10^{41}\,\lunits/(\sfrunits)]$} \\'

       printf, lun, '\cline{3-7}'

       printf, lun, '\multicolumn{2}{c}{} & '
       printf, lun, '\colhead{$P_{25}$} & '
       printf, lun, '\colhead{$P_{50}$} & '
       printf, lun, '\colhead{$P_{75}$} & '
       printf, lun, '\colhead{$\langle R\rangle$} & '
       printf, lun, '\colhead{$\sigma_{R}$}'
       printf, lun, '}'
       printf, lun, '\startdata'

       for i = 0L, ndoit-1L do printf, lun, $
         running.binctr[doit[i]], ' & ', interpol(xabs,x,running.binctr[doit[i]]), ' & ', $
         running.sigy25[doit[i]], ' & ', running.medy[doit[i]], ' & ', $
         running.sigy75[doit[i]], ' & ', running.meany[doit[i]], ' & ', $
         running.stddev[doit[i]], ' \\', $
         format='(F5.2,A3,F6.2,A3,F6.3,A3,F6.3,A3,F6.3,A3,F6.3,A3,F6.3,A3)'

       printf, lun, '\enddata'
       printf, lun, '%\tablenotetext{a}{}'
       printf, lun, '\tablecomments{The columns labeled $P_{25}$, $P_{50}$, and $P_{75}$ '+$
         'give the $25$, $50$ (median), and $75$ percentile of the $\log\,[\sfr/\loiiobs]$ '+$
         'distribution, respectively, in bins of $0.5$~dex in luminosity.  $\langle R\rangle$ '+$
         'and $\sigma_{R}$ give the mean and standard-deviation of the distribution in each bin.}'
       printf, lun, '%\tablerefs{}'
       printf, lun, '\end{deluxetable}'
       free_lun, lun

    endif

; assume a constant ratio at low luminosity    
;
;   xbestfit = findgen((LBlo-min(xfit))/0.01+1L)*0.01+min(xfit)
;   constant = mean(running.medy[lolum])
;   constant_err = mean(running.sigy[lolum])
;   ybestfit = xbestfit*0.0+constant
;   
;   djs_oplot, xbestfit, ybestfit, thick=postthick2, line=0 ;, color='orange'
;
;   splog, 'L(B)/L(B)_sun vs SFR(Ha)/L([O II])_obs coefficients: '
;   niceprint, 10^constant/1D40, constant_err*10^constant*alog(10.0)/1D40
    
; fit the high-luminosity end
    
;   xfitaxis = findgen((LBhi-LBlo)/0.05+1L)*0.05 + LBlo
;   yfitaxis = interpol(running.medy[doit],running.binctr[doit],xfitaxis)
;   yerraxis = sqrt(interpol(running.sigy[doit]^2,running.binctr[doit],xfitaxis))

;   xfitaxis = running.binctr[doit]
;   yfitaxis = running.medy[doit]
;   yerraxis = running.sigy[doit]
    
;   coeff = im_linefit(xfitaxis-LBlo,yfitaxis,yerr=yerraxis,coeff_fixed=[0,0],$
;     coeff_guess=[1.0,0.25],coeff_err=coeff_err)
;   coeff = linfit(xfitaxis-LBlo,yfitaxis,measure_errors=yerraxis,sigma=coeff_err)
;   coeff = linfit(xfitaxis-LBlo,yfitaxis,measure_errors=yerraxis,sigma=coeff_err)
;   coeff = im_linefit(xfitaxis-LBlo,yfitaxis,coeff_fixed=[1,0],$
;     coeff_guess=[constant,0.25],coeff_err=coeff_err)

;   coeff = linfit(running.binctr[doit],running.medy[doit])
;   sixlin, running.binctr[doit], running.medy[doit], a, siga, b, sigb
;   sixlin, xfitaxis, yfitaxis, a, siga, b, sigb
;   coeff = [a[2],b[2]]

;   splog, 'L(B)/L(B)_sun vs SFR(Ha)/L([O II])_obs coefficients: '
;   niceprint, LBlocut^(-coeff[1])*10^(coeff[0]-40.0), coeff[1]
;
;   xbestfit = findgen((LBhi-LBlo)/0.01+1L)*0.01+LBlo
;   ybestfit = poly(xbestfit-LBlo,coeff)
;   djs_oplot, xbestfit, ybestfit, thick=postthick2, line=0;, color='orange'

; compute the residual scatter

;   hilum = where(xfit gt lblo)
;   resid = yfit[hilum] - poly(xfit[hilum]-lblo,coeff)
;   jj = im_stats(resid,/verbose)

    im_openclose, postscript=postscript, /close    

; ------------------------------------------------------------
; L(B) versus SFR(Ha)/L(Hb)_obs - Integrated
; ------------------------------------------------------------

    psname = 'lb_vs_sfr_ha_lhb_obs'
    im_openclose, pspath+psname, postscript=postscript, xsize=8.5, ysize=8.8, encapsulated=encapsulated

    pagemaker, nx=1, ny=1, xspace=0, yspace=0.0, width=6.8, height=6.8, $
      xmargin=[1.4,0.3], ymargin=[0.9,1.1], xpage=8.5, ypage=8.8, $
      position=pos, /normal

    if keyword_set(blackwhite) then begin
       notlirgscolor = 'gray'
       lirgscolor = 'dark gray'
       nolirdatacolor = 'dark gray'
       medcolor = 'black'
    endif else begin
       notlirgscolor = 'gray'
       lirgscolor = 'medium orchid'
       nolirdatacolor = 'dodger blue'
       medcolor = 'black'
    endelse

    medbin = 0.5
    minpts = 3
    minx = 7.0 ; -3.0
    
; Atlas - SF
    
    indx = where((atlasdust.b_lum_obs gt -900) and (atlasnodust.sfr_h_alpha gt -900.0) and $
      (atlasdust.h_beta_ew_uncor[0] ge ewcut),nindx)

    x = atlasdust[indx].b_lum_obs; - alog10(LBnorm)
    xerr = atlasdust[indx].b_lum_obs_err
    xabs = atlasdust[indx].m_b_obs

    hb = atlasdust[indx].h_beta_lum[0] + alog10(lsun) - alog10(elumnorm)
    hb_err = atlasdust[indx].h_beta_lum[1]

    sfr_ha = atlasnodust[indx].sfr_h_alpha
    sfr_ha_err = atlasnodust[indx].sfr_h_alpha_err

    y = sfr_ha - hb
    yerr = sqrt(sfr_ha_err^2 + hb_err^2)
    
; Atlas - AGN
    
    indx_agn = where((atlasdust_agn.b_lum_obs gt -900) and (atlasnodust_agn.sfr_h_alpha gt -900.0) and $
      (atlasdust_agn.h_beta_ew_uncor[0] ge ewcut),nindx_agn)

    x_agn = atlasdust_agn[indx_agn].b_lum_obs; - alog10(LBnorm)
    xerr_agn = atlasdust_agn[indx_agn].b_lum_obs_err

    hb_agn = atlasdust_agn[indx_agn].h_beta_lum[0] + alog10(lsun) - alog10(elumnorm)
    hb_err_agn = atlasdust_agn[indx_agn].h_beta_lum[1]

    sfr_ha_agn = atlasnodust_agn[indx_agn].sfr_h_alpha
    sfr_ha_err_agn = atlasnodust_agn[indx_agn].sfr_h_alpha_err

    y_agn = sfr_ha_agn - hb_agn
    yerr_agn = sqrt(sfr_ha_err_agn^2 + hb_err_agn^2)
    
; NFGS - SF
    
    indxnfgs = where((nfgsdust.b_lum_obs gt -900.0) and (nfgsnodust.sfr_h_alpha gt -900.0) and $
      (nfgsdust.h_beta_ew_uncor[0] ge ewcut),nindxnfgs)

    xnfgs = nfgsdust[indxnfgs].b_lum_obs; - alog10(LBnorm)
    xerrnfgs = nfgsdust[indxnfgs].b_lum_obs_err
    
    hbnfgs = nfgsdust[indxnfgs].h_beta_lum[0] + alog10(lsun) - alog10(elumnorm)
    hbnfgs_err = nfgsdust[indxnfgs].h_beta_lum[1]

    sfrnfgs_ha = nfgsnodust[indxnfgs].sfr_h_alpha
    sfrnfgs_ha_err = nfgsnodust[indxnfgs].sfr_h_alpha_err

    ynfgs = sfrnfgs_ha - hbnfgs
    yerrnfgs = sqrt(sfrnfgs_ha_err^2 + hbnfgs_err^2)

; NFGS - AGN
    
    indxnfgs_agn = where((nfgsdust_agn.b_lum_obs gt -900.0) and (nfgsnodust_agn.sfr_h_alpha gt -900.0) and $
      (nfgsdust_agn.h_beta_ew_uncor[0] ge ewcut),nindxnfgs_agn)

    xnfgs_agn = nfgsdust_agn[indxnfgs_agn].b_lum_obs; - alog10(LBnorm)
    xerrnfgs_agn = nfgsdust_agn[indxnfgs_agn].b_lum_obs_err
    
    hbnfgs_agn = nfgsdust_agn[indxnfgs_agn].h_beta_lum[0] + alog10(lsun) - alog10(elumnorm)
    hbnfgs_err_agn = nfgsdust_agn[indxnfgs_agn].h_beta_lum[1]

    sfrnfgs_ha_agn = nfgsnodust_agn[indxnfgs_agn].sfr_h_alpha
    sfrnfgs_ha_err_agn = nfgsnodust_agn[indxnfgs_agn].sfr_h_alpha_err

    ynfgs_agn = sfrnfgs_ha_agn - hbnfgs_agn
    yerrnfgs_agn = sqrt(sfrnfgs_ha_err_agn^2 + hbnfgs_err_agn^2)

; combine the two samples - SF

    xbig = [x,xnfgs] & xerrbig = [xerr,xerrnfgs]
    ybig = [y,ynfgs] & yerrbig = [yerr,yerrnfgs]
    lir = [atlasdust[indx].ir_lum,nfgsdust[indxnfgs].ir_lum]
    gal = [atlasdust[indx].galaxy,nfgsdust[indxnfgs].galaxy]

    lirgs = where((lir gt -900.0) and (lir gt 11.0))
    notlirgs = where((lir gt -900.0) and (lir lt 11.0))
    nolirdata = where((lir lt -900.0))

;   w = where(xbig[nolirdata] gt 9.5)
;   niceprint, gal[nolirdata[w]], xbig[nolirdata[w]], ybig[nolirdata[w]], lir[nolirdata[w]]
    
; combine the two samples - AGN

    xbig_agn = [x_agn,xnfgs_agn] & xerrbig_agn = [xerr_agn,xerrnfgs_agn]
    ybig_agn = [y_agn,ynfgs_agn] & yerrbig_agn = [yerr_agn,yerrnfgs_agn]
    lir_agn = [atlasdust_agn[indx_agn].ir_lum,nfgsdust_agn[indxnfgs_agn].ir_lum]

    lirgs_agn = where((lir_agn gt -900.0) and (lir_agn gt 11.0))
    notlirgs_agn = where((lir_agn gt -900.0) and (lir_agn lt 11.0))
    nolirdata_agn = where((lir_agn lt -900.0))
    
; make the plot

    xtitle = 'log L(B) [L(B)'+sunsymbol()+']'
    ytitle = 'log [10^{41} \psi/L(H\beta)_{obs}] [erg s^{-1}/'+sfr_units()+']'

    xrange = LBrange; - alog10(LBnorm)
    yrange = alog10([1.8D-41,3D-40]) - alog10(sfrnorm) + alog10(elumnorm)

; SF galaxies    
    
    atlas1d_lineplot, xbig[notlirgs], ybig[notlirgs], xerrbig[notlirgs], yerrbig[notlirgs], $
      plottype=1, postscript=postscript, xtitle=xtitle, ytitle=ytitle, xrange=xrange, $
      xstyle=11, yrange=yrange, legendtype=0, /right, /top, yfrac=1.5, xfrac=1.5, charsize=charsize_8, $
      atlaspsize=1.1, position=pos[*,0], atlascolor=notlirgscolor, blackwhite=blackwhite
    axis, /xaxis, xrange=interpol(xabs,x,!x.crange), xthick=postthick, xstyle=1, $
      charsize=charsize_8, charthick=postthick
    im_xyouts_title, xtitle='M_{B} [mag]', charsize=charsize_8, charthick=postthick

    im_symbols, 108, psize=1.4, /fill, color=fsc_color(lirgscolor,!d.table_size-5)
    djs_oplot, xbig[lirgs], ybig[lirgs], ps=8

    im_symbols, 105, psize=1.4, /fill, color=fsc_color(nolirdatacolor,!d.table_size-6)
    djs_oplot, xbig[nolirdata], ybig[nolirdata], ps=8

; AGN galaxies    
    
    atlas1d_lineplot, xbig_agn[notlirgs_agn], ybig_agn[notlirgs_agn], $
      xerrbig_agn[notlirgs_agn], yerrbig_agn[notlirgs_agn], blackwhite=blackwhite, $
      plottype=1, postscript=postscript, xtitle=xtitle, ytitle=ytitle, xrange=xrange, $
      yrange=yrange, legendtype=0, /right, /top, yfrac=1.5, xfrac=1.5, charsize=charsize_8, $
      atlaspsize=1.0, position=pos[*,0], atlasfill=0, /overplot, thick=postthick, atlascolor=notlirgscolor

    im_symbols, 108, psize=1.3, fill=0, color=fsc_color(lirgscolor,!d.table_size-5), thick=postthick
    djs_oplot, xbig_agn[lirgs_agn], ybig_agn[lirgs_agn], ps=8

    im_symbols, 105, psize=1.2, fill=0, color=fsc_color(nolirdatacolor,!d.table_size-6), thick=postthick
    djs_oplot, xbig_agn[nolirdata_agn], ybig_agn[nolirdata_agn], ps=8

    label = ['L(IR) < 10^{11} L'+sunsymbol(),'L(IR) > 10^{11} L'+sunsymbol(),$
      'Undetected with IRAS']
    
;   im_legend, textoidl(label), /left, /top, box=0, charsize=charsize_8, $
;     charthick=postthick, psym=[106,108,105], /fill, symsize=1.3, $
;     spacing=1.7, color=djs_icolor(['blue','red','dark green'])

; fit only to certain data

;   xfit = [xbig,xbig_agn]
;   yfit = [ybig,ybig_agn]
    xfit = xbig
    yfit = ybig
;   xfit = [xbig[notlirgs],xbig[nolirdata]]
;   yfit = [ybig[notlirgs],ybig[nolirdata]]
    
; 9.0 is M_B = -17.07 and 8.7 is M_B = -16.3; 8.57 is -16

    LBlocut = 10.0^6.0 ; 10.0^8.7
    LBhicut = 10.0^11.0

    LBlo = alog10(LBlocut)>min(xfit)
    LBhi = alog10(LBhicut)<max(xfit)
;   LBlo = alog10(LBlocut/LBnorm)>min(xfit)
;   LBhi = alog10(LBhicut/LBnorm)<max(xfit)

; overlay the running median    

    running = im_medxbin(xfit,yfit,medbin,minx=minx,minpts=minpts,/verbose)

    doit = where((running.binctr gt LBlo) and (running.binctr le LBhi),ndoit)
    lolum = where((running.binctr lt LBlo))

    im_symbols, 108, psize=3.0, thick=postthick, color=medcolor;, /fill
;   oploterror, running.binctr[doit], running.medy[doit], $
;     running.sigy[doit], ps=8, errthick=(postthick-2L)>2L
    oploterror, running.binctr[doit], running.medy[doit], running.binsz/2.0+fltarr(ndoit), $
      running.sigy75[doit]-running.medy[doit], ps=8, errthick=(postthick)>2L, /hibar, $
      color=djs_icolor(medcolor), errcolor=djs_icolor(medcolor)
    oploterror, running.binctr[doit], running.medy[doit], running.binsz/2.0+fltarr(ndoit), $
      running.medy[doit]-running.sigy25[doit], ps=8, errthick=(postthick)>2L, /lobar, $
      color=djs_icolor(medcolor), errcolor=djs_icolor(medcolor)
;   oploterror, running.binctr[doit], running.medy[doit], running.binsz/2.0+fltarr(ndoit), $
;     running.sigyup[doit]-running.medy[doit], ps=8, errthick=(postthick-1L)>2L, /hibar, $
;     color=djs_icolor(medcolor), errcolor=djs_icolor(medcolor)
;   oploterror, running.binctr[doit], running.medy[doit], running.binsz/2.0+fltarr(ndoit), $
;     running.medy[doit]-running.sigylo[doit], ps=8, errthick=(postthick-1L)>2L, /lobar, $
;     color=djs_icolor(medcolor), errcolor=djs_icolor(medcolor)
;   oploterror, running.binctr[doit], running.medy[doit], $
;     running.sigyup[doit]-running.medy[doit], ps=8, errthick=(postthick-1L)>2L, /hibar
;   oploterror, running.binctr[doit], running.medy[doit], $
;     running.medy[doit]-running.sigylo[doit], ps=8, errthick=(postthick-1L)>2L, /lobar

    if keyword_set(paper) then begin
    
; generate a binary FITS table with the statistics

       hb_sfrfile = 'hb_sfr.fits'

       hb_sfr = replicate(template_sfr,ndoit)
       hb_sfr.loglb  = running.binctr[doit]
       hb_sfr.mb     = interpol(xabs,x,running.binctr[doit])
       hb_sfr.p25    = running.sigy25[doit]
       hb_sfr.p50    = running.medy[doit]
       hb_sfr.p75    = running.sigy75[doit]
       hb_sfr.mean   = running.meany[doit]
       hb_sfr.stddev = running.stddev[doit]
       
       mwrfits, hb_sfr, sfrspath+hb_sfrfile, /create

; generate a latex table with the statistics

       openw, lun, latexpath+'hb_sfr.tex', /get_lun
       printf, lun, '\begin{deluxetable}{ccccccc}'
       printf, lun, '\tabletypesize{\small}'
       printf, lun, '\tablecolumns{7}'
       printf, lun, '\tablecaption{\hblam{} Star-Formation Rate Conversion Factors \label{table:hb_sfr}}'
       printf, lun, '\tablewidth{0in}'
       printf, lun, '\tablehead{'
       printf, lun, '\colhead{$\log\,\lb$} & '
       printf, lun, '\colhead{\mb} & '
       printf, lun, '\multicolumn{5}{c}{$\log\,[\sfr/\lhbobs]$} \\'

       printf, lun, '\colhead{[\lbsun]} & '
       printf, lun, '\colhead{[mag]} & '
       printf, lun, '\multicolumn{5}{c}{$[10^{41}\,\lunits/(\sfrunits)]$} \\'

       printf, lun, '\cline{3-7}'

       printf, lun, '\multicolumn{2}{c}{} & '
       printf, lun, '\colhead{$P_{25}$} & '
       printf, lun, '\colhead{$P_{50}$} & '
       printf, lun, '\colhead{$P_{75}$} & '
       printf, lun, '\colhead{$\langle R\rangle$} & '
       printf, lun, '\colhead{$\sigma_{R}$}'
       printf, lun, '}'
       printf, lun, '\startdata'

       for i = 0L, ndoit-1L do printf, lun, $
         running.binctr[doit[i]], ' & ', interpol(xabs,x,running.binctr[doit[i]]), ' & ', $
         running.sigy25[doit[i]], ' & ', running.medy[doit[i]], ' & ', $
         running.sigy75[doit[i]], ' & ', running.meany[doit[i]], ' & ', $
         running.stddev[doit[i]], ' \\', $
         format='(F5.2,A3,F6.2,A3,F5.3,A3,F5.3,A3,F5.3,A3,F5.3,A3,F5.3,A3)'

       printf, lun, '\enddata'
       printf, lun, '%\tablenotetext{a}{}'
       printf, lun, '\tablecomments{The columns labeled $P_{25}$, $P_{50}$, and $P_{75}$ '+$
         'give the $25$, $50$ (median), and $75$ percentile of the $\sfr/\lhbobs$ '+$
         'distribution, respectively, in bins of $0.5$~dex in luminosity.  $\langle R\rangle$ '+$
         'and $\sigma_{R}$ give the mean and standard-deviation of the distribution in each bin.}'
       printf, lun, '%\tablerefs{}'
       printf, lun, '\end{deluxetable}'
       free_lun, lun

    endif

; assume a constant ratio at low luminosity    

;   xbestfit = findgen((LBlo-min(xfit))/0.01+1L)*0.01+min(xfit)
;   constant = mean(running.medy[lolum])
;   constant_err = mean(running.sigy[lolum])
;   ybestfit = xbestfit*0.0+constant
    
;   djs_oplot, xbestfit, ybestfit, thick=postthick2, line=0 ;, color='orange'

;   splog, 'log L(B)/L(B)_sun vs log [SFR(Ha)/L(Hb)_obs] low-luminosity constant: ', constant, constant_err
;   splog, 'L(B)/L(B)_sun vs SFR(Ha)/L(Hb)_obs coefficients: '
;   niceprint, 10^constant/1D40, constant_err*10^constant*alog(10.0)/1D40
    
; now fit the high-luminosity end; constrain the intercept
    
    xfitaxis = findgen((LBhi-LBlo)/0.05+1L)*0.05 + LBlo
    yfitaxis = interpol(running.medy[doit],running.binctr[doit],xfitaxis)
    yerraxis = sqrt(interpol(running.sigy[doit]^2,running.binctr[doit],xfitaxis))

;   coeff = im_linefit(xfitaxis-LBlo,yfitaxis,yerr=yerraxis,coeff_fixed=[1,0],$
;     coeff_guess=[constant,0.25],coeff_err=coeff_err)
;   coeff = im_linefit(xfitaxis-LBlo,yfitaxis,coeff_fixed=[1,0],$
;     coeff_guess=[constant,0.25],coeff_err=coeff_err)

;   coeff = linfit(running.binctr[doit],running.medy[doit])
;   sixlin, running.binctr[doit], running.medy[doit], a, siga, b, sigb
;   sixlin, xfitaxis, yfitaxis, a, siga, b, sigb
;   coeff = [a[2],b[2]]

;   splog, 'log L(B)/L(B)_sun vs log [SFR(Ha)/L(Hb)_obs] coefficients: '
;   niceprint, coeff, coeff_err
;   splog, 'L(B)/L(B)_sun vs SFR(Ha)/L(Hb)_obs coefficients: '
;   niceprint, LBlocut^(-coeff[1])*10^(coeff[0]-40.0), coeff[1]

;   xbestfit = findgen((LBhi-LBlo)/0.01+1L)*0.01+LBlo
;   ybestfit = poly(xbestfit-LBlo,coeff)
;   djs_oplot, xbestfit, ybestfit, thick=postthick2, line=0;, color='orange'

; compute the residual scatter

;   hilum = where(xfit gt lblo)
;   resid = yfit[hilum] - poly(xfit[hilum]-lblo,coeff)
;   jj = im_stats(resid,/verbose)

;   lum = lsun*10^atlasdust[indx].h_beta_lum[0]
;   blum = 10^atlasdust[indx].b_lum_obs
;   hilum = where(blum gt LBlocut)
;   hbsfr = 2.45D-43*lum[hilum]*blum[hilum]^0.23
    
    im_openclose, postscript=postscript, /close    

; ------------------------------------------------------------
; L(B) vs 12+log(O/H) [Integrated]
; ------------------------------------------------------------

    psname = 'lb_vs_12oh_oiiinii_niiha'
    im_openclose, pspath+psname, postscript=postscript, xsize=8.5, ysize=8.3, encapsulated=encapsulated

; the size of this figure is meant to match lb_vs_ehbha.eps

    pagemaker, nx=1, ny=1, xspace=0.0, yspace=0.0, width=6.3, height=6.3, $
      xmargin=[1.1,1.1], ymargin=[0.9,1.1], xpage=8.5, ypage=8.3, $
      position=pos, /normal

    indx = where((atlasdust.b_lum_obs gt -900.0) and (atlasdust.zstrong_12oh_oiiinii_niiha gt -900.0),nindx)

    x = atlasdust[indx].b_lum_obs
    xerr = atlasdust[indx].b_lum_obs_err
    xabs = atlasdust[indx].m_b_obs
    
    y = atlasdust[indx].zstrong_12oh_oiiinii_niiha
    yerr = atlasdust[indx].zstrong_12oh_oiiinii_niiha_err

    indxnfgs = where((nfgsdust.b_lum_obs gt -900) and (nfgsdust.zstrong_12oh_oiiinii_niiha gt -900.0),nindxnfgs)

    xnfgs = nfgsdust[indxnfgs].b_lum_obs
    xerrnfgs = nfgsdust[indxnfgs].b_lum_obs_err
    
    ynfgs = nfgsdust[indxnfgs].zstrong_12oh_oiiinii_niiha
    yerrnfgs = nfgsdust[indxnfgs].zstrong_12oh_oiiinii_niiha_err

    xtitle = 'log L(B) [L(B)'+sunsymbol()+']'
    ytitle = '12 + log (O/H)'

    xrange = LBrange
    yrange = ohrange2
    
    atlas1d_lineplot, x, y, xerr, yerr, plottype=3, postscript=postscript, ysty=3, $
      ytitle=ytitle, xtitle=xtitle, xrange=xrange, yrange=yrange, legendtype=0, $
      xstyle=11, /right, /top, xnfgs=xnfgs, ynfgs=ynfgs, xerrnfgs=xerrnfgs, blackwhite=blackwhite, $
      yerrnfgs=yerrnfgs, position=pos[*,0], charsize=charsize_9, atlaspsize=1.1, nfgspsize=1.4
    axis, /xaxis, xrange=interpol(xabs,x,!x.crange), xthick=postthick, xstyle=1, $
      charsize=charsize_9, charthick=postthick
    im_xyouts_title, xtitle='M_{B} [mag]', charsize=charsize_9, charthick=postthick

    im_openclose, postscript=postscript, /close

; ------------------------------------------------------------
; L(B) vs E(Hb-Ha) [Integrated]
; ------------------------------------------------------------

    psname = 'lb_vs_ehbha'
    im_openclose, pspath+psname, postscript=postscript, xsize=8.5, ysize=8.3, encapsulated=encapsulated

    pagemaker, nx=1, ny=1, xspace=0, yspace=0.0, width=6.3, height=6.3, $
      xmargin=[1.1,1.1], ymargin=[0.9,1.1], xpage=8.5, ypage=8.3, $
      position=pos, /normal

    indx = where((atlasdust.b_lum_obs gt -900.0) and (atlasnodust.ehbha_err gt 0.0),nindx)

    x = atlasdust[indx].b_lum_obs
    xerr = atlasdust[indx].b_lum_obs_err
    xabs = atlasdust[indx].m_b_obs
    
    y = atlasnodust[indx].ehbha
    yerr = atlasnodust[indx].ehbha_err
    yebv = y/(k_lambda(4861,/odonnell)-k_lambda(6563,/odonnell))

    indxnfgs = where((nfgsdust.b_lum_obs gt -900) and (nfgsnodust.ehbha_err gt 0.0),nindxnfgs)

    xnfgs = nfgsdust[indxnfgs].b_lum_obs
    xerrnfgs = nfgsdust[indxnfgs].b_lum_obs_err
    
    ynfgs = nfgsnodust[indxnfgs].ehbha
    yerrnfgs = nfgsnodust[indxnfgs].ehbha_err

    xtitle = 'log L(B) [L(B)'+sunsymbol()+']'
    ytitle = 'E(H\beta-H\alpha) [mag]'

    xrange = LBrange
    yrange = ehbharange
    
    atlas1d_lineplot, x, y, xerr, yerr, plottype=3, postscript=postscript, ystyle=11, $
      ytitle=ytitle, xtitle=xtitle, xrange=xrange, yrange=yrange, legendtype=0, $
      xstyle=11, /right, /top, xnfgs=xnfgs, ynfgs=ynfgs, xerrnfgs=xerrnfgs, blackwhite=blackwhite, $
      yerrnfgs=yerrnfgs, position=pos[*,0], charsize=charsize_9, atlaspsize=1.1, nfgspsize=1.4
    axis, /xaxis, xrange=interpol(xabs,x,!x.crange), xthick=postthick, xstyle=1, $
      charsize=charsize_9, charthick=postthick
    axis, /yaxis, yrange=interpol(yebv,y,!y.crange), ythick=postthick, ystyle=1, $
      charsize=charsize_9, charthick=postthick
    im_xyouts_title, xtitle='M_{B} [mag]', charsize=charsize_9, charthick=postthick
    im_xyouts_title, ytitle='E(B-V) [mag]', charsize=charsize_9, charthick=postthick, xspacing=10.0

    im_openclose, postscript=postscript, /close

; ------------------------------------------------------------
; High-Redshift L([O II]) vs E(B-V)
; ------------------------------------------------------------

    psname = 'highz_loii_vs_ebv'
    im_openclose, pspath+psname, postscript=postscript, xsize=8.5, ysize=8.0, encapsulated=encapsulated

    pagemaker, nx=1, ny=1, xspace=0, yspace=0.0, width=6.6, height=6.6, $
      xmargin=[1.5,0.4], ymargin=[0.3,1.1], xpage=8.5, ypage=8.0, $
      position=pos, /normal

; SF galaxies    

    indx = where((atlasnodust.ebv_hahb_err gt 0.0) and (atlasnodust.oii_3727_lum[0] gt -900),nindx)
    
    x = atlasnodust[indx].oii_3727_lum[0] + alog10(lsun)
    xerr = atlasnodust[indx].oii_3727_lum[1]
    y = atlasnodust[indx].ebv_hahb
    yerr = atlasnodust[indx].ebv_hahb_err

    indxnfgs = where((nfgsnodust.ebv_hahb_err gt 0.0) and (nfgsnodust.oii_3727_lum[0] gt -900),nindxnfgs)

    xnfgs = nfgsnodust[indxnfgs].oii_3727_lum[0] + alog10(lsun)
    xerrnfgs = nfgsnodust[indxnfgs].oii_3727_lum[1]

    ynfgs = nfgsnodust[indxnfgs].ebv_hahb
    yerrnfgs = nfgsnodust[indxnfgs].ebv_hahb_err

    xbig = [x,xnfgs]
    ybig = [y,ynfgs]
    xerrbig = [xerr,xerrnfgs]
    yerrbig = [yerr,yerrnfgs]

; AGN

    indx_agn = where((atlasnodust_agn.ebv_hahb_err gt 0.0) and (atlasnodust_agn.oii_3727_lum[0] gt -900),nindx_agn)
    
    x_agn = atlasnodust_agn[indx_agn].oii_3727_lum[0] + alog10(lsun)
    xerr_agn = atlasnodust_agn[indx_agn].oii_3727_lum[1]
    y_agn = atlasnodust_agn[indx_agn].ebv_hahb
    yerr_agn = atlasnodust_agn[indx_agn].ebv_hahb_err

    indxnfgs_agn = where((nfgsnodust_agn.ebv_hahb_err gt 0.0) and $
      (nfgsnodust_agn.oii_3727_lum[0] gt -900),nindxnfgs_agn)

    xnfgs_agn = nfgsnodust_agn[indxnfgs_agn].oii_3727_lum[0] + alog10(lsun)
    xerrnfgs_agn = nfgsnodust_agn[indxnfgs_agn].oii_3727_lum[1]

    ynfgs_agn = nfgsnodust_agn[indxnfgs_agn].ebv_hahb
    yerrnfgs_agn = nfgsnodust_agn[indxnfgs_agn].ebv_hahb_err

    xbig_agn = [x_agn,xnfgs_agn]
    ybig_agn = [y_agn,ynfgs_agn]
    xerrbig_agn = [xerr_agn,xerrnfgs_agn]
    yerrbig_agn = [yerr_agn,yerrnfgs_agn]

; average error bar
    
    avgerrlocal = djs_median([yerrbig,yerrbig_agn])

; ####################
; Maier (2005)
; ####################
    
    indxm05 = where((m05nodust.ebv_hahb_err gt 0.0) and (m05nodust.oii_3727_lum[0] gt -900),nindxm05)

    xm05 = m05nodust[indxm05].oii_3727_lum[0] + alog10(lsun)
    xerrm05 = m05nodust[indxm05].oii_3727_lum[1]
    ym05 = m05nodust[indxm05].ebv_hahb
    yerrm05 = m05nodust[indxm05].ebv_hahb_err

    avgerrm05 = djs_median(yerrm05)
    splog, 'Maier: redshift = ', minmax(m05[indxm05].z_obj), median(m05[indxm05].z_obj)

; ####################
    
    xtitle = 'log L([O II])_{cor} [erg s^{-1}]'
    ytitle = 'E(B-V) [mag]'

    xrange = Loiikewleyrange
    yrange = ehbharange

    djs_plot, [0], [0], /nodata, xthick=postthick, ythick=postthick, charsize=charsize_8, $
      charthick=postthick, thick=postthick, xtitle=xtitle, ytitle=ytitle, $
      ymargin=[4,3], xstyle=3, ystyle=3, xrange=xrange, yrange=yrange, position=pos[*,0]

    im_symbols, localsym, psize=localpsize, /fill, color=djs_icolor(localcolor)
    djs_oplot, xbig, ybig, ps=8, color=localcolor
;   oploterror, xbig, ybig, xerrbig, yerrbig, ps=8, thick=postthick, errthick=highz_errthick, $
;     /nohat, color=djs_icolor(localcolor), errcolor=djs_icolor(localcolor)

    im_symbols, localsym, psize=localpsize, fill=0, color=djs_icolor(localcolor), thick=postthick
    djs_oplot, xbig_agn, ybig_agn, ps=8, color=localcolor
;   oploterror, xbig_agn, ybig_agn, xerrbig_agn, yerrbig_agn, ps=8, thick=postthick, errthick=highz_errthick, $
;     /nohat, color=djs_icolor(localcolor), errcolor=djs_icolor(localcolor)

    im_symbols, m05sym, psize=m05psize, /fill, color=djs_icolor(m05color), thick=postthick
;   djs_oplot, xm05, ym05, ps=8, color=m05color, thick=postthick
    oploterror, xm05, ym05, xerrm05, yerrm05, ps=8, thick=postthick, errthick=highz_errthick, $
      /nohat, color=djs_icolor(m05color), errcolor=djs_icolor(m05color)

    djs_oplot, !x.crange, [0,0], line=0, thick=postthick

; overplot the Kewley calibration

    xaxis = findgen((xrange[1]-xrange[0])/0.01+1)*0.01+xrange[0]
    yaxis = 0.174*xaxis - 6.84
    chop = where(yaxis gt 0.0)
    djs_oplot, xaxis[chop], yaxis[chop], line=2, thick=postthick;, color='dark green'

; legend

    label = ['Local Sample - Star-Forming','Local Sample - AGN','Maier et al. (2005)']
    psym = [localsym,localsym,m05sym]
    fill = [1,0,1]
    color = djs_icolor([localcolor,localcolor,m05color])
    postthick1 = postthick
    im_legend, textoidl(label), /left, /top, box=0, charsize=charsize_2, $
      charthick=postthick1, psym=psym, fill=fill, symsize=1.5, $
      spacing=1.8, thick=postthick1

    im_openclose, postscript=postscript, /close    

; ------------------------------------------------------------
; L(B) vs 12+log(O/H) [Integrated + SDSS]
; ------------------------------------------------------------

    psname = 'lb_vs_12oh_oiiinii_niiha_2panel'
    im_openclose, pspath+psname, postscript=postscript, xsize=8.5, ysize=5.5, encapsulated=encapsulated

    pagemaker, nx=2, ny=1, height=3.5, width=3.5, xmargin=[1.2,0.3], $
      ymargin=[0.9,1.1], xspace=0, yspace=0, xpage=8.5, ypage=5.5, $
      position=pos, /normal

; Integrated    
    
    indx = where((atlasdust.b_lum_obs gt -900.0) and (atlasdust.zstrong_12oh_oiiinii_niiha gt -900.0),nindx)

    x = atlasdust[indx].b_lum_obs
    xerr = atlasdust[indx].b_lum_obs_err
    xabs = atlasdust[indx].m_b_obs
    
    y = atlasdust[indx].zstrong_12oh_oiiinii_niiha
    yerr = atlasdust[indx].zstrong_12oh_oiiinii_niiha_err

    indxnfgs = where((nfgsdust.b_lum_obs gt -900) and (nfgsdust.zstrong_12oh_oiiinii_niiha gt -900.0),nindxnfgs)

    xnfgs = nfgsdust[indxnfgs].b_lum_obs
    xerrnfgs = nfgsdust[indxnfgs].b_lum_obs_err
    
    ynfgs = nfgsdust[indxnfgs].zstrong_12oh_oiiinii_niiha
    yerrnfgs = nfgsdust[indxnfgs].zstrong_12oh_oiiinii_niiha_err

    xtitle = 'log L(B) [L(B)'+sunsymbol()+']'
    ytitle = '12 + log (O/H)'

    xrange = LBrange
    yrange = ohrange2
    
    atlas1d_lineplot, x, y, xerr, yerr, plottype=3, postscript=postscript, ysty=3, $
      ytitle=ytitle, xtitle=xtitle, xrange=xrange, yrange=yrange, legendtype=0, $
      xstyle=11, /right, /top, xnfgs=xnfgs, ynfgs=ynfgs, xerrnfgs=xerrnfgs, blackwhite=blackwhite, $
      yerrnfgs=yerrnfgs, position=pos[*,0], charsize=charsize_6
    axis, /xaxis, xrange=interpol(xabs,x,!x.crange), xthick=postthick, xstyle=1, $
      charsize=charsize_6, charthick=postthick
    im_xyouts_title, xtitle='M_{B} [mag]', charsize=charsize_6, charthick=postthick

; SDSS    
    
    indx = where((sdssancillary.b_lum gt -900) and (sdssdust.zstrong_12oh_oiiinii_niiha gt -900.0),nindx)

    x = sdssancillary[indx].b_lum
    xabs = sdssancillary[indx].m_b
    xerr = sdssancillary[indx].b_lum_err
    
    y = sdssdust[indx].zstrong_12oh_oiiinii_niiha
    yerr = sdssdust[indx].zstrong_12oh_oiiinii_niiha_err

    sdss_lineplot, x, y, xerr, yerr, plottype=3, postscript=postscript, $
      ytitle='', xtitle=xtitle, xrange=xrange, yrange=yrange, legendtype=0, $
      xstyle=11, /right, /top, position=pos[*,1], charsize=charsize_6, /noerase, $
      ytickname=replicate(' ',10)
    axis, /xaxis, xrange=interpol(xabs,x,!x.crange), xthick=postthick, xstyle=1, $
      charsize=charsize_6, charthick=postthick
    im_xyouts_title, xtitle='M_{B} [mag]', charsize=charsize_6, charthick=postthick

    im_openclose, postscript=postscript, /close

; ------------------------------------------------------------
; 2-panel - SFR(IR) vs SFR(Ha)/SFR(IR)
; ------------------------------------------------------------
    
    psname = 'sfr_ir_vs_sfr_hair_2panel'
    im_openclose, pspath+psname, postscript=postscript, xsize=8.5, ysize=7.5, encapsulated=encapsulated

    pagemaker, nx=1, ny=2, height=[3.0,3.0], width=7.0, xmargin=[1.2,0.3], $
      ymargin=[0.4,1.1], xspace=0, yspace=0, xpage=8.5, ypage=7.5, $
      position=pos, /normal

    indx = where((atlasnodust.ir_flux gt -900) and (atlasnodust.ebv_hahb_err gt 0.0),nindx)

    x = atlasdust[indx].ir_lum + alog10(lsun) + alog10(irconst)
    xerr = atlasdust[indx].ir_lum_err
    x = atlasdust[indx].ir_lum + alog10(lsun) + alog10(irconst)
    
    lir = atlasdust[indx].ir_flux ; [erg/s/cm2]
    lir_err = sqrt(atlasdust[indx].ir_flux_err^2 + (atlasdust[indx].ir_flux*0.15)^2)

    indxnfgs = where((nfgsdust.ir_flux gt -900) and (nfgsnodust.ebv_hahb_err gt 0.0),nindxnfgs)

    xnfgs = nfgsnodust[indxnfgs].ir_lum + alog10(lsun) + alog10(irconst)
    xerrnfgs = nfgsnodust[indxnfgs].ir_lum_err
    
    lirnfgs = nfgsdust[indxnfgs].ir_flux
    lirnfgs_err = sqrt(nfgsdust[indxnfgs].ir_flux_err^2 + (nfgsdust[indxnfgs].ir_flux*0.15)^2)
    
; flag the LIRGS and ULIRGS

    lirgs = where(atlasdust[indx].ir_lum gt 11.0,comp=notlirgs)
    lirgsnfgs = where(nfgsdust[indxnfgs].ir_lum gt 11.0,notlirgsnfgs)

; plotting variables    
    
    xrange = sfrIRrange
    yrange = sfrHaIRrange

    xtitle = 'log \psi(IR) [M'+sunsymbol()+' yr^{-1}]'
    ytitle = 'log \psi(H\alpha)/\psi(IR)'

    lhalir = alog10(4.5D-44/7.9D-42)

; ##########################################################
; Panel 1: No reddening correction
; ##########################################################

    ha = atlasdust[indx].h_alpha[0]
    ha_err = atlasdust[indx].h_alpha[1]
    
    y = alog10(ha/lir) - lhalir
    yerr = im_compute_error(ha,ha_err,lir,lir_err,/log)

    hanfgs = nfgsdust[indxnfgs].h_alpha[0]
    hanfgs_err = nfgsdust[indxnfgs].h_alpha[1]
    
    ynfgs = alog10(hanfgs/lirnfgs) - lhalir
    yerrnfgs = im_compute_error(hanfgs,hanfgs_err,lirnfgs,lirnfgs_err,/log)

    stats = im_stats([y,ynfgs],/verbose)
    xstr = strtrim(string(stats.median_rej,format='(F12.2)'),2)+' ('+$
      strtrim(string(stats.mean_rej,format='(F12.2)'),2)+$
      '\pm'+strtrim(string(stats.sigma_rej,format='(F12.2)'),2)+')'
;   xstr = strtrim(string(stats.median,format='(F12.2)'),2)+'\pm'+strtrim(string(stats.sig68mean,format='(F12.2)'),2)

    ytitle = 'log \psi(H\alpha)/\psi(IR)'

    atlas1d_lineplot, x, y, xerr, yerr, plottype=3, postscript=postscript, $
      xtitle='', ytitle=ytitle, xrange=xrange, yrange=yrange, legendtype=0, $
      charsize=charsize_7, position=pos[*,0], xnfgs=xnfgs, ynfgs=ynfgs, $
      xerrnfgs=xerrnfgs, yerrnfgs=yerrnfgs, xtickname=replicate(' ',10), $
      yfrac=8, blackwhite=blackwhite
    djs_oplot, !x.crange, [0,0], line=0, thick=postthick
;   djs_oplot, lir_bell, lratio_bell, line=2, thick=postthick

;   plotsym, 0, 1.6, thick=(postthick-3L)>1L
;   djs_oplot, x[lirgs], y[lirgs], ps=8
;   djs_oplot, xnfgs[lirgsnfgs], ynfgs[lirgsnfgs], ps=8

    legend, '(a)', /left, /top, box=0, charsize=charsize_3, charthick=postthick
    legend, textoidl('Observed H\alpha'), /left, /bottom, box=0, charsize=charsize_3, charthick=postthick
    legend, textoidl(xstr), /right, /top, box=0, charsize=charsize_3, charthick=postthick
;   label = ['A(H\alpha) = 0']
;   legend, textoidl(label), /left, /bottom, box=0, charsize=charsize_3, $
;     charthick=postthick

    w = where([x,xnfgs] ge 9.0 and [x,xnfgs] le 11.0)
;   junk = im_stats(([y,ynfgs])[w],/verbose)

; ############################################################
; Panel 2: Good Reddening Correction
; ############################################################

    ha = atlasnodust[indx].h_alpha[0]
    ha_err = atlasnodust[indx].h_alpha[1]

    y = alog10(ha/lir) - lhalir
    yerr = im_compute_error(ha,ha_err,lir,lir_err,/log)

    hanfgs = nfgsnodust[indxnfgs].h_alpha[0]
    hanfgs_err = nfgsnodust[indxnfgs].h_alpha[1]

    ynfgs = alog10(hanfgs/lirnfgs) - lhalir
    yerrnfgs = im_compute_error(hanfgs,hanfgs_err,lirnfgs,lirnfgs_err,/log)

    stats = im_stats([y,ynfgs],/verbose,/no_head)
    xstr = strtrim(string(stats.median_rej,format='(F12.2)'),2)+' ('+$
      strtrim(string(stats.mean_rej,format='(F12.2)'),2)+$
      '\pm'+strtrim(string(stats.sigma_rej,format='(F12.2)'),2)+')'
;   xstr = strtrim(string(stats.median,format='(F12.2)'),2)+'\pm'+strtrim(string(stats.sig68mean,format='(F12.2)'),2)

    ytitle = 'log \psi(H\alpha)/\psi(IR)'

    atlas1d_lineplot, x, y, xerr, yerr, plottype=3, postscript=postscript, $
      xtitle=xtitle, ytitle=ytitle, xrange=xrange, yrange=yrange, legendtype=0, $
      charsize=charsize_7, position=pos[*,1], /noerase, xnfgs=xnfgs, ynfgs=ynfgs, $
      xerrnfgs=xerrnfgs, yerrnfgs=yerrnfgs, blackwhite=blackwhite
    djs_oplot, !x.crange, [0,0], line=0, thick=postthick
;   djs_oplot, lir_bell, lratio_bell, line=2, thick=postthick

;   plotsym, 0, 1.6, thick=(postthick-3L)>1L
;   djs_oplot, x[lirgs], y[lirgs], ps=8
;   djs_oplot, xnfgs[lirgsnfgs], ynfgs[lirgsnfgs], ps=8

    legend, '(b)', /left, /top, box=0, charsize=charsize_3, charthick=postthick
    legend, textoidl('Extinction-Corrected H\alpha'), /left, /bottom, box=0, charsize=charsize_3, charthick=postthick
    legend, textoidl(xstr), /right, /top, box=0, charsize=charsize_3, charthick=postthick
;   label = ['Individual A(H\alpha)']
;   legend, textoidl(label), /left, /bottom, box=0, charsize=charsize_3, $
;     charthick=postthick

    w = where([x,xnfgs] ge 9.0 and [x,xnfgs] le 11.0)
;   junk = im_stats(([y,ynfgs])[w],/verbose)

; overplot the prediction by Hirashita et al. (2003):

    eta = 0.4
    epsilon = 0.5
    frac = 0.6

;   hbi_lir = 10^(findgen((12.5-7.0)/0.05+1)*0.05+7.0) ; Hirashita, Buat, & Inoue
;   hbi_lha = lsun*10^(findgen((9.0-5.0)/0.05+1)*0.05+5.0)

    hbi_lir = 10^atlasdust[indx].ir_lum
    srt = sort(hbi_lir)
    hbi_lir = hbi_lir[srt]
    hbi_lha = lsun*10^atlasdust[indx].h_alpha_lum[0]
    hbi_lha = hbi_lha[srt]
    
    ir_sfr = 1.79E-10*(1-eta) / (0.13 - 0.085*frac + 0.87*epsilon) * hbi_lir
    ha_sfr = 7.89D-42 / frac * hbi_lha

;   oplot, alog10(ir_sfr), alog10(ha_sfr/ir_sfr), line=2, thick=postthick
    
    im_openclose, postscript=postscript, /close

; ------------------------------------------------------------
; L(B) versus SFR(Ha)/L(Hb)_obs - SDSS
; ------------------------------------------------------------

    psname = 'sdss_lb_vs_sfr_ha_lhb_obs'
    im_openclose, pspath+psname, postscript=postscript, xsize=8.5, ysize=8.8, encapsulated=encapsulated

    pagemaker, nx=1, ny=1, xspace=0, yspace=0.0, width=6.8, height=6.8, $
      xmargin=[1.4,0.3], ymargin=[0.9,1.1], xpage=8.5, ypage=8.8, $
      position=pos, /normal

    medbin = 0.5
    minpts = 3
    minx = 8.5 ; -3.0
    medcolor = 'orange'
    
    indx = where((sdssancillary.b_lum gt -900) and (sdssnodust.sfr_h_alpha gt -900.0) and $
      (sdssdust.h_beta_ew_uncor[0] ge ewcut),nindx)

    x = sdssancillary[indx].b_lum; - alog10(LBnorm)
    xerr = sdssancillary[indx].b_lum_err
    xabs = sdssancillary[indx].m_b

    hb = sdssdust[indx].h_beta_lum[0] + alog10(lsun) - alog10(elumnorm)
    hb_err = sdssdust[indx].h_beta_lum[1]

    sfr_ha = sdssnodust[indx].sfr_h_alpha
    sfr_ha_err = sdssnodust[indx].sfr_h_alpha_err

    y = sfr_ha - hb
    yerr = sqrt(sfr_ha_err^2 + hb_err^2)
    
; combine the two samples

    xbig = x & xerrbig = xerr
    ybig = y & yerrbig = yerr

; make the plot

    xtitle = 'log L(B) [L(B)'+sunsymbol()+']'
    ytitle = 'log [10^{41} \psi/L(H\beta)_{obs}] [erg s^{-1}/'+sfr_units()+']'

    xrange = LBrange; - alog10(LBnorm)
    yrange = alog10([1.8D-41,3D-40]) - alog10(sfrnorm) + alog10(elumnorm)

    sdss_lineplot, xbig, ybig, xerrbig, yerrbig, $
      plottype=1, postscript=postscript, xtitle=xtitle, ytitle=ytitle, xrange=xrange, $
      xstyle=11, yrange=yrange, legendtype=0, /right, /top, yfrac=1.5, xfrac=1.5, charsize=charsize_8, $
      position=pos[*,0]
    axis, /xaxis, xrange=interpol(xabs,x,!x.crange), xthick=postthick, xstyle=1, $
      charsize=charsize_8, charthick=postthick
    im_xyouts_title, xtitle='M_{B} [mag]', charsize=charsize_8, charthick=postthick

; fit only to certain data

    xfit = xbig
    yfit = ybig
    
; 9.0 is M_B = -17.07 and 8.7 is M_B = -16.3; 8.57 is -16

    LBlocut = 10.0^9.0 ; 10.0^8.7
    LBhicut = 10.0^11.0

    LBlo = alog10(LBlocut)>min(xfit)
    LBhi = alog10(LBhicut)<max(xfit)
;   LBlo = alog10(LBlocut/LBnorm)>min(xfit)
;   LBhi = alog10(LBhicut/LBnorm)<max(xfit)

; overlay the running median    

    running = im_medxbin(xfit,yfit,medbin,minx=minx,minpts=minpts,/verbose)

    doit = where((running.binctr gt LBlo) and (running.binctr le LBhi),ndoit)
    lolum = where((running.binctr lt LBlo))

    im_symbols, 108, psize=3.0, thick=postthick, color=medcolor;, /fill
;   oploterror, running.binctr[doit], running.medy[doit], $
;     running.sigy[doit], ps=8, errthick=(postthick-2L)>2L
    oploterror, running.binctr[doit], running.medy[doit], running.binsz/2.0+fltarr(ndoit), $
      running.sigy75[doit]-running.medy[doit], ps=8, errthick=(postthick)>2L, /hibar, $
      color=djs_icolor(medcolor), errcolor=djs_icolor(medcolor)
    oploterror, running.binctr[doit], running.medy[doit], running.binsz/2.0+fltarr(ndoit), $
      running.medy[doit]-running.sigy25[doit], ps=8, errthick=(postthick)>2L, /lobar, $
      color=djs_icolor(medcolor), errcolor=djs_icolor(medcolor)
;   oploterror, running.binctr[doit], running.medy[doit], running.binsz/2.0+fltarr(ndoit), $
;     running.sigyup[doit]-running.medy[doit], ps=8, errthick=(postthick-1L)>2L, /hibar, $
;     color=djs_icolor(medcolor), errcolor=djs_icolor(medcolor)
;   oploterror, running.binctr[doit], running.medy[doit], running.binsz/2.0+fltarr(ndoit), $
;     running.medy[doit]-running.sigylo[doit], ps=8, errthick=(postthick-1L)>2L, /lobar, $
;     color=djs_icolor(medcolor), errcolor=djs_icolor(medcolor)
;   oploterror, running.binctr[doit], running.medy[doit], $
;     running.sigyup[doit]-running.medy[doit], ps=8, errthick=(postthick-1L)>2L, /hibar
;   oploterror, running.binctr[doit], running.medy[doit], $
;     running.medy[doit]-running.sigylo[doit], ps=8, errthick=(postthick-1L)>2L, /lobar

    if keyword_set(paper) then begin
    
; generate a binary FITS table with the statistics

       hb_sfrfile = 'sdss_hb_sfr.fits'

       hb_sfr = replicate(template_sfr,ndoit)
       hb_sfr.loglb  = running.binctr[doit]
       hb_sfr.mb     = interpol(xabs,x,running.binctr[doit])
       hb_sfr.p25    = running.sigy25[doit]
       hb_sfr.p50    = running.medy[doit]
       hb_sfr.p75    = running.sigy75[doit]
       hb_sfr.mean   = running.meany[doit]
       hb_sfr.stddev = running.stddev[doit]
       
       mwrfits, hb_sfr, sfrspath+hb_sfrfile, /create

; generate a latex table with the statistics

       openw, lun, latexpath+'sdss_hb_sfr.tex', /get_lun
       printf, lun, '\begin{deluxetable}{ccccccc}'
       printf, lun, '\tabletypesize{\small}'
       printf, lun, '\tablecolumns{7}'
       printf, lun, '\tablecaption{SDSS \hblam{} Star-Formation Rate Conversion Factors \label{table:hb_sfr}}'
       printf, lun, '\tablewidth{0in}'
       printf, lun, '\tablehead{'
       printf, lun, '\colhead{$\log\,\lb$} & '
       printf, lun, '\colhead{\mb} & '
       printf, lun, '\multicolumn{5}{c}{$\log\,[\sfr/\lhbobs]$} \\'

       printf, lun, '\colhead{[\lbsun]} & '
       printf, lun, '\colhead{[mag]} & '
       printf, lun, '\multicolumn{5}{c}{$[10^{41}\,\lunits/(\sfrunits)]$} \\'

       printf, lun, '\cline{3-7}'

       printf, lun, '\multicolumn{2}{c}{} & '
       printf, lun, '\colhead{$P_{25}$} & '
       printf, lun, '\colhead{$P_{50}$} & '
       printf, lun, '\colhead{$P_{75}$} & '
       printf, lun, '\colhead{$\langle R\rangle$} & '
       printf, lun, '\colhead{$\sigma_{R}$}'
       printf, lun, '}'
       printf, lun, '\startdata'

       for i = 0L, ndoit-1L do printf, lun, $
         running.binctr[doit[i]], ' & ', interpol(xabs,x,running.binctr[doit[i]]), ' & ', $
         running.sigy25[doit[i]], ' & ', running.medy[doit[i]], ' & ', $
         running.sigy75[doit[i]], ' & ', running.meany[doit[i]], ' & ', $
         running.stddev[doit[i]], ' \\', $
         format='(F5.2,A3,F6.2,A3,F5.3,A3,F5.3,A3,F5.3,A3,F5.3,A3,F5.3,A3)'

       printf, lun, '\enddata'
       printf, lun, '%\tablenotetext{a}{}'
       printf, lun, '\tablecomments{The columns labeled $P_{25}$, $P_{50}$, and $P_{75}$ '+$
         'give the $25$, $50$ (median), and $75$ percentile of the $\log\,[\sfr/\lhbobs]$ '+$
         'distribution, respectively, in bins of $0.5$~dex in luminosity.  $\langle R\rangle$ '+$
         'and $\sigma_{R}$ give the mean and standard-deviation of the distribution in each bin.}'
       printf, lun, '%\tablerefs{}'
       printf, lun, '\end{deluxetable}'
       free_lun, lun

    endif
    
    im_openclose, postscript=postscript, /close    

; ------------------------------------------------------------
; 3-panel - L(B) vs [O III]/Ha - Integrated
; ------------------------------------------------------------

    psname = 'lb_vs_oiiiha_3panel'
    im_openclose, pspath+psname, postscript=postscript, xsize=8.5, ysize=9.5, encapsulated=encapsulated

    pagemaker, nx=1, ny=3, height=[2.5,2.5,2.5], width=7.0, xmargin=[1.1,0.4], $
      ymargin=[0.9,1.1], xspace=0, yspace=0, xpage=8.5, ypage=9.5, $
      position=pos, /normal

    medbin = 0.5
    minpts = 10
    minx = 7.75

    indx = where((atlasdust.oiii_5007[0]/atlasdust.oiii_5007[1] gt snrcut) and $
      (atlasnodust.b_lum_obs gt -900.0) and (atlasnodust.ehbha_err gt 0.0),nindx)

    x = atlasnodust[indx].b_lum_obs
    xerr = atlasnodust[indx].b_lum_obs_err
    xabs = atlasnodust[indx].m_b_obs

    indxnfgs = where((nfgsdust.oiii_5007[0]/nfgsdust.oiii_5007[1] gt snrcut) and $
      (nfgsnodust.b_lum_obs gt -900.0) and (nfgsnodust.ehbha_err gt 0.0),nindxnfgs)

    xnfgs = nfgsnodust[indxnfgs].b_lum_obs
    xerrnfgs = nfgsnodust[indxnfgs].b_lum_obs_err

    xtitle = 'log L(B) [L(B)'+sunsymbol()+']'
;   xtitle = 'log [L(B)/L(B,'+im_sunsymbol()+')]'
    ytitle = 'log ([O III]/H\alpha)'

    xrange = LBrange
    yrange = oiiihacorrange

; ##########################################################
; Panel 1: [O III]_obs, Ha_obs
; ##########################################################

    y1 = atlasdust[indx].oiii_5007[0]
    y1err = atlasdust[indx].oiii_5007[1]

    y2 = atlasdust[indx].h_alpha[0]
    y2err = atlasdust[indx].h_alpha[1]

    y = alog10(y1/y2)
    yerr = im_compute_error(y1,y1err,y2,y2err,/log)

; NFGS    
    
    y1nfgs = nfgsdust[indxnfgs].oiii_5007[0]
    y1errnfgs = nfgsdust[indxnfgs].oiii_5007[1]

    y2nfgs = nfgsdust[indxnfgs].h_alpha[0]
    y2errnfgs = nfgsdust[indxnfgs].h_alpha[1]

    ynfgs = alog10(y1nfgs/y2nfgs)
    yerrnfgs = im_compute_error(y1nfgs,y1errnfgs,y2nfgs,y2errnfgs,/log)

    stats = im_stats([y,ynfgs],/verbose,/baremin)
    xstr = strtrim(string(stats.median_rej,format='(F12.2)'),2)+' ('+$
      strtrim(string(stats.mean_rej,format='(F12.2)'),2)+$
      '\pm'+strtrim(string(stats.sigma_rej,format='(F12.2)'),2)+')'
;   xstr = strtrim(string(stats.median,format='(F12.2)'),2)+'\pm'+strtrim(string(stats.sig68mean,format='(F12.2)'),2)
    rcor = r_correlate([x,xnfgs],[y,ynfgs],zd=zd,probd=probd)
    splog, 'Spearman rank test: ', rcor, probd, zd

    ytitle = 'log ([O III]/H\alpha)_{obs}'

    atlas1d_lineplot, x, y, xerr, yerr, plottype=3, postscript=postscript, $
      ytitle=ytitle, xrange=xrange, yrange=yrange, legendtype=0, $
      xstyle=11, /right, /top, xnfgs=xnfgs, ynfgs=ynfgs, xerrnfgs=xerrnfgs, $
      yerrnfgs=yerrnfgs, position=pos[*,0], charsize=charsize_4, xtickname=replicate(' ',10)
    axis, /xaxis, xrange=interpol(xabs,x,!x.crange), xthick=postthick, xstyle=1, $
      charsize=charsize_4, charthick=postthick
    im_xyouts_title, xtitle='M_{B} [mag]', charsize=charsize_4, charthick=postthick

;   print & running = im_medxbin([x,xnfgs],[y,ynfgs],medbin,minx=minx,minpts=minpts,/verbose) & print
;   oplot, running.binctr, running.medy, thick=postthick, line=0
;   oplot, running.binctr, running.medy+running.sigy, thick=postthick, line=2
;   oplot, running.binctr, running.medy-running.sigy, thick=postthick, line=2
;   oploterror, running.binctr, running.medy, replicate(running.binsz/2.0,running.nbins), $
;     running.sigy, ps=-3, thick=postthick, errthick=postthick, /nohat, line=0, errstyle=0
;   niceprint, running.binctr, running.medy, running.sigy

    legend, '(a)', /left, /top, box=0, charsize=charsize_3, charthick=postthick
    legend, textoidl(xstr), /right, /top, box=0, charsize=charsize_3, charthick=postthick

    w = where([x,xnfgs] gt 9.5)
;   junk = im_stats(([y,ynfgs])[w],/ver)

; ##########################################################
; Panel 2: [O III]_cor, Ha_cor
; ##########################################################

    y1 = atlasnodust[indx].oiii_5007[0]
    y1err = atlasnodust[indx].oiii_5007[1]
 
    y2 = atlasnodust[indx].h_alpha[0]
    y2err = atlasnodust[indx].h_alpha[1]
 
    y = alog10(y1/y2)
    yerr = im_compute_error(y1,y1err,y2,y2err,/log)
 
    y1nfgs = nfgsnodust[indxnfgs].oiii_5007[0]
    y1errnfgs = nfgsnodust[indxnfgs].oiii_5007[1]
 
    y2nfgs = nfgsnodust[indxnfgs].h_alpha[0]
    y2errnfgs = nfgsnodust[indxnfgs].h_alpha[1]
 
    ynfgs = alog10(y1nfgs/y2nfgs)
    yerrnfgs = im_compute_error(y1nfgs,y1errnfgs,y2nfgs,y2errnfgs,/log)

;   y1 = atlasnodust[indx].oiii_5007_lum[0]
;   y1err = atlasnodust[indx].oiii_5007_lum[1]
;
;   y2 = atlasnodust[indx].h_alpha_lum[0]
;   y2err = atlasnodust[indx].h_alpha_lum[1]
;
;   y = y1 - y2 - loghaconst - hasfrconstoffset
;   yerr = sqrt(y1err^2 + y2err^2)
;
;   y1nfgs = nfgsnodust[indxnfgs].oiii_5007_lum[0]
;   y1errnfgs = nfgsnodust[indxnfgs].oiii_5007_lum[1]
;
;   y2nfgs = nfgsnodust[indxnfgs].h_alpha_lum[0]
;   y2errnfgs = nfgsnodust[indxnfgs].h_alpha_lum[1]
;
;   ynfgs = y1nfgs - y2nfgs - loghaconst - hasfrconstoffset
;   yerrnfgs = sqrt(y1errnfgs^2 + y2errnfgs^2)
    
    stats = im_stats([y,ynfgs],/verbose,/baremin,/no_head)
    xstr = strtrim(string(stats.median_rej,format='(F12.2)'),2)+' ('+$
      strtrim(string(stats.mean_rej,format='(F12.2)'),2)+$
      '\pm'+strtrim(string(stats.sigma_rej,format='(F12.2)'),2)+')'
;   xstr = strtrim(string(stats.median,format='(F12.2)'),2)+'\pm'+strtrim(string(stats.sig68mean,format='(F12.2)'),2)
    rcor = r_correlate([x,xnfgs],[y,ynfgs],zd=zd,probd=probd)
    splog, 'Spearman rank test: ', rcor, probd, zd

    ytitle = 'log ([O III]/H\alpha)_{cor}'
;   ytitle = 'log [L([O III])_{cor}/\psi(H\alpha)]-'+string(hasfrconstoffset,format='(I0)')
;   yrange2 = sfroiiharange - hasfrconstoffset

    atlas1d_lineplot, x, y, xerr, yerr, plottype=3, postscript=postscript, $
      xtitle='', ytitle=ytitle, xrange=xrange, yrange=yrange, legendtype=0, /noerase, $
      xstyle=3, /right, /top, xtickname=replicate(' ',10), $
      xnfgs=xnfgs, ynfgs=ynfgs, xerrnfgs=xerrnfgs, $
      yerrnfgs=yerrnfgs, position=pos[*,1], charsize=charsize_4
    
    print & running = im_medxbin([x,xnfgs],[y,ynfgs],medbin,minx=minx,minpts=minpts,/verbose) & print
;   oplot, running.binctr, running.medy, thick=postthick, line=0
;   oplot, running.binctr, running.medy+running.sigy, thick=postthick, line=2
;   oplot, running.binctr, running.medy-running.sigy, thick=postthick, line=2
;   oploterror, running.binctr, running.medy, replicate(running.binsz/2.0,running.nbins), $
;     running.sigy, ps=-3, thick=postthick, errthick=postthick, /nohat, line=0, errstyle=0

    legend, '(b)', /left, /top, box=0, charsize=charsize_3, charthick=postthick
    legend, textoidl(xstr), /right, /top, box=0, charsize=charsize_3, charthick=postthick

    w = where([x,xnfgs] gt 9.5)
;   junk = im_stats(([y,ynfgs])[w],/ver)

; ##########################################################
; Panel 3: [O III]_obs, SFR(Ha)
; ##########################################################

    y1 = atlasdust[indx].oiii_5007_lum[0]
    y1err = atlasdust[indx].oiii_5007_lum[1]

    y2 = atlasnodust[indx].h_alpha_lum[0]
    y2err = atlasnodust[indx].h_alpha_lum[1]

    y = y1 - y2 - loghaconst - hasfrconstoffset
    yerr = sqrt(y1err^2 + y2err^2)

; NFGS    
    
    y1nfgs = nfgsdust[indxnfgs].oiii_5007_lum[0]
    y1errnfgs = nfgsdust[indxnfgs].oiii_5007_lum[1]

    y2nfgs = nfgsnodust[indxnfgs].h_alpha_lum[0]
    y2errnfgs = nfgsnodust[indxnfgs].h_alpha_lum[1]

    ynfgs = y1nfgs - y2nfgs - loghaconst - hasfrconstoffset
    yerrnfgs = sqrt(y1errnfgs^2 + y2errnfgs^2)
    
;   y1 = atlasdust[indx].oiii_5007[0]
;   y1err = atlasdust[indx].oiii_5007[1]
;
;   y2 = atlasnodust[indx].h_alpha[0]
;   y2err = atlasnodust[indx].h_alpha[1]
;
;   y = alog10(y1/y2)
;   yerr = im_compute_error(y1,y1err,y2,y2err,/log)
;
;   y1nfgs = nfgsdust[indxnfgs].oiii_5007[0]
;   y1errnfgs = nfgsdust[indxnfgs].oiii_5007[1]
;
;   y2nfgs = nfgsnodust[indxnfgs].h_alpha[0]
;   y2errnfgs = nfgsnodust[indxnfgs].h_alpha[1]
;
;   ynfgs = alog10(y1nfgs/y2nfgs)
;   yerrnfgs = im_compute_error(y1nfgs,y1errnfgs,y2nfgs,y2errnfgs,/log)

    stats = im_stats([y,ynfgs],/verbose,/baremin,/no_head)
    xstr = strtrim(string(stats.median_rej,format='(F12.2)'),2)+' ('+$
      strtrim(string(stats.mean_rej,format='(F12.2)'),2)+$
      '\pm'+strtrim(string(stats.sigma_rej,format='(F12.2)'),2)+')'
;   xstr = strtrim(string(stats.sig68mean,format='(F12.2)'),2)
;   xstr = strtrim(string(stats.median,format='(F12.2)'),2)+'\pm'+strtrim(string(stats.sig68mean,format='(F12.2)'),2)
    rcor = r_correlate([x,xnfgs],[y,ynfgs],zd=zd,probd=probd)
    splog, 'Spearman rank test: ', rcor, probd, zd

;   ytitle = 'log ([O III]_{obs}/H\alpha_{cor})'
    ytitle = 'log [10^{-41} L([O III])_{obs}/\psi(H\alpha)]'
    yrange2 = sfroiiharange - hasfrconstoffset

    atlas1d_lineplot, x, y, xerr, yerr, plottype=3, postscript=postscript, $
      xtitle=xtitle, ytitle=ytitle, xrange=xrange, yrange=yrange2, legendtype=0, /noerase, $
      xstyle=3, /right, /top, $
      xnfgs=xnfgs, ynfgs=ynfgs, xerrnfgs=xerrnfgs, $
      yerrnfgs=yerrnfgs, position=pos[*,2], charsize=charsize_4
;   djs_oplot, !x.crange, alog10(1.0)*[1,1], line=0, thick=(postthick-4)>2

    print & running = im_medxbin([x,xnfgs],[y,ynfgs],medbin,minx=minx,minpts=minpts,/verbose) & print
;   oplot, running.binctr, running.medy, thick=postthick, line=0
;   oplot, running.binctr, running.medy+running.sigy, thick=postthick, line=2
;   oplot, running.binctr, running.medy-running.sigy, thick=postthick, line=2
;   oploterror, running.binctr, running.medy, replicate(running.binsz/2.0,running.nbins), $
;     running.sigy, ps=-3, thick=postthick, errthick=postthick, /nohat, line=0, errstyle=0

    legend, '(c)', /left, /top, box=0, charsize=charsize_3, charthick=postthick
    legend, textoidl(xstr), /right, /top, box=0, charsize=charsize_3, charthick=postthick

;   w = where([x,xnfgs] lt 9)
    w2 = where([x,xnfgs] le 9.5)
    w3 = where([x,xnfgs] gt 9.5)
;   junk = im_stats(([y,ynfgs])[w],/ver)
    junk = im_stats(([y,ynfgs])[w2],/ver)
    junk = im_stats(([y,ynfgs])[w3],/ver)

;   w = where([x,xnfgs] ge 10.25 and [x,xnfgs] le 10.75)
;   junk = im_stats(([y,ynfgs])[w],/verbose)

    im_openclose, postscript=postscript, /close

; ------------------------------------------------------------
; 3-panel - L(B) vs [O III]/Ha - SDSS
; ------------------------------------------------------------

    psname = 'sdss_lb_vs_oiiiha_3panel'
    im_openclose, pspath+psname, postscript=postscript, xsize=8.5, ysize=9.5, encapsulated=encapsulated

    pagemaker, nx=1, ny=3, height=[2.5,2.5,2.5], width=7.0, xmargin=[1.1,0.4], $
      ymargin=[0.9,1.1], xspace=0, yspace=0, xpage=8.5, ypage=9.5, $
      position=pos, /normal

    indx = where((sdssdust.oiii_5007[0]/sdssdust.oiii_5007[1] gt snrcut) and $
      (sdssancillary.b_lum gt -900.0) and (sdssnodust.ehbha_err gt 0),nindx)

    x = sdssancillary[indx].b_lum
    xerr = sdssancillary[indx].b_lum_err
    xabs = sdssancillary[indx].m_b

    xtitle = 'log L(B) [L(B)'+sunsymbol()+']'
    ytitle = 'log ([O III]/H\alpha)'

    xrange = LBrange
    yrange = oiiihacorrange

; ##########################################################
; Panel 1: [O III]_obs, Ha_obs
; ##########################################################

    y1 = sdssdust[indx].oiii_5007[0]
    y1err = sdssdust[indx].oiii_5007[1]

    y2 = sdssdust[indx].h_alpha[0]
    y2err = sdssdust[indx].h_alpha[1]

    y = alog10(y1/y2)
    yerr = im_compute_error(y1,y1err,y2,y2err,/log)

    stats = im_stats(y,/verbose,/baremin)
    xstr = strtrim(string(stats.median_rej,format='(F12.2)'),2)+' ('+$
      strtrim(string(stats.mean_rej,format='(F12.2)'),2)+$
      '\pm'+strtrim(string(stats.sigma_rej,format='(F12.2)'),2)+')'
;   xstr = strtrim(string(stats.median,format='(F12.2)'),2)+'\pm'+strtrim(string(stats.sig68mean,format='(F12.2)'),2)

    ytitle = 'log ([O III]/H\alpha)_{obs}'

    sdss_lineplot, x, y, xerr, yerr, plottype=3, postscript=postscript, $
      ytitle=ytitle, xrange=xrange, yrange=yrange, legendtype=0, $
      xstyle=11, /right, /top, position=pos[*,0], charsize=charsize_4, xtickname=replicate(' ',10)
    axis, /xaxis, xrange=interpol(xabs,x,!x.crange), xthick=postthick, xstyle=1, $
      charsize=charsize_4, charthick=postthick
    im_xyouts_title, xtitle='M_{B} [mag]', charsize=charsize_4, charthick=postthick
;   xyouts, mean(!x.crange), !y.crange[1]+0.16*(!y.crange[1]-!y.crange[0]), textoidl('M_{B} [mag]'), $
;     align=0.5, charsize=charsize_4, charthick=postthick

;   print & running = im_medxbin([x,xnfgs],[y,ynfgs],medbin,minx=minx,minpts=minpts,/verbose) & print
;   oploterror, running.binctr, running.medy, replicate(running.binsz/2.0,running.nbins), $
;     running.sigy, ps=-3, thick=postthick, errthick=postthick, /nohat, line=0, errstyle=0

    legend, '(a)', /left, /top, box=0, charsize=charsize_3, charthick=postthick
    legend, textoidl(xstr), /right, /top, box=0, charsize=charsize_3, charthick=postthick

; ##########################################################
; Panel 2: [O III]_cor, Ha_cor
; ##########################################################

    y1 = sdssnodust[indx].oiii_5007[0]
    y1err = sdssnodust[indx].oiii_5007[1]

    y2 = sdssnodust[indx].h_alpha[0]
    y2err = sdssnodust[indx].h_alpha[1]

    y = alog10(y1/y2)
    yerr = im_compute_error(y1,y1err,y2,y2err,/log)

    stats = im_stats(y,/verbose,/baremin,/no_head)
    xstr = strtrim(string(stats.median_rej,format='(F12.2)'),2)+' ('+$
      strtrim(string(stats.mean_rej,format='(F12.2)'),2)+$
      '\pm'+strtrim(string(stats.sigma_rej,format='(F12.2)'),2)+')'
;   xstr = strtrim(string(stats.median,format='(F12.2)'),2)+'\pm'+strtrim(string(stats.sig68mean,format='(F12.2)'),2)

    ytitle = 'log ([O III]/H\alpha)_{cor}'
;   ytitle = 'log [L([O III])_{cor}/\psi(H\alpha)]-'+string(hasfrconstoffset,format='(I0)')
;   yrange2 = sfroiiharange - hasfrconstoffset

    sdss_lineplot, x, y, xerr, yerr, plottype=3, postscript=postscript, $
      xtitle='', ytitle=ytitle, xrange=xrange, yrange=yrange, $
      legendtype=0, /noerase, xtickname=replicate(' ',10), $
      xstyle=3, /right, /top, position=pos[*,1], charsize=charsize_4
    
    print & running = im_medxbin([x,xnfgs],[y,ynfgs],medbin,minx=minx,minpts=minpts,/verbose) & print
;   oploterror, running.binctr, running.medy, replicate(running.binsz/2.0,running.nbins), $
;     running.sigy, ps=-3, thick=postthick, errthick=postthick, /nohat, line=0, errstyle=0

    legend, '(b)', /left, /top, box=0, charsize=charsize_3, charthick=postthick
    legend, textoidl(xstr), /right, /top, box=0, charsize=charsize_3, charthick=postthick

    w = where(x ge 10.25 and x le 10.75)
;   junk = im_stats(y[w],/verbose)
    
; ##########################################################
; Panel 3: [O III]_obs, SFR(Ha)
; ##########################################################

    y1 = sdssdust[indx].oiii_5007_lum[0]
    y1err = sdssdust[indx].oiii_5007_lum[1]

    y2 = sdssnodust[indx].h_alpha_lum[0]
    y2err = sdssnodust[indx].h_alpha_lum[1]

    y = y1 - y2 - loghaconst - hasfrconstoffset
    yerr = sqrt(y1err^2 + y2err^2)

    stats = im_stats(y,/verbose,/baremin,/no_head)
    xstr = strtrim(string(stats.median_rej,format='(F12.2)'),2)+' ('+$
      strtrim(string(stats.mean_rej,format='(F12.2)'),2)+$
      '\pm'+strtrim(string(stats.sigma_rej,format='(F12.2)'),2)+')'
;   xstr = strtrim(string(stats.median,format='(F12.2)'),2)+'\pm'+strtrim(string(stats.sig68mean,format='(F12.2)'),2)

;   ytitle = 'log ([O III]_{obs}/H\alpha_{cor})'
    ytitle = 'log [10^{-41} L([O III])_{obs}/\psi(H\alpha)]'
    yrange2 = sfroiiharange - hasfrconstoffset

    sdss_lineplot, x, y, xerr, yerr, plottype=3, postscript=postscript, $
      xtitle=xtitle, ytitle=ytitle, xrange=xrange, yrange=yrange2, legendtype=0, /noerase, $
      xstyle=3, /right, /top, position=pos[*,2], charsize=charsize_4
;   djs_oplot, !x.crange, alog10(1.0)*[1,1], line=0, thick=(postthick-4)>2
    
    print & running = im_medxbin([x,xnfgs],[y,ynfgs],medbin,minx=minx,minpts=minpts,/verbose) & print
;   oploterror, running.binctr, running.medy, replicate(running.binsz/2.0,running.nbins), $
;     running.sigy, ps=-3, thick=postthick, errthick=postthick, /nohat, line=0, errstyle=0

    legend, '(c)', /left, /top, box=0, charsize=charsize_3, charthick=postthick
    legend, textoidl(xstr), /right, /top, box=0, charsize=charsize_3, charthick=postthick

    im_openclose, postscript=postscript, /close

; ------------------------------------------------------------
; E(Hb-Ha) vs [O III]/Ha [Integrated + SDSS]
; ------------------------------------------------------------

    psname = 'ehbha_vs_oiiiha_2panel'
    im_openclose, pspath+psname, postscript=postscript, xsize=8.5, ysize=5.5, encapsulated=encapsulated

    pagemaker, nx=2, ny=1, height=3.5, width=3.5, xmargin=[1.2,0.3], $
      ymargin=[0.9,1.1], xspace=0, yspace=0, xpage=8.5, ypage=5.5, $
      position=pos, /normal

    medbin = 0.05
    minpts = 10
    minx = 0.025

; Integrated    
    
    indx = where((atlasdust.oiii_5007[0]/atlasdust.oiii_5007[1] gt snrcut) and $
      (atlasnodust.ehbha_err gt 0.0),nindx)

    x = atlasnodust[indx].ehbha
    xerr = atlasnodust[indx].ehbha_err
    xebv = x/(k_lambda(4861,/odonnell)-k_lambda(6563,/odonnell))

    indxnfgs = where((nfgsdust.oiii_5007[0]/nfgsdust.oiii_5007[1] gt snrcut) and $
      (nfgsnodust.ehbha_err gt 0.0),nindxnfgs)

    xnfgs = nfgsnodust[indxnfgs].ehbha
    xerrnfgs = nfgsnodust[indxnfgs].ehbha_err

    xtitle = 'E(H\beta-H\alpha) [mag]'
    ytitle = 'log ([O III] \lambda5007/H\alpha)_{obs}'

    xrange = ehbharange
    yrange = oiiiharange
    
    y1 = atlasdust[indx].oiii_5007[0]
    y1err = atlasdust[indx].oiii_5007[1]

    y2 = atlasdust[indx].h_alpha[0]
    y2err = atlasdust[indx].h_alpha[1]

    y = alog10(y1/y2)
    yerr = im_compute_error(y1,y1err,y2,y2err,/log)

    y1nfgs = nfgsdust[indxnfgs].oiii_5007[0]
    y1errnfgs = nfgsdust[indxnfgs].oiii_5007[1]

    y2nfgs = nfgsdust[indxnfgs].h_alpha[0]
    y2errnfgs = nfgsdust[indxnfgs].h_alpha[1]

    ynfgs = alog10(y1nfgs/y2nfgs)
    yerrnfgs = im_compute_error(y1nfgs,y1errnfgs,y2nfgs,y2errnfgs,/log)

    stats = im_stats([y,ynfgs],/verbose,/baremin)
    rcor = r_correlate([x,xnfgs],[y,ynfgs],zd=zd,probd=probd)
;   splog, 'Spearman rank test: ', rcor, probd, zd

    atlas1d_lineplot, x, y, xerr, yerr, plottype=3, postscript=postscript, ysty=1, $
      ytitle=ytitle, xtitle=xtitle, xrange=xrange, yrange=yrange, legendtype=0, $
      xstyle=11, /right, /top, xnfgs=xnfgs, ynfgs=ynfgs, xerrnfgs=xerrnfgs, $
      yerrnfgs=yerrnfgs, position=pos[*,0], charsize=charsize_6
    axis, /xaxis, xrange=interpol(xebv,x,!x.crange), xthick=postthick, xstyle=1, $
      charsize=charsize_6, charthick=postthick
    im_xyouts_title, xtitle='E(B-V) [mag]', charsize=charsize_6, charthick=postthick

;   legend, '(a)', /left, /top, box=0, charsize=charsize_2, charthick=postthick

    print & running = im_medxbin([x,xnfgs],[y,ynfgs],medbin,minx=minx,minpts=minpts,/verbose) & print

; overplot reddening relations

    ratio = 0.0

    x_ebv = findgen((1.3-(-0.1))/0.1)*0.1+(-0.1)
    yodonnell =  -0.4*x_ebv*(k_lambda(5007.0,/odonnell)-k_lambda(6563.0,/odonnell)) + ratio
    ycharlot = -0.4*x_ebv*(k_lambda(5007,/charlot)-k_lambda(6563,/charlot)) + ratio

    djs_oplot, x_ebv, yodonnell, line=0, thick=postthick
    djs_oplot, x_ebv, ycharlot, line=2, thick=postthick

    label = ["O'Donnell (1994)",'Charlot & Fall (2000)']
    linestyle = [0,2]
    legend, label, /right, /top, box=0, linestyle=linestyle, $
      charsize=charsize_3, charthick=postthick, thick=postthick

; SDSS    
    
    indx = where((sdssdust.oiii_5007[0]/sdssdust.oiii_5007[1] gt snrcut) and $
      (sdssnodust.ehbha_err gt 0.0),nindx)

    x = sdssnodust[indx].ehbha
    xerr = sdssnodust[indx].ehbha_err
    xebv = x/(k_lambda(4861,/odonnell)-k_lambda(6563,/odonnell))

    y1 = sdssdust[indx].oiii_5007[0]
    y1err = sdssdust[indx].oiii_5007[1]

    y2 = sdssdust[indx].h_alpha[0]
    y2err = sdssdust[indx].h_alpha[1]

    y = alog10(y1/y2)
    yerr = im_compute_error(y1,y1err,y2,y2err,/log)

    stats = im_stats(y,/verbose,/baremin,/no_head)
    rcor = r_correlate(x,y,zd=zd,probd=probd)
;   splog, 'Spearman rank test: ', rcor, probd, zd

    sdss_lineplot, x, y, xerr, yerr, plottype=3, postscript=postscript, $
      ytitle='', xtitle=xtitle, xrange=xrange, yrange=yrange, legendtype=0, $
      xstyle=11, /right, /top, position=pos[*,1], charsize=charsize_6, /noerase, $
      ytickname=replicate(' ',10)
    axis, /xaxis, xrange=interpol(xebv,x,!x.crange), xthick=postthick, xstyle=1, $
      charsize=charsize_6, charthick=postthick
    im_xyouts_title, xtitle='E(B-V) [mag]', charsize=charsize_6, charthick=postthick

;   legend, '(b)', /left, /top, box=0, charsize=charsize_2, charthick=postthick

    print & running = im_medxbin(x,y,medbin,minx=minx,minpts=minpts,/verbose) & print

; overplot reddening relations

    ratio = 0.0

    x_ebv = findgen((1.3-(-0.1))/0.1)*0.1+(-0.1)
    yodonnell =  -0.4*x_ebv*(k_lambda(5007.0,/odonnell)-k_lambda(6563.0,/odonnell)) + ratio
    ycharlot = -0.4*x_ebv*(k_lambda(5007,/charlot)-k_lambda(6563,/charlot)) + ratio

    djs_oplot, x_ebv, yodonnell, line=0, thick=postthick
    djs_oplot, x_ebv, ycharlot, line=2, thick=postthick

    label = ["O'Donnell (1994)",'Charlot & Fall (2000)']
    linestyle = [0,2]
;   legend, label, /right, /top, box=0, linestyle=linestyle, $
;     charsize=charsize_2, charthick=postthick, thick=postthick

    im_openclose, postscript=postscript, /close

; ------------------------------------------------------------
; 12+log (O/H) [Te] vs [O II]/Ha
; ------------------------------------------------------------

    psname = '12oh_te_vs_oiiha'
    im_openclose, pspath+psname, postscript=postscript, xsize=8.5, ysize=8.0, encapsulated=encapsulated

    pagemaker, nx=1, ny=1, xspace=0, yspace=0.0, width=6.6, height=6.6, $
      xmargin=[1.5,0.4], ymargin=[0.3,1.1], xpage=8.5, ypage=8.0, $
      position=pos, /normal

    if (n_elements(hiiall) eq 0L) then hiiall = read_hii_regions(/nosdss,/nokiss)
    
    indx = where((hiiall.oii_h_alpha gt -900.0) and (hiiall.ZT_log12oh gt -900.0),nindx)
    
    x = hiiall[indx].ZT_log12oh
    xerr = hiiall[indx].ZT_log12oh_err

    y = hiiall[indx].oii_h_alpha
    yerr = hiiall[indx].oii_h_alpha_err

;   xrange = interpol(xOHgood,xgood,reverse(oiiiniirange))
    xrange = hiiohrange
    yrange = oiiharange

;   xtitle = '12 + log (O/H) [T_{e}]'
    xtitle = '12 + log (O/H)_{Te}'
    ytitle = 'log ([O II]/H\alpha)_{cor}'

    zindicators_lineplot, x, y, xerr, yerr, plottype=2, postscript=postscript, $
      xtitle=xtitle, ytitle=ytitle, xrange=xrange, yrange=yrange, legendtype=0, $
      position=pos[*,0], charsize=charsize_8, /right, /top, hiipsize=0.7, $
      hiicolor='purple', /errorleft

; overplot the Kewley grids

    if keyword_set(postscript) then thisthick = 6.0 else thickthick = 2.0
    
    plot_kewley_grids, plotnumber=17, model=3, labeltype=4, /log12oh, /overplot, $
      /noZgrid, thick=3.0, Ugridcolor='dark green', Ulinestyle=0, $
      postscript=postscript, Ulabelcolor='', Ucharsize=1.0, charthick=thisthick
    xyouts, 9.08, -0.18, 'log U', align=0.5, /data, charsize=1.1, charthick=thisthick

    r12 = 8.15
    r23 = 8.7
    
    oplot, r12*[1,1], [-1.25,0.75], line=1, thick=postthick
    oplot, r23*[1,1], [-1.25,0.75], line=1, thick=postthick

    xyouts, 7.7, 0.65, 'R1', align=0.5, /data, charsize=charsize_8, charthick=postthick
    xyouts, 8.4, 0.65, 'R2', align=0.5, /data, charsize=charsize_8, charthick=postthick
    xyouts, 9.0, 0.65, 'R3', align=0.5, /data, charsize=charsize_8, charthick=postthick

    im_openclose, postscript=postscript, /close    

; ------------------------------------------------------------
; ([O III]/Hb)/([N II]/Ha) vs E(Hb-Ha) [Integrated + SDSS]
; ------------------------------------------------------------

    psname = '12oh_oiiinii_niiha_vs_ehbha_2panel'
    im_openclose, pspath+psname, postscript=postscript, xsize=8.5, ysize=4.65, encapsulated=encapsulated

    pagemaker, nx=2, ny=1, height=3.15, width=3.15, xmargin=[1.1,1.1], $
      ymargin=[0.5,1.0], xspace=0, yspace=0, xpage=8.5, ypage=4.65, $
      position=pos, /normal

    xtitle = '12 + log (O/H)'
    ytitle = 'E(H\beta-H\alpha) [mag]'

    xrange = ohrange2
    yrange = ehbharange

; Integrated    
    
    indx = where((atlasdust.zstrong_12oh_oiiinii_niiha gt -900) and (atlasnodust.ehbha_err gt 0),nindx)

    x = atlasdust[indx].zstrong_12oh_oiiinii_niiha
    xerr = sqrt(atlasdust[indx].zstrong_12oh_oiiinii_niiha_err^2 + 0.1^2)

    y = atlasnodust[indx].ehbha
    yerr = atlasnodust[indx].ehbha_err
    yebv = y/(k_lambda(4861,/odonnell)-k_lambda(6563,/odonnell))

    indxnfgs = where((nfgsdust.zstrong_12oh_oiiinii_niiha gt -900) and (nfgsnodust.ehbha_err gt 0.0),nindxnfgs)

    xnfgs = nfgsdust[indxnfgs].zstrong_12oh_oiiinii_niiha
    xerrnfgs = sqrt(nfgsdust[indxnfgs].zstrong_12oh_oiiinii_niiha_err^2 + 0.1^2)

    ynfgs = nfgsnodust[indxnfgs].ehbha
    yerrnfgs = nfgsnodust[indxnfgs].ehbha_err

    xbig = [x,xnfgs]
    ybig = [y,ynfgs]

    atlas1d_lineplot, x, y, xerr, yerr, plottype=3, postscript=postscript, $n
      xtitle=xtitle, ytitle=ytitle, xrange=xrange, yrange=yrange, legendtype=0, $
      /left, /top, position=pos[*,0], charsize=charsize_4, $
      xnfgs=xnfgs, ynfgs=ynfgs, xerrnfgs=xerrnfgs, yerrnfgs=yerrnfgs

; SDSS    

    indx = where((sdssdust.zstrong_12oh_oiiinii_niiha gt -900) and (sdssnodust.ehbha_err gt 0.0),nindx)

    x = sdssdust[indx].zstrong_12oh_oiiinii_niiha
    xerr = sqrt(sdssdust[indx].zstrong_12oh_oiiinii_niiha_err^2 + 0.1^2)

    y = sdssnodust[indx].ehbha
    yerr = sdssnodust[indx].ehbha_err
    yebv = y/(k_lambda(4861,/odonnell)-k_lambda(6563,/odonnell))

    sdss_lineplot, x, y, xerr, yerr, plottype=3, postscript=postscript, $
      xtitle=xtitle, ytitle='', xrange=xrange, yrange=yrange, legendtype=0, $
      ystyle=11, xmargin=[8,6], /left, /top, position=pos[*,1], $
      ytickname=replicate(' ',10), /noerase, charsize=charsize_4
    axis, /yaxis, yrange=interpol(yebv,y,!y.crange), ythick=postthick, $
      charthick=postthick, charsize=charsize_4, ystyle=1
    im_xyouts_title, ytitle='E(B-V) [mag]', charsize=charsize_4, charthick=postthick

    im_openclose, postscript=postscript, /close    

; ------------------------------------------------------------
; R23 vs [O II]/Ha_cor [Integrated + SDSS]
; ------------------------------------------------------------

    psname = 'r23_vs_oiiha_2panel'
    im_openclose, pspath+psname, postscript=postscript, xsize=8.5, ysize=5.0, encapsulated=encapsulated

    pagemaker, nx=2, ny=1, height=3.5, width=3.5, xmargin=[1.2,0.3], $
      ymargin=[0.4,1.1], xspace=0, yspace=0, xpage=8.5, ypage=5.0, $
      position=pos, /normal

    xtitle = 'log (R_{23})'
    ytitle = 'log ([O II]/H\alpha)_{cor}'

    xrange = R23range
    yrange = oiiharange

; Integrated    

    lineratio, atlasnodust, ['OII_3727','OIII_4959','OIII_5007'], 'H_BETA', $
      'OII_3727', 'H_ALPHA', x, xerr, y, yerr, index=indx, nindex=nindx, snrcut=snrcut, $
      xsyserr=syserr, ysyserr=syserr
    
    lineratio, nfgsnodust, ['OII_3727','OIII_4959','OIII_5007'], 'H_BETA', $
      'OII_3727', 'H_ALPHA', xnfgs, xerrnfgs, ynfgs, yerrnfgs, $
      index=indxnfgs, snrcut=snrcut, xsyserr=nfgssyserr, ysyserr=nfgssyserr

    xbig = [x,xnfgs]
    ybig = [y,ynfgs]

; HII regions

    good = where((hii.oii_h_alpha gt -900.0) and (hii.R23 gt -900.0))
    xregion = hii[good].R23 & xerrregion = hii[good].R23_err
    yregion = hii[good].oii_h_alpha & yerrregion = hii[good].oii_h_alpha_err

    atlas1d_lineplot, x, y, xerr, yerr, plottype=3, postscript=postscript, $
      xtitle=xtitle, ytitle=ytitle, xrange=xrange, yrange=yrange, legendtype=0, $
      /left, /top, position=pos[*,0], charsize=charsize_6, $
      xregion=xregion, yregion=yregion, xerrregion=xerrregion, yerrregion=yerrregion, $
      xnfgs=xnfgs, ynfgs=ynfgs, xerrnfgs=xerrnfgs, yerrnfgs=yerrnfgs

; show that this is basically y = x + constant.  take the intercept to
; be log 2.86 and assume a mean [O III]/[O II] ratio of OIII_OII_RATIO

    oiii_oii_ratio = 0.3
    
    ymean = dkboxstats(ybig,xwidth=10,boxstat='mean')
    ysig = dkboxstats(ybig,xwidth=10,boxstat='sigma',sigrej=2.0)
;   djs_oplot, xbig, ymean+ysig, line=2
;   djs_oplot, xbig, ymean-ysig, line=2
    
    xaxis = findgen(((+2.0)-(-1.0))/0.1)*0.1+(-1.0)
    yline = poly(xaxis,[-alog10(2.86)-alog10(1+oiii_oii_ratio),1.0])
    oplot, xaxis, yline, linestyle=2, thick=postthick

    fitstr = textoidl('log [O II]/H\alpha = log R_{23} - log 2.86 - log (1+'+$
      string(oiii_oii_ratio,format='(F4.2)')+')')
;   legend, fitstr, /left, /bottom, box=0, charsize=charsize, charthick=postthick, $
;     clear=postscript;, line=2

; SDSS

    lineratio, sdssnodust, ['OII_3727','OIII_4959','OIII_5007'], 'H_BETA', $
      'OII_3727', 'H_ALPHA', x, xerr, y, yerr, index=indx, nindex=nindx, snrcut=snrcut
    
    sdss_lineplot, x, y, xerr, yerr, plottype=3, postscript=postscript, $
      xtitle=xtitle, ytitle='', xrange=xrange, yrange=yrange, legendtype=0, $
      /left, /top, position=pos[*,1], $
      ytickname=replicate(' ',10), /noerase, charsize=charsize_6
    oplot, xaxis, yline, linestyle=2, thick=postthick

    im_openclose, postscript=postscript, /close    

; ------------------------------------------------------------
; 3-panel - L(B) vs Hb/Ha - Integrated
; ------------------------------------------------------------
    
    psname = 'mass_vs_hbha_3panel'
    im_openclose, pspath+psname, postscript=postscript, xsize=8.5, ysize=9.5, encapsulated=encapsulated

    pagemaker, nx=1, ny=3, height=[2.5,2.5,2.5], width=7.0, xmargin=[1.1,0.4], $
      ymargin=[0.9,1.1], xspace=0, yspace=0, xpage=8.5, ypage=9.5, $
      position=pos, /normal

    medbin = 0.5
    minpts = 10
    minx = 7.75
    
    indx = where((atlasdust.mass_bv_b gt -900) and (atlasnodust.ebv_hahb_err gt 0.0) and $
      (atlasdust.h_beta_ew_uncor[0] ge ewcut),nindx)

    ewstatcor = 4.0*atlasdust[indx].babs_h_beta_continuum[0]
    
    x = atlasdust[indx].mass_bv_b
    xerr = atlasdust[indx].mass_bv_b_err
    xabs = atlasdust[indx].m_b_obs
    
    indxnfgs = where((nfgsnodust.mass_bv_b gt -900.0) and (nfgsnodust.ebv_hahb_err gt 0.0) and $
      (nfgsdust.h_beta_ew_uncor[0] ge ewcut),nindxnfgs)

    ewstatcornfgs = 4.0*nfgsdust[indxnfgs].babs_h_beta_continuum[0]

    xnfgs = nfgsdust[indxnfgs].mass_bv_b
    xerrnfgs = nfgsdust[indxnfgs].mass_bv_b_err
    
    xrange = massrange
    yrange = hbharange; + alog10(HaHb)

    xtitle = 'log M [M'+sunsymbol()+']'
    ytitle = 'log (H\beta/H\alpha)'

; ##########################################################
; Panel 1: No absorption correction, A(Hb)=0, A(Ha)=0
; ##########################################################

    hb = atlasdust[indx].h_beta_uncor[0]
    hb_err = atlasdust[indx].h_beta_uncor[1]
    
    ha = atlasdust[indx].h_alpha[0]
    ha_err = atlasdust[indx].h_alpha[1]
    
    y = alog10(hb/ha); + alog10(HaHb)
    yerr = im_compute_error(hb,hb_err,ha,ha_err,/log)

    hbnfgs = nfgsdust[indxnfgs].h_beta_uncor[0]
    hbnfgs_err = nfgsdust[indxnfgs].h_beta_uncor[1]
    
    hanfgs = nfgsdust[indxnfgs].h_alpha[0]
    hanfgs_err = nfgsdust[indxnfgs].h_alpha[1]
    
    ynfgs = alog10(hbnfgs/hanfgs); + alog10(HaHb)
    yerrnfgs = im_compute_error(hbnfgs,hbnfgs_err,hanfgs,hanfgs_err,/log)

    stats = im_stats([y,ynfgs],/verbose,/baremin)
    xstr = strtrim(string(stats.median_rej,format='(F12.2)'),2)+' ('+$
      strtrim(string(stats.mean_rej,format='(F12.2)'),2)+$
      '\pm'+strtrim(string(stats.sigma_rej,format='(F12.2)'),2)+')'
;   xstr = strtrim(string(stats.median,format='(F12.2)'),2)+'\pm'+strtrim(string(stats.sig68mean,format='(F12.2)'),2)
    
    ytitle = 'log (H\beta/H\alpha)_{obs}'
;   ytitle = 'log (H\beta/H\alpha)_{obs}+'+string(alog10(HaHb),format='(F4.2)')

    atlas1d_lineplot, x, y, xerr, yerr, plottype=3, postscript=postscript, $
      xtitle='', ytitle=ytitle, xrange=xrange, yrange=yrange, legendtype=0, $
      charsize=charsize_5, xstyle=11, xtickname=replicate(' ',10), position=pos[*,0], $
      xnfgs=xnfgs, ynfgs=ynfgs, xerrnfgs=xerrnfgs, yerrnfgs=yerrnfgs
    axis, /xaxis, xrange=interpol(xabs,x,!x.crange), xthick=postthick, xstyle=1, $
      charsize=charsize_5, charthick=postthick
    im_xyouts_title, xtitle='M_{B} [mag]', charsize=charsize_5, charthick=postthick
;   djs_oplot, !x.crange, alog10(1.0)*[1,1], line=0, thick=postthick
    djs_oplot, !x.crange, alog10(1.0/HaHb)*[1,1], line=2, thick=postthick

;   print & running = im_medxbin([x,xnfgs],[y,ynfgs],medbin,minx=minx,minpts=minpts,/verbose) & print
;   oplot, running.binctr, running.medy, thick=postthick, line=0
;   oplot, running.binctr, running.medy+running.sigy, thick=postthick, line=2
;   oplot, running.binctr, running.medy-running.sigy, thick=postthick, line=2
;   oploterror, running.binctr, running.medy, replicate(running.binsz/2.0,running.nbins), $
;     running.sigy, ps=-3, thick=postthick, errthick=postthick, /nohat, line=0, errstyle=0

;   legend, textoidl('(a) EW(H\beta) (uncorrected) > '+string(ewcut,format='(I0)'))+' '+angstrom(), $
;     /left, /top, box=0, charsize=charsize_3, charthick=postthick
    legend, '(a)', /left, /top, box=0, charsize=charsize_3, charthick=postthick
    legend, textoidl(xstr), /right, /top, box=0, charsize=charsize_3, charthick=postthick

    label = ['No Stellar Absorption Correction']
    legend, textoidl(label), /left, /bottom, box=0, charsize=charsize_3, $
      charthick=postthick

; ############################################################
; Panel 2: Absorption-Corrected, A(Hb)=0, A(Ha)=0
; ############################################################

    hb = atlasdust[indx].h_beta[0]
    hb_err = atlasdust[indx].h_beta[1]
    
    ha = atlasdust[indx].h_alpha[0]
    ha_err = atlasdust[indx].h_alpha[1]
    
    y = alog10(hb/ha); + alog10(HaHb)
    yerr = im_compute_error(hb,hb_err,ha,ha_err,/log)

    hbnfgs = nfgsdust[indxnfgs].h_beta[0]
    hbnfgs_err = nfgsdust[indxnfgs].h_beta[1]
    
    hanfgs = nfgsdust[indxnfgs].h_alpha[0]
    hanfgs_err = nfgsdust[indxnfgs].h_alpha[1]
    
    ynfgs = alog10(hbnfgs/hanfgs); + alog10(HaHb)
    yerrnfgs = im_compute_error(hbnfgs,hbnfgs_err,hanfgs,hanfgs_err,/log)

    stats = im_stats([y,ynfgs],/verbose,/baremin,/no_head)
    xstr = strtrim(string(stats.median_rej,format='(F12.2)'),2)+' ('+$
      strtrim(string(stats.mean_rej,format='(F12.2)'),2)+$
      '\pm'+strtrim(string(stats.sigma_rej,format='(F12.2)'),2)+')'
;   xstr = strtrim(string(stats.median,format='(F12.2)'),2)+'\pm'+strtrim(string(stats.sig68mean,format='(F12.2)'),2)

;   ytitle = 'log (H\beta/H\alpha)_{obs}+'+string(alog10(HaHb),format='(F4.2)')
    ytitle = 'log (H\beta/H\alpha)_{obs}'

    atlas1d_lineplot, x, y, xerr, yerr, plottype=3, postscript=postscript, $
      xtitle='', ytitle=ytitle, xrange=xrange, yrange=yrange, legendtype=0, $
      charsize=charsize_5, xtickname=replicate(' ',10), $
      position=pos[*,1], /noerase, xnfgs=xnfgs, ynfgs=ynfgs, xerrnfgs=xerrnfgs, yerrnfgs=yerrnfgs
;   djs_oplot, !x.crange, alog10(1.0)*[1,1], line=0, thick=postthick
    djs_oplot, !x.crange, alog10(1.0/HaHb)*[1,1], line=2, thick=postthick

;   print & running = im_medxbin([x,xnfgs],[y,ynfgs],medbin,minx=minx,minpts=minpts,/verbose) & print
;   oplot, running.binctr, running.medy, thick=postthick, line=0
;   oplot, running.binctr, running.medy+running.sigy, thick=postthick, line=2
;   oplot, running.binctr, running.medy-running.sigy, thick=postthick, line=2
;   oploterror, running.binctr, running.medy, replicate(running.binsz/2.0,running.nbins), $
;     running.sigy, ps=-3, thick=postthick, errthick=postthick, /nohat, line=0, errstyle=0

    legend, '(b)', /left, /top, box=0, charsize=charsize_3, charthick=postthick
    legend, textoidl(xstr), /right, /top, box=0, charsize=charsize_3, charthick=postthick

    label = ['Stellar Absorption Corrected']
    legend, textoidl(label), /left, /bottom, box=0, charsize=charsize_3, $
      charthick=postthick

; ############################################################
; Panel 3: Absorption-Corrected, A(Hb)=0, Individual A(Ha) 
; ############################################################

;   hb = alog10((atlasdust[indx].h_beta_uncor[0]+ewstatcor)*4*!pi*((atlasdust[indx].distance*3.086D24)^2)/lsun)
;   hb = alog10(atlasdust[indx].h_beta_uncor[0]*4*!pi*((atlasdust[indx].distance*3.086D24)^2)/lsun)
    hb = atlasdust[indx].h_beta_lum[0]
    hb_err = atlasdust[indx].h_beta_lum[1]

    ha = atlasnodust[indx].h_alpha_lum[0]
    ha_err = atlasnodust[indx].h_alpha_lum[1]

    y = hb - ha - loghaconst - hasfrconstoffset
    yerr = sqrt(ha_err^2 + hb_err^2)

;   hbnfgs = alog10((nfgsdust[indxnfgs].h_beta_uncor[0]+ewstatcornfgs)*4*!pi*((nfgsdust[indxnfgs].distance*3.086D24)^2)/lsun)
;   hbnfgs = alog10(nfgsdust[indxnfgs].h_beta_uncor[0]*4*!pi*((nfgsdust[indxnfgs].distance*3.086D24)^2)/lsun)
    hbnfgs = nfgsdust[indxnfgs].h_beta_lum[0]
    hbnfgs_err = nfgsdust[indxnfgs].h_beta_lum[1]

    hanfgs = nfgsnodust[indxnfgs].h_alpha_lum[0]
    hanfgs_err = nfgsnodust[indxnfgs].h_alpha_lum[1]

    ynfgs = hbnfgs - hanfgs - loghaconst - hasfrconstoffset
    yerrnfgs = sqrt(hanfgs_err^2 + hbnfgs_err^2)
       
    stats = im_stats([y,ynfgs],/verbose,/baremin,/no_head)
    xstr = strtrim(string(stats.median_rej,format='(F12.2)'),2)+' ('+$
      strtrim(string(stats.mean_rej,format='(F12.2)'),2)+$
      '\pm'+strtrim(string(stats.sigma_rej,format='(F12.2)'),2)+')'
;   xstr = strtrim(string(stats.median,format='(F12.2)'),2)+'\pm'+strtrim(string(stats.sig68mean,format='(F12.2)'),2)
        
;   ytitle = 'log [L(H\beta)_{obs}/\psi]-'+string(hasfrconstoffset,format='(I0)')
    ytitle = 'log [10^{-41} L(H\beta)_{obs}/\psi(H\alpha)]'
    yrange = sfrhbharange - hasfrconstoffset

    atlas1d_lineplot, x, y, xerr, yerr, plottype=3, postscript=postscript, $
      xtitle=xtitle, ytitle=ytitle, xrange=xrange, yrange=yrange, legendtype=0, $
      charsize=charsize_5, position=pos[*,2], /noerase, $
      xnfgs=xnfgs, ynfgs=ynfgs, xerrnfgs=xerrnfgs, yerrnfgs=yerrnfgs
;   djs_oplot, !x.crange, -(alog10(HaHb)+loghaconst+hasfrconstoffset)*[1,1], line=0, thick=(postthick-4)>2

    print & running = im_medxbin([x,xnfgs],[y,ynfgs],medbin,minx=minx,minpts=minpts,/verbose) & print
;   oplot, running.binctr, running.medy, thick=postthick, line=0
;   oplot, running.binctr, running.medy+running.sigy, thick=postthick, line=2
;   oplot, running.binctr, running.medy-running.sigy, thick=postthick, line=2
;   oploterror, running.binctr, running.medy, replicate(running.binsz/2.0,running.nbins), $
;     running.sigy, ps=-3, thick=postthick, errthick=postthick, /nohat, line=0, errstyle=0

    legend, '(c)', /left, /top, box=0, charsize=charsize_3, charthick=postthick
    legend, textoidl(xstr), /right, /top, box=0, charsize=charsize_3, charthick=postthick

    label = ['Stellar Absorption Corrected']
    legend, textoidl(label), /left, /bottom, box=0, charsize=charsize_3, $
      charthick=postthick

    im_openclose, postscript=postscript, /close

; ------------------------------------------------------------
; L(B) versus SFR(Ha)/L(U)_obs - Integrated
; ------------------------------------------------------------
    
    psname = 'lb_vs_sfr_ha_lu_obs'
    im_openclose, pspath+psname, postscript=postscript, xsize=8.5, ysize=8.0, encapsulated=encapsulated

    pagemaker, nx=1, ny=1, xspace=0, yspace=0.0, width=6.8, height=6.8, $
      xmargin=[1.4,0.3], ymargin=[0.9,1.1], xpage=8.5, ypage=8.8, $
      position=pos, /normal

    if keyword_set(blackwhite) then begin
       notlirgscolor = 'gray'
       lirgscolor = 'dark gray'
       nolirdatacolor = 'dark gray'
       medcolor = 'black'
    endif else begin
       notlirgscolor = 'gray'
       lirgscolor = 'medium orchid'
       nolirdatacolor = 'dodger blue'
       medcolor = 'black'
    endelse

    medbin = 1.0
    minpts = 5
    minx = 8.0
    
; Atlas - SF
    
    indx = where((atlasdust.b_lum_obs gt -900) and (atlasnodust.sfr_h_alpha gt -900.0) and $
      (atlasdust.u_lum_obs gt -900.0),nindx)

    x = atlasdust[indx].b_lum_obs; - alog10(LBnorm)
    xerr = atlasdust[indx].b_lum_obs_err
    xabs = atlasdust[indx].m_b_obs

    uband = atlasdust[indx].u_lum_obs + alog10(lusun) - alog10(Ulumnorm)
    uband_err = atlasdust[indx].u_lum_obs_err

    sfr_ha = atlasnodust[indx].sfr_h_alpha
    sfr_ha_err = atlasnodust[indx].sfr_h_alpha_err

    y = sfr_ha - uband
    yerr = sqrt(sfr_ha_err^2 + uband_err^2)

; Atlas - AGN
    
    indx_agn = where((atlasdust_agn.b_lum_obs gt -900) and (atlasnodust_agn.sfr_h_alpha gt -900.0) and $
      (atlasdust_agn.u_lum_obs gt -900.0),nindx_agn)

    x_agn = atlasdust_agn[indx_agn].b_lum_obs; - alog10(LBnorm)
    xerr_agn = atlasdust_agn[indx_agn].b_lum_obs_err

    uband_agn = atlasdust_agn[indx_agn].u_lum_obs + alog10(lusun) - alog10(Ulumnorm)
    uband_err_agn = atlasdust_agn[indx_agn].u_lum_obs_err

    sfr_ha_agn = atlasnodust_agn[indx_agn].sfr_h_alpha
    sfr_ha_err_agn = atlasnodust_agn[indx_agn].sfr_h_alpha_err

    y_agn = sfr_ha_agn - uband_agn
    yerr_agn = sqrt(sfr_ha_err_agn^2 + uband_err_agn^2)
    
; NFGS - SF
    
    indxnfgs = where((nfgsdust.b_lum_obs gt -900.0) and (nfgsnodust.sfr_h_alpha gt -900.0) and $
      (nfgsdust.u_lum_obs gt -900.0),nindxnfgs)

    xnfgs = nfgsdust[indxnfgs].b_lum_obs; - alog10(LBnorm)
    xerrnfgs = nfgsdust[indxnfgs].b_lum_obs_err
    
    ubandnfgs = nfgsdust[indxnfgs].u_lum_obs + alog10(lusun) - alog10(Ulumnorm)
    ubandnfgs_err = nfgsdust[indxnfgs].u_lum_obs_err

    sfrnfgs_ha = nfgsnodust[indxnfgs].sfr_h_alpha
    sfrnfgs_ha_err = nfgsnodust[indxnfgs].sfr_h_alpha_err

    ynfgs = sfrnfgs_ha - ubandnfgs
    yerrnfgs = sqrt(sfrnfgs_ha_err^2 + ubandnfgs_err^2)

; NFGS - AGN
    
    indxnfgs_agn = where((nfgsdust_agn.b_lum_obs gt -900.0) and (nfgsnodust_agn.sfr_h_alpha gt -900.0) and $
      (nfgsdust_agn.u_lum_obs gt -900.0),nindxnfgs_agn)

    xnfgs_agn = nfgsdust_agn[indxnfgs_agn].b_lum_obs; - alog10(LBnorm)
    xerrnfgs_agn = nfgsdust_agn[indxnfgs_agn].b_lum_obs_err
    
    ubandnfgs_agn = nfgsdust_agn[indxnfgs_agn].u_lum_obs + alog10(lusun) - alog10(Ulumnorm)
    ubandnfgs_err_agn = nfgsdust_agn[indxnfgs_agn].u_lum_obs_err

    sfrnfgs_ha_agn = nfgsnodust_agn[indxnfgs_agn].sfr_h_alpha
    sfrnfgs_ha_err_agn = nfgsnodust_agn[indxnfgs_agn].sfr_h_alpha_err

    ynfgs_agn = sfrnfgs_ha_agn - ubandnfgs_agn
    yerrnfgs_agn = sqrt(sfrnfgs_ha_err_agn^2 + ubandnfgs_err_agn^2)

; combine the two samples - SF

    xbig = [x,xnfgs] & xerrbig = [xerr,xerrnfgs]
    ybig = [y,ynfgs] & yerrbig = [yerr,yerrnfgs]
    lir = [atlasdust[indx].ir_lum,nfgsdust[indxnfgs].ir_lum]

    lirgs = where((lir gt -900.0) and (lir gt 11.0))
    notlirgs = where((lir gt -900.0) and (lir lt 11.0))
    nolirdata = where((lir lt -900.0))
    
; combine the two samples - AGN

    xbig_agn = [x_agn,xnfgs_agn] & xerrbig_agn = [xerr_agn,xerrnfgs_agn]
    ybig_agn = [y_agn,ynfgs_agn] & yerrbig_agn = [yerr_agn,yerrnfgs_agn]
    lir_agn = [atlasdust_agn[indx_agn].ir_lum,nfgsdust_agn[indxnfgs_agn].ir_lum]

    lirgs_agn = where((lir_agn gt -900.0) and (lir_agn gt 11.0))
    notlirgs_agn = where((lir_agn gt -900.0) and (lir_agn lt 11.0))
    nolirdata_agn = where((lir_agn lt -900.0))
    
; make the plot

    xtitle = 'log L(B) [L(B)'+sunsymbol()+']'
    ytitle = 'log [10^{42} \psi(H\alpha)/L(U)_{obs}] [erg s^{-1}/'+sfr_units()+']'

    xrange = LBrange; - alog10(LBnorm)
    yrange = alog10([1D-45,1D-41]) - alog10(sfrnorm) + alog10(Ulumnorm)
   
; SF galaxies    
    
    atlas1d_lineplot, xbig[notlirgs], ybig[notlirgs], xerrbig[notlirgs], yerrbig[notlirgs], $
      plottype=1, postscript=postscript, xtitle=xtitle, ytitle=ytitle, xrange=xrange, $
      xstyle=11, yrange=yrange, legendtype=0, /right, /top, yfrac=1.5, xfrac=1.5, charsize=charsize_8, $
      atlaspsize=1.3, position=pos[*,0], atlascolor=notlirgscolor
    axis, /xaxis, xrange=interpol(xabs,x,!x.crange), xthick=postthick, xstyle=1, $
      charsize=charsize_8, charthick=postthick
    im_xyouts_title, xtitle='M_{B} [mag]', charsize=charsize_8, charthick=postthick

    im_symbols, 108, psize=1.4, /fill, color=fsc_color(lirgscolor,!d.table_size-5)
    djs_oplot, xbig[lirgs], ybig[lirgs], ps=8

    im_symbols, 105, psize=1.4, /fill, color=fsc_color(nolirdatacolor,!d.table_size-6)
    djs_oplot, xbig[nolirdata], ybig[nolirdata], ps=8

; AGN galaxies    
    
    atlas1d_lineplot, xbig_agn[notlirgs_agn], ybig_agn[notlirgs_agn], xerrbig_agn[notlirgs_agn], yerrbig_agn[notlirgs_agn], $
      plottype=1, postscript=postscript, xtitle=xtitle, ytitle=ytitle, xrange=xrange, $
      yrange=yrange, legendtype=0, /right, /top, yfrac=1.5, xfrac=1.5, charsize=charsize_8, $
      atlaspsize=1.1, position=pos[*,0], atlasfill=0, /overplot, thick=postthick, atlascolor=notlirgscolor

    im_symbols, 108, psize=1.3, fill=0, color=fsc_color(lirgscolor,!d.table_size-5), thick=postthick
    djs_oplot, xbig_agn[lirgs_agn], ybig_agn[lirgs_agn], ps=8

    im_symbols, 105, psize=1.2, fill=0, color=fsc_color(nolirdatacolor,!d.table_size-6), thick=postthick
    djs_oplot, xbig_agn[nolirdata_agn], ybig_agn[nolirdata_agn], ps=8

; fit only to certain data

;   xfit = [xbig,xbig_agn]
;   yfit = [ybig,ybig_agn]
    xfit = xbig
    yfit = ybig
;   xfit = [xbig[notlirgs],xbig[nolirdata]]
;   yfit = [ybig[notlirgs],ybig[nolirdata]]
    
; 9.0 is M_B = -17.07 and 8.7 is M_B = -16.3; 8.57 is -16

    LBlocut = 10.0^6.0 ; 10.0^8.7
    LBhicut = 10.0^11.0

    LBlo = alog10(LBlocut)>min(xfit)
    LBhi = alog10(LBhicut)<max(xfit)
;   LBlo = alog10(LBlocut/LBnorm)>min(xfit)
;   LBhi = alog10(LBhicut/LBnorm)<max(xfit)

; overlay the running median    

    running = im_medxbin(xfit,yfit,medbin,minx=minx,minpts=minpts,/verbose)

    doit = where((running.binctr gt LBlo) and (running.binctr le LBhi),ndoit)
    lolum = where((running.binctr lt LBlo))

    im_symbols, 108, psize=3.0, thick=postthick, color=medcolor;, /fill
;   oploterror, running.binctr[doit], running.medy[doit], $
;     running.sigy[doit], ps=8, errthick=(postthick-2L)>2L
    oploterror, running.binctr[doit], running.medy[doit], running.binsz/2.0+fltarr(ndoit), $
      running.sigy75[doit]-running.medy[doit], ps=8, errthick=(postthick-1L)>2L, /hibar, $
      color=djs_icolor(medcolor), errcolor=djs_icolor(medcolor)
    oploterror, running.binctr[doit], running.medy[doit], running.binsz/2.0+fltarr(ndoit), $
      running.medy[doit]-running.sigy25[doit], ps=8, errthick=(postthick-1L)>2L, /lobar, $
      color=djs_icolor(medcolor), errcolor=djs_icolor(medcolor)
;   oploterror, running.binctr[doit], running.medy[doit], running.binsz/2.0+fltarr(ndoit), $
;     running.sigyup[doit]-running.medy[doit], ps=8, errthick=(postthick-1L)>2L, /hibar, $
;     color=djs_icolor(medcolor), errcolor=djs_icolor(medcolor)
;   oploterror, running.binctr[doit], running.medy[doit], running.binsz/2.0+fltarr(ndoit), $
;     running.medy[doit]-running.sigylo[doit], ps=8, errthick=(postthick-1L)>2L, /lobar, $
;     color=djs_icolor(medcolor), errcolor=djs_icolor(medcolor)
;   oploterror, running.binctr[doit], running.medy[doit], $
;     running.sigyup[doit]-running.medy[doit], ps=8, errthick=(postthick-1L)>2L, /hibar
;   oploterror, running.binctr[doit], running.medy[doit], $
;     running.medy[doit]-running.sigylo[doit], ps=8, errthick=(postthick-1L)>2L, /lobar

    if keyword_set(paper) then begin
    
; generate a binary FITS table with the statistics

       uband_sfrfile = 'uband_sfr.fits'

       uband_sfr = replicate(template_sfr,ndoit)
       uband_sfr.loglb  = running.binctr[doit]
       uband_sfr.mb     = interpol(xabs,x,running.binctr[doit])
       uband_sfr.p25    = running.sigy25[doit]
       uband_sfr.p50    = running.medy[doit]
       uband_sfr.p75    = running.sigy75[doit]
       uband_sfr.mean   = running.meany[doit]
       uband_sfr.stddev = running.stddev[doit]
       
       mwrfits, uband_sfr, sfrspath+uband_sfrfile, /create

; generate a latex table with the statistics

       openw, lun, latexpath+'uband_sfr.tex', /get_lun
       printf, lun, '\begin{deluxetable}{ccccccc}'
       printf, lun, '\tabletypesize{\small}'
       printf, lun, '\tablecolumns{7}'
       printf, lun, '\tablecaption{U-band Star-Formation Rate Conversion Factors \label{table:uband_sfr}}'
       printf, lun, '\tablewidth{0in}'
       printf, lun, '\tablehead{'
       printf, lun, '\colhead{$\log\,\lb$} & '
       printf, lun, '\colhead{\mb} & '
       printf, lun, '\multicolumn{5}{c}{$\log\,[\sfr/\luobs]$} \\'

       printf, lun, '\colhead{[\lbsun]} & '
       printf, lun, '\colhead{[mag]} & '
       printf, lun, '\multicolumn{5}{c}{$[10^{42}\,\lunits/(\sfrunits)]$} \\'

       printf, lun, '\cline{3-7}'

       printf, lun, '\multicolumn{2}{c}{} & '
       printf, lun, '\colhead{$P_{25}$} & '
       printf, lun, '\colhead{$P_{50}$} & '
       printf, lun, '\colhead{$P_{75}$} & '
       printf, lun, '\colhead{$\langle R\rangle$} & '
       printf, lun, '\colhead{$\sigma_{R}$}'
       printf, lun, '}'
       printf, lun, '\startdata'

       for i = 0L, ndoit-1L do printf, lun, $
         running.binctr[doit[i]], ' & ', interpol(xabs,x,running.binctr[doit[i]]), ' & ', $
         running.sigy25[doit[i]], ' & ', running.medy[doit[i]], ' & ', $
         running.sigy75[doit[i]], ' & ', running.meany[doit[i]], ' & ', $
         running.stddev[doit[i]], ' \\', $
         format='(F5.2,A3,F6.2,A3,F6.3,A3,F6.3,A3,F6.3,A3,F6.3,A3,F6.3,A3)'

       printf, lun, '\enddata'
       printf, lun, '%\tablenotetext{a}{}'
       printf, lun, '\tablecomments{The columns labeled $P_{25}$, $P_{50}$, and $P_{75}$ '+$
         'give the $25$, $50$ (median), and $75$ percentile of the $\sfr/\luobs$ '+$
         'distribution, respectively, in bins of $0.5$~dex in luminosity.  $\langle R\rangle$ '+$
         'and $\sigma_{R}$ give the mean and standard-deviation of the distribution in each bin.}'
       printf, lun, '%\tablerefs{}'
       printf, lun, '\end{deluxetable}'
       free_lun, lun

    endif

    im_openclose, postscript=postscript, /close    

; --------------------------------------------------    
; GENERATE THE WEB PAGE
; --------------------------------------------------    
    
    if keyword_set(postscript) and keyword_set(encapsulated) and (not keyword_set(blackwhite)) then begin
       im_ps2html, htmlbase, html_path=html_path, cleanpng=0, npscols=3, _extra=extra
    endif

; --------------------------------------------------    
; WRITE OUT PLOTS FOR THE PAPER
; --------------------------------------------------    

    if keyword_set(paper) then begin

       splog, 'Writing paper plots to '+paperpath+'.'
       paperplots = [$
         'niiha_vs_oiiihb',$
         'histogram_properties',$
         'sdss_histogram_properties',$
         'sdss_ebv_hahb_vs_ebv_hbhg',$
         'lb_vs_sfr_hair_2panel',$
         'lb_vs_hbha_3panel',$
         'sdss_lb_vs_hbha_3panel',$
         'lb_vs_oiiha_3panel',$
         'sdss_lb_vs_oiiha_3panel',$
         'ehbha_vs_oiiha_2panel',$
         '12oh_oiiinii_niiha_vs_oiiha_2panel',$
         'lb_vs_loiii_sfr_ha_2panel',$
         'd4000_vs_lu_lha_2panel',$
         'sdss_d4000_vs_lu_lha_2panel',$
         'lb_vs_ehbha',$
         'lb_vs_12oh_oiiinii_niiha',$
         'lb_vs_sfr_ha_lhb_obs',$
         'sfr_ha_vs_sfr_hb_2panel',$         
         'lb_vs_sfr_ha_loii_obs',$
         'sfr_ha_vs_sfr_oii_4panel',$
         'highz_lb_vs_hahb',$
         'highz_lb_vs_oiiha',$
         'highz_lb_vs_oiiioii'$
         ]+'.eps'

       for k = 0L, n_elements(paperplots)-1L do $
         spawn, ['/bin/cp -f '+html_path+htmlbase+'/'+paperplots[k]+' '+paperpath], /sh

    endif

stop

; OUT OF DATE

;; ------------------------------------------------------------
;; 3-panel - D4000 vs L(U)/L(Ha) - Integrated
;; ------------------------------------------------------------
;
;    psname = 'd4000_vs_lu_lha_3panel'
;    im_openclose, pspath+psname, postscript=postscript, xsize=7.0, ysize=9.0, encapsulated=encapsulated
;
;    pagemaker, nx=1, ny=3, height=[2.5,2.5,2.5], width=5.5, xmargin=[1.1,0.4], $
;      ymargin=[0.4,1.1], xspace=0, yspace=0, xpage=7.0, ypage=9.0, $
;      position=pos, /normal
;
;    dust_fraction = 0.5
;    
;    indx = where((atlasdust.d4000_narrow_model[1] gt 0.0) and (atlasdust.u_obs gt -900) and $
;      (atlasnodust.ebv_hahb_err gt 0),nindx)
;
;    x = atlasdust[indx].d4000_narrow_model[0]
;    xerr = atlasdust[indx].d4000_narrow_model[1]
;    xlum = atlasdust[indx].b_lum_obs
;
;    indxnfgs = where((nfgsdust.d4000_narrow_model[1] gt 0.0) and (nfgsdust.u_obs gt -900) and $
;      (nfgsnodust.ebv_hahb_err gt 0),nindxnfgs)
;
;    xnfgs = nfgsdust[indxnfgs].d4000_narrow_model[0]
;    xerrnfgs = nfgsdust[indxnfgs].d4000_narrow_model[1]
;    xlumnfgs = nfgsdust[indxnfgs].b_lum_obs
;
;    xtitle = 'D_{n}(4000)'
;    ytitle = 'log [L(U)/L(H\alpha)] '
;
;    xrange = D4000range2
;    yrange = Uhacorrange
;    
;; ##########################################################
;; Panel 1: U observed, Ha observed
;; ##########################################################
;
;    y1 = Uconstant*10^(-0.4*atlasdust[indx].synth_u_obs)
;    y1err = 0.4*atlasdust[indx].synth_u_obs_err*alog(10.0)*y1
;    
;    y2 = atlasdust[indx].h_alpha[0]
;    y2err = atlasdust[indx].h_alpha[1]
;
;    y = alog10(y1/y2)
;    yerr = im_compute_error(y1,y1err,y2,y2err,/log)
;
;    y1nfgs = Uconstant*10^(-0.4*nfgsdust[indxnfgs].synth_u_obs)
;    y1errnfgs = 0.4*nfgsdust[indxnfgs].synth_u_err*alog(10.0)*y1
;    
;    y2nfgs = nfgsdust[indxnfgs].h_alpha[0]
;    y2errnfgs = nfgsdust[indxnfgs].h_alpha[1]
;
;    ynfgs = alog10(y1nfgs/y2nfgs)
;    yerrnfgs = im_compute_error(y1nfgs,y1errnfgs,y2nfgs,y2errnfgs,/log)
;    
;    stats = im_stats([y,ynfgs],/verbose,/baremin)
;    xstr = strtrim(string(stats.median_rej,format='(F12.2)'),2)+' ('+$
;      strtrim(string(stats.mean_rej,format='(F12.2)'),2)+$
;      '\pm'+strtrim(string(stats.sigma_rej,format='(F12.2)'),2)+')'
;    rcor = r_correlate([x,xnfgs],[y,ynfgs],zd=zd,probd=probd)
;    splog, 'Spearman rank test: ', rcor, probd, zd
;
;    ytitle = 'log [L(U)/L(H\alpha)]_{obs} '
;
;    atlas1d_lineplot, x, y, xerr, yerr, plottype=3, postscript=postscript, $
;      ytitle=ytitle, xrange=xrange, yrange=yrange, legendtype=0, $
;      /right, /top, xnfgs=xnfgs, ynfgs=ynfgs, xerrnfgs=xerrnfgs, blackwhite=blackwhite, $
;      yerrnfgs=yerrnfgs, position=pos[*,0], charsize=charsize_5, xtickname=replicate(' ',10)
;
;    running = im_medxbin([x,xnfgs],[y,ynfgs],0.1,minpts=10)
;
;    legend, '(a)', /left, /top, box=0, charsize=charsize_3, charthick=postthick
;    legend, textoidl(xstr), /right, /top, box=0, charsize=charsize_3, charthick=postthick
;
;; ##########################################################
;; Panel 2: U corrected, Ha corrected
;; ##########################################################
;
;;   AU = atlasnodust[indx].ebv_hahb*k_lambda(Uinfo.weff,/charlot)*dust_fraction
;;   AU_err = atlasnodust[indx].ebv_hahb_err*k_lambda(Uinfo.weff,/charlot)*dust_fraction
;;   y1 = Uconstant*10^(-0.4*atlasdust[indx].synth_u_obs)*10^(0.4*AU)
;;   y1err = sqrt((0.4*atlasdust[indx].synth_u_obs_err*alog(10.0)*y1*10^(0.4*AU))^2 + $
;;     (y1*alog(10.0)*AU_err)^2)
;;   y2 = atlasnodust[indx].h_alpha[0]
;;   y2err = atlasnodust[indx].h_alpha[1]
;;   y = alog10(y1/y2)
;;   yerr = im_compute_error(y1,y1err,y2,y2err,/log)
;
;    AU = atlasnodust[indx].ebv_hahb*kl_u*dust_fraction
;    u_dust = Uconstant*10^(-0.4*atlasdust[indx].synth_u_obs)
;    u_dust_err = atlasdust[indx].synth_u_obs_err*u_dust*alog(10.0)/2.5
;    u_nodust = u_dust*10^(0.4*AU)
;
;    hb_dust = atlasdust[indx].h_beta[0]
;    hb_dust_err = atlasdust[indx].h_beta[1]
;
;    ha_dust = atlasdust[indx].h_alpha[0]
;    ha_dust_err = atlasdust[indx].h_alpha[1]
;    ha_nodust = atlasnodust[indx].h_alpha[0]
;
;    y1 = u_nodust/ha_nodust
;    y1err = sqrt((y1*(1+bigbprime_u-bigbprime_ha)*(ha_dust_err/ha_dust))^2 + $
;      (y1*(bigbprime_u-bigbprime_ha)*(hb_dust_err/hb_dust))^2 + (y1*u_dust_err/u_dust)^2)
;    y = alog10(y1)
;    yerr = y1err/y1/alog(10.0)    
;
;; NFGS    
;    
;;   y1nfgs = Uconstant*10^(-0.4*nfgsdust[indxnfgs].synth_u_obs)*10^(0.4*AUnfgs)
;;   y1errnfgs = sqrt((0.4*nfgsdust[indxnfgs].synth_u_obs_err*alog(10.0)*y1nfgs*10^(0.4*AUnfgs))^2 + $
;;     (y1nfgs*alog(10.0)*AUnfgs_err)^2)
;;   y2nfgs = nfgsnodust[indxnfgs].h_alpha[0]
;;   y2errnfgs = nfgsnodust[indxnfgs].h_alpha[1]
;;   ynfgs = alog10(y1nfgs/y2nfgs)
;;   yerrnfgs = im_compute_error(y1nfgs,y1errnfgs,y2nfgs,y2errnfgs,/log)
;    
;    AUnfgs = nfgsnodust[indxnfgs].ebv_hahb*kl_u*dust_fraction
;    unfgs_dust = Uconstant*10^(-0.4*nfgsdust[indxnfgs].synth_u_obs)
;    unfgs_dust_err = nfgsdust[indxnfgs].synth_u_obs_err*unfgs_dust*alog(10.0)/2.5
;    unfgs_nodust = unfgs_dust*10^(0.4*AUnfgs)
;
;    hbnfgs_dust = nfgsdust[indxnfgs].h_beta[0]
;    hbnfgs_dust_err = nfgsdust[indxnfgs].h_beta[1]
;
;    hanfgs_dust = nfgsdust[indxnfgs].h_alpha[0]
;    hanfgs_dust_err = nfgsdust[indxnfgs].h_alpha[1]
;    hanfgs_nodust = nfgsnodust[indxnfgs].h_alpha[0]
;
;    y1nfgs = unfgs_nodust/hanfgs_nodust
;    y1errnfgs = sqrt((y1nfgs*(1+bigbprime_u-bigbprime_ha)*(hanfgs_dust_err/hanfgs_dust))^2 + $
;      (y1nfgs*(bigbprime_u-bigbprime_ha)*(hbnfgs_dust_err/hbnfgs_dust))^2 + (y1nfgs*unfgs_dust_err/unfgs_dust)^2)
;    
;    ynfgs = alog10(y1nfgs)
;    yerrnfgs = y1errnfgs/y1nfgs/alog(10.0)
;    
;; combine the samples and give the SFR conversion ratio
;
;    ratio = loghaconst+[y,ynfgs]
;    stats = im_stats(ratio)
;    constant = 10^stats.median
;    constant_err = constant*alog(10.0)*stats.sig68mean
;
;    splog, 'SFR(Ha)/L(U)_cor coefficients: '
;    print, constant, constant_err
;
;; now move on with your life    
;    
;    stats = im_stats([y,ynfgs],/verbose,/baremin,/no_head)
;    xstr = strtrim(string(stats.median_rej,format='(F12.2)'),2)+' ('+$
;      strtrim(string(stats.mean_rej,format='(F12.2)'),2)+$
;      '\pm'+strtrim(string(stats.sigma_rej,format='(F12.2)'),2)+')'
;    rcor = r_correlate([x,xnfgs],[y,ynfgs],zd=zd,probd=probd)
;    splog, 'Spearman rank test: ', rcor, probd, zd
;
;    ytitle = 'log [L(U)/L(H\alpha)]_{cor}'
;
;    atlas1d_lineplot, x, y, xerr, yerr, plottype=3, postscript=postscript, $
;      xtitle='', ytitle=ytitle, xrange=xrange, yrange=yrange, legendtype=0, $
;      /noerase, xstyle=3, /right, /top, xtickname=replicate(' ',10), $
;      xnfgs=xnfgs, ynfgs=ynfgs, xerrnfgs=xerrnfgs, blackwhite=blackwhite, $
;      yerrnfgs=yerrnfgs, position=pos[*,1], charsize=charsize_5
;    
;    running = im_medxbin([x,xnfgs],[y,ynfgs],0.1,minpts=10)
;
;    legend, '(b)', /left, /top, box=0, charsize=charsize_3, charthick=postthick
;    legend, textoidl(xstr), /right, /top, box=0, charsize=charsize_3, charthick=postthick
;    legend, textoidl('E(B-V)_{stars}/E(B-V)_{gas} = 0.5'), /left, /bottom, box=0, charsize=charsize_3, charthick=postthick
;;   legend, textoidl('(b) E(B-V)_{stars}/E(B-V)_{gas} = 0.5'), /left, /top, box=0, charsize=charsize_3, charthick=postthick
;
;; overlay the running median    
;
;    running = im_medxbin([x,xnfgs],[y,ynfgs],0.05,minpts=3)
;
;    splog, 'Scatter about the running median:'
;    runningy = interpol(running.medy,running.binctr,[x,xnfgs])
;    stats = im_stats([y,ynfgs]-runningy,/verbose,/baremin,/no_head)
;
;; ##########################################################
;; Panel 3: U observed, Ha corrected
;; ##########################################################
;
;    u_dust = Uconstant*10^(-0.4*atlasdust[indx].synth_u_obs)
;    u_dust_err = atlasdust[indx].synth_u_obs_err*u_dust*alog(10.0)/2.5
;
;    hb_dust = atlasdust[indx].h_beta[0]
;    hb_dust_err = atlasdust[indx].h_beta[1]
;
;    ha_dust = atlasdust[indx].h_alpha[0]
;    ha_dust_err = atlasdust[indx].h_alpha[1]
;    ha_nodust = atlasnodust[indx].h_alpha[0]
;    ha_nodust_err = sqrt((ha_nodust*(1-bigbprime_ha)*ha_dust_err/ha_dust)^2 + (bigbprime_ha*ha_nodust*hb_dust_err/hb_dust)^2)
;
;    y = alog10(u_dust/ha_nodust) - loghaconst - Uhasfrconstoffset
;    yerr = im_compute_error(u_dust,u_dust_err,ha_nodust,ha_nodust_err,/log)
;    
;; the code below was wrong!!!  the U-band luminosity was normalized to
;; the U-band solar luminosity, while the H-alpha luminosity was
;; normalized to the bolometric luminosity!
;
;;   cor = (atlasdust[indx].synth_u-atlasdust[indx].synth_u_obs) > 0.0 ; <-- NOTE!
;;   y1 = atlasdust[indx].synth_u_lum + 0.4*(mbolsun-uinfo.solarmags)
;;;  y1 = atlasdust[indx].synth_u_lum - 0.4*cor + 0.4*(mbolsun-uinfo.solarmags)
;;   y1err = atlasdust[indx].synth_u_lum_err
;;   y2 = atlasnodust[indx].h_alpha_lum[0]
;;   y2err = atlasnodust[indx].h_alpha_lum[1]
;;   y = y1 - y2 - loghaconst - Uhasfrconstoffset
;;   yerr = sqrt(y1err^2 + y2err^2)
;
;;;  y2 = atlasdust[indx].synth_u_lum
;;;  y3 = atlasdust[indx].synth_u_lum + 0.4*(mbolsun-uinfo.solarmags)
;;;  y4 = 0.4*(mbolsun-atlasdust[indx].synth_M_u)   
;    
;; NFGS    
;    
;    unfgs_dust = Uconstant*10^(-0.4*nfgsdust[indxnfgs].synth_u_obs)
;    unfgs_dust_err = nfgsdust[indxnfgs].synth_u_obs_err*unfgs_dust*alog(10.0)/2.5
;
;    hbnfgs_dust = nfgsdust[indxnfgs].h_beta[0]
;    hbnfgs_dust_err = nfgsdust[indxnfgs].h_beta[1]
;
;    hanfgs_dust = nfgsdust[indxnfgs].h_alpha[0]
;    hanfgs_dust_err = nfgsdust[indxnfgs].h_alpha[1]
;    hanfgs_nodust = nfgsnodust[indxnfgs].h_alpha[0]
;    hanfgs_nodust_err = sqrt((hanfgs_nodust*(1-bigbprime_ha)*hanfgs_dust_err/hanfgs_dust)^2 + $
;      (bigbprime_ha*hanfgs_nodust*hbnfgs_dust_err/hbnfgs_dust)^2)
;
;    ynfgs = alog10(unfgs_dust/hanfgs_nodust) - loghaconst - Uhasfrconstoffset
;    yerrnfgs = im_compute_error(unfgs_dust,unfgs_dust_err,hanfgs_nodust,hanfgs_nodust_err,/log)
;
;;   cor = (nfgsdust[indxnfgs].synth_u-nfgsdust[indxnfgs].synth_u_obs) > 0.0 ; <-- NOTE!
;;   y1nfgs = nfgsdust[indxnfgs].synth_u_lum - 0.4*cor + 0.4*(mbolsun-uinfo.solarmags)
;;   y1errnfgs = nfgsdust[indxnfgs].synth_u_lum_err
;;   y2nfgs = nfgsnodust[indxnfgs].h_alpha_lum[0]
;;   y2errnfgs = nfgsnodust[indxnfgs].h_alpha_lum[1]
;;   ynfgs = y1nfgs - y2nfgs - loghaconst - Uhasfrconstoffset
;;   yerrnfgs = sqrt(y1errnfgs^2 + y2errnfgs^2)
;    
;; stats
;    
;    stats = im_stats([y,ynfgs],/verbose)
;    xstr = strtrim(string(stats.median_rej,format='(F12.2)'),2)+' ('+$
;      strtrim(string(stats.mean_rej,format='(F12.2)'),2)+$
;      '\pm'+strtrim(string(stats.sigma_rej,format='(F12.2)'),2)+')'
;    rcor = r_correlate([x,xnfgs],[y,ynfgs],zd=zd,probd=probd)
;    splog, 'Spearman rank test: ', rcor, probd, zd
;
;; derive the SFR conversion factor for the full sample
;
;    constant = 10.0^(-(stats.median_rej+Uhasfrconstoffset))
;    constant_err = constant*alog(10.0)*stats.sigma_rej
;    
;    splog, 'SFR(U) coefficients: '
;    print, constant, constant_err
;
;; derive the SFR conversion factor for just luminous galaxies 
;
;    w = where([xlum,xlumnfgs] gt 9.5)
;    lumstats = im_stats(([y,ynfgs])[w])
;    
;    constant2 = 10.0^(-(lumstats.median_rej+Uhasfrconstoffset))
;    constant2_err = constant*alog(10.0)*lumstats.sigma_rej
;    
;    splog, 'SFR(U) coefficients for log LB > 9.5: '
;    print, constant2, constant2_err
;
;; now move on with your life    
;
;    ytitle = 'log [10^{-43} L(U)_{obs}/\psi(H\alpha)]'
;    yrange2 = sfrUharange - Uhasfrconstoffset
;
;    atlas1d_lineplot, x, y, xerr, yerr, plottype=3, postscript=postscript, $
;      xrange=xrange, yrange=yrange2, legendtype=0, /noerase, xtitle=xtitle, $
;      /right, /top, xnfgs=xnfgs, ynfgs=ynfgs, xerrnfgs=xerrnfgs, $
;      yerrnfgs=yerrnfgs, position=pos[*,2], charsize=charsize_5, ytitle=ytitle, blackwhite=blackwhite
;
;    running = im_medxbin([x,xnfgs],[y,ynfgs],0.1,minpts=10)
;
;; overplot the Cram et al. SFR calibration
;
;;   djs_oplot, !x.crange, -(alog10(HaHb)+loghaconst)*[1,1], line=0, thick=(postthick-4)>2
;
;    legend, '(c)', /left, /top, box=0, charsize=charsize_3, charthick=postthick
;    legend, textoidl(xstr), /right, /top, box=0, charsize=charsize_3, charthick=postthick
;
;    im_openclose, postscript=postscript, /close
;
;stop
;    
;; ------------------------------------------------------------
;; 3-panel - D4000 vs L(U)/L(Ha) - SDSS
;; ------------------------------------------------------------
;
;    psname = 'sdss_d4000_vs_lu_lha_3panel'
;    im_openclose, pspath+psname, postscript=postscript, xsize=7.0, ysize=9.0, encapsulated=encapsulated
;
;    pagemaker, nx=1, ny=3, height=[2.5,2.5,2.5], width=5.5, xmargin=[1.1,0.4], $
;      ymargin=[0.4,1.1], xspace=0, yspace=0, xpage=7.0, ypage=9.0, $
;      position=pos, /normal
;
;    dust_fraction = 0.5
;
;    indx = where((sdssdust.model_d4000_narrow[1] gt 0.0) and (sdssancillary.fiber_u gt -900) and $
;      finite(sdssancillary.fiber_u) and (sdssnodust.ebv_hahb_err gt 0),nindx)
;
;    x = sdssdust[indx].model_d4000_narrow[0]
;    xerr = sdssdust[indx].model_d4000_narrow[1]
;
;    xtitle = 'D_{n}(4000)'
;    ytitle = 'log [L(U)/L(H\alpha)] '
;
;    xrange = D4000range2
;    yrange = Uhacorrange
;    
;; ##########################################################
;; Panel 1: U observed, Ha observed
;; ##########################################################
;
;    y1 = Uconstant*10^(-0.4*sdssancillary[indx].fiber_U)
;    y1err = 0.4*sdssancillary[indx].fiber_U_err*alog(10.0)*y1
;    
;    y2 = sdssdust[indx].h_alpha[0]
;    y2err = sdssdust[indx].h_alpha[1]
;
;    y = alog10(y1/y2)
;    yerr = im_compute_error(y1,y1err,y2,y2err,/log)
;
;    stats = im_stats(y,/verbose,/baremin)
;    xstr = strtrim(string(stats.median_rej,format='(F12.2)'),2)+' ('+$
;      strtrim(string(stats.mean_rej,format='(F12.2)'),2)+$
;      '\pm'+strtrim(string(stats.sigma_rej,format='(F12.2)'),2)+')'
;    rcor = r_correlate(x,y,zd=zd,probd=probd)
;    splog, 'Spearman rank test: ', rcor, probd, zd
;
;    ytitle = 'log [L(U)/L(H\alpha)]_{obs} '
;
;    sdss_lineplot, x, y, xerr, yerr, plottype=3, postscript=postscript, $
;      ytitle=ytitle, xrange=xrange, yrange=yrange, legendtype=0, $
;      /right, /top, position=pos[*,0], charsize=charsize_5, xtickname=replicate(' ',10)
;
;    running = im_medxbin(x,y,0.1,minpts=100)
;
;    legend, '(a)', /left, /top, box=0, charsize=charsize_3, charthick=postthick
;    legend, textoidl(xstr), /right, /top, box=0, charsize=charsize_3, charthick=postthick
;
;; ##########################################################
;; Panel 2: U corrected, Ha corrected
;; ##########################################################
;
;    AU = sdssnodust[indx].ebv_hahb*k_lambda(Uinfo.weff,/charlot)*dust_fraction
;    AU_err = sdssnodust[indx].ebv_hahb_err*k_lambda(Uinfo.weff,/charlot)*dust_fraction
;
;    y1 = Uconstant*10^(-0.4*sdssancillary[indx].fiber_U)*10^(0.4*AU)
;    y1err = sqrt((0.4*sdssancillary[indx].fiber_U_err*alog(10.0)*y1*10^(0.4*AU))^2 + $
;      (y1*alog(10.0)*AU_err)^2)
; 
;    y2 = sdssnodust[indx].h_alpha[0]
;    y2err = sdssnodust[indx].h_alpha[1]
; 
;    y = alog10(y1/y2)
;    yerr = im_compute_error(y1,y1err,y2,y2err,/log)
;
;    stats = im_stats(y,/verbose,/baremin,/no_head)
;    xstr = strtrim(string(stats.median_rej,format='(F12.2)'),2)+' ('+$
;      strtrim(string(stats.mean_rej,format='(F12.2)'),2)+$
;      '\pm'+strtrim(string(stats.sigma_rej,format='(F12.2)'),2)+')'
;    rcor = r_correlate(x,y,zd=zd,probd=probd)
;    splog, 'Spearman rank test: ', rcor, probd, zd
;
;    ytitle = 'log [L(U)/L(H\alpha)]_{cor}'
;
;    sdss_lineplot, x, y, xerr, yerr, plottype=3, postscript=postscript, xtickname=replicate(' ',10), $
;      xtitle='', ytitle=ytitle, xrange=xrange, yrange=yrange, legendtype=0, /noerase, $
;      xstyle=3, /right, /top, position=pos[*,1], charsize=charsize_5
;    
;    running = im_medxbin(x,y,0.1,minpts=100)
;
;;   legend, textoidl('(b) E(B-V)_{stars}/E(B-V)_{gas} = 0.5'), /left, /top, box=0, charsize=charsize_3, charthick=postthick
;    legend, '(b)', /left, /top, box=0, charsize=charsize_3, charthick=postthick
;    legend, textoidl(xstr), /right, /top, box=0, charsize=charsize_3, charthick=postthick
;    legend, textoidl('E(B-V)_{stars}/E(B-V)_{gas} = 0.5'), /left, /bottom, box=0, charsize=charsize_3, charthick=postthick
;
;;   label = textoidl(['A_{stars}/A_{gas} = 0.5'])
;;   legend, label, /left, /bottom, box=0, charsize=charsize_3, charthick=postthick
;
;; overlay the running median    
;
;    running = im_medxbin(x,y,0.05,minpts=100)
;
;; ##########################################################
;; Panel 3: U observed, Ha corrected
;; ##########################################################
;
;    y1 = sdssancillary[indx].fiber_U_lum
;    y1err = sdssancillary[indx].fiber_U_lum_err
;    
;    y2 = sdssnodust[indx].h_alpha_lum[0]
;    y2err = sdssnodust[indx].h_alpha_lum[1]
;
;    y = y1 - y2 - loghaconst - Uhasfrconstoffset
;    yerr = sqrt(y1err^2 + y2err^2)
;
;    stats = im_stats(y,/verbose)
;    xstr = strtrim(string(stats.median_rej,format='(F12.2)'),2)+' ('+$
;      strtrim(string(stats.mean_rej,format='(F12.2)'),2)+$
;      '\pm'+strtrim(string(stats.sigma_rej,format='(F12.2)'),2)+')'
;    rcor = r_correlate(x,y,zd=zd,probd=probd)
;    splog, 'Spearman rank test: ', rcor, probd, zd
;
;    ytitle = 'log [10^{-43} L(U)_{obs}/\psi(H\alpha)]'
;    yrange2 = sfrUharange - Uhasfrconstoffset
;
;    sdss_lineplot, x, y, xerr, yerr, plottype=3, postscript=postscript, xtitle=xtitle, $
;      xrange=xrange, yrange=yrange2, legendtype=0, /noerase, $
;      /right, /top, position=pos[*,2], charsize=charsize_5, ytitle=ytitle
;    
;    running = im_medxbin(x,y,0.1,minpts=100)
;
;    legend, '(c)', /left, /top, box=0, charsize=charsize_3, charthick=postthick
;    legend, textoidl(xstr), /right, /top, box=0, charsize=charsize_3, charthick=postthick
;
;    im_openclose, postscript=postscript, /close

;; ------------------------------------------------------------
;; D4000 vs [L(U)/L(Ha)]_cor - Integrated for different values of eta 
;; ------------------------------------------------------------
;
;    psname = 'd4000_vs_lulhacor_10panel'
;    im_openclose, pspath+psname, postscript=postscript, xsize=8.5, ypage=11.5, encapsulated=encapsulated
;
;    pagemaker, nx=1, ny=10, width=6.0, xmargin=[1.1,0.4], height=replicate(1.0,10), $
;      ymargin=[0.4,1.1], xspace=0, yspace=0, xpage=8.5, ypage=11.5, $
;      position=pos, /normal
;
;    dust_fraction = findgen(10)/10+0.1
;    ndust = n_elements(dust_fraction)
;    
;    indx = where((atlasdust.d4000_narrow_model[1] gt 0.0) and (atlasdust.u_obs gt -900) and $
;      (atlasnodust.ebv_hahb_err gt 0),nindx)
;
;    x = atlasdust[indx].d4000_narrow_model[0]
;    xerr = atlasdust[indx].d4000_narrow_model[1]
;    xlum = atlasdust[indx].b_lum_obs
;
;    indxnfgs = where((nfgsdust.d4000_narrow_model[1] gt 0.0) and (nfgsdust.u_obs gt -900) and $
;      (nfgsnodust.ebv_hahb_err gt 0),nindxnfgs)
;
;    xnfgs = nfgsdust[indxnfgs].d4000_narrow_model[0]
;    xerrnfgs = nfgsdust[indxnfgs].d4000_narrow_model[1]
;    xlumnfgs = nfgsdust[indxnfgs].b_lum_obs
;
;    xrange = D4000range2
;    yrange = Uhacorrange
;
;    xtickname = replicate(' ',10)
;    
;    for idust = 0L, ndust-1L do begin
;
;       if (idust eq ndust-1L) then begin
;          delvarx, xtickname
;          xtitle = 'D_{n}(4000)'          
;       endif else xtitle = ''
;       
;       if (idust eq 1L) or (idust eq (ndust-1L)/2L) or (idust eq (ndust-2L)) then begin
;          ytitle = 'log [L(U)/L(H\alpha)]_{cor}'
;       endif else ytitle = ''
;       
;       AU = atlasnodust[indx].ebv_hahb*k_lambda(Uinfo.weff,/charlot)*dust_fraction[idust]
;       AU_err = atlasnodust[indx].ebv_hahb_err*k_lambda(Uinfo.weff,/charlot)*dust_fraction[idust]
;
;       y1 = Uconstant*10^(-0.4*atlasdust[indx].synth_u_obs)*10^(0.4*AU)
;       y1err = sqrt((0.4*atlasdust[indx].synth_u_obs_err*alog(10.0)*y1*10^(0.4*AU))^2 + $
;                    (y1*alog(10.0)*AU_err)^2)
;       
;       y2 = atlasnodust[indx].h_alpha[0]
;       y2err = atlasnodust[indx].h_alpha[1]
;       
;       y = alog10(y1/y2)
;       yerr = im_compute_error(y1,y1err,y2,y2err,/log)
;
;; NFGS    
;       
;       AUnfgs = nfgsnodust[indxnfgs].ebv_hahb*k_lambda(Uinfo.weff,/charlot)*dust_fraction[idust]
;       AUnfgs_err = nfgsnodust[indxnfgs].ebv_hahb_err*k_lambda(Uinfo.weff,/charlot)*dust_fraction[idust]
;
;       y1nfgs = Uconstant*10^(-0.4*nfgsdust[indxnfgs].synth_u_obs)*10^(0.4*AUnfgs)
;       y1errnfgs = sqrt((0.4*nfgsdust[indxnfgs].synth_u_obs_err*alog(10.0)*y1nfgs*10^(0.4*AUnfgs))^2 + $
;                        (y1nfgs*alog(10.0)*AUnfgs_err)^2)
;       
;       y2nfgs = nfgsnodust[indxnfgs].h_alpha[0]
;       y2errnfgs = nfgsnodust[indxnfgs].h_alpha[1]
;       
;       ynfgs = alog10(y1nfgs/y2nfgs)
;       yerrnfgs = im_compute_error(y1nfgs,y1errnfgs,y2nfgs,y2errnfgs,/log)
;       
;       stats = im_stats([y,ynfgs],/verbose,/baremin,/no_head)
;       xstr = strtrim(string(stats.median_rej,format='(F12.2)'),2)+' ('+$
;         strtrim(string(stats.mean_rej,format='(F12.2)'),2)+$
;         '\pm'+strtrim(string(stats.sigma_rej,format='(F12.2)'),2)+')'
;       rcor = r_correlate([x,xnfgs],[y,ynfgs],zd=zd,probd=probd)
;
;       atlas1d_lineplot, x, y, xerr, yerr, plottype=1, postscript=postscript, $
;         xtitle=xtitle, ytitle=ytitle, xrange=xrange, yrange=yrange, legendtype=0, $
;         xstyle=3, /right, /top, xtickname=xtickname, $
;         xnfgs=xnfgs, ynfgs=ynfgs, xerrnfgs=xerrnfgs, blackwhite=blackwhite, $
;         yerrnfgs=yerrnfgs, position=pos[*,idust], charsize=charsize_0, $
;         noerase=(idust gt 0L)
;       
;       running = im_medxbin([x,xnfgs],[y,ynfgs],0.1,minpts=10)
;
;       legend, textoidl(xstr), /right, /bottom, box=0, charsize=0.8, charthick=postthick3
;       legend, string(dust_fraction[idust],format='(F3.1)'), $
;         /left, /bottom, box=0, charsize=charsize_0, charthick=postthick3
;
;; overlay the running median    
;
;;      running = im_medxbin([x,xnfgs],[y,ynfgs],0.05,minpts=3)
;
;       splog, 'Scatter about the running median:'
;       runningy = interpol(running.medy,running.binctr,[x,xnfgs])
;       stats = im_stats([y,ynfgs]-runningy,/verbose,/baremin,/no_head)
;
;    endfor
;       
;    im_openclose, postscript=postscript, /close

; ------------------------------------------------------------
; SFR(Ha) vs various SFR(U) calibrations - SDSS
; ------------------------------------------------------------
    
    psname = 'sdss_sfr_ha_vs_sfr_uband_2panel'
    im_openclose, pspath+psname, postscript=postscript, xsize=8.5, ysize=8.0, encapsulated=encapsulated

    pagemaker, nx=1, ny=1, xspace=0, yspace=0.0, width=6.8, height=6.8, $
      xmargin=[1.4,0.3], ymargin=[0.9,1.1], xpage=8.5, ypage=8.8, $
      position=pos, /normal

    xrange = sfrharange
    yrange = residrange_Usfrs

    xtitle = 'log \psi(H\alpha) ['+sfr_units()+']'
    ytitle = 'log \psi(U)/\psi(H\alpha)'
    
    indx = where((sdssancillary.fiber_u_lum gt -900.0) and (sdssancillary.b_lum gt -900) and $
      (sdssnodust.sfr_h_alpha gt -900.0),nindx)

    x = sdssnodust[indx].sfr_h_alpha
    xerr = sdssnodust[indx].sfr_h_alpha_err
    
    y1 = x
    y1err = xerr
    
    lubandobs = sdssancillary[indx].fiber_u_lum[0] + alog10(lusun)
    loglb = sdssancillary[indx].b_lum

    y2 = alog10(uband_sfr(loglb,lubandobs))    
    y2err = y2*0.0

    resid = y2-y1
    resid_err = sqrt(y1err^2+y2err^2)

; make the plot    

    stats = im_stats(resid,/verbose,/baremin)
    xstr = strtrim(string(stats.median_rej,format='(F12.2)'),2)+' ('+$
      strtrim(string(stats.mean_rej,format='(F12.2)'),2)+$
      '\pm'+strtrim(string(stats.sigma_rej,format='(F12.2)'),2)+')'
    
    sdss_lineplot, x, resid, xerr, resid_err, plottype=1, postscript=postscript, $
      xtitle=xtitle, ytitle=ytitle, xrange=xrange, yrange=yrange, legendtype=0, $
      charsize=charsize_7, position=pos[*,0]
    djs_oplot, !x.crange, [0,0], line=2, thick=postthick

    legend, textoidl(xstr), /right, /top, box=0, charsize=charsize_3, charthick=postthick

    im_openclose, postscript=postscript, /close    

stop
    
; ------------------------------------------------------------
; L(B) versus SFR(Ha)/L(U)_obs - SDSS
; ------------------------------------------------------------
    
    psname = 'sdss_lb_vs_sfr_ha_lu_obs'
    im_openclose, pspath+psname, postscript=postscript, xsize=8.5, ysize=8.0, encapsulated=encapsulated

    pagemaker, nx=1, ny=1, xspace=0, yspace=0.0, width=6.8, height=6.8, $
      xmargin=[1.4,0.3], ymargin=[0.9,1.1], xpage=8.5, ypage=8.8, $
      position=pos, /normal

    if keyword_set(blackwhite) then begin
       medcolor = 'grey'
    endif else begin
       medcolor = 'red'
    endelse

    medbin = 1.0
    minpts = 5
    minx = 9.0
    
    indx = where((sdssancillary.b_lum gt -900) and (sdssnodust.sfr_h_alpha gt -900.0) and $
      (sdssancillary.fiber_u_lum gt -900.0),nindx)

    x = sdssancillary[indx].b_lum; - alog10(LBnorm)
    xerr = sdssancillary[indx].b_lum_err
    xabs = sdssancillary[indx].m_b

    uband = sdssancillary[indx].fiber_u_lum + alog10(lusun) - alog10(Ulumnorm)
    uband_err = sdssancillary[indx].fiber_u_lum_err

    sfr_ha = sdssnodust[indx].sfr_h_alpha
    sfr_ha_err = sdssnodust[indx].sfr_h_alpha_err

    y = sfr_ha - uband
    yerr = sqrt(sfr_ha_err^2 + uband_err^2)

; make the plot

    xtitle = 'log L(B) [L(B)'+sunsymbol()+']'
    ytitle = 'log [10^{42} \psi(H\alpha)/L(U)_{obs}] [erg s^{-1}/'+sfr_units()+']'

    xrange = LBrange; - alog10(LBnorm)
    yrange = alog10([1D-45,1D-41]) - alog10(sfrnorm) + alog10(Ulumnorm)
   
    sdss_lineplot, x, y, xerr, yerr, plottype=1, postscript=postscript, $
      xtitle=xtitle, ytitle=ytitle, xrange=xrange, $
      xstyle=11, yrange=yrange, legendtype=0, /right, /top, yfrac=1.5, xfrac=1.5, charsize=charsize_8, $
      position=pos[*,0]
    axis, /xaxis, xrange=interpol(xabs,x,!x.crange), xthick=postthick, xstyle=1, $
      charsize=charsize_8, charthick=postthick
    im_xyouts_title, xtitle='M_{B} [mag]', charsize=charsize_8, charthick=postthick

; fit only to certain data

    xfit = x
    yfit = y
    
; 9.0 is M_B = -17.07 and 8.7 is M_B = -16.3; 8.57 is -16

    LBlocut = 10.0^6.0 ; 10.0^8.7
    LBhicut = 10.0^11.0

    LBlo = alog10(LBlocut)>min(xfit)
    LBhi = alog10(LBhicut)<max(xfit)

; overlay the running median    

    running = im_medxbin(xfit,yfit,medbin,minx=minx,minpts=minpts,/verbose)

    doit = where((running.binctr gt LBlo) and (running.binctr le LBhi),ndoit)
    lolum = where((running.binctr lt LBlo))

    im_symbols, 108, psize=3.0, thick=postthick, color=medcolor;, /fill
    oploterror, running.binctr[doit], running.medy[doit], running.binsz/2.0+fltarr(ndoit), $
      running.sigy75[doit]-running.medy[doit], ps=8, errthick=(postthick-1L)>2L, /hibar, $
      color=djs_icolor(medcolor), errcolor=djs_icolor(medcolor)
    oploterror, running.binctr[doit], running.medy[doit], running.binsz/2.0+fltarr(ndoit), $
      running.medy[doit]-running.sigy25[doit], ps=8, errthick=(postthick-1L)>2L, /lobar, $
      color=djs_icolor(medcolor), errcolor=djs_icolor(medcolor)

    im_openclose, postscript=postscript, /close    

; ------------------------------------------------------------
; SFR(Ha) vs various SFR(U) calibrations - Integrated
; ------------------------------------------------------------
    
    psname = 'sfr_ha_vs_sfr_uband_2panel'
    im_openclose, pspath+psname, postscript=postscript, xsize=7.0, ysize=7.5, encapsulated=encapsulated

    pagemaker, nx=1, ny=2, height=3.0*[1,1], width=5.5, xmargin=[1.1,0.4], $
      ymargin=[0.4,1.1], xspace=0, yspace=0, xpage=7.0, ypage=7.5, $
      position=pos, /normal

    xrange = sfrharange
    yrange = residrange_Usfrs

    xtitle = 'log \psi(H\alpha) ['+sfr_units()+']'
    ytitle = 'log \psi(U)/\psi(H\alpha)'
    
; ##########################################################
; Panel 1
; ##########################################################

; SF galaxies
    
    indx = where((sfrs.sfr_ha gt -900) and (sfrs.sfr_uband_cram gt -900),nindx)

    x = sfrs[indx].sfr_ha
    xerr = sfrs[indx].sfr_ha_err
    
    y1 = sfrs[indx].sfr_ha
    y1err = sfrs[indx].sfr_ha_err
    
    y2 = sfrs[indx].sfr_uband_cram
    y2err = sfrs[indx].sfr_uband_cram_err

    resid = y2-y1
    resid_err = sqrt(y1err^2+y2err^2)

; AGN
    
    indxagn = where((sfrs_agn.sfr_ha gt -900) and (sfrs_agn.sfr_uband_cram gt -900),nindxagn)

    xagn = sfrs_agn[indxagn].sfr_ha
    xerragn = sfrs_agn[indxagn].sfr_ha_err
    
    y1agn = sfrs_agn[indxagn].sfr_ha
    y1erragn = sfrs_agn[indxagn].sfr_ha_err
    
    y2agn = sfrs_agn[indxagn].sfr_uband_cram
    y2erragn = sfrs_agn[indxagn].sfr_uband_cram_err

    residagn = y2agn-y1agn
    residagn_err = sqrt(y1erragn^2+y2erragn^2)

; make the plot    

    stats = im_stats(resid,/verbose,/baremin)
    xstr = strtrim(string(stats.median_rej,format='(F12.2)'),2)+' ('+$
      strtrim(string(stats.mean_rej,format='(F12.2)'),2)+$
      '\pm'+strtrim(string(stats.sigma_rej,format='(F12.2)'),2)+')'
;   xstr = strtrim(string(stats.median,format='(F12.2)'),2)+'\pm'+strtrim(string(stats.sig68mean,format='(F12.2)'),2)
    
    atlas1d_lineplot, x, resid, xerr, resid_err, plottype=1, postscript=postscript, $
      xtitle='', ytitle=ytitle, xrange=xrange, yrange=yrange, legendtype=0, $
      charsize=charsize_7, xtickname=replicate(' ',10), position=pos[*,0], $
      xminor=3, yminor=3, atlaspsize=0.8, ytickinterval=1.0
    djs_oplot, !x.crange, [0,0], line=2, thick=postthick

; now overplot the AGN

    atlas1d_lineplot, xagn, residagn, xerragn, residagn_err, $
      postscript=postscript, atlaspsize=0.8, atlasfill=0, /overplot, atlascolor='red'

    legend, '(a)', /left, /top, box=0, charsize=charsize_3, $
      charthick=postthick, clear=keyword_set(postscript)
    legend, 'Cram et al. (1998)', /left, /bottom, box=0, charsize=charsize_3, $
      charthick=postthick, clear=keyword_set(postscript)
    legend, textoidl(xstr), /right, /top, box=0, charsize=charsize_3, charthick=postthick

; ##########################################################
; Panel 2
; ##########################################################

; SF galaxies
    
    indx = where((sfrs.sfr_ha gt -900) and (sfrs.sfr_uband_best gt -900),nindx)

    x = sfrs[indx].sfr_ha
    xerr = sfrs[indx].sfr_ha_err
    
    y1 = sfrs[indx].sfr_ha
    y1err = sfrs[indx].sfr_ha_err
    
    y2 = sfrs[indx].sfr_uband_best
    y2err = sfrs[indx].sfr_uband_best_err

    resid = y2-y1
    resid_err = sqrt(y1err^2+y2err^2)

; AGN
    
    indxagn = where((sfrs_agn.sfr_ha gt -900) and (sfrs_agn.sfr_uband_best gt -900),nindxagn)

    xagn = sfrs_agn[indxagn].sfr_ha
    xerragn = sfrs_agn[indxagn].sfr_ha_err
    
    y1agn = sfrs_agn[indxagn].sfr_ha
    y1erragn = sfrs_agn[indxagn].sfr_ha_err
    
    y2agn = sfrs_agn[indxagn].sfr_uband_best
    y2erragn = sfrs_agn[indxagn].sfr_uband_best_err

    residagn = y2agn-y1agn
    residagn_err = sqrt(y1erragn^2+y2erragn^2)

; make the plot    

    stats = im_stats(resid,/verbose,/baremin)
    xstr = strtrim(string(stats.median_rej,format='(F12.2)'),2)+' ('+$
      strtrim(string(stats.mean_rej,format='(F12.2)'),2)+$
      '\pm'+strtrim(string(stats.sigma_rej,format='(F12.2)'),2)+')'
;   xstr = strtrim(string(stats.median,format='(F12.2)'),2)+'\pm'+strtrim(string(stats.sig68mean,format='(F12.2)'),2)
    
    atlas1d_lineplot, x, resid, xerr, resid_err, plottype=1, postscript=postscript, $
      xtitle=xtitle, ytitle=ytitle, xrange=xrange, yrange=yrange, legendtype=0, $
      charsize=charsize_7, position=pos[*,1], /noerase, $
      xminor=3, yminor=3, atlaspsize=0.8, ytickinterval=1.0
    djs_oplot, !x.crange, [0,0], line=2, thick=postthick

; now overplot the AGN

    atlas1d_lineplot, xagn, residagn, xerragn, residagn_err, $
      postscript=postscript, atlaspsize=0.8, atlasfill=0, /overplot, atlascolor='red'

    legend, '(b)', /left, /top, box=0, charsize=charsize_3, $
      charthick=postthick, clear=keyword_set(postscript)
    legend, 'This Paper', /left, /bottom, box=0, charsize=charsize_3, $
      charthick=postthick, clear=keyword_set(postscript)
    legend, textoidl(xstr), /right, /top, box=0, charsize=charsize_3, charthick=postthick

    im_openclose, postscript=postscript, /close    

stop    

return
end    

