pro sfrs_talkplots, atlasdust, atlasnodust, nfgsdust, nfgsnodust, atlasdust_agn, $
  atlasnodust_agn, nfgsdust_agn, nfgsnodust_agn, sdssdust, sdssnodust, $
  sdssancillary, hii, postscript=postscript, cleanpng=cleanpng, $
  _extra=extra

; sfrs_talkplots,atlasdust,atlasnodust,nfgsdust,nfgsnodust,atlasdust_agn,atlasnodust_agn,nfgsdust_agn,nfgsnodust_agn,sdssdust,sdssnodust,sdssancillary,hii

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
    
    htmlbase = 'sfrs_talk'

    html_path = atlas_path(/web)+'analysis/'
    pspath = atlas_path(/web)+'analysis/'+htmlbase+'/'

    if keyword_set(postscript) then begin
       postthick = 8.0 
       postthick2 = 12.0
       postthick3 = 5.0
       postthick4 = 10.0
       symthick = 5.0
       highz_errthick = 4.0
       talkcolor = 'white'
       talk = 1L
    endif else begin
       blackwhite = 0L
       postthick = 2.0
       postthick2 = 2.0
       postthick3 = 2.0
       postthick4 = 2.0
       symthick = 2.0
       highz_errthick = 1.0
       im_window, 0, /square
       talkcolor = ''
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

    if keyword_set(blackwhite) then localcolor = 'medium gray' else localcolor = 'white'
    localsym = 108
    localpsize = 0.7
    
    if keyword_set(blackwhite) then localagncolor = 'medium gray' else localagncolor = 'white'
    localagnsym = 108
    localagnpsize = 1.0

    if keyword_set(blackwhite) then g99color = 'dark gray' else g99color = 'yellow'
    g99sym = 122
    g99psize = 1.8

    if keyword_set(blackwhite) then t02color = 'dark gray' else t02color = 'blue'
    t02sym = 104
    t02psize = 1.5

    if keyword_set(blackwhite) then h02color = 'dark gray' else h02color = 'green'
    h02sym = 106
    h02psize = 1.2

    if keyword_set(blackwhite) then liang04color = 'dark gray' else liang04color = 'yellow'
    liang04sym = 106
    liang04psize = 1.1

    if keyword_set(blackwhite) then lilly03color = 'dark gray' else lilly03color = 'green'
    lilly03sym = 122
    lilly03psize = 1.4

    if keyword_set(blackwhite) then k04color = 'dark gray' else k04color = 'white'
    k04sym = '104'
    k04psize = 0.8

    if keyword_set(blackwhite) then m05color = 'dark gray' else m05color = 'red'
    m05sym = 105
    m05psize = 1.4
    
    if keyword_set(blackwhite) then shap05color = 'dark gray' else shap05color = 'green'
    shap05sym = 122
    shap05psize = 1.9
    
    if keyword_set(blackwhite) then sava05color = 'dark gray' else sava05color = 'dodger blue'
    sava05sym = 115
    sava05psize = 2.0
    
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
; BEGIN TALK PLOTS
; ##################################################

; ------------------------------------------------------------
; wavelength versus redshift
; ------------------------------------------------------------

    psname = 'wavelength_vs_redshift'
    im_openclose, pspath+psname, postscript=postscript, xsize=8.5, ysize=8.5, encapsulated=encapsulated

    pagemaker, nx=1, ny=1, xspace=0, yspace=0.0, width=7.1, height=7.1, $
      xmargin=[1.1,0.3], ymargin=[0.3,1.1], xpage=8.5, ypage=8.5, $
      position=pos, /normal

    polyfill, [1,1,0,0,1], [1,0,0,1,1], /normal, color=djs_icolor('black')

    zmax = 5.2
    zmin = 0.0
    dz = 0.1
    zgrid = findgen((zmax-zmin)/dz+1)*dz+zmin

    lamgrid = findgen((26E4-3500.0)/1.0+1)*1.0+3500.0
    lamgrid = lamgrid / 1E4

;   lamrest = [3727.0,4861.0,6563.0] / 1E4 ; OII, H-beta, OIII, Ha
;   linename = ['[O II]','Hb','Ha']
;   linestyle = [0,3,5]

    lamrest = [3727.0,4340.0,4861.0,6563.0] / 1E4 ; OII, H-beta, OIII, Ha
    linename = ['[O II]','Hg','Hb','Ha']
    linestyle = [0,1,3,5]

;   linestyle = [0,3,5]
;   lamrest = [3727.0,4861.0,5007.0,6563.0] / 1E4 ; OII, H-beta, OIII, Ha
;   linename = ['[O II]','Hb','[O III]','Ha']
;   linestyle = [0,2,3,5]
    
;   zgrid = alog10(1+zgrid)
;   lamgrid = alog10(lamgrid)
;   lamrest = alog10(lamrest)

;   xtitle = 'log (\lambda) ['+angstrom()+']' 
    xtitle = textoidl('Observed Wavelength [\mu'+'m]')
;   xtitle = 'Wavelength ['+angstrom()+']'
    ytitle = ' Redshift' ; 'Redshift'
;   ytitle = 'log (1+z)' ; 'Redshift'

    xrange = [1600.0,26E3] / 1E4
;   xrange = [3000.0,30E3]
;   xrange = alog10([3000.0,30E3])
;   xrange = [min(lamrest)+min(zgrid),max(lamrest)+max(zgrid)]
    yrange = minmax(zgrid)
    
    plot, [0], [0], /nodata, xrange=xrange, yrange=yrange, /noerase, $
      charsize=2.0, charthick=postthick, xthick=postthick, ythick=postthick, $
      xsty=1, ysty=1, xtitle=xtitle, ytitle=ytitle, position=pos[*,0], color=djs_icolor(talkcolor);, /xlog, /ylog, $
;     ytickinterval=0.01, ytickname=['0','1','2','3','4','5'], yminor=50 ;, yticks=4
    
; fill the optical wavelength range
    
    polyfill, [3500.0,3500.0,9500,9500] / 1E4, [yrange[0],yrange[1],$
      yrange[1],yrange[0]], color=djs_icolor('orange'), /fill;/line_fill, $
;     spacing=0.1, orientation=135
    polyfill, [3500.0,3500.0,9500,9500] / 1E4, [yrange[0],yrange[1],$
      yrange[1],yrange[0]], color=djs_icolor('orange'), /fill;/line_fill, $
;     spacing=0.1, orientation=45
    xyouts, 6500.0/1E4, 4.5, 'Optical', charsize=2.0, charthick=postthick, align=0.5, /data, color=djs_icolor(talkcolor)

    polyfill, [11E3,11E3,14E3,14E3] / 1E4, [yrange[0],yrange[1],$
      yrange[1],yrange[0]], color=djs_icolor('light blue'), /fill;/line_fill, $
;     spacing=0.1, orientation=135
    polyfill, [11E3,11E3,14E3,14E3] / 1E4, [yrange[0],yrange[1],$
      yrange[1],yrange[0]], color=djs_icolor('light blue'), /fill;/line_fill, $
;     spacing=0.1, orientation=45
    xyouts, 12.5E3/1E4, 4.5, 'J', charsize=2.0, charthick=postthick, align=0.5, /data, color=djs_icolor(talkcolor)

    polyfill, [15E3,15E3,18E3,18E3] / 1E4, [yrange[0],yrange[1],$
      yrange[1],yrange[0]], color=djs_icolor('dark green'), /fill;/line_fill, $
;     spacing=0.1, orientation=135
    polyfill, [15E3,15E3,18E3,18E3] / 1E4, [yrange[0],yrange[1],$
      yrange[1],yrange[0]], color=djs_icolor('dark green'), /fill;/line_fill, $
;     spacing=0.1, orientation=45
    xyouts, 16.5E3/1E4, 4.5, 'H', charsize=2.0, charthick=postthick, align=0.5, /data, color=djs_icolor(talkcolor)

    polyfill, [20E3,20E3,23E3,23E3] / 1E4, [yrange[0],yrange[1],$
      yrange[1],yrange[0]], color=djs_icolor('red'), /fill;/line_fill, $
;     spacing=0.1, orientation=135
    polyfill, [20E3,20E3,23E3,23E3] / 1E4, [yrange[0],yrange[1],$
      yrange[1],yrange[0]], color=djs_icolor('red'), /fill;/line_fill, $
;     spacing=0.1, orientation=45
    xyouts, 21.5E3/1E4, 4.5, 'K', charsize=2.0, charthick=postthick, align=0.5, /data, color=djs_icolor(talkcolor)

    for i = 0L, n_elements(lamrest)-1L do djs_oplot, lamrest[i]*(1+zgrid), $
      zgrid, line=linestyle[i], thick=postthick, color=talkcolor

; overplot the tick marks    
    
    plot, [0], [0], /nodata, xrange=xrange, yrange=yrange, /noerase, $
      charsize=2.0, charthick=postthick, xthick=postthick, ythick=postthick, $
      xsty=1, ysty=1, xtitle='', ytitle='', position=pos[*,0], color=djs_icolor(talkcolor)
    
; legend

;   legend, textoidl(['[O II]','H\beta','H\alpha']), /right, $
    legend, textoidl(['[O II]','H\gamma','H\beta','H\alpha']), /right, $
;   legend, textoidl(['[O II]','H\beta','H\alpha']), /right, $
;   legend, textoidl(['[O II]','H\beta','[O III]','H\alpha']), /right, $
      /bottom, charsize=1.2, charthick=postthick, line=linestyle, thick=postthick, $
      clear=keyword_set(postscript), spacing=1.5, box=0

; when are each of the emission lines within the atmospheric windows?

    for iline = 0L, n_elements(linename)-1L do begin
       print, linename[iline]+':'
       print, ' Optical: <'+strtrim(string(interpol(zgrid,lamrest[iline]*(1+zgrid),0.95),format='(F12.1)'),2)
       print, ' J-band : '+strtrim(string(interpol(zgrid,lamrest[iline]*(1+zgrid),1.1),format='(F12.1)'),2)+'-'+$
         strtrim(string(interpol(zgrid,lamrest[iline]*(1+zgrid),1.4),format='(F12.1)'),2)
       print, ' H-band : '+strtrim(string(interpol(zgrid,lamrest[iline]*(1+zgrid),1.5),format='(F12.1)'),2)+'-'+$
         strtrim(string(interpol(zgrid,lamrest[iline]*(1+zgrid),1.8),format='(F12.1)'),2)
       print, ' K-band : '+strtrim(string(interpol(zgrid,lamrest[iline]*(1+zgrid),2.0),format='(F12.1)'),2)+'-'+$
         strtrim(string(interpol(zgrid,lamrest[iline]*(1+zgrid),2.3),format='(F12.1)'),2)
       print
    endfor
    
    im_openclose, postscript=postscript, /close    

; ------------------------------------------------------------
; L(B) versus SFR(Ha)/L([O II])_obs
; ------------------------------------------------------------
    
    psname = 'lb_vs_sfr_ha_loii_obs_nomedian'
    im_openclose, pspath+psname, postscript=postscript, xsize=8.5, ysize=8.8, encapsulated=encapsulated

    polyfill, [1,1,0,0,1], [1,0,0,1,1], /normal, color=djs_icolor('black')

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
       if keyword_set(talk) then medcolor = 'red' else medcolor = ''
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
      plottype=1, postscript=postscript, xtitle=xtitle, ytitle=ytitle, xrange=xrange, /noerase, $
      xstyle=11, yrange=yrange, legendtype=0, /right, /top, yfrac=1.5, xfrac=1.5, charsize=charsize_8, $
      atlaspsize=1.1, position=pos[*,0], atlascolor=notlirgscolor, blackwhite=blackwhite, talk=talk
    axis, /xaxis, xrange=interpol(xabs,x,!x.crange), xthick=postthick, xstyle=1, $
      charsize=charsize_8, charthick=postthick, color=djs_icolor(talkcolor)
    im_xyouts_title, xtitle='M_{B} [mag]', charsize=charsize_8, charthick=postthick, color=djs_icolor(talkcolor)

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

    im_symbols, 108, psize=1.3, fill=0, color=fsc_color(lirgscolor,!d.table_size-5), thick=postthick4
    djs_oplot, xbig_agn[lirgs_agn], ybig_agn[lirgs_agn], ps=8

    im_symbols, 105, psize=1.2, fill=0, color=fsc_color(nolirdatacolor,!d.table_size-6), thick=postthick4
    djs_oplot, xbig_agn[nolirdata_agn], ybig_agn[nolirdata_agn], ps=8

    label = ['L(IR) < 10^{11} L'+sunsymbol(),'L(IR) > 10^{11} L'+sunsymbol(),$
      'Undetected with IRAS']

; fit only to certain data

    xfit = xbig
    yfit = ybig
    
; 9.0 is M_B = -17.07 and 8.7 is M_B = -16.3; 8.57 is -16

    LBlocut = 10.0^6.0  ; 10^7.5
    LBhicut = 10.0^11.0 ; 10^10.8

    LBlo = alog10(LBlocut)>min(xfit)
    LBhi = alog10(LBhicut)<max(xfit)
    
; overlay the running median    

    running = im_medxbin(xfit,yfit,medbin,minx=minx,minpts=minpts,/verbose)

    doit = where((running.binctr gt LBlo) and (running.binctr le LBhi),ndoit)
    lolum = where((running.binctr lt LBlo))

;   im_symbols, 108, psize=3.0, thick=postthick, color=medcolor;, /fill
;   oploterror, running.binctr[doit], running.medy[doit], running.binsz/2.0+fltarr(ndoit), $
;     running.sigy75[doit]-running.medy[doit], ps=8, errthick=(postthick)>2L, /hibar, $
;     color=djs_icolor(medcolor), errcolor=djs_icolor(medcolor)
;   oploterror, running.binctr[doit], running.medy[doit], running.binsz/2.0+fltarr(ndoit), $
;     running.medy[doit]-running.sigy25[doit], ps=8, errthick=(postthick)>2L, /lobar, $
;     color=djs_icolor(medcolor), errcolor=djs_icolor(medcolor)

    im_openclose, postscript=postscript, /close    

; ------------------------------------------------------------
; L(B) versus SFR(Ha)/L([O II])_obs
; ------------------------------------------------------------
    
    psname = 'lb_vs_sfr_ha_loii_obs_withmedian'
    im_openclose, pspath+psname, postscript=postscript, xsize=8.5, ysize=8.8, encapsulated=encapsulated

    polyfill, [1,1,0,0,1], [1,0,0,1,1], /normal, color=djs_icolor('black')

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
       if keyword_set(talk) then medcolor = 'red' else medcolor = ''
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
      plottype=1, postscript=postscript, xtitle=xtitle, ytitle=ytitle, xrange=xrange, /noerase, $
      xstyle=11, yrange=yrange, legendtype=0, /right, /top, yfrac=1.5, xfrac=1.5, charsize=charsize_8, $
      atlaspsize=1.1, position=pos[*,0], atlascolor=notlirgscolor, blackwhite=blackwhite, talk=talk
    axis, /xaxis, xrange=interpol(xabs,x,!x.crange), xthick=postthick, xstyle=1, $
      charsize=charsize_8, charthick=postthick, color=djs_icolor(talkcolor)
    im_xyouts_title, xtitle='M_{B} [mag]', charsize=charsize_8, charthick=postthick, color=djs_icolor(talkcolor)

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

    im_symbols, 108, psize=1.3, fill=0, color=fsc_color(lirgscolor,!d.table_size-5), thick=postthick4
    djs_oplot, xbig_agn[lirgs_agn], ybig_agn[lirgs_agn], ps=8

    im_symbols, 105, psize=1.2, fill=0, color=fsc_color(nolirdatacolor,!d.table_size-6), thick=postthick4
    djs_oplot, xbig_agn[nolirdata_agn], ybig_agn[nolirdata_agn], ps=8

    label = ['L(IR) < 10^{11} L'+sunsymbol(),'L(IR) > 10^{11} L'+sunsymbol(),$
      'Undetected with IRAS']

; fit only to certain data

    xfit = xbig
    yfit = ybig
    
; 9.0 is M_B = -17.07 and 8.7 is M_B = -16.3; 8.57 is -16

    LBlocut = 10.0^6.0  ; 10^7.5
    LBhicut = 10.0^11.0 ; 10^10.8

    LBlo = alog10(LBlocut)>min(xfit)
    LBhi = alog10(LBhicut)<max(xfit)
    
; overlay the running median    

    running = im_medxbin(xfit,yfit,medbin,minx=minx,minpts=minpts,/verbose)

    doit = where((running.binctr gt LBlo) and (running.binctr le LBhi),ndoit)
    lolum = where((running.binctr lt LBlo))

    im_symbols, 108, psize=3.0, thick=postthick, color=medcolor;, /fill
    oploterror, running.binctr[doit], running.medy[doit], running.binsz/2.0+fltarr(ndoit), $
      running.sigy75[doit]-running.medy[doit], ps=8, errthick=(postthick)>2L, /hibar, $
      color=djs_icolor(medcolor), errcolor=djs_icolor(medcolor)
    oploterror, running.binctr[doit], running.medy[doit], running.binsz/2.0+fltarr(ndoit), $
      running.medy[doit]-running.sigy25[doit], ps=8, errthick=(postthick)>2L, /lobar, $
      color=djs_icolor(medcolor), errcolor=djs_icolor(medcolor)
    
; overplot the K98 transformation assuming one magnitude of extinction
; at H-alpha

;   djs_oplot, !x.crange, alog10(1.4*10^(0.4*1.0))*[1,1], line=2, thick=postthick4, color=djs_icolor(talkcolor)
    
    im_openclose, postscript=postscript, /close    

; ------------------------------------------------------------
; L(B) versus SFR(Ha)/L([O II])_obs
; ------------------------------------------------------------
    
    psname = 'lb_vs_sfr_ha_loii_obs_withmedian_withk98'
    im_openclose, pspath+psname, postscript=postscript, xsize=8.5, ysize=8.8, encapsulated=encapsulated

    polyfill, [1,1,0,0,1], [1,0,0,1,1], /normal, color=djs_icolor('black')

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
       if keyword_set(talk) then medcolor = 'red' else medcolor = ''
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
      plottype=1, postscript=postscript, xtitle=xtitle, ytitle=ytitle, xrange=xrange, /noerase, $
      xstyle=11, yrange=yrange, legendtype=0, /right, /top, yfrac=1.5, xfrac=1.5, charsize=charsize_8, $
      atlaspsize=1.1, position=pos[*,0], atlascolor=notlirgscolor, blackwhite=blackwhite, talk=talk
    axis, /xaxis, xrange=interpol(xabs,x,!x.crange), xthick=postthick, xstyle=1, $
      charsize=charsize_8, charthick=postthick, color=djs_icolor(talkcolor)
    im_xyouts_title, xtitle='M_{B} [mag]', charsize=charsize_8, charthick=postthick, color=djs_icolor(talkcolor)

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

    im_symbols, 108, psize=1.3, fill=0, color=fsc_color(lirgscolor,!d.table_size-5), thick=postthick4
    djs_oplot, xbig_agn[lirgs_agn], ybig_agn[lirgs_agn], ps=8

    im_symbols, 105, psize=1.2, fill=0, color=fsc_color(nolirdatacolor,!d.table_size-6), thick=postthick4
    djs_oplot, xbig_agn[nolirdata_agn], ybig_agn[nolirdata_agn], ps=8

    label = ['L(IR) < 10^{11} L'+sunsymbol(),'L(IR) > 10^{11} L'+sunsymbol(),$
      'Undetected with IRAS']

; fit only to certain data

    xfit = xbig
    yfit = ybig
    
; 9.0 is M_B = -17.07 and 8.7 is M_B = -16.3; 8.57 is -16

    LBlocut = 10.0^6.0  ; 10^7.5
    LBhicut = 10.0^11.0 ; 10^10.8

    LBlo = alog10(LBlocut)>min(xfit)
    LBhi = alog10(LBhicut)<max(xfit)
    
; overlay the running median    

    running = im_medxbin(xfit,yfit,medbin,minx=minx,minpts=minpts,/verbose)

    doit = where((running.binctr gt LBlo) and (running.binctr le LBhi),ndoit)
    lolum = where((running.binctr lt LBlo))

    im_symbols, 108, psize=3.0, thick=postthick, color=medcolor;, /fill
    oploterror, running.binctr[doit], running.medy[doit], running.binsz/2.0+fltarr(ndoit), $
      running.sigy75[doit]-running.medy[doit], ps=8, errthick=(postthick)>2L, /hibar, $
      color=djs_icolor(medcolor), errcolor=djs_icolor(medcolor)
    oploterror, running.binctr[doit], running.medy[doit], running.binsz/2.0+fltarr(ndoit), $
      running.medy[doit]-running.sigy25[doit], ps=8, errthick=(postthick)>2L, /lobar, $
      color=djs_icolor(medcolor), errcolor=djs_icolor(medcolor)
    
; overplot the K98 transformation assuming one magnitude of extinction
; at H-alpha

    djs_oplot, !x.crange, alog10(1.4*10^(0.4*1.0))*[1,1], line=2, thick=postthick4, color=djs_icolor(talkcolor)
    
    im_openclose, postscript=postscript, /close    

; ------------------------------------------------------------
; 2-panel - SFR(IR) vs SFR(Ha)/SFR(IR)
; ------------------------------------------------------------
    
    psname = 'sfr_ir_vs_sfr_hair_haobs'
    im_openclose, pspath+psname, postscript=postscript, xsize=8.5, ysize=9.0, encapsulated=encapsulated

    polyfill, [1,1,0,0,1], [1,0,0,1,1], /normal, color=djs_icolor('black')

    pagemaker, nx=1, ny=2, xspace=0, yspace=0.0, width=6.5, height=[5.0,2.5], $
      xmargin=[1.5,0.5], ymargin=[0.4,1.1], xpage=8.5, ypage=9.0, $
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
    yrange = xrange
    yrange2 = sfrHaIRrange2

    xtitle = 'log \psi(IR) [M'+sunsymbol()+' yr^{-1}]'
    ytitle = 'log \psi(H\alpha) [M'+sunsymbol()+' yr^{-1}]'
;   ytitle2 = 'log \Delta\psi' 
    ytitle2 = 'log [\psi(H\alpha)/\psi(IR)]'

    lhalir = alog10(4.5D-44/7.9D-42)

; ##########################################################
; Panel 1: No reddening correction
; ##########################################################

    ha = atlasdust[indx].h_alpha[0]
    ha_err = atlasdust[indx].h_alpha[1]
    
    y = alog10(ha/lir) - lhalir
    yerr = im_compute_error(ha,ha_err,lir,lir_err,/log)

    sfr_ha = alog10(haconst*ha*4.0*!dpi*(atlasdust[indx].distance*3.826D24)^2)
    sfr_ir = alog10(irconst*lir*4.0*!dpi*(atlasdust[indx].distance*3.826D24)^2)
    
    hanfgs = nfgsdust[indxnfgs].h_alpha[0]
    hanfgs_err = nfgsdust[indxnfgs].h_alpha[1]
    
    ynfgs = alog10(hanfgs/lirnfgs) - lhalir
    yerrnfgs = im_compute_error(hanfgs,hanfgs_err,lirnfgs,lirnfgs_err,/log)

    sfr_hanfgs = alog10(haconst*hanfgs*4.0*!dpi*(nfgsdust[indxnfgs].distance*3.826D24)^2)
    sfr_irnfgs = alog10(irconst*lirnfgs*4.0*!dpi*(nfgsdust[indxnfgs].distance*3.826D24)^2)
    
    stats = im_stats([y,ynfgs],/verbose)
    xstr = strtrim(string(stats.median_rej,format='(F12.2)'),2)+' ('+$
      strtrim(string(stats.mean_rej,format='(F12.2)'),2)+$
      '\pm'+strtrim(string(stats.sigma_rej,format='(F12.2)'),2)+')'

    atlas1d_lineplot, sfr_ir, sfr_ha, xerr, yerr, plottype=1, postscript=postscript, $
      xtitle='', ytitle=ytitle, xrange=xrange, yrange=yrange, legendtype=0, $
      charsize=singlecharsize, position=pos[*,0], xnfgs=sfr_irnfgs, ynfgs=sfr_hanfgs, $
      xerrnfgs=xerrnfgs, yerrnfgs=yerrnfgs, xtickname=replicate(' ',10), $
      yfrac=8, blackwhite=blackwhite, talk=talk
    djs_oplot, !x.crange, !y.crange, line=0, thick=postthick4, color=djs_icolor(talkcolor)

    legend, textoidl('Observed H\alpha'), /left, /top, box=0, $
      charsize=singlecharsize, charthick=postthick, textcolor=djs_icolor(talkcolor)

    atlas1d_lineplot, x, y, xerr, yerr, plottype=1, postscript=postscript, $
      xtitle=xtitle, ytitle=ytitle2, xrange=xrange, yrange=yrange2, legendtype=0, $
      charsize=singlecharsize, position=pos[*,1], xnfgs=xnfgs, ynfgs=ynfgs, $
      xerrnfgs=xerrnfgs, yerrnfgs=yerrnfgs, /noerase, $
      yfrac=8, blackwhite=blackwhite, talk=talk
    djs_oplot, !x.crange, [0,0], line=0, thick=postthick4, color=djs_icolor(talkcolor)

    im_openclose, postscript=postscript, /close
    
; ############################################################
; Panel 2: Good Reddening Correction
; ############################################################

    psname = 'sfr_ir_vs_sfr_hair_hacor'
    im_openclose, pspath+psname, postscript=postscript, xsize=8.5, ysize=9.0, encapsulated=encapsulated

    polyfill, [1,1,0,0,1], [1,0,0,1,1], /normal, color=djs_icolor('black')

    pagemaker, nx=1, ny=2, xspace=0, yspace=0.0, width=6.5, height=[5.0,2.5], $
      xmargin=[1.5,0.5], ymargin=[0.4,1.1], xpage=8.5, ypage=9.0, $
      position=pos, /normal

    ha = atlasnodust[indx].h_alpha[0]
    ha_err = atlasnodust[indx].h_alpha[1]

    y = alog10(ha/lir) - lhalir
    yerr = im_compute_error(ha,ha_err,lir,lir_err,/log)

    sfr_ha = alog10(haconst*ha*4.0*!dpi*(atlasdust[indx].distance*3.826D24)^2)
    sfr_ir = alog10(irconst*lir*4.0*!dpi*(atlasdust[indx].distance*3.826D24)^2)
    
    hanfgs = nfgsnodust[indxnfgs].h_alpha[0]
    hanfgs_err = nfgsnodust[indxnfgs].h_alpha[1]

    ynfgs = alog10(hanfgs/lirnfgs) - lhalir
    yerrnfgs = im_compute_error(hanfgs,hanfgs_err,lirnfgs,lirnfgs_err,/log)

    sfr_hanfgs = alog10(haconst*hanfgs*4.0*!dpi*(nfgsdust[indxnfgs].distance*3.826D24)^2)
    sfr_irnfgs = alog10(irconst*lirnfgs*4.0*!dpi*(nfgsdust[indxnfgs].distance*3.826D24)^2)
    
    stats = im_stats([y,ynfgs],/verbose,/no_head)
    xstr = strtrim(string(stats.median_rej,format='(F12.2)'),2)+' ('+$
      strtrim(string(stats.mean_rej,format='(F12.2)'),2)+$
      '\pm'+strtrim(string(stats.sigma_rej,format='(F12.2)'),2)+')'

    atlas1d_lineplot, sfr_ir, sfr_ha, xerr, yerr, plottype=1, postscript=postscript, $
      xtitle='', ytitle=ytitle, xrange=xrange, yrange=yrange, legendtype=0, $
      charsize=singlecharsize, position=pos[*,0], xnfgs=sfr_irnfgs, ynfgs=sfr_hanfgs, $
      xerrnfgs=xerrnfgs, yerrnfgs=yerrnfgs, xtickname=replicate(' ',10), $
      yfrac=8, blackwhite=blackwhite, talk=talk
    djs_oplot, !x.crange, !y.crange, line=0, thick=postthick4, color=djs_icolor(talkcolor)

    legend, textoidl('Extinction-Corrected H\alpha'), /left, /top, box=0, $
      charsize=singlecharsize, charthick=postthick, textcolor=djs_icolor(talkcolor)

    atlas1d_lineplot, x, y, xerr, yerr, plottype=1, postscript=postscript, $
      xtitle=xtitle, ytitle=ytitle2, xrange=xrange, yrange=yrange2, legendtype=0, $
      charsize=singlecharsize, position=pos[*,1], xnfgs=xnfgs, ynfgs=ynfgs, $
      xerrnfgs=xerrnfgs, yerrnfgs=yerrnfgs, /noerase, $
      yfrac=8, blackwhite=blackwhite, talk=talk
    djs_oplot, !x.crange, [0,0], line=0, thick=postthick4, color=djs_icolor(talkcolor)

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

    atlas1d_lineplot, x, y, xerr, yerr, plottype=1, postscript=postscript, $
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
    
    sdss_lineplot, x, y, xerr, yerr, plottype=1, postscript=postscript, $
      xtitle=xtitle, ytitle='', xrange=xrange, yrange=yrange, legendtype=0, $
      /left, /top, position=pos[*,1], $
      ytickname=replicate(' ',10), /noerase, charsize=charsize_6
    oplot, xaxis, yline, linestyle=2, thick=postthick

    im_openclose, postscript=postscript, /close    

; ------------------------------------------------------------
; L(B) versus SFR(Ha)/L(Hb)_obs - Integrated
; ------------------------------------------------------------

    psname = 'lb_vs_sfr_ha_lhb_obs_nomedian'
    im_openclose, pspath+psname, postscript=postscript, xsize=8.5, ysize=8.8, encapsulated=encapsulated

    polyfill, [1,1,0,0,1], [1,0,0,1,1], /normal, color=djs_icolor('black')

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
       if keyword_set(talk) then medcolor = 'red' else medcolor = ''
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
      plottype=1, postscript=postscript, xtitle=xtitle, ytitle=ytitle, xrange=xrange, /noerase, $
      xstyle=11, yrange=yrange, legendtype=0, /right, /top, yfrac=1.5, xfrac=1.5, charsize=charsize_8, $
      atlaspsize=1.1, position=pos[*,0], atlascolor=notlirgscolor, blackwhite=blackwhite, talk=talk
    axis, /xaxis, xrange=interpol(xabs,x,!x.crange), xthick=postthick, xstyle=1, $
      charsize=charsize_8, charthick=postthick, color=djs_icolor(talkcolor)
    im_xyouts_title, xtitle='M_{B} [mag]', charsize=charsize_8, charthick=postthick, color=djs_icolor(talkcolor)

    im_symbols, 108, psize=1.4, /fill, color=fsc_color(lirgscolor,!d.table_size-5)
    djs_oplot, xbig[lirgs], ybig[lirgs], ps=8

    im_symbols, 105, psize=1.4, /fill, color=fsc_color(nolirdatacolor,!d.table_size-6)
    djs_oplot, xbig[nolirdata], ybig[nolirdata], ps=8

; AGN galaxies    
    
    atlas1d_lineplot, xbig_agn[notlirgs_agn], ybig_agn[notlirgs_agn], $
      xerrbig_agn[notlirgs_agn], yerrbig_agn[notlirgs_agn], blackwhite=blackwhite, talk=talk, $
      plottype=1, postscript=postscript, xtitle=xtitle, ytitle=ytitle, xrange=xrange, $
      yrange=yrange, legendtype=0, /right, /top, yfrac=1.5, xfrac=1.5, charsize=charsize_8, $
      atlaspsize=1.0, position=pos[*,0], atlasfill=0, /overplot, thick=postthick, atlascolor=notlirgscolor

    im_symbols, 108, psize=1.3, fill=0, color=fsc_color(lirgscolor,!d.table_size-5), thick=postthick4
    djs_oplot, xbig_agn[lirgs_agn], ybig_agn[lirgs_agn], ps=8

    im_symbols, 105, psize=1.2, fill=0, color=fsc_color(nolirdatacolor,!d.table_size-6), thick=postthick4
    djs_oplot, xbig_agn[nolirdata_agn], ybig_agn[nolirdata_agn], ps=8

    label = ['L(IR) < 10^{11} L'+sunsymbol(),'L(IR) > 10^{11} L'+sunsymbol(),$
      'Undetected with IRAS']

; fit only to certain data

    xfit = xbig
    yfit = ybig
    
; 9.0 is M_B = -17.07 and 8.7 is M_B = -16.3; 8.57 is -16

    LBlocut = 10.0^6.0 ; 10.0^8.7
    LBhicut = 10.0^11.0

    LBlo = alog10(LBlocut)>min(xfit)
    LBhi = alog10(LBhicut)<max(xfit)

; overlay the running median    

    running = im_medxbin(xfit,yfit,medbin,minx=minx,minpts=minpts,/verbose)

    doit = where((running.binctr gt LBlo) and (running.binctr le LBhi),ndoit)
    lolum = where((running.binctr lt LBlo))

;   im_symbols, 108, psize=3.0, thick=postthick, color=medcolor;, /fill
;   oploterror, running.binctr[doit], running.medy[doit], running.binsz/2.0+fltarr(ndoit), $
;     running.sigy75[doit]-running.medy[doit], ps=8, errthick=(postthick)>2L, /hibar, $
;     color=djs_icolor(medcolor), errcolor=djs_icolor(medcolor)
;   oploterror, running.binctr[doit], running.medy[doit], running.binsz/2.0+fltarr(ndoit), $
;     running.medy[doit]-running.sigy25[doit], ps=8, errthick=(postthick)>2L, /lobar, $
;     color=djs_icolor(medcolor), errcolor=djs_icolor(medcolor)

; now fit the high-luminosity end; constrain the intercept
    
    xfitaxis = findgen((LBhi-LBlo)/0.05+1L)*0.05 + LBlo
    yfitaxis = interpol(running.medy[doit],running.binctr[doit],xfitaxis)
    yerraxis = sqrt(interpol(running.sigy[doit]^2,running.binctr[doit],xfitaxis))

    im_openclose, postscript=postscript, /close    

; ------------------------------------------------------------
; L(B) versus SFR(Ha)/L(Hb)_obs - Integrated
; ------------------------------------------------------------

    psname = 'lb_vs_sfr_ha_lhb_obs_withmedian'
    im_openclose, pspath+psname, postscript=postscript, xsize=8.5, ysize=8.8, encapsulated=encapsulated

    polyfill, [1,1,0,0,1], [1,0,0,1,1], /normal, color=djs_icolor('black')

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
       if keyword_set(talk) then medcolor = 'red' else medcolor = ''
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
      plottype=1, postscript=postscript, xtitle=xtitle, ytitle=ytitle, xrange=xrange, /noerase, $
      xstyle=11, yrange=yrange, legendtype=0, /right, /top, yfrac=1.5, xfrac=1.5, charsize=charsize_8, $
      atlaspsize=1.1, position=pos[*,0], atlascolor=notlirgscolor, blackwhite=blackwhite, talk=talk
    axis, /xaxis, xrange=interpol(xabs,x,!x.crange), xthick=postthick, xstyle=1, $
      charsize=charsize_8, charthick=postthick, color=djs_icolor(talkcolor)
    im_xyouts_title, xtitle='M_{B} [mag]', charsize=charsize_8, charthick=postthick, color=djs_icolor(talkcolor)

    im_symbols, 108, psize=1.4, /fill, color=fsc_color(lirgscolor,!d.table_size-5)
    djs_oplot, xbig[lirgs], ybig[lirgs], ps=8

    im_symbols, 105, psize=1.4, /fill, color=fsc_color(nolirdatacolor,!d.table_size-6)
    djs_oplot, xbig[nolirdata], ybig[nolirdata], ps=8

; AGN galaxies    
    
    atlas1d_lineplot, xbig_agn[notlirgs_agn], ybig_agn[notlirgs_agn], $
      xerrbig_agn[notlirgs_agn], yerrbig_agn[notlirgs_agn], blackwhite=blackwhite, talk=talk, $
      plottype=1, postscript=postscript, xtitle=xtitle, ytitle=ytitle, xrange=xrange, $
      yrange=yrange, legendtype=0, /right, /top, yfrac=1.5, xfrac=1.5, charsize=charsize_8, $
      atlaspsize=1.0, position=pos[*,0], atlasfill=0, /overplot, thick=postthick, atlascolor=notlirgscolor

    im_symbols, 108, psize=1.3, fill=0, color=fsc_color(lirgscolor,!d.table_size-5), thick=postthick4
    djs_oplot, xbig_agn[lirgs_agn], ybig_agn[lirgs_agn], ps=8

    im_symbols, 105, psize=1.2, fill=0, color=fsc_color(nolirdatacolor,!d.table_size-6), thick=postthick4
    djs_oplot, xbig_agn[nolirdata_agn], ybig_agn[nolirdata_agn], ps=8

    label = ['L(IR) < 10^{11} L'+sunsymbol(),'L(IR) > 10^{11} L'+sunsymbol(),$
      'Undetected with IRAS']

; fit only to certain data

    xfit = xbig
    yfit = ybig
    
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
      running.sigy75[doit]-running.medy[doit], ps=8, errthick=(postthick)>2L, /hibar, $
      color=djs_icolor(medcolor), errcolor=djs_icolor(medcolor)
    oploterror, running.binctr[doit], running.medy[doit], running.binsz/2.0+fltarr(ndoit), $
      running.medy[doit]-running.sigy25[doit], ps=8, errthick=(postthick)>2L, /lobar, $
      color=djs_icolor(medcolor), errcolor=djs_icolor(medcolor)

; now fit the high-luminosity end; constrain the intercept
    
    xfitaxis = findgen((LBhi-LBlo)/0.05+1L)*0.05 + LBlo
    yfitaxis = interpol(running.medy[doit],running.binctr[doit],xfitaxis)
    yerraxis = sqrt(interpol(running.sigy[doit]^2,running.binctr[doit],xfitaxis))

    im_openclose, postscript=postscript, /close    

; ------------------------------------------------------------
; High-Redshift L(B) vs [O III]/[O II]
; ------------------------------------------------------------

    psname = 'highz_lb_vs_oiiioii'
    im_openclose, pspath+psname, postscript=postscript, xsize=8.5, ysize=8.0, encapsulated=encapsulated

    polyfill, [1,1,0,0,1], [1,0,0,1,1], /normal, color=djs_icolor('black')

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
    ytitle = textoidl('log ([O III] \lambda5007/[O II])_{obs}')

    xrange = LBrange
    yrange = oiiioiirange2

    plot, [0], [0], /nodata, xthick=postthick, ythick=postthick, charsize=charsize_8, $
      charthick=postthick, thick=postthick, xtitle=xtitle, ytitle=ytitle, /noerase, color=djs_icolor(talkcolor), $
      ymargin=[4,3], xstyle=11, ystyle=11, xrange=xrange, yrange=yrange, position=pos[*,0]
    axis, /yaxis, yrange=interpol(U,sbgrids[3,*].oiii_5007_oii,!y.crange), ythick=postthick, ystyle=1, $
      charsize=charsize_8, charthick=postthick, color=djs_icolor(talkcolor)
    axis, /xaxis, xrange=interpol(xabs,x,!x.crange), xthick=postthick, xstyle=1, $
      charsize=charsize_8, charthick=postthick, color=djs_icolor(talkcolor)
    im_xyouts_title, xtitle='M_{B} [mag]', charsize=charsize_8, charthick=postthick, color=djs_icolor(talkcolor)
    im_xyouts_title, ytitle='log U', charsize=charsize_8, charthick=postthick, xspacing=11.0, color=djs_icolor(talkcolor)

    im_symbols, localsym, psize=localpsize, /fill, color=fsc_color(localcolor,!d.table_size-151)
    djs_oplot, xbig, ybig, ps=8, color=fsc_color(localcolor,!d.table_size-151)

    im_symbols, localagnsym, psize=localagnpsize, fill=0, color=fsc_color(localagncolor,!d.table_size-151), thick=symthick
    djs_oplot, xbig_agn, ybig_agn, ps=8, color=fsc_color(localagncolor,!d.table_size-151), thick=symthick

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
    psym = [localsym,localagnsym,m05sym,sava05sym,liang04sym,lilly03sym]
    fill = [1,0,1,1,1,1]
    color = fsc_color([localcolor,localagncolor,m05color,sava05color,liang04color,lilly03color],!d.table_size-[151,180,152,153,154,155])

    postthick1 = postthick
    im_legend, textoidl(label), /left, /bottom, box=0, charsize=charsize_2, $
      charthick=postthick1, psym=psym, fill=fill, symsize=1.3, color=color, $
      spacing=1.8, thick=postthick1, textcolor=djs_icolor(replicate(talkcolor,6))
    
    im_openclose, postscript=postscript, /close    

; ------------------------------------------------------------
; High-Redshift L(B) vs [O II]/Ha
; ------------------------------------------------------------

    psname = 'highz_lb_vs_oiiha'
    im_openclose, pspath+psname, postscript=postscript, xsize=8.5, ysize=8.8, encapsulated=encapsulated

    polyfill, [1,1,0,0,1], [1,0,0,1,1], /normal, color=djs_icolor('black')

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
    ytitle = textoidl('log ([O II]/H\alpha)_{obs}')

    xrange = LBrange
    yrange = [-1.6,1.0] ; oiihacorrange

    plot, [0], [0], /nodata, xthick=postthick, ythick=postthick, charsize=charsize_8, $
      charthick=postthick, thick=postthick, xtitle=xtitle, ytitle=ytitle, color=djs_icolor(talkcolor), $
      ymargin=[4,3], xstyle=11, ystyle=3, xrange=xrange, yrange=yrange, position=pos[*,0], /noerase
    axis, /xaxis, xrange=interpol(xabs,x,!x.crange), xthick=postthick, xstyle=1, $
      charsize=charsize_8, charthick=postthick, color=djs_icolor(talkcolor)
    im_xyouts_title, xtitle='M_{B} [mag]', charsize=charsize_8, charthick=postthick, color=djs_icolor(talkcolor)

    im_symbols, localsym, psize=localpsize, /fill, color=fsc_color(localcolor,!d.table_size-101)
    djs_oplot, xbig, ybig, ps=8, color=fsc_color(localcolor,!d.table_size-101)

    im_symbols, localagnsym, psize=localagnpsize, fill=0, color=fsc_color(localagncolor,!d.table_size-101), thick=symthick
    djs_oplot, xbig_agn, ybig_agn, ps=8, color=fsc_color(localagncolor,!d.table_size-101)

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
    psym = [localsym,localagnsym,m05sym,t02sym,h02sym,g99sym]
    fill = [1,0,1,1,1,1]
    color = fsc_color([localcolor,localagncolor,m05color,t02color,h02color,g99color],!d.table_size-[101,101,102,103,104,105])
    postthick1 = postthick
    im_legend, textoidl(label), /left, /bottom, box=0, charsize=charsize_2, $
      charthick=postthick1, psym=psym, fill=fill, color=color, $
      spacing=1.8, thick=postthick1, textcolor=djs_icolor(replicate(talkcolor,6))

    im_openclose, postscript=postscript, /close    

; ------------------------------------------------------------
; High-Redshift L(B) vs Hb/Ha
; ------------------------------------------------------------

    psname = 'highz_lb_vs_hahb'
    im_openclose, pspath+psname, postscript=postscript, xsize=8.5, ysize=8.3, encapsulated=encapsulated

    polyfill, [1,1,0,0,1], [1,0,0,1,1], /normal, color=djs_icolor('black')

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
    ytitle = textoidl('log (H\alpha/H\beta)_{obs}')

    xrange = LBrange
    yrange = hahbrange

    plot, [0], [0], /nodata, xthick=postthick, ythick=postthick, charsize=charsize_8, $
      charthick=postthick, thick=postthick, xtitle=xtitle, ytitle=ytitle, color=djs_icolor(talkcolor), $
      ymargin=[4,3], xstyle=11, ystyle=11, xrange=xrange, yrange=yrange, position=pos[*,0], /noerase
    axis, /xaxis, xrange=interpol(xabs,x,!x.crange), xthick=postthick, xstyle=1, $
      charsize=charsize_8, charthick=postthick, color=djs_icolor(talkcolor)
    axis, /yaxis, yrange=interpol(yebv,y,!y.crange), ythick=postthick, ystyle=1, $
      charsize=charsize_8, charthick=postthick, color=djs_icolor(talkcolor)
    im_xyouts_title, xtitle='M_{B} [mag]', charsize=charsize_8, charthick=postthick, color=djs_icolor(talkcolor)
    im_xyouts_title, ytitle='E(B-V) [mag]', charsize=charsize_8, charthick=postthick, xspacing=9.0, color=djs_icolor(talkcolor)

    im_symbols, localsym, psize=localpsize, /fill, color=fsc_color(localcolor,!d.table_size-75)
    djs_oplot, xbig, ybig, ps=8, color=fsc_color(localcolor,!d.table_size-75)

    im_symbols, localagnsym, psize=localagnpsize, fill=0, color=fsc_color(localagncolor,!d.table_size-75), thick=symthick
    djs_oplot, xbig_agn, ybig_agn, ps=8, color=fsc_color(localagncolor,!d.table_size-75)

; put liang on bottom because darkest

    im_symbols, liang04sym, psize=liang04psize, /fill, color=fsc_color(liang04color,!d.table_size-79), thick=postthick
    djs_oplot, xliang04, yliang04, ps=8, color=fsc_color(liang04color,!d.table_size-79), thick=postthick

    im_symbols, m05sym, psize=m05psize, /fill, color=fsc_color(m05color,!d.table_size-76), thick=postthick
    djs_oplot, xm05, ym05, ps=8, color=fsc_color(m05color,!d.table_size-76), thick=postthick

    im_symbols, shap05sym, psize=shap05psize, /fill, color=fsc_color(shap05color,!d.table_size-77), thick=postthick
    djs_oplot, xshap05, yshap05, ps=8, color=fsc_color(shap05color,!d.table_size-77), thick=postthick

    im_symbols, sava05sym, psize=sava05psize, /fill, color=fsc_color(sava05color,!d.table_size-78), thick=postthick
    djs_oplot, xsava05, ysava05, ps=8, color=fsc_color(sava05color,!d.table_size-78), thick=postthick

    djs_oplot, !x.crange, alog10(HaHb)*[1,1], line=0, thick=postthick, color=djs_icolor(talkcolor)

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
    psym = [localsym,localagnsym,m05sym,shap05sym,sava05sym,liang04sym]
    fill = [1,0,1,1,1,1]
    color = fsc_color([localcolor,localagncolor,m05color,shap05color,sava05color,liang04color],!d.table_size-[75,75,76,77,78,79])
    postthick1 = postthick
    im_legend, textoidl(label), /left, /bottom, box=0, charsize=charsize_2, $
      charthick=postthick1, psym=psym, fill=fill, symsize=1.3, color=color, $
      spacing=1.8, thick=postthick1, textcolor=djs_icolor(replicate(talkcolor,6))

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
    
    atlas1d_lineplot, x, y, xerr, yerr, plottype=1, postscript=postscript, ysty=3, $
      ytitle=ytitle, xtitle=xtitle, xrange=xrange, yrange=yrange, legendtype=0, $
      xstyle=11, /right, /top, xnfgs=xnfgs, ynfgs=ynfgs, xerrnfgs=xerrnfgs, blackwhite=blackwhite, $
      yerrnfgs=yerrnfgs, position=pos[*,0], charsize=charsize_9, atlaspsize=1.1, nfgspsize=1.4
    axis, /xaxis, xrange=interpol(xabs,x,!x.crange), xthick=postthick, xstyle=1, $
      charsize=charsize_9, charthick=postthick
    im_xyouts_title, xtitle='M_{B} [mag]', charsize=charsize_9, charthick=postthick

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
    
    atlas1d_lineplot, x, y, xerr, yerr, plottype=1, postscript=postscript, ysty=3, $
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

    sdss_lineplot, x, y, xerr, yerr, plottype=1, postscript=postscript, $
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

    atlas1d_lineplot, x, y, xerr, yerr, plottype=1, postscript=postscript, $
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

    sdss_lineplot, x, y, xerr, yerr, plottype=1, postscript=postscript, $
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
; 12+log(O/H) vs ([O II]/Ha)_cor [Integrated + SDSS]
; ------------------------------------------------------------

    psname = '12oh_oiiinii_niiha_vs_oiiha_2panel_nogrids_nolines'
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

    atlas1d_lineplot, x, y, xerr, yerr, plottype=1, postscript=postscript, $n
      xtitle=xtitle, ytitle=ytitle, xrange=xrange, yrange=yrange, legendtype=0, $
      /left, /top, position=pos[*,0], charsize=charsize_5, yfrac=1.5, /errorleft, $
      xregion=xregion, yregion=yregion, xerrregion=xerrregion, yerrregion=yerrregion, $
      xnfgs=xnfgs, ynfgs=ynfgs, xerrnfgs=xerrnfgs, yerrnfgs=yerrnfgs, blackwhite=blackwhite

    if keyword_set(postscript) then thisthick = 6.0 else thickthick = 2.0

;   plot_kewley_grids, plotnumber=17, model=3, labeltype=4, /log12oh, /overplot, $
;     /noZgrid, thick=3.0, Ugridcolor=ugridcolor, Ulinestyle=0, $
;     postscript=postscript, Ulabelcolor='', Ucharsize=0.8, charthick=thisthick
;   xyouts, 9.1, -0.18, 'log U', align=0.5, /data, charsize=0.9, charthick=thisthick

; overlay the three metallicity regions

    r12 = 8.15
    r23 = 8.7

;   oplot, r12*[1,1], [-1.25,0.75], line=1, thick=postthick
;   oplot, r23*[1,1], [-1.25,0.75], line=1, thick=postthick
 
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
    
    sdss_lineplot, x, y, xerr, yerr, plottype=1, postscript=postscript, $
      xtitle=xtitle, ytitle='', xrange=xrange, yrange=yrange, legendtype=0, $
      /left, /top, position=pos[*,1], yfrac=1.5, /errorleft, $
      ytickname=replicate(' ',10), /noerase, charsize=charsize_5;, $
;     xregion=xregion, yregion=yregion, xerrregion=xerrregion, yerrregion=yerrregion

    if keyword_set(postscript) then thisthick = 6.0 else thickthick = 2.0

;   plot_kewley_grids, plotnumber=17, model=3, labeltype=4, /log12oh, /overplot, $
;     /noZgrid, thick=3.0, Ugridcolor=ugridcolor, Ulinestyle=0, $
;     postscript=postscript, Ulabelcolor='', Ucharsize=0.8, charthick=thisthick
;   xyouts, 9.1, -0.18, 'log U', align=0.5, /data, charsize=0.9, charthick=thisthick

; overlay the three metallicity regions

    r12 = 8.15
    r23 = 8.7

;   oplot, r12*[1,1], [-1.25,0.75], line=1, thick=postthick
;   oplot, r23*[1,1], [-1.25,0.75], line=1, thick=postthick
 
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
; 12+log(O/H) vs ([O II]/Ha)_cor [Integrated + SDSS]
; ------------------------------------------------------------

    psname = '12oh_oiiinii_niiha_vs_oiiha_2panel_nogrids_withlines'
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

    atlas1d_lineplot, x, y, xerr, yerr, plottype=1, postscript=postscript, $n
      xtitle=xtitle, ytitle=ytitle, xrange=xrange, yrange=yrange, legendtype=0, $
      /left, /top, position=pos[*,0], charsize=charsize_5, yfrac=1.5, /errorleft, $
      xregion=xregion, yregion=yregion, xerrregion=xerrregion, yerrregion=yerrregion, $
      xnfgs=xnfgs, ynfgs=ynfgs, xerrnfgs=xerrnfgs, yerrnfgs=yerrnfgs, blackwhite=blackwhite

    if keyword_set(postscript) then thisthick = 6.0 else thickthick = 2.0

;   plot_kewley_grids, plotnumber=17, model=3, labeltype=4, /log12oh, /overplot, $
;     /noZgrid, thick=3.0, Ugridcolor=ugridcolor, Ulinestyle=0, $
;     postscript=postscript, Ulabelcolor='', Ucharsize=0.8, charthick=thisthick
;   xyouts, 9.1, -0.18, 'log U', align=0.5, /data, charsize=0.9, charthick=thisthick

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
    
    sdss_lineplot, x, y, xerr, yerr, plottype=1, postscript=postscript, $
      xtitle=xtitle, ytitle='', xrange=xrange, yrange=yrange, legendtype=0, $
      /left, /top, position=pos[*,1], yfrac=1.5, /errorleft, $
      ytickname=replicate(' ',10), /noerase, charsize=charsize_5;, $
;     xregion=xregion, yregion=yregion, xerrregion=xerrregion, yerrregion=yerrregion

    if keyword_set(postscript) then thisthick = 6.0 else thickthick = 2.0

;   plot_kewley_grids, plotnumber=17, model=3, labeltype=4, /log12oh, /overplot, $
;     /noZgrid, thick=3.0, Ugridcolor=ugridcolor, Ulinestyle=0, $
;     postscript=postscript, Ulabelcolor='', Ucharsize=0.8, charthick=thisthick
;   xyouts, 9.1, -0.18, 'log U', align=0.5, /data, charsize=0.9, charthick=thisthick

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
; 12+log(O/H) vs ([O II]/Ha)_cor [Integrated + SDSS]
; ------------------------------------------------------------

    psname = '12oh_oiiinii_niiha_vs_oiiha_2panel_grids_withlines'
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

    atlas1d_lineplot, x, y, xerr, yerr, plottype=1, postscript=postscript, $n
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
    
    sdss_lineplot, x, y, xerr, yerr, plottype=1, postscript=postscript, $
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

    atlas1d_lineplot, x, y, xerr, yerr, plottype=1, postscript=postscript, $n
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

    sdss_lineplot, x, y, xerr, yerr, plottype=1, postscript=postscript, $
      xtitle=xtitle, ytitle='', xrange=xrange, yrange=yrange, legendtype=0, $
      ystyle=11, xmargin=[8,6], /left, /top, position=pos[*,1], $
      ytickname=replicate(' ',10), /noerase, charsize=charsize_4
    axis, /yaxis, yrange=interpol(yebv,y,!y.crange), ythick=postthick, $
      charthick=postthick, charsize=charsize_4, ystyle=1
    im_xyouts_title, ytitle='E(B-V) [mag]', charsize=charsize_4, charthick=postthick

    im_openclose, postscript=postscript, /close    

; ---------------------------------------------------------------------------    
; B magnitude distribution - Integrated + SDSS
; ---------------------------------------------------------------------------    

    psname = 'sdss_histogram_bmag'
    im_openclose, pspath+psname, postscript=postscript, xsize=8.5, ysize=8.8, encapsulated=encapsulated

    polyfill, [1,1,0,0,1], [1,0,0,1,1], /normal, color=djs_icolor('black')

    pagemaker, nx=1, ny=1, xspace=0, yspace=0.0, width=6.8, height=6.8, $
      xmargin=[1.3,0.4], ymargin=[0.9,1.1], xpage=8.5, ypage=8.8, $
      position=pos, /normal

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

    plot, [0], [0], /nodata, xthick=postthick, ythick=postthick, charsize=singlecharsize, $
      charthick=postthick, thick=postthick, xtitle=textoidl(xtitle), ytitle='Fraction', /noerase, $
      xstyle=11, ystyle=1, xrange=xrange, yrange=yrange, position=pos[*,0], color=djs_icolor(talkcolor)
    axis, /xaxis, xrange=interpol(xabs,x,!x.crange), xthick=postthick, xstyle=1, $
      charsize=singlecharsize, charthick=postthick, color=djs_icolor(talkcolor)
    im_xyouts_title, xtitle='M_{B} [mag]', charsize=singlecharsize, charthick=postthick, color=djs_icolor(talkcolor)

    im_plothist, x, bin=binsize, /overplot, thick=postthick4, /fraction, /halfbin, color=djs_icolor('green')
    im_plothist, xbig, bin=binsize, thick=postthick4, line=0, /overplot, /fraction, /halfbin, color=djs_icolor('red')

    legend, ['Integrated Sample','SDSS'], textcolor=djs_icolor([talkcolor,talkcolor]), $
      /left, /top, box=0, charsize=1.7, charthick=postthick4, line=[0,0], $
      color=djs_icolor(['red','green']), thick=postthick4, spacing=2.0
    
    im_openclose, postscript=postscript, /close    

; ---------------------------------------------------------------------------    
; A(Ha) distribution - Integrated + SDSS
; ---------------------------------------------------------------------------    

    psname = 'sdss_histogram_aha'
    im_openclose, pspath+psname, postscript=postscript, xsize=8.5, ysize=8.8, encapsulated=encapsulated

    polyfill, [1,1,0,0,1], [1,0,0,1,1], /normal, color=djs_icolor('black')

    pagemaker, nx=1, ny=1, xspace=0, yspace=0.0, width=6.8, height=6.8, $
      xmargin=[1.3,0.4], ymargin=[0.9,1.1], xpage=8.5, ypage=8.8, $
      position=pos, /normal

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

    plot, [0], [0], /nodata, xthick=postthick, ythick=postthick, charsize=singlecharsize, $
      charthick=postthick, thick=postthick, xtitle=textoidl(xtitle), ytitle='Fraction', $
      xstyle=11, ystyle=1, xrange=xrange, yrange=yrange, position=pos[*,0], /noerase, color=djs_icolor(talkcolor)
    axis, /xaxis, xrange=interpol(xebv,x,!x.crange), xthick=postthick, xstyle=1, $
      charsize=singlecharsize, charthick=postthick, color=djs_icolor(talkcolor)
    im_xyouts_title, xtitle='E(B-V) [mag]', charsize=singlecharsize, charthick=postthick, color=djs_icolor(talkcolor)

    im_plothist, x, bin=binsize, /overplot, thick=postthick4, /fraction, /halfbin, color=djs_icolor('green')
    im_plothist, xbig, bin=binsize, thick=postthick4, line=0, /overplot, /fraction, /halfbin, color=djs_icolor('red')

    legend, ['Integrated Sample','SDSS'], textcolor=djs_icolor([talkcolor,talkcolor]), $
      /right, /top, box=0, charsize=1.7, charthick=postthick4, line=[0,0], $
      color=djs_icolor(['red','green']), thick=postthick4, spacing=2.0
    
    im_openclose, postscript=postscript, /close    

; ---------------------------------------------------------------------------    
; 12+log(O/H) distribution - Integrated + SDSS
; ---------------------------------------------------------------------------    

    psname = 'sdss_histogram_log12oh'
    im_openclose, pspath+psname, postscript=postscript, xsize=8.5, ysize=8.3, encapsulated=encapsulated

    polyfill, [1,1,0,0,1], [1,0,0,1,1], /normal, color=djs_icolor('black')

    pagemaker, nx=1, ny=1, xspace=0, yspace=0.0, width=6.8, height=6.8, $
      xmargin=[1.3,0.4], ymargin=[0.4,1.1], xpage=8.5, ypage=8.3, $
      position=pos, /normal

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
    
    plot, [0], [0], /nodata, xthick=postthick, ythick=postthick, charsize=singlecharsize, $
      charthick=postthick, thick=postthick, xtitle=textoidl(xtitle), ytitle='Fraction', $
      xstyle=1, ystyle=1, xrange=xrange, yrange=yrange, position=pos[*,0], /noerase, color=djs_icolor(talkcolor)

    im_plothist, x, bin=binsize, /overplot, thick=postthick4, /fraction, /halfbin, color=djs_icolor('green')
    im_plothist, xbig, bin=binsize, thick=postthick4, line=0, /overplot, /fraction, /halfbin, color=djs_icolor('red')

    legend, ['Integrated Sample','SDSS'], textcolor=djs_icolor([talkcolor,talkcolor]), $
      /left, /top, box=0, charsize=1.7, charthick=postthick4, line=[0,0], $
      color=djs_icolor(['red','green']), thick=postthick4, spacing=2.0
    
    im_openclose, postscript=postscript, /close    

; ---------------------------------------------------------------------------    
; [O III]/[O II] distribution - Integrated + SDSS
; ---------------------------------------------------------------------------    

    psname = 'sdss_histogram_oiiioii'
    im_openclose, pspath+psname, postscript=postscript, xsize=8.5, ysize=8.3, encapsulated=encapsulated

    polyfill, [1,1,0,0,1], [1,0,0,1,1], /normal, color=djs_icolor('black')

    pagemaker, nx=1, ny=1, xspace=0, yspace=0.0, width=6.8, height=6.8, $
      xmargin=[1.3,0.4], ymargin=[0.4,1.1], xpage=8.5, ypage=8.3, $
      position=pos, /normal

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

    plot, [0], [0], /nodata, xthick=postthick, ythick=postthick, charsize=singlecharsize, $
      charthick=postthick, thick=postthick, xtitle=textoidl(xtitle), ytitle='Fraction', $
      xstyle=1, ystyle=1, xrange=xrange, yrange=yrange, position=pos[*,0], /noerase, color=djs_icolor(talkcolor)

    im_plothist, x, bin=binsize, /overplot, thick=postthick4, /fraction, /halfbin, color=djs_icolor('green')
    im_plothist, xbig, bin=binsize, thick=postthick4, line=0, /overplot, /fraction, /halfbin, color=djs_icolor('red')

    legend, ['Integrated Sample','SDSS'], textcolor=djs_icolor([talkcolor,talkcolor]), $
      /left, /top, box=0, charsize=1.7, charthick=postthick4, line=[0,0], $
      color=djs_icolor(['red','green']), thick=postthick4, spacing=2.0
    
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
    
    atlas1d_lineplot, x, y, xerr, yerr, plottype=1, postscript=postscript, ysty=3, $
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

    sdss_lineplot, x, y, xerr, yerr, plottype=1, postscript=postscript, $
      ytitle='', xtitle=xtitle, xrange=xrange, yrange=yrange, legendtype=0, $
      xstyle=11, /right, /top, position=pos[*,1], charsize=charsize_6, /noerase, $
      ytickname=replicate(' ',10)
    axis, /xaxis, xrange=interpol(xabs,x,!x.crange), xthick=postthick, xstyle=1, $
      charsize=charsize_6, charthick=postthick
    im_xyouts_title, xtitle='M_{B} [mag]', charsize=charsize_6, charthick=postthick

    im_openclose, postscript=postscript, /close

; --------------------------------------------------    
; GENERATE THE WEB PAGE
; --------------------------------------------------    
    
    if keyword_set(postscript) then begin
       im_ps2html, htmlbase, html_path=html_path, cleanpng=0, npscols=3, _extra=extra, /psfind
    endif

stop    
    
return
end    

