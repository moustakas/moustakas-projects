; OLD PLOTS, SOME USEABLE - jm05sep17uofa

; ------------------------------------------------------------
; L(B) versus SFR(Ha)/L(Hb)_obs - Integrated
; ------------------------------------------------------------

    psname = 'mass_vs_sfr_ha_lhb_obs'
    im_openclose, pspath+psname, postscript=postscript, xsize=8.5, ysize=8.8, /encapsulated

    pagemaker, nx=1, ny=1, xspace=0, yspace=0.0, width=6.8, height=6.8, $
      xmargin=[1.4,0.3], ymargin=[0.9,1.1], xpage=8.5, ypage=8.8, $
      position=pos, /normal

    medbin = 0.5
    minpts = 3
    minx = 7.0 ; -3.0
    medcolor = 'dark blue'
    
; Atlas - SF
    
    indx = where((atlasdust.mass_bv_b gt -900) and (atlasnodust.sfr_h_alpha gt -900.0) and $
      (atlasdust.h_beta_ew_uncor[0] ge ewcut),nindx)

    x = atlasdust[indx].mass_bv_b; - alog10(LBnorm)
    xerr = atlasdust[indx].mass_bv_b_err
    xabs = atlasdust[indx].m_b_obs

    hb = atlasdust[indx].h_beta_lum[0] + alog10(lsun) - alog10(elumnorm)
    hb_err = atlasdust[indx].h_beta_lum[1]

    sfr_ha = atlasnodust[indx].sfr_h_alpha
    sfr_ha_err = atlasnodust[indx].sfr_h_alpha_err

    y = sfr_ha - hb
    yerr = sqrt(sfr_ha_err^2 + hb_err^2)
    
; Atlas - AGN
    
    indx_agn = where((atlasdust_agn.mass_bv_b gt -900) and (atlasnodust_agn.sfr_h_alpha gt -900.0) and $
      (atlasdust_agn.h_beta_ew_uncor[0] ge ewcut),nindx_agn)

    x_agn = atlasdust_agn[indx_agn].mass_bv_b; - alog10(LBnorm)
    xerr_agn = atlasdust_agn[indx_agn].mass_bv_b_err

    hb_agn = atlasdust_agn[indx_agn].h_beta_lum[0] + alog10(lsun) - alog10(elumnorm)
    hb_err_agn = atlasdust_agn[indx_agn].h_beta_lum[1]

    sfr_ha_agn = atlasnodust_agn[indx_agn].sfr_h_alpha
    sfr_ha_err_agn = atlasnodust_agn[indx_agn].sfr_h_alpha_err

    y_agn = sfr_ha_agn - hb_agn
    yerr_agn = sqrt(sfr_ha_err_agn^2 + hb_err_agn^2)
    
; NFGS - SF
    
    indxnfgs = where((nfgsdust.mass_bv_b gt -900.0) and (nfgsnodust.sfr_h_alpha gt -900.0) and $
      (nfgsdust.h_beta_ew_uncor[0] ge ewcut),nindxnfgs)

    xnfgs = nfgsdust[indxnfgs].mass_bv_b; - alog10(LBnorm)
    xerrnfgs = nfgsdust[indxnfgs].mass_bv_b_err
    
    hbnfgs = nfgsdust[indxnfgs].h_beta_lum[0] + alog10(lsun) - alog10(elumnorm)
    hbnfgs_err = nfgsdust[indxnfgs].h_beta_lum[1]

    sfrnfgs_ha = nfgsnodust[indxnfgs].sfr_h_alpha
    sfrnfgs_ha_err = nfgsnodust[indxnfgs].sfr_h_alpha_err

    ynfgs = sfrnfgs_ha - hbnfgs
    yerrnfgs = sqrt(sfrnfgs_ha_err^2 + hbnfgs_err^2)

; NFGS - AGN
    
    indxnfgs_agn = where((nfgsdust_agn.mass_bv_b gt -900.0) and (nfgsnodust_agn.sfr_h_alpha gt -900.0) and $
      (nfgsdust_agn.h_beta_ew_uncor[0] ge ewcut),nindxnfgs_agn)

    xnfgs_agn = nfgsdust_agn[indxnfgs_agn].mass_bv_b; - alog10(LBnorm)
    xerrnfgs_agn = nfgsdust_agn[indxnfgs_agn].mass_bv_b_err
    
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

    xtitle = 'log M [M'+sunsymbol()+']'
    ytitle = 'log [10^{41} \psi/L(H\beta)_{obs}] [erg s^{-1}/'+sfr_units()+']'

    xrange = massrange
    yrange = alog10([1.8D-41,3D-40]) - alog10(sfrnorm) + alog10(elumnorm)

; SF galaxies    
    
    atlas1d_lineplot, xbig[notlirgs], ybig[notlirgs], xerrbig[notlirgs], yerrbig[notlirgs], $
      plottype=1, postscript=postscript, xtitle=xtitle, ytitle=ytitle, xrange=xrange, $
      yrange=yrange, legendtype=0, /right, /top, yfrac=1.5, xfrac=1.5, charsize=charsize_8, $
      atlaspsize=1.3, position=pos[*,0], atlascolor='light green'

    im_symbols, 108, psize=1.6, /fill, color='magenta'
    djs_oplot, xbig[lirgs], ybig[lirgs], ps=8

    im_symbols, 105, psize=1.5, /fill, color='light blue'
    djs_oplot, xbig[nolirdata], ybig[nolirdata], ps=8

; AGN galaxies    
    
    atlas1d_lineplot, xbig_agn[notlirgs_agn], ybig_agn[notlirgs_agn], xerrbig_agn[notlirgs_agn], yerrbig_agn[notlirgs_agn], $
      plottype=1, postscript=postscript, xtitle=xtitle, ytitle=ytitle, xrange=xrange, $
      yrange=yrange, legendtype=0, /right, /top, yfrac=1.5, xfrac=1.5, charsize=charsize_8, $
      atlaspsize=1.1, position=pos[*,0], atlasfill=0, /overplot, thick=postthick, atlascolor='light green'

    im_symbols, 108, psize=1.4, fill=0, color='magenta', thick=postthick
    djs_oplot, xbig_agn[lirgs_agn], ybig_agn[lirgs_agn], ps=8

    im_symbols, 105, psize=1.3, fill=0, color='light blue', thick=postthick
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
    oploterror, running.binctr[doit], running.medy[doit], running.binsz/2.0+fltarr(ndoit), $
      running.sigy75[doit]-running.medy[doit], ps=8, errthick=(postthick-1L)>2L, /hibar, $
      color=djs_icolor(medcolor), errcolor=djs_icolor(medcolor)
    oploterror, running.binctr[doit], running.medy[doit], running.binsz/2.0+fltarr(ndoit), $
      running.medy[doit]-running.sigy25[doit], ps=8, errthick=(postthick-1L)>2L, /lobar, $
      color=djs_icolor(medcolor), errcolor=djs_icolor(medcolor)
    
    im_openclose, postscript=postscript, /close    

; ------------------------------------------------------------
; SFR(Ha) vs various SFR(Hb) calibrations
; ------------------------------------------------------------
    
    psname = 'sfr_ha_vs_sfr_hb_3panel'
    im_openclose, pspath+psname, postscript=postscript, /encapsulated, xsize=7.0, ysize=9.0

    pagemaker, nx=1, ny=3, height=[2.5,2.5,2.5], width=5.5, xmargin=[1.1,0.4], $
      ymargin=[0.4,1.1], xspace=0, yspace=0, xpage=7.0, ypage=9.0, $
      position=pos, /normal

    xrange = sfrharange
    yrange = residrange_hbsfrs

    xtitle = 'log \psi(H\alpha)'
    ytitle = 'log \psi(H\beta)/\psi(H\alpha)'
    
; ##########################################################
; Panel 1
; ##########################################################

; Atlas    
    
    indx = where((atlas_sfrs.sfr_ha gt -900) and (atlas_sfrs.sfr_hb_uncor gt -900),nindx)

    x = atlas_sfrs[indx].sfr_ha
    xerr = atlas_sfrs[indx].sfr_ha_err
    
    y1 = atlas_sfrs[indx].sfr_ha
    y1err = atlas_sfrs[indx].sfr_ha_err
    
    y2 = atlas_sfrs[indx].sfr_hb_uncor
    y2err = atlas_sfrs[indx].sfr_hb_uncor_err

    resid = y2-y1
    resid_err = sqrt(y1err^2+y2err^2)

; NFGS    
    
    indxnfgs = where((nfgs_sfrs.sfr_ha gt -900) and (nfgs_sfrs.sfr_hb_uncor gt -900),nindxnfgs)

    xnfgs = nfgs_sfrs[indxnfgs].sfr_ha
    xerrnfgs = nfgs_sfrs[indxnfgs].sfr_ha_err
    
    y1nfgs = nfgs_sfrs[indxnfgs].sfr_ha
    y1errnfgs = nfgs_sfrs[indxnfgs].sfr_ha_err
    
    y2nfgs = nfgs_sfrs[indxnfgs].sfr_hb_uncor
    y2errnfgs = nfgs_sfrs[indxnfgs].sfr_hb_uncor_err

    residnfgs = y2nfgs-y1nfgs
    residnfgs_err = sqrt(y1errnfgs^2+y2errnfgs^2)

; combine the two samples

    xbig = [x,xnfgs]
    xerrbig = [xerr,xerrnfgs]
    residbig = [resid,residnfgs]
    residerrbig = [resid_err,residnfgs_err]
    
; make the plot    

    stats = im_stats(residbig,/verbose,/baremin)
    xstr = strtrim(string(stats.median_rej,format='(F12.2)'),2)+' ('+$
      strtrim(string(stats.mean_rej,format='(F12.2)'),2)+$
      '\pm'+strtrim(string(stats.sigma_rej,format='(F12.2)'),2)+')'
;   xstr = strtrim(string(stats.median,format='(F12.2)'),2)+'\pm'+strtrim(string(stats.sig68mean,format='(F12.2)'),2)
    
    atlas1d_lineplot, xbig, residbig, xerrbig, residerrbig, plottype=3, postscript=postscript, $
      xtitle='', ytitle='', xrange=xrange, yrange=yrange, legendtype=0, $
      charsize=charsize_7, xtickname=replicate(' ',10), position=pos[*,0], $
      xminor=3, yminor=3, atlaspsize=0.8, ytickinterval=1.0;, nfgspsize=0.8, $
;     xnfgs=xnfgs, ynfgs=residnfgs, xerrnfgs=xerrnfgs, yerrnfgs=residnfgs_err
    djs_oplot, !x.crange, [0,0], line=2, thick=postthick
    
    legend, textoidl('(a) Observed L(H\beta)'), /left, /top, box=0, charsize=charsize_3, $
      charthick=postthick, clear=keyword_set(postscript)
    legend, textoidl(xstr), /right, /top, box=0, charsize=charsize_3, charthick=postthick

; ##########################################################
; Panel 2
; ##########################################################

; Atlas    
    
    indx = where((atlas_sfrs.sfr_ha gt -900) and (atlas_sfrs.sfr_hb_hbhg gt -900),nindx)

    x = atlas_sfrs[indx].sfr_ha
    xerr = atlas_sfrs[indx].sfr_ha_err
    
    y1 = atlas_sfrs[indx].sfr_ha
    y1err = atlas_sfrs[indx].sfr_ha_err
    
    y2 = atlas_sfrs[indx].sfr_hb_hbhg
    y2err = atlas_sfrs[indx].sfr_hb_hbhg_err

    resid = y2-y1
    resid_err = sqrt(y1err^2+y2err^2)

; NFGS    
    
    indxnfgs = where((nfgs_sfrs.sfr_ha gt -900) and (nfgs_sfrs.sfr_hb_hbhg gt -900),nindxnfgs)

    xnfgs = nfgs_sfrs[indxnfgs].sfr_ha
    xerrnfgs = nfgs_sfrs[indxnfgs].sfr_ha_err
    
    y1nfgs = nfgs_sfrs[indxnfgs].sfr_ha
    y1errnfgs = nfgs_sfrs[indxnfgs].sfr_ha_err
    
    y2nfgs = nfgs_sfrs[indxnfgs].sfr_hb_hbhg
    y2errnfgs = nfgs_sfrs[indxnfgs].sfr_hb_hbhg_err

    residnfgs = y2nfgs-y1nfgs
    residnfgs_err = sqrt(y1errnfgs^2+y2errnfgs^2)

; combine the two samples

    xbig = [x,xnfgs]
    xerrbig = [xerr,xerrnfgs]
    residbig = [resid,residnfgs]
    residerrbig = [resid_err,residnfgs_err]
    
; make the plot    

    stats = im_stats(residbig,/verbose,/baremin,/no_head)
    xstr = strtrim(string(stats.median_rej,format='(F12.2)'),2)+' ('+$
      strtrim(string(stats.mean_rej,format='(F12.2)'),2)+$
      '\pm'+strtrim(string(stats.sigma_rej,format='(F12.2)'),2)+')'
;   xstr = strtrim(string(stats.median,format='(F12.2)'),2)+'\pm'+strtrim(string(stats.sig68mean,format='(F12.2)'),2)
    
    atlas1d_lineplot, xbig, residbig, xerrbig, residerrbig, plottype=3, postscript=postscript, $
      xtitle='', ytitle=ytitle, xrange=xrange, yrange=yrange, legendtype=0, $
      charsize=charsize_7, xtickname=replicate(' ',10), position=pos[*,1], $
      xminor=3, yminor=3, /noerase, atlaspsize=0.8, ytickinterval=1.0;, nfgspsize=0.8, $
;     xnfgs=xnfgs, ynfgs=residnfgs, xerrnfgs=xerrnfgs, yerrnfgs=residnfgs_err
    djs_oplot, !x.crange, [0,0], line=2, thick=postthick
    
    legend, textoidl('(b) L(H\beta) corrected using H\beta/H\gamma'), /left, $
      /top, box=0, charsize=charsize_3, charthick=postthick, clear=keyword_set(postscript)
    legend, textoidl(xstr), /right, /top, box=0, charsize=charsize_3, charthick=postthick

; ##########################################################
; Panel 3
; ##########################################################

; Atlas    
    
    indx = where((atlas_sfrs.sfr_ha gt -900) and (atlas_sfrs.sfr_hb_best gt -900),nindx)

    x = atlas_sfrs[indx].sfr_ha
    xerr = atlas_sfrs[indx].sfr_ha_err
    
    y1 = atlas_sfrs[indx].sfr_ha
    y1err = atlas_sfrs[indx].sfr_ha_err
    
    y2 = atlas_sfrs[indx].sfr_hb_best
    y2err = atlas_sfrs[indx].sfr_hb_best_err

    resid = y2-y1
    resid_err = sqrt(y1err^2+y2err^2)

; NFGS    
    
    indxnfgs = where((nfgs_sfrs.sfr_ha gt -900) and (nfgs_sfrs.sfr_hb_best gt -900),nindxnfgs)

    xnfgs = nfgs_sfrs[indxnfgs].sfr_ha
    xerrnfgs = nfgs_sfrs[indxnfgs].sfr_ha_err
    
    y1nfgs = nfgs_sfrs[indxnfgs].sfr_ha
    y1errnfgs = nfgs_sfrs[indxnfgs].sfr_ha_err
    
    y2nfgs = nfgs_sfrs[indxnfgs].sfr_hb_best
    y2errnfgs = nfgs_sfrs[indxnfgs].sfr_hb_best_err

    residnfgs = y2nfgs-y1nfgs
    residnfgs_err = sqrt(y1errnfgs^2+y2errnfgs^2)

; combine the two samples

    xbig = [x,xnfgs]
    xerrbig = [xerr,xerrnfgs]
    residbig = [resid,residnfgs]
    residerrbig = [resid_err,residnfgs_err]
    
; make the plot    

    stats = im_stats(residbig,/verbose,/baremin,/no_head)
    xstr = strtrim(string(stats.median_rej,format='(F12.2)'),2)+' ('+$
      strtrim(string(stats.mean_rej,format='(F12.2)'),2)+$
      '\pm'+strtrim(string(stats.sigma_rej,format='(F12.2)'),2)+')'
;   xstr = strtrim(string(stats.median,format='(F12.2)'),2)+'\pm'+strtrim(string(stats.sig68mean,format='(F12.2)'),2)

;   yrange = [-1,1]
    
    atlas1d_lineplot, xbig, residbig, xerrbig, residerrbig, plottype=3, postscript=postscript, $
      xtitle=xtitle, ytitle='', xrange=xrange, yrange=yrange, legendtype=0, $
      charsize=charsize_7, position=pos[*,2], $
      xminor=3, yminor=3, /noerase, atlaspsize=0.8, ytickinterval=1.0;, nfgspsize=0.8, $
;     xnfgs=xnfgs, ynfgs=residnfgs, xerrnfgs=xerrnfgs, yerrnfgs=residnfgs_err
    djs_oplot, !x.crange, [0,0], line=2, thick=postthick
    
    legend, '(c) This Paper', /left, /top, box=0, charsize=charsize_3, $
      charthick=postthick, clear=keyword_set(postscript)
    legend, textoidl(xstr), /right, /top, box=0, charsize=charsize_3, charthick=postthick

    im_openclose, postscript=postscript, /close    

; ------------------------------------------------------------
; SFR(Ha) vs various SFR([O II]) calibrations - residual plots
; ------------------------------------------------------------

    psname = 'sfr_ha_vs_sfr_oii_6panel'
    im_openclose, pspath+psname, postscript=postscript, xsize=8.5, ysize=8.75, /encapsulated

    pagemaker, nx=2, ny=3, yspace=0, xspace=0, width=3.15*[1,1], height=2.4*[1,1,1], $
      xmargin=[1.1,1.1], ymargin=[0.45,1.1], xpage=8.5, ypage=8.75, position=pos, /normal

    xtitle = 'log \psi(H\alpha) ['+sfr_units()+']'
    ytitle = 'log \psi([O II])/\psi(H\alpha)'
    
    xrange = sfrharange
    yrange = residrange_sfrs
    
; ##########################################################
; Panel 1
; ##########################################################

; Atlas    
    
    indx = where((atlas_sfrs.sfr_ha gt -900) and (atlas_sfrs.sfr_oii_uncor gt -900),nindx)

    x = atlas_sfrs[indx].sfr_ha
    xerr = atlas_sfrs[indx].sfr_ha_err
    
    y = atlas_sfrs[indx].sfr_oii_uncor
    yerr = atlas_sfrs[indx].sfr_oii_uncor_err

    resid = y-x
    resid_err = sqrt(yerr^2+xerr^2)

; NFGS    
    
    indxnfgs = where((nfgs_sfrs.sfr_ha gt -900) and (nfgs_sfrs.sfr_oii_uncor gt -900),nindxnfgs)

    xnfgs = nfgs_sfrs[indxnfgs].sfr_ha
    xerrnfgs = nfgs_sfrs[indxnfgs].sfr_ha_err
    
    ynfgs = nfgs_sfrs[indxnfgs].sfr_oii_uncor
    yerrnfgs = nfgs_sfrs[indxnfgs].sfr_oii_uncor_err

    residnfgs = ynfgs-xnfgs
    residnfgs_err = sqrt(yerrnfgs^2+xerrnfgs^2)

; combine the two samples

    xbig = [x,xnfgs]
    xerrbig = [xerr,xerrnfgs]
    residbig = [resid,residnfgs]
    residerrbig = [resid_err,residnfgs_err]
    
; make the plot    

    stats = im_stats(residbig,/verbose,/baremin)
    xstr = strtrim(string(stats.median_rej,format='(F12.2)'),2)+' ('+$
      strtrim(string(stats.mean_rej,format='(F12.2)'),2)+$
      '\pm'+strtrim(string(stats.sigma_rej,format='(F12.2)'),2)+')'
;   xstr = strtrim(string(stats.median,format='(F12.2)'),2)+'\pm'+strtrim(string(stats.sig68mean,format='(F12.2)'),2)
    
    atlas1d_lineplot, xbig, residbig, xerrbig, residerrbig, plottype=3, postscript=postscript, $
;   atlas1d_lineplot, x, resid, xerr, resid_err, plottype=3, postscript=postscript, $
      xtitle='', ytitle='', xrange=xrange, yrange=yrange, legendtype=0, $
      charsize=charsize_5, xtickname=replicate(' ',10), position=pos[*,0], $
      xminor=3, yminor=3, atlaspsize=0.5;, nfgspsize=0.8, $
;     xnfgs=xnfgs, ynfgs=residnfgs, xerrnfgs=xerrnfgs, yerrnfgs=residnfgs_err
    djs_oplot, !x.crange, [0,0], line=0, thick=postthick
    
    legend, '(a) '+textoidl('L([O II])_{obs}'), /left, /top, box=0, charsize=charsize_3, $
      charthick=postthick, clear=keyword_set(postscript)
    legend, textoidl(xstr), /left, /bottom, box=0, charsize=charsize_3, $
      charthick=postthick, clear=keyword_set(postscript)

; ##########################################################
; Panel 2
; ##########################################################

; Atlas    
    
    indx = where((atlas_sfrs.sfr_ha gt -900) and (atlas_sfrs.sfr_oii_k98 gt -900),nindx)

    x = atlas_sfrs[indx].sfr_ha
    xerr = atlas_sfrs[indx].sfr_ha_err
    
    y = atlas_sfrs[indx].sfr_oii_k98
    yerr = atlas_sfrs[indx].sfr_oii_k98_err

    resid = y-x
    resid_err = sqrt(yerr^2+xerr^2)

; NFGS    
    
    indxnfgs = where((nfgs_sfrs.sfr_ha gt -900) and (nfgs_sfrs.sfr_oii_k98 gt -900),nindxnfgs)

    xnfgs = nfgs_sfrs[indxnfgs].sfr_ha
    xerrnfgs = nfgs_sfrs[indxnfgs].sfr_ha_err
    
    ynfgs = nfgs_sfrs[indxnfgs].sfr_oii_k98
    yerrnfgs = nfgs_sfrs[indxnfgs].sfr_oii_k98_err

    residnfgs = ynfgs-xnfgs
    residnfgs_err = sqrt(yerrnfgs^2+xerrnfgs^2)

; combine the two samples

    xbig = [x,xnfgs]
    xerrbig = [xerr,xerrnfgs]
    residbig = [resid,residnfgs]
    residerrbig = [resid_err,residnfgs_err]
    
; make the plot    

    stats = im_stats(residbig,/verbose,/baremin,/no_head)
    xstr = strtrim(string(stats.median_rej,format='(F12.2)'),2)+' ('+$
      strtrim(string(stats.mean_rej,format='(F12.2)'),2)+$
      '\pm'+strtrim(string(stats.sigma_rej,format='(F12.2)'),2)+')'
;   xstr = strtrim(string(stats.median,format='(F12.2)'),2)+'\pm'+strtrim(string(stats.sig68mean,format='(F12.2)'),2)
    
    atlas1d_lineplot, xbig, residbig, xerrbig, residerrbig, plottype=3, postscript=postscript, $
;   atlas1d_lineplot, x, resid, xerr, resid_err, plottype=3, postscript=postscript, $
      xtitle='', ytitle='', xrange=xrange, yrange=yrange, legendtype=0, $
      charsize=charsize_5, xtickname=replicate(' ',10), position=pos[*,1], $
      ytickname=replicate(' ',10), xminor=3, yminor=3, ystyle=11, /noerase, $
      atlaspsize=0.5;, nfgspsize=0.8, xnfgs=xnfgs, ynfgs=residnfgs, xerrnfgs=xerrnfgs, yerrnfgs=residnfgs_err
    axis, /yaxis, yrange=yrange, ythick=postthick, ytitle='', $
      charthick=postthick, charsize=charsize_5, ysty=1, yminor=3
    im_xyouts_title, ytitle='', charsize=charsize_5, charthick=postthick
    djs_oplot, !x.crange, [0,0], line=0, thick=postthick
    
    legend, '(b) K98', /left, /top, box=0, charsize=charsize_3, $
      charthick=postthick, clear=keyword_set(postscript)
    legend, textoidl(xstr), /left, /bottom, box=0, charsize=charsize_3, $
      charthick=postthick, clear=keyword_set(postscript)

; ##########################################################
; Panel 3
; ##########################################################

; Atlas    
    
    indx = where((atlas_sfrs.sfr_ha gt -900) and (atlas_sfrs.sfr_oii_k04 gt -900),nindx)

    x = atlas_sfrs[indx].sfr_ha
    xerr = atlas_sfrs[indx].sfr_ha_err
    
    y = atlas_sfrs[indx].sfr_oii_k04
    yerr = atlas_sfrs[indx].sfr_oii_k04_err

    resid = y-x
    resid_err = sqrt(yerr^2+xerr^2)

; NFGS    
    
    indxnfgs = where((nfgs_sfrs.sfr_ha gt -900) and (nfgs_sfrs.sfr_oii_k04 gt -900),nindxnfgs)

    xnfgs = nfgs_sfrs[indxnfgs].sfr_ha
    xerrnfgs = nfgs_sfrs[indxnfgs].sfr_ha_err
    
    ynfgs = nfgs_sfrs[indxnfgs].sfr_oii_k04
    yerrnfgs = nfgs_sfrs[indxnfgs].sfr_oii_k04_err

    residnfgs = ynfgs-xnfgs
    residnfgs_err = sqrt(yerrnfgs^2+xerrnfgs^2)

; combine the two samples

    xbig = [x,xnfgs]
    xerrbig = [xerr,xerrnfgs]
    residbig = [resid,residnfgs]
    residerrbig = [resid_err,residnfgs_err]
    
; make the plot    

    stats = im_stats(residbig,/verbose,/baremin,/no_head)
    xstr = strtrim(string(stats.median_rej,format='(F12.2)'),2)+' ('+$
      strtrim(string(stats.mean_rej,format='(F12.2)'),2)+$
      '\pm'+strtrim(string(stats.sigma_rej,format='(F12.2)'),2)+')'
;   xstr = strtrim(string(stats.median,format='(F12.2)'),2)+'\pm'+strtrim(string(stats.sig68mean,format='(F12.2)'),2)
    
    atlas1d_lineplot, xbig, residbig, xerrbig, residerrbig, plottype=3, postscript=postscript, $
;   atlas1d_lineplot, x, resid, xerr, resid_err, plottype=3, postscript=postscript, $
      xtitle='', ytitle=ytitle, xrange=xrange, yrange=yrange, legendtype=0, $
      charsize=charsize_5, xtickname=replicate(' ',10), position=pos[*,2], $
      xminor=3, yminor=3, /noerase, atlaspsize=0.5;, nfgspsize=0.8, $
;     xnfgs=xnfgs, ynfgs=residnfgs, xerrnfgs=xerrnfgs, yerrnfgs=residnfgs_err
    djs_oplot, !x.crange, [0,0], line=0, thick=postthick
    
    legend, '(c) KGJ04', /left, /top, box=0, charsize=charsize_3, $
      charthick=postthick, clear=keyword_set(postscript)
    legend, textoidl(xstr), /left, /bottom, box=0, charsize=charsize_3, $
      charthick=postthick, clear=keyword_set(postscript)

; ##########################################################
; Panel 4
; ##########################################################

; Atlas    
    
    indx = where((atlas_sfrs.sfr_ha gt -900) and (atlas_sfrs.sfr_oii_hbhg gt -900),nindx)

    x = atlas_sfrs[indx].sfr_ha
    xerr = atlas_sfrs[indx].sfr_ha_err
    
    y = atlas_sfrs[indx].sfr_oii_hbhg
    yerr = atlas_sfrs[indx].sfr_oii_hbhg_err

    resid = y-x
    resid_err = sqrt(yerr^2+xerr^2)

; NFGS    
    
    indxnfgs = where((nfgs_sfrs.sfr_ha gt -900) and (nfgs_sfrs.sfr_oii_hbhg gt -900),nindxnfgs)

    xnfgs = nfgs_sfrs[indxnfgs].sfr_ha
    xerrnfgs = nfgs_sfrs[indxnfgs].sfr_ha_err
    
    ynfgs = nfgs_sfrs[indxnfgs].sfr_oii_hbhg
    yerrnfgs = nfgs_sfrs[indxnfgs].sfr_oii_hbhg_err

    residnfgs = ynfgs-xnfgs
    residnfgs_err = sqrt(yerrnfgs^2+xerrnfgs^2)

; combine the two samples

    xbig = [x,xnfgs]
    xerrbig = [xerr,xerrnfgs]
    residbig = [resid,residnfgs]
    residerrbig = [resid_err,residnfgs_err]
    
; make the plot    

    stats = im_stats(residbig,/verbose,/baremin,/no_head)
    xstr = strtrim(string(stats.median_rej,format='(F12.2)'),2)+' ('+$
      strtrim(string(stats.mean_rej,format='(F12.2)'),2)+$
      '\pm'+strtrim(string(stats.sigma_rej,format='(F12.2)'),2)+')'
;   xstr = strtrim(string(stats.median,format='(F12.2)'),2)+'\pm'+strtrim(string(stats.sig68mean,format='(F12.2)'),2)
    
    atlas1d_lineplot, xbig, residbig, xerrbig, residerrbig, plottype=3, postscript=postscript, $
;   atlas1d_lineplot, x, resid, xerr, resid_err, plottype=3, postscript=postscript, $
      xtitle='', ytitle='', xrange=xrange, yrange=yrange, legendtype=0, $
      charsize=charsize_5, xtickname=replicate(' ',10), position=pos[*,3], $
      ytickname=replicate(' ',10), xminor=3, yminor=3, ystyle=11, /noerase, $
      atlaspsize=0.5;, nfgspsize=0.8, xnfgs=xnfgs, ynfgs=residnfgs, xerrnfgs=xerrnfgs, yerrnfgs=residnfgs_err
    axis, /yaxis, yrange=yrange, ythick=postthick, ytitle='', $
      charthick=postthick, charsize=charsize_5, ysty=1, yminor=3
    im_xyouts_title, ytitle=ytitle, charsize=charsize_5, charthick=postthick
    djs_oplot, !x.crange, [0,0], line=0, thick=postthick
    
    legend, textoidl('(d) H\beta/H\gamma'), /left, /top, box=0, charsize=charsize_3, $
      charthick=postthick, clear=keyword_set(postscript)
    legend, textoidl(xstr), /left, /bottom, box=0, charsize=charsize_3, $
      charthick=postthick, clear=keyword_set(postscript)

; ##########################################################
; Panel 5
; ##########################################################

; Atlas    
    
    indx = where((atlas_sfrs.sfr_ha gt -900) and (atlas_sfrs.sfr_oii_best gt -900),nindx)

    x = atlas_sfrs[indx].sfr_ha
    xerr = atlas_sfrs[indx].sfr_ha_err
    
    y = atlas_sfrs[indx].sfr_oii_best
    yerr = atlas_sfrs[indx].sfr_oii_best_err

    resid = y-x
    resid_err = sqrt(yerr^2+xerr^2)

; NFGS    
    
    indxnfgs = where((nfgs_sfrs.sfr_ha gt -900) and (nfgs_sfrs.sfr_oii_best gt -900),nindxnfgs)

    xnfgs = nfgs_sfrs[indxnfgs].sfr_ha
    xerrnfgs = nfgs_sfrs[indxnfgs].sfr_ha_err
    
    ynfgs = nfgs_sfrs[indxnfgs].sfr_oii_best
    yerrnfgs = nfgs_sfrs[indxnfgs].sfr_oii_best_err

    residnfgs = ynfgs-xnfgs
    residnfgs_err = sqrt(yerrnfgs^2+xerrnfgs^2)

; combine the two samples

    xbig = [x,xnfgs]
    xerrbig = [xerr,xerrnfgs]
    residbig = [resid,residnfgs]
    residerrbig = [resid_err,residnfgs_err]
    
; make the plot    

    stats = im_stats(residbig,/verbose,/baremin,/no_head)
    xstr = strtrim(string(stats.median_rej,format='(F12.2)'),2)+' ('+$
      strtrim(string(stats.mean_rej,format='(F12.2)'),2)+$
      '\pm'+strtrim(string(stats.sigma_rej,format='(F12.2)'),2)+')'
;   xstr = strtrim(string(stats.median,format='(F12.2)'),2)+'\pm'+strtrim(string(stats.sig68mean,format='(F12.2)'),2)
    
    atlas1d_lineplot, xbig, residbig, xerrbig, residerrbig, plottype=3, postscript=postscript, $
;   atlas1d_lineplot, x, resid, xerr, resid_err, plottype=3, postscript=postscript, $
      xtitle=xtitle, ytitle='', xrange=xrange, yrange=yrange, legendtype=0, $
      charsize=charsize_5, position=pos[*,4], $
      xminor=3, yminor=3, /noerase, atlaspsize=0.5;, nfgspsize=0.8, $
;     xnfgs=xnfgs, ynfgs=residnfgs, xerrnfgs=xerrnfgs, yerrnfgs=residnfgs_err
    djs_oplot, !x.crange, [0,0], line=0, thick=postthick
    
    legend, '(e) This Paper', /left, /top, box=0, charsize=charsize_3, $
      charthick=postthick, clear=keyword_set(postscript)
    legend, textoidl(xstr), /left, /bottom, box=0, charsize=charsize_3, $
      charthick=postthick, clear=keyword_set(postscript)

; ##########################################################
; Panel 6
; ##########################################################

; Atlas    
    
    indx = where((atlas_sfrs.sfr_ha gt -900) and (atlas_sfrs.sfr_oii_hahb gt -900),nindx)

    x = atlas_sfrs[indx].sfr_ha
    xerr = atlas_sfrs[indx].sfr_ha_err
    
    y = atlas_sfrs[indx].sfr_oii_hahb
    yerr = atlas_sfrs[indx].sfr_oii_hahb_err

    resid = y-x
    resid_err = sqrt(yerr^2+xerr^2)

; NFGS    
    
    indxnfgs = where((nfgs_sfrs.sfr_ha gt -900) and (nfgs_sfrs.sfr_oii_hahb gt -900),nindxnfgs)

    xnfgs = nfgs_sfrs[indxnfgs].sfr_ha
    xerrnfgs = nfgs_sfrs[indxnfgs].sfr_ha_err
    
    ynfgs = nfgs_sfrs[indxnfgs].sfr_oii_hahb
    yerrnfgs = nfgs_sfrs[indxnfgs].sfr_oii_hahb_err

    residnfgs = ynfgs-xnfgs
    residnfgs_err = sqrt(yerrnfgs^2+xerrnfgs^2)

; combine the two samples

    xbig = [x,xnfgs]
    xerrbig = [xerr,xerrnfgs]
    residbig = [resid,residnfgs]
    residerrbig = [resid_err,residnfgs_err]
    
; make the plot    

    stats = im_stats(residbig,/verbose,/baremin,/no_head)
    xstr = strtrim(string(stats.median_rej,format='(F12.2)'),2)+' ('+$
      strtrim(string(stats.mean_rej,format='(F12.2)'),2)+$
      '\pm'+strtrim(string(stats.sigma_rej,format='(F12.2)'),2)+')'
;   xstr = strtrim(string(stats.median,format='(F12.2)'),2)+'\pm'+strtrim(string(stats.sig68mean,format='(F12.2)'),2)
    
    atlas1d_lineplot, xbig, residbig, xerrbig, residerrbig, plottype=3, postscript=postscript, $
;   atlas1d_lineplot, x, resid, xerr, resid_err, plottype=3, postscript=postscript, $
      xtitle=xtitle, ytitle='', xrange=xrange, yrange=yrange, legendtype=0, $
      charsize=charsize_5, position=pos[*,5], $
      ytickname=replicate(' ',10), xminor=3, yminor=3, ystyle=11, /noerase, $
      atlaspsize=0.5;, nfgspsize=0.8, xnfgs=xnfgs, ynfgs=residnfgs, xerrnfgs=xerrnfgs, yerrnfgs=residnfgs_err
    axis, /yaxis, yrange=yrange, ythick=postthick, ytitle='', $
      charthick=postthick, charsize=charsize_5, ysty=1, yminor=3
    im_xyouts_title, ytitle='', charsize=charsize_5, charthick=postthick
    djs_oplot, !x.crange, [0,0], line=0, thick=postthick
    
    legend, textoidl('(f) H\alpha/H\beta Correction'), /left, /top, box=0, charsize=charsize_3, $
      charthick=postthick, clear=keyword_set(postscript)
    legend, textoidl(xstr), /left, /bottom, box=0, charsize=charsize_3, $
      charthick=postthick, clear=keyword_set(postscript)

    im_openclose, postscript=postscript, /close    

; ------------------------------------------------------------
; 4-panel empirical reddening indicators - ATLAS/NFGS
; ------------------------------------------------------------

    psname = 'ebv_empirical_4panel'
    im_openclose, pspath+psname, postscript=postscript, /encapsulated, xsize=8.5, ysize=9.0

    pagemaker, nx=2, ny=2, xspace=0.0, yspace=1.0, width=3.05*[1,1], height=3.05*[1,1], $
      xmargin=[1.2,1.2], ymargin=[0.8,1.1], xpage=8.5, ypage=9.0, position=pos, /normal

    ebvrange = [-0.05,0.85]
    yrange = ebvrange
    
    ytitle = 'E(B-V) from H\alpha/H\beta'

    hgsnrcut = 10.0
    
; ############################################################
; Panel 1:     
; ############################################################

    indx = where((atlasnodust.ebv_hahb_err gt 0.0) and (atlasnodust.ebv_hbhg_err gt 0.0) and $
      (atlasdust.h_gamma[0]/atlasdust.h_gamma[1] gt hgsnrcut),nindx)
;   indx = where((atlasnodust.ebv_hahb_err gt 0.0) and (atlasnodust.ebv_hbhg_err gt 0.0),nindx)
    
    x = atlasnodust[indx].ebv_hbhg
    xerr = atlasnodust[indx].ebv_hbhg_err
    
    y = atlasnodust[indx].ebv_hahb
    yerr = atlasnodust[indx].ebv_hahb_err

    indxnfgs = where((nfgsnodust.ebv_hahb_err gt 0.0) and (nfgsnodust.ebv_hbhg_err gt 0.0) and $
      (nfgsdust.h_gamma[0]/nfgsdust.h_gamma[1] gt hgsnrcut),nindxnfgs)
;   indxnfgs = where((nfgsnodust.ebv_hahb_err gt 0.0) and (nfgsnodust.ebv_hbhg_err gt 0.0),nindxnfgs)
    
    xnfgs = nfgsnodust[indxnfgs].ebv_hbhg
    xerrnfgs = nfgsnodust[indxnfgs].ebv_hbhg_err
    
    ynfgs = nfgsnodust[indxnfgs].ebv_hahb
    yerrnfgs = nfgsnodust[indxnfgs].ebv_hahb_err

    residuals = [x,xnfgs]-[y,ynfgs]
    stats = im_stats(residuals,/verbose)
    
    xtitle = 'E(B-V) from H\beta/H\gamma'
    
;   xrange = [-0.1,1.6]
    xrange = ebvrange

    rcor = r_correlate([x,xnfgs],[y,ynfgs],zd=zd,probd=probd)
    splog, 'Spearman rank test: ', rcor, probd, zd

    atlas1d_lineplot, x, y, xerr, yerr, plottype=1, postscript=postscript, $
      xtitle=xtitle, ytitle=ytitle, xrange=xrange, yrange=yrange, legendtype=0, $
      /left, /top, yfrac=1.2, xfrac=1.3, position=pos[*,0], charsize=charsize_5, $
      xnfgs=xnfgs, ynfgs=ynfgs, xerrnfgs=xerrnfgs, yerrnfgs=yerrnfgs
;   djs_oplot, !x.crange, !y.crange, line=0, thick=(postthick-4)>2
    djs_oplot, [0,!x.crange[1]], [0,!y.crange[1]], line=0, thick=(postthick-2)>2    
;   legend, '(a)', /left, /top, box=0, charsize=charsize_3, charthick=postthick
    legend, textoidl('(a) S/N (H\gamma) > '+string(hgsnrcut,format='(I0)')), $
      /left, /top, box=0, charsize=charsize_3, charthick=postthick

;   label = ['S/N [H\gamma] > 10']
;   legend, textoidl(label), /left, /bottom, box=0, charsize=charsize_3, $
;     charthick=postthick
    
    
; ############################################################
; Panel 2:     
; ############################################################

    indx = where((atlasnodust.ebv_hahb_err gt 0.0) and (atlasdust.b_lum_obs gt -900.0),nindx)

    x = atlasdust[indx].b_lum_obs
    xerr = atlasdust[indx].b_lum_obs_err
    xabs = atlasdust[indx].m_b_obs

    y = atlasnodust[indx].ebv_hahb
    yerr = atlasnodust[indx].ebv_hahb_err

    indxnfgs = where((nfgsnodust.ebv_hahb_err gt 0.0) and (nfgsdust.b_lum_obs gt -900.0),nindxnfgs)

    xnfgs = nfgsdust[indxnfgs].b_lum_obs
    xerrnfgs = nfgsdust[indxnfgs].b_lum_obs_err

    ynfgs = nfgsnodust[indxnfgs].ebv_hahb
    yerrnfgs = nfgsnodust[indxnfgs].ebv_hahb_err
    
    xtitle = 'log L(B) [L(B)'+sunsymbol()+']'
    xrange = LBrange

    atlas1d_lineplot, x, y, xerr, yerr, plottype=1, postscript=postscript, $
      xtitle=xtitle, ytitle='', xrange=xrange, yrange=yrange, legendtype=0, $
      xstyle=11, ymargin=[4,3], xmargin=[8,6], position=pos[*,1], $
      /left, /top, /noerase, charsize=charsize_5, ystyle=11, $
      ytickname=replicate(' ',10), $
      xnfgs=xnfgs, ynfgs=ynfgs, xerrnfgs=xerrnfgs, yerrnfgs=yerrnfgs
    axis, /xaxis, xrange=interpol(xabs,x,!x.crange), xthick=postthick, xstyle=1, $
      charsize=charsize_5, charthick=postthick
    im_xyouts_title, xtitle='M_{B} [mag]', charsize=charsize_5, charthick=postthick
    axis, /yaxis, yrange=yrange, ythick=postthick, ytitle='', $
      charthick=postthick, charsize=charsize_5, ysty=1
    im_xyouts_title, ytitle=ytitle, charsize=charsize_5, charthick=postthick
    legend, '(b)', /left, /top, box=0, charsize=charsize_5, charthick=postthick

; ############################################################
; Panel 3:     
; ############################################################
    
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

    xtitle = 'log L([O II]) [erg s^{-1}]'
    xrange = Loiikewleyrange

    atlas1d_lineplot, x, y, xerr, yerr, plottype=1, postscript=postscript, $
      xtitle=xtitle, ytitle=ytitle, xrange=xrange, yrange=yrange, legendtype=0, $
      charsize=charsize_5, position=pos[*,2], /noerase, xminor=3, $
      /left, /top, xnfgs=xnfgs, ynfgs=ynfgs, xerrnfgs=xerrnfgs, yerrnfgs=yerrnfgs
    legend, '(c)', /left, /top, box=0, charsize=charsize_5, charthick=postthick

; overplot the Kewley relation    

    xaxis = findgen((xrange[1]-xrange[0])/0.01+1)*0.01+xrange[0]
    yaxis = 0.174*xaxis - 6.84
    chop = where(yaxis gt 0.0)

    djs_oplot, xaxis[chop], yaxis[chop], line=2, thick=(postthick-2)>2;, color='dark green'
;   legend, ['KGJ04'], line=0, /right, /top, charsize=charsize, $
;     charthick=postthick, box=0, thick=postthick, clear=postscript

; ############################################################
; Panel 4:     
; ############################################################

    cut = where(atlasnodust.ebv_hahb_err gt 0.0)
    lineratio, atlasnodust[cut], 'OII_3727_EW', '', '', '', x, xerr, $
      dum1, dum2, index=indx, nindex=nindx, snrcut=snrcut, $
      xsyserr=syserr, ysyserr=syserr

    y = atlasnodust[cut[indx]].ebv_hahb
    yerr = atlasnodust[cut[indx]].ebv_hahb_err

    cutnfgs = where(nfgsnodust.ebv_hahb_err gt 0.0)
    lineratio, nfgsnodust[cutnfgs], 'OII_3727_EW', '', '', '', xnfgs, xerrnfgs, $
      dum1, dum2, index=indxnfgs, nindex=nindxnfgs, snrcut=snrcut, $
      xsyserr=nfgssyserr, ysyserr=nfgssyserr

    ynfgs = nfgsnodust[cutnfgs[indxnfgs]].ebv_hahb
    yerrnfgs = nfgsnodust[cutnfgs[indxnfgs]].ebv_hahb_err
    
    xtitle = 'log EW([O II]) ['+angstrom()+']'
    xrange = ewoiirange

    atlas1d_lineplot, x, y, xerr, yerr, plottype=1, postscript=postscript, $
      xtitle=xtitle, ytitle='', xrange=xrange, yrange=yrange, legendtype=0, $
      ystyle=11, xmargin=[8,6], ymargin=[4,3], /noerase, /right, /top, $
      charsize=charsize_5, ytickname=replicate(' ',10), position=pos[*,3], $
      xnfgs=xnfgs, ynfgs=ynfgs, xerrnfgs=xerrnfgs, yerrnfgs=yerrnfgs
    axis, /yaxis, yrange=yrange, ythick=postthick, ytitle='', $
      charthick=postthick, charsize=charsize_5, ysty=1
    im_xyouts_title, ytitle=ytitle, charsize=charsize_5, charthick=postthick
    legend, '(d)', /left, /top, box=0, charsize=charsize_3, charthick=postthick

    im_openclose, postscript=postscript, /close    
    

; ------------------------------------------------------------
; ([O III]/Hb)/([N II]/Ha) vs E(Hb-Ha) [Integrated + SDSS]
; ------------------------------------------------------------

    psname = 'oiiinii_vs_ehbha_2panel'
    im_openclose, pspath+psname, postscript=postscript, /encapsulated, xsize=8.5, ysize=5.15

    pagemaker, nx=2, ny=1, height=3.15, width=3.15, xmargin=[1.1,1.1], $
      ymargin=[0.9,1.1], xspace=0, yspace=0, xpage=8.5, ypage=5.15, $
      position=pos, /normal

    xtitle = 'log {([O III]/H\beta)/([N II]/H\alpha)}'
    ytitle = 'E(H\beta-H\alpha) [mag]'

    xrange = reverse(oiiiniirange)
    yrange = ehbharange

; Integrated    
    
    indx = where((atlasdust.zstrong_oiiinii gt -900) and (atlasnodust.ehbha_err gt 0),nindx)

    x = atlasdust[indx].zstrong_oiiinii
    xerr = atlasdust[indx].zstrong_oiiinii_err

    y = atlasnodust[indx].ehbha
    yerr = atlasnodust[indx].ehbha_err
    yebv = y/(k_lambda(4861,/odonnell)-k_lambda(6563,/odonnell))

    OHgood = where(atlasdust[indx].zstrong_12oh_oiiinii_pettini gt -900)
    xgood = x[OHgood]
    xOHgood = atlasdust[indx[OHgood]].zstrong_12oh_oiiinii_pettini

    indxnfgs = where((nfgsdust.zstrong_oiiinii gt -900) and (nfgsnodust.ehbha_err gt 0.0),nindxnfgs)

    xnfgs = nfgsdust[indxnfgs].zstrong_oiiinii
    xerrnfgs = nfgsdust[indxnfgs].zstrong_oiiinii_err

    ynfgs = nfgsnodust[indxnfgs].ehbha
    yerrnfgs = nfgsnodust[indxnfgs].ehbha_err

    xbig = [x,xnfgs]
    ybig = [y,ynfgs]

    atlas1d_lineplot, x, y, xerr, yerr, plottype=3, postscript=postscript, $n
      xtitle=xtitle, ytitle=ytitle, xrange=xrange, yrange=yrange, legendtype=0, $
      xstyle=11, ymargin=[4,3], /left, /top, position=pos[*,0], $
      charsize=charsize_3, /xreverse, /errorleft, xfrac=8.0, xnfgs=xnfgs, ynfgs=ynfgs, $
      xerrnfgs=xerrnfgs, yerrnfgs=yerrnfgs
    axis, /xaxis, xrange=interpol(xOHgood,xgood,!x.crange), xthick=postthick, xstyle=1, $
      charsize=charsize_3, charthick=postthick
    im_xyouts_title, xtitle='12 + log (O/H)', charsize=charsize_3, charthick=postthick

; SDSS    

    indx = where((sdssdust.zstrong_oiiinii gt -900) and (sdssnodust.ehbha_err gt 0.0),nindx)

    x = sdssdust[indx].zstrong_oiiinii
    xerr = sdssdust[indx].zstrong_oiiinii_err

    y = sdssnodust[indx].ehbha
    yerr = sdssnodust[indx].ehbha_err
    yebv = y/(k_lambda(4861,/odonnell)-k_lambda(6563,/odonnell))

    OHgood = where(sdssdust[indx].zstrong_12oh_oiiinii_pettini gt -900)
    xgood = x[OHgood]
    xOHgood = sdssdust[indx[OHgood]].zstrong_12oh_oiiinii_pettini

    sdss_lineplot, x, y, xerr, yerr, plottype=3, postscript=postscript, $
      xtitle=xtitle, ytitle='', xrange=xrange, yrange=yrange, legendtype=0, $
      xstyle=11, ymargin=[4,3], ystyle=11, xmargin=[8,6], /left, /top, /xreverse, $
      /errorleft, xfrac=8.0, position=pos[*,1], $
      ytickname=replicate(' ',10), /noerase, charsize=charsize_3
    axis, /xaxis, xrange=interpol(xOHgood,xgood,!x.crange), xthick=postthick, xstyle=1, $
      charsize=charsize_3, charthick=postthick
    im_xyouts_title, xtitle='12 + log (O/H)', charsize=charsize_3, charthick=postthick
    axis, /yaxis, yrange=interpol(yebv,y,!y.crange), ythick=postthick, $
      charthick=postthick, charsize=charsize_3, ystyle=1
    im_xyouts_title, ytitle='E(B-V) [mag]', charsize=charsize_3, charthick=postthick

;   sdss_lineplot, sdssdust[indx].zstrong_obj, x, x*0.0, y*0.0, $
;     xrange=[0,0.25], yrange=xrange, xtitle='z/z_{max}', ytitle=xtitle, /yreverse
;   sdss_lineplot, sdssdust[indx].zstrong_obj/sdssdust[indx].zzmax, x, x*0.0, y*0.0, $
;     xrange=[0,1.1], yrange=xrange, xtitle='z/z_{max}', ytitle=xtitle, /yreverse

    im_openclose, postscript=postscript, /close    

; ------------------------------------------------------------
; L(B) vs various SFR(Hb) calibrations
; ------------------------------------------------------------
    
    psname = 'lb_vs_sfr_hb_3panel'
    im_openclose, pspath+psname, postscript=postscript, /encapsulated, xsize=7.0, ysize=9.0

    pagemaker, nx=1, ny=3, height=[2.5,2.5,2.5], width=5.5, xmargin=[1.1,0.4], $
      ymargin=[0.4,1.1], xspace=0, yspace=0, xpage=7.0, ypage=9.0, $
      position=pos, /normal

    xrange = lbrange
    yrange = residrange_hbsfrs

    xtitle = 'log L(B) [L(B)'+sunsymbol()+']'
    ytitle = 'log \psi(H\beta)/\psi(H\alpha)'
    
; ##########################################################
; Panel 1
; ##########################################################

; Atlas    
    
    indx = where((atlas_sfrs.sfr_ha gt -900) and (atlas_sfrs.sfr_hb_uncor gt -900),nindx)

    x = atlasdust[indx].b_lum_obs
    xerr = atlasdust[indx].b_lum_obs_err
    
    y1 = atlas_sfrs[indx].sfr_ha
    y1err = atlas_sfrs[indx].sfr_ha_err
    
    y2 = atlas_sfrs[indx].sfr_hb_uncor
    y2err = atlas_sfrs[indx].sfr_hb_uncor_err

    resid = y2-y1
    resid_err = sqrt(y1err^2+y2err^2)

; NFGS    
    
    indxnfgs = where((nfgs_sfrs.sfr_ha gt -900) and (nfgs_sfrs.sfr_hb_uncor gt -900),nindxnfgs)

    xnfgs = nfgsdust[indxnfgs].b_lum_obs
    xerrnfgs = nfgsdust[indxnfgs].b_lum_obs_err
    
    y1nfgs = nfgs_sfrs[indxnfgs].sfr_ha
    y1errnfgs = nfgs_sfrs[indxnfgs].sfr_ha_err
    
    y2nfgs = nfgs_sfrs[indxnfgs].sfr_hb_uncor
    y2errnfgs = nfgs_sfrs[indxnfgs].sfr_hb_uncor_err

    residnfgs = y2nfgs-y1nfgs
    residnfgs_err = sqrt(y1errnfgs^2+y2errnfgs^2)

; combine the two samples

    xbig = [x,xnfgs]
    xerrbig = [xerr,xerrnfgs]
    residbig = [resid,residnfgs]
    residerrbig = [resid_err,residnfgs_err]
    
; make the plot    

    stats = im_stats(residbig,/verbose,/baremin)
    xstr = strtrim(string(stats.median_rej,format='(F12.2)'),2)+' ('+$
      strtrim(string(stats.mean_rej,format='(F12.2)'),2)+$
      '\pm'+strtrim(string(stats.sigma_rej,format='(F12.2)'),2)+')'
;   xstr = strtrim(string(stats.median,format='(F12.2)'),2)+'\pm'+strtrim(string(stats.sig68mean,format='(F12.2)'),2)
    
    atlas1d_lineplot, xbig, residbig, xerrbig, residerrbig, plottype=3, postscript=postscript, $
      xtitle='', ytitle='', xrange=xrange, yrange=yrange, legendtype=0, $
      charsize=charsize_7, xtickname=replicate(' ',10), position=pos[*,0], $
      xminor=3, yminor=3, atlaspsize=0.8, ytickinterval=1.0;, nfgspsize=0.8, $
;     xnfgs=xnfgs, ynfgs=residnfgs, xerrnfgs=xerrnfgs, yerrnfgs=residnfgs_err
    djs_oplot, !x.crange, [0,0], line=2, thick=postthick
    
    legend, textoidl('(a) Observed L(H\beta)'), /left, /top, box=0, charsize=charsize_3, charthick=postthick
    legend, textoidl(xstr), /right, /top, box=0, charsize=charsize_3, charthick=postthick

; ##########################################################
; Panel 2
; ##########################################################

; Atlas    
    
    indx = where((atlas_sfrs.sfr_ha gt -900) and (atlas_sfrs.sfr_hb_hbhg gt -900),nindx)

    x = atlasdust[indx].b_lum_obs
    xerr = atlasdust[indx].b_lum_obs_err
    
    y1 = atlas_sfrs[indx].sfr_ha
    y1err = atlas_sfrs[indx].sfr_ha_err
    
    y2 = atlas_sfrs[indx].sfr_hb_hbhg
    y2err = atlas_sfrs[indx].sfr_hb_hbhg_err

    resid = y2-y1
    resid_err = sqrt(y1err^2+y2err^2)

; NFGS    
    
    indxnfgs = where((nfgs_sfrs.sfr_ha gt -900) and (nfgs_sfrs.sfr_hb_hbhg gt -900),nindxnfgs)

    xnfgs = nfgsdust[indxnfgs].b_lum_obs
    xerrnfgs = nfgsdust[indxnfgs].b_lum_obs_err
    
    y1nfgs = nfgs_sfrs[indxnfgs].sfr_ha
    y1errnfgs = nfgs_sfrs[indxnfgs].sfr_ha_err
    
    y2nfgs = nfgs_sfrs[indxnfgs].sfr_hb_hbhg
    y2errnfgs = nfgs_sfrs[indxnfgs].sfr_hb_hbhg_err

    residnfgs = y2nfgs-y1nfgs
    residnfgs_err = sqrt(y1errnfgs^2+y2errnfgs^2)

; combine the two samples

    xbig = [x,xnfgs]
    xerrbig = [xerr,xerrnfgs]
    residbig = [resid,residnfgs]
    residerrbig = [resid_err,residnfgs_err]
    
; make the plot    

    stats = im_stats(residbig,/verbose,/baremin,/no_head)
    xstr = strtrim(string(stats.median_rej,format='(F12.2)'),2)+' ('+$
      strtrim(string(stats.mean_rej,format='(F12.2)'),2)+$
      '\pm'+strtrim(string(stats.sigma_rej,format='(F12.2)'),2)+')'
;   xstr = strtrim(string(stats.median,format='(F12.2)'),2)+'\pm'+strtrim(string(stats.sig68mean,format='(F12.2)'),2)
    
    atlas1d_lineplot, xbig, residbig, xerrbig, residerrbig, plottype=3, postscript=postscript, $
      xtitle='', ytitle=ytitle, xrange=xrange, yrange=yrange, legendtype=0, $
      charsize=charsize_7, xtickname=replicate(' ',10), position=pos[*,1], $
      xminor=3, yminor=3, /noerase, atlaspsize=0.8, ytickinterval=1.0;, nfgspsize=0.8, $
;     xnfgs=xnfgs, ynfgs=residnfgs, xerrnfgs=xerrnfgs, yerrnfgs=residnfgs_err
    djs_oplot, !x.crange, [0,0], line=2, thick=postthick
    
    legend, textoidl('(b) L(H\beta) corrected using H\beta/H\gamma'), /left, $
      /top, box=0, charsize=charsize_3, charthick=postthick
    legend, textoidl(xstr), /right, /top, box=0, charsize=charsize_3, charthick=postthick

; ##########################################################
; Panel 3
; ##########################################################

; Atlas    
    
    indx = where((atlas_sfrs.sfr_ha gt -900) and (atlas_sfrs.sfr_hb_best gt -900),nindx)

    x = atlasdust[indx].b_lum_obs
    xerr = atlasdust[indx].b_lum_obs_err
    
    y1 = atlas_sfrs[indx].sfr_ha
    y1err = atlas_sfrs[indx].sfr_ha_err
    
    y2 = atlas_sfrs[indx].sfr_hb_best
    y2err = atlas_sfrs[indx].sfr_hb_best_err

    resid = y2-y1
    resid_err = sqrt(y1err^2+y2err^2)

; NFGS    
    
    indxnfgs = where((nfgs_sfrs.sfr_ha gt -900) and (nfgs_sfrs.sfr_hb_best gt -900),nindxnfgs)

    xnfgs = nfgsdust[indxnfgs].b_lum_obs
    xerrnfgs = nfgsdust[indxnfgs].b_lum_obs_err
    
    y1nfgs = nfgs_sfrs[indxnfgs].sfr_ha
    y1errnfgs = nfgs_sfrs[indxnfgs].sfr_ha_err
    
    y2nfgs = nfgs_sfrs[indxnfgs].sfr_hb_best
    y2errnfgs = nfgs_sfrs[indxnfgs].sfr_hb_best_err

    residnfgs = y2nfgs-y1nfgs
    residnfgs_err = sqrt(y1errnfgs^2+y2errnfgs^2)

; combine the two samples

    xbig = [x,xnfgs]
    xerrbig = [xerr,xerrnfgs]
    residbig = [resid,residnfgs]
    residerrbig = [resid_err,residnfgs_err]
    
; make the plot    

    stats = im_stats(residbig,/verbose,/baremin,/no_head)
    xstr = strtrim(string(stats.median_rej,format='(F12.2)'),2)+' ('+$
      strtrim(string(stats.mean_rej,format='(F12.2)'),2)+$
      '\pm'+strtrim(string(stats.sigma_rej,format='(F12.2)'),2)+')'
;   xstr = strtrim(string(stats.median,format='(F12.2)'),2)+'\pm'+strtrim(string(stats.sig68mean,format='(F12.2)'),2)

;   yrange = [-1,1]
    
    atlas1d_lineplot, xbig, residbig, xerrbig, residerrbig, plottype=3, postscript=postscript, $
      xtitle=xtitle, ytitle='', xrange=xrange, yrange=yrange, legendtype=0, $
      charsize=charsize_7, position=pos[*,2], $
      xminor=3, yminor=3, /noerase, atlaspsize=0.8, ytickinterval=1.0;, nfgspsize=0.8, $
;     xnfgs=xnfgs, ynfgs=residnfgs, xerrnfgs=xerrnfgs, yerrnfgs=residnfgs_err
    djs_oplot, !x.crange, [0,0], line=2, thick=postthick
    
    legend, '(c) This Paper', /left, /top, box=0, charsize=charsize_3, charthick=postthick
    legend, textoidl(xstr), /right, /top, box=0, charsize=charsize_3, charthick=postthick

    im_openclose, postscript=postscript, /close    

; ------------------------------------------------------------
; L(Ha)/L(IR) vs A(Ha)
; ------------------------------------------------------------

    psname = 'LHa_LIR_vs_AHa'
    im_openclose, pspath+psname, postscript=postscript, /encapsulated

; Atlas    
    
    indx = where((atlasnodust.ir_flux gt -900) and (atlasnodust.ebv_hahb_err gt 0.0),nindx)

    y = atlasnodust[indx].ebv_hahb*k_lambda(6563.0,/odonnell)
    yerr = atlasnodust[indx].ebv_hahb_err*k_lambda(6563.0,/odonnell)
    
    lir = atlasdust[indx].ir_flux ; [erg/s/cm2]
    lir_err = sqrt(atlasdust[indx].ir_flux_err^2 + (atlasdust[indx].ir_flux*0.15)^2)

    ha = atlasdust[indx].h_alpha[0]
    ha_err = atlasdust[indx].h_alpha[1]

    x = alog10(ha/lir)
    xerr = im_compute_error(ha,ha_err,lir,lir_err,/log)

; NFGS    
    
    indxnfgs = where((nfgsdust.ir_flux gt -900) and (nfgsnodust.ebv_hahb_err gt 0.0),nindxnfgs)

    ynfgs = nfgsnodust[indxnfgs].ebv_hahb*k_lambda(6563.0,/odonnell)
    yerrnfgs = nfgsnodust[indxnfgs].ebv_hahb_err*k_lambda(6563.0,/odonnell)
    
    lirnfgs = nfgsdust[indxnfgs].ir_flux
    lirnfgs_err = sqrt(nfgsdust[indxnfgs].ir_flux_err^2 + (nfgsdust[indxnfgs].ir_flux*0.15)^2)
    
    hanfgs = nfgsdust[indxnfgs].h_alpha[0]
    hanfgs_err = nfgsdust[indxnfgs].h_alpha[1]

    xnfgs = alog10(hanfgs/lirnfgs)
    xerrnfgs = im_compute_error(hanfgs,hanfgs_err,lirnfgs,lirnfgs_err,/log)

    xrange = haLIRrange
    yrange = AHarange

    xtitle = 'log [L(H\alpha)_{obs}/L(IR)]'
    ytitle = 'A(H\alpha)'

;   atlas1d_lineplot, x, y, xerr, yerr, plottype=3, postscript=postscript, $
;     xtitle=xtitle, ytitle=ytitle, xrange=xrange, yrange=yrange, legendtype=0, $
;     xnfgs=xnfgs, ynfgs=ynfgs, xerrnfgs=xerrnfgs, yerrnfgs=yerrnfgs, yfrac=3.4

    atlas1d_lineplot, x, y, xerr, yerr, postscript=postscript, /nodata, $
      xtitle=xtitle, ytitle=ytitle, xrange=xrange, yrange=yrange

; overlay the theoretical flux ratios; This will compute the
; attenuation at 1600 A for a galaxy with a measured F(IR)/F(1600)
; flux ratio of 10.0, an age of 10.0, a metallicity of 0.0 (solar), an
; a_d (fraction of ionizing photons absorbed by dust) of 0.5, and
; constant star formation (fit_info_c).  The output of the above
; command would be:

;   read_many_fit, frmodels, sf_type='c' ; burst SES models
    read_many_fit, frmodels, sf_type='c' ; constant SES models

    fr = findgen((!x.crange[1]-!x.crange[0])/0.01+1)*0.01+!x.crange[0]
    nfr = n_elements(fr)

; compute the attenuation at Ha assuming various FR parameters; assume
; constant star formation for ages ranging from 1-100 Myr; vary a_d
; from 0.0 to 0.5; and vary the metallicity over [-0.4,0.0,0.4] =
; (0.4, 1.0, and 2.5 times solar).

    ads = [0.0,0.25]
;   ads = [0.01,0.25,0.5]
    nads = n_elements(ads)

    mets = [-0.4,0.0,0.4]
;   metsmax = 0.4 & metsmin = -0.4 & dmets = 0.1
;   mets = findgen((metsmax-metsmin)/dmets+1)*dmets+metsmin
    nmets = n_elements(mets)

    ages = [10.0,50.0,100.0]
;   agesmax = 100.0 & agesmin = 10.0 & dages = 10.0
;   ages = findgen((agesmax-agesmin)/dages+1)*dages+agesmin
    nages = n_elements(ages)
    
    AHa = fltarr(nads,nmets,nages,nfr)

    for iads = 0L, nads-1L do begin
       for imets = 0L, nmets-1L do begin
          for iages = 0L, nages-1L do begin
             for ifr = 0L, nfr-1L do begin
                AHa[iads,imets,iages,ifr] = fr2atten(10^(-fr[ifr]),'Ha',6563,$
                  ages[iages],mets[imets],ads[iads],frmodels,/silent)
             endfor
          endfor
       endfor
    endfor

    nAHa = nads*nmets*nages
    AHa = reform(aha,nAHa,nfr)

    bigfr = rebin(reform(fr,1,nfr),nAHa,nfr)

    left = min(aha[*,0],leftindx)
    right = max(aha[*,0],rightindx)

    polyx = [reform(bigfr[leftindx,*]),reverse(reform(bigfr[rightindx,*]))]
    polyy = [reform(aha[leftindx,*]),reverse(reform(aha[rightindx,*]))]

;   oplot, bigfr[leftindx,*], aha[leftindx,*], line=0, thick=2
;   oplot, bigfr[rightindx,*], aha[rightindx,*], line=0, thick=2

;   polyx = [bigfr[leftindx,*],bigfr[rightindx,*]]
;   polyy = [aha[leftindx,*],aha[rightindx,*]]
;   polyx = [reform(bigfr[leftindx,*]),reform(bigfr[rightindx,*])]
;   polyy = [reform(aha[leftindx,*]),reform(aha[rightindx,*])]

;   loadct, 0, /silent
;   device, decompose=0
;   polyfill, polyx, polyy, color=125, noclip=0, linestyle=3, /line_fill, $
;     spacing=0.01, clip=[!x.crange[0],!y.crange[0],!x.crange[1],!y.crange[1]]
;   device, /decompose

    polyfill, polyx, polyy, color=djs_icolor('grey'), noclip=0, $
      clip=[!x.crange[0],!y.crange[0],!x.crange[1],!y.crange[1]]

;   polyfill, polyx, polyy, color=djs_icolor('grey'), noclip=0, linestyle=3, /line_fill, $
;     spacing=0.025, orientation=45, clip=[!x.crange[0],!y.crange[0],!x.crange[1],!y.crange[1]]

; now overplot the data

    atlas1d_lineplot, x, y, xerr, yerr, plottype=3, postscript=postscript, $
      legendtype=0, yfrac=3.4, xnfgs=xnfgs, ynfgs=ynfgs, xerrnfgs=xerrnfgs, $
      yerrnfgs=yerrnfgs, /overplot
    
    im_openclose, postscript=postscript, /close    

; ------------------------------------------------------------
; L(IR)/L(UV) vs A(Ha)
; ------------------------------------------------------------
    
    psname = 'LIR_LUV_vs_AHa'
    im_openclose, pspath+psname, postscript=postscript, /encapsulated

; Atlas    
    
    indx = where((atlasnodust.ir_flux gt -900) and (atlasnodust.fuv_flux gt -900) and $
      (atlasnodust.ebv_hahb_err gt 0.0),nindx)

    lir = atlasdust[indx].ir_flux ; [erg/s/cm2]
    lir_err = sqrt(atlasdust[indx].ir_flux_err^2 + (atlasdust[indx].ir_flux*0.15)^2)

    luv = atlasdust[indx].fuv_flux
    luv_err = atlasdust[indx].fuv_flux_err

    x = alog10(lir/luv)
    xerr = im_compute_error(lir,lir_err,luv,luv_err,/log)
    
    y = atlasnodust[indx].ebv_hahb*k_lambda(6563.0,/odonnell)
    yerr = atlasnodust[indx].ebv_hahb_err*k_lambda(6563.0,/odonnell)
    
; NFGS    
    
    indxnfgs = where((nfgsdust.ir_flux gt -900) and (nfgsdust.fuv_flux gt -900) and $
      (nfgsnodust.ebv_hahb_err gt 0.0),nindxnfgs)

    lirnfgs = nfgsdust[indxnfgs].ir_flux
    lirnfgs_err = sqrt(nfgsdust[indxnfgs].ir_flux_err^2 + (nfgsdust[indxnfgs].ir_flux*0.15)^2)

    luvnfgs = nfgsdust[indxnfgs].fuv_flux
    luvnfgs_err = nfgsdust[indxnfgs].fuv_flux_err

    xnfgs = alog10(lirnfgs/luvnfgs)
    xerrnfgs = im_compute_error(lirnfgs,lirnfgs_err,luvnfgs,luvnfgs_err,/log)
    
    ynfgs = nfgsnodust[indxnfgs].ebv_hahb*k_lambda(6563.0,/odonnell)
    yerrnfgs = nfgsnodust[indxnfgs].ebv_hahb_err*k_lambda(6563.0,/odonnell)
    
    xrange = LIRLUVrange
    yrange = AHarange

    xtitle = 'log [L(IR)/L(UV)]'
    ytitle = 'A(H\alpha)'

    atlas1d_lineplot, x, y, xerr, yerr, plottype=3, postscript=postscript, $
      xtitle=xtitle, ytitle=ytitle, xrange=xrange, yrange=yrange, legendtype=0, $
      xnfgs=xnfgs, ynfgs=ynfgs, xerrnfgs=xerrnfgs, yerrnfgs=yerrnfgs

    im_openclose, postscript=postscript, /close    

; ------------------------------------------------------------
; E(B-V) [Ha/Hb] vs E(B-V) [Hb/Hg]
; ------------------------------------------------------------

    psname = 'sdss_ebv_hahb_vs_ebv_hbhg'
    im_openclose, pspath+psname, postscript=postscript, /encapsulated
    
    indx = where((sdssnodust.ebv_hahb_err gt 0.0) and (sdssnodust.ebv_hbhg_err gt 0.0),nindx)
    
    x = sdssnodust[indx].ebv_hahb
    xerr = sdssnodust[indx].ebv_hahb_err

    y = sdssnodust[indx].ebv_hbhg
    yerr = sdssnodust[indx].ebv_hbhg_err

    indxatlas = where((atlasnodust.ebv_hahb_err gt 0.0) and (atlasnodust.ebv_hbhg_err gt 0.0),nindxatlas)
    
    xatlas = atlasnodust[indxatlas].ebv_hahb
    xerratlas = atlasnodust[indxatlas].ebv_hahb_err

    yatlas = atlasnodust[indxatlas].ebv_hbhg
    yerratlas = atlasnodust[indxatlas].ebv_hbhg_err
    
    indxnfgs = where((nfgsnodust.ebv_hahb_err gt 0.0) and (nfgsnodust.ebv_hbhg_err gt 0.0),nindxnfgs)
    
    xnfgs = nfgsnodust[indxnfgs].ebv_hahb
    xerrnfgs = nfgsnodust[indxnfgs].ebv_hahb_err

    ynfgs = nfgsnodust[indxnfgs].ebv_hbhg
    yerrnfgs = nfgsnodust[indxnfgs].ebv_hbhg_err
    
    xtitle = 'E(B-V) [H\alpha/H\beta]'
    ytitle = 'E(B-V) [H\beta/H\gamma]'
    
    xrange = [0,1.5]
    yrange = xrange

    sdss_lineplot, x, y, xerr, yerr, plottype=3, postscript=postscript, $
      xtitle=xtitle, ytitle=ytitle, xrange=xrange, yrange=yrange, legendtype=0, $
      /left, /top, yfrac=1.3, $
      xatlas=xatlas, yatlas=yatlas, xerratlas=xerratlas, yerratlas=yerratlas, $
      xnfgs=xnfgs, ynfgs=ynfgs, xerrnfgs=xerrnfgs, yerrnfgs=yerrnfgs
    djs_oplot, !x.crange, !y.crange, line=0, thick=postthick

    im_openclose, postscript=postscript, /close

; ---------------------------------------------------------------------------
; Hg Measured vs Hg Predicted
; ---------------------------------------------------------------------------

    psname = 'sdss_hg_measured_vs_hg_predicted'
    im_openclose, pspath+psname, postscript=postscript, /encapsulated

    lineratio, sdssdust, '', '', 'H_GAMMA', 'H_GAMMA_DUST_PREDICT', $
      dum1, dum2, y, yerr, index=indx, nindex=nindx, snrcut=0.0

    x = alog10(sdssdust[indx].h_gamma[0]/sdssdust[indx].h_gamma[1])
    xerr = x*0.0
    
    lineratio, atlasdust, '', '', 'H_GAMMA', 'H_GAMMA_DUST_PREDICT', $
      dum1, dum2, yatlas, yerratlas, index=indxatlas, nindex=nindxatlas, snrcut=0.0

    xatlas = alog10(atlasdust[indxatlas].h_gamma[0]/atlasdust[indxatlas].h_gamma[1])
    xerratlas = xatlas*0.0
    
    lineratio, nfgsdust, '', '', 'H_GAMMA', 'H_GAMMA_DUST_PREDICT', $
      dum1, dum2, ynfgs, yerrnfgs, index=indxnfgs, nindex=nindxnfgs, snrcut=0.0

    xnfgs = alog10(nfgsdust[indxnfgs].h_gamma[0]/nfgsdust[indxnfgs].h_gamma[1])
    xerrnfgs = xnfgs*0.0
       
    stats = im_stats(y,/verbose)

    xrange = alog10([1.0,300])
    yrange = [-0.7,0.7]

    xtitle = 'log (S/N) [H\gamma]'
    ytitle = 'log (H\gamma_{measured}/H\gamma_{predicted})'
    
    sdss_lineplot, x, y, xerr, yerr, plottype=3, postscript=postscript, $
      xtitle=xtitle, ytitle=ytitle, xrange=xrange, yrange=yrange, legendtype=0, $
      /right, /top, $
      xatlas=xatlas, yatlas=yatlas, xerratlas=xerratlas, yerratlas=yerratlas, $
      xnfgs=xnfgs, ynfgs=ynfgs, xerrnfgs=xerrnfgs, yerrnfgs=yerrnfgs
    djs_oplot, !x.crange, [0.0,0.0], line=0, thick=postthick

    im_openclose, postscript=postscript, /close

; ------------------------------------------------------------
; 3-panel - EW(Hb) vs Hb/Ha - Integrated
; ------------------------------------------------------------
    
    psname = 'ewhb_vs_hbha_multipanel'
    im_openclose, pspath+psname, postscript=postscript, /encapsulated, xsize=8.5, ysize=8.0

    pagemaker, nx=1, ny=3, height=[2.0,2.0,2.0], width=5.5, xmargin=[2.0,1.0], $
      ymargin=[0.9,1.1], xspace=0, yspace=0, xpage=8.5, ypage=8.0, $
      position=pos, /normal

    indx = where((atlasnodust.ebv_hahb_err gt 0.0) and (atlasdust.h_beta_ew_uncor[0] ge ewcut),nindx)

    x = alog10(atlasdust[indx].h_beta_ew[0])
    xerr = atlasdust[indx].h_beta_ew[1]/atlasdust[indx].h_beta_ew[0]/alog(10.0)
    
    indxnfgs = where((nfgsnodust.ebv_hahb_err gt 0.0) and (nfgsdust.h_beta_ew_uncor[0] ge ewcut),nindxnfgs)

    xnfgs = alog10(nfgsnodust[indxnfgs].h_beta_ew[0])
    xerrnfgs = nfgsnodust[indxnfgs].h_beta_ew[1]/nfgsnodust[indxnfgs].h_beta_ew[0]/alog(10.0)
    
    xrange = ewhbrange
    yrange = hbharange

    xtitle = 'log [EW(H\beta)] ['+angstrom()+']'
    ytitle = 'log (H\beta/H\alpha)'

; ##########################################################
; Panel 1: No absorption correction, A(Hb)=0, A(Ha)=0
; ##########################################################

    hb = atlasdust[indx].h_beta_uncor[0]
    hb_err = atlasdust[indx].h_beta_uncor[1]
    
    ha = atlasdust[indx].h_alpha[0]
    ha_err = atlasdust[indx].h_alpha[1]
    
    y = alog10(hb/ha)
    yerr = im_compute_error(hb,hb_err,ha,ha_err,/log)

    hbnfgs = nfgsdust[indxnfgs].h_beta_uncor[0]
    hbnfgs_err = nfgsdust[indxnfgs].h_beta_uncor[1]
    
    hanfgs = nfgsdust[indxnfgs].h_alpha[0]
    hanfgs_err = nfgsdust[indxnfgs].h_alpha[1]
    
    ynfgs = alog10(hbnfgs/hanfgs)
    yerrnfgs = im_compute_error(hbnfgs,hbnfgs_err,hanfgs,hanfgs_err,/log)

    stats = im_stats([y,ynfgs],/verbose)
    xstr = strtrim(string(stats.median_rej,format='(F12.2)'),2)+' ('+$
      strtrim(string(stats.mean_rej,format='(F12.2)'),2)+$
      '\pm'+strtrim(string(stats.sigma_rej,format='(F12.2)'),2)+')'
;   xstr = strtrim(string(stats.median,format='(F12.2)'),2)+'\pm'+strtrim(string(stats.sig68mean,format='(F12.2)'),2)

    ytitle = 'log (H\beta/H\alpha)_{obs}'

    atlas1d_lineplot, x, y, xerr, yerr, plottype=3, postscript=postscript, $
      xtitle='', ytitle=ytitle, xrange=xrange, yrange=yrange, legendtype=0, $
      charsize=charsize_3, xtickname=replicate(' ',10), position=pos[*,0], $
      xnfgs=xnfgs, ynfgs=ynfgs, xerrnfgs=xerrnfgs, yerrnfgs=yerrnfgs
    djs_oplot, !x.crange, alog10(1.0/HaHb)*[1,1], line=2, thick=postthick

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
    
    y = alog10(hb/ha)
    yerr = im_compute_error(hb,hb_err,ha,ha_err,/log)

    hbnfgs = nfgsdust[indxnfgs].h_beta[0]
    hbnfgs_err = nfgsdust[indxnfgs].h_beta[1]
    
    hanfgs = nfgsdust[indxnfgs].h_alpha[0]
    hanfgs_err = nfgsdust[indxnfgs].h_alpha[1]
    
    ynfgs = alog10(hbnfgs/hanfgs)
    yerrnfgs = im_compute_error(hbnfgs,hbnfgs_err,hanfgs,hanfgs_err,/log)

    stats = im_stats([y,ynfgs],/verbose)
    xstr = strtrim(string(stats.median_rej,format='(F12.2)'),2)+' ('+$
      strtrim(string(stats.mean_rej,format='(F12.2)'),2)+$
      '\pm'+strtrim(string(stats.sigma_rej,format='(F12.2)'),2)+')'
;   xstr = strtrim(string(stats.median,format='(F12.2)'),2)+'\pm'+strtrim(string(stats.sig68mean,format='(F12.2)'),2)

    ytitle = 'log (H\beta/H\alpha)_{obs}'

    atlas1d_lineplot, x, y, xerr, yerr, plottype=3, postscript=postscript, $
      xtitle='', ytitle=ytitle, xrange=xrange, yrange=yrange, legendtype=0, $
      charsize=charsize_3, xtickname=replicate(' ',10), $
      position=pos[*,1], /noerase, xnfgs=xnfgs, ynfgs=ynfgs, xerrnfgs=xerrnfgs, yerrnfgs=yerrnfgs
    djs_oplot, !x.crange, alog10(1.0/HaHb)*[1,1], line=2, thick=postthick

    legend, '(b)', /left, /top, box=0, charsize=charsize_3, charthick=postthick
    legend, textoidl(xstr), /right, /top, box=0, charsize=charsize_3, charthick=postthick

    label = ['Stellar Absorption Corrected']
    legend, textoidl(label), /left, /bottom, box=0, charsize=charsize_3, $
      charthick=postthick

; ############################################################
; Panel 3: Absorption-Corrected, A(Hb)=0, Individual A(Ha) 
; ############################################################

    hb = atlasdust[indx].h_beta[0];*10^(0.4*Ahb)
    hb_err = atlasdust[indx].h_beta[1];*10^(0.4*Ahb)+atlasdust[indx].h_beta[0]*alog(10.0)*10^(0.4*Ahb)*Ahb_err

    ha = atlasnodust[indx].h_alpha[0]
    ha_err = atlasnodust[indx].h_alpha[1]

    y = alog10(hb/ha)
    yerr = im_compute_error(hb,hb_err,ha,ha_err,/log)

    hbnfgs = nfgsdust[indxnfgs].h_beta[0];*10^(0.4*Ahb)
    hbnfgs_err = nfgsdust[indxnfgs].h_beta[1];*10^(0.4*Ahb)+nfgsdust[indxnfgs].h_beta[0]*alog(10.0)*10^(0.4*Ahb)*Ahb_err
    hanfgs = nfgsnodust[indxnfgs].h_alpha[0]
    hanfgs_err = nfgsnodust[indxnfgs].h_alpha[1]

    ynfgs = alog10(hbnfgs/hanfgs)
    yerrnfgs = im_compute_error(hbnfgs,hbnfgs_err,hanfgs,hanfgs_err,/log)
       
    stats = im_stats([y,ynfgs],/verbose)
    xstr = strtrim(string(stats.median_rej,format='(F12.2)'),2)+' ('+$
      strtrim(string(stats.mean_rej,format='(F12.2)'),2)+$
      '\pm'+strtrim(string(stats.sigma_rej,format='(F12.2)'),2)+')'
;   xstr = strtrim(string(stats.median,format='(F12.2)'),2)+'\pm'+strtrim(string(stats.sig68mean,format='(F12.2)'),2)
        
    ytitle = 'log (H\beta_{obs}/H\alpha_{cor})'

    atlas1d_lineplot, x, y, xerr, yerr, plottype=3, postscript=postscript, $
      xtitle=xtitle, ytitle=ytitle, xrange=xrange, yrange=yrange, legendtype=0, $
      charsize=charsize_3, position=pos[*,2], /noerase, $
      xnfgs=xnfgs, ynfgs=ynfgs, xerrnfgs=xerrnfgs, yerrnfgs=yerrnfgs
    djs_oplot, !x.crange, alog10(1.0/HaHb)*[1,1], line=2, thick=postthick

    legend, '(c)', /left, /top, box=0, charsize=charsize_3, charthick=postthick
    legend, textoidl(xstr), /right, /top, box=0, charsize=charsize_3, charthick=postthick

    label = ['Stellar Absorption Corrected']
    legend, textoidl(label), /left, /bottom, box=0, charsize=charsize_3, $
      charthick=postthick

    im_openclose, postscript=postscript, /close

; ------------------------------------------------------------
; L(IR) vs A(Ha) - Integrated
; ------------------------------------------------------------

    psname = 'LIR_vs_AHa'
    im_openclose, pspath+psname, postscript=postscript, /encapsulated

    indx = where((atlasdust.ir_lum gt -900.0) and (atlasnodust.ebv_hahb_err gt 0.0),nindx)

    kha = k_lambda(6563,/odonnel)    

    x = atlasdust[indx].ir_lum
    xerr = atlasdust[indx].ir_lum_err

    y = atlasnodust[indx].ebv_hahb*kha
    yerr = atlasnodust[indx].ebv_hahb_err*kha

    indxnfgs = where((nfgsdust.ir_lum gt -900.0) and (nfgsnodust.ehbha_err gt 0.0),nindxnfgs)

    xnfgs = nfgsdust[indxnfgs].ir_lum
    xerrnfgs = nfgsdust[indxnfgs].ir_lum_err

    ynfgs = nfgsnodust[indxnfgs].ebv_hahb*kha
    yerrnfgs = nfgsnodust[indxnfgs].ebv_hahb_err*kha
    
    xtitle = 'log [L(IR)/L'+sunsymbol()+']'
    ytitle = 'A(H\alpha)'

    xrange = LIRrange
    yrange = AHarange

    atlas1d_lineplot, x, y, xerr, yerr, plottype=3, postscript=postscript, $
      xtitle=xtitle, ytitle=ytitle, xrange=xrange, yrange=yrange, legendtype=0, $
      /left, /top, xnfgs=xnfgs, ynfgs=ynfgs, xerrnfgs=xerrnfgs, yerrnfgs=yerrnfgs

    im_openclose, postscript=postscript, /close    
    
; ------------------------------------------------------------
; 2-panel - L(B) vs L(Ha)/L(IR)
; ------------------------------------------------------------
    
    psname = 'LB_vs_LHa_LIR_2panel'
    im_openclose, pspath+psname, postscript=postscript, /encapsulated, xsize=7.0, ysize=8.2

    pagemaker, nx=1, ny=2, height=[3.0,3.0], width=5.5, xmargin=[1.1,0.4], $
      ymargin=[1.1,1.1], xspace=0, yspace=0, xpage=7.0, ypage=8.2, $
      position=pos, /normal

    indx = where((atlasdust.b_lum_obs gt -900) and (atlasnodust.ir_flux gt -900) and $
      (atlasnodust.ebv_hahb_err gt 0.0),nindx)

    x = atlasdust[indx].b_lum_obs
    xerr = atlasdust[indx].b_lum_obs_err
    xabs = atlasdust[indx].m_b_obs
    
    lir = atlasdust[indx].ir_flux ; [erg/s/cm2]
    lir_err = sqrt(atlasdust[indx].ir_flux_err^2 + (atlasdust[indx].ir_flux*0.15)^2)

    indxnfgs = where((nfgsdust.b_lum_obs gt -900) and (nfgsdust.ir_flux gt -900) and $
      (nfgsnodust.ebv_hahb_err gt 0.0),nindxnfgs)

    xnfgs = nfgsnodust[indxnfgs].b_lum_obs
    xerrnfgs = nfgsnodust[indxnfgs].b_lum_obs_err
    
    lirnfgs = nfgsdust[indxnfgs].ir_flux
    lirnfgs_err = sqrt(nfgsdust[indxnfgs].ir_flux_err^2 + (nfgsdust[indxnfgs].ir_flux*0.15)^2)
    
    xrange = LBrange
    yrange = haLIRrange

    xtitle = 'log L(B) [L(B)'+sunsymbol()+']'
    ytitle = 'log [L(H\alpha)/L(IR)]'

    lhalir = alog10(4.5D-44/7.9D-42)
    
; ##########################################################
; Panel 1: No reddening correction
; ##########################################################

    ha = atlasdust[indx].h_alpha[0]
    ha_err = atlasdust[indx].h_alpha[1]
    
    y = alog10(ha/lir)
    yerr = im_compute_error(ha,ha_err,lir,lir_err,/log)

    hanfgs = nfgsdust[indxnfgs].h_alpha[0]
    hanfgs_err = nfgsdust[indxnfgs].h_alpha[1]
    
    ynfgs = alog10(hanfgs/lirnfgs)
    yerrnfgs = im_compute_error(hanfgs,hanfgs_err,lirnfgs,lirnfgs_err,/log)

    stats = im_stats([y,ynfgs],/verbose)
    xstr = strtrim(string(stats.median_rej,format='(F12.2)'),2)+' ('+$
      strtrim(string(stats.mean_rej,format='(F12.2)'),2)+$
      '\pm'+strtrim(string(stats.sigma_rej,format='(F12.2)'),2)+')'
;   xstr = strtrim(string(stats.median,format='(F12.2)'),2)+'\pm'+strtrim(string(stats.sig68mean,format='(F12.2)'),2)

    ytitle = 'log [L(H\alpha)_{obs}/L(IR)]'

    atlas1d_lineplot, x, y, xerr, yerr, plottype=3, postscript=postscript, $
      xtitle='', ytitle=ytitle, xrange=xrange, yrange=yrange, legendtype=0, $
      charsize=charsize, xstyle=11, position=pos[*,0], xnfgs=xnfgs, ynfgs=ynfgs, $
      xerrnfgs=xerrnfgs, yerrnfgs=yerrnfgs, xtickname=replicate(' ',10), yminor=3
    djs_oplot, !x.crange, lhalir*[1,1], line=0, thick=postthick
    axis, /xaxis, xrange=interpol(xabs,x,!x.crange), xthick=postthick, xstyle=1, $
      charsize=charsize, charthick=postthick
    xyouts, mean(!x.crange), !y.crange[1]+0.14*(!y.crange[1]-!y.crange[0]), textoidl('M_{B} [mag]'), $
      align=0.5, charsize=charsize, charthick=postthick

    legend, '(a)', /left, /top, box=0, charsize=charsize, charthick=postthick

    legend, textoidl(xstr), /right, /top, box=0, charsize=charsize_3, charthick=postthick
;   label = ['A(H\alpha) = 0']
;   legend, textoidl(label), /left, /bottom, box=0, charsize=charsize_3, $
;     charthick=postthick

; ############################################################
; Panel 2: Good Reddening Correction
; ############################################################

    ha = atlasnodust[indx].h_alpha[0]
    ha_err = atlasnodust[indx].h_alpha[1]

    y = alog10(ha/lir)
    yerr = im_compute_error(ha,ha_err,lir,lir_err,/log)

    hanfgs = nfgsnodust[indxnfgs].h_alpha[0]
    hanfgs_err = nfgsnodust[indxnfgs].h_alpha[1]

    ynfgs = alog10(hanfgs/lirnfgs)
    yerrnfgs = im_compute_error(hanfgs,hanfgs_err,lirnfgs,lirnfgs_err,/log)

    stats = im_stats([y,ynfgs],/verbose,/no_head)
    xstr = strtrim(string(stats.median_rej,format='(F12.2)'),2)+' ('+$
      strtrim(string(stats.mean_rej,format='(F12.2)'),2)+$
      '\pm'+strtrim(string(stats.sigma_rej,format='(F12.2)'),2)+')'
;   xstr = strtrim(string(stats.median,format='(F12.2)'),2)+'\pm'+strtrim(string(stats.sig68mean,format='(F12.2)'),2)

    ytitle = 'log [L(H\alpha)_{cor}/L(IR)]'

    atlas1d_lineplot, x, y, xerr, yerr, plottype=3, postscript=postscript, $
      xtitle=xtitle, ytitle=ytitle, xrange=xrange, yrange=yrange, legendtype=0, $
      charsize=charsize, position=pos[*,1], /noerase, xnfgs=xnfgs, ynfgs=ynfgs, $
      xerrnfgs=xerrnfgs, yerrnfgs=yerrnfgs, yminor=3
    djs_oplot, !x.crange, lhalir*[1,1], line=2, thick=postthick

    legend, '(b)', /left, /top, box=0, charsize=charsize, charthick=postthick

    legend, textoidl(xstr), /right, /top, box=0, charsize=charsize_3, charthick=postthick
;   label = ['Individual A(H\alpha)']
;   legend, textoidl(label), /left, /bottom, box=0, charsize=charsize_3, $
;     charthick=postthick

    im_openclose, postscript=postscript, /close

; ------------------------------------------------------------
; 2-panel - L(IR) vs L(Ha)/L(IR)+L(UV)
; ------------------------------------------------------------
    
    psname = 'LIR_vs_LHa_LIR_LUV'
    im_openclose, pspath+psname, postscript=postscript, /encapsulated, $
      xsize=8.5, ysize=7.5

    pagemaker, nx=1, ny=2, height=[3.0,3.0], width=5.5, xmargin=[2.0,1.0], $
      ymargin=[0.4,1.1], xspace=0, yspace=0, xpage=8.5, ypage=7.5, $
      position=pos, /normal

    indx = where((atlasnodust.ir_flux gt -900) and (atlasnodust.fuv_flux gt -900) and $
      (atlasnodust.ebv_hahb_err gt 0.0),nindx)

    x = atlasdust[indx].ir_lum
    xerr = atlasdust[indx].ir_lum_err

    lir = atlasdust[indx].ir_flux         ; [erg/s/cm2]
    lir_err = sqrt(atlasdust[indx].ir_flux_err^2 + (atlasdust[indx].ir_flux*0.15)^2)

    luv = atlasdust[indx].fuv_flux         ; [erg/s/cm2]
    luv_err = atlasdust[indx].fuv_flux_err ; [erg/s/cm2]

    lir_luv = lir + luv
    lir_luv_err = sqrt(lir_err^2 + luv_err^2)

    indxnfgs = where((nfgsdust.ir_flux gt -900) and (nfgsdust.fuv_flux gt -900) and $
      (nfgsnodust.ebv_hahb_err gt 0.0),nindxnfgs)

    xnfgs = nfgsnodust[indxnfgs].ir_lum
    xerrnfgs = nfgsnodust[indxnfgs].ir_lum_err
    
    lirnfgs = nfgsdust[indxnfgs].ir_flux ; [erg/s/cm2]
    lirnfgs_err = sqrt(nfgsdust[indxnfgs].ir_flux_err^2 + (nfgsdust[indxnfgs].ir_flux*0.15)^2)

    luvnfgs = nfgsdust[indxnfgs].fuv_flux ; [erg/s/cm2]
    luvnfgs_err = nfgsdust[indxnfgs].fuv_flux_err ; [erg/s/cm2]

    lir_luvnfgs = lirnfgs + luvnfgs
    lir_luvnfgs_err = sqrt(lirnfgs_err^2 + luvnfgs_err^2)

    xrange = LIRrange
    yrange = haLIRrange

    xtitle = 'log [L(IR)/L'+sunsymbol()+']'
    ytitle = 'log {L(H\alpha)/[L(IR)+L(UV)]}'

    lhalir = alog10(4.5D-44/7.9D-42)
    
; ##########################################################
; Panel 1: No reddening correction
; ##########################################################

    ha = atlasdust[indx].h_alpha[0]
    ha_err = atlasdust[indx].h_alpha[1]
    
    y = alog10(ha/lir_luv)
    yerr = im_compute_error(ha,ha_err,lir_luv,lir_luv_err,/log)

    hanfgs = nfgsdust[indxnfgs].h_alpha[0]
    hanfgs_err = nfgsdust[indxnfgs].h_alpha[1]
    
    ynfgs = alog10(hanfgs/lir_luvnfgs)
    yerrnfgs = im_compute_error(hanfgs,hanfgs_err,lir_luvnfgs,lir_luvnfgs_err,/log)

    stats = im_stats([y,ynfgs],/verbose)
    xstr = strtrim(string(stats.median_rej,format='(F12.2)'),2)+' ('+$
      strtrim(string(stats.mean_rej,format='(F12.2)'),2)+$
      '\pm'+strtrim(string(stats.sigma_rej,format='(F12.2)'),2)+')'
;   xstr = strtrim(string(stats.median,format='(F12.2)'),2)+'\pm'+strtrim(string(stats.sig68mean,format='(F12.2)'),2)

    ytitle = 'log {L(H\alpha)_{obs}/[L(IR)+L(UV)]}'

    atlas1d_lineplot, x, y, xerr, yerr, plottype=3, postscript=postscript, $
      xtitle='', ytitle=ytitle, xrange=xrange, yrange=yrange, legendtype=0, $
      charsize=charsize_3, position=pos[*,0], xnfgs=xnfgs, ynfgs=ynfgs, $
      xerrnfgs=xerrnfgs, yerrnfgs=yerrnfgs, xtickname=replicate(' ',10), yminor=3, $
      yfrac=8
    djs_oplot, !x.crange, lhalir*[1,1], line=2, thick=postthick

    legend, '(a)', /left, /top, box=0, charsize=charsize_3, charthick=postthick
    legend, textoidl(xstr), /right, /top, box=0, charsize=charsize_3, charthick=postthick

; ############################################################
; Panel 2: Good Reddening Correction
; ############################################################

    ha = atlasnodust[indx].h_alpha[0]
    ha_err = atlasnodust[indx].h_alpha[1]

    y = alog10(ha/lir_luv)
    yerr = im_compute_error(ha,ha_err,lir_luv,lir_luv_err,/log)

    hanfgs = nfgsnodust[indxnfgs].h_alpha[0]
    hanfgs_err = nfgsnodust[indxnfgs].h_alpha[1]

    ynfgs = alog10(hanfgs/lir_luvnfgs)
    yerrnfgs = im_compute_error(hanfgs,hanfgs_err,lir_luvnfgs,lir_luvnfgs_err,/log)

    stats = im_stats([y,ynfgs],/verbose,/no_head)
    xstr = strtrim(string(stats.median_rej,format='(F12.2)'),2)+' ('+$
      strtrim(string(stats.mean_rej,format='(F12.2)'),2)+$
      '\pm'+strtrim(string(stats.sigma_rej,format='(F12.2)'),2)+')'
;   xstr = strtrim(string(stats.median,format='(F12.2)'),2)+'\pm'+strtrim(string(stats.sig68mean,format='(F12.2)'),2)

    ytitle = 'log {L(H\alpha)_{cor}/[L(IR)+L(UV)]}'

    atlas1d_lineplot, x, y, xerr, yerr, plottype=3, postscript=postscript, $
      xtitle=xtitle, ytitle=ytitle, xrange=xrange, yrange=yrange, legendtype=0, $
      charsize=charsize_3, position=pos[*,1], /noerase, xnfgs=xnfgs, ynfgs=ynfgs, $
      xerrnfgs=xerrnfgs, yerrnfgs=yerrnfgs, yminor=3
    djs_oplot, !x.crange, lhalir*[1,1], line=2, thick=postthick

    legend, '(b)', /left, /top, box=0, charsize=charsize_3, charthick=postthick
    legend, textoidl(xstr), /right, /top, box=0, charsize=charsize_3, charthick=postthick
    
    im_openclose, postscript=postscript, /close

; ------------------------------------------------------------
; 2-panel - L(B) vs Hb/Ha - Integrated - Hb/Hg corrected
; ------------------------------------------------------------
    
    psname = 'LB_vs_hbha_hbhg_corrected'
    im_openclose, pspath+psname, postscript=postscript, /encapsulated, xsize=8.5, ysize=7.5

    pagemaker, nx=1, ny=2, height=[3.0,3.0], width=5.5, xmargin=[2.0,1.0], $
      ymargin=[0.4,1.1], xspace=0, yspace=0, xpage=8.5, ypage=7.5, $
      position=pos, /normal

    indx = where((atlasdust.b_lum_obs gt -900) and (atlasnodust.ebv_hahb_err gt 0.0) and $
      (atlasnodust.ebv_hbhg_err gt 0.0) and (atlasdust.h_beta_ew_uncor[0] ge ewcut),nindx)

    x = atlasdust[indx].b_lum_obs
    xabs = atlasdust[indx].m_b_obs
    xerr = atlasdust[indx].b_lum_obs_err
    
    indxnfgs = where((nfgsdust.b_lum_obs gt -900) and (nfgsnodust.ebv_hahb_err gt 0.0) and $
      (nfgsnodust.ebv_hbhg_err gt 0.0) and (nfgsdust.h_beta_ew_uncor[0] ge ewcut),nindxnfgs)

    xnfgs = nfgsnodust[indxnfgs].b_lum_obs
    xerrnfgs = nfgsnodust[indxnfgs].b_lum_obs_err
    
    xrange = LBrange
    yrange = [-2.5,2.1]

    xtitle = 'log L(B) [L(B)'+sunsymbol()+']'
    ytitle = 'log (H\beta/H\alpha)'

; ############################################################
; Panel 1: Absorption-Corrected, Individual A(Ha), Hb corrected using Hb/Hg
; ############################################################

    kl = k_lambda(4861,/odonnell)
    Ahb = atlasnodust[indx].ebv_hbhg*kl
    Ahb_err = atlasnodust[indx].ebv_hbhg_err*kl    
    hb = atlasdust[indx].h_beta[0]*10^(0.4*Ahb)
    hb_err = atlasdust[indx].h_beta[1]*10^(0.4*Ahb)+atlasdust[indx].h_beta[0]*alog(10.0)*10^(0.4*Ahb)*Ahb_err

    ha = atlasnodust[indx].h_alpha[0]
    ha_err = atlasnodust[indx].h_alpha[1]

    y = alog10(hb/ha)
    yerr = im_compute_error(hb,hb_err,ha,ha_err,/log)

    Ahbnfgs = nfgsnodust[indxnfgs].ebv_hbhg*kl
    Ahbnfgs_err = nfgsnodust[indxnfgs].ebv_hbhg_err*kl

    hbnfgs = nfgsdust[indxnfgs].h_beta[0]*10^(0.4*Ahbnfgs)
    hbnfgs_err = nfgsdust[indxnfgs].h_beta[1]*10^(0.4*Ahbnfgs) + $
      nfgsdust[indxnfgs].h_beta[0]*alog(10.0)*10^(0.4*Ahbnfgs)*Ahbnfgs_err

    hanfgs = nfgsnodust[indxnfgs].h_alpha[0]
    hanfgs_err = nfgsnodust[indxnfgs].h_alpha[1]

    ynfgs = alog10(hbnfgs/hanfgs)
    yerrnfgs = im_compute_error(hbnfgs,hbnfgs_err,hanfgs,hanfgs_err,/log)
       
    stats = im_stats([y,ynfgs],/verbose)
    xstr = strtrim(string(stats.median_rej,format='(F12.2)'),2)+' ('+$
      strtrim(string(stats.mean_rej,format='(F12.2)'),2)+$
      '\pm'+strtrim(string(stats.sigma_rej,format='(F12.2)'),2)+')'
;   xstr = strtrim(string(stats.median,format='(F12.2)'),2)+'\pm'+strtrim(string(stats.sig68mean,format='(F12.2)'),2)
        
    ytitle = 'log (H\beta/H\alpha)_{cor}'

    atlas1d_lineplot, x, y, xerr, yerr, plottype=3, postscript=postscript, $
      xtitle='', ytitle=ytitle, xrange=xrange, yrange=yrange, legendtype=0, yfrac=1.5, $
      charsize=charsize_3, position=pos[*,0], xtickname=replicate(' ',10), $
      xnfgs=xnfgs, ynfgs=ynfgs, xerrnfgs=xerrnfgs, yerrnfgs=yerrnfgs
    djs_oplot, !x.crange, alog10(1.0/HaHb)*[1,1], line=2, thick=postthick

    legend, '(a)', /left, /top, box=0, charsize=charsize_3, charthick=postthick
    legend, textoidl(xstr), /right, /top, box=0, charsize=charsize_3, charthick=postthick

    label = ['H\beta corrected using H\beta/H\gamma']
    legend, textoidl(label), /left, /bottom, box=0, charsize=charsize_3, $
      charthick=postthick

; ############################################################
; Panel 2: Absorption-Corrected, Individual A(Ha), Hb corrected using
; Hb/Hg, S/N cut on H-gamma
; ############################################################

    hgsnrcut = where(atlasdust[indx].h_gamma[0]/atlasdust[indx].h_gamma[1] gt 10.0)
    hgsnrcutnfgs = where(nfgsdust[indxnfgs].h_gamma[0]/nfgsdust[indxnfgs].h_gamma[1] gt 10.0)
    
    kl = k_lambda(4861,/odonnell)
    Ahb = atlasnodust[indx[hgsnrcut]].ebv_hbhg*kl
    Ahb_err = atlasnodust[indx[hgsnrcut]].ebv_hbhg_err*kl
    
    hb = atlasdust[indx[hgsnrcut]].h_beta[0]*10^(0.4*Ahb)
    hb_err = atlasdust[indx[hgsnrcut]].h_beta[1]*10^(0.4*Ahb) + $
      atlasdust[indx[hgsnrcut]].h_beta[0]*alog(10.0)*10^(0.4*Ahb)*Ahb_err

    ha = atlasnodust[indx[hgsnrcut]].h_alpha[0]
    ha_err = atlasnodust[indx[hgsnrcut]].h_alpha[1]

    y = alog10(hb/ha)
    yerr = im_compute_error(hb,hb_err,ha,ha_err,/log)

    Ahbnfgs = nfgsnodust[indxnfgs[hgsnrcutnfgs]].ebv_hbhg*kl
    Ahbnfgs_err = nfgsnodust[indxnfgs[hgsnrcutnfgs]].ebv_hbhg_err*kl

    hbnfgs = nfgsdust[indxnfgs[hgsnrcutnfgs]].h_beta[0]*10^(0.4*Ahbnfgs)
    hbnfgs_err = nfgsdust[indxnfgs[hgsnrcutnfgs]].h_beta[1]*10^(0.4*Ahbnfgs) + $
      nfgsdust[indxnfgs[hgsnrcutnfgs]].h_beta[0]*alog(10.0)*10^(0.4*Ahbnfgs)*Ahbnfgs_err

    hanfgs = nfgsnodust[indxnfgs[hgsnrcutnfgs]].h_alpha[0]
    hanfgs_err = nfgsnodust[indxnfgs[hgsnrcutnfgs]].h_alpha[1]

    ynfgs = alog10(hbnfgs/hanfgs)
    yerrnfgs = im_compute_error(hbnfgs,hbnfgs_err,hanfgs,hanfgs_err,/log)
       
    stats = im_stats([y,ynfgs],/verbose)
    xstr = strtrim(string(stats.median_rej,format='(F12.2)'),2)+' ('+$
      strtrim(string(stats.mean_rej,format='(F12.2)'),2)+$
      '\pm'+strtrim(string(stats.sigma_rej,format='(F12.2)'),2)+')'
;   xstr = strtrim(string(stats.median,format='(F12.2)'),2)+'\pm'+strtrim(string(stats.sig68mean,format='(F12.2)'),2)
        
    ytitle = 'log (H\beta/H\alpha)_{cor}'

    atlas1d_lineplot, x, y, xerr, yerr, plottype=3, postscript=postscript, $
      xtitle=xtitle, ytitle=ytitle, xrange=xrange, yrange=yrange, legendtype=0, yfrac=1.5, $
      charsize=charsize_3, position=pos[*,1], /noerase, $
      xnfgs=xnfgs, ynfgs=ynfgs, xerrnfgs=xerrnfgs, yerrnfgs=yerrnfgs
    djs_oplot, !x.crange, alog10(1.0/HaHb)*[1,1], line=2, thick=postthick

    legend, '(b)', /left, /top, box=0, charsize=charsize_3, charthick=postthick
    legend, textoidl(xstr), /right, /top, box=0, charsize=charsize_3, charthick=postthick

    label = ['H\beta corrected using H\beta/H\gamma','S/N [H\gamma] > 10']
    legend, textoidl(label), /left, /bottom, box=0, charsize=charsize_3, $
      charthick=postthick

    im_openclose, postscript=postscript, /close

; ------------------------------------------------------------
; compare the distributions of 60/100-micron ratios
; ------------------------------------------------------------
    
    psname = 'histogram_60_100'
    im_openclose, pspath+psname, postscript=postscript, /encapsulated

    indx = where(atlasdust.iras_60 gt -900.0 and atlasdust.iras_100 gt -900)
    indxnfgs = where(nfgsdust.iras_60 gt -900.0 and nfgsdust.iras_100 gt -900)

    x = alog10(atlasdust[indx].iras_60/atlasdust[indx].iras_100)
    xnfgs = alog10(nfgsdust[indxnfgs].iras_60/nfgsdust[indxnfgs].iras_100)

    stats = im_stats(x)
    nfgsstats = im_stats(xnfgs)

    xstr = ['Atlas: '+strtrim(string(stats.median,format='(F12.2)'),2)+'\pm' + $
      strtrim(string(stats.sig68mean,format='(F12.2)'),2),$
      'NFGS: '+strtrim(string(nfgsstats.median,format='(F12.2)'),2)+'\pm' + $
      strtrim(string(nfgsstats.sig68mean,format='(F12.2)'),2)]

    binsize = 0.05
    
    xbig = [x,xnfgs]
    plothist, xbig, bin=binsize, xbigbin, ybigbin, /noplot, /halfbin
    yrange = minmax(ybigbin)*[1.0,1.1]

    xrange = [-1,0.3]
    xtitle = 'log [S_{\nu}(60)/S_{\nu}(100)]'
    
    djs_plot, [0], [0], /nodata, xthick=postthick, ythick=postthick, charsize=2.0, $
      charthick=postthick, thick=postthick, xtitle=textoidl(xtitle), ytitle='Number', $
      ystyle=1, xrange=xrange, yrange=yrange, xstyle=1
    plothist, x, bin=binsize, /overplot, color=djs_icolor('blue'), thick=postthick, /halfbin
    plothist, xnfgs, bin=binsize, thick=postthick, line=0, /halfbin, $
      /overplot, /fill, /fline, forientation=45, fcolor=djs_icolor('red'), $
      color=djs_icolor('red'), fspacing=0.05
    plothist, xbig, bin=binsize, thick=postthick, line=2, /overplot, /halfbin
    legend, textoidl(xstr), /right, /top, box=0, charsize=charsize_3, charthick=postthick

    im_openclose, postscript=postscript, /close
    
; ------------------------------------------------------------
; SFR(Ha) versus SFR(Hb) - uncorrected
; ------------------------------------------------------------
    
    psname = 'sfr_ha_vs_sfr_hb_uncorrected'
    im_openclose, pspath+psname, postscript=postscript, /encapsulated, xsize=8.5, ysize=8.5

    pagemaker, nx=1, ny=2, xspace=0, yspace=0, xmargin=[2.0,0.5], $
      ymargin=[0.4,1.1], width=5.5, height=[4.5,2.5], xpage=8.5, ypage=8.5, $
      /normal, position=pos

    indx = where((atlas_sfrs.sfr_ha gt -900) and (atlas_sfrs.sfr_hb_uncor gt -900),nindx)
    
    x = atlas_sfrs[indx].sfr_ha
    xerr = atlas_sfrs[indx].sfr_ha_err

    y = atlas_sfrs[indx].sfr_hb_uncor
    yerr = atlas_sfrs[indx].sfr_hb_uncor_err

    resid = y-x
    resid_err = sqrt(yerr^2+xerr^2)
    
; NFGS    
    
    indxnfgs = where((nfgs_sfrs.sfr_ha gt -900) and (nfgs_sfrs.sfr_hb_uncor gt -900),nindxnfgs)

    xnfgs = nfgs_sfrs[indxnfgs].sfr_ha
    xerrnfgs = nfgs_sfrs[indxnfgs].sfr_ha_err

    ynfgs = nfgs_sfrs[indxnfgs].sfr_hb_uncor
    yerrnfgs = nfgs_sfrs[indxnfgs].sfr_hb_uncor_err

    residnfgs = ynfgs-xnfgs
    residnfgs_err = sqrt(yerrnfgs^2+xerrnfgs^2)
    
; make the plot    

    stats = im_stats([resid,residnfgs],/verbose)
    xstr = strtrim(string(stats.median_rej,format='(F12.2)'),2)+' ('+$
      strtrim(string(stats.mean_rej,format='(F12.2)'),2)+$
      '\pm'+strtrim(string(stats.sigma_rej,format='(F12.2)'),2)+')'
;   xstr = strtrim(string(stats.median,format='(F12.2)'),2)+'\pm'+strtrim(string(stats.sig68mean,format='(F12.2)'),2)
    
    xrange = sfrharange
    yrange = sfrharange
    yrange2 = residrange_sfrs

    xtitle = 'log \psi(H\alpha)'
    ytitle = 'log \psi(H\beta)'
    ytitle2 = 'log [\psi(H\beta)/\psi(H\alpha)]'

    atlas1d_lineplot, x, y, xerr, yerr, plottype=3, postscript=postscript, $
      xtitle='', ytitle=ytitle, xrange=xrange, yrange=yrange, legendtype=0, $
      charsize=charsize_5, position=pos[*,0], xtickname=replicate(' ',10), $
      xnfgs=xnfgs, ynfgs=ynfgs, xerrnfgs=xerrnfgs, yerrnfgs=yerrnfgs
    djs_oplot, !x.crange, !y.crange, line=0, thick=postthick

    label = [$
      'H\beta Uncorrected for Reddening',$
      'Observed EW(H\beta) > '+string(hbewcut_sfrs,format='(I0)')+' '+angstrom()$
      ]
    legend, textoidl(label), /left, /top, box=0, charsize=charsize_3, charthick=postthick, $
      spacing=1.8

    atlas1d_lineplot, x, resid, xerr, resid_err, plottype=3, postscript=postscript, $
      xtitle=xtitle, ytitle=ytitle2, xrange=xrange, yrange=yrange2, legendtype=0, $
      charsize=charsize_5, position=pos[*,1], /noerase, $
      xnfgs=xnfgs, ynfgs=residnfgs, xerrnfgs=xerrnfgs, yerrnfgs=residnfgs_err
    djs_oplot, !x.crange, [0,0], line=0, thick=postthick

    legend, textoidl(xstr), /right, /top, box=0, charsize=charsize_3, charthick=postthick

    im_openclose, postscript=postscript, /close

; ------------------------------------------------------------
; SFR(Ha) versus SFR(Hb) - Hb/Hg corrected
; ------------------------------------------------------------
    
    psname = 'sfr_ha_vs_sfr_hb_hbhg_corrected'
    im_openclose, pspath+psname, postscript=postscript, /encapsulated, xsize=8.5, ysize=8.5

    pagemaker, nx=1, ny=2, xspace=0, yspace=0, xmargin=[2.0,0.5], $
      ymargin=[0.4,1.1], width=5.5, height=[4.5,2.5], xpage=8.5, ypage=8.5, $
      /normal, position=pos

    indx = where((atlas_sfrs.sfr_ha gt -900) and (atlas_sfrs.sfr_hb_hbhg gt -900),nindx)
    
    x = atlas_sfrs[indx].sfr_ha
    xerr = atlas_sfrs[indx].sfr_ha_err

    y = atlas_sfrs[indx].sfr_hb_hbhg
    yerr = atlas_sfrs[indx].sfr_hb_hbhg_err

    resid = y-x
    resid_err = sqrt(yerr^2+xerr^2)
    
; NFGS    
    
    indxnfgs = where((nfgs_sfrs.sfr_ha gt -900) and (nfgs_sfrs.sfr_hb_hbhg gt -900),nindxnfgs)

    xnfgs = nfgs_sfrs[indxnfgs].sfr_ha
    xerrnfgs = nfgs_sfrs[indxnfgs].sfr_ha_err

    ynfgs = nfgs_sfrs[indxnfgs].sfr_hb_hbhg
    yerrnfgs = nfgs_sfrs[indxnfgs].sfr_hb_hbhg_err

    residnfgs = ynfgs-xnfgs
    residnfgs_err = sqrt(yerrnfgs^2+xerrnfgs^2)
    
; make the plot    

    stats = im_stats([resid,residnfgs],/verbose)
    xstr = strtrim(string(stats.median_rej,format='(F12.2)'),2)+' ('+$
      strtrim(string(stats.mean_rej,format='(F12.2)'),2)+$
      '\pm'+strtrim(string(stats.sigma_rej,format='(F12.2)'),2)+')'
;   xstr = strtrim(string(stats.median,format='(F12.2)'),2)+'\pm'+strtrim(string(stats.sig68mean,format='(F12.2)'),2)
    
    xrange = sfrharange
    yrange = sfrharange
    yrange2 = residrange_sfrs

    xtitle = 'log \psi(H\alpha)'
    ytitle = 'log \psi(H\beta)'
    ytitle2 = 'log [\psi(H\beta)/\psi(H\alpha)]'

    atlas1d_lineplot, x, y, xerr, yerr, plottype=3, postscript=postscript, $
      xtitle='', ytitle=ytitle, xrange=xrange, yrange=yrange, legendtype=0, $
      charsize=charsize_5, position=pos[*,0], xtickname=replicate(' ',10), $
      xnfgs=xnfgs, ynfgs=ynfgs, xerrnfgs=xerrnfgs, yerrnfgs=yerrnfgs
    djs_oplot, !x.crange, !y.crange, line=0, thick=postthick

    label = [$
      'H\beta corrected using H\beta/H\gamma ratio',$
      'Observed EW(H\beta) > '+string(hbewcut_sfrs,format='(I0)')+' '+angstrom(),$
      'Signal-to-Noise [H\gamma] > '+string(hgsnrcut_sfrs,format='(I0)')$
      ]
    legend, textoidl(label), /left, /top, box=0, charsize=charsize_3, charthick=postthick, $
      spacing=1.8

    atlas1d_lineplot, x, resid, xerr, resid_err, plottype=3, postscript=postscript, $
      xtitle=xtitle, ytitle=ytitle2, xrange=xrange, yrange=yrange2, legendtype=0, $
      charsize=charsize_5, position=pos[*,1], /noerase, $
      xnfgs=xnfgs, ynfgs=residnfgs, xerrnfgs=xerrnfgs, yerrnfgs=residnfgs_err
    djs_oplot, !x.crange, [0,0], line=0, thick=postthick

    legend, textoidl(xstr), /right, /top, box=0, charsize=charsize_3, charthick=postthick

    im_openclose, postscript=postscript, /close

; ------------------------------------------------------------
; SFR(Ha) versus SFR([O II]) - uncorrected
; ------------------------------------------------------------
    
    psname = 'sfr_ha_vs_sfr_oii_uncorrected'
    im_openclose, pspath+psname, postscript=postscript, /encapsulated, xsize=8.5, ysize=8.5

    pagemaker, nx=1, ny=2, xspace=0, yspace=0, xmargin=[2.0,0.5], $
      ymargin=[0.4,1.1], width=5.5, height=[4.5,2.5], xpage=8.5, ypage=8.5, $
      /normal, position=pos

    indx = where((atlas_sfrs.sfr_ha gt -900) and (atlas_sfrs.sfr_oii_uncor gt -900),nindx)
    
    x = atlas_sfrs[indx].sfr_ha
    xerr = atlas_sfrs[indx].sfr_ha_err

    y = atlas_sfrs[indx].sfr_oii_uncor
    yerr = atlas_sfrs[indx].sfr_oii_uncor_err

    resid = y-x
    resid_err = sqrt(yerr^2+xerr^2)
    
; NFGS    
    
    indxnfgs = where((nfgs_sfrs.sfr_ha gt -900) and (nfgs_sfrs.sfr_oii_uncor gt -900),nindxnfgs)

    xnfgs = nfgs_sfrs[indxnfgs].sfr_ha
    xerrnfgs = nfgs_sfrs[indxnfgs].sfr_ha_err

    ynfgs = nfgs_sfrs[indxnfgs].sfr_oii_uncor
    yerrnfgs = nfgs_sfrs[indxnfgs].sfr_oii_uncor_err

    residnfgs = ynfgs-xnfgs
    residnfgs_err = sqrt(yerrnfgs^2+xerrnfgs^2)
    
; make the plot    

    stats = im_stats([resid,residnfgs],/verbose)
    xstr = strtrim(string(stats.median_rej,format='(F12.2)'),2)+' ('+$
      strtrim(string(stats.mean_rej,format='(F12.2)'),2)+$
      '\pm'+strtrim(string(stats.sigma_rej,format='(F12.2)'),2)+')'
;   xstr = strtrim(string(stats.median,format='(F12.2)'),2)+'\pm'+strtrim(string(stats.sig68mean,format='(F12.2)'),2)
    
    xrange = sfrharange
    yrange = sfrharange
    yrange2 = residrange_sfrs

    xtitle = 'log \psi(H\alpha)'
    ytitle = 'log \psi([O II])'
    ytitle2 = 'log [\psi(H\beta)/\psi(H\alpha)]'

    atlas1d_lineplot, x, y, xerr, yerr, plottype=3, postscript=postscript, $
      xtitle='', ytitle=ytitle, xrange=xrange, yrange=yrange, legendtype=0, $
      charsize=charsize_5, position=pos[*,0], xtickname=replicate(' ',10), $
      xnfgs=xnfgs, ynfgs=ynfgs, xerrnfgs=xerrnfgs, yerrnfgs=yerrnfgs
    djs_oplot, !x.crange, !y.crange, line=0, thick=postthick

    label = [$
      '[O II] Uncorrected for Reddening'$
      ]
    legend, textoidl(label), /left, /top, box=0, charsize=charsize_3, charthick=postthick, $
      spacing=1.8

    atlas1d_lineplot, x, resid, xerr, resid_err, plottype=3, postscript=postscript, $
      xtitle=xtitle, ytitle=ytitle2, xrange=xrange, yrange=yrange2, legendtype=0, $
      charsize=charsize_5, position=pos[*,1], /noerase, $
      xnfgs=xnfgs, ynfgs=residnfgs, xerrnfgs=xerrnfgs, yerrnfgs=residnfgs_err
    djs_oplot, !x.crange, [0,0], line=0, thick=postthick

    legend, textoidl(xstr), /right, /top, box=0, charsize=charsize_3, charthick=postthick

    im_openclose, postscript=postscript, /close
    
; ------------------------------------------------------------
; SFR(Ha) versus SFR([O II]) - Hb/Hg corrected
; ------------------------------------------------------------
    
    psname = 'sfr_ha_vs_sfr_oii_hbhg_corrected'
    im_openclose, pspath+psname, postscript=postscript, /encapsulated, xsize=8.5, ysize=8.5

    pagemaker, nx=1, ny=2, xspace=0, yspace=0, xmargin=[2.0,0.5], $
      ymargin=[0.4,1.1], width=5.5, height=[4.5,2.5], xpage=8.5, ypage=8.5, $
      /normal, position=pos

    indx = where((atlas_sfrs.sfr_ha gt -900) and (atlas_sfrs.sfr_oii_hbhg gt -900),nindx)
    
    x = atlas_sfrs[indx].sfr_ha
    xerr = atlas_sfrs[indx].sfr_ha_err

    y = atlas_sfrs[indx].sfr_oii_hbhg
    yerr = atlas_sfrs[indx].sfr_oii_hbhg_err

    resid = y-x
    resid_err = sqrt(yerr^2+xerr^2)
    
; NFGS    
    
    indxnfgs = where((nfgs_sfrs.sfr_ha gt -900) and (nfgs_sfrs.sfr_oii_hbhg gt -900),nindxnfgs)

    xnfgs = nfgs_sfrs[indxnfgs].sfr_ha
    xerrnfgs = nfgs_sfrs[indxnfgs].sfr_ha_err

    ynfgs = nfgs_sfrs[indxnfgs].sfr_oii_hbhg
    yerrnfgs = nfgs_sfrs[indxnfgs].sfr_oii_hbhg_err

    residnfgs = ynfgs-xnfgs
    residnfgs_err = sqrt(yerrnfgs^2+xerrnfgs^2)
    
; make the plot    

    stats = im_stats([resid,residnfgs],/verbose)
    xstr = strtrim(string(stats.median_rej,format='(F12.2)'),2)+' ('+$
      strtrim(string(stats.mean_rej,format='(F12.2)'),2)+$
      '\pm'+strtrim(string(stats.sigma_rej,format='(F12.2)'),2)+')'
;   xstr = strtrim(string(stats.median,format='(F12.2)'),2)+'\pm'+strtrim(string(stats.sig68mean,format='(F12.2)'),2)
    
    xrange = sfrharange
    yrange = sfrharange
    yrange2 = residrange_sfrs

    xtitle = 'log \psi(H\alpha)'
    ytitle = 'log \psi([O II])'
    ytitle2 = 'log [\psi(H\beta)/\psi(H\alpha)]'

    atlas1d_lineplot, x, y, xerr, yerr, plottype=3, postscript=postscript, $
      xtitle='', ytitle=ytitle, xrange=xrange, yrange=yrange, legendtype=0, $
      charsize=charsize_5, position=pos[*,0], xtickname=replicate(' ',10), $
      xnfgs=xnfgs, ynfgs=ynfgs, xerrnfgs=xerrnfgs, yerrnfgs=yerrnfgs
    djs_oplot, !x.crange, !y.crange, line=0, thick=postthick

    label = [$
      '[O II] corrected using H\beta/H\gamma ratio',$
      'Signal-to-Noise [H\gamma] > '+string(hgsnrcut_sfrs,format='(I0)')$
      ]
    legend, textoidl(label), /left, /top, box=0, charsize=charsize_3, charthick=postthick, $
      spacing=1.8

    atlas1d_lineplot, x, resid, xerr, resid_err, plottype=3, postscript=postscript, $
      xtitle=xtitle, ytitle=ytitle2, xrange=xrange, yrange=yrange2, legendtype=0, $
      charsize=charsize_5, position=pos[*,1], /noerase, $
      xnfgs=xnfgs, ynfgs=residnfgs, xerrnfgs=xerrnfgs, yerrnfgs=residnfgs_err
    djs_oplot, !x.crange, [0,0], line=0, thick=postthick

    legend, textoidl(xstr), /right, /top, box=0, charsize=charsize_3, charthick=postthick

    im_openclose, postscript=postscript, /close

; ------------------------------------------------------------
; inter-compare all the various SFR([O II]) calibrations and SFR(Ha)
; ------------------------------------------------------------

    psname = 'sfr_oii_intercomparison'
    im_openclose, pspath+psname, postscript=postscript, /encapsulated

    pagemaker, nx=3, ny=3, position=pos, /normal, xspace=0.0, $
      xmargin=[1.3,0.2], ymargin=[0.2,1.1], yspace=0.0

    xrange = [6.9,9.5]
    yrange = xrange
    charsize = 1.3

; --------------------------------------------------    
; Panel 1 - SFR(Ha) vs SFR([O II]) [K04]
; --------------------------------------------------    

    indx = where((atlas_sfrs.sfr_ha gt -900) and (atlas_sfrs.sfr_oii_k04 gt -900),nindx)
    
    x = atlas_sfrs[indx].sfr_ha
    xerr = atlas_sfrs[indx].sfr_ha_err

    y = atlas_sfrs[indx].sfr_oii_k04
    yerr = atlas_sfrs[indx].sfr_oii_k04_err

    resid = y-x
    resid_err = sqrt(yerr^2+xerr^2)
    
; NFGS    
    
    indxnfgs = where((nfgs_sfrs.sfr_ha gt -900) and (nfgs_sfrs.sfr_oii_k04 gt -900),nindxnfgs)

    xnfgs = nfgs_sfrs[indxnfgs].sfr_ha
    xerrnfgs = nfgs_sfrs[indxnfgs].sfr_ha_err

    ynfgs = nfgs_sfrs[indxnfgs].sfr_oii_k04
    yerrnfgs = nfgs_sfrs[indxnfgs].sfr_oii_k04_err

    residnfgs = ynfgs-xnfgs
    residnfgs_err = sqrt(yerrnfgs^2+xerrnfgs^2)
    
; make the plot    

    stats = im_stats([resid,residnfgs],/verbose)
    xstr = strtrim(string(stats.median_rej,format='(F12.2)'),2)+' ('+$
      strtrim(string(stats.mean_rej,format='(F12.2)'),2)+$
      '\pm'+strtrim(string(stats.sigma_rej,format='(F12.2)'),2)+')'
;   xstr = strtrim(string(stats.median,format='(F12.2)'),2)+'\pm'+strtrim(string(stats.sig68mean,format='(F12.2)'),2)
    
    xrange = sfrharange
    yrange = sfrharange
    yrange2 = residrange_sfrs

    xtitle = 'log \psi(H\alpha)'
    ytitle = 'log \psi([O II]) [K04]'
    ytitle2 = 'log [\psi(H\beta)/\psi(H\alpha)]'

    atlas1d_lineplot, x, y, xerr, yerr, plottype=1, postscript=postscript, $
      xtitle='', ytitle=ytitle, xrange=xrange, yrange=yrange, legendtype=0, $
      charsize=charsize_5, position=pos[*,0], xtickname=replicate(' ',10), $
      xnfgs=xnfgs, ynfgs=ynfgs, xerrnfgs=xerrnfgs, yerrnfgs=yerrnfgs, xminor=3, $
      yminor=3
    djs_oplot, !x.crange, !y.crange, line=0, thick=postthick

    legend, '(a)', /left, /top, box=0, charsize=charsize_5, charthick=postthick
    legend, textoidl(xstr), /right, /bottom, box=0, charsize=charsize_3, charthick=postthick

; --------------------------------------------------    
; Panel 2 - SFR([O II]) [Hb/Hg] vs SFR([O II]) [K04]
; --------------------------------------------------    

    indx = where((atlas_sfrs.sfr_oii_hbhg gt -900) and (atlas_sfrs.sfr_oii_k04 gt -900),nindx)
    
    x = atlas_sfrs[indx].sfr_oii_hbhg
    xerr = atlas_sfrs[indx].sfr_oii_hbhg_err

    y = atlas_sfrs[indx].sfr_oii_k04
    yerr = atlas_sfrs[indx].sfr_oii_k04_err

    resid = y-x
    resid_err = sqrt(yerr^2+xerr^2)
    
; NFGS    
    
    indxnfgs = where((nfgs_sfrs.sfr_oii_hbhg gt -900) and (nfgs_sfrs.sfr_oii_k04 gt -900),nindxnfgs)

    xnfgs = nfgs_sfrs[indxnfgs].sfr_oii_hbhg
    xerrnfgs = nfgs_sfrs[indxnfgs].sfr_oii_hbhg_err

    ynfgs = nfgs_sfrs[indxnfgs].sfr_oii_k04
    yerrnfgs = nfgs_sfrs[indxnfgs].sfr_oii_k04_err

    residnfgs = ynfgs-xnfgs
    residnfgs_err = sqrt(yerrnfgs^2+xerrnfgs^2)
    
; make the plot    

    stats = im_stats([resid,residnfgs],/verbose)
    xstr = strtrim(string(stats.median_rej,format='(F12.2)'),2)+' ('+$
      strtrim(string(stats.mean_rej,format='(F12.2)'),2)+$
      '\pm'+strtrim(string(stats.sigma_rej,format='(F12.2)'),2)+')'
;   xstr = strtrim(string(stats.median,format='(F12.2)'),2)+'\pm'+strtrim(string(stats.sig68mean,format='(F12.2)'),2)
    
    xrange = sfrharange
    yrange = sfrharange
    yrange2 = residrange_sfrs

    xtitle = 'log \psi([O II]) [H\beta/H\gamma]'
    ytitle = 'log \psi([O II]) [K04]'
    ytitle2 = 'log [\psi(H\beta)/\psi(H\alpha)]'

    atlas1d_lineplot, x, y, xerr, yerr, plottype=1, postscript=postscript, $
      xtitle='', ytitle='', xrange=xrange, yrange=yrange, legendtype=0, $
      charsize=charsize_5, position=pos[*,1], xtickname=replicate(' ',10), $
      /noerase, ytickname=replicate(' ',10), $
      xnfgs=xnfgs, ynfgs=ynfgs, xerrnfgs=xerrnfgs, yerrnfgs=yerrnfgs, xminor=3, $
      yminor=3
    djs_oplot, !x.crange, !y.crange, line=0, thick=postthick

    legend, '(b)', /left, /top, box=0, charsize=charsize_5, charthick=postthick
    legend, textoidl(xstr), /right, /bottom, box=0, charsize=charsize_3, charthick=postthick

; --------------------------------------------------    
; Panel 3 - SFR([O II]) [K98] vs SFR([O II]) [K04]
; --------------------------------------------------    

    indx = where((atlas_sfrs.sfr_oii_k98 gt -900) and (atlas_sfrs.sfr_oii_k04 gt -900),nindx)
    
    x = atlas_sfrs[indx].sfr_oii_k98
    xerr = atlas_sfrs[indx].sfr_oii_k98_err

    y = atlas_sfrs[indx].sfr_oii_k04
    yerr = atlas_sfrs[indx].sfr_oii_k04_err

    resid = y-x
    resid_err = sqrt(yerr^2+xerr^2)
    
; NFGS    
    
    indxnfgs = where((nfgs_sfrs.sfr_oii_k98 gt -900) and (nfgs_sfrs.sfr_oii_k04 gt -900),nindxnfgs)

    xnfgs = nfgs_sfrs[indxnfgs].sfr_oii_k98
    xerrnfgs = nfgs_sfrs[indxnfgs].sfr_oii_k98_err

    ynfgs = nfgs_sfrs[indxnfgs].sfr_oii_k04
    yerrnfgs = nfgs_sfrs[indxnfgs].sfr_oii_k04_err

    residnfgs = ynfgs-xnfgs
    residnfgs_err = sqrt(yerrnfgs^2+xerrnfgs^2)
    
; make the plot    

    stats = im_stats([resid,residnfgs],/verbose)
    xstr = strtrim(string(stats.median_rej,format='(F12.2)'),2)+' ('+$
      strtrim(string(stats.mean_rej,format='(F12.2)'),2)+$
      '\pm'+strtrim(string(stats.sigma_rej,format='(F12.2)'),2)+')'
;   xstr = strtrim(string(stats.median,format='(F12.2)'),2)+'\pm'+strtrim(string(stats.sig68mean,format='(F12.2)'),2)
    
    xrange = sfrharange
    yrange = sfrharange
    yrange2 = residrange_sfrs

    xtitle = 'log \psi([O II]) [K98]'
    ytitle = 'log \psi([O II]) [K04]'
    ytitle2 = 'log [\psi(H\beta)/\psi(H\alpha)]'

    atlas1d_lineplot, x, y, xerr, yerr, plottype=1, postscript=postscript, $
      xtitle=xtitle, ytitle='', xrange=xrange, yrange=yrange, legendtype=0, $
      charsize=charsize_5, position=pos[*,2], /noerase, ytickname=replicate(' ',10), $
      xnfgs=xnfgs, ynfgs=ynfgs, xerrnfgs=xerrnfgs, yerrnfgs=yerrnfgs, xminor=3, $
      yminor=3
    djs_oplot, !x.crange, !y.crange, line=0, thick=postthick

    legend, '(c)', /left, /top, box=0, charsize=charsize_5, charthick=postthick
    legend, textoidl(xstr), /right, /bottom, box=0, charsize=charsize_3, charthick=postthick

; --------------------------------------------------    
; Panel 4 - SFR(Ha) vs SFR([O II]) [K98]
; --------------------------------------------------    

    indx = where((atlas_sfrs.sfr_ha gt -900) and (atlas_sfrs.sfr_oii_k98 gt -900),nindx)
    
    x = atlas_sfrs[indx].sfr_ha
    xerr = atlas_sfrs[indx].sfr_ha_err

    y = atlas_sfrs[indx].sfr_oii_k98
    yerr = atlas_sfrs[indx].sfr_oii_k98_err

    resid = y-x
    resid_err = sqrt(yerr^2+xerr^2)
    
; NFGS    
    
    indxnfgs = where((nfgs_sfrs.sfr_ha gt -900) and (nfgs_sfrs.sfr_oii_k98 gt -900),nindxnfgs)

    xnfgs = nfgs_sfrs[indxnfgs].sfr_ha
    xerrnfgs = nfgs_sfrs[indxnfgs].sfr_ha_err

    ynfgs = nfgs_sfrs[indxnfgs].sfr_oii_k98
    yerrnfgs = nfgs_sfrs[indxnfgs].sfr_oii_k98_err

    residnfgs = ynfgs-xnfgs
    residnfgs_err = sqrt(yerrnfgs^2+xerrnfgs^2)
    
; make the plot    

    stats = im_stats([resid,residnfgs],/verbose)
    xstr = strtrim(string(stats.median_rej,format='(F12.2)'),2)+' ('+$
      strtrim(string(stats.mean_rej,format='(F12.2)'),2)+$
      '\pm'+strtrim(string(stats.sigma_rej,format='(F12.2)'),2)+')'
;   xstr = strtrim(string(stats.median,format='(F12.2)'),2)+'\pm'+strtrim(string(stats.sig68mean,format='(F12.2)'),2)
    
    xrange = sfrharange
    yrange = sfrharange
    yrange2 = residrange_sfrs

    xtitle = 'log \psi(H\alpha)'
    ytitle = 'log \psi([O II]) [K98]'
    ytitle2 = 'log [\psi(H\beta)/\psi(H\alpha)]'

    atlas1d_lineplot, x, y, xerr, yerr, plottype=1, postscript=postscript, $
      xtitle='', ytitle=ytitle, xrange=xrange, yrange=yrange, legendtype=0, $
      charsize=charsize_5, position=pos[*,3], /noerase, xtickname=replicate(' ',10), $
      xnfgs=xnfgs, ynfgs=ynfgs, xerrnfgs=xerrnfgs, yerrnfgs=yerrnfgs, xminor=3, $
      yminor=3
    djs_oplot, !x.crange, !y.crange, line=0, thick=postthick

    legend, '(d)', /left, /top, box=0, charsize=charsize_5, charthick=postthick
    legend, textoidl(xstr), /right, /bottom, box=0, charsize=charsize_3, charthick=postthick

; --------------------------------------------------    
; Panel 5 - SFR([O II]) [Hb/Hg] vs SFR([O II]) [K98]
; --------------------------------------------------    

    indx = where((atlas_sfrs.sfr_oii_hbhg gt -900) and (atlas_sfrs.sfr_oii_k98 gt -900),nindx)
    
    x = atlas_sfrs[indx].sfr_oii_hbhg
    xerr = atlas_sfrs[indx].sfr_oii_hbhg_err

    y = atlas_sfrs[indx].sfr_oii_k98
    yerr = atlas_sfrs[indx].sfr_oii_k98_err

    resid = y-x
    resid_err = sqrt(yerr^2+xerr^2)
    
; NFGS    
    
    indxnfgs = where((nfgs_sfrs.sfr_oii_hbhg gt -900) and (nfgs_sfrs.sfr_oii_k98 gt -900),nindxnfgs)

    xnfgs = nfgs_sfrs[indxnfgs].sfr_oii_hbhg
    xerrnfgs = nfgs_sfrs[indxnfgs].sfr_oii_hbhg_err

    ynfgs = nfgs_sfrs[indxnfgs].sfr_oii_k98
    yerrnfgs = nfgs_sfrs[indxnfgs].sfr_oii_k98_err

    residnfgs = ynfgs-xnfgs
    residnfgs_err = sqrt(yerrnfgs^2+xerrnfgs^2)
    
; make the plot    

    stats = im_stats([resid,residnfgs],/verbose)
    xstr = strtrim(string(stats.median_rej,format='(F12.2)'),2)+' ('+$
      strtrim(string(stats.mean_rej,format='(F12.2)'),2)+$
      '\pm'+strtrim(string(stats.sigma_rej,format='(F12.2)'),2)+')'
;   xstr = strtrim(string(stats.median,format='(F12.2)'),2)+'\pm'+strtrim(string(stats.sig68mean,format='(F12.2)'),2)
    
    xrange = sfrharange
    yrange = sfrharange
    yrange2 = residrange_sfrs

    xtitle = 'log \psi([O II]) [H\beta/H\gamma]'
    ytitle = 'log \psi([O II]) [K98]'
    ytitle2 = 'log [\psi(H\beta)/\psi(H\alpha)]'

    atlas1d_lineplot, x, y, xerr, yerr, plottype=1, postscript=postscript, $
      xtitle=xtitle, ytitle='', xrange=xrange, yrange=yrange, legendtype=0, $
      charsize=charsize_5, position=pos[*,4], /noerase, ytickname=replicate(' ',10), $
      xnfgs=xnfgs, ynfgs=ynfgs, xerrnfgs=xerrnfgs, yerrnfgs=yerrnfgs, xminor=3, $
      yminor=3
    djs_oplot, !x.crange, !y.crange, line=0, thick=postthick

    legend, '(e)', /left, /top, box=0, charsize=charsize_5, charthick=postthick
    legend, textoidl(xstr), /right, /bottom, box=0, charsize=charsize_3, charthick=postthick

; --------------------------------------------------    
; Panel 6 - Blank
; --------------------------------------------------    

; --------------------------------------------------    
; Panel 7 - SFR(Ha) vs SFR([O II]) [Hb/Hg]
; --------------------------------------------------    

    indx = where((atlas_sfrs.sfr_ha gt -900) and (atlas_sfrs.sfr_oii_hbhg gt -900),nindx)
    
    x = atlas_sfrs[indx].sfr_ha
    xerr = atlas_sfrs[indx].sfr_ha_err

    y = atlas_sfrs[indx].sfr_oii_hbhg
    yerr = atlas_sfrs[indx].sfr_oii_hbhg_err

    resid = y-x
    resid_err = sqrt(yerr^2+xerr^2)
    
; NFGS    
    
    indxnfgs = where((nfgs_sfrs.sfr_ha gt -900) and (nfgs_sfrs.sfr_oii_hbhg gt -900),nindxnfgs)

    xnfgs = nfgs_sfrs[indxnfgs].sfr_ha
    xerrnfgs = nfgs_sfrs[indxnfgs].sfr_ha_err

    ynfgs = nfgs_sfrs[indxnfgs].sfr_oii_hbhg
    yerrnfgs = nfgs_sfrs[indxnfgs].sfr_oii_hbhg_err

    residnfgs = ynfgs-xnfgs
    residnfgs_err = sqrt(yerrnfgs^2+xerrnfgs^2)
    
; make the plot    

    stats = im_stats([resid,residnfgs],/verbose)
    xstr = strtrim(string(stats.median_rej,format='(F12.2)'),2)+' ('+$
      strtrim(string(stats.mean_rej,format='(F12.2)'),2)+$
      '\pm'+strtrim(string(stats.sigma_rej,format='(F12.2)'),2)+')'
;   xstr = strtrim(string(stats.median,format='(F12.2)'),2)+'\pm'+strtrim(string(stats.sig68mean,format='(F12.2)'),2)
    
    xrange = sfrharange
    yrange = sfrharange
    yrange2 = residrange_sfrs

    xtitle = 'log \psi(H\alpha)'
    ytitle = 'log \psi([O II]) [H\beta/H\gamma]'
    ytitle2 = 'log [\psi(H\beta)/\psi(H\alpha)]'

    atlas1d_lineplot, x, y, xerr, yerr, plottype=1, postscript=postscript, $
      xtitle=xtitle, ytitle=ytitle, xrange=xrange, yrange=yrange, legendtype=0, $
      charsize=charsize_5, position=pos[*,6], /noerase, $
      xnfgs=xnfgs, ynfgs=ynfgs, xerrnfgs=xerrnfgs, yerrnfgs=yerrnfgs, xminor=3, $
      yminor=3
    djs_oplot, !x.crange, !y.crange, line=0, thick=postthick

    legend, '(f)', /left, /top, box=0, charsize=charsize_5, charthick=postthick
    legend, textoidl(xstr), /right, /bottom, box=0, charsize=charsize_3, charthick=postthick

; --------------------------------------------------    
; Panel 8 - Blank
; --------------------------------------------------    

; --------------------------------------------------    
; Panel 9 - Blank
; --------------------------------------------------    

    im_openclose, postscript=postscript, /close    

; ------------------------------------------------------------
; SFR(Ha) vs various SFR(Hb) calibrations
; ------------------------------------------------------------
    
    psname = 'sfr_ha_vs_sfr_hb_multipanel'
    im_openclose, pspath+psname, postscript=postscript, /encapsulated, xsize=8.5, ysize=9.5

    pagemaker, nx=1, ny=3, height=[2.5,2.5,2.5], width=5.5, xmargin=[2.0,1.0], $
      ymargin=[0.9,1.1], xspace=0, yspace=0, xpage=8.5, ypage=9.5, $
      position=pos, /normal

; ##########################################################
; Panel 1
; ##########################################################

; Atlas    
    
    indx = where((atlas_sfrs.sfr_ha gt -900) and (atlas_sfrs.sfr_hb_uncor gt -900),nindx)

    x = atlas_sfrs[indx].sfr_ha
    xerr = atlas_sfrs[indx].sfr_ha_err
    
    y = atlas_sfrs[indx].sfr_hb_uncor
    yerr = atlas_sfrs[indx].sfr_hb_uncor_err

    resid = y-x
    resid_err = sqrt(yerr^2+xerr^2)

; NFGS    
    
    indxnfgs = where((nfgs_sfrs.sfr_ha gt -900) and (nfgs_sfrs.sfr_hb_uncor gt -900),nindxnfgs)

    xnfgs = nfgs_sfrs[indxnfgs].sfr_ha
    xerrnfgs = nfgs_sfrs[indxnfgs].sfr_ha_err
    
    ynfgs = nfgs_sfrs[indxnfgs].sfr_hb_uncor
    yerrnfgs = nfgs_sfrs[indxnfgs].sfr_hb_uncor_err

    residnfgs = ynfgs-xnfgs
    residnfgs_err = sqrt(yerrnfgs^2+xerrnfgs^2)

; make the plot    

    stats = im_stats([resid,residnfgs],/verbose)
    xstr = strtrim(string(stats.median_rej,format='(F12.2)'),2)+' ('+$
      strtrim(string(stats.mean_rej,format='(F12.2)'),2)+$
      '\pm'+strtrim(string(stats.sigma_rej,format='(F12.2)'),2)+')'
;   xstr = strtrim(string(stats.median,format='(F12.2)'),2)+'\pm'+strtrim(string(stats.sig68mean,format='(F12.2)'),2)
    
    xtitle = 'log \psi(H\alpha)'
    ytitle = 'log \psi(H\beta) [Observed]'

    xrange = sfrharange
    yrange = xrange
    
    atlas1d_lineplot, x, y, xerr, yerr, plottype=1, postscript=postscript, $
      xtitle='', ytitle=ytitle, xrange=xrange, yrange=yrange, legendtype=0, $
      charsize=charsize, xtickname=replicate(' ',10), position=pos[*,0], $
      xminor=3, yminor=3, $
      xnfgs=xnfgs, ynfgs=ynfgs, xerrnfgs=xerrnfgs, yerrnfgs=yerrnfgs
    djs_oplot, !x.crange, !y.crange, line=0, thick=postthick
    
    legend, '(a)', /left, /top, box=0, charsize=charsize_5, charthick=postthick
    legend, textoidl(xstr), /right, /bottom, box=0, charsize=charsize_3, charthick=postthick

; ##########################################################
; Panel 2
; ##########################################################

; Atlas    
    
    indx = where((atlas_sfrs.sfr_ha gt -900) and (atlas_sfrs.sfr_hb_hbhg gt -900),nindx)

    x = atlas_sfrs[indx].sfr_ha
    xerr = atlas_sfrs[indx].sfr_ha_err
    
    y = atlas_sfrs[indx].sfr_hb_hbhg
    yerr = atlas_sfrs[indx].sfr_hb_hbhg_err

    resid = y-x
    resid_err = sqrt(yerr^2+xerr^2)

; NFGS    
    
    indxnfgs = where((nfgs_sfrs.sfr_ha gt -900) and (nfgs_sfrs.sfr_hb_hbhg gt -900),nindxnfgs)

    xnfgs = nfgs_sfrs[indxnfgs].sfr_ha
    xerrnfgs = nfgs_sfrs[indxnfgs].sfr_ha_err
    
    ynfgs = nfgs_sfrs[indxnfgs].sfr_hb_hbhg
    yerrnfgs = nfgs_sfrs[indxnfgs].sfr_hb_hbhg_err

    residnfgs = ynfgs-xnfgs
    residnfgs_err = sqrt(yerrnfgs^2+xerrnfgs^2)

; make the plot    

    stats = im_stats([resid,residnfgs],/verbose)
    xstr = strtrim(string(stats.median_rej,format='(F12.2)'),2)+' ('+$
      strtrim(string(stats.mean_rej,format='(F12.2)'),2)+$
      '\pm'+strtrim(string(stats.sigma_rej,format='(F12.2)'),2)+')'
;   xstr = strtrim(string(stats.median,format='(F12.2)'),2)+'\pm'+strtrim(string(stats.sig68mean,format='(F12.2)'),2)
    
    xtitle = 'log \psi(H\alpha)'
    ytitle = 'log \psi(H\beta) [H\beta/H\gamma]'

    xrange = sfrharange
    yrange = xrange
    
    atlas1d_lineplot, x, y, xerr, yerr, plottype=1, postscript=postscript, $
      xtitle='', ytitle=ytitle, xrange=xrange, yrange=yrange, legendtype=0, $
      charsize=charsize, xtickname=replicate(' ',10), position=pos[*,1], $
      xminor=3, yminor=3, /noerase, $
      xnfgs=xnfgs, ynfgs=ynfgs, xerrnfgs=xerrnfgs, yerrnfgs=yerrnfgs
    djs_oplot, !x.crange, !y.crange, line=0, thick=postthick
    
    legend, '(b)', /left, /top, box=0, charsize=charsize_5, charthick=postthick
    legend, textoidl(xstr), /right, /bottom, box=0, charsize=charsize_3, charthick=postthick

; ##########################################################
; Panel 3
; ##########################################################

; Atlas    
    
    indx = where((atlas_sfrs.sfr_ha gt -900) and (atlas_sfrs.sfr_hb_best gt -900),nindx)

    x = atlas_sfrs[indx].sfr_ha
    xerr = atlas_sfrs[indx].sfr_ha_err
    
    y = atlas_sfrs[indx].sfr_hb_best
    yerr = atlas_sfrs[indx].sfr_hb_best_err

    resid = y-x
    resid_err = sqrt(yerr^2+xerr^2)

; NFGS    
    
    indxnfgs = where((nfgs_sfrs.sfr_ha gt -900) and (nfgs_sfrs.sfr_hb_best gt -900),nindxnfgs)

    xnfgs = nfgs_sfrs[indxnfgs].sfr_ha
    xerrnfgs = nfgs_sfrs[indxnfgs].sfr_ha_err
    
    ynfgs = nfgs_sfrs[indxnfgs].sfr_hb_best
    yerrnfgs = nfgs_sfrs[indxnfgs].sfr_hb_best_err

    residnfgs = ynfgs-xnfgs
    residnfgs_err = sqrt(yerrnfgs^2+xerrnfgs^2)

; make the plot    

    stats = im_stats([resid,residnfgs],/verbose)
    xstr = strtrim(string(stats.median_rej,format='(F12.2)'),2)+' ('+$
      strtrim(string(stats.mean_rej,format='(F12.2)'),2)+$
      '\pm'+strtrim(string(stats.sigma_rej,format='(F12.2)'),2)+')'
;   xstr = strtrim(string(stats.median,format='(F12.2)'),2)+'\pm'+strtrim(string(stats.sig68mean,format='(F12.2)'),2)
    
    xtitle = 'log \psi(H\alpha)'
    ytitle = 'log \psi(H\beta) [This Paper]'

    xrange = sfrharange
    yrange = xrange
    
    atlas1d_lineplot, x, y, xerr, yerr, plottype=1, postscript=postscript, $
      xtitle=xtitle, ytitle=ytitle, xrange=xrange, yrange=yrange, legendtype=0, $
      charsize=charsize, position=pos[*,2], $
      xminor=3, yminor=3, /noerase, $
      xnfgs=xnfgs, ynfgs=ynfgs, xerrnfgs=xerrnfgs, yerrnfgs=yerrnfgs
    djs_oplot, !x.crange, !y.crange, line=0, thick=postthick
    
    legend, '(c)', /left, /top, box=0, charsize=charsize_5, charthick=postthick
    legend, textoidl(xstr), /right, /bottom, box=0, charsize=charsize_3, charthick=postthick
    
    im_openclose, postscript=postscript, /close    

; ------------------------------------------------------------
; SFR(Ha) vs various SFR([O II]) calibrations
; ------------------------------------------------------------

    psname = 'sfr_ha_vs_sfr_oii_multipanel'
    im_openclose, pspath+psname, postscript=postscript, xsize=8.5, ysize=11.0, /encapsulated

    pagemaker, nx=2, ny=3, yspace=0, xspace=0, width=3.15*[1,1], height=3.15*[1,1,1], $
      xmargin=[1.1,1.1], ymargin=[0.45,1.1], xpage=8.5, ypage=11.0, position=pos, /normal

; ##########################################################
; Panel 1
; ##########################################################

; Atlas    
    
    indx = where((atlas_sfrs.sfr_ha gt -900) and (atlas_sfrs.sfr_oii_uncor gt -900),nindx)

    x = atlas_sfrs[indx].sfr_ha
    xerr = atlas_sfrs[indx].sfr_ha_err
    
    y = atlas_sfrs[indx].sfr_oii_uncor
    yerr = atlas_sfrs[indx].sfr_oii_uncor_err

    resid = y-x
    resid_err = sqrt(yerr^2+xerr^2)

; NFGS    
    
    indxnfgs = where((nfgs_sfrs.sfr_ha gt -900) and (nfgs_sfrs.sfr_oii_uncor gt -900),nindxnfgs)

    xnfgs = nfgs_sfrs[indxnfgs].sfr_ha
    xerrnfgs = nfgs_sfrs[indxnfgs].sfr_ha_err
    
    ynfgs = nfgs_sfrs[indxnfgs].sfr_oii_uncor
    yerrnfgs = nfgs_sfrs[indxnfgs].sfr_oii_uncor_err

    residnfgs = ynfgs-xnfgs
    residnfgs_err = sqrt(yerrnfgs^2+xerrnfgs^2)

; make the plot    

    stats = im_stats([resid,residnfgs],/verbose)
    xstr = strtrim(string(stats.median_rej,format='(F12.2)'),2)+' ('+$
      strtrim(string(stats.mean_rej,format='(F12.2)'),2)+$
      '\pm'+strtrim(string(stats.sigma_rej,format='(F12.2)'),2)+')'
;   xstr = strtrim(string(stats.median,format='(F12.2)'),2)+'\pm'+strtrim(string(stats.sig68mean,format='(F12.2)'),2)
    
    xtitle = 'log \psi(H\alpha)'
    ytitle = 'log \psi([O II]) [Observed]'

    xrange = sfrharange
    yrange = xrange
    
    atlas1d_lineplot, x, y, xerr, yerr, plottype=1, postscript=postscript, $
      xtitle='', ytitle=ytitle, xrange=xrange, yrange=yrange, legendtype=0, $
      charsize=charsize, xtickname=replicate(' ',10), position=pos[*,0], $
      xminor=3, yminor=3, $
      xnfgs=xnfgs, ynfgs=ynfgs, xerrnfgs=xerrnfgs, yerrnfgs=yerrnfgs
    djs_oplot, !x.crange, !y.crange, line=0, thick=postthick
    
    legend, '(a)', /left, /top, box=0, charsize=charsize_5, charthick=postthick
    legend, textoidl(xstr), /right, /bottom, box=0, charsize=charsize_3, charthick=postthick

; ##########################################################
; Panel 2
; ##########################################################

; Atlas    
    
    indx = where((atlas_sfrs.sfr_ha gt -900) and (atlas_sfrs.sfr_oii_k98 gt -900),nindx)

    x = atlas_sfrs[indx].sfr_ha
    xerr = atlas_sfrs[indx].sfr_ha_err
    
    y = atlas_sfrs[indx].sfr_oii_k98
    yerr = atlas_sfrs[indx].sfr_oii_k98_err

    resid = y-x
    resid_err = sqrt(yerr^2+xerr^2)

; NFGS    
    
    indxnfgs = where((nfgs_sfrs.sfr_ha gt -900) and (nfgs_sfrs.sfr_oii_k98 gt -900),nindxnfgs)

    xnfgs = nfgs_sfrs[indxnfgs].sfr_ha
    xerrnfgs = nfgs_sfrs[indxnfgs].sfr_ha_err
    
    ynfgs = nfgs_sfrs[indxnfgs].sfr_oii_k98
    yerrnfgs = nfgs_sfrs[indxnfgs].sfr_oii_k98_err

    residnfgs = ynfgs-xnfgs
    residnfgs_err = sqrt(yerrnfgs^2+xerrnfgs^2)

; make the plot    

    stats = im_stats([resid,residnfgs],/verbose)
    xstr = strtrim(string(stats.median_rej,format='(F12.2)'),2)+' ('+$
      strtrim(string(stats.mean_rej,format='(F12.2)'),2)+$
      '\pm'+strtrim(string(stats.sigma_rej,format='(F12.2)'),2)+')'
;   xstr = strtrim(string(stats.median,format='(F12.2)'),2)+'\pm'+strtrim(string(stats.sig68mean,format='(F12.2)'),2)
    
    xtitle = 'log \psi(H\alpha)'
    ytitle = 'log \psi([O II]) [K98]'

    xrange = sfrharange
    yrange = xrange
    
    atlas1d_lineplot, x, y, xerr, yerr, plottype=1, postscript=postscript, $
      xtitle='', ytitle='', xrange=xrange, yrange=yrange, legendtype=0, $
      charsize=charsize, xtickname=replicate(' ',10), position=pos[*,1], $
      ytickname=replicate(' ',10), xminor=3, yminor=3, ystyle=11, /noerase, $
      xnfgs=xnfgs, ynfgs=ynfgs, xerrnfgs=xerrnfgs, yerrnfgs=yerrnfgs
    axis, /yaxis, yrange=yrange, ythick=postthick, ytitle=textoidl(ytitle), $
      charthick=postthick, charsize=charsize, ysty=1
    djs_oplot, !x.crange, !y.crange, line=0, thick=postthick
    
    legend, '(b)', /left, /top, box=0, charsize=charsize_5, charthick=postthick
    legend, textoidl(xstr), /right, /bottom, box=0, charsize=charsize_3, charthick=postthick

; ##########################################################
; Panel 3
; ##########################################################

; Atlas    
    
    indx = where((atlas_sfrs.sfr_ha gt -900) and (atlas_sfrs.sfr_oii_k04 gt -900),nindx)

    x = atlas_sfrs[indx].sfr_ha
    xerr = atlas_sfrs[indx].sfr_ha_err
    
    y = atlas_sfrs[indx].sfr_oii_k04
    yerr = atlas_sfrs[indx].sfr_oii_k04_err

    resid = y-x
    resid_err = sqrt(yerr^2+xerr^2)

; NFGS    
    
    indxnfgs = where((nfgs_sfrs.sfr_ha gt -900) and (nfgs_sfrs.sfr_oii_k04 gt -900),nindxnfgs)

    xnfgs = nfgs_sfrs[indxnfgs].sfr_ha
    xerrnfgs = nfgs_sfrs[indxnfgs].sfr_ha_err
    
    ynfgs = nfgs_sfrs[indxnfgs].sfr_oii_k04
    yerrnfgs = nfgs_sfrs[indxnfgs].sfr_oii_k04_err

    residnfgs = ynfgs-xnfgs
    residnfgs_err = sqrt(yerrnfgs^2+xerrnfgs^2)

; make the plot    

    stats = im_stats([resid,residnfgs],/verbose)
    xstr = strtrim(string(stats.median_rej,format='(F12.2)'),2)+' ('+$
      strtrim(string(stats.mean_rej,format='(F12.2)'),2)+$
      '\pm'+strtrim(string(stats.sigma_rej,format='(F12.2)'),2)+')'
;   xstr = strtrim(string(stats.median,format='(F12.2)'),2)+'\pm'+strtrim(string(stats.sig68mean,format='(F12.2)'),2)
    
    xtitle = 'log \psi(H\alpha)'
    ytitle = 'log \psi([O II]) [K04]'

    xrange = sfrharange
    yrange = xrange
    
    atlas1d_lineplot, x, y, xerr, yerr, plottype=1, postscript=postscript, $
      xtitle='', ytitle=ytitle, xrange=xrange, yrange=yrange, legendtype=0, $
      charsize=charsize, xtickname=replicate(' ',10), position=pos[*,2], $
      xminor=3, yminor=3, /noerase, $
      xnfgs=xnfgs, ynfgs=ynfgs, xerrnfgs=xerrnfgs, yerrnfgs=yerrnfgs
    djs_oplot, !x.crange, !y.crange, line=0, thick=postthick
    
    legend, '(c)', /left, /top, box=0, charsize=charsize_5, charthick=postthick
    legend, textoidl(xstr), /right, /bottom, box=0, charsize=charsize_3, charthick=postthick

; ##########################################################
; Panel 4
; ##########################################################

; Atlas    
    
    indx = where((atlas_sfrs.sfr_ha gt -900) and (atlas_sfrs.sfr_oii_hbhg gt -900),nindx)

    x = atlas_sfrs[indx].sfr_ha
    xerr = atlas_sfrs[indx].sfr_ha_err
    
    y = atlas_sfrs[indx].sfr_oii_hbhg
    yerr = atlas_sfrs[indx].sfr_oii_hbhg_err

    resid = y-x
    resid_err = sqrt(yerr^2+xerr^2)

; NFGS    
    
    indxnfgs = where((nfgs_sfrs.sfr_ha gt -900) and (nfgs_sfrs.sfr_oii_hbhg gt -900),nindxnfgs)

    xnfgs = nfgs_sfrs[indxnfgs].sfr_ha
    xerrnfgs = nfgs_sfrs[indxnfgs].sfr_ha_err
    
    ynfgs = nfgs_sfrs[indxnfgs].sfr_oii_hbhg
    yerrnfgs = nfgs_sfrs[indxnfgs].sfr_oii_hbhg_err

    residnfgs = ynfgs-xnfgs
    residnfgs_err = sqrt(yerrnfgs^2+xerrnfgs^2)

; make the plot    

    stats = im_stats([resid,residnfgs],/verbose)
    xstr = strtrim(string(stats.median_rej,format='(F12.2)'),2)+' ('+$
      strtrim(string(stats.mean_rej,format='(F12.2)'),2)+$
      '\pm'+strtrim(string(stats.sigma_rej,format='(F12.2)'),2)+')'
;   xstr = strtrim(string(stats.median,format='(F12.2)'),2)+'\pm'+strtrim(string(stats.sig68mean,format='(F12.2)'),2)
    
    xtitle = 'log \psi(H\alpha)'
    ytitle = 'log \psi([O II]) [H\beta/H\gamma]'

    xrange = sfrharange
    yrange = xrange
    
    atlas1d_lineplot, x, y, xerr, yerr, plottype=1, postscript=postscript, $
      xtitle='', ytitle='', xrange=xrange, yrange=yrange, legendtype=0, $
      charsize=charsize, xtickname=replicate(' ',10), position=pos[*,3], $
      ytickname=replicate(' ',10), xminor=3, yminor=3, ystyle=11, /noerase, $
      xnfgs=xnfgs, ynfgs=ynfgs, xerrnfgs=xerrnfgs, yerrnfgs=yerrnfgs
    axis, /yaxis, yrange=yrange, ythick=postthick, ytitle=textoidl(ytitle), $
      charthick=postthick, charsize=charsize, ysty=1
    djs_oplot, !x.crange, !y.crange, line=0, thick=postthick
    
    legend, '(d)', /left, /top, box=0, charsize=charsize_5, charthick=postthick
    legend, textoidl(xstr), /right, /bottom, box=0, charsize=charsize_3, charthick=postthick

; ##########################################################
; Panel 5
; ##########################################################

; Atlas    
    
    indx = where((atlas_sfrs.sfr_ha gt -900) and (atlas_sfrs.sfr_oii_best gt -900),nindx)

    x = atlas_sfrs[indx].sfr_ha
    xerr = atlas_sfrs[indx].sfr_ha_err
    
    y = atlas_sfrs[indx].sfr_oii_best
    yerr = atlas_sfrs[indx].sfr_oii_best_err

    resid = y-x
    resid_err = sqrt(yerr^2+xerr^2)

; NFGS    
    
    indxnfgs = where((nfgs_sfrs.sfr_ha gt -900) and (nfgs_sfrs.sfr_oii_best gt -900),nindxnfgs)

    xnfgs = nfgs_sfrs[indxnfgs].sfr_ha
    xerrnfgs = nfgs_sfrs[indxnfgs].sfr_ha_err
    
    ynfgs = nfgs_sfrs[indxnfgs].sfr_oii_best
    yerrnfgs = nfgs_sfrs[indxnfgs].sfr_oii_best_err

    residnfgs = ynfgs-xnfgs
    residnfgs_err = sqrt(yerrnfgs^2+xerrnfgs^2)

; make the plot    

    stats = im_stats([resid,residnfgs],/verbose)
    xstr = strtrim(string(stats.median_rej,format='(F12.2)'),2)+' ('+$
      strtrim(string(stats.mean_rej,format='(F12.2)'),2)+$
      '\pm'+strtrim(string(stats.sigma_rej,format='(F12.2)'),2)+')'
;   xstr = strtrim(string(stats.median,format='(F12.2)'),2)+'\pm'+strtrim(string(stats.sig68mean,format='(F12.2)'),2)
    
    xtitle = 'log \psi(H\alpha)'
    ytitle = 'log \psi([O II]) [This Paper]'

    xrange = sfrharange
    yrange = xrange
    
    atlas1d_lineplot, x, y, xerr, yerr, plottype=1, postscript=postscript, $
      xtitle=xtitle, ytitle=ytitle, xrange=xrange, yrange=yrange, legendtype=0, $
      charsize=charsize, position=pos[*,4], $
      xminor=3, yminor=3, /noerase, $
      xnfgs=xnfgs, ynfgs=ynfgs, xerrnfgs=xerrnfgs, yerrnfgs=yerrnfgs
    djs_oplot, !x.crange, !y.crange, line=0, thick=postthick
    
    legend, '(e)', /left, /top, box=0, charsize=charsize_5, charthick=postthick
    legend, textoidl(xstr), /right, /bottom, box=0, charsize=charsize_3, charthick=postthick

; ##########################################################
; Panel 6
; ##########################################################

; Atlas    
    
    indx = where((atlas_sfrs.sfr_ha gt -900) and (atlas_sfrs.sfr_oii_hahb gt -900),nindx)

    x = atlas_sfrs[indx].sfr_ha
    xerr = atlas_sfrs[indx].sfr_ha_err
    
    y = atlas_sfrs[indx].sfr_oii_hahb
    yerr = atlas_sfrs[indx].sfr_oii_hahb_err

    resid = y-x
    resid_err = sqrt(yerr^2+xerr^2)

; NFGS    
    
    indxnfgs = where((nfgs_sfrs.sfr_ha gt -900) and (nfgs_sfrs.sfr_oii_hahb gt -900),nindxnfgs)

    xnfgs = nfgs_sfrs[indxnfgs].sfr_ha
    xerrnfgs = nfgs_sfrs[indxnfgs].sfr_ha_err
    
    ynfgs = nfgs_sfrs[indxnfgs].sfr_oii_hahb
    yerrnfgs = nfgs_sfrs[indxnfgs].sfr_oii_hahb_err

    residnfgs = ynfgs-xnfgs
    residnfgs_err = sqrt(yerrnfgs^2+xerrnfgs^2)

; make the plot    

    stats = im_stats([resid,residnfgs],/verbose)
    xstr = strtrim(string(stats.median_rej,format='(F12.2)'),2)+' ('+$
      strtrim(string(stats.mean_rej,format='(F12.2)'),2)+$
      '\pm'+strtrim(string(stats.sigma_rej,format='(F12.2)'),2)+')'
;   xstr = strtrim(string(stats.median,format='(F12.2)'),2)+'\pm'+strtrim(string(stats.sig68mean,format='(F12.2)'),2)
    
    xtitle = 'log \psi(H\alpha)'
    ytitle = 'log \psi([O II]) [H\alpha/H\beta]'

    xrange = sfrharange
    yrange = xrange
    
    atlas1d_lineplot, x, y, xerr, yerr, plottype=1, postscript=postscript, $
      xtitle=xtitle, ytitle='', xrange=xrange, yrange=yrange, legendtype=0, $
      charsize=charsize, position=pos[*,5], $
      ytickname=replicate(' ',10), xminor=3, yminor=3, ystyle=11, /noerase, $
      xnfgs=xnfgs, ynfgs=ynfgs, xerrnfgs=xerrnfgs, yerrnfgs=yerrnfgs
    axis, /yaxis, yrange=yrange, ythick=postthick, ytitle=textoidl(ytitle), $
      charthick=postthick, charsize=charsize, ysty=1
    djs_oplot, !x.crange, !y.crange, line=0, thick=postthick
    
    legend, '(f)', /left, /top, box=0, charsize=charsize_5, charthick=postthick
    legend, textoidl(xstr), /right, /bottom, box=0, charsize=charsize_3, charthick=postthick

    im_openclose, postscript=postscript, /close    

; ------------------------------------------------------------
; visualize the various SFR calibrations
; ------------------------------------------------------------
    
    psname = 'visualize_sfrs'
    im_openclose, pspath+psname, postscript=postscript, /encapsulated

;   pagemaker, nx=1, ny=1, ypage=4.9, xspace=0, xmargin=[1.1,0.2], $
;     height=3.08, ymargin=[0.75,1.0], position=pos, /normal
    
    xmax = -39 & xmin = -44.3 & dx = 0.005
    xaxis = findgen((xmax-xmin)/dx+1)*dx+xmin
    xrange = minmax(xaxis)
    
    yrange = [0,0.02]
;   yrange = [0,0.06]

    xtitle = 'log [\psi(H\alpha)/L(\lambda)]'
    ytitle = 'Relative Probability'

    djs_plot, [0], [0], /nodata, xrange=xrange, yrange=yrange, xtitle=xtitle, $
      ytitle=ytitle, charsize=1.6, charthick=postthick, xthick=postthick, $
      ythick=postthick, xsty=3, ysty=1, ytickname=replicate(' ',10)
    
; --------------------------------------------------
; SFR(Ha)/L(Hb)
; --------------------------------------------------

    med = stats.sfr_LHbobs_Abs_median
    sig = stats.sfr_LHbobs_Abs_error

    yaxis = 1.0/sqrt(2*!pi)*sig
    
    yaxis = gauss1(xaxis,[med,sig,1.0])
    yaxis = yaxis/total(yaxis)

    plotit = where(yaxis gt 1E-4)
    djs_oplot, xaxis[plotit], yaxis[plotit], line=1, thick=postthick+3, color='dark green'

; --------------------------------------------------
; SFR(Ha)/L([O II]) observed
; --------------------------------------------------

    med = stats.sfr_Loiiobs_median
    sig = stats.sfr_Loiiobs_error
    
    yaxis = gauss1(xaxis,[med,sig,1.0])
    yaxis = yaxis/total(yaxis)

    plotit = where(yaxis gt 1E-4)
    djs_oplot, xaxis[plotit], yaxis[plotit], line=0, thick=postthick+3, color='red'

; --------------------------------------------------
; SFR(Ha)/L([O II]) corrected
; --------------------------------------------------

    med = stats.sfr_Loiicor_median
    sig = stats.sfr_Loiicor_error
    
    yaxis = gauss1(xaxis,[med,sig,1.0])
    yaxis = yaxis/total(yaxis)

    plotit = where(yaxis gt 1E-4)
    djs_oplot, xaxis[plotit], yaxis[plotit], line=2, thick=postthick+3, color='red'

; --------------------------------------------------
; SFR(Ha)/L([O II]) corrected - low Z
; --------------------------------------------------

    med = stats.sfr_Loiicor_lowZ_median
    sig = stats.sfr_Loiicor_lowZ_error
    
    yaxis = gauss1(xaxis,[med,sig,1.0])
    yaxis = yaxis/total(yaxis)

;   plotit = where(yaxis gt 1E-4)
;   djs_oplot, xaxis[plotit], yaxis[plotit], line=0, thick=postthick, color='blue'

; --------------------------------------------------
; SFR(Ha)/L([O II]) corrected - mid Z
; --------------------------------------------------

    med = stats.sfr_Loiicor_midZ_median
    sig = stats.sfr_Loiicor_midZ_error
    
    yaxis = gauss1(xaxis,[med,sig,1.0])
    yaxis = yaxis/total(yaxis)

;   plotit = where(yaxis gt 1E-4)
;   djs_oplot, xaxis[plotit], yaxis[plotit], line=0, thick=postthick, color='yellow'

; --------------------------------------------------
; SFR(Ha)/L([O II]) corrected - high Z
; --------------------------------------------------

    med = stats.sfr_Loiicor_highZ_median
    sig = stats.sfr_Loiicor_highZ_error
    
    yaxis = gauss1(xaxis,[med,sig,1.0])
    yaxis = yaxis/total(yaxis)

;   plotit = where(yaxis gt 1E-4)
;   djs_oplot, xaxis[plotit], yaxis[plotit], line=0, thick=postthick, color='purple'

; --------------------------------------------------
; SFR(Ha)/L(U) observed
; --------------------------------------------------

    med = stats.sfr_LUobs_median
    sig = stats.sfr_LUobs_error
    
    yaxis = gauss1(xaxis,[med,sig,1.0])
    yaxis = yaxis/total(yaxis)

    plotit = where(yaxis gt 1E-4)
    djs_oplot, xaxis[plotit], yaxis[plotit], line=0, thick=postthick, color='purple'

; --------------------------------------------------
; SFR(Ha)/L(U) corrected
; --------------------------------------------------

    med = stats.sfr_LUcor_median
    sig = stats.sfr_LUcor_error
    
    yaxis = gauss1(xaxis,[med,sig,1.0])
    yaxis = yaxis/total(yaxis)

    plotit = where(yaxis gt 1E-4)
    djs_oplot, xaxis[plotit], yaxis[plotit], line=2, thick=postthick, color='purple'

; --------------------------------------------------
; legend
; --------------------------------------------------

    label = [$
      'H\beta_{obs}', $
      '[O II]_{obs}', $
      '[O II]_{cor}', $
      'L(U)_{obs}', $
      'L(U)_{cor}' $
      ]
    color = [$
      'dark green', $
      'red', $
      'red', $
      'purple', $
      'purple' $
      ]
    line = [$
      1, $
      0, $
      2, $
      0, $
      2]
    
    legend, textoidl(label), /left, /top, box=0, charsize=1.8, $
      charthick=postthick, textcolor=djs_icolor(color), color=djs_icolor(color), $
      thick=postthick, line=line, spacing=2

; overplot a white line to cover up the color

;   djs_oplot, !x.crange, [0,0], thick=postthick, line=0, color='black'
    
    im_openclose, postscript=postscript, /close

; ------------------------------------------------------------
; EW(Ha) vs EW(Hb) - [Integrated + SDSS]
; ------------------------------------------------------------

    psname = 'ewha_vs_ewhb_multipanel'
    im_openclose, pspath+psname, postscript=postscript, /encapsulated, xsize=8.5, ysize=5.2

    pagemaker, nx=2, ny=1, height=3.5, width=3.5, xmargin=[1.2,0.3], $
      ymargin=[0.7,1.0], xspace=0, yspace=0, xpage=8.5, ypage=5.2, $
      position=pos, /normal

    xtitle = 'log EW(H\alpha)  ['+angstrom()+']' 
    ytitle = 'log EW(H\beta)  ['+angstrom()+']'
    
    xrange = ewharange
    yrange = ewhbrange

; Integrated    
    
    lineratio, atlasdust, 'H_ALPHA_EW', '', 'H_BETA_EW', '', x, xerr, $
      y, yerr, index=indx, nindex=nindx, snrcut=snrcut, $
      xsyserr=syserr, ysyserr=syserr

    lineratio, nfgsdust, 'H_ALPHA_EW', '', 'H_BETA_EW', '', xnfgs, xerrnfgs, $
      ynfgs, yerrnfgs, index=indxnfgs, nindex=nindxnfgs, snrcut=snrcut, $
      xsyserr=nfgssyserr, ysyserr=nfgssyserr

    xbig = [x,xnfgs]
    ybig = [y,ynfgs]
    
; the error bar is tiny!    
    
    atlas1d_lineplot, x, y, xerr, yerr, plottype=1, postscript=postscript, $
      xtitle=xtitle, ytitle=ytitle, xrange=xrange, yrange=yrange, legendtype=0, $
      /right, /bottom, position=pos[*,0], charsize=charsize_2, $
      xnfgs=xnfgs, ynfgs=ynfgs, xerrnfgs=xerrnfgs, yerrnfgs=yerrnfgs
    legend, '(a)', /left, /top, box=0, charsize=charsize_2, charthick=postthick

; overplot a bi-sector fit
    
    sixlin, xbig, ybig, a, siga, b, sigb
    coeff = [a[2],b[2]] ; Ordinary Least Squares Bisector
    splog, 'EW(Ha) vs EW(Hb) - Integrated: ', coeff
    
    axis = findgen((3.0-0.5)/0.01+1)*0.01+0.5
    yfit = poly(axis,coeff)    
;   djs_oplot, axis, yfit, thick=postthick, line=0
;   djs_oplot, [0,3], [0,3], thick=postthick-3, line=2

    fitstr = textoidl('log EW(H\beta) = '+string(coeff[0],format='(F5.2)')+' +'+$
      string(coeff[1],format='(F6.3)')+' log EW(H\alpha)')
;   legend, fitstr, /left, /top, box=0, charsize=1.6, charthick=postthick

    splog, 'Statistics about the bisector fit:'
    stats = im_stats(ybig-poly(xbig,coeff),/verbose)
    
; SDSS    
    
    lineratio, sdssdust, 'H_ALPHA_EW', '', 'H_BETA_EW', '', x, xerr, $
      y, yerr, index=indx, nindex=nindx, snrcut=snrcut

    sixlin, x, y, a, siga, b, sigb
    coeff = [a[2],b[2]] ; Ordinary Least Squares Bisector
    splog, 'EW(Ha) vs EW(Hb) - SDSS: ', coeff

    splog, 'Statistics about the bisector fit:'
    stats = im_stats(y-poly(x,coeff),/verbose)
    
    sdss_lineplot, x, y, xerr, yerr, plottype=1, postscript=postscript, $
      xtitle=xtitle, ytitle='', xrange=xrange, yrange=yrange, legendtype=0, $
      /right, /bottom, position=pos[*,1], ytickname=replicate(' ',10), $
      charsize=charsize_2, /noerase
    legend, '(b)', /left, /top, box=0, charsize=charsize_2, charthick=postthick

    im_openclose, postscript=postscript, /close        
    
; ------------------------------------------------------------
; 2-panel - L(IR) vs L(Ha)/L(IR)
; ------------------------------------------------------------
    
    psname = 'LIR_vs_LHa_LIR_multipanel'
    im_openclose, pspath+psname, postscript=postscript, /encapsulated, xsize=8.5, ysize=7.5

    pagemaker, nx=1, ny=2, height=[3.0,3.0], width=5.5, xmargin=[2.0,1.0], $
      ymargin=[0.4,1.1], xspace=0, yspace=0, xpage=8.5, ypage=7.5, $
      position=pos, /normal

    indx = where((atlasnodust.ir_flux gt -900) and (atlasnodust.ebv_hahb_err gt 0.0),nindx)

    x = atlasdust[indx].ir_lum
    xerr = atlasdust[indx].ir_lum_err
    
    lir = atlasdust[indx].ir_flux ; [erg/s/cm2]
    lir_err = sqrt(atlasdust[indx].ir_flux_err^2 + (atlasdust[indx].ir_flux*0.15)^2)

    indxnfgs = where((nfgsdust.ir_flux gt -900) and (nfgsnodust.ebv_hahb_err gt 0.0),nindxnfgs)

    xnfgs = nfgsnodust[indxnfgs].ir_lum
    xerrnfgs = nfgsnodust[indxnfgs].ir_lum_err
    
    lirnfgs = nfgsdust[indxnfgs].ir_flux
    lirnfgs_err = sqrt(nfgsdust[indxnfgs].ir_flux_err^2 + (nfgsdust[indxnfgs].ir_flux*0.15)^2)
    
    xrange = LIRrange
    yrange = haLIRrange

    xtitle = 'log [L(IR)/L'+sunsymbol()+']'
    ytitle = 'log [L(H\alpha)/L(IR)]'

    lhalir = alog10(4.5D-44/7.9D-42)

    lir_bell = findgen((xrange[1]-xrange[0])/0.01+1)*0.01+xrange[0]
    sfrir_bell = lir_bell*0.0D
    hi = where(lir_bell gt 11.0,comp=lo)
    sfrir_bell[hi] = 1.57D-10*(1+sqrt(1E9/10^lir_bell[hi]))/lsun
    sfrir_bell[lo] = 1.17D-10*(1+sqrt(1E9/10^lir_bell[lo]))/lsun
    lratio_bell = alog10(sfrir_bell) - alog10(7.9D-42)

; ##########################################################
; Panel 1: No reddening correction
; ##########################################################

    ha = atlasdust[indx].h_alpha[0]
    ha_err = atlasdust[indx].h_alpha[1]
    
    y = alog10(ha/lir)
    yerr = im_compute_error(ha,ha_err,lir,lir_err,/log)

    hanfgs = nfgsdust[indxnfgs].h_alpha[0]
    hanfgs_err = nfgsdust[indxnfgs].h_alpha[1]
    
    ynfgs = alog10(hanfgs/lirnfgs)
    yerrnfgs = im_compute_error(hanfgs,hanfgs_err,lirnfgs,lirnfgs_err,/log)

    stats = im_stats([y,ynfgs],/verbose)
    xstr = strtrim(string(stats.median_rej,format='(F12.2)'),2)+' ('+$
      strtrim(string(stats.mean_rej,format='(F12.2)'),2)+$
      '\pm'+strtrim(string(stats.sigma_rej,format='(F12.2)'),2)+')'
;   xstr = strtrim(string(stats.median,format='(F12.2)'),2)+'\pm'+strtrim(string(stats.sig68mean,format='(F12.2)'),2)

    ytitle = 'log [L(H\alpha)_{obs}/L(IR)]'

    atlas1d_lineplot, x, y, xerr, yerr, plottype=3, postscript=postscript, $
      xtitle='', ytitle=ytitle, xrange=xrange, yrange=yrange, legendtype=0, $
      charsize=charsize_2, position=pos[*,0], xnfgs=xnfgs, ynfgs=ynfgs, $
      xerrnfgs=xerrnfgs, yerrnfgs=yerrnfgs, xtickname=replicate(' ',10), yminor=3, $
      yfrac=8
    djs_oplot, !x.crange, lhalir*[1,1], line=2, thick=postthick
;   djs_oplot, lir_bell, lratio_bell, line=2, thick=postthick

    legend, '(a)', /left, /top, box=0, charsize=charsize_2, charthick=postthick

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

    y = alog10(ha/lir)
    yerr = im_compute_error(ha,ha_err,lir,lir_err,/log)

    hanfgs = nfgsnodust[indxnfgs].h_alpha[0]
    hanfgs_err = nfgsnodust[indxnfgs].h_alpha[1]

    ynfgs = alog10(hanfgs/lirnfgs)
    yerrnfgs = im_compute_error(hanfgs,hanfgs_err,lirnfgs,lirnfgs_err,/log)

    stats = im_stats([y,ynfgs],/verbose,/no_head)
    xstr = strtrim(string(stats.median_rej,format='(F12.2)'),2)+' ('+$
      strtrim(string(stats.mean_rej,format='(F12.2)'),2)+$
      '\pm'+strtrim(string(stats.sigma_rej,format='(F12.2)'),2)+')'
;   xstr = strtrim(string(stats.median,format='(F12.2)'),2)+'\pm'+strtrim(string(stats.sig68mean,format='(F12.2)'),2)

    ytitle = 'log [L(H\alpha)_{cor}/L(IR)]'

    atlas1d_lineplot, x, y, xerr, yerr, plottype=3, postscript=postscript, $
      xtitle=xtitle, ytitle=ytitle, xrange=xrange, yrange=yrange, legendtype=0, $
      charsize=charsize_2, position=pos[*,1], /noerase, xnfgs=xnfgs, ynfgs=ynfgs, $
      xerrnfgs=xerrnfgs, yerrnfgs=yerrnfgs, yminor=3
      
    djs_oplot, !x.crange, lhalir*[1,1], line=2, thick=postthick
;   djs_oplot, lir_bell, lratio_bell, line=2, thick=postthick

    legend, '(b)', /left, /top, box=0, charsize=charsize_2, charthick=postthick

    legend, textoidl(xstr), /right, /top, box=0, charsize=charsize_3, charthick=postthick
;   label = ['Individual A(H\alpha)']
;   legend, textoidl(label), /left, /bottom, box=0, charsize=charsize_3, $
;     charthick=postthick

    w = where([x,xnfgs] ge 9.0 and [x,xnfgs] le 11.0)
;   junk = im_stats(([y,ynfgs])[w],/verbose)

    im_openclose, postscript=postscript, /close

; ------------------------------------------------------------
; L(B) versus SFR(Ha)/L(Hb)_obs
; ------------------------------------------------------------
    
    psname = 'LB_vs_sfr_ha_LHb_obs'
    im_openclose, pspath+psname, postscript=postscript, xsize=8.5, ysize=8.0, /encapsulated

    pagemaker, nx=1, ny=1, xspace=0, yspace=0.0, width=6.0, height=6.0, $
      xmargin=[2.0,0.5], ymargin=[0.9,1.1], xpage=8.5, ypage=8.0, $
      position=pos, /normal

; Atlas    
    
    indx = where((atlasdust.b_lum_obs gt -900) and (atlasnodust.ebv_hahb_err gt 0.0) and $
      (atlasdust.h_beta_ew_uncor[0] ge ewcut),nindx)

    x = atlasdust[indx].b_lum_obs
    xerr = atlasdust[indx].b_lum_obs_err
    xabs = atlasdust[indx].m_b_obs

    hb = atlasdust[indx].h_beta_lum[0]
    hb_err = atlasdust[indx].h_beta_lum[1]

    ha = atlasnodust[indx].h_alpha_lum[0]
    ha_err = atlasnodust[indx].h_alpha_lum[1]

    y = ha + loghaconst - hb
    yerr = sqrt(ha_err^2 + hb_err^2)

; linearize    
    
;   y = 10^y
;   yerr = alog(10)*yerr*y
    
; NFGS    
    
    indxnfgs = where((nfgsnodust.b_lum_obs gt -900.0) and (nfgsnodust.ebv_hahb_err gt 0.0) and $
      (nfgsdust.h_beta_ew_uncor[0] ge ewcut),nindxnfgs)

    xnfgs = nfgsdust[indxnfgs].b_lum_obs
    xerrnfgs = nfgsdust[indxnfgs].b_lum_obs_err
    
    hbnfgs = nfgsdust[indxnfgs].h_beta_lum[0]
    hbnfgs_err = nfgsdust[indxnfgs].h_beta_lum[1]

    hanfgs = nfgsnodust[indxnfgs].h_alpha_lum[0]
    hanfgs_err = nfgsnodust[indxnfgs].h_alpha_lum[1]

    ynfgs = hanfgs + loghaconst - hbnfgs
    yerrnfgs = sqrt(hanfgs_err^2 + hbnfgs_err^2)

; linearize    
    
;   ynfgs = 10^ynfgs
;   yerrnfgs = alog(10)*yerrnfgs*ynfgs
    
; make the plot

    xtitle = 'log L(B) [L(B)'+sunsymbol()+']'
;   ytitle = '\psi(H\alpha) / L(H\beta)_{obs} [10^{-40} M'+sunsymbol()+' yr^{-1}/erg s^{-1}]'
    ytitle = 'log [\psi(H\alpha) / L(H\beta)_{obs}]  (M'+sunsymbol()+' yr^{-1}/erg s^{-1})'


    scale = 1.0
;   scale = 1D40
    
    xrange = LBrange
    yrange = scale*[1.9D-41,3D-40]
    yrange = alog10([1.9D-41,3D-40])
    
    atlas1d_lineplot, x, scale*y, xerr, scale*yerr, plottype=1, postscript=postscript, $
      xtitle=xtitle, ytitle=ytitle, xrange=xrange, yrange=yrange, legendtype=0, $
      /right, /top, yfrac=1.5, xfrac=1.5, charsize=1.8, xstyle=11, ylog=(scale ne 1.0), $
      position=pos[*,0], xnfgs=xnfgs, ynfgs=scale*ynfgs, xerrnfgs=xerrnfgs, yerrnfgs=scale*yerrnfgs
    axis, /xaxis, xrange=interpol(xabs,x,!x.crange), xthick=postthick, xstyle=1, $
      charsize=1.8, charthick=postthick
    xyouts, (pos[2]-pos[0])/2.0+pos[0], pos[3]+0.08*(pos[3]-pos[1]), $
      textoidl('M_{B} [mag]'), align=0.5, charsize=1.8, charthick=postthick, /normal
    
; flag the LIRGS and ULIRGS

    lirgs = where(atlasdust[indx].ir_lum gt 11.0,comp=notlirgs)
    lirgsnfgs = where(nfgsdust[indxnfgs].ir_lum gt 11.0,notlirgsnfgs)

    plotsym, 0, 1.6, thick=(postthick-3L)>1L
    djs_oplot, x[lirgs], scale*y[lirgs], ps=8
    djs_oplot, xnfgs[lirgsnfgs], scale*ynfgs[lirgsnfgs], ps=8
    
    xbig = [x[notlirgs],xnfgs[notlirgsnfgs]]
    ybig = [y[notlirgs],ynfgs[notlirgsnfgs]]
    ybigerr = [yerr[notlirgs],yerrnfgs[notlirgsnfgs]]

; overlay the running median    

    running = im_medxbin(xbig,ybig,0.3,minpts=5)
;   oploterror, running.binctr, scale*running.medy, replicate(running.binsz/2.0,running.nbins), $
;     scale*running.sigy, ps=-3, thick=(postthick-2L)>2L, errthick=(postthick-2L)>2L, /nohat, $
;     line=0, errstyle=0

; do the fit in log-linear space; constrain the fit to asymptote to
; the median ratio below 10^9 L_sun, and only fit above 10^9 L_sun;
; subtract the first luminosity element to ensure that the intercept
; of the function equals the median below 10^9 L_sun

    LBcut = 8.77 ; 8.7 ; 9.0 is M_B = -17.07 and 8.7 is M_B = -16.3; 8.57 is -16
    indx = where(running.binctr gt LBcut)
    yesfit = where(xbig gt LBcut,comp=nofit)

;   med = median(ybig[nofit])
    
    xfitaxis = findgen((max(running.binctr[indx])-LBcut)/0.05+1L)*0.05 + LBcut
    yfitaxis = interpol(running.medy[indx],running.binctr[indx],xfitaxis)

    sixlin, xfitaxis, scale*yfitaxis, a, siga, b, sigb
    coeff = [a[2],b[2]]/scale
    splog, 'log L(B)/L(B)_sun vs log [SFR(Ha)/L(Hb)_obs] coefficients: ', coeff

;   order = 1L
;   coeff_guess = [scale*med,1.0D]
;   coeff_fixed = [1B,0B]
;
;   coeff = im_polyfit(xbig[yesfit]-xstart,scale*ybig[yesfit],order,$;yerr=scale*ybigerr[yesfit],$
;     coeff_guess=coeff_guess,coeff_fixed=coeff_fixed,coeff_err=coeff_err,$
;     yfit=polyfit,chi2=chi2)
;
;   coeff = coeff/scale & coeff_err = coeff_err/scale & polyfit = polyfit/scale
;
;   xaxis = findgen((11.2-7.3)/0.01+1L)*0.01+7.3
;   yfit = xaxis*0.0
;   no = where(xaxis le LBcut,comp=yes)
;
;   yfit[no] = med
;   yfit[yes] = interpol(polyfit,xbig[yesfit],xaxis[yes])
    
    xaxis = findgen((max(xbig)-LBcut)/0.01+1L)*0.01+LBcut
    yfit = poly(xaxis,coeff)

    djs_oplot, xaxis, scale*yfit, thick=postthick+2, line=0, color='dark green'

    im_openclose, postscript=postscript, /close    

; ------------------------------------------------------------
; L(B) versus SFR(Ha)/L([O II])_obs
; ------------------------------------------------------------
    
    psname = 'LB_vs_sfr_ha_Loii_obs'
    im_openclose, pspath+psname, postscript=postscript, xsize=8.5, ysize=8.0, /encapsulated

    pagemaker, nx=1, ny=1, xspace=0, yspace=0.0, width=6.0, height=6.0, $
      xmargin=[2.0,0.5], ymargin=[0.9,1.1], xpage=8.5, ypage=8.0, $
      position=pos, /normal

; Atlas    
    
    indx = where((atlasdust.oii_3727[0]/atlasdust.oii_3727[1] gt snrcut) and $
      (atlasnodust.b_lum_obs gt -900.0) and (atlasnodust.ehbha_err gt 0.0),nindx)

    x = atlasdust[indx].b_lum_obs
    xerr = atlasdust[indx].b_lum_obs_err
    xabs = atlasdust[indx].m_b_obs

    oii = atlasdust[indx].oii_3727_lum[0]
    oii_err = atlasdust[indx].oii_3727_lum[1]

    ha = atlasnodust[indx].h_alpha_lum[0]
    ha_err = atlasnodust[indx].h_alpha_lum[1]

    y = ha + loghaconst - oii
    yerr = sqrt(ha_err^2 + oii_err^2)

; linearize    
    
;   y = 10^y
;   yerr = alog(10)*yerr*y
    
; NFGS    
    
    indxnfgs = where((nfgsdust.oii_3727[0]/nfgsdust.oii_3727[1] gt snrcut) and $
      (nfgsnodust.b_lum_obs gt -900.0) and (nfgsnodust.ehbha_err gt 0.0),nindxnfgs)

    xnfgs = nfgsdust[indxnfgs].b_lum_obs
    xerrnfgs = nfgsdust[indxnfgs].b_lum_obs_err
    
    oiinfgs = nfgsdust[indxnfgs].oii_3727_lum[0]
    oiinfgs_err = nfgsdust[indxnfgs].oii_3727_lum[1]

    hanfgs = nfgsnodust[indxnfgs].h_alpha_lum[0]
    hanfgs_err = nfgsnodust[indxnfgs].h_alpha_lum[1]

    ynfgs = hanfgs + loghaconst - oiinfgs
    yerrnfgs = sqrt(hanfgs_err^2 + oiinfgs_err^2)

; linearize    
    
;   ynfgs = 10^ynfgs
;   yerrnfgs = alog(10)*yerrnfgs*ynfgs
    
; make the plot

    xtitle = 'log L(B) [L(B)'+sunsymbol()+']'
    ytitle = 'log [\psi(H\alpha) / L([O II])_{obs}] (M'+sunsymbol()+' yr^{-1}/erg s^{-1})'
;   ytitle = '\psi(H\alpha) / L([O II])_{obs}  [10^{-40} M'+sunsymbol()+' yr^{-1}/erg s^{-1}]'

    scale = 1.0
;   scale = 1D40
    
    xrange = LBrange
;   yrange = scale*[2.0D-42,7D-40]
    yrange = alog10([2.0D-42,7D-40])
    
    atlas1d_lineplot, x, scale*y, xerr, scale*yerr, plottype=1, postscript=postscript, $
      xtitle=xtitle, ytitle=ytitle, xrange=xrange, yrange=yrange, legendtype=0, $
      /right, /top, yfrac=1.5, xfrac=1.5, charsize=1.8, xstyle=11, ylog=(scale ne 1.0), $
      position=pos[*,0], xnfgs=xnfgs, ynfgs=scale*ynfgs, xerrnfgs=xerrnfgs, yerrnfgs=scale*yerrnfgs
    axis, /xaxis, xrange=interpol(xabs,x,!x.crange), xthick=postthick, xstyle=1, $
      charsize=1.8, charthick=postthick
    xyouts, (pos[2]-pos[0])/2.0+pos[0], pos[3]+0.08*(pos[3]-pos[1]), $
      textoidl('M_{B} [mag]'), align=0.5, charsize=1.8, charthick=postthick, /normal

; flag the LIRGS and ULIRGS

    lirgs = where(atlasdust[indx].ir_lum gt 11.0,comp=notlirgs)
    lirgsnfgs = where(nfgsdust[indxnfgs].ir_lum gt 11.0,notlirgsnfgs)

    plotsym, 0, 1.6, thick=(postthick-3L)>1L
    djs_oplot, x[lirgs], scale*y[lirgs], ps=8
    djs_oplot, xnfgs[lirgsnfgs], scale*ynfgs[lirgsnfgs], ps=8
    
    xbig = [x,xnfgs]
    ybig = [y,ynfgs]
    ybigerr = [yerr,yerrnfgs]

;   xbig = [x[notlirgs],xnfgs[notlirgsnfgs]]
;   ybig = [y[notlirgs],ynfgs[notlirgsnfgs]]
;   ybigerr = [yerr[notlirgs],yerrnfgs[notlirgsnfgs]]

; overlay the running median    

    running = im_medxbin(xbig,ybig,0.45,minpts=15,maxx=10.8)
;   oploterror, running.binctr, scale*running.medy, replicate(running.binsz/2.0,running.nbins), $
;     scale*running.sigyup, ps=-3, thick=(postthick-2L)>2L, errthick=(postthick-2L)>2L, /nohat, $
;     line=0, errstyle=0, /hibar
;   oploterror, running.binctr, scale*running.medy, replicate(running.binsz/2.0,running.nbins), $
;     scale*running.sigylo, ps=-3, thick=(postthick-2L)>2L, errthick=(postthick-2L)>2L, /nohat, $
;     line=0, errstyle=0, /lobar
;   oplot, running.binctr, scale*running.medy, ps=-3, thick=(postthick-2L)>2L, line=0
;   oploterror, running.binctr, scale*running.medy, replicate(running.binsz/2.0,running.nbins), $
;     0.0*scale*running.sigy, ps=-3, thick=(postthick-2L)>2L, $
;     errthick=(postthick-2L)>2L, /nohat, line=0, errstyle=0
    
; do the fit in log-linear space to the running median; only fit above
; 10^8 L_sun; over-interpolate the running median; L_B = 8.77 is M_B =
; -16.5 

    LBcut = 8.77 ; min(xbig)
    indx = where(running.binctr gt LBcut)
    yesfit = where(xbig gt LBcut,comp=nofit)

    xfitaxis = findgen((max(running.binctr[indx])-LBcut)/0.05+1L)*0.05 + LBcut
    yfitaxis = interpol(running.medy[indx],running.binctr[indx],xfitaxis)

    sixlin, xfitaxis, scale*yfitaxis, a, siga, b, sigb
    coeff = [a[2],b[2]]/scale
    splog, 'log L(B)/L(B)_sun vs log [SFR(Ha)/L(OII)_obs] coefficients: ', coeff

;   sixlin, xbig[yesfit], scale*ybig[yesfit], a, siga, b, sigb
;   coeff = [a[2],b[2]]/scale

;   xfitaxis = findgen((max(running.binctr[indx])-LBcut)/0.05+1L)*0.05 + LBcut
;   yfitaxis = interpol(running.medy[indx],running.binctr[indx],xfitaxis)
;   coeff = im_polyfit(xfitaxis-LBcut,scale*yfitaxis,order,yfit=polyfit)
 
;   order = 1L
;   coeff = robust_poly_fit(xbig-LBcut,scale*ybig,order,polyfit)
;   coeff = im_polyfit(xbig-LBcut,scale*ybig,order,yfit=polyfit)
;   coeff = coeff/scale & polyfit = polyfit/scale
 
    xaxis = findgen((max(xbig)-LBcut)/0.01+1L)*0.01+LBcut
;   yfit = poly(xaxis-LBcut,coeff)
    yfit = poly(xaxis,coeff)

    djs_oplot, xaxis, scale*yfit, thick=postthick+2, line=0, color='dark green'

    im_openclose, postscript=postscript, /close    

; ------------------------------------------------------------
; SFR([O II]) [Kewley] vs SFR(Ha)
; ------------------------------------------------------------

    psname = 'sfr_compare_oii_kewley04_ha'
    im_openclose, pspath+psname, postscript=postscript, /encapsulated

; integrated sample    
    
    indx = where((atlasdust.oii_3727_lum[0] gt -900) and (atlasnodust.ebv_hahb_err gt 0.0) and $
      (atlasdust.oii_3727[0]/atlasdust.oii_3727[1] gt snrcut) and $
      (atlasdust.oiii_5007[0]/atlasdust.oiii_5007[1] gt snrcut) and $
      (atlasdust.h_beta[0]/atlasdust.h_beta[1] gt snrcut),nindx)

    indxnfgs = where((nfgsdust.oii_3727_lum[0] gt -900) and (nfgsnodust.ebv_hahb_err gt 0.0) and $
      (nfgsdust.oii_3727[0]/nfgsdust.oii_3727[1] gt snrcut) and $
      (nfgsdust.oiii_5007[0]/nfgsdust.oiii_5007[1] gt snrcut) and $
      (nfgsdust.h_beta[0]/nfgsdust.h_beta[1] gt snrcut),nindxnfgs)

    indxsdss = where((sdssdust.oii_3727_lum[0] gt -900) and (sdssnodust.ebv_hahb_err gt 0.0) and $
      (sdssdust.oii_3727[0]/sdssdust.oii_3727[1] gt snrcut) and $
      (sdssdust.oiii_5007[0]/sdssdust.oiii_5007[1] gt snrcut) and $
      (sdssdust.h_beta[0]/sdssdust.h_beta[1] gt snrcut),nindxsdss)

    oii3727  = [ [atlasdust[indx].oii_3727], [nfgsdust[indxnfgs].oii_3727] ]
    oiii4959 = [ [atlasdust[indx].oiii_4959], [nfgsdust[indxnfgs].oiii_4959] ]
    oiii5007 = [ [atlasdust[indx].oiii_5007], [nfgsdust[indxnfgs].oiii_5007] ]
    halpha   = [ [atlasdust[indx].h_alpha], [nfgsdust[indxnfgs].h_alpha] ]

    lumoii_o = lsun*10^[atlasdust[indx].oii_3727_lum[0],nfgsdust[indxnfgs].oii_3727_lum[0]]
    sfr_ha = [atlasnodust[indx].sfr_h_alpha[0],nfgsnodust[indxnfgs].sfr_h_alpha[0]]
;   sfr_ha = 10^[atlasnodust[indx].sfr_h_alpha[0],nfgsnodust[indxnfgs].sfr_h_alpha[0]]

    ebv_true = [atlasnodust[indx].ebv_hahb,nfgsnodust[indxnfgs].ebv_hahb]
    logoh_true = [atlasnodust[indx].zstrong_12oh_oiiinii_pettini,nfgsnodust[indxnfgs].zstrong_12oh_oiiinii_pettini]

; SDSS

    indxsdss = where((sdssdust.oii_3727_lum[0] gt -900) and (sdssnodust.ebv_hahb_err gt 0.0) and $
      (sdssdust.oii_3727[0]/sdssdust.oii_3727[1] gt snrcut) and $
      (sdssdust.oiii_5007[0]/sdssdust.oiii_5007[1] gt snrcut) and $
      (sdssdust.h_beta[0]/sdssdust.h_beta[1] gt snrcut),nindxsdss)

    oii3727_sdss  = sdssdust[indxsdss].oii_3727
    oiii4959_sdss = sdssdust[indxsdss].oiii_4959
    oiii5007_sdss = sdssdust[indxsdss].oiii_5007
    halpha_sdss   = sdssdust[indxsdss].h_alpha

    lumoii_o_sdss = lsun*10^sdssdust[indxsdss].oii_3727_lum[0]
    sfr_ha_sdss = sdssnodust[indxsdss].sfr_h_alpha[0]
;   sfr_ha_sdss = 10^sdssnodust[indxsdss].sfr_h_alpha[0]

; first, compute the intrinsic luminosity from her eq. (18)

    lumoii_i = 3.11D-20*lumoii_o^1.495
    lumoii_i_sdss = 3.11D-20*lumoii_o_sdss^1.495

; second, estimate E(B-V) from eq. (16)

    ebv = (0.174 * alog10(lumoii_i) - 6.84) > 0.0
    ebv_sdss = (0.174 * alog10(lumoii_i_sdss) - 6.84) > 0.0

;   plot, ebv_true, ebv, xrange=[-0.1,1], yrange=[-0.1,1], xsty=3, ysty=3, ps=4
;   oplot, !x.crange, !y.crange
    
; third, correct R23 for the E(B-V) estimated above

; integrated sample    
    
    line = {$
      linename:       ['OII_3727','H_BETA','OIII_4959','OIII_5007','H_ALPHA'], $
      oii_3727:       [0.0D,0.0D], $
      oii_3727_wave:  atlasdust[0].oii_3727_wave, $
      h_beta:         [0.0D,0.0D], $
      h_beta_wave:    atlasdust[0].h_beta_wave, $
      oiii_4959:      [0.0D,0.0D], $
      oiii_4959_wave: atlasdust[0].oiii_4959_wave, $
      oiii_5007:      [0.0D,0.0D], $
      oiii_5007_wave: atlasdust[0].oiii_5007_wave, $
      h_alpha:        [0.0D,0.0D], $
      h_alpha_wave:   atlasdust[0].h_alpha_wave}
        
    line = replicate(line,nindx+nindxnfgs)

    hbeta = halpha
    hbeta[0,*] = halpha[0,*]/(return_tbalmer(/HaHb)*10^(ebv/2.18417))
    
    line.oii_3727  = oii3727
    line.h_beta    = hbeta
    line.oiii_4959 = oiii4959
    line.oiii_5007 = oiii5007
    line.h_alpha   = halpha

    linenodust = iunred_linedust(line,snrcut=3.0,/silent,/nopropagate)

; SDSS    
    
    line_sdss = {$
      linename:       ['OII_3727','H_BETA','OIII_4959','OIII_5007','H_ALPHA'], $
      oii_3727:       [0.0D,0.0D], $
      oii_3727_wave:  sdssdust[0].oii_3727_wave, $
      h_beta:         [0.0D,0.0D], $
      h_beta_wave:    sdssdust[0].h_beta_wave, $
      oiii_4959:      [0.0D,0.0D], $
      oiii_4959_wave: sdssdust[0].oiii_4959_wave, $
      oiii_5007:      [0.0D,0.0D], $
      oiii_5007_wave: sdssdust[0].oiii_5007_wave, $
      h_alpha:        [0.0D,0.0D], $
      h_alpha_wave:   sdssdust[0].h_alpha_wave}
        
    line_sdss = replicate(line_sdss,nindxsdss)

    hbeta_sdss = halpha_sdss
    hbeta_sdss[0,*] = halpha_sdss[0,*]/(return_tbalmer(/HaHb)*10^(ebv_sdss/2.18417))
    
    line_sdss.oii_3727  = oii3727_sdss
    line_sdss.h_beta    = hbeta_sdss
    line_sdss.oiii_4959 = oiii4959_sdss
    line_sdss.oiii_5007 = oiii5007_sdss
    line_sdss.h_alpha   = halpha_sdss

    linenodust_sdss = iunred_linedust(line_sdss,snrcut=3.0,/silent,/nopropagate)
    
; fourth, use eq. (11) to estimate the abundance from the Z94 diagnostic

    abund = im_abundance(linenodust,snrcut=3.0,/notstrict)
    logoh = abund.zstrong_12oh_zkh94

    abund_sdss = im_abundance(linenodust_sdss,snrcut=3.0,/notstrict)
    logoh_sdss = abund_sdss.zstrong_12oh_zkh94

; fifth, use eq. (10) to compute the reddening and abundance corrected
; SFR([O II],Z) estimate

    sfr_oii = alog10(7.9D-42*lumoii_i/(-1.75*logoh+16.73))
    sfr_oii_theory = alog10(7.9D-42*lumoii_i/(-1857.24+612.693*logoh-67.0264*logoh^2+2.43209*logoh^3))
    sfr_ratio = sfr_oii - sfr_ha

    stats = im_stats(sfr_ratio)
    xstr = strtrim(string(stats.median_rej,format='(F12.2)'),2)+' ('+$
      strtrim(string(stats.mean_rej,format='(F12.2)'),2)+$
      '\pm'+strtrim(string(stats.sigma_rej,format='(F12.2)'),2)+')'
;   xstr = strtrim(string(stats.median,format='(F12.2)'),2)+'\pm'+strtrim(string(stats.sig68mean,format='(F12.2)'),2)

    sfr_oii_sdss = alog10(7.9D-42*lumoii_i_sdss/(-1.75*logoh_sdss+16.73)) ; [M_sun/yr]
    sfr_ratio_sdss = sfr_oii_sdss - sfr_ha_sdss

    stats_sdss = im_stats(sfr_ratio_sdss)
    xstr_sdss = strtrim(string(stats_sdss.median,format='(F12.2)'),2)+'\pm'+$
      strtrim(string(stats_sdss.sig68mean,format='(F12.2)'),2)
    
    xrange = [min(sfr_ha)<min(sfr_oii),max(sfr_ha)>max(sfr_oii)]
    yrange = xrange

    xrange2 = xrange
    yrange2 = [-1.4,1.4] ; max(abs(sfr_ratio))*[-1.1,1.1]
    
    xtitle = 'log [\psi(H\alpha)]'
    ytitle = 'log {\psi([O II])} [Kewley]'
    ytitle2 = 'log {\psi(H\alpha)/\psi([O II])}'
    
    pagemaker, nx=1, ny=2, yspace=0, xspace=0, xmargin=[1.9,0.5], $
      ymargin=[0.3,1.3], position=pos, height=[6.0,3.4], /normal

    sdss_lineplot, sfr_ha_sdss, sfr_oii_sdss, sfr_ha_sdss*0.0, sfr_oii_sdss*0.0, $
      plottype=3, postscript=postscript, ytitle=ytitle, xrange=xrange, yrange=yrange, $
      legendtype=0, position=pos[*,0], charsize=1.6, xtickname=replicate(' ',10), $
      xatlas=sfr_ha, yatlas=sfr_oii, xerratlas=sfr_ha*0.0, yerratlas=sfr_oii*0.0, $
      atlaspsize=0.5
    djs_oplot, !x.crange, !y.crange, line=0, thick=postthick
    
    sdss_lineplot, sfr_ha_sdss, sfr_ratio_sdss, sfr_ha_sdss*0.0, sfr_ratio_sdss*0.0, $
      plottype=3, postscript=postscript, ytitle=ytitle2, xrange=xrange2, yrange=yrange2, $
      legendtype=0, position=pos[*,1], charsize=1.6, xtitle=xtitle, $
      xatlas=sfr_ha, yatlas=sfr_ratio, xerratlas=sfr_ha*0.0, yerratlas=sfr_ratio*0.0, $
      atlaspsize=0.5, /noerase
    djs_oplot, !x.crange, [0,0], line=0, thick=postthick

    legend, textoidl(['SDSS Sample '+xstr_sdss,'Integrated Sample '+xstr]), $
      /right, /top, box=0, charsize=charsize_3, charthick=postthick, spacing=1.5

;   plotsym, 0, 1, /fill
;   djs_plot, sfr_ha, sfr_oii, ps=8, xsty=3, ysty=3, xrange=xrange, $
;     yrange=yrange, xthick=postthick, ythick=postthick, charsize=1.7, $
;     charthick=postthick, xtitle='', ytitle=ytitle, /xlog, /ylog, $
;     position=pos[*,0], xtickname=replicate(' ',10), color='dark green'
;   djs_oplot, 10^!x.crange, 10^!y.crange, line=0, thick=postthick
;
;   djs_plot, sfr_ha, sfr_ratio, ps=8, xsty=3, ysty=3, xrange=xrange2, $
;     yrange=yrange2, xthick=postthick, ythick=postthick, charsize=1.7, $
;     charthick=postthick, xtitle=xtitle, ytitle=ytitle2, /xlog, $
;     position=pos[*,1], /noerase, color='dark green'
;   djs_oplot, 10^!x.crange, [0,0], line=0, thick=postthick
;   legend, textoidl(xstr), /left, /top, box=0, charsize=charsize_3, charthick=postthick

    im_openclose, postscript=postscript, /close    

; ------------------------------------------------------------
; SFR([O II]) [Kennicutt] vs SFR(Ha)
; ------------------------------------------------------------

    psname = 'sfr_compare_oii_kenn98_ha'
    im_openclose, pspath+psname, postscript=postscript, /encapsulated

    indx = where((atlasdust.oii_3727_lum[0] gt -900) and (atlasnodust.sfr_h_alpha[0] gt -900.0),nindx)
    indxnfgs = where((nfgsdust.oii_3727_lum[0] gt -900) and (nfgsnodust.sfr_h_alpha[0] gt -900.0),nindxnfgs)
    indxsdss = where((sdssdust.oii_3727_lum[0] gt -900) and (sdssnodust.sfr_h_alpha[0] gt -900.0),nindxsdss)

    lum_oii = lsun*10^[atlasdust[indx].oii_3727_lum[0],nfgsdust[indxnfgs].oii_3727_lum[0]]
    sfr_oii = alog10(1.4D-41*lum_oii)
    sfr_ha = [atlasnodust[indx].sfr_h_alpha[0],nfgsnodust[indxnfgs].sfr_h_alpha[0]]

    lum_oii_sdss = lsun*10^sdssdust[indxsdss].oii_3727_lum[0]
    sfr_oii_sdss = alog10(1.4D-41*lum_oii_sdss)
    sfr_ha_sdss = sdssnodust[indxsdss].sfr_h_alpha[0]

    sfr_ratio = sfr_oii - sfr_ha
    stats = im_stats(sfr_ratio)
    xstr = strtrim(string(stats.median_rej,format='(F12.2)'),2)+' ('+$
      strtrim(string(stats.mean_rej,format='(F12.2)'),2)+$
      '\pm'+strtrim(string(stats.sigma_rej,format='(F12.2)'),2)+')'
;   xstr = strtrim(string(stats.median,format='(F12.2)'),2)+'\pm'+strtrim(string(stats.sig68mean,format='(F12.2)'),2)
    
    sfr_ratio_sdss = sfr_oii_sdss - sfr_ha_sdss
    stats_sdss = im_stats(sfr_ratio_sdss)
    xstr_sdss = strtrim(string(stats_sdss.median,format='(F12.2)'),2)+'\pm'+$
      strtrim(string(stats_sdss.sig68mean,format='(F12.2)'),2)
    
    xrange = [min(sfr_ha)<min(sfr_oii),max(sfr_ha)>max(sfr_oii)]
    yrange = xrange

    xrange2 = xrange
    yrange2 = max(abs(sfr_ratio))*[-1.1,1.1]
    
    xtitle = 'log [\psi(H\alpha)]'
    ytitle = 'log {\psi([O II])} [Kennicutt 1998]'
    ytitle2 = 'log {\psi(H\alpha)/\psi([O II])}'
    
    pagemaker, nx=1, ny=2, yspace=0, xspace=0, xmargin=[1.9,0.5], $
      ymargin=[0.3,1.3], position=pos, height=[6.0,3.4], /normal

    sdss_lineplot, sfr_ha_sdss, sfr_oii_sdss, sfr_ha_sdss*0.0, sfr_oii_sdss*0.0, $
      plottype=3, postscript=postscript, ytitle=ytitle, xrange=xrange, yrange=yrange, $
      legendtype=0, position=pos[*,0], charsize=1.6, xtickname=replicate(' ',10), $
      xatlas=sfr_ha, yatlas=sfr_oii, xerratlas=sfr_ha*0.0, yerratlas=sfr_oii*0.0, $
      atlaspsize=0.5
    djs_oplot, !x.crange, !y.crange, line=0, thick=postthick
    
    sdss_lineplot, sfr_ha_sdss, sfr_ratio_sdss, sfr_ha_sdss*0.0, sfr_ratio_sdss*0.0, $
      plottype=3, postscript=postscript, ytitle=ytitle2, xrange=xrange2, yrange=yrange2, $
      legendtype=0, position=pos[*,1], charsize=1.6, xtitle=xtitle, $
      xatlas=sfr_ha, yatlas=sfr_ratio, xerratlas=sfr_ha*0.0, yerratlas=sfr_ratio*0.0, $
      atlaspsize=0.5, /noerase
    djs_oplot, !x.crange, [0,0], line=0, thick=postthick

    legend, textoidl(['SDSS Sample '+xstr_sdss,'Integrated Sample '+xstr]), $
      /right, /top, box=0, charsize=charsize_3, charthick=postthick, spacing=1.5

;   plotsym, 0, 1, /fill
;   djs_plot, sfr_ha, sfr_oii, ps=8, xsty=3, ysty=3, xrange=xrange, $
;     yrange=yrange, xthick=postthick, ythick=postthick, charsize=1.7, $
;     charthick=postthick, xtitle='', ytitle=ytitle, /xlog, /ylog, $
;     position=pos[*,0], xtickname=replicate(' ',10), color='dark green'
;   djs_oplot, 10^!x.crange, 10^!y.crange, line=0, thick=postthick
;
;   djs_plot, sfr_ha, sfr_ratio, ps=8, xsty=3, ysty=3, xrange=xrange2, $
;     yrange=yrange2, xthick=postthick, ythick=postthick, charsize=1.7, $
;     charthick=postthick, xtitle=xtitle, ytitle=ytitle2, /xlog, $
;     position=pos[*,1], /noerase, color='dark green'
;   djs_oplot, 10^!x.crange, [0,0], line=0, thick=postthick
;   legend, textoidl(xstr), /left, /top, box=0, charsize=charsize_3, charthick=postthick

    im_openclose, postscript=postscript, /close    

; ------------------------------------------------------------
; L([O II]) vs E(B-V) - integrated
; ------------------------------------------------------------

    psname = 'LOII_vs_ebv'
    im_openclose, pspath+psname, postscript=postscript, /encapsulated

    indx = where((atlasnodust.ebv_hahb_err gt 0.0) and (atlasnodust.oii_3727_lum[0] gt -900),nindx)

    x = atlasnodust[indx].oii_3727_lum[0] + alog10(lsun)
    xerr = atlasnodust[indx].oii_3727_lum[1]

    y = atlasnodust[indx].ebv_hahb
    yerr = atlasnodust[indx].ebv_hahb_err
    
    indxnfgs = where((nfgsnodust.ebv_hahb_err gt 0.0) and $
      (nfgsnodust.oii_3727_lum[0] gt -900),nindxnfgs)

    xnfgs = nfgsnodust[indxnfgs].oii_3727_lum[0] + alog10(lsun)
    xerrnfgs = nfgsnodust[indxnfgs].oii_3727_lum[1]

    ynfgs = nfgsnodust[indxnfgs].ebv_hahb
    yerrnfgs = nfgsnodust[indxnfgs].ebv_hahb_err

    xtitle = 'log L([O II]) [erg s^{-1}]'
    ytitle = 'E(B-V) [mag]'

    xrange = Loiikewleyrange
    yrange = ehbharange

    atlas1d_lineplot, x, y, xerr, yerr, plottype=3, postscript=postscript, $
      xtitle=xtitle, ytitle=ytitle, xrange=xrange, yrange=yrange, legendtype=0, $
      /left, /top, xnfgs=xnfgs, ynfgs=ynfgs, xerrnfgs=xerrnfgs, yerrnfgs=yerrnfgs

; overplot the Kewley relation    

    xaxis = findgen((xrange[1]-xrange[0])/0.01+1)*0.01+xrange[0]
    yaxis = 0.174*xaxis - 6.84
    chop = where(yaxis gt 0.0)
    djs_oplot, xaxis[chop], yaxis[chop], line=0, thick=postthick;, color='dark green'

    legend, ['KGJ04'], line=0, /right, /top, charsize=charsize, $
      charthick=postthick, box=0, thick=postthick, clear=postscript

    im_openclose, postscript=postscript, /close    

; ------------------------------------------------------------
; L([O II]) vs E(B-V) - SDSS
; ------------------------------------------------------------

    psname = 'sdss_LOII_vs_ebv'
    im_openclose, pspath+psname, postscript=postscript, /encapsulated

    indx = where((sdssnodust.ebv_hahb_err gt 0.0) and (sdssnodust.oii_3727_lum[0] gt -900.0),nindx)

    x = sdssnodust[indx].oii_3727_lum[0] + alog10(lsun)
    xerr = sdssnodust[indx].oii_3727_lum[1]

    y = sdssnodust[indx].ebv_hahb
    yerr = sdssnodust[indx].ebv_hahb_err
    
    xtitle = 'log L([O II]) [erg s^{-1}]'
    ytitle = 'E(B-V) [mag]'

    xrange = Loiikewleyrange
    yrange = ehbharange

    sdss_lineplot, x, y, xerr, yerr, plottype=3, postscript=postscript, $
      xtitle=xtitle, ytitle=ytitle, xrange=xrange, yrange=yrange, legendtype=0, $
      /left, /top

; overplot the Kewley relation    

    xaxis = findgen((xrange[1]-xrange[0])/0.01+1)*0.01+xrange[0]
    yaxis = 0.174*xaxis - 6.84
    chop = where(yaxis gt 0.0)
    djs_oplot, xaxis[chop], yaxis[chop], line=0, thick=postthick;, color='dark green'

;   legend, ['KGJ04'], line=0, /right, /top, charsize=1.5, $
;     charthick=postthick, box=0, color=djs_icolor(['dark green']), $
;     thick=postthick, clear=postscript

    im_openclose, postscript=postscript, /close    

; ------------------------------------------------------------
; L(B) vs E(Hb-Ha) - Integrated
; ------------------------------------------------------------

    psname = 'LB_vs_ehbha'
    im_openclose, pspath+psname, postscript=postscript, /encapsulated

    indx = where((atlasdust.b_lum_obs gt -900.0) and (atlasnodust.ehbha_err gt 0.0),nindx)

    x = atlasdust[indx].b_lum_obs
    xerr = atlasdust[indx].b_lum_obs_err
    xabs = atlasdust[indx].m_b_obs

    y = atlasnodust[indx].ehbha
    yerr = atlasnodust[indx].ehbha_err
    yebv = y/(k_lambda(4861,/odonnell)-k_lambda(6563,/odonnell))

    indxnfgs = where((nfgsdust.b_lum_obs gt -900.0) and (nfgsnodust.ehbha_err gt 0.0),nindxnfgs)

    xnfgs = nfgsdust[indxnfgs].b_lum_obs
    xerrnfgs = nfgsdust[indxnfgs].b_lum_obs_err

    ynfgs = nfgsnodust[indxnfgs].ehbha
    yerrnfgs = nfgsnodust[indxnfgs].ehbha_err
    
    xtitle = 'log L(B) [L(B)'+sunsymbol()+']'
    ytitle = 'E(H\beta-H\alpha) [mag]'

    xrange = LBrange
    yrange = ehbharange

    atlas1d_lineplot, x, y, xerr, yerr, plottype=3, postscript=postscript, $
      xtitle=xtitle, ytitle=ytitle, xrange=xrange, yrange=yrange, legendtype=0, $
      xstyle=11, ymargin=[4,3], ystyle=11, xmargin=[8,6], /left, /top, $
      xnfgs=xnfgs, ynfgs=ynfgs, xerrnfgs=xerrnfgs, yerrnfgs=yerrnfgs
    axis, /xaxis, xrange=interpol(xabs,x,!x.crange), xthick=postthick, xstyle=1, $
      charsize=singlecharsize, charthick=postthick
    xyouts, mean(!x.crange), !y.crange[1]+0.08*(!y.crange[1]-!y.crange[0]), textoidl('M_{B} [mag]'), $
      align=0.5, charsize=singlecharsize, charthick=postthick
    axis, /yaxis, yrange=interpol(yebv,y,!y.crange), ythick=postthick, $
      charthick=postthick, charsize=2.0, ytitle='E(B-V) [mag]', ysty=1
    
    im_openclose, postscript=postscript, /close    
    
; ------------------------------------------------------------
; L(B) vs E(Hb-Ha) - SDSS
; ------------------------------------------------------------

    psname = 'sdss_LB_vs_ehbha'
    im_openclose, pspath+psname, postscript=postscript, /encapsulated

    indx = where((sdssdust.b_lum gt -900.0) and (sdssnodust.ehbha_err gt 0.0),nindx)

    x = sdssdust[indx].b_lum
    xerr = sdssdust[indx].b_lum_err
    xabs = sdssdust[indx].m_b

    y = sdssnodust[indx].ehbha
    yerr = sdssnodust[indx].ehbha_err
    yebv = y/(k_lambda(4861,/odonnell)-k_lambda(6563,/odonnell))

    xtitle = 'log L(B) [L(B)'+sunsymbol()+']'
    ytitle = 'E(H\beta-H\alpha) [mag]'

    xrange = LBrange
    yrange = ehbharange

    sdss_lineplot, x, y, xerr, yerr, plottype=3, postscript=postscript, $
      xtitle=xtitle, ytitle=ytitle, xrange=xrange, yrange=yrange, legendtype=0, $
      xstyle=11, ymargin=[4,3], ystyle=11, xmargin=[8,6], /left, /top
    axis, /xaxis, xrange=interpol(xabs,x,!x.crange), xthick=postthick, xstyle=1, $
      charsize=singlecharsize, charthick=postthick
    xyouts, mean(!x.crange), !y.crange[1]+0.08*(!y.crange[1]-!y.crange[0]), textoidl('M_{B} [mag]'), $
      align=0.5, charsize=singlecharsize, charthick=postthick
    axis, /yaxis, yrange=interpol(yebv,y,!y.crange), ythick=postthick, $
      charthick=postthick, charsize=2.0, ytitle='E(B-V) [mag]', ysty=1
    
    im_openclose, postscript=postscript, /close    
    
; ------------------------------------------------------------
; EW([O II]) vs E(Hb-Ha)
; ------------------------------------------------------------

    psname = 'ewoii_vs_ehbha'
    im_openclose, pspath+psname, postscript=postscript, /encapsulated

    lineratio, atlasnodust, 'OII_3727_EW', '', '', '', x, xerr, $
      dum1, dum2, index=indx, nindex=nindx, snrcut=snrcut, $
      xsyserr=syserr, ysyserr=syserr

    y = atlasnodust[indx].ehbha
    yerr = atlasnodust[indx].ehbha_err
    yebv = y/(k_lambda(4861,/odonnell)-k_lambda(6563,/odonnell))

    lineratio, nfgsnodust, 'OII_3727_EW', '', '', '', xnfgs, xerrnfgs, $
      dum1, dum2, index=indxnfgs, nindex=nindxnfgs, snrcut=snrcut, $
      xsyserr=nfgssyserr, ysyserr=nfgssyserr
    ynfgs = nfgsnodust[indxnfgs].ehbha
    yerrnfgs = nfgsnodust[indxnfgs].ehbha_err
    
    xtitle = 'log EW([O II])  ['+angstrom()+']'
    ytitle = 'E(H\beta-H\alpha) [mag]'

    xrange = ewoiirange
    yrange = ehbharange

    atlas1d_lineplot, x, y, xerr, yerr, plottype=3, postscript=postscript, $
      xtitle=xtitle, ytitle=ytitle, xrange=xrange, yrange=yrange, legendtype=0, $
      ystyle=11, xmargin=[8,6], ymargin=[4,3], /right, /top, xnfgs=xnfgs, ynfgs=ynfgs, $
      xerrnfgs=xerrnfgs, yerrnfgs=yerrnfgs
    axis, /yaxis, yrange=interpol(yebv,y,!y.crange), ythick=postthick, $
      charthick=postthick, charsize=2.0, ytitle = 'E(B-V) [mag]', ysty=1

    im_openclose, postscript=postscript, /close        

; ------------------------------------------------------------
; EW([O II]) vs E(Hb-Ha)
; ------------------------------------------------------------

    psname = 'sdss_ewoii_vs_ehbha'
    im_openclose, pspath+psname, postscript=postscript, /encapsulated

    lineratio, sdssnodust, 'OII_3727_EW', '', '', '', x, xerr, $
      dum1, dum2, index=indx, nindex=nindx, snrcut=snrcut

    y = sdssnodust[indx].ehbha
    yerr = sdssnodust[indx].ehbha_err
    yebv = y/(k_lambda(4861,/odonnell)-k_lambda(6563,/odonnell))

    xtitle = 'log EW([O II])  ['+angstrom()+']'
    ytitle = 'E(H\beta-H\alpha) [mag]'

    xrange = ewoiirange
    yrange = ehbharange

    sdss_lineplot, x, y, xerr, yerr, plottype=3, postscript=postscript, $
      xtitle=xtitle, ytitle=ytitle, xrange=xrange, yrange=yrange, legendtype=0, $
      ystyle=11, xmargin=[8,6], ymargin=[4,3], /right, /top
    axis, /yaxis, yrange=interpol(yebv,y,!y.crange), ythick=postthick, $
      charthick=postthick, charsize=2.0, ytitle = 'E(B-V) [mag]', ysty=1

    im_openclose, postscript=postscript, /close        

; ------------------------------------------------------------
; Mass vs SFR - Fun plot for Janice - [Integrated + SDSS]
; ------------------------------------------------------------

    psname = 'mass_vs_sfr'
    im_openclose, cwd()+psname, postscript=postscript, xsize=8.5, ysize=5.0

    pagemaker, nx=2, ny=1, height=3.5, width=3.5, xmar=[1.2,0.3], $
      ymar=[0.5,1.0], xpage=8.5, ypage=5.0, pos=pos, xspace=0, /normal

; Integrated    
    
    indx = where((atlasnodust.sfr_h_alpha_err gt -900) and (atlasnodust.mass_bv_b gt -900),nindx)

    x = atlasnodust[indx].mass_bv_b
    xerr = atlasnodust[indx].mass_bv_b_err
    y = atlasnodust[indx].sfr_h_alpha
    yerr = atlasnodust[indx].sfr_h_alpha_err

    indxnfgs = where((nfgsnodust.sfr_h_alpha_err gt -900) and (nfgsnodust.mass_bv_b gt -900.0),nindxnfgs)

    xnfgs = nfgsnodust[indxnfgs].mass_bv_b
    xerrnfgs = nfgsnodust[indxnfgs].mass_bv_b_err
    ynfgs = nfgsnodust[indxnfgs].sfr_h_alpha
    yerrnfgs = nfgsnodust[indxnfgs].sfr_h_alpha_err

    xtitle = 'log (M / M'+sunsymbol()+')'
    ytitle = '\psi (H\alpha) [M'+sunsymbol()+' yr^{-1}]'

    xrange = [5.8,12.5]
    yrange = [-3.2,3.0]
    
    atlas1d_lineplot, x, y, xerr, yerr, plottype=3, postscript=postscript, $
      ytitle=ytitle, xtitle=xtitle, xrange=xrange, yrange=yrange, legendtype=0, $
      /right, /top, xnfgs=xnfgs, ynfgs=ynfgs, xerrnfgs=xerrnfgs, $
      yerrnfgs=yerrnfgs, position=pos[*,0], charsize=charsize_2

; SDSS    
    
;   indx = where((sdssnodust.sfr_h_alpha_err gt -900) and (sdssnodust.kauffmann_mass gt 0.0),nindx)
    indx = where((sdssnodust.sfr_h_alpha_err gt -900) and (sdssnodust.mass_bv_b gt -900),nindx)

;   x = sdssnodust[indx].kauffmann_mass
;   xerr = x*0.0
    x = sdssnodust[indx].mass_bv_b
    xerr = sdssnodust[indx].mass_bv_b_err
    y = sdssnodust[indx].sfr_h_alpha
    yerr = sdssnodust[indx].sfr_h_alpha_err

    sdss_lineplot, x, y, xerr, yerr, plottype=3, postscript=postscript, $
      ytitle='', xtitle=xtitle, xrange=xrange, yrange=yrange, legendtype=0, $
      /right, /top, position=pos[*,1], charsize=charsize_2, /noerase, $
      ytickname=replicate(' ',10)

    im_openclose, postscript=postscript, /close
    
; ------------------------------------------------------------
; 3-panel - L(B) vs L(U)/L(Ha) - SDSS
; ------------------------------------------------------------

    psname = 'sdss_LB_vs_LU_LHa_multipanel'
    im_openclose, pspath+psname, postscript=postscript, /encapsulated

    indx = where((sdssdust.b_lum gt -900.0) and (sdssdust.fiber_u gt -900) and $
      (sdssnodust.ebv_hahb_err gt 0),nindx)

    x = sdssdust[indx].b_lum
    xerr = sdssdust[indx].b_lum_err
    xabs = sdssdust[indx].m_b

    xtitle = 'log L(B) [L(B)'+sunsymbol()+']'
    ytitle = 'log [L(U)/L(H\alpha)] '

    xrange = LBrange
    yrange = Uhacorrange
    
    pagemaker, nx=1, ny=3, yspace=0, xspace=0, xmargin=[2.0,1.0], $
      ymargin=[1.1,1.1], width=5.5, position=pos, /normal

; ##########################################################
; Panel 1: U observed, Ha observed
; ##########################################################

    y1 = Uconstant*10^(-0.4*sdssdust[indx].fiber_U)
    y1err = 0.4*sdssdust[indx].fiber_U_err*alog(10.0)*y1
    
    y2 = sdssdust[indx].h_alpha[0]
    y2err = sdssdust[indx].h_alpha[1]

    y = alog10(y1/y2)
    yerr = im_compute_error(y1,y1err,y2,y2err,/log)

    stats = im_stats(y,/verbose)
    xstr = strtrim(string(stats.median,format='(F12.2)'),2)+'\pm'+strtrim(string(stats.sig68mean,format='(F12.2)'),2)
    rcor = r_correlate(x,y,zd=zd,probd=probd)
    splog, 'Spearman rank test: ', rcor, probd

    ytitle = 'log [L(U)/L(H\alpha)]_{obs} '

    sdss_lineplot, x, y, xerr, yerr, plottype=3, postscript=postscript, $
      ytitle=ytitle, xrange=xrange, yrange=yrange, legendtype=0, xminor=3, yminor=3, $
      /right, /top, position=pos[*,0], charsize=charsize_5, xtickname=replicate(' ',10)

    legend, '(a)', /left, /top, box=0, charsize=charsize_5, charthick=postthick
    legend, textoidl(xstr), /right, /top, box=0, charsize=charsize_3, charthick=postthick

; ##########################################################
; Panel 2: U observed, Ha corrected
; ##########################################################

    y1 = Uconstant*10^(-0.4*sdssdust[indx].fiber_U)
    y1err = 0.4*sdssdust[indx].fiber_U_err*alog(10.0)*y1
    
    y2 = sdssnodust[indx].h_alpha[0]
    y2err = sdssnodust[indx].h_alpha[1]

    y = alog10(y1/y2)
    yerr = im_compute_error(y1,y1err,y2,y2err,/log)

    stats = im_stats(y,/verbose)
    xstr = strtrim(string(stats.median,format='(F12.2)'),2)+'\pm'+strtrim(string(stats.sig68mean,format='(F12.2)'),2)
    rcor = r_correlate(x,y,zd=zd,probd=probd)
    splog, 'Spearman rank test: ', rcor, probd

    ytitle = 'log [L(U)_{obs}/L(H\alpha)_{cor}] '

    sdss_lineplot, x, y, xerr, yerr, plottype=3, postscript=postscript, $
      xrange=xrange, yrange=yrange, legendtype=0, /noerase, xtickname=replicate(' ',10), $
      /right, /top, position=pos[*,1], charsize=charsize_5, ytitle=ytitle
    
    legend, '(b)', /left, /top, box=0, charsize=charsize_5, charthick=postthick
    legend, textoidl(xstr), /right, /top, box=0, charsize=charsize_3, charthick=postthick

; ##########################################################
; Panel 3: U corrected, Ha corrected
; ##########################################################

    dust_fraction = 0.3
    
    AU = sdssnodust[indx].ebv_hahb*k_lambda(Uinfo.weff,/charlot)*dust_fraction
    AU_err = sdssnodust[indx].ebv_hahb_err*k_lambda(Uinfo.weff,/charlot)*dust_fraction

    y1 = Uconstant*10^(-0.4*sdssdust[indx].fiber_U)*10^(0.4*AU)
    y1err = sqrt((0.4*sdssdust[indx].fiber_U_err*alog(10.0)*y1*10^(0.4*AU))^2 + $
      (y1*alog(10.0)*AU_err)^2)

    y2 = sdssnodust[indx].h_alpha[0]
    y2err = sdssnodust[indx].h_alpha[1]

    y = alog10(y1/y2)
    yerr = im_compute_error(y1,y1err,y2,y2err,/log)

    stats = im_stats(y,/verbose)
    xstr = strtrim(string(stats.median,format='(F12.2)'),2)+'\pm'+strtrim(string(stats.sig68mean,format='(F12.2)'),2)
    rcor = r_correlate(x,y,zd=zd,probd=probd)
    splog, 'Spearman rank test: ', rcor, probd

    ytitle = 'log [L(U)/L(H\alpha)]_{cor}'

    sdss_lineplot, x, y, xerr, yerr, plottype=3, postscript=postscript, $
      xtitle=xtitle, ytitle=ytitle, xrange=xrange, yrange=yrange, legendtype=0, /noerase, $
      xstyle=3, /right, /top, position=pos[*,2], charsize=charsize_5, xminor=3, yminor=3
    
    legend, '(c)', /left, /top, box=0, charsize=charsize_5, charthick=postthick
    legend, textoidl(xstr), /right, /top, box=0, charsize=charsize_3, charthick=postthick

    im_openclose, postscript=postscript, /close

; ------------------------------------------------------------
; 3-panel - L(B) vs L(U)/L(Ha) - Integrated
; ------------------------------------------------------------

    psname = 'LB_vs_LU_LHa_multipanel'
    im_openclose, pspath+psname, postscript=postscript, /encapsulated

    pagemaker, nx=1, ny=3, yspace=0, xspace=0, xmargin=[2.0,1.0], $
      ymargin=[1.1,1.1], width=5.5, position=pos, /normal

    indx = where((atlasdust.b_lum_obs gt -900.0) and (atlasdust.u_obs gt -900) and $
      (atlasnodust.ebv_hahb_err gt 0),nindx)

    x = atlasdust[indx].b_lum_obs
    xerr = atlasdust[indx].b_lum_obs_err
    xabs = atlasnodust[indx].m_b_obs

    indxnfgs = where((nfgsdust.b_lum_obs gt -900.0) and (nfgsdust.u_obs gt -900) and $
      (nfgsnodust.ebv_hahb_err gt 0),nindxnfgs)

    xnfgs = nfgsdust[indxnfgs].b_lum_obs
    xerrnfgs = nfgsdust[indxnfgs].b_lum_obs_err

    xtitle = 'log L(B) [L(B)'+sunsymbol()+']'
    ytitle = 'log [L(U)/L(H\alpha)] '

    xrange = LBrange
    yrange = Uhacorrange
    
; ##########################################################
; Panel 1: U observed, Ha observed
; ##########################################################

;   y1 = atlasdust[indx].u_lum_obs
;   y1err = atlasdust[indx].u_lum_obs_err
;   y2 = atlasnodust[indx].h_alpha_lum[0]
;   y2err = atlasnodust[indx].h_alpha_lum[1]
;   y = y1 - y2
;   yerr = sqrt(y1err^2 + y2err^2)
;   y1nfgs = nfgsdust[indxnfgs].u_lum_obs
;   y1errnfgs = nfgsdust[indxnfgs].u_lum_obs_err
;   y2nfgs = nfgsdust[indxnfgs].h_alpha_lum[0]
;   y2errnfgs = nfgsdust[indxnfgs].h_alpha_lum[1]
;   ynfgs = y1nfgs - y2nfgs
;   yerrnfgs = sqrt(y1errnfgs^2 + y2errnfgs^2)

    y1 = Uconstant*10^(-0.4*atlasdust[indx].u_obs)
    y1err = 0.4*atlasdust[indx].u_obs_err*alog(10.0)*y1
    
    y2 = atlasdust[indx].h_alpha[0]
    y2err = atlasdust[indx].h_alpha[1]

    y = alog10(y1/y2)
    yerr = im_compute_error(y1,y1err,y2,y2err,/log)

    y1nfgs = Uconstant*10^(-0.4*nfgsdust[indxnfgs].u_obs)
    y1errnfgs = 0.4*nfgsdust[indxnfgs].u_obs_err*alog(10.0)*y1
    
    y2nfgs = nfgsdust[indxnfgs].h_alpha[0]
    y2errnfgs = nfgsdust[indxnfgs].h_alpha[1]

    ynfgs = alog10(y1nfgs/y2nfgs)
    yerrnfgs = im_compute_error(y1nfgs,y1errnfgs,y2nfgs,y2errnfgs,/log)
    
    stats = im_stats([y,ynfgs],/verbose)
    xstr = strtrim(string(stats.median,format='(F12.2)'),2)+'\pm'+strtrim(string(stats.sig68mean,format='(F12.2)'),2)
    rcor = r_correlate([x,xnfgs],[y,ynfgs],zd=zd,probd=probd)
    splog, 'Spearman rank test: ', rcor, probd

    ytitle = 'log [L(U)/L(H\alpha)]_{obs} '
;   ytitle = 'log [L(U)_{obs}/L(H\alpha)_{obs}] '

    atlas1d_lineplot, x, y, xerr, yerr, plottype=3, postscript=postscript, $
      ytitle=ytitle, xrange=xrange, yrange=yrange, legendtype=0, xminor=3, yminor=3, $
      /right, /top, xnfgs=xnfgs, ynfgs=ynfgs, xerrnfgs=xerrnfgs, $
      yerrnfgs=yerrnfgs, position=pos[*,0], charsize=charsize_5, xtickname=replicate(' ',10)

    legend, '(a)', /left, /top, box=0, charsize=charsize_5, charthick=postthick
    legend, textoidl(xstr), /right, /top, box=0, charsize=charsize_3, charthick=postthick

;   label = ['A(U) = 0','A(H\alpha) = 0']
;   legend, textoidl(label), /left, /bottom, box=0, charsize=charsize_3, $
;     charthick=postthick

; ##########################################################
; Panel 2: U observed, Ha corrected
; ##########################################################

    y1 = Uconstant*10^(-0.4*atlasdust[indx].u_obs)
    y1err = 0.4*atlasdust[indx].u_obs_err*alog(10.0)*y1
    
    y2 = atlasnodust[indx].h_alpha[0]
    y2err = atlasnodust[indx].h_alpha[1]

    y = alog10(y1/y2)
    yerr = im_compute_error(y1,y1err,y2,y2err,/log)

    y1nfgs = Uconstant*10^(-0.4*nfgsdust[indxnfgs].u_obs)
    y1errnfgs = 0.4*nfgsdust[indxnfgs].u_obs_err*alog(10.0)*y1
    
    y2nfgs = nfgsnodust[indxnfgs].h_alpha[0]
    y2errnfgs = nfgsnodust[indxnfgs].h_alpha[1]

    ynfgs = alog10(y1nfgs/y2nfgs)
    yerrnfgs = im_compute_error(y1nfgs,y1errnfgs,y2nfgs,y2errnfgs,/log)
    
    stats = im_stats([y,ynfgs],/verbose)
    xstr = strtrim(string(stats.median,format='(F12.2)'),2)+'\pm'+strtrim(string(stats.sig68mean,format='(F12.2)'),2)
    rcor = r_correlate([x,xnfgs],[y,ynfgs],zd=zd,probd=probd)
    splog, 'Spearman rank test: ', rcor, probd

    ytitle = 'log [L(U)_{obs}/L(H\alpha)_{cor}] '

    atlas1d_lineplot, x, y, xerr, yerr, plottype=3, postscript=postscript, $
      xrange=xrange, yrange=yrange, legendtype=0, /noerase, xtickname=replicate(' ',10), $
      /right, /top, xnfgs=xnfgs, ynfgs=ynfgs, xerrnfgs=xerrnfgs, xminor=3, yminor=3, $
      yerrnfgs=yerrnfgs, position=pos[*,1], charsize=charsize_5, ytitle=ytitle
    
    legend, '(b)', /left, /top, box=0, charsize=charsize_5, charthick=postthick
    legend, textoidl(xstr), /right, /top, box=0, charsize=charsize_3, charthick=postthick

;   label = ['A(U) = 0','Individual A(H\alpha)']
;   legend, textoidl(label), /left, /bottom, box=0, charsize=charsize_3, $
;     charthick=postthick

; ##########################################################
; Panel 3: U corrected, Ha corrected
; ##########################################################

    dust_fraction = 0.3
    
    AU = atlasnodust[indx].ebv_hahb*k_lambda(Uinfo.weff,/charlot)*dust_fraction
    AU_err = atlasnodust[indx].ebv_hahb_err*k_lambda(Uinfo.weff,/charlot)*dust_fraction

    y1 = Uconstant*10^(-0.4*atlasdust[indx].u_obs)*10^(0.4*AU)
    y1err = sqrt((0.4*atlasdust[indx].u_obs_err*alog(10.0)*y1*10^(0.4*AU))^2 + $
      (y1*alog(10.0)*AU_err)^2)

    y2 = atlasnodust[indx].h_alpha[0]
    y2err = atlasnodust[indx].h_alpha[1]

    y = alog10(y1/y2)
    yerr = im_compute_error(y1,y1err,y2,y2err,/log)

; NFGS    
    
    AUnfgs = nfgsnodust[indxnfgs].ebv_hahb*k_lambda(Uinfo.weff,/charlot)*dust_fraction
    AUnfgs_err = nfgsnodust[indxnfgs].ebv_hahb_err*k_lambda(Uinfo.weff,/charlot)*dust_fraction

    y1nfgs = Uconstant*10^(-0.4*nfgsdust[indxnfgs].u_obs)*10^(0.4*AUnfgs)
    y1errnfgs = sqrt((0.4*nfgsdust[indxnfgs].u_obs_err*alog(10.0)*y1nfgs*10^(0.4*AUnfgs))^2 + $
      (y1nfgs*alog(10.0)*AUnfgs_err)^2)

    y2nfgs = nfgsnodust[indxnfgs].h_alpha[0]
    y2errnfgs = nfgsnodust[indxnfgs].h_alpha[1]

    ynfgs = alog10(y1nfgs/y2nfgs)
    yerrnfgs = im_compute_error(y1nfgs,y1errnfgs,y2nfgs,y2errnfgs,/log)
    
    stats = im_stats([y,ynfgs],/verbose)
    xstr = strtrim(string(stats.median,format='(F12.2)'),2)+'\pm'+strtrim(string(stats.sig68mean,format='(F12.2)'),2)
    rcor = r_correlate([x,xnfgs],[y,ynfgs],zd=zd,probd=probd)
    splog, 'Spearman rank test: ', rcor, probd

    ytitle = 'log [L(U)/L(H\alpha)]_{cor}'

    atlas1d_lineplot, x, y, xerr, yerr, plottype=3, postscript=postscript, $
      xtitle=xtitle, ytitle=ytitle, xrange=xrange, yrange=yrange, legendtype=0, $
      /noerase, xstyle=3, /right, /top, xnfgs=xnfgs, ynfgs=ynfgs, xerrnfgs=xerrnfgs, $
      yerrnfgs=yerrnfgs, position=pos[*,2], charsize=charsize_5, xminor=3, yminor=3
    
    legend, '(c)', /left, /top, box=0, charsize=charsize_5, charthick=postthick
    legend, textoidl(xstr), /right, /top, box=0, charsize=charsize_3, charthick=postthick

;   label = ['Individual A(U)','Individual A(H\alpha)']
;   legend, textoidl(label), /left, /bottom, box=0, charsize=charsize_3, $
;     charthick=postthick

    im_openclose, postscript=postscript, /close

; ------------------------------------------------------------
; 12+log(O/H) [OIIINII and NIIHA] vs [O III]_cor/Ha_cor [Integrated + SDSS]
; ------------------------------------------------------------

    psname = '12oh_OIIINII_NIIHA_vs_oiiiha_multipanel'
    im_openclose, pspath+psname, postscript=postscript, /encapsulated, xsize=8.5, ysize=5.7

    pagemaker, nx=2, ny=1, height=3.5, width=3.5, xmargin=[1.2,0.3], $
      ymargin=[0.5,1.0], xspace=0, yspace=0, xpage=8.5, ypage=5.0, $
      position=pos, /normal

    xtitle = '12 + log (O/H)'
    ytitle = 'log ([O III] \lambda5007/H\alpha)_{cor}'

    xrange = hiiohrange
    yrange = oiiiharange

; Integrated    
    
    indx1 = where(atlasdust.zstrong_12oh_oiiinii_pettini gt -900,nindx1)
    lineratio, atlasnodust[indx1], '', '', 'OIII_5007', 'H_ALPHA', $
      dum1, dum2, y, yerr, index=indx2, nindex=nindx2, snrcut=snrcut, $
      xsyserr=syserr, ysyserr=syserr
    x = atlasdust[indx1[indx2]].zstrong_12oh_oiiinii_pettini
    xerr = atlasdust[indx1[indx2]].zstrong_12oh_oiiinii_pettini_err

    indx3 = where((atlasdust.zstrong_12oh_niiha_denicolo gt -900) and (atlasdust.zstrong_12oh_oiiinii_pettini lt -900),nindx3)
    if (nindx3 ne 0L) then begin
       lineratio, atlasnodust[indx3], '', '', 'OIII_5007', 'H_ALPHA', $
         dum1, dum2, y2, y2err, index=indx4, nindex=nindx4, snrcut=snrcut, $
      xsyserr=syserr, ysyserr=syserr
       y = [y,y2]
       yerr = [yerr,y2err]
       x = [x,atlasdust[indx3[indx4]].zstrong_12oh_niiha_denicolo]
       xerr = [xerr,atlasdust[indx3[indx4]].zstrong_12oh_niiha_denicolo_err]
    endif

;   niceprint, x, xerr, y, yerr

    cutnfgs = where(nfgsdust.zstrong_12oh_oiiinii_pettini gt -900,ncutnfgs)
    lineratio, nfgsnodust[cutnfgs], '', '', 'OIII_5007', 'H_ALPHA', $
      dum1, dum2, ynfgs, yerrnfgs, index=indx2nfgs, nindex=nindx2nfgs, snrcut=snrcut, $
      xsyserr=nfgssyserr, ysyserr=nfgssyserr
    xnfgs = nfgsdust[cutnfgs[indx2nfgs]].zstrong_12oh_oiiinii_pettini
    xerrnfgs = nfgsdust[cutnfgs[indx2nfgs]].zstrong_12oh_oiiinii_pettini_err

    indx3nfgs = where((nfgsdust.zstrong_12oh_niiha_denicolo gt -900) and (nfgsdust.zstrong_12oh_oiiinii_pettini lt -900),nindx3nfgs)
    if (nindx3nfgs ne 0L) then begin
       lineratio, nfgsnodust[indx3], '', '', 'OIII_5007', 'H_ALPHA', $
         dum1, dum2, y2nfgs, y2errnfgs, index=indx4nfgs, nindex=nindx4nfgs, snrcut=snrcut, $
      xsyserr=nfgssyserr, ysyserr=nfgssyserr
       ynfgs = [ynfgs,y2nfgs]
       yerrnfgs = [yerrnfgs,y2errnfgs]
       xnfgs = [xnfgs,nfgsdust[indx3nfgs[indx4nfgs]].zstrong_12oh_niiha_denicolo]
       xerrnfgs = [xerrnfgs,nfgsdust[indx3nfgs[indx4nfgs]].zstrong_12oh_niiha_denicolo_err]
    endif

    xbig = [x,xnfgs]
    ybig = [y,ynfgs]

; HII regions    
    
    indx1 = where((hii.zstrong_12oh_oiiinii_pettini gt -900) and (hii.oiii_h_alpha gt -900),nindx1)
    yregion = hii[indx1].oiii_h_alpha
    yerrregion = hii[indx1].oiii_h_alpha_err
    xregion = hii[indx1].zstrong_12oh_oiiinii_pettini
    xerrregion = hii[indx1].zstrong_12oh_oiiinii_pettini_err

    indx2 = where((hii.zstrong_12oh_niiha_denicolo gt -900) and (hii.zstrong_12oh_oiiinii_pettini lt -900) and $
      (hii.oiii_h_alpha gt -900),nindx2)
    if (nindx2 ne 0L) then begin
       yregion = [yregion,hii[indx2].oiii_h_alpha]
       yerrregion = [yerrregion,hii[indx2].oiii_h_alpha_err]
       xregion = [xregion,hii[indx2].zstrong_12oh_niiha_denicolo]
       xerrregion = [xerrregion,hii[indx2].zstrong_12oh_niiha_denicolo_err]
    endif
       
    atlas1d_lineplot, x, y, xerr, yerr, plottype=3, postscript=postscript, $n
      xtitle=xtitle, ytitle=ytitle, xrange=xrange, yrange=yrange, legendtype=0, $
      /left, /top, position=pos[*,0], charsize=charsize_5, yfrac=1.5, /errorleft, $
      xregion=xregion, yregion=yregion, xerrregion=xerrregion, yerrregion=yerrregion, $
      xnfgs=xnfgs, ynfgs=ynfgs, xerrnfgs=xerrnfgs, yerrnfgs=yerrnfgs

    if keyword_set(postscript) then thisthick = 6.0 else thickthick = 2.0

    plot_kewley_grids, plotnumber=17, model=3, labeltype=4, /log12oh, /overplot, $
      /noZgrid, thick=3.0, Ugridcolor='dark green', Ulinestyle=0, $
      postscript=postscript, Ulabelcolor='black', Ucharsize=0.8, charthick=thisthick
    xyouts, 9.1, -0.18, 'log U', align=0.5, /data, charsize=0.9, charthick=thisthick

; overlay the three metallicity regions

    r12 = 8.15
    r23 = 8.7

    oplot, r12*[1,1], [-1.25,0.75], line=1, thick=postthick
    oplot, r23*[1,1], [-1.25,0.75], line=1, thick=postthick
 
    xyouts, 7.7, 0.65, 'R1', align=0.5, /data, charsize=charsize_5, charthick=postthick
    xyouts, 8.4, 0.65, 'R2', align=0.5, /data, charsize=charsize_5, charthick=postthick
    xyouts, 9.0, 0.65, 'R3', align=0.5, /data, charsize=charsize_5, charthick=postthick

    region1 = where(xbig lt r12,nregion1)
    splog, 'Region 1, 2, 3: '
    stats = im_stats(ybig[region1],sigrej=3.0,/verbose)

    regioniiha = where((xbig gt r12) and (xbig lt r23),nregioniiha)
    stats = im_stats(ybig[regioniiha],sigrej=3.0,/no_head,/verbose)

    region3 = where(xbig gt r23,nregion3)
    stats = im_stats(ybig[region3],sigrej=3.0,/no_head,/verbose)

; SDSS    

    indx1 = where(sdssdust.zstrong_12oh_oiiinii_pettini gt -900,nindx1)
    lineratio, sdssnodust[indx1], '', '', 'OIII_5007', 'H_ALPHA', $
      dum1, dum2, y, yerr, index=indx2, nindex=nindx2, snrcut=snrcut
    x = sdssdust[indx1[indx2]].zstrong_12oh_oiiinii_pettini
    xerr = sdssdust[indx1[indx2]].zstrong_12oh_oiiinii_pettini_err

    indx3 = where((sdssdust.zstrong_12oh_niiha_denicolo gt -900) and (sdssdust.zstrong_12oh_oiiinii_pettini lt -900),nindx3)
    if (nindx3 ne 0L) then begin
       lineratio, sdssnodust[indx3], '', '', 'OIII_5007', 'H_ALPHA', $
         dum1, dum2, y2, y2err, index=indx4, nindex=nindx4, snrcut=snrcut
       y = [y,y2]
       yerr = [yerr,y2err]
       x = [x,sdssdust[indx3[indx4]].zstrong_12oh_niiha_denicolo]
       xerr = [xerr,sdssdust[indx3[indx4]].zstrong_12oh_niiha_denicolo_err]
    endif

    xbig = x
    ybig = y
    
    sdss_lineplot, x, y, xerr, yerr, plottype=3, postscript=postscript, $
      xtitle=xtitle, ytitle='', xrange=xrange, yrange=yrange, legendtype=0, $
      /left, /top, position=pos[*,1], yfrac=1.5, /errorleft, $
      ytickname=replicate(' ',10), /noerase, charsize=charsize_5;, $
;     xregion=xregion, yregion=yregion, xerrregion=xerrregion, yerrregion=yerrregion

    if keyword_set(postscript) then thisthick = 6.0 else thickthick = 2.0

    plot_kewley_grids, plotnumber=17, model=3, labeltype=4, /log12oh, /overplot, $
      /noZgrid, thick=3.0, Ugridcolor='dark green', Ulinestyle=0, $
      postscript=postscript, Ulabelcolor='black', Ucharsize=0.8, charthick=thisthick
    xyouts, 9.1, -0.18, 'log U', align=0.5, /data, charsize=0.9, charthick=thisthick

; overlay the three metallicity regions

    r12 = 8.14
    r23 = 8.68

    oplot, r12*[1,1], [-1.25,0.75], line=1, thick=postthick
    oplot, r23*[1,1], [-1.25,0.75], line=1, thick=postthick
 
    xyouts, 7.7, 0.65, 'R1', align=0.5, /data, charsize=charsize_5, charthick=postthick
    xyouts, 8.4, 0.65, 'R2', align=0.5, /data, charsize=charsize_5, charthick=postthick
    xyouts, 9.0, 0.65, 'R3', align=0.5, /data, charsize=charsize_5, charthick=postthick

    splog, 'Region 1, 2, 3: '
    region1 = where(xbig lt r12,nregion1)
    if (nregion1 ne 0L) then stats = im_stats(ybig[region1],sigrej=3.0,/verbose)

    regioniiha = where((xbig gt r12) and (xbig lt r23),nregioniiha)
    if (nregioniiha ne 0L) then stats = im_stats(ybig[regioniiha],sigrej=3.0,/no_head,/verbose)

    region3 = where(xbig gt r23,nregion3)
    if (nregion3 ne 0L) then stats = im_stats(ybig[region3],sigrej=3.0,/no_head,/verbose)

    im_openclose, postscript=postscript, /close    

; ------------------------------------------------------------
; OIIINII vs [O II]_cor/Ha_cor [Integrated + SDSS]
; ------------------------------------------------------------

    psname = 'OIIINII_vs_oiiha_multipanel'
    im_openclose, pspath+psname, postscript=postscript, /encapsulated, xsize=8.5, ysize=5.7

    pagemaker, nx=2, ny=1, height=3.5, width=3.5, xmargin=[1.2,0.3], $
      ymargin=[1.2,1.0], xspace=0, yspace=0, xpage=8.5, ypage=5.7, $
      position=pos, /normal

    xtitle = 'log {([O III]/H\beta)/([N II]/H\alpha)}_{obs}'
;   xtitle = 'log {([O III]/H\beta)/([N II]/H\alpha)}_{cor}'
;   xtitle = 'log {([O III] \lambda5007/H\beta)/([N II] \lambda6584/H\alpha)}'
    ytitle = 'log ([O II]/H\alpha)_{cor}'

    xrange = (hiiohrange-8.73)/(-0.32)
;   xrange = reverse(oiiiniirange)
    yrange = oiiharange

; Integrated    
    
    cut = where(atlasdust.zstrong_oiiinii gt -900)
    lineratio, atlasnodust[cut], '', '', 'OII_3727', 'H_ALPHA', $
      dum1, dum2, y, yerr, index=indx, nindex=nindx, snrcut=snrcut, $
      xsyserr=syserr, ysyserr=syserr
    x = atlasdust[cut[indx]].zstrong_oiiinii
    xerr = atlasdust[cut[indx]].zstrong_oiiinii_err

    OHgood = where(atlasdust[cut[indx]].zstrong_12oh_oiiinii_pettini gt -900)
    xgood = x[OHgood]
    xOHgood = atlasdust[cut[indx[OHgood]]].zstrong_12oh_oiiinii_pettini

    cutnfgs = where(nfgsdust.zstrong_oiiinii gt -900)
    lineratio, nfgsnodust[cut], '', '', 'OII_3727', 'H_ALPHA', $
      dum1, dum2, ynfgs, yerrnfgs, index=indxnfgs, snrcut=snrcut, $
      xsyserr=nfgssyserr, ysyserr=nfgssyserr
    xnfgs = nfgsdust[cutnfgs[indxnfgs]].zstrong_oiiinii
    xerrnfgs = nfgsdust[cutnfgs[indxnfgs]].zstrong_oiiinii_err

    xbig = [x,xnfgs]
    ybig = [y,ynfgs]

    good = where((hii.zstrong_oiiinii gt -900.0) and (hii.oii_h_alpha gt -900.0))
    xregion = hii[good].zstrong_oiiinii & xerrregion = hii[good].zstrong_oiiinii_err
    yregion = hii[good].oii_h_alpha & yerrregion = hii[good].oii_h_alpha_err
       
    atlas1d_lineplot, x, y, xerr, yerr, plottype=3, postscript=postscript, $n
      xtitle=xtitle, ytitle=ytitle, xrange=xrange, yrange=yrange, legendtype=0, $
      xstyle=11, ymargin=[4,3], /left, /top, position=pos[*,0], charsize=charsize_5, $
      /xreverse, xfrac=10, /errorleft, $
      xregion=xregion, yregion=yregion, xerrregion=xerrregion, yerrregion=yerrregion, $
      xnfgs=xnfgs, ynfgs=ynfgs, xerrnfgs=xerrnfgs, yerrnfgs=yerrnfgs
    
; overlay the three metallicity regions

    OH12 = 8.14 & OH23 = 8.68
    r12 = interpol(xgood,xOHgood,OH12) & r23 = interpol(xgood,xOHgood,OH23)
    OH12 = interpol(xOHgood,xgood,r12) & OH23 = interpol(xOHgood,xgood,r23)
    Z12 = 10^(OH12 - Zsun_new) & Z23 = 10^(OH23 - Zsun_new)    

    splog, 'OIIINII = '+string(r12,format='(F4.1)')+' --> 12+log(O/H) = '+$
      string(OH12,format='(F4.2)')+', Z/Z_sun = '+string(Z12,format='(F4.2)')
    splog, 'OIIINII = '+string(r23,format='(F4.1)')+' --> 12+log(O/H) = '+$
      string(OH23,format='(F4.2)')+', Z/Z_sun = '+string(Z23,format='(F4.2)')
    
    oplot, r12*[1,1], [-1.25,0.75], line=1, thick=postthick
    oplot, r23*[1,1], [-1.25,0.75], line=1, thick=postthick
 
    xyouts, -0.84, 0.65, 'R3', align=0.5, /data, charsize=charsize_5, charthick=postthick
    xyouts,  1.03, 0.65, 'R2', align=0.5, /data, charsize=charsize_5, charthick=postthick
    xyouts,  3.22, 0.65, 'R1', align=0.5, /data, charsize=charsize_5, charthick=postthick

    splog, 'Region 1, 2, 3: '
    region1 = where(xbig gt r12,nregion1)
    if (nregion1 ne 0L) then stats = im_stats(ybig[region1],/verbose)

    regioniiha = where((xbig lt r12) and (xbig gt r23),nregioniiha)
    if (nregioniiha ne 0L) then stats = im_stats(ybig[regioniiha],/verbose)

    region3 = where(xbig lt r23,nregion3)    
    if (nregion3 ne 0L) then stats = im_stats(ybig[region3],/verbose)

; now label the upper abscissa

    axis, /xaxis, xrange=interpol(xOHgood,xgood,!x.crange), xthick=postthick, xstyle=1, $
      charsize=charsize_5, charthick=postthick, /save
    xyouts, mean(!x.crange), !y.crange[1]+0.10*(!y.crange[1]-!y.crange[0]), textoidl('12 + log (O/H)'), $
      align=0.5, charsize=charsize_5, charthick=postthick
    
    if keyword_set(postscript) then thisthick = 6.0 else thickthick = 2.0
    
    plot_kewley_grids, plotnumber=17, model=3, labeltype=4, /log12oh, /overplot, $
      /noZgrid, thick=3.0, Ugridcolor='dark green', Ulinestyle=0, $
      postscript=postscript, Ulabelcolor='black', Ucharsize=0.8, charthick=thisthick
    xyouts, 9.1, -0.18, 'log U', align=0.5, /data, charsize=0.9, charthick=thisthick

; SDSS    

    cut = where(sdssdust.zstrong_oiiinii gt -900)
    lineratio, sdssnodust[cut], '', '', 'OII_3727', 'H_ALPHA', $
      dum1, dum2, y, yerr, index=indx, nindex=nindx, snrcut=snrcut
    x = sdssdust[cut[indx]].zstrong_oiiinii
    xerr = sdssdust[cut[indx]].zstrong_oiiinii_err

    OHgood = where(sdssdust[cut[indx]].zstrong_12oh_oiiinii_pettini gt -900)
    xgood = x[OHgood]
    xOHgood = sdssdust[cut[indx[OHgood]]].zstrong_12oh_oiiinii_pettini

    sdss_lineplot, x, y, xerr, yerr, plottype=3, postscript=postscript, $
      xtitle=xtitle, ytitle='', xrange=xrange, yrange=yrange, legendtype=0, $
      xstyle=11, ymargin=[4,3], /left, /top, /xreverse, xfrac=10, position=pos[*,1], $
      ytickname=replicate(' ',10), /noerase, charsize=charsize_5, /errorleft ;, $
;     xregion=xregion, yregion=yregion, xerrregion=xerrregion, yerrregion=yerrregion

; overlay the three metallicity regions

    OH12 = 8.14 & OH23 = 8.68
    r12 = interpol(xgood,xOHgood,OH12) & r23 = interpol(xgood,xOHgood,OH23)
    OH12 = interpol(xOHgood,xgood,r12) & OH23 = interpol(xOHgood,xgood,r23)
    Z12 = 10^(OH12 - Zsun_new) & Z23 = 10^(OH23 - Zsun_new)    

    splog, 'OIIINII = '+string(r12,format='(F4.1)')+' --> 12+log(O/H) = '+$
      string(OH12,format='(F4.2)')+', Z/Z_sun = '+string(Z12,format='(F4.2)')
    splog, 'OIIINII = '+string(r23,format='(F4.1)')+' --> 12+log(O/H) = '+$
      string(OH23,format='(F4.2)')+', Z/Z_sun = '+string(Z23,format='(F4.2)')
    
    oplot, r12*[1,1], [-1.25,0.75], line=1, thick=postthick
    oplot, r23*[1,1], [-1.25,0.75], line=1, thick=postthick
 
    xyouts, -0.84, 0.65, 'R3', align=0.5, /data, charsize=charsize_5, charthick=postthick
    xyouts,  1.03, 0.65, 'R2', align=0.5, /data, charsize=charsize_5, charthick=postthick
    xyouts,  3.22, 0.65, 'R1', align=0.5, /data, charsize=charsize_5, charthick=postthick

    splog, 'Region 1, 2, 3: '
    region1 = where(x gt r12,nregion1)
    if (nregion1 ne 0L) then stats = im_stats(y[region1],/verbose)

    regioniiha = where((x lt r12) and (x gt r23),nregioniiha)
    if (nregioniiha ne 0L) then stats = im_stats(y[regioniiha],/verbose)

    region3 = where(x lt r23,nregion3)    
    if (nregion3 ne 0L) then stats = im_stats(y[region3],/verbose)

; now label the upper abscissa

    axis, /xaxis, xrange=interpol(xOHgood,xgood,!x.crange), xthick=postthick, xstyle=1, $
      charsize=charsize_5, charthick=postthick, /save
    xyouts, mean(!x.crange), !y.crange[1]+0.10*(!y.crange[1]-!y.crange[0]), textoidl('12 + log (O/H)'), $
      align=0.5, charsize=charsize_5, charthick=postthick

    if keyword_set(postscript) then thisthick = 6.0 else thickthick = 2.0
    
    plot_kewley_grids, plotnumber=17, model=3, labeltype=4, /log12oh, /overplot, $
      /noZgrid, thick=3.0, Ugridcolor='dark green', Ulinestyle=0, $
      postscript=postscript, Ulabelcolor='black', Ucharsize=0.8, charthick=thisthick
    xyouts, 9.1, -0.18, 'log U', align=0.5, /data, charsize=0.9, charthick=thisthick
    
    im_openclose, postscript=postscript, /close    

; ------------------------------------------------------------
; 12+log(O/H) [OIIINII] vs [O II]_cor/Ha_cor [Integrated + SDSS]
; ------------------------------------------------------------

    psname = '12oh_OIIINII_vs_oiiha_multipanel'
    im_openclose, pspath+psname, postscript=postscript, /encapsulated, xsize=8.5, ysize=5.7

    pagemaker, nx=2, ny=1, height=3.5, width=3.5, xmargin=[1.2,0.3], $
      ymargin=[0.5,1.0], xspace=0, yspace=0, xpage=8.5, ypage=5.0, $
      position=pos, /normal

    xtitle = '12 + log (O/H)'
    ytitle = 'log ([O II]/H\alpha)_{cor}'

    xrange = hiiohrange
    yrange = oiiharange

; Integrated    
    
    indx1 = where(atlasdust.zstrong_12oh_oiiinii_pettini gt -900,nindx1)
    lineratio, atlasnodust[indx1], '', '', 'OII_3727', 'H_ALPHA', $
      dum1, dum2, y, yerr, index=indx2, nindex=nindx2, snrcut=snrcut, $
      xsyserr=syserr, ysyserr=syserr
    x = atlasdust[indx1[indx2]].zstrong_12oh_oiiinii_pettini
    xerr = atlasdust[indx1[indx2]].zstrong_12oh_oiiinii_pettini_err

    indx1nfgs = where(nfgsdust.zstrong_12oh_oiiinii_pettini gt -900,nindx1nfgs)
    lineratio, nfgsnodust[indx1nfgs], '', '', 'OII_3727', 'H_ALPHA', $
      dum1, dum2, ynfgs, yerrnfgs, index=indx2nfgs, nindex=nindx2nfgs, snrcut=snrcut, $
      xsyserr=nfgssyserr, ysyserr=nfgssyserr
    xnfgs = nfgsdust[indx1nfgs[indx2nfgs]].zstrong_12oh_oiiinii_pettini
    xerrnfgs = nfgsdust[indx1nfgs[indx2nfgs]].zstrong_12oh_oiiinii_pettini_err

    xbig = [x,xnfgs]
    ybig = [y,ynfgs]

; HII regions    
    
    indx1 = where((hii.zstrong_12oh_oiiinii_pettini gt -900) and (hii.oii_h_alpha gt -900),nindx1)
    yregion = hii[indx1].oii_h_alpha
    yerrregion = hii[indx1].oii_h_alpha_err
    xregion = hii[indx1].zstrong_12oh_oiiinii_pettini
    xerrregion = hii[indx1].zstrong_12oh_oiiinii_pettini_err

    atlas1d_lineplot, x, y, xerr, yerr, plottype=3, postscript=postscript, $n
      xtitle=xtitle, ytitle=ytitle, xrange=xrange, yrange=yrange, legendtype=0, $
      /left, /top, position=pos[*,0], charsize=charsize_5, yfrac=1.5, /errorleft, $
      xregion=xregion, yregion=yregion, xerrregion=xerrregion, yerrregion=yerrregion, $
      xnfgs=xnfgs, ynfgs=ynfgs, xerrnfgs=xerrnfgs, yerrnfgs=yerrnfgs

    if keyword_set(postscript) then thisthick = 6.0 else thickthick = 2.0

    plot_kewley_grids, plotnumber=17, model=3, labeltype=4, /log12oh, /overplot, $
      /noZgrid, thick=3.0, Ugridcolor='dark green', Ulinestyle=0, $
      postscript=postscript, Ulabelcolor='black', Ucharsize=0.8, charthick=thisthick
    xyouts, 9.1, -0.18, 'log U', align=0.5, /data, charsize=0.9, charthick=thisthick

; overlay the three metallicity regions

    r12 = 8.14
    r23 = 8.68

    oplot, r12*[1,1], [-1.25,0.75], line=1, thick=postthick
    oplot, r23*[1,1], [-1.25,0.75], line=1, thick=postthick
 
    xyouts, 7.7, 0.65, 'R1', align=0.5, /data, charsize=charsize_5, charthick=postthick
    xyouts, 8.4, 0.65, 'R2', align=0.5, /data, charsize=charsize_5, charthick=postthick
    xyouts, 9.0, 0.65, 'R3', align=0.5, /data, charsize=charsize_5, charthick=postthick

    splog, 'Region 1, 2, 3: '
    region1 = where(xbig lt r12,nregion1)
    if (nregion1 gt 2L) then stats = im_stats(ybig[region1],sigrej=3.0,/verbose)

    regioniiha = where((xbig gt r12) and (xbig lt r23),nregioniiha)
    if (nregioniiha gt 2L) then stats = im_stats(ybig[regioniiha],sigrej=3.0,/no_head,/verbose)

    region3 = where(xbig gt r23,nregion3)
    if (nregion3 gt 2L) then stats = im_stats(ybig[region3],sigrej=3.0,/no_head,/verbose)

; SDSS    

    indx1 = where(sdssdust.zstrong_12oh_oiiinii_pettini gt -900,nindx1)
    lineratio, sdssnodust[indx1], '', '', 'OII_3727', 'H_ALPHA', $
      dum1, dum2, y, yerr, index=indx2, nindex=nindx2, snrcut=snrcut
    x = sdssdust[indx1[indx2]].zstrong_12oh_oiiinii_pettini
    xerr = sdssdust[indx1[indx2]].zstrong_12oh_oiiinii_pettini_err

    xbig = x
    ybig = y
    
    sdss_lineplot, x, y, xerr, yerr, plottype=3, postscript=postscript, $
      xtitle=xtitle, ytitle='', xrange=xrange, yrange=yrange, legendtype=0, $
      /left, /top, position=pos[*,1], yfrac=1.5, /errorleft, $
      ytickname=replicate(' ',10), /noerase, charsize=charsize_5;, $
;     xregion=xregion, yregion=yregion, xerrregion=xerrregion, yerrregion=yerrregion

    if keyword_set(postscript) then thisthick = 6.0 else thickthick = 2.0

    plot_kewley_grids, plotnumber=17, model=3, labeltype=4, /log12oh, /overplot, $
      /noZgrid, thick=3.0, Ugridcolor='dark green', Ulinestyle=0, $
      postscript=postscript, Ulabelcolor='black', Ucharsize=0.8, charthick=thisthick
    xyouts, 9.1, -0.18, 'log U', align=0.5, /data, charsize=0.9, charthick=thisthick

; overlay the three metallicity regions

    r12 = 8.14
    r23 = 8.68

    oplot, r12*[1,1], [-1.25,0.75], line=1, thick=postthick
    oplot, r23*[1,1], [-1.25,0.75], line=1, thick=postthick
 
    xyouts, 7.7, 0.65, 'R1', align=0.5, /data, charsize=charsize_5, charthick=postthick
    xyouts, 8.4, 0.65, 'R2', align=0.5, /data, charsize=charsize_5, charthick=postthick
    xyouts, 9.0, 0.65, 'R3', align=0.5, /data, charsize=charsize_5, charthick=postthick

    splog, 'Region 1, 2, 3: '
    region1 = where(xbig lt r12,nregion1)
    if (nregion1 gt 2L) then stats = im_stats(ybig[region1],sigrej=3.0,/verbose)

    regioniiha = where((xbig gt r12) and (xbig lt r23),nregioniiha)
    if (nregioniiha gt 2L) then stats = im_stats(ybig[regioniiha],sigrej=3.0,/no_head,/verbose)

    region3 = where(xbig gt r23,nregion3)
    if (nregion3 gt 2L) then stats = im_stats(ybig[region3],sigrej=3.0,/no_head,/verbose)

    im_openclose, postscript=postscript, /close    

; ------------------------------------------------------------
; 3-panel - L(B) vs [O II]/Ha - Integrated
; ------------------------------------------------------------

    psname = 'LB_vs_oiiha_multipanel'
    im_openclose, pspath+psname, postscript=postscript, /encapsulated

    pagemaker, nx=1, ny=3, yspace=0, xspace=0, xmargin=[2.0,1.0], $
      ymargin=[1.1,1.1], width=5.5, position=pos, /normal

    indx = where((atlasdust.oii_3727[0]/atlasdust.oii_3727[1] gt snrcut) and $
      (atlasnodust.b_lum_obs gt -900.0) and (atlasnodust.ehbha_err gt 0.0),nindx)

    x = atlasnodust[indx].b_lum_obs
    xerr = atlasnodust[indx].b_lum_obs_err
    xabs = atlasnodust[indx].m_b_obs

    indxnfgs = where((nfgsdust.oii_3727[0]/nfgsdust.oii_3727[1] gt snrcut) and $
      (nfgsnodust.b_lum_obs gt -900.0) and (nfgsnodust.ehbha_err gt 0.0),nindxnfgs)

    xnfgs = nfgsnodust[indxnfgs].b_lum_obs
    xerrnfgs = nfgsnodust[indxnfgs].b_lum_obs_err

    Aha = 1.0
    Aha_err = 0.1
    
    Aoii = Aha*k_lambda(3727,/odonnel)/k_lambda(6563,/odonnel)
    Aoii_err = Aha_err*k_lambda(3727,/odonnel)/k_lambda(6563,/odonnel)

    xtitle = 'log L(B) [L(B)'+sunsymbol()+']'
;   xtitle = 'log [L(B)/L(B,'+im_sunsymbol()+')]'
    ytitle = 'log ([O II]/H\alpha)'

    xrange = LBrange
    yrange = oiihacorrange

; ##########################################################
; Panel 1: A([O II])=0, A(Ha)=0
; ##########################################################

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

    stats = im_stats([y,ynfgs],/verbose)
    xstr = strtrim(string(stats.median,format='(F12.2)'),2)+'\pm'+strtrim(string(stats.sig68mean,format='(F12.2)'),2)
    rcor = r_correlate([x,xnfgs],[y,ynfgs],zd=zd,probd=probd)
    splog, 'Spearman rank test: ', rcor, probd

    ytitle = 'log ([O II]/H\alpha)_{obs}'

    atlas1d_lineplot, x, y, xerr, yerr, plottype=3, postscript=postscript, $
      ytitle=ytitle, xrange=xrange, yrange=yrange, legendtype=0, $
      xstyle=11, /right, /top, xnfgs=xnfgs, ynfgs=ynfgs, xerrnfgs=xerrnfgs, $
      yerrnfgs=yerrnfgs, position=pos[*,0], charsize=charsize_3, xtickname=replicate(' ',10)
    axis, /xaxis, xrange=interpol(xabs,x,!x.crange), xthick=postthick, xstyle=1, $
      charsize=charsize_3, charthick=postthick
    xyouts, mean(!x.crange), !y.crange[1]+0.16*(!y.crange[1]-!y.crange[0]), textoidl('M_{B} [mag]'), $
      align=0.5, charsize=charsize_3, charthick=postthick

    legend, '(a)', /left, /top, box=0, charsize=charsize_3, charthick=postthick
    legend, textoidl(xstr), /right, /top, box=0, charsize=charsize_3, charthick=postthick

; ##########################################################
; Panel 2: A([O II])=0, Individual A(Ha)
; ##########################################################

    y1 = atlasdust[indx].oii_3727[0]
    y1err = atlasdust[indx].oii_3727[1]

    y2 = atlasnodust[indx].h_alpha[0]
    y2err = atlasnodust[indx].h_alpha[1]

    y = alog10(y1/y2)
    yerr = im_compute_error(y1,y1err,y2,y2err,/log)

    y1nfgs = nfgsdust[indxnfgs].oii_3727[0]
    y1errnfgs = nfgsdust[indxnfgs].oii_3727[1]

    y2nfgs = nfgsnodust[indxnfgs].h_alpha[0]
    y2errnfgs = nfgsnodust[indxnfgs].h_alpha[1]

    ynfgs = alog10(y1nfgs/y2nfgs)
    yerrnfgs = im_compute_error(y1nfgs,y1errnfgs,y2nfgs,y2errnfgs,/log)

    stats = im_stats([y,ynfgs],/verbose)
    xstr = strtrim(string(stats.median,format='(F12.2)'),2)+'\pm'+strtrim(string(stats.sig68mean,format='(F12.2)'),2)
    rcor = r_correlate([x,xnfgs],[y,ynfgs],zd=zd,probd=probd)
    splog, 'Spearman rank test: ', rcor, probd

    ytitle = 'log ([O II]_{obs}/H\alpha_{cor})'

    atlas1d_lineplot, x, y, xerr, yerr, plottype=3, postscript=postscript, $
      xtitle='', ytitle=ytitle, xrange=xrange, yrange=yrange, legendtype=0, /noerase, $
      xstyle=3, /right, /top, xtickname=replicate(' ',10), $
      xnfgs=xnfgs, ynfgs=ynfgs, xerrnfgs=xerrnfgs, $
      yerrnfgs=yerrnfgs, position=pos[*,1], charsize=charsize_3
    
    legend, '(b)', /left, /top, box=0, charsize=charsize_3, charthick=postthick
    legend, textoidl(xstr), /right, /top, box=0, charsize=charsize_3, charthick=postthick

    w = where([x,xnfgs] lt 9)
    w2 = where([x,xnfgs] gt 10 and [x,xnfgs] lt 11)
;   junk = im_stats(([y,ynfgs])[w],/ver)
;   junk = im_stats(([y,ynfgs])[w2],/ver)

; ##########################################################
; Panel 3: Individual A([O II]), Individual A(Ha)
; ##########################################################

    y1 = atlasnodust[indx].oii_3727[0]
    y1err = atlasnodust[indx].oii_3727[1]

    y2 = atlasnodust[indx].h_alpha[0]
    y2err = atlasnodust[indx].h_alpha[1]

    y = alog10(y1/y2)
    yerr = im_compute_error(y1,y1err,y2,y2err,/log)

    y1nfgs = nfgsnodust[indxnfgs].oii_3727[0]
    y1errnfgs = nfgsnodust[indxnfgs].oii_3727[1]

    y2nfgs = nfgsnodust[indxnfgs].h_alpha[0]
    y2errnfgs = nfgsnodust[indxnfgs].h_alpha[1]

    ynfgs = alog10(y1nfgs/y2nfgs)
    yerrnfgs = im_compute_error(y1nfgs,y1errnfgs,y2nfgs,y2errnfgs,/log)

    stats = im_stats([y,ynfgs],/verbose)
    xstr = strtrim(string(stats.median,format='(F12.2)'),2)+'\pm'+strtrim(string(stats.sig68mean,format='(F12.2)'),2)
    rcor = r_correlate([x,xnfgs],[y,ynfgs],zd=zd,probd=probd)
    splog, 'Spearman rank test: ', rcor, probd

    ytitle = 'log ([O II]/H\alpha)_{cor}'

    atlas1d_lineplot, x, y, xerr, yerr, plottype=3, postscript=postscript, $
      xtitle=xtitle, ytitle=ytitle, xrange=xrange, yrange=yrange, legendtype=0, /noerase, $
      xstyle=3, /right, /top, xnfgs=xnfgs, ynfgs=ynfgs, xerrnfgs=xerrnfgs, $
      yerrnfgs=yerrnfgs, position=pos[*,2], charsize=charsize_3
    
    legend, '(c)', /left, /top, box=0, charsize=charsize_3, charthick=postthick
    legend, textoidl(xstr), /right, /top, box=0, charsize=charsize_3, charthick=postthick

    im_openclose, postscript=postscript, /close
    
; ------------------------------------------------------------
; visualize the H-alpha SFR calibrations
; ------------------------------------------------------------
    
    psname = 'visualize_sfr_ha'
    im_openclose, pspath+psname, postscript=postscript, /encapsulated

    irconst = 4.5D-44          ; K98 conversion L(IR) --> SFR(IR)
    haconst = 7.9-42           ; K98 conversion L(IR) --> SFR(IR)

; --------------------------------------------------
; L(Ha)/L(IR) individual extinction correction
; --------------------------------------------------

    log_ha_ir_median = -2.31 ; from the data
    
    ha_obs_median = alog10(irconst) - log_ha_ir_median ; L(Ha) --> SFR(IR)
    ha_obs_sigma = 0.19                              ; transformation scatter 

    xmax = ha_obs_median + 2 & xmin = ha_obs_median - 2 & dx = 0.01
    xaxis = findgen((xmax-xmin)/dx+1)*dx+xmin
    xrange = minmax(xaxis)
    
    yaxis = gauss1(xaxis,[ha_obs_median,ha_obs_sigma,1.0])
    yaxis = yaxis/total(yaxis)

    yrange = [0,max(yaxis)*1.15]

    xtitle = 'log [\psi(IR)/L(H\alpha)]'
    ytitle = 'Probability'

    djs_plot, [0], [0], /nodata, xrange=xrange, yrange=yrange, xtitle=xtitle, $
      ytitle=ytitle, charsize=2.0, charthick=postthick, xthick=postthick, $
      ythick=postthick, xsty=3, ysty=3
    
    djs_oplot, xaxis, yaxis, line=0, thick=postthick, color='dark green'

; --------------------------------------------------
; L(Ha)/L(IR) mean extinction correction
; --------------------------------------------------

    log_ha_ir_median = -2.19 ; from the data
    
    ha_obs_median = alog10(irconst) - log_ha_ir_median ; L(Ha) --> SFR(IR)
    ha_obs_sigma = 0.38                              ; transformation scatter 

    yaxis = gauss1(xaxis,[ha_obs_median,ha_obs_sigma,1.0])
    yaxis = yaxis/total(yaxis)

    djs_oplot, xaxis, yaxis, line=2, thick=postthick, color='red'

; --------------------------------------------------
; L(Ha)/L(IR) observed
; --------------------------------------------------

    log_ha_ir_median = -2.59 ; from the data
    
    ha_obs_median = alog10(irconst) - log_ha_ir_median ; L(Ha) --> SFR(IR)
    ha_obs_sigma = 0.38                              ; transformation scatter 

    yaxis = gauss1(xaxis,[ha_obs_median,ha_obs_sigma,1.0])
    yaxis = yaxis/total(yaxis)

    djs_oplot, xaxis, yaxis, line=3, thick=postthick, color='dark blue'

; --------------------------------------------------
; legend
; --------------------------------------------------

    label = ['A(H\alpha) = 0','Mean A(H\alpha)','Individual A(H\alpha)']
    color = ['dark blue','red','dark green']
    line = [3,2,1]
    legend, textoidl(label), /right, /top, box=0, charsize=1.5, $
      charthick=postthick, textcolor=djs_icolor(color), color=djs_icolor(color), $
      thick=postthick, line=line, spacing=0.07
    
    im_openclose, postscript=postscript, /close

; ------------------------------------------------------------
; ionization parameter vs [O II]_cor/Ha_cor
; ------------------------------------------------------------

; how does the [O II]/Ha ratio in our data vary with ionization
; parameter?  it forms a very narrow sequence.  if we can estimate the
; ionization parameter from [O III]/[O II] and then overlay the data
; on this grid...
    
    psname = 'logU_vs_oiiha_cor'
    im_openclose, pspath+psname, postscript=postscript, /encapsulated

    lineratio, atlasnodust, '', '', 'OII_3727', 'H_ALPHA', $
      dum1, dum2, y, yerr, index=indx, nindex=nindx, snrcut=snrcut, $
      xsyserr=syserr, ysyserr=syserr

    lineratio, nfgsnodust, '', '', 'OII_3727', 'H_ALPHA', $
      dum1, dum2, ynfgs, yerrnfgs, index=indxnfgs, snrcut=snrcut, $
      xsyserr=nfgssyserr, ysyserr=nfgssyserr
    ybig = [y,ynfgs]
    
; compute statistics in the y axis

    stats = im_stats(ybig,sigrej=3.0)

    xvert = [max(U),min(U),min(U),max(U)] ; ionization parameter range
    sigma1 = [$
      stats.median-1.0*stats.sigma,stats.median-1.0*stats.sigma,$
      stats.median+1.0*stats.sigma,stats.median+1.0*stats.sigma]
    sigma2 = [$
      stats.median-2.0*stats.sigma,stats.median-2.0*stats.sigma,$
      stats.median+2.0*stats.sigma,stats.median+2.0*stats.sigma]
    sigma3 = [$
      stats.median-3.0*stats.sigma,stats.median-3.0*stats.sigma,$
      stats.median+3.0*stats.sigma,stats.median+3.0*stats.sigma]

    xtitle = 'log U'
    ytitle = 'log ([O II]/H\alpha)_{cor}'

    xrange = [-4.0,-1.8]
    yrange = oiiharange

    plot_kewley_grids, plotnumber=15, model=3, labeltype=3, xrange=xrange, $
      yrange=yrange, /noUgrid, thick=postthick, Zgridcolor='dark red', $
      postscript=postscript, ymargin=[4,3]
    plot_kewley_grids, plotnumber=15, model=8, labeltype=2, /overplot, /noUgrid, $
      thick=postthick, Zgridcolor='purple', Zlinestyle=2, $
      postscript=postscript
    djs_plotlimitbox, xvert[0:1], [sigma2[[0,2]]], line=0, color='navy', $
      thick=postthick
    polyfill, xvert, sigma1, /data, color=djs_icolor('grey'), $
      linestyle=0, /line_fill, orientation=45, spacing=0.1
    legend, ['Model Starburst','Model H II Region'], line=[0,2], /left, /top, $
      charsize=charsize_3, charthick=postthick, box=0, color=djs_icolor(['dark red','purple']), $
      thick=postthick, clear=postscript

    im_openclose, postscript=postscript, /close    

; ------------------------------------------------------------
; logZ vs [O II]_cor/Ha_cor
; ------------------------------------------------------------

; how does the [O II]/Ha ratio in our data vary with abundance?
    
    psname = 'logZ_vs_oiiha_cor'
    im_openclose, pspath+psname, postscript=postscript, /encapsulated

    lineratio, atlasnodust, '', '', 'OII_3727', 'H_ALPHA', $
      dum1, dum2, y, yerr, index=indx, nindex=nindx, snrcut=snrcut, $
      xsyserr=syserr, ysyserr=syserr

    lineratio, nfgsnodust, '', '', 'OII_3727', 'H_ALPHA', $
      dum1, dum2, ynfgs, yerrnfgs, index=indxnfgs, snrcut=snrcut, $
      xsyserr=nfgssyserr, ysyserr=nfgssyserr
    ybig = [y,ynfgs]
    
; compute statistics in the y axis

    stats = im_stats(ybig,sigrej=3.0)

    xvert = alog10([min(Z),max(Z),max(Z),min(Z)]) ; metallicity range
    sigma1 = [$
      stats.median-1.0*stats.sigma,stats.median-1.0*stats.sigma,$
      stats.median+1.0*stats.sigma,stats.median+1.0*stats.sigma]
    sigma2 = [$
      stats.median-2.0*stats.sigma,stats.median-2.0*stats.sigma,$
      stats.median+2.0*stats.sigma,stats.median+2.0*stats.sigma]
    sigma3 = [$
      stats.median-3.0*stats.sigma,stats.median-3.0*stats.sigma,$
      stats.median+3.0*stats.sigma,stats.median+3.0*stats.sigma]

    xtitle = 'Z/Z'+sunsymbol()
    ytitle = 'log ([O II]/H\alpha)_{cor}'

    xrange = [-1.5,0.7]
    yrange = oiiharange

    xoh = alog10(Z)+Zsun_new
    coeff = [21.21,-2.29]
    good = where(xoh gt 8.4)
    oii_ha_kewley = alog10(coeff[0]+coeff[1]*xoh)
    xkewley = xoh-Zsun_new
    
    plot_kewley_grids, plotnumber=17, model=3, labeltype=4, xrange=xrange, $
      yrange=yrange, /noZgrid, thick=postthick, Ugridcolor='dark red', $
      xstyle=11, ymargin=[4,3], postscript=postscript, Umax=-2.5, Ulinestyle=0

;   djs_oplot, xkewley[good], oii_ha_kewley[good], line=1, thick=4.0, color='red'

    plot_kewley_grids, plotnumber=17, model=8, labeltype=5, /overplot, /noZgrid, $
      thick=postthick, Ugridcolor='purple', Ulinestyle=2,$
      postscript=postscript
    djs_plotlimitbox, xvert[0:1], [sigma2[[0,2]]], line=0, color='navy', $
      thick=postthick
    polyfill, xvert, sigma1, /data, color=djs_icolor('grey'), $
      linestyle=0, /line_fill, orientation=45, spacing=0.1
;   legend, ['Model Starburst','Model H II Region'], line=[0,2], /left, /top, $
;     charsize=charsize_3, charthick=postthick, box=0, color=djs_icolor(['dark red','purple']), $
;     thick=postthick, clear=postscript

    axis, /xaxis, xrange=interpol(xoh,alog10(Z),!x.crange), xthick=postthick, xstyle=1, $
      charsize=singlecharsize, charthick=postthick, /save
    xyouts, mean(!x.crange), !y.crange[1]+0.07*(!y.crange[1]-!y.crange[0]), textoidl('12 + log (O/H)'), $
      align=0.5, charsize=singlecharsize, charthick=postthick

    r12 = 8.14
    r23 = 8.68
    
    oplot, r12*[1,1], [-1.25,0.75], line=1, thick=postthick
    oplot, r23*[1,1], [-1.25,0.75], line=1, thick=postthick

    xyouts, 7.7, 0.65, 'region 1', align=0.5, /data, charsize=charsize, charthick=postthick
    xyouts, 8.4, 0.65, 'region 2', align=0.5, /data, charsize=charsize, charthick=postthick
    xyouts, 9.0, 0.65, 'region 3', align=0.5, /data, charsize=charsize, charthick=postthick

    im_openclose, postscript=postscript, /close    

; ------------------------------------------------------------
; OIIINII vs [O II]_cor/Ha_cor
; ------------------------------------------------------------

    psname = 'OIIINII_vs_oiiha_cor'
    im_openclose, pspath+psname, postscript=postscript, /encapsulated

    lineratio, atlasnodust, '', '', 'OII_3727', 'H_ALPHA', $
      dum1, dum2, y, yerr, index=indx, nindex=nindx, snrcut=snrcut, $
      xsyserr=syserr, ysyserr=syserr
    x = atlasnodust[indx].zstrong_oiiinii
    xerr = atlasnodust[indx].zstrong_oiiinii_err

    OHgood = where(atlasnodust[indx].zstrong_12oh_oiiinii_pettini gt -900)
    xgood = x[OHgood]
    xOHgood = atlasnodust[indx[OHgood]].zstrong_12oh_oiiinii_pettini

    xbig = x
    ybig = y
    
    lineratio, nfgsnodust, '', '', 'OII_3727', 'H_ALPHA', $
      dum1, dum2, ynfgs, yerrnfgs, index=indxnfgs, snrcut=snrcut, $
      xsyserr=nfgssyserr, ysyserr=nfgssyserr
    xnfgs = nfgsnodust[indxnfgs].zstrong_oiiinii
    xerrnfgs = nfgsnodust[indxnfgs].zstrong_oiiinii_err

    xbig = [xbig,xnfgs]
    ybig = [ybig,ynfgs]

    xtitle = 'log {([O III] \lambda5007/H\beta)/([N II] \lambda6584/H\alpha)}'
    ytitle = 'log ([O II]/H\alpha)_{cor}'

    xrange = reverse(oiiiniirange)
    yrange = oiiharange

    good = where((hii.zstrong_oiiinii gt -900.0) and (hii.oii_h_alpha gt -900.0))
    xregion = hii[good].zstrong_oiiinii & xerrregion = hii[good].zstrong_oiiinii_err
    yregion = hii[good].oii_h_alpha & yerrregion = hii[good].oii_h_alpha_err
       
    atlas1d_lineplot, x, y, xerr, yerr, plottype=3, postscript=postscript, $n
      xtitle=xtitle, ytitle=ytitle, xrange=xrange, yrange=yrange, legendtype=0, $
      xstyle=11, ymargin=[4,3], /left, /top, xregion=xregion, yregion=yregion, xerrregion=xerrregion, $
      yerrregion=yerrregion, xnfgs=xnfgs, ynfgs=ynfgs, xerrnfgs=xerrnfgs, yerrnfgs=yerrnfgs
    axis, /xaxis, xrange=interpol(xOHgood,xgood,!x.crange), xthick=postthick, xstyle=1, $
      charsize=singlecharsize, charthick=postthick
    xyouts, mean(!x.crange), !y.crange[1]+0.07*(!y.crange[1]-!y.crange[0]), textoidl('12 + log (O/H)'), $
      align=0.5, charsize=singlecharsize, charthick=postthick

; overlay the three metallicity regions defined in the NIIHA-OII/Ha plot

    OH12 = 8.14 & OH23 = 8.68
    r12 = interpol(xgood,xOHgood,OH12) & r23 = interpol(xgood,xOHgood,OH23)
    OH12 = interpol(xOHgood,xgood,r12) & OH23 = interpol(xOHgood,xgood,r23)
    Z12 = 10^(OH12 - Zsun_new) & Z23 = 10^(OH23 - Zsun_new)    

    splog, 'OIIINII = '+string(r12,format='(F4.1)')+' --> 12+log(O/H) = '+$
      string(OH12,format='(F4.2)')+', Z/Z_sun = '+string(Z12,format='(F4.2)')
    splog, 'OIIINII = '+string(r23,format='(F4.1)')+' --> 12+log(O/H) = '+$
      string(OH23,format='(F4.2)')+', Z/Z_sun = '+string(Z23,format='(F4.2)')
    
    oplot, r12*[1,1], [-1.25,0.75], line=1, thick=postthick
    oplot, r23*[1,1], [-1.25,0.75], line=1, thick=postthick
 
    xyouts, -0.5, 0.65, 'region 3', align=0.5, /data, charsize=charsize, charthick=postthick
    xyouts,  1.0, 0.65, 'region 2', align=0.5, /data, charsize=charsize, charthick=postthick
    xyouts,  2.5, 0.65, 'region 1', align=0.5, /data, charsize=charsize, charthick=postthick

    region1 = where(xbig gt r12,nregion1)
    splog, 'Region 1, 2, 3: '
    stats = im_stats(ybig[region1],sigrej=3.0,/verbose)

    regioniiha = where((xbig lt r12) and (xbig gt r23),nregioniiha)
    stats = im_stats(ybig[regioniiha],sigrej=3.0,/no_head,/verbose)

    region3 = where(xbig lt r23,nregion3)
    stats = im_stats(ybig[region3],sigrej=3.0,/no_head,/verbose)

    im_openclose, postscript=postscript, /close    

; ------------------------------------------------------------
; OIIINII vs [O II]_cor/Ha_cor
; ------------------------------------------------------------

    psname = 'sdss_OIIINII_vs_oiiha_cor'
    im_openclose, pspath+psname, postscript=postscript, /encapsulated

    lineratio, sdssnodust, '', '', 'OII_3727', 'H_ALPHA', $
      dum1, dum2, y, yerr, index=indx, nindex=nindx, snrcut=snrcut
    x = sdssnodust[indx].zstrong_oiiinii
    xerr = sdssnodust[indx].zstrong_oiiinii_err

    OHgood = where(sdssnodust[indx].zstrong_12oh_oiiinii_pettini gt -900)
    xgood = x[OHgood]
    xOHgood = sdssnodust[indx[OHgood]].zstrong_12oh_oiiinii_pettini

    xbig = x
    ybig = y

    xtitle = 'log {([O III] \lambda5007/H\beta)/([N II] \lambda6584/H\alpha)}'
    ytitle = 'log ([O II]/H\alpha)_{cor}'

    xrange = reverse(oiiiniirange)
    yrange = oiiharange

    sdss_lineplot, x, y, xerr, yerr, plottype=3, postscript=postscript, $
      xtitle=xtitle, ytitle=ytitle, xrange=xrange, yrange=yrange, legendtype=0, $
      xstyle=11, ymargin=[4,3], /left, /top, /xreverse
    axis, /xaxis, xrange=interpol(xOHgood,xgood,!x.crange), xthick=postthick, xstyle=1, $
      charsize=singlecharsize, charthick=postthick
    xyouts, mean(!x.crange), !y.crange[1]+0.07*(!y.crange[1]-!y.crange[0]), textoidl('12 + log (O/H)'), $
      align=0.5, charsize=singlecharsize, charthick=postthick

; overlay the three metallicity regions defined in the NIIHA-[O II]/Ha plot

    OH12 = 8.14 & OH23 = 8.68
    r12 = interpol(xgood,xOHgood,OH12) & r23 = interpol(xgood,xOHgood,OH23)
    OH12 = interpol(xOHgood,xgood,r12) & OH23 = interpol(xOHgood,xgood,r23)
    Z12 = 10^(OH12 - Zsun_new) & Z23 = 10^(OH23 - Zsun_new)    

    splog, 'OIIINII = '+string(r12,format='(F4.1)')+' --> 12+log(O/H) = '+$
      string(OH12,format='(F4.2)')+', Z/Z_sun = '+string(Z12,format='(F4.2)')
    splog, 'OIIINII = '+string(r23,format='(F4.1)')+' --> 12+log(O/H) = '+$
      string(OH23,format='(F4.2)')+', Z/Z_sun = '+string(Z23,format='(F4.2)')
    
    oplot, r12*[1,1], [-1.25,0.75], line=1, thick=postthick
    oplot, r23*[1,1], [-1.25,0.75], line=1, thick=postthick
 
    xyouts, -0.5, 0.65, 'region 3', align=0.5, /data, charsize=1.5, charthick=postthick
    xyouts,  1.0, 0.65, 'region 2', align=0.5, /data, charsize=1.5, charthick=postthick
    xyouts,  2.5, 0.65, 'region 1', align=0.5, /data, charsize=1.5, charthick=postthick

    region1 = where(xbig gt r12,nregion1)
    if (nregion1 ne 0L) then begin
       splog, 'Region 1, 2, 3: '
       stats = im_stats(ybig[region1],sigrej=3.0,/verbose)
    endif

    if (nregioniiha ne 0L) then begin
       regioniiha = where((xbig lt r12) and (xbig gt r23),nregioniiha)
       stats = im_stats(ybig[regioniiha],sigrej=3.0,/no_head,/verbose)
    endif

    if (nregion3 ne 0L) then begin
       region3 = where(xbig lt r23,nregion3)
       stats = im_stats(ybig[region3],sigrej=3.0,/no_head,/verbose)
    endif

    im_openclose, postscript=postscript, /close    

; ------------------------------------------------------------
; L(IR)/L(U) vs A(Ha)
; ------------------------------------------------------------
    
    psname = 'LIR_LU_vs_AHa'
    im_openclose, pspath+psname, postscript=postscript, /encapsulated

    indx = where((atlasdust.ir_lum gt -900) and (atlasdust.u_lum_obs gt -900) and $
      (atlasnodust.ebv_hahb_err gt 0.0),nindx)

    x = atlasdust[indx].ir_lum-atlasdust[indx].u_lum_obs
    xerr = sqrt(atlasdust[indx].ir_lum_err^2 + atlasdust[indx].u_lum_obs_err^2)
    y = atlasnodust[indx].ebv_hahb*k_lambda(6563.0,/odonnell)
    yerr = atlasnodust[indx].ebv_hahb_err*k_lambda(6563.0,/odonnell)
    
    indxnfgs = where((nfgsdust.ir_lum gt -900) and (nfgsdust.u_lum_obs gt -900) and $
      (nfgsnodust.ebv_hahb_err gt 0.0),nindxnfgs)

    xnfgs = nfgsdust[indxnfgs].ir_lum-nfgsdust[indxnfgs].u_lum_obs
    xerrnfgs = sqrt(nfgsdust[indxnfgs].ir_lum_err^2 + nfgsdust[indxnfgs].u_lum_obs_err^2)
    ynfgs = nfgsnodust[indxnfgs].ebv_hahb*k_lambda(6563.0,/odonnell)
    yerrnfgs = nfgsnodust[indxnfgs].ebv_hahb_err*k_lambda(6563.0,/odonnell)

    xrange = LIRLUrange
    yrange = AHarange

    xtitle = 'log [L(IR)/L(U)]'
    ytitle = 'A(H\alpha)'

    atlas1d_lineplot, x, y, xerr, yerr, plottype=3, postscript=postscript, $
      xtitle=xtitle, ytitle=ytitle, xrange=xrange, yrange=yrange, legendtype=0, $
      xnfgs=xnfgs, ynfgs=ynfgs, xerrnfgs=xerrnfgs, yerrnfgs=yerrnfgs

; overplot the expected relation assuming *all* the IR luminosity is
; due to absorption of the U-band flux    

    xaxis = findgen((xrange[1]-xrange[0])/0.01+1)*0.01+xrange[0]
    yaxis = 0.9208*alog(1+10^xaxis)
    
    djs_oplot, xaxis, 0.5*yaxis, line=0, thick=postthick
    djs_oplot, xaxis, 0.75*yaxis, line=2, thick=postthick
    djs_oplot, xaxis, 0.25*yaxis, line=2, thick=postthick

    im_openclose, postscript=postscript, /close

; ------------------------------------------------------------
; 2-panel [O II] plot - E(Hb-Ha) vs [O II]/Ha
; ------------------------------------------------------------

    psname = 'ehbha_vs_oiiha_multipanel'
    im_openclose, pspath+psname, postscript=postscript, /encapsulated

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
    yrange = oiihacorrange
    
    pagemaker, nx=1, ny=2, yspace=0, xspace=0, xmargin=[1.25,0.35], $
      ymargin=[1.1,1.1], position=pos, /normal

; ##########################################################
; Panel 1: [O II] observed, Ha observed
; ##########################################################

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

    stats = im_stats([y,ynfgs],/verbose)
    rcor = r_correlate([x,xnfgs],[y,ynfgs],zd=zd,probd=probd)
    splog, 'Spearman rank for E(Hb-Ha) vs [O II]/Ha observed: ', rcor, probd

;   ytitle = 'log ([O II]/H\alpha)_{obs}'

    atlas1d_lineplot, x, y, xerr, yerr, plottype=3, postscript=postscript, $
      ytitle=ytitle, xrange=xrange, yrange=yrange, legendtype=0, xminor=3, yminor=3, $
      xstyle=11, /right, /top, xnfgs=xnfgs, ynfgs=ynfgs, xerrnfgs=xerrnfgs, $
      yerrnfgs=yerrnfgs, position=pos[*,0], charsize=charsize, xtickname=replicate(' ',10)
    axis, /xaxis, xrange=interpol(xebv,x,!x.crange), xthick=postthick, xstyle=1, $
      charsize=charsize, charthick=postthick
    xyouts, mean(!x.crange), !y.crange[1]+0.14*(!y.crange[1]-!y.crange[0]), textoidl('E(B-V) [mag]'), $
      align=0.5, charsize=charsize, charthick=postthick

    legend, '(a)', /left, /top, box=0, charsize=charsize, charthick=postthick

    label = ['A([O II]) = 0','A(H\alpha) = 0']
    legend, textoidl(label), /left, /bottom, box=0, charsize=charsize_3, charthick=postthick

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

; ##########################################################
; Panel 2: [O II] corrected, Ha corrected
; ##########################################################

    y1 = atlasnodust[indx].oii_3727[0]
    y1err = atlasnodust[indx].oii_3727[1]

    y2 = atlasnodust[indx].h_alpha[0]
    y2err = atlasnodust[indx].h_alpha[1]

    y = alog10(y1/y2)
    yerr = im_compute_error(y1,y1err,y2,y2err,/log)

    y1nfgs = nfgsnodust[indxnfgs].oii_3727[0]
    y1errnfgs = nfgsnodust[indxnfgs].oii_3727[1]

    y2nfgs = nfgsnodust[indxnfgs].h_alpha[0]
    y2errnfgs = nfgsnodust[indxnfgs].h_alpha[1]

    ynfgs = alog10(y1nfgs/y2nfgs)
    yerrnfgs = im_compute_error(y1nfgs,y1errnfgs,y2nfgs,y2errnfgs,/log)

    stats = im_stats([y,ynfgs],/verbose)
    rcor = r_correlate([x,xnfgs],[y,ynfgs],zd=zd,probd=probd)
    splog, 'Spearman rank for E(Hb-Ha) vs [O II]/Ha corrected: ', rcor, probd

;   ytitle = 'log ([O II]/H\alpha)_{cor}'

    atlas1d_lineplot, x, y, xerr, yerr, plottype=3, postscript=postscript, $
      xtitle=xtitle, ytitle=ytitle, xrange=xrange, yrange=yrange, legendtype=0, /noerase, $
      xstyle=3, /right, /top, xnfgs=xnfgs, ynfgs=ynfgs, xerrnfgs=xerrnfgs, $
      yerrnfgs=yerrnfgs, position=pos[*,1], charsize=charsize, xminor=3, yminor=3
    
    legend, '(b)', /left, /top, box=0, charsize=charsize, charthick=postthick

    label = ['Individual A([O II])','Individual A(H\alpha)']
    legend, textoidl(label), /left, /bottom, box=0, charsize=charsize_3, charthick=postthick

    im_openclose, postscript=postscript, /close

; ------------------------------------------------------------
; 2-panel - E(Hb-Ha) vs [O II]/Ha
; ------------------------------------------------------------

    psname = 'sdss_ehbha_vs_oiiha_multipanel'
    im_openclose, pspath+psname, postscript=postscript, /encapsulated

    indx = where((sdssdust.oii_3727[0]/sdssdust.oii_3727[1] gt snrcut) and $
;     (sdssnodust.oii_3727[0]/sdssnodust.oii_3727[1] gt snrcut) and $
;     (sdssnodust.h_alpha[0]/sdssnodust.h_alpha[1] gt snrcut) and $
      (sdssnodust.ehbha_err gt 0.0),nindx)

    x = sdssnodust[indx].ehbha
    xerr = sdssnodust[indx].ehbha_err
    xebv = x/(k_lambda(4861,/odonnell)-k_lambda(6563,/odonnell))

    xtitle = 'E(H\beta-H\alpha) [mag]'
    ytitle = 'log ([O II]/H\alpha)'

    xrange = ehbharange
    yrange = oiihacorrange
    
    pagemaker, nx=1, ny=2, yspace=0, xspace=0, xmargin=[1.25,0.35], $
      ymargin=[1.1,1.1], position=pos, /normal

; ##########################################################
; Panel 1: [O II] observed, Ha observed
; ##########################################################

    y1 = sdssdust[indx].oii_3727[0]
    y1err = sdssdust[indx].oii_3727[1]

    y2 = sdssdust[indx].h_alpha[0]
    y2err = sdssdust[indx].h_alpha[1]

    y = alog10(y1/y2)
    yerr = im_compute_error(y1,y1err,y2,y2err,/log)

    stats = im_stats(y,/verbose)
    rcor = r_correlate(x,y,zd=zd,probd=probd)
    splog, 'Spearman rank for E(Hb-Ha) vs [O II]/Ha observed: ', rcor, probd

;   ytitle = 'log ([O II]/H\alpha)_{obs}'

    sdss_lineplot, x, y, xerr, yerr, plottype=3, postscript=postscript, $
      ytitle=ytitle, xrange=xrange, yrange=yrange, legendtype=0, xminor=3, yminor=3, $
      xstyle=11, /right, /top, position=pos[*,0], charsize=charsize, xtickname=replicate(' ',10)
    axis, /xaxis, xrange=interpol(xebv,x,!x.crange), xthick=postthick, xstyle=1, $
      charsize=charsize, charthick=postthick
    xyouts, mean(!x.crange), !y.crange[1]+0.13*(!y.crange[1]-!y.crange[0]), textoidl('E(B-V) [mag]'), $
      align=0.5, charsize=charsize, charthick=postthick

    legend, '(a)', /left, /top, box=0, charsize=charsize, charthick=postthick

    label = ['A([O II]) = 0','A(H\alpha) = 0']
    legend, textoidl(label), /left, /bottom, box=0, charsize=charsize_3, charthick=postthick

; overplot reddening relations

    ratio = 0.0

    x_ebv = findgen((1.3-(-0.1))/0.1)*0.1+(-0.1)
    yodonnell =  -0.4*x_ebv*(k_lambda(3727.0,/odonnell)-k_lambda(6563.0,/odonnell)) + ratio
    ysmc = -0.4*x_ebv*(k_lambda(3727,/smc)-k_lambda(6563,/smc)) + ratio

    djs_oplot, x_ebv, yodonnell, line=0, thick=postthick
    djs_oplot, x_ebv, ysmc, line=2, thick=postthick

;   label = ['Milky Way','SMC Bar']
;   linestyle = [0,2]
;   legend, label, /right, /top, box=0, linestyle=linestyle, $
;     charsize=charsize_3, charthick=postthick, thick=postthick, clear=postscript

; ##########################################################
; Panel 2: [O II] corrected, Ha corrected
; ##########################################################

    y1 = sdssnodust[indx].oii_3727[0]
    y1err = sdssnodust[indx].oii_3727[1]

    y2 = sdssnodust[indx].h_alpha[0]
    y2err = sdssnodust[indx].h_alpha[1]

    y = alog10(y1/y2)
    yerr = im_compute_error(y1,y1err,y2,y2err,/log)

    stats = im_stats(y,/verbose)
    rcor = r_correlate(x,y,zd=zd,probd=probd)
    splog, 'Spearman rank for E(Hb-Ha) vs [O II]/Ha corrected: ', rcor, probd

;   ytitle = 'log ([O II]/H\alpha)_{cor}'

    sdss_lineplot, x, y, xerr, yerr, plottype=3, postscript=postscript, $
      xtitle=xtitle, ytitle=ytitle, xrange=xrange, yrange=yrange, legendtype=0, /noerase, $
      xstyle=3, /right, /top, position=pos[*,1], charsize=charsize, xminor=3, yminor=3
    
    legend, '(b)', /left, /top, box=0, charsize=charsize, charthick=postthick

    label = ['Individual A([O II])','Individual A(H\alpha)']
    legend, textoidl(label), /left, /bottom, box=0, charsize=charsize_3, charthick=postthick

    im_openclose, postscript=postscript, /close

; ------------------------------------------------------------
; 4-panel [O II] plot - L(B) vs [O II]/Ha
; ------------------------------------------------------------

    psname = 'LB_vs_oiiha_multipanel'
    im_openclose, pspath+psname, postscript=postscript, /encapsulated

    indx = where((atlasdust.oii_3727[0]/atlasdust.oii_3727[1] gt snrcut) and $
      (atlasnodust.b_lum_obs gt -900.0) and (atlasnodust.ehbha_err gt 0.0),nindx)

    x = atlasnodust[indx].b_lum_obs
    xerr = atlasnodust[indx].b_lum_obs_err
    xabs = atlasnodust[indx].m_b_obs

    indxnfgs = where((nfgsdust.oii_3727[0]/nfgsdust.oii_3727[1] gt snrcut) and $
      (nfgsnodust.b_lum_obs gt -900.0) and (nfgsnodust.ehbha_err gt 0.0),nindxnfgs)

    xnfgs = nfgsnodust[indxnfgs].b_lum_obs
    xerrnfgs = nfgsnodust[indxnfgs].b_lum_obs_err

    Aha = 1.0
    Aha_err = 0.1
    
    Aoii = Aha*k_lambda(3727,/odonnel)/k_lambda(6563,/odonnel)
    Aoii_err = Aha_err*k_lambda(3727,/odonnel)/k_lambda(6563,/odonnel)

    xtitle = 'log L(B) [L(B)'+sunsymbol()+']'
    ytitle = 'log ([O II]/H\alpha)'

    xrange = LBrange
    yrange = oiihacorrange

    pagemaker, nx=2, ny=2, yspace=0, xspace=0, xmargin=[1.25,0.35], $
      ymargin=[0.8,1.1], position=pos, /normal

; ##########################################################
; Panel 1: A([O II])=0, A(Ha)=0
; ##########################################################

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

    stats = im_stats([y,ynfgs],/verbose)
    xstr = strtrim(string(stats.median,format='(F12.2)'),2)+'\pm'+strtrim(string(stats.sig68mean,format='(F12.2)'),2)
    rcor = r_correlate([x,xnfgs],[y,ynfgs],zd=zd,probd=probd)
    splog, 'Spearman rank for E(Hb-Ha) Panel a: ', rcor, probd

    atlas1d_lineplot, x, y, xerr, yerr, plottype=3, postscript=postscript, $
      ytitle=ytitle, xrange=xrange, yrange=yrange, legendtype=0, $
      xstyle=11, /right, /top, xnfgs=xnfgs, ynfgs=ynfgs, xerrnfgs=xerrnfgs, $
      yerrnfgs=yerrnfgs, position=pos[*,0], charsize=charsize, xtickname=replicate(' ',10)
    axis, /xaxis, xrange=interpol(xabs,x,!x.crange), xthick=postthick, xstyle=1, $
      charsize=charsize, charthick=postthick
    xyouts, mean(!x.crange), !y.crange[1]+0.12*(!y.crange[1]-!y.crange[0]), textoidl('M_{B} [mag]'), $
      align=0.5, charsize=charsize, charthick=postthick
;   djs_oplot, !x.crange, alog10(8.9D-42/2D-41)*[1,1], line=2, thick=postthick

    legend, '(a)', /left, /top, box=0, charsize=charsize, charthick=postthick
    legend, textoidl(xstr), /right, /top, box=0, charsize=charsize_3, charthick=postthick

    label = ['A([O II]) = 0','A(H\alpha) = 0']
    legend, textoidl(label), /left, /bottom, box=0, charsize=charsize_3, $
      charthick=postthick

; ##########################################################
; Panel 2: A([O II])=0, Mean A(Ha)
; ##########################################################

    y1 = atlasdust[indx].oii_3727[0]
    y1err = atlasdust[indx].oii_3727[1]
;   y1 = atlasdust[indx].oii_3727[0]*10^(0.4*Aoii)
;   y1err = atlasdust[indx].oii_3727[1]*10^(0.4*Aoii) + $
;     atlasdust[indx].oii_3727[0]*alog(10.0)*10^(0.4*Aoii)*Aoii_err

    y2 = atlasdust[indx].h_alpha[0]*10^(0.4*Aha)
    y2err = atlasdust[indx].h_alpha[1]*10^(0.4*Aha) + $
      atlasdust[indx].h_alpha[0]*alog(10.0)*10^(0.4*Aha)*Aha_err

    y = alog10(y1/y2)
    yerr = im_compute_error(y1,y1err,y2,y2err,/log)

    y1nfgs = nfgsdust[indxnfgs].oii_3727[0]
    y1errnfgs = nfgsdust[indxnfgs].oii_3727[1]
;   y1nfgs = nfgsdust[indxnfgs].oii_3727[0]*10^(0.4*Aoii)
;   y1errnfgs = nfgsdust[indxnfgs].oii_3727[1]*10^(0.4*Aoii) + $
;     nfgsdust[indxnfgs].oii_3727[0]*alog(10.0)*10^(0.4*Aoii)*Aoii_err

    y2nfgs = nfgsdust[indxnfgs].h_alpha[0]*10^(0.4*Aha)
    y2errnfgs = nfgsdust[indxnfgs].h_alpha[1]*10^(0.4*Aha) + $
      nfgsdust[indxnfgs].h_alpha[0]*alog(10.0)*10^(0.4*Aha)*Aha_err

    ynfgs = alog10(y1nfgs/y2nfgs)
    yerrnfgs = im_compute_error(y1nfgs,y1errnfgs,y2nfgs,y2errnfgs,/log)

    stats = im_stats([y,ynfgs],/verbose)
    xstr = strtrim(string(stats.median,format='(F12.2)'),2)+'\pm'+strtrim(string(stats.sig68mean,format='(F12.2)'),2)
    rcor = r_correlate([x,xnfgs],[y,ynfgs],zd=zd,probd=probd)
    splog, 'Spearman rank for E(Hb-Ha) Panel b: ', rcor, probd

    atlas1d_lineplot, x, y, xerr, yerr, plottype=3, postscript=postscript, $
      xrange=xrange, yrange=yrange, legendtype=0, /noerase, xtickname=replicate(' ',10), $
      xstyle=11, /right, /top, xnfgs=xnfgs, ynfgs=ynfgs, xerrnfgs=xerrnfgs, $
      yerrnfgs=yerrnfgs, position=pos[*,1], charsize=charsize, ytickname=replicate(' ',10)
    axis, /xaxis, xrange=interpol(xabs,x,!x.crange), xthick=postthick, xstyle=1, $
      charsize=charsize, charthick=postthick
    xyouts, mean(!x.crange), !y.crange[1]+0.12*(!y.crange[1]-!y.crange[0]), textoidl('M_{B} [mag]'), $
      align=0.5, charsize=charsize, charthick=postthick
    
    legend, '(b)', /left, /top, box=0, charsize=charsize, charthick=postthick
    legend, textoidl(xstr), /right, /top, box=0, charsize=charsize_3, charthick=postthick

    label = ['A([O II]) = 0','Mean A(H\alpha)']
;   label = ['Mean A([O II])','Mean A(H\alpha)']
    legend, textoidl(label), /left, /bottom, box=0, charsize=charsize_3, $
      charthick=postthick

; ##########################################################
; Panel 3: A([O II])=0, Individual A(Ha)
; ##########################################################

    y1 = atlasdust[indx].oii_3727[0]
    y1err = atlasdust[indx].oii_3727[1]

    y2 = atlasnodust[indx].h_alpha[0]
    y2err = atlasnodust[indx].h_alpha[1]

    y = alog10(y1/y2)
    yerr = im_compute_error(y1,y1err,y2,y2err,/log)

    y1nfgs = nfgsdust[indxnfgs].oii_3727[0]
    y1errnfgs = nfgsdust[indxnfgs].oii_3727[1]

    y2nfgs = nfgsnodust[indxnfgs].h_alpha[0]
    y2errnfgs = nfgsnodust[indxnfgs].h_alpha[1]

    ynfgs = alog10(y1nfgs/y2nfgs)
    yerrnfgs = im_compute_error(y1nfgs,y1errnfgs,y2nfgs,y2errnfgs,/log)

    stats = im_stats([y,ynfgs],/verbose)
    xstr = strtrim(string(stats.median,format='(F12.2)'),2)+'\pm'+strtrim(string(stats.sig68mean,format='(F12.2)'),2)
    rcor = r_correlate([x,xnfgs],[y,ynfgs],zd=zd,probd=probd)
    splog, 'Spearman rank for E(Hb-Ha) Panel c: ', rcor, probd

    atlas1d_lineplot, x, y, xerr, yerr, plottype=3, postscript=postscript, $
      xtitle=xtitle, ytitle=ytitle, xrange=xrange, yrange=yrange, legendtype=0, /noerase, $
      xstyle=3, /right, /top, xnfgs=xnfgs, ynfgs=ynfgs, xerrnfgs=xerrnfgs, $
      yerrnfgs=yerrnfgs, position=pos[*,2], charsize=charsize
    
    legend, '(c)', /left, /top, box=0, charsize=charsize, charthick=postthick
    legend, textoidl(xstr), /right, /top, box=0, charsize=charsize_3, charthick=postthick

    label = ['A([O II]) = 0','Individual A(H\alpha)']
    legend, textoidl(label), /left, /bottom, box=0, charsize=charsize_3, $
      charthick=postthick

; ##########################################################
; Panel 4: Individual A([O II]), Individual A(Ha)
; ##########################################################

    y1 = atlasnodust[indx].oii_3727[0]
    y1err = atlasnodust[indx].oii_3727[1]

    y2 = atlasnodust[indx].h_alpha[0]
    y2err = atlasnodust[indx].h_alpha[1]

    y = alog10(y1/y2)
    yerr = im_compute_error(y1,y1err,y2,y2err,/log)

    y1nfgs = nfgsnodust[indxnfgs].oii_3727[0]
    y1errnfgs = nfgsnodust[indxnfgs].oii_3727[1]

    y2nfgs = nfgsnodust[indxnfgs].h_alpha[0]
    y2errnfgs = nfgsnodust[indxnfgs].h_alpha[1]

    ynfgs = alog10(y1nfgs/y2nfgs)
    yerrnfgs = im_compute_error(y1nfgs,y1errnfgs,y2nfgs,y2errnfgs,/log)

    stats = im_stats([y,ynfgs],/verbose)
    xstr = strtrim(string(stats.median,format='(F12.2)'),2)+'\pm'+strtrim(string(stats.sig68mean,format='(F12.2)'),2)
    rcor = r_correlate([x,xnfgs],[y,ynfgs],zd=zd,probd=probd)
    splog, 'Spearman rank for E(Hb-Ha) Panel d: ', rcor, probd

    atlas1d_lineplot, x, y, xerr, yerr, plottype=3, postscript=postscript, $
      xtitle=xtitle, xrange=xrange, yrange=yrange, legendtype=0, /noerase, $
      xstyle=3, /right, /top, xnfgs=xnfgs, ynfgs=ynfgs, xerrnfgs=xerrnfgs, $
      yerrnfgs=yerrnfgs, position=pos[*,3], charsize=charsize, ytickname=replicate(' ',10)
    
    legend, '(d)', /left, /top, box=0, charsize=charsize, charthick=postthick
    legend, textoidl(xstr), /right, /top, box=0, charsize=charsize_3, charthick=postthick

    label = ['Individual A([O II])','Individual A(H\alpha)']
    legend, textoidl(label), /left, /bottom, box=0, charsize=charsize_3, $
      charthick=postthick

    im_openclose, postscript=postscript, /close
    
; ------------------------------------------------------------
; 4-panel [O II] plot - L(B) vs [O II]/Ha
; ------------------------------------------------------------

    psname = 'sdss_LB_vs_oiiha_multipanel'
    im_openclose, pspath+psname, postscript=postscript, /encapsulated

    indx = where((sdssdust.oii_3727[0]/sdssdust.oii_3727[1] gt snrcut) and $
      (sdssdust.b_lum gt -900.0) and (sdssnodust.ehbha_err gt 0),nindx)

    x = sdssdust[indx].b_lum
    xerr = sdssdust[indx].b_lum_err
    xabs = sdssdust[indx].m_b

    Aha = 1.0
    Aha_err = 0.1
    
    Aoii = Aha*k_lambda(3727,/odonnel)/k_lambda(6563,/odonnel)
    Aoii_err = Aha_err*k_lambda(3727,/odonnel)/k_lambda(6563,/odonnel)

    xtitle = 'log L(B) [L(B)'+sunsymbol()+']'
    ytitle = 'log ([O II]/H\alpha)'

    xrange = LBrange
    yrange = oiihacorrange

    pagemaker, nx=2, ny=2, yspace=0, xspace=0, xmargin=[1.25,0.35], $
      ymargin=[0.8,1.1], position=pos, /normal

; ##########################################################
; Panel 1: A([O II])=0, A(Ha)=0
; ##########################################################

    y1 = sdssdust[indx].oii_3727[0]
    y1err = sdssdust[indx].oii_3727[1]

    y2 = sdssdust[indx].h_alpha[0]
    y2err = sdssdust[indx].h_alpha[1]

    y = alog10(y1/y2)
    yerr = im_compute_error(y1,y1err,y2,y2err,/log)

    stats = im_stats(y,/verbose)
    xstr = strtrim(string(stats.median,format='(F12.2)'),2)+'\pm'+strtrim(string(stats.sig68mean,format='(F12.2)'),2)

    sdss_lineplot, x, y, xerr, yerr, plottype=3, postscript=postscript, $
      ytitle=ytitle, xrange=xrange, yrange=yrange, legendtype=0, $
      xstyle=11, /right, /top, position=pos[*,0], charsize=charsize, xtickname=replicate(' ',10)
    axis, /xaxis, xrange=interpol(xabs,x,!x.crange), xthick=postthick, xstyle=1, $
      charsize=charsize, charthick=postthick
    xyouts, mean(!x.crange), !y.crange[1]+0.12*(!y.crange[1]-!y.crange[0]), textoidl('M_{B} [mag]'), $
      align=0.5, charsize=charsize, charthick=postthick
;   djs_oplot, !x.crange, alog10(8.9D-42/2D-41)*[1,1], line=2, thick=postthick

    legend, '(a)', /left, /top, box=0, charsize=charsize, charthick=postthick
    legend, textoidl(xstr), /right, /top, box=0, charsize=charsize_3, charthick=postthick

    label = ['A([O II]) = 0','A(H\alpha) = 0']
    legend, textoidl(label), /left, /bottom, box=0, charsize=charsize_3, $
      charthick=postthick

; ##########################################################
; Panel 2: A([O II])=0, Mean A(Ha)
; ##########################################################

    y1 = sdssdust[indx].oii_3727[0]
    y1err = sdssdust[indx].oii_3727[1]

    y2 = sdssdust[indx].h_alpha[0]*10^(0.4*Aha)
    y2err = sdssdust[indx].h_alpha[1]*10^(0.4*Aha) + $
      sdssdust[indx].h_alpha[0]*alog(10.0)*10^(0.4*Aha)*Aha_err

    y = alog10(y1/y2)
    yerr = im_compute_error(y1,y1err,y2,y2err,/log)

    stats = im_stats(y,/verbose)
    xstr = strtrim(string(stats.median,format='(F12.2)'),2)+'\pm'+strtrim(string(stats.sig68mean,format='(F12.2)'),2)

    sdss_lineplot, x, y, xerr, yerr, plottype=3, postscript=postscript, $
      xrange=xrange, yrange=yrange, legendtype=0, /noerase, xtickname=replicate(' ',10), $
      xstyle=11, /right, /top, position=pos[*,1], charsize=charsize, ytickname=replicate(' ',10)
    axis, /xaxis, xrange=interpol(xabs,x,!x.crange), xthick=postthick, xstyle=1, $
      charsize=charsize, charthick=postthick
    xyouts, mean(!x.crange), !y.crange[1]+0.12*(!y.crange[1]-!y.crange[0]), textoidl('M_{B} [mag]'), $
      align=0.5, charsize=charsize, charthick=postthick
    
    legend, '(b)', /left, /top, box=0, charsize=charsize, charthick=postthick
    legend, textoidl(xstr), /right, /top, box=0, charsize=charsize_3, charthick=postthick

    label = ['A([O II]) = 0','Mean A(H\alpha)']
    legend, textoidl(label), /left, /bottom, box=0, charsize=charsize_3, $
      charthick=postthick

; ##########################################################
; Panel 3: A([O II])=0, Individual A(Ha)
; ##########################################################

    y1 = sdssdust[indx].oii_3727[0]
    y1err = sdssdust[indx].oii_3727[1]

    y2 = sdssnodust[indx].h_alpha[0]
    y2err = sdssnodust[indx].h_alpha[1]

    y = alog10(y1/y2)
    yerr = im_compute_error(y1,y1err,y2,y2err,/log)
    
    stats = im_stats(y,/verbose)
    xstr = strtrim(string(stats.median,format='(F12.2)'),2)+'\pm'+strtrim(string(stats.sig68mean,format='(F12.2)'),2)

    sdss_lineplot, x, y, xerr, yerr, plottype=3, postscript=postscript, $
      xtitle=xtitle, ytitle=ytitle, xrange=xrange, yrange=yrange, legendtype=0, /noerase, $
      xstyle=3, /right, /top, position=pos[*,2], charsize=charsize
    
    legend, '(c)', /left, /top, box=0, charsize=charsize, charthick=postthick
    legend, textoidl(xstr), /right, /top, box=0, charsize=charsize_3, charthick=postthick

    label = ['A([O II]) = 0','Individual A(H\alpha)']
    legend, textoidl(label), /left, /bottom, box=0, charsize=charsize_3, $
      charthick=postthick

; ##########################################################
; Panel 4: Individual A([O II]), Individual A(Ha)
; ##########################################################

    y1 = sdssnodust[indx].oii_3727[0]
    y1err = sdssnodust[indx].oii_3727[1]

    y2 = sdssnodust[indx].h_alpha[0]
    y2err = sdssnodust[indx].h_alpha[1]

    y = alog10(y1/y2)
    yerr = im_compute_error(y1,y1err,y2,y2err,/log)

    stats = im_stats(y,/verbose)
    xstr = strtrim(string(stats.median,format='(F12.2)'),2)+'\pm'+strtrim(string(stats.sig68mean,format='(F12.2)'),2)

    sdss_lineplot, x, y, xerr, yerr, plottype=3, postscript=postscript, $
      xtitle=xtitle, xrange=xrange, yrange=yrange, legendtype=0, /noerase, $
      xstyle=3, /right, /top, position=pos[*,3], charsize=charsize, ytickname=replicate(' ',10)
    
    legend, '(d)', /left, /top, box=0, charsize=charsize, charthick=postthick
    legend, textoidl(xstr), /right, /top, box=0, charsize=charsize_3, charthick=postthick

    label = ['Individual A([O II])','Individual A(H\alpha)']
    legend, textoidl(label), /left, /bottom, box=0, charsize=charsize_3, $
      charthick=postthick

    im_openclose, postscript=postscript, /close
    
; ------------------------------------------------------------
; 4-panel H-beta plot - L(B) vs Hb/Ha
; ------------------------------------------------------------
    
    psname = 'LB_vs_hbha_multipanel'
    im_openclose, pspath+psname, postscript=postscript, /encapsulated

    indx = where((atlasdust.b_lum_obs gt -900) and (atlasnodust.ebv_hahb_err gt 0.0) and $
      (atlasdust.h_beta_ew_uncor[0] ge ewcut),nindx)

    x = atlasdust[indx].b_lum_obs
    xabs = atlasdust[indx].m_b_obs
    xerr = atlasdust[indx].b_lum_obs_err
    
    indxnfgs = where((nfgsdust.b_lum_obs gt -900) and (nfgsnodust.ebv_hahb_err gt 0.0) and $
      (nfgsdust.h_beta_ew_uncor[0] ge ewcut),nindxnfgs)

    xnfgs = nfgsnodust[indxnfgs].b_lum_obs
    xerrnfgs = nfgsnodust[indxnfgs].b_lum_obs_err
    
    Aha = 1.0
    Aha_err = 0.1

    Ahb = Aha*k_lambda(4861,/odonnell)/k_lambda(6563,/odonnell)
    Ahb_err = 0.1*k_lambda(4861,/odonnell)/k_lambda(6563,/odonnell)

    xrange = LBrange
    yrange = [-1.9,-0.1]

    xtitle = 'log L(B) [L(B)'+sunsymbol()+']'
    ytitle = 'log (H\beta/H\alpha)'

    pagemaker, nx=2, ny=2, yspace=0, xspace=0, xmargin=[1.25,0.35], $
      ymargin=[0.8,1.1], position=pos, /normal

; ##########################################################
; Panel 1: No absorption correction, A(Hb)=0, A(Ha)=0
; ##########################################################

    hb = atlasdust[indx].h_beta_uncor[0]
    hb_err = atlasdust[indx].h_beta_uncor[1]
    
    ha = atlasdust[indx].h_alpha[0]
    ha_err = atlasdust[indx].h_alpha[1]
    
    y = alog10(hb/ha)
    yerr = im_compute_error(hb,hb_err,ha,ha_err,/log)

    hbnfgs = nfgsdust[indxnfgs].h_beta_uncor[0]
    hbnfgs_err = nfgsdust[indxnfgs].h_beta_uncor[1]
    
    hanfgs = nfgsdust[indxnfgs].h_alpha[0]
    hanfgs_err = nfgsdust[indxnfgs].h_alpha[1]
    
    ynfgs = alog10(hbnfgs/hanfgs)
    yerrnfgs = im_compute_error(hbnfgs,hbnfgs_err,hanfgs,hanfgs_err,/log)

    stats = im_stats([y,ynfgs],/verbose)
    xstr = strtrim(string(stats.median,format='(F12.2)'),2)+'\pm'+strtrim(string(stats.sig68mean,format='(F12.2)'),2)

    atlas1d_lineplot, x, y, xerr, yerr, plottype=3, postscript=postscript, $
      xtitle='', ytitle=ytitle, xrange=xrange, yrange=yrange, legendtype=0, $
      charsize=charsize, xstyle=11, xtickname=replicate(' ',10), position=pos[*,0], $
      xnfgs=xnfgs, ynfgs=ynfgs, xerrnfgs=xerrnfgs, yerrnfgs=yerrnfgs
    axis, /xaxis, xrange=interpol(xabs,x,!x.crange), xthick=postthick, xstyle=1, $
      charsize=charsize, charthick=postthick
    xyouts, mean(!x.crange), !y.crange[1]+0.11*(!y.crange[1]-!y.crange[0]), textoidl('M_{B} [mag]'), $
      align=0.5, charsize=charsize, charthick=postthick
    djs_oplot, !x.crange, alog10(1.0/HaHb)*[1,1], line=2, thick=postthick

    legend, '(a)', /left, /top, box=0, charsize=charsize, charthick=postthick
    legend, textoidl(xstr), /right, /top, box=0, charsize=charsize_3, charthick=postthick

    label = ['No Absorption Correction','A(H\beta) = 0, A(H\alpha) = 0']
    legend, textoidl(label), /left, /bottom, box=0, charsize=charsize_3, $
      charthick=postthick

; ############################################################
; Panel 2: Absorption-Corrected, A(Hb)=0, A(Ha)=0
; ############################################################

    hb = atlasdust[indx].h_beta[0]
    hb_err = atlasdust[indx].h_beta[1]
    
    ha = atlasdust[indx].h_alpha[0]
    ha_err = atlasdust[indx].h_alpha[1]
    
    y = alog10(hb/ha)
    yerr = im_compute_error(hb,hb_err,ha,ha_err,/log)

    hbnfgs = nfgsdust[indxnfgs].h_beta[0]
    hbnfgs_err = nfgsdust[indxnfgs].h_beta[1]
    
    hanfgs = nfgsdust[indxnfgs].h_alpha[0]
    hanfgs_err = nfgsdust[indxnfgs].h_alpha[1]
    
    ynfgs = alog10(hbnfgs/hanfgs)
    yerrnfgs = im_compute_error(hbnfgs,hbnfgs_err,hanfgs,hanfgs_err,/log)

    stats = im_stats([y,ynfgs],/verbose)
    xstr = strtrim(string(stats.median,format='(F12.2)'),2)+'\pm'+strtrim(string(stats.sig68mean,format='(F12.2)'),2)

    atlas1d_lineplot, x, y, xerr, yerr, plottype=3, postscript=postscript, $
      xtitle='', ytitle='', xrange=xrange, yrange=yrange, legendtype=0, $
      charsize=charsize, xstyle=11, xtickname=replicate(' ',10), ytickname=replicate(' ',10), $
      position=pos[*,1], /noerase, xnfgs=xnfgs, ynfgs=ynfgs, xerrnfgs=xerrnfgs, yerrnfgs=yerrnfgs
    axis, /xaxis, xrange=interpol(xabs,x,!x.crange), xthick=postthick, xstyle=1, $
      charsize=charsize, charthick=postthick
    xyouts, mean(!x.crange), !y.crange[1]+0.11*(!y.crange[1]-!y.crange[0]), textoidl('M_{B} [mag]'), $
      align=0.5, charsize=charsize, charthick=postthick
    djs_oplot, !x.crange, alog10(1.0/HaHb)*[1,1], line=2, thick=postthick

    legend, '(b)', /left, /top, box=0, charsize=charsize, charthick=postthick

    label = ['Absorption-Corrected','A(H\beta) = 0, A(H\alpha) = 0']
    legend, textoidl(xstr), /right, /top, box=0, charsize=charsize_3, charthick=postthick
    legend, textoidl(label), /left, /bottom, box=0, charsize=charsize_3, $
      charthick=postthick

; ############################################################
; Panel 3: Absorption-Corrected, A(Hb)=0, Mean A(Ha)
; ############################################################

    hb = atlasdust[indx].h_beta[0];*10^(0.4*Ahb)
    hb_err = atlasdust[indx].h_beta[1];*10^(0.4*Ahb)+atlasdust[indx].h_beta[0]*alog(10.0)*10^(0.4*Ahb)*Ahb_err

    ha = atlasdust[indx].h_alpha[0]*10^(0.4*Aha)
    ha_err = atlasdust[indx].h_alpha[1]*10^(0.4*Aha) + $
      atlasdust[indx].h_alpha[0]*alog(10.0)*10^(0.4*Aha)*Aha_err

    y = alog10(hb/ha)
    yerr = im_compute_error(hb,hb_err,ha,ha_err,/log)

    hbnfgs = nfgsdust[indxnfgs].h_beta[0];*10^(0.4*Ahb)
    hbnfgs_err = nfgsdust[indxnfgs].h_beta[1];*10^(0.4*Ahb)+nfgsdust[indxnfgs].h_beta[0]*alog(10.0)*10^(0.4*Ahb)*Ahb_err

    hanfgs = nfgsdust[indxnfgs].h_alpha[0]*10^(0.4*Aha)
    hanfgs_err = nfgsdust[indxnfgs].h_alpha[1]*10^(0.4*Aha) + $
      nfgsdust[indxnfgs].h_alpha[0]*alog(10.0)*10^(0.4*Aha)*Aha_err

    ynfgs = alog10(hbnfgs/hanfgs)
    yerrnfgs = im_compute_error(hbnfgs,hbnfgs_err,hanfgs,hanfgs_err,/log)

    stats = im_stats([y,ynfgs],/verbose)
    xstr = strtrim(string(stats.median,format='(F12.2)'),2)+'\pm'+strtrim(string(stats.sig68mean,format='(F12.2)'),2)
    
    atlas1d_lineplot, x, y, xerr, yerr, plottype=3, postscript=postscript, $
      xtitle=xtitle, ytitle=ytitle, xrange=xrange, yrange=yrange, legendtype=0, $
      charsize=charsize, position=pos[*,2], /noerase, xnfgs=xnfgs, ynfgs=ynfgs, $
      xerrnfgs=xerrnfgs, yerrnfgs=yerrnfgs
    djs_oplot, !x.crange, alog10(1.0/HaHb)*[1,1], line=2, thick=postthick

    legend, '(c)', /left, /top, box=0, charsize=charsize, charthick=postthick
    legend, textoidl(xstr), /right, /top, box=0, charsize=charsize_3, charthick=postthick

    label = ['Absorption-Corrected','A(H\beta) = 0, Mean A(H\alpha)']
;   label = ['Absorption-Corrected','Mean A(H\beta), Mean A(H\alpha)']
    legend, textoidl(label), /left, /bottom, box=0, charsize=charsize_3, $
      charthick=postthick
    
; ############################################################
; Panel 4: Absorption-Corrected, A(Hb)=0, Mean A(Ha)
; ############################################################

    hb = atlasdust[indx].h_beta[0];*10^(0.4*Ahb)
    hb_err = atlasdust[indx].h_beta[1];*10^(0.4*Ahb)+atlasdust[indx].h_beta[0]*alog(10.0)*10^(0.4*Ahb)*Ahb_err

    ha = atlasnodust[indx].h_alpha[0]
    ha_err = atlasnodust[indx].h_alpha[1]

    y = alog10(hb/ha)
    yerr = im_compute_error(hb,hb_err,ha,ha_err,/log)

    hbnfgs = nfgsdust[indxnfgs].h_beta[0];*10^(0.4*Ahb)
    hbnfgs_err = nfgsdust[indxnfgs].h_beta[1];*10^(0.4*Ahb)+nfgsdust[indxnfgs].h_beta[0]*alog(10.0)*10^(0.4*Ahb)*Ahb_err

    hanfgs = nfgsnodust[indxnfgs].h_alpha[0]
    hanfgs_err = nfgsnodust[indxnfgs].h_alpha[1]

    ynfgs = alog10(hbnfgs/hanfgs)
    yerrnfgs = im_compute_error(hbnfgs,hbnfgs_err,hanfgs,hanfgs_err,/log)
       
    stats = im_stats([y,ynfgs],/verbose)
    xstr = strtrim(string(stats.median,format='(F12.2)'),2)+'\pm'+strtrim(string(stats.sig68mean,format='(F12.2)'),2)
        
    atlas1d_lineplot, x, y, xerr, yerr, plottype=3, postscript=postscript, $
      xtitle=xtitle, ytitle='', xrange=xrange, yrange=yrange, legendtype=0, $
      charsize=charsize, ytickname=replicate(' ',10), position=pos[*,3], /noerase, $
      xnfgs=xnfgs, ynfgs=ynfgs, xerrnfgs=xerrnfgs, yerrnfgs=yerrnfgs
    djs_oplot, !x.crange, alog10(1.0/HaHb)*[1,1], line=2, thick=postthick

    legend, '(d)', /left, /top, box=0, charsize=charsize, charthick=postthick
    legend, textoidl(xstr), /right, /top, box=0, charsize=charsize_3, charthick=postthick

    label = ['Absorption-Corrected','A(H\beta) = 0, Individual A(H\alpha)']
;   label = ['Absorption-Corrected','Mean A(H\beta), Individual A(H\alpha)']
    legend, textoidl(label), /left, /bottom, box=0, charsize=charsize_3, $
      charthick=postthick

    im_openclose, postscript=postscript, /close

; ------------------------------------------------------------
; 4-panel H-beta plot - L(B) vs Hb/Ha
; ------------------------------------------------------------
    
    psname = 'sdss_LB_vs_hbha_multipanel'
    im_openclose, pspath+psname, postscript=postscript, /encapsulated

    indx = where((sdssdust.b_lum gt -900) and (sdssnodust.ebv_hahb_err gt 0.0) and $
      (sdssdust.h_beta_ew_uncor[0] ge ewcut),nindx)

    x = sdssdust[indx].b_lum
    xabs = sdssdust[indx].m_b
    xerr = sdssdust[indx].b_lum_err
    
    Aha = 1.0
    Aha_err = 0.1

    Ahb = Aha*k_lambda(4861,/odonnell)/k_lambda(6563,/odonnell)
    Ahb_err = 0.1*k_lambda(4861,/odonnell)/k_lambda(6563,/odonnell)

    xrange = LBrange
    yrange = [-1.9,-0.1]

    xtitle = 'log L(B) [L(B)'+sunsymbol()+']'
    ytitle = 'log (H\beta/H\alpha)'

    pagemaker, nx=2, ny=2, yspace=0, xspace=0, xmargin=[1.25,0.35], $
      ymargin=[0.8,1.1], position=pos, /normal

; ##########################################################
; Panel 1: No absorption correction, A(Hb)=0, A(Ha)=0
; ##########################################################

    hb = sdssdust[indx].h_beta_uncor[0]
    hb_err = sdssdust[indx].h_beta_uncor[1]
    
    ha = sdssdust[indx].h_alpha[0]
    ha_err = sdssdust[indx].h_alpha[1]
    
    y = alog10(hb/ha)
    yerr = im_compute_error(hb,hb_err,ha,ha_err,/log)

    stats = im_stats(y,/verbose)
    xstr = strtrim(string(stats.median,format='(F12.2)'),2)+'\pm'+strtrim(string(stats.sig68mean,format='(F12.2)'),2)

    sdss_lineplot, x, y, xerr, yerr, plottype=3, postscript=postscript, $
      xtitle='', ytitle=ytitle, xrange=xrange, yrange=yrange, legendtype=0, $
      charsize=charsize, xstyle=11, xtickname=replicate(' ',10), position=pos[*,0]
    axis, /xaxis, xrange=interpol(xabs,x,!x.crange), xthick=postthick, xstyle=1, $
      charsize=charsize, charthick=postthick
    xyouts, mean(!x.crange), !y.crange[1]+0.11*(!y.crange[1]-!y.crange[0]), textoidl('M_{B} [mag]'), $
      align=0.5, charsize=charsize, charthick=postthick
    djs_oplot, !x.crange, alog10(1.0/HaHb)*[1,1], line=2, thick=postthick

    legend, '(a)', /left, /top, box=0, charsize=charsize, charthick=postthick
    legend, textoidl(xstr), /right, /top, box=0, charsize=charsize_3, charthick=postthick

    label = ['No Absorption Correction','A(H\beta) = 0, A(H\alpha) = 0']
    legend, textoidl(label), /left, /bottom, box=0, charsize=charsize_3, $
      charthick=postthick

; ############################################################
; Panel 2: Absorption-Corrected, A(Hb)=0, A(Ha)=0
; ############################################################

    hb = sdssdust[indx].h_beta[0]
    hb_err = sdssdust[indx].h_beta[1]
    
    ha = sdssdust[indx].h_alpha[0]
    ha_err = sdssdust[indx].h_alpha[1]
    
    y = alog10(hb/ha)
    yerr = im_compute_error(hb,hb_err,ha,ha_err,/log)

    stats = im_stats(y,/verbose)
    xstr = strtrim(string(stats.median,format='(F12.2)'),2)+'\pm'+strtrim(string(stats.sig68mean,format='(F12.2)'),2)

    sdss_lineplot, x, y, xerr, yerr, plottype=3, postscript=postscript, $
      xtitle='', ytitle='', xrange=xrange, yrange=yrange, legendtype=0, $
      charsize=charsize, xstyle=11, xtickname=replicate(' ',10), ytickname=replicate(' ',10), $
      position=pos[*,1], /noerase
    axis, /xaxis, xrange=interpol(xabs,x,!x.crange), xthick=postthick, xstyle=1, $
      charsize=charsize, charthick=postthick
    xyouts, mean(!x.crange), !y.crange[1]+0.11*(!y.crange[1]-!y.crange[0]), textoidl('M_{B} [mag]'), $
      align=0.5, charsize=charsize, charthick=postthick
    djs_oplot, !x.crange, alog10(1.0/HaHb)*[1,1], line=2, thick=postthick

    legend, '(b)', /left, /top, box=0, charsize=charsize, charthick=postthick
    legend, textoidl(xstr), /right, /top, box=0, charsize=charsize_3, charthick=postthick

    label = ['Absorption-Corrected','A(H\beta) = 0, A(H\alpha) = 0']
    legend, textoidl(label), /left, /bottom, box=0, charsize=charsize_3, $
      charthick=postthick
    
; ############################################################
; Panel 3: Absorption-Corrected, Mean A(Hb), Mean A(Ha)
; ############################################################

    hb = sdssdust[indx].h_beta[0];*10^(0.4*Ahb)
    hb_err = sdssdust[indx].h_beta[1];*10^(0.4*Ahb)+sdssdust[indx].h_beta[0]*alog(10.0)*10^(0.4*Ahb)*Ahb_err

    ha = sdssdust[indx].h_alpha[0]*10^(0.4*Aha)
    ha_err = sdssdust[indx].h_alpha[1]*10^(0.4*Aha) + $
      sdssdust[indx].h_alpha[0]*alog(10.0)*10^(0.4*Aha)*Aha_err

    y = alog10(hb/ha)
    yerr = im_compute_error(hb,hb_err,ha,ha_err,/log)

    stats = im_stats(y,/verbose)
    xstr = strtrim(string(stats.median,format='(F12.2)'),2)+'\pm'+strtrim(string(stats.sig68mean,format='(F12.2)'),2)
    
    sdss_lineplot, x, y, xerr, yerr, plottype=3, postscript=postscript, $
      xtitle=xtitle, ytitle=ytitle, xrange=xrange, yrange=yrange, legendtype=0, $
      charsize=charsize, position=pos[*,2], /noerase
    djs_oplot, !x.crange, alog10(1.0/HaHb)*[1,1], line=2, thick=postthick

    legend, '(c)', /left, /top, box=0, charsize=charsize, charthick=postthick
    legend, textoidl(xstr), /right, /top, box=0, charsize=charsize_3, charthick=postthick

    label = ['Absorption-Corrected','A(H\beta) = 0, Mean A(H\alpha)']
    legend, textoidl(label), /left, /bottom, box=0, charsize=charsize_3, $
      charthick=postthick
    
; ###############################################################
; Panel 4: Absorption-Corrected, A(Hb)=0, Individual A(Ha)
; ###############################################################

    hb = sdssdust[indx].h_beta[0];*10^(0.4*Ahb)
    hb_err = sdssdust[indx].h_beta[1];*10^(0.4*Ahb)+sdssdust[indx].h_beta[0]*alog(10.0)*10^(0.4*Ahb)*Ahb_err

    ha = sdssnodust[indx].h_alpha[0]
    ha_err = sdssnodust[indx].h_alpha[1]

    y = alog10(hb/ha)
    yerr = im_compute_error(hb,hb_err,ha,ha_err,/log)

    stats = im_stats(y,/verbose)
    xstr = strtrim(string(stats.median,format='(F12.2)'),2)+'\pm'+strtrim(string(stats.sig68mean,format='(F12.2)'),2)

    sdss_lineplot, x, y, xerr, yerr, plottype=3, postscript=postscript, $
      xtitle=xtitle, ytitle='', xrange=xrange, yrange=yrange, legendtype=0, $
      charsize=charsize, ytickname=replicate(' ',10), position=pos[*,3], /noerase
    djs_oplot, !x.crange, alog10(1.0/HaHb)*[1,1], line=2, thick=postthick

    legend, '(d)', /left, /top, box=0, charsize=charsize, charthick=postthick
    legend, textoidl(xstr), /right, /top, box=0, charsize=charsize_3, charthick=postthick

    label = ['Absorption-Corrected','A(H\beta) = 0, Individual A(H\alpha)']
    legend, textoidl(label), /left, /bottom, box=0, charsize=charsize_3, $
      charthick=postthick

    im_openclose, postscript=postscript, /close

; ------------------------------------------------------------
; L(IR)/L([O II]) vs A(Ha)
; ------------------------------------------------------------
    
    psname = 'LIR_LOII_vs_AHa'
    im_openclose, pspath+psname, postscript=postscript, /encapsulated

    indx = where((atlasdust.ir_flux gt -900) and (atlasdust.oii_3727[0]/atlasdust.oii_3727[1] gt snrcut) and $
      (atlasnodust.ebv_hahb_err gt 0.0),nindx)

    x1 = atlasdust[indx].ir_flux
    x1err = atlasdust[indx].ir_flux_err
    x2 = atlasdust[indx].oii_3727[0]
    x2err = atlasdust[indx].oii_3727[1]

    x = alog10(x1/x2)
    xerr = im_compute_error(x1,x1err,x2,x2err,/log)

    y = atlasnodust[indx].ebv_hahb*k_lambda(6563.0,/odonnell)
    yerr = atlasnodust[indx].ebv_hahb_err*k_lambda(6563.0,/odonnell)
    
    indxnfgs = where((nfgsdust.ir_flux gt -900) and (nfgsdust.oii_3727[0]/nfgsdust.oii_3727[1] gt snrcut) and $
      (nfgsnodust.ebv_hahb_err gt 0.0),nindxnfgs)

    x1nfgs = nfgsdust[indxnfgs].ir_flux
    x1errnfgs = nfgsdust[indxnfgs].ir_flux_err
    x2nfgs = nfgsdust[indxnfgs].oii_3727[0]
    x2errnfgs = nfgsdust[indxnfgs].oii_3727[1]
    
    xnfgs = alog10(x1nfgs/x2nfgs)
    xerrnfgs = im_compute_error(x1nfgs,x1errnfgs,x2nfgs,x2errnfgs,/log)

    ynfgs = nfgsnodust[indxnfgs].ebv_hahb*k_lambda(6563.0,/odonnell)
    yerrnfgs = nfgsnodust[indxnfgs].ebv_hahb_err*k_lambda(6563.0,/odonnell)

    xrange = LIRLOIIrange
    yrange = AHarange

    xtitle = 'log (L_{IR}/L_{[O II], obs})'
    ytitle = 'A(H\alpha)'

    atlas1d_lineplot, x, y, xerr, yerr, plottype=3, postscript=postscript, $
      xtitle=xtitle, ytitle=ytitle, xrange=xrange, yrange=yrange, legendtype=0, $
      xnfgs=xnfgs, ynfgs=ynfgs, xerrnfgs=xerrnfgs, yerrnfgs=yerrnfgs

; overplot the expected relation assuming *all* the IR luminosity is
; due to absorption of the U-band flux    

    xaxis = findgen((xrange[1]-xrange[0])/0.01+1)*0.01+xrange[0]
    yaxis = 0.9208*alog(1+10^xaxis)
    
    djs_oplot, xaxis, 0.5*yaxis, line=0, thick=postthick
    djs_oplot, xaxis, 0.75*yaxis, line=2, thick=postthick
    djs_oplot, xaxis, 0.25*yaxis, line=2, thick=postthick
    
    im_openclose, postscript=postscript, /close

; ------------------------------------------------------------
; 2-panel H-alpha plot - L(B) vs L(Ha)/L(FIR)
; ------------------------------------------------------------
    
    psname = 'LB_vs_LHa_LFIR_multipanel'
    im_openclose, pspath+psname, postscript=postscript, /encapsulated

    indx = where((atlasdust.b_lum_obs gt -900) and (atlasnodust.fir_flux gt -900) and $
      (atlasnodust.ebv_hahb_err gt 0.0),nindx)

    x = atlasdust[indx].b_lum_obs
    xerr = atlasdust[indx].b_lum_obs_err
    xabs = atlasdust[indx].m_b_obs
    
    lir = atlasdust[indx].fir_flux ; [erg/s/cm2]
    lir_err = atlasdust[indx].fir_flux_err

    indxnfgs = where((nfgsdust.b_lum_obs gt -900) and (nfgsdust.fir_flux gt -900) and $
      (nfgsnodust.ebv_hahb_err gt 0.0),nindxnfgs)

    xnfgs = nfgsnodust[indxnfgs].b_lum_obs
    xerrnfgs = nfgsnodust[indxnfgs].b_lum_obs_err
    
    lirnfgs = nfgsdust[indxnfgs].fir_flux
    lirnfgs_err = nfgsdust[indxnfgs].fir_flux_err
    
    xrange = LBrange
    yrange = haLIRrange

    xtitle = 'log L(B) [L(B)'+sunsymbol()+']'
    ytitle = 'log [L(H\alpha)/L(FIR)]'

    lhalir = alog10(4.5D-44/7.9D-42)
    
    pagemaker, nx=1, ny=2, yspace=0, xspace=0, xmargin=[1.25,0.35], $
      ymargin=[1.1,1.1], position=pos, /normal

; ##########################################################
; Panel 1: No reddening correction
; ##########################################################

    ha = atlasdust[indx].h_alpha[0]
    ha_err = atlasdust[indx].h_alpha[1]
    
    y = alog10(ha/lir)
    yerr = im_compute_error(ha,ha_err,lir,lir_err,/log)

    hanfgs = nfgsdust[indxnfgs].h_alpha[0]
    hanfgs_err = nfgsdust[indxnfgs].h_alpha[1]
    
    ynfgs = alog10(hanfgs/lirnfgs)
    yerrnfgs = im_compute_error(hanfgs,hanfgs_err,lirnfgs,lirnfgs_err,/log)

    stats = im_stats([y,ynfgs],/verbose)
    xstr = strtrim(string(stats.median,format='(F12.2)'),2)+'\pm'+strtrim(string(stats.sig68mean,format='(F12.2)'),2)

    atlas1d_lineplot, x, y, xerr, yerr, plottype=3, postscript=postscript, $
      xtitle='', ytitle=ytitle, xrange=xrange, yrange=yrange, legendtype=0, $
      charsize=charsize, xstyle=11, position=pos[*,0], xnfgs=xnfgs, ynfgs=ynfgs, $
      xerrnfgs=xerrnfgs, yerrnfgs=yerrnfgs, xtickname=replicate(' ',10), yminor=3
;   djs_oplot, !x.crange, lhalir*[1,1], line=2, thick=postthick
    axis, /xaxis, xrange=interpol(xabs,x,!x.crange), xthick=postthick, xstyle=1, $
      charsize=charsize, charthick=postthick
    xyouts, mean(!x.crange), !y.crange[1]+0.14*(!y.crange[1]-!y.crange[0]), textoidl('M_{B} [mag]'), $
      align=0.5, charsize=charsize, charthick=postthick

    legend, '(a)', /left, /top, box=0, charsize=charsize, charthick=postthick

    label = ['A(H\alpha) = 0']
    legend, textoidl(xstr), /right, /top, box=0, charsize=charsize_3, charthick=postthick
    legend, textoidl(label), /left, /bottom, box=0, charsize=charsize_3, $
      charthick=postthick

; ############################################################
; Panel 2: Good Reddening Correction
; ############################################################

    ha = atlasnodust[indx].h_alpha[0]
    ha_err = atlasnodust[indx].h_alpha[1]

    y = alog10(ha/lir)
    yerr = im_compute_error(ha,ha_err,lir,lir_err,/log)

    hanfgs = nfgsnodust[indxnfgs].h_alpha[0]
    hanfgs_err = nfgsnodust[indxnfgs].h_alpha[1]

    ynfgs = alog10(hanfgs/lirnfgs)
    yerrnfgs = im_compute_error(hanfgs,hanfgs_err,lirnfgs,lirnfgs_err,/log)

    stats = im_stats([y,ynfgs],/verbose,/no_head)
    xstr = strtrim(string(stats.median,format='(F12.2)'),2)+'\pm'+strtrim(string(stats.sig68mean,format='(F12.2)'),2)

    atlas1d_lineplot, x, y, xerr, yerr, plottype=3, postscript=postscript, $
      xtitle=xtitle, ytitle=ytitle, xrange=xrange, yrange=yrange, legendtype=0, $
      charsize=charsize, position=pos[*,1], /noerase, xnfgs=xnfgs, ynfgs=ynfgs, $
      xerrnfgs=xerrnfgs, yerrnfgs=yerrnfgs, yminor=3
;   djs_oplot, !x.crange, lhalir*[1,1], line=2, thick=postthick

    legend, '(c)', /left, /top, box=0, charsize=charsize, charthick=postthick

    label = ['Individual A(H\alpha)']
    legend, textoidl(xstr), /right, /top, box=0, charsize=charsize_3, charthick=postthick
    legend, textoidl(label), /left, /bottom, box=0, charsize=charsize_3, $
      charthick=postthick

    im_openclose, postscript=postscript, /close

; ------------------------------------------------------------
; [N II]/Ha vs [O II]_cor/Ha_cor
; ------------------------------------------------------------

    psname = 'niiha_vs_oiiha_cor'
    im_openclose, pspath+psname, postscript=postscript, /encapsulated

    lineratio, atlasnodust, 'NII_6584', 'H_ALPHA', 'OII_3727', 'H_ALPHA', $
      x, xerr, y, yerr, index=indx, nindex=nindx, snrcut=snrcut, $
      xsyserr=syserr, ysyserr=syserr

    OHgood = where(atlasnodust[indx].zstrong_12oh_niiha gt -900)
    xgood = x[OHgood]
    xOHgood = atlasnodust[indx[OHgood]].zstrong_12oh_niiha
    
    xbig = x
    ybig = y

    if keyword_set(kenn92) then begin
       lineratio, kenn92nodust, 'NII_6584', 'H_ALPHA', 'OII_3727', 'H_ALPHA', $
         xkenn92, xerrkenn92, ykenn92, yerrkenn92, index=indxkenn92, snrcut=snrcut
    endif    

    lineratio, nfgsnodust, 'NII_6584', 'H_ALPHA', 'OII_3727', 'H_ALPHA', $
      xnfgs, xerrnfgs, ynfgs, yerrnfgs, index=indxnfgs, snrcut=snrcut, $
      xsyserr=nfgssyserr, ysyserr=nfgssyserr
    xbig = [xbig,xnfgs]
    ybig = [ybig,ynfgs]

    xtitle = 'log ([N II] \lambda6584/H\alpha)_{cor}'
    ytitle = 'log ([O II]/H\alpha)_{cor}'

    xrange = niiharangenoagn
    yrange = oiiharange

; HII regions

    good = where((hii.nii_6584_h_alpha gt -900.0) and (hii.oii_h_alpha gt -900.0))
    xregion = hii[good].nii_6584_h_alpha & xerrregion = hii[good].nii_6584_h_alpha_err
    yregion = hii[good].oii_h_alpha & yerrregion = hii[good].oii_h_alpha_err

    atlas1d_lineplot, x, y, xerr, yerr, plottype=3, postscript=postscript, $
      xtitle=xtitle, ytitle=ytitle, xrange=xrange, yrange=yrange, legendtype=0, $
      xstyle=11, ymargin=[4,3], /right, /top, $
      xregion=xregion, yregion=yregion, xerrregion=xerrregion, yerrregion=yerrregion, $
;     xkenn92=xkenn92, ykenn92=ykenn92, xerrkenn92=xerrkenn92, yerrkenn92=yerrkenn92, $
      xnfgs=xnfgs, ynfgs=ynfgs, xerrnfgs=xerrnfgs, yerrnfgs=yerrnfgs
    axis, /xaxis, xrange=interpol(xOHgood,xgood,!x.crange), xthick=postthick, xstyle=1, $
      charsize=singlecharsize, charthick=postthick
    xyouts, mean(!x.crange), !y.crange[1]+0.07*(!y.crange[1]-!y.crange[0]), textoidl('12 + log (O/H)'), $
      align=0.5, charsize=singlecharsize, charthick=postthick
    
; overlay the three metallicity regions of interest

    r12 = -1.3 & r23 = -0.6
    OH12 = interpol(xOHgood,xgood,r12) & OH23 = interpol(xOHgood,xgood,r23)
    Z12 = 10^(OH12 - Zsun_new) & Z23 = 10^(OH23 - Zsun_new)    
    
    splog, '[N II]/Ha = '+string(r12,format='(F4.1)')+' --> 12+log(O/H) = '+$
      string(OH12,format='(F4.2)')+', Z/Z_sun = '+string(Z12,format='(F4.2)')
    splog, '[N II]/Ha = '+string(r23,format='(F4.1)')+' --> 12+log(O/H) = '+$
      string(OH23,format='(F4.2)')+', Z/Z_sun = '+string(Z23,format='(F4.2)')
    
    oplot, r12*[1,1], [-1.25,0.75], line=1, thick=postthick
    oplot, r23*[1,1], [-1.25,0.75], line=1, thick=postthick
 
    xyouts, -1.75, 0.65, 'region 1', align=0.5, /data, charsize=charsize, charthick=postthick
    xyouts, -0.95, 0.65, 'region 2', align=0.5, /data, charsize=charsize, charthick=postthick
    xyouts, -0.25, 0.65, 'region 3', align=0.5, /data, charsize=charsize, charthick=postthick

    region1 = where(xbig lt r12,nregion1)
    splog, 'Region 1, 2, 3: '
    stats = im_stats(ybig[region1],sigrej=3.0,/verbose)

    regioniiha = where((xbig gt r12) and (xbig lt r23),nregioniiha)
    stats = im_stats(ybig[regioniiha],sigrej=3.0,/no_head,/verbose)

    region3 = where(xbig gt r23,nregion3)
    stats = im_stats(ybig[region3],sigrej=3.0,/no_head,/verbose)

;   if keyword_set(kewley_grids) then begin
;      plot_kewley_grids, plotnumber=11, model=3, labeltype=4, /noZgrid, $
;        /overplot, Umax=-2.5, Umin=-3.5, Zmax=2.9, Zmin=0.15, postscript=postscript
;      plot_kewley_grids, plotnumber=11, model=3, labeltype=2, /noUgrid, $
;        /overplot, Umax=-2.5, Umin=-3.5, Zmax=2.9, Zmin=0.15, postscript=postscript
;   endif

    im_openclose, postscript=postscript, /close    

; ------------------------------------------------------------
; [N II]/Ha vs [O II]_cor/Ha_cor
; ------------------------------------------------------------

    psname = 'sdss_niiha_vs_oiiha_cor'
    im_openclose, pspath+psname, postscript=postscript, /encapsulated

    lineratio, sdssnodust, 'NII_6584', 'H_ALPHA', 'OII_3727', 'H_ALPHA', $
      x, xerr, y, yerr, index=indx, nindex=nindx, snrcut=snrcut

    OHgood = where(sdssnodust[indx].zstrong_12oh_niiha gt -900)
    xgood = x[OHgood]
    xOHgood = sdssnodust[indx[OHgood]].zstrong_12oh_niiha
    
    xbig = x
    ybig = y

    xtitle = 'log ([N II] \lambda6584/H\alpha)_{cor}'
    ytitle = 'log ([O II]/H\alpha)_{cor}'

    xrange = niiharangenoagn
    yrange = oiiharange

    sdss_lineplot, x, y, xerr, yerr, plottype=3, postscript=postscript, $
      xtitle=xtitle, ytitle=ytitle, xrange=xrange, yrange=yrange, legendtype=0, $
      xstyle=11, ymargin=[4,3], /right, /top
    axis, /xaxis, xrange=interpol(xOHgood,xgood,!x.crange), xthick=postthick, xstyle=1, $
      charsize=singlecharsize, charthick=postthick
    xyouts, mean(!x.crange), !y.crange[1]+0.07*(!y.crange[1]-!y.crange[0]), textoidl('12 + log (O/H)'), $
      align=0.5, charsize=singlecharsize, charthick=postthick

; overlay the three metallicity regions of interest

    r12 = -1.3 & r23 = -0.6
    OH12 = interpol(xOHgood,xgood,r12) & OH23 = interpol(xOHgood,xgood,r23)
    Z12 = 10^(OH12 - Zsun_new) & Z23 = 10^(OH23 - Zsun_new)    
    
    splog, '[N II]/Ha = '+string(r12,format='(F4.1)')+' --> 12+log(O/H) = '+$
      string(OH12,format='(F4.2)')+', Z/Z_sun = '+string(Z12,format='(F4.2)')
    splog, '[N II]/Ha = '+string(r23,format='(F4.1)')+' --> 12+log(O/H) = '+$
      string(OH23,format='(F4.2)')+', Z/Z_sun = '+string(Z23,format='(F4.2)')
    
    oplot, r12*[1,1], [-1.25,0.75], line=1, thick=postthick
    oplot, r23*[1,1], [-1.25,0.75], line=1, thick=postthick
 
    xyouts, -1.75, 0.65, 'region 1', align=0.5, /data, charsize=1.5, charthick=postthick
    xyouts, -0.95, 0.65, 'region 2', align=0.5, /data, charsize=1.5, charthick=postthick
    xyouts, -0.25, 0.65, 'region 3', align=0.5, /data, charsize=1.5, charthick=postthick

    region1 = where(xbig lt r12,nregion1)
    splog, 'Region 1, 2, 3: '
    stats = im_stats(ybig[region1],sigrej=3.0,/verbose)

    regioniiha = where((xbig gt r12) and (xbig lt r23),nregioniiha)
    stats = im_stats(ybig[regioniiha],sigrej=3.0,/no_head,/verbose)

    region3 = where(xbig gt r23,nregion3)
    stats = im_stats(ybig[region3],sigrej=3.0,/no_head,/verbose)

    im_openclose, postscript=postscript, /close    

; ------------------------------------------------------------
; 3-panel - L(B) vs L(Ha)/L(IR) - Integrated
; ------------------------------------------------------------
    
    psname = 'LB_vs_LHa_LIR_multipanel'
    im_openclose, pspath+psname, postscript=postscript, /encapsulated

    indx = where((atlasdust.b_lum_obs gt -900) and (atlasnodust.ir_flux gt -900) and $
      (atlasnodust.ebv_hahb_err gt 0.0),nindx)

    x = atlasdust[indx].b_lum_obs
    xerr = atlasdust[indx].b_lum_obs_err
    xabs = atlasdust[indx].m_b_obs
    
    lir = atlasdust[indx].ir_flux ; [erg/s/cm2]
    lir_err = sqrt(atlasdust[indx].ir_flux_err^2 + (atlasdust[indx].ir_flux*0.15)^2)

    indxnfgs = where((nfgsdust.b_lum_obs gt -900) and (nfgsdust.ir_flux gt -900) and $
      (nfgsnodust.ebv_hahb_err gt 0.0),nindxnfgs)

    xnfgs = nfgsnodust[indxnfgs].b_lum_obs
    xerrnfgs = nfgsnodust[indxnfgs].b_lum_obs_err
    
    lirnfgs = nfgsdust[indxnfgs].ir_flux
    lirnfgs_err = sqrt(nfgsdust[indxnfgs].ir_flux_err^2 + (nfgsdust[indxnfgs].ir_flux*0.15)^2)
    
    xrange = LBrange
    yrange = haLIRrange

    xtitle = 'log L(B) [L(B)'+sunsymbol()+']'
    ytitle = 'log [L(H\alpha)/L(IR)]'

    lhalir = alog10(4.5D-44/7.9D-42)
    
; don't mess with the width and xmargin: these values were selected to
; center the plot on the page    

    pagemaker, nx=1, ny=3, yspace=0, xspace=0, xmargin=[1.1,0.2], $
      ymargin=[0.75,1.0], position=pos, /normal

; ##########################################################
; Panel 1: No reddening correction
; ##########################################################

    ha = atlasdust[indx].h_alpha[0]
    ha_err = atlasdust[indx].h_alpha[1]
    
    y = alog10(ha/lir)
    yerr = im_compute_error(ha,ha_err,lir,lir_err,/log)

    hanfgs = nfgsdust[indxnfgs].h_alpha[0]
    hanfgs_err = nfgsdust[indxnfgs].h_alpha[1]
    
    ynfgs = alog10(hanfgs/lirnfgs)
    yerrnfgs = im_compute_error(hanfgs,hanfgs_err,lirnfgs,lirnfgs_err,/log)

    stats = im_stats([y,ynfgs],/verbose)
    xstr = strtrim(string(stats.median,format='(F12.2)'),2)+'\pm'+strtrim(string(stats.sig68mean,format='(F12.2)'),2)

    atlas1d_lineplot, x, y, xerr, yerr, plottype=3, postscript=postscript, $
      xtitle='', ytitle='', xrange=xrange, yrange=yrange, legendtype=0, $
      charsize=charsize, xstyle=11, position=pos[*,0], xnfgs=xnfgs, ynfgs=ynfgs, $
      xerrnfgs=xerrnfgs, yerrnfgs=yerrnfgs, xtickname=replicate(' ',10), yminor=3
    djs_oplot, !x.crange, lhalir*[1,1], line=2, thick=postthick
    axis, /xaxis, xrange=interpol(xabs,x,!x.crange), xthick=postthick, xstyle=1, $
      charsize=charsize, charthick=postthick
    xyouts, mean(!x.crange), !y.crange[1]+0.14*(!y.crange[1]-!y.crange[0]), textoidl('M_{B} [mag]'), $
      align=0.5, charsize=charsize, charthick=postthick

    legend, '(a)', /left, /top, box=0, charsize=charsize, charthick=postthick

    label = ['A(H\alpha) = 0']
    legend, textoidl(xstr), /right, /top, box=0, charsize=charsize_3, charthick=postthick
    legend, textoidl(label), /left, /bottom, box=0, charsize=charsize_3, $
      charthick=postthick

; ############################################################
; Panel 2: A(Ha) = 1
; ############################################################

    Aha = 1.0
    Aha_err = 0.1
    
    ha = atlasdust[indx].h_alpha[0]*10^(0.4*Aha)
    ha_err = atlasdust[indx].h_alpha[1]*10^(0.4*Aha) + $
      atlasdust[indx].h_alpha[0]*alog(10.0)*10^(0.4*Aha)*Aha_err

    y = alog10(ha/lir)
    yerr = im_compute_error(ha,ha_err,lir,lir_err,/log)

    hanfgs = nfgsdust[indxnfgs].h_alpha[0]*10^(0.4*Aha)
    hanfgs_err = nfgsdust[indxnfgs].h_alpha[1]*10^(0.4*Aha) + $
      nfgsdust[indxnfgs].h_alpha[0]*alog(10.0)*10^(0.4*Aha)*Aha_err

    ynfgs = alog10(hanfgs/lirnfgs)
    yerrnfgs = im_compute_error(hanfgs,hanfgs_err,lirnfgs,lirnfgs_err,/log)

    stats = im_stats([y,ynfgs],/verbose,/no_head)
    xstr = strtrim(string(stats.median,format='(F12.2)'),2)+'\pm'+strtrim(string(stats.sig68mean,format='(F12.2)'),2)

    atlas1d_lineplot, x, y, xerr, yerr, plottype=3, postscript=postscript, $
      xtitle='', ytitle=ytitle, xrange=xrange, yrange=yrange, legendtype=0, $
      charsize=charsize, xtickname=replicate(' ',10), yminor=3, $
      position=pos[*,1], /noerase, xnfgs=xnfgs, ynfgs=ynfgs, xerrnfgs=xerrnfgs, yerrnfgs=yerrnfgs
    djs_oplot, !x.crange, lhalir*[1,1], line=2, thick=postthick

    legend, '(b)', /left, /top, box=0, charsize=charsize, charthick=postthick

    label = ['Mean A(H\alpha)']
    legend, textoidl(xstr), /right, /top, box=0, charsize=charsize_3, charthick=postthick
    legend, textoidl(label), /left, /bottom, box=0, charsize=charsize_3, $
      charthick=postthick
    
; ############################################################
; Panel 3: Good Reddening Correction
; ############################################################

    ha = atlasnodust[indx].h_alpha[0]
    ha_err = atlasnodust[indx].h_alpha[1]

    y = alog10(ha/lir)
    yerr = im_compute_error(ha,ha_err,lir,lir_err,/log)

    hanfgs = nfgsnodust[indxnfgs].h_alpha[0]
    hanfgs_err = nfgsnodust[indxnfgs].h_alpha[1]

    ynfgs = alog10(hanfgs/lirnfgs)
    yerrnfgs = im_compute_error(hanfgs,hanfgs_err,lirnfgs,lirnfgs_err,/log)

    stats = im_stats([y,ynfgs],/verbose,/no_head)
    xstr = strtrim(string(stats.median,format='(F12.2)'),2)+'\pm'+strtrim(string(stats.sig68mean,format='(F12.2)'),2)

    atlas1d_lineplot, x, y, xerr, yerr, plottype=3, postscript=postscript, $
      xtitle=xtitle, ytitle='', xrange=xrange, yrange=yrange, legendtype=0, $
      charsize=charsize, position=pos[*,2], /noerase, xnfgs=xnfgs, ynfgs=ynfgs, $
      xerrnfgs=xerrnfgs, yerrnfgs=yerrnfgs, yminor=3
    djs_oplot, !x.crange, lhalir*[1,1], line=2, thick=postthick

    legend, '(c)', /left, /top, box=0, charsize=charsize, charthick=postthick

    label = ['Individual A(H\alpha)']
    legend, textoidl(xstr), /right, /top, box=0, charsize=charsize_3, charthick=postthick
    legend, textoidl(label), /left, /bottom, box=0, charsize=charsize_3, $
      charthick=postthick
    
    im_openclose, postscript=postscript, /close

; ------------------------------------------------------------
; 3-panel [O II] plot - E(Hb-Ha) vs [O II]/Ha
; ------------------------------------------------------------

    psname = 'ehbha_vs_oiiha_multipanel'
    im_openclose, pspath+psname, postscript=postscript, /encapsulated

    indx = where((atlasdust.oii_3727[0]/atlasdust.oii_3727[1] gt snrcut) and $
      (atlasnodust.oii_3727[0]/atlasnodust.oii_3727[1] gt snrcut) and $
      (atlasdust.h_alpha[0]/atlasdust.h_alpha[1] gt snrcut) and $
      (atlasnodust.h_alpha[0]/atlasnodust.h_alpha[1] gt snrcut) and $
      (atlasnodust.ehbha_err gt 0.0),nindx)

    x = atlasnodust[indx].ehbha
    xerr = atlasnodust[indx].ehbha_err
    xebv = x/(k_lambda(4861,/odonnell)-k_lambda(6563,/odonnell))

    indxnfgs = where((nfgsdust.oii_3727[0]/nfgsdust.oii_3727[1] gt snrcut) and $
      (nfgsnodust.oii_3727[0]/nfgsnodust.oii_3727[1] gt snrcut) and $
      (nfgsdust.h_alpha[0]/nfgsdust.h_alpha[1] gt snrcut) and $
      (nfgsnodust.h_alpha[0]/nfgsnodust.h_alpha[1] gt snrcut) and $
      (nfgsnodust.ehbha_err gt 0.0),nindxnfgs)

    xnfgs = nfgsnodust[indxnfgs].ehbha
    xerrnfgs = nfgsnodust[indxnfgs].ehbha_err

    xtitle = 'E(H\beta-H\alpha) [mag]'
    ytitle = 'log ([O II]/H\alpha)'

    xrange = ehbharange
    yrange = oiihacorrange
    
; don't mess with the width and xmargin: these values were selected to
; center the plot on the page    
    
    pagemaker, nx=1, ny=3, yspace=0, xspace=0, xmargin=[2.0,1.0], $
      ymargin=[1.1,1.1], width=5.5, position=pos, /normal

; ##########################################################
; Panel 1: [O II] observed, Ha observed
; ##########################################################

    y1 = atlasdust[indx].oii_3727[0]
    y1err = atlasdust[indx].oii_3727[1]

    y2 = atlasdust[indx].h_alpha[0]
    y2err = atlasdust[indx].h_alpha[1]

    y = alog10(y1/y2)
    yerr = im_compute_error(y1,y1err,y2,y2err,/log)

    xbig = x
    ybig = y
    
    y1nfgs = nfgsdust[indxnfgs].oii_3727[0]
    y1errnfgs = nfgsdust[indxnfgs].oii_3727[1]

    y2nfgs = nfgsdust[indxnfgs].h_alpha[0]
    y2errnfgs = nfgsdust[indxnfgs].h_alpha[1]

    ynfgs = alog10(y1nfgs/y2nfgs)
    yerrnfgs = im_compute_error(y1nfgs,y1errnfgs,y2nfgs,y2errnfgs,/log)

    ybig = [ybig,ynfgs]
    xbig = [xbig,xnfgs]

    stats = im_stats(ybig,/verbose)
    rcor = r_correlate(xbig,ybig,zd=zd,probd=probd)
    splog, 'Spearman rank for E(Hb-Ha) vs [O II]/Ha observed: ', rcor, probd

    atlas1d_lineplot, x, y, xerr, yerr, plottype=3, postscript=postscript, $
      ytitle='', xrange=xrange, yrange=yrange, legendtype=0, xminor=3, yminor=3, $
      xstyle=11, /right, /top, xnfgs=xnfgs, ynfgs=ynfgs, xerrnfgs=xerrnfgs, $
      yerrnfgs=yerrnfgs, position=pos[*,0], charsize=charsize, xtickname=replicate(' ',10)
    axis, /xaxis, xrange=interpol(xebv,x,!x.crange), xthick=postthick, xstyle=1, $
      charsize=charsize, charthick=postthick
    xyouts, mean(!x.crange), !y.crange[1]+0.14*(!y.crange[1]-!y.crange[0]), textoidl('E(B-V) [mag]'), $
      align=0.5, charsize=charsize, charthick=postthick

    legend, '(a)', /left, /top, box=0, charsize=charsize, charthick=postthick

    label = ['A([O II]) = 0','A(H\alpha) = 0']
    legend, textoidl(label), /left, /bottom, box=0, charsize=charsize_3, charthick=postthick

; overplot reddening relations

    ratio = 0.0

    x_ebv = findgen((1.1-(-0.5))/0.1)*0.1+(-0.5)
    yodonnell =  -0.4*x_ebv*(k_lambda(3727.0,/odonnell)-k_lambda(6563.0,/odonnell)) + ratio
    ysmc = -0.4*x_ebv*(k_lambda(3727,/smc)-k_lambda(6563,/smc)) + ratio

    djs_oplot, x_ebv, yodonnell, line=0, thick=postthick
    djs_oplot, x_ebv, ysmc, line=2, thick=postthick

    label = ['Milky Way','SMC Bar']
    linestyle = [0,2]
    legend, label, /right, /top, box=0, linestyle=linestyle, $
      charsize=charsize_3, charthick=postthick, thick=postthick, clear=postscript

; ##########################################################
; Panel 2: [O II] observed, Ha corrected
; ##########################################################

    y1 = atlasdust[indx].oii_3727[0]
    y1err = atlasdust[indx].oii_3727[1]

    y2 = atlasnodust[indx].h_alpha[0]
    y2err = atlasnodust[indx].h_alpha[1]

    y = alog10(y1/y2)
    yerr = im_compute_error(y1,y1err,y2,y2err,/log)

    xbig = x
    ybig = y
    
    y1nfgs = nfgsdust[indxnfgs].oii_3727[0]
    y1errnfgs = nfgsdust[indxnfgs].oii_3727[1]

    y2nfgs = nfgsnodust[indxnfgs].h_alpha[0]
    y2errnfgs = nfgsnodust[indxnfgs].h_alpha[1]

    ynfgs = alog10(y1nfgs/y2nfgs)
    yerrnfgs = im_compute_error(y1nfgs,y1errnfgs,y2nfgs,y2errnfgs,/log)

    ybig = [ybig,ynfgs]
    xbig = [xbig,xnfgs]

    stats = im_stats(ybig,/verbose)

    atlas1d_lineplot, x, y, xerr, yerr, plottype=3, postscript=postscript, $
      xrange=xrange, yrange=yrange, legendtype=0, /noerase, xtickname=replicate(' ',10), $
      /right, /top, xnfgs=xnfgs, ynfgs=ynfgs, xerrnfgs=xerrnfgs, xminor=3, yminor=3, $
      yerrnfgs=yerrnfgs, position=pos[*,1], charsize=charsize, ytitle=ytitle
    
    legend, '(b)', /left, /top, box=0, charsize=charsize, charthick=postthick

    label = ['A([O II]) = 0','Individual A(H\alpha)']
    legend, textoidl(label), /left, /bottom, box=0, charsize=charsize_3, charthick=postthick

; overplot a histogram

    binsize = 0.1
    plothist, ybig, bin=binsize, xbin, ybin, /noplot, /halfbin

    yhistrange = minmax(ybin)*[1.0,1.05]

; ##########################################################
; Panel 3: [O II] corrected, Ha corrected
; ##########################################################

    y1 = atlasnodust[indx].oii_3727[0]
    y1err = atlasnodust[indx].oii_3727[1]

    y2 = atlasnodust[indx].h_alpha[0]
    y2err = atlasnodust[indx].h_alpha[1]

    y = alog10(y1/y2)
    yerr = im_compute_error(y1,y1err,y2,y2err,/log)

    xbig = x
    ybig = y
    
    y1nfgs = nfgsnodust[indxnfgs].oii_3727[0]
    y1errnfgs = nfgsnodust[indxnfgs].oii_3727[1]

    y2nfgs = nfgsnodust[indxnfgs].h_alpha[0]
    y2errnfgs = nfgsnodust[indxnfgs].h_alpha[1]

    ynfgs = alog10(y1nfgs/y2nfgs)
    yerrnfgs = im_compute_error(y1nfgs,y1errnfgs,y2nfgs,y2errnfgs,/log)

    ybig = [ybig,ynfgs]
    xbig = [xbig,xnfgs]

    stats = im_stats(ybig,/verbose)
    rcor = r_correlate(xbig,ybig,zd=zd,probd=probd)
    splog, 'Spearman rank for E(Hb-Ha) vs [O II]/Ha corrected: ', rcor, probd

    atlas1d_lineplot, x, y, xerr, yerr, plottype=3, postscript=postscript, $
      xtitle=xtitle, xrange=xrange, yrange=yrange, legendtype=0, /noerase, $
      xstyle=3, /right, /top, xnfgs=xnfgs, ynfgs=ynfgs, xerrnfgs=xerrnfgs, $
      yerrnfgs=yerrnfgs, position=pos[*,2], charsize=charsize, xminor=3, yminor=3
    
    legend, '(c)', /left, /top, box=0, charsize=charsize, charthick=postthick

    label = ['Individual A([O II])','Individual A(H\alpha)']
    legend, textoidl(label), /left, /bottom, box=0, charsize=charsize_3, charthick=postthick

    im_openclose, postscript=postscript, /close

; ------------------------------------------------------------
; 3-panel - E(Hb-Ha) vs [O II]/Ha
; ------------------------------------------------------------

    psname = 'sdss_ehbha_vs_oiiha_multipanel'
    im_openclose, pspath+psname, postscript=postscript, /encapsulated

    indx = where((sdssdust.oii_3727[0]/sdssdust.oii_3727[1] gt snrcut) and $
      (sdssnodust.oii_3727[0]/sdssnodust.oii_3727[1] gt snrcut) and $
      (sdssdust.h_alpha[0]/sdssdust.h_alpha[1] gt snrcut) and $
      (sdssnodust.h_alpha[0]/sdssnodust.h_alpha[1] gt snrcut) and $
      (sdssnodust.ehbha_err gt 0.0),nindx)

    x = sdssnodust[indx].ehbha
    xerr = sdssnodust[indx].ehbha_err
    xebv = x/(k_lambda(4861,/odonnell)-k_lambda(6563,/odonnell))

    xtitle = 'E(H\beta-H\alpha) [mag]'
    ytitle = 'log ([O II]/H\alpha)'

    xrange = ehbharange
    yrange = oiihacorrange
    
; don't mess with the width and xmargin: these values were selected to
; center the plot on the page    
    
    pagemaker, nx=1, ny=3, yspace=0, xspace=0, xmargin=[2.0,1.0], $
      ymargin=[1.1,1.1], width=5.5, position=pos, /normal

; ##########################################################
; Panel 1: [O II] observed, Ha observed
; ##########################################################

    y1 = sdssdust[indx].oii_3727[0]
    y1err = sdssdust[indx].oii_3727[1]

    y2 = sdssdust[indx].h_alpha[0]
    y2err = sdssdust[indx].h_alpha[1]

    y = alog10(y1/y2)
    yerr = im_compute_error(y1,y1err,y2,y2err,/log)

    xbig = x
    ybig = y
    
    stats = im_stats(ybig,/verbose)
    rcor = r_correlate(xbig,ybig,zd=zd,probd=probd)
    splog, 'Spearman rank for E(Hb-Ha) vs [O II]/Ha observed: ', rcor, probd

    sdss_lineplot, x, y, xerr, yerr, plottype=3, postscript=postscript, $
      ytitle='', xrange=xrange, yrange=yrange, legendtype=0, xminor=3, yminor=3, $
      xstyle=11, /right, /top, position=pos[*,0], charsize=charsize, xtickname=replicate(' ',10)
    axis, /xaxis, xrange=interpol(xebv,x,!x.crange), xthick=postthick, xstyle=1, $
      charsize=charsize, charthick=postthick
    xyouts, mean(!x.crange), !y.crange[1]+0.14*(!y.crange[1]-!y.crange[0]), textoidl('E(B-V) [mag]'), $
      align=0.5, charsize=charsize, charthick=postthick

    legend, '(a)', /left, /top, box=0, charsize=charsize, charthick=postthick

    label = ['A([O II]) = 0','A(H\alpha) = 0']
    legend, textoidl(label), /left, /bottom, box=0, charsize=charsize_3, charthick=postthick

; overplot reddening relations

    ratio = 0.0

    x_ebv = findgen((1.1-(-0.5))/0.1)*0.1+(-0.5)
    yodonnell =  -0.4*x_ebv*(k_lambda(3727.0,/odonnell)-k_lambda(6563.0,/odonnell)) + ratio
    ysmc = -0.4*x_ebv*(k_lambda(3727,/smc)-k_lambda(6563,/smc)) + ratio

    djs_oplot, x_ebv, yodonnell, line=0, thick=postthick
    djs_oplot, x_ebv, ysmc, line=2, thick=postthick

;   label = ['Milky Way','SMC Bar']
;   linestyle = [0,2]
;   legend, label, /right, /top, box=0, linestyle=linestyle, $
;     charsize=charsize_3, charthick=postthick, thick=postthick, clear=postscript

; ##########################################################
; Panel 2: [O II] observed, Ha corrected
; ##########################################################

    y1 = sdssdust[indx].oii_3727[0]
    y1err = sdssdust[indx].oii_3727[1]

    y2 = sdssnodust[indx].h_alpha[0]
    y2err = sdssnodust[indx].h_alpha[1]

    y = alog10(y1/y2)
    yerr = im_compute_error(y1,y1err,y2,y2err,/log)

    xbig = x
    ybig = y
    
    stats = im_stats(ybig,/verbose)

    sdss_lineplot, x, y, xerr, yerr, plottype=3, postscript=postscript, $
      xrange=xrange, yrange=yrange, legendtype=0, /noerase, xtickname=replicate(' ',10), $
      /right, /top, position=pos[*,1], charsize=charsize, ytitle=ytitle
    
    legend, '(b)', /left, /top, box=0, charsize=charsize, charthick=postthick

    label = ['A([O II]) = 0','Individual A(H\alpha)']
    legend, textoidl(label), /left, /bottom, box=0, charsize=charsize_3, charthick=postthick

; ##########################################################
; Panel 3: [O II] corrected, Ha corrected
; ##########################################################

    y1 = sdssnodust[indx].oii_3727[0]
    y1err = sdssnodust[indx].oii_3727[1]

    y2 = sdssnodust[indx].h_alpha[0]
    y2err = sdssnodust[indx].h_alpha[1]

    y = alog10(y1/y2)
    yerr = im_compute_error(y1,y1err,y2,y2err,/log)

    xbig = x
    ybig = y
    
    stats = im_stats(ybig,/verbose)
    rcor = r_correlate(xbig,ybig,zd=zd,probd=probd)
    splog, 'Spearman rank for E(Hb-Ha) vs [O II]/Ha corrected: ', rcor, probd

    sdss_lineplot, x, y, xerr, yerr, plottype=3, postscript=postscript, $
      xtitle=xtitle, xrange=xrange, yrange=yrange, legendtype=0, /noerase, $
      xstyle=3, /right, /top, position=pos[*,2], charsize=charsize, xminor=3, yminor=3
    
    legend, '(c)', /left, /top, box=0, charsize=1.3, charthick=postthick

    label = ['Individual A([O II])','Individual A(H\alpha)']
    legend, textoidl(label), /left, /bottom, box=0, charsize=charsize_3, charthick=postthick

    im_openclose, postscript=postscript, /close

; ------------------------------------------------------------
; 6-panel H-beta plot - L(B) vs Hb/Ha
; ------------------------------------------------------------
    
    psname = 'LB_vs_hbha_multipanel'
    im_openclose, pspath+psname, postscript=postscript, /encapsulated

    indx = where((atlasdust.b_lum_obs gt -900) and (atlasnodust.ebv_hahb_err gt 0.0) and $
      (atlasdust.h_beta_ew_uncor[0] ge ewcut),nindx)

    x = atlasdust[indx].b_lum_obs
    xabs = atlasdust[indx].m_b_obs
    xerr = atlasdust[indx].b_lum_obs_err
    
    indxnfgs = where((nfgsdust.b_lum_obs gt -900) and (nfgsnodust.ebv_hahb_err gt 0.0) and $
      (nfgsdust.h_beta_ew_uncor[0] ge ewcut),nindxnfgs)

    xnfgs = nfgsnodust[indxnfgs].b_lum_obs
    xerrnfgs = nfgsnodust[indxnfgs].b_lum_obs_err
    
    Aha = 1.0
    Aha_err = 0.1

    Ahb = Aha*k_lambda(4861,/odonnell)/k_lambda(6563,/odonnell)
    Ahb_err = 0.1*k_lambda(4861,/odonnell)/k_lambda(6563,/odonnell)

    xrange = LBrange
    yrange = [-2.4,0.8]

    xtitle = 'log L(B) [L(B)'+sunsymbol()+']'
    ytitle = 'log (H\beta/H\alpha)'

    pagemaker, nx=2, ny=3, yspace=0, xspace=0, xmargin=[1.25,0.35], $
      ymargin=[0.85,1.05], position=pos, /normal

; ##########################################################
; Panel 1: No absorption correction, A(Hb)=0, A(Ha)=0
; ##########################################################

    hb = atlasdust[indx].h_beta_uncor[0]
    hb_err = atlasdust[indx].h_beta_uncor[1]
    
    ha = atlasdust[indx].h_alpha[0]
    ha_err = atlasdust[indx].h_alpha[1]
    
    y = alog10(hb/ha)
    yerr = im_compute_error(hb,hb_err,ha,ha_err,/log)

    hbnfgs = nfgsdust[indxnfgs].h_beta_uncor[0]
    hbnfgs_err = nfgsdust[indxnfgs].h_beta_uncor[1]
    
    hanfgs = nfgsdust[indxnfgs].h_alpha[0]
    hanfgs_err = nfgsdust[indxnfgs].h_alpha[1]
    
    ynfgs = alog10(hbnfgs/hanfgs)
    yerrnfgs = im_compute_error(hbnfgs,hbnfgs_err,hanfgs,hanfgs_err,/log)

    stats = im_stats([y,ynfgs],/verbose)
    xstr = strtrim(string(stats.median,format='(F12.2)'),2)+'\pm'+strtrim(string(stats.sig68mean,format='(F12.2)'),2)

    atlas1d_lineplot, x, y, xerr, yerr, plottype=3, postscript=postscript, $
      xtitle='', ytitle=ytitle, xrange=xrange, yrange=yrange, legendtype=0, $
      charsize=charsize, xstyle=11, xtickname=replicate(' ',10), position=pos[*,0], $
      xnfgs=xnfgs, ynfgs=ynfgs, xerrnfgs=xerrnfgs, yerrnfgs=yerrnfgs
    axis, /xaxis, xrange=interpol(xabs,x,!x.crange), xthick=postthick, xstyle=1, $
      charsize=charsize, charthick=postthick
    xyouts, mean(!x.crange), !y.crange[1]+0.16*(!y.crange[1]-!y.crange[0]), textoidl('M_{B} [mag]'), $
      align=0.5, charsize=charsize, charthick=postthick
    djs_oplot, !x.crange, alog10(1.0/HaHb)*[1,1], line=2, thick=postthick

    legend, '(a)', /left, /top, box=0, charsize=charsize, charthick=postthick
    legend, textoidl(xstr), /right, /top, box=0, charsize=charsize_3, charthick=postthick

    label = ['No Absorption Correction','A(H\beta) = 0, A(H\alpha) = 0']
    legend, textoidl(label), /left, /bottom, box=0, charsize=charsize_3, $
      charthick=postthick

; ############################################################
; Panel 2: Absorption-Corrected, A(Hb)=0, A(Ha)=0
; ############################################################

    hb = atlasdust[indx].h_beta[0]
    hb_err = atlasdust[indx].h_beta[1]
    
    ha = atlasdust[indx].h_alpha[0]
    ha_err = atlasdust[indx].h_alpha[1]
    
    y = alog10(hb/ha)
    yerr = im_compute_error(hb,hb_err,ha,ha_err,/log)

    hbnfgs = nfgsdust[indxnfgs].h_beta[0]
    hbnfgs_err = nfgsdust[indxnfgs].h_beta[1]
    
    hanfgs = nfgsdust[indxnfgs].h_alpha[0]
    hanfgs_err = nfgsdust[indxnfgs].h_alpha[1]
    
    ynfgs = alog10(hbnfgs/hanfgs)
    yerrnfgs = im_compute_error(hbnfgs,hbnfgs_err,hanfgs,hanfgs_err,/log)

    stats = im_stats([y,ynfgs],/verbose)
    xstr = strtrim(string(stats.median,format='(F12.2)'),2)+'\pm'+strtrim(string(stats.sig68mean,format='(F12.2)'),2)

    atlas1d_lineplot, x, y, xerr, yerr, plottype=3, postscript=postscript, $
      xtitle='', ytitle='', xrange=xrange, yrange=yrange, legendtype=0, $
      charsize=charsize, xstyle=11, xtickname=replicate(' ',10), ytickname=replicate(' ',10), $
      position=pos[*,1], /noerase, xnfgs=xnfgs, ynfgs=ynfgs, xerrnfgs=xerrnfgs, yerrnfgs=yerrnfgs
    axis, /xaxis, xrange=interpol(xabs,x,!x.crange), xthick=postthick, xstyle=1, $
      charsize=charsize, charthick=postthick
    xyouts, mean(!x.crange), !y.crange[1]+0.16*(!y.crange[1]-!y.crange[0]), textoidl('M_{B} [mag]'), $
      align=0.5, charsize=charsize, charthick=postthick
    djs_oplot, !x.crange, alog10(1.0/HaHb)*[1,1], line=2, thick=postthick

    legend, '(b)', /left, /top, box=0, charsize=charsize, charthick=postthick

    label = ['Absorption-Corrected','A(H\beta) = 0, A(H\alpha) = 0']
    legend, textoidl(xstr), /right, /top, box=0, charsize=charsize_3, charthick=postthick
    legend, textoidl(label), /left, /bottom, box=0, charsize=charsize_3, $
      charthick=postthick

; ############################################################
; Panel 3: Absorption-Corrected, A(Hb)=0, Mean A(Ha)
; ############################################################

    hb = atlasdust[indx].h_beta[0]
    hb_err = atlasdust[indx].h_beta[1]
    
    ha = atlasdust[indx].h_alpha[0]*10^(0.4*Aha)
    ha_err = atlasdust[indx].h_alpha[1]*10^(0.4*Aha) + $
      atlasdust[indx].h_alpha[0]*alog(10.0)*10^(0.4*Aha)*Aha_err

    y = alog10(hb/ha)
    yerr = im_compute_error(hb,hb_err,ha,ha_err,/log)

    hbnfgs = nfgsdust[indxnfgs].h_beta[0]
    hbnfgs_err = nfgsdust[indxnfgs].h_beta[1]
    
    hanfgs = nfgsdust[indxnfgs].h_alpha[0]*10^(0.4*Aha)
    hanfgs_err = nfgsdust[indxnfgs].h_alpha[1]*10^(0.4*Aha) + $
      nfgsdust[indxnfgs].h_alpha[0]*alog(10.0)*10^(0.4*Aha)*Aha_err

    ynfgs = alog10(hbnfgs/hanfgs)
    yerrnfgs = im_compute_error(hbnfgs,hbnfgs_err,hanfgs,hanfgs_err,/log)

    stats = im_stats([y,ynfgs],/verbose)
    xstr = strtrim(string(stats.median,format='(F12.2)'),2)+'\pm'+strtrim(string(stats.sig68mean,format='(F12.2)'),2)
    
    atlas1d_lineplot, x, y, xerr, yerr, plottype=3, postscript=postscript, $
      xtitle='', ytitle=ytitle, xrange=xrange, yrange=yrange, legendtype=0, $
      charsize=charsize, position=pos[*,2], /noerase, xnfgs=xnfgs, ynfgs=ynfgs, $
      xerrnfgs=xerrnfgs, yerrnfgs=yerrnfgs, xtickname=replicate(' ',10)
    djs_oplot, !x.crange, alog10(1.0/HaHb)*[1,1], line=2, thick=postthick

    legend, '(c)', /left, /top, box=0, charsize=charsize, charthick=postthick

    label = ['Absorption-Corrected','A(H\beta) = 0, Mean A(H\alpha)']
    legend, textoidl(xstr), /right, /top, box=0, charsize=charsize_3, charthick=postthick
    legend, textoidl(label), /left, /bottom, box=0, charsize=charsize_3, $
      charthick=postthick
    
; ############################################################
; Panel 4: Absorption-Corrected, Mean A(Hb), Mean A(Ha)
; ############################################################

    hb = atlasdust[indx].h_beta[0]*10^(0.4*Ahb)
    hb_err = atlasdust[indx].h_beta[1]*10^(0.4*Ahb) + $
      atlasdust[indx].h_beta[0]*alog(10.0)*10^(0.4*Ahb)*Ahb_err

    ha = atlasdust[indx].h_alpha[0]*10^(0.4*Aha)
    ha_err = atlasdust[indx].h_alpha[1]*10^(0.4*Aha) + $
      atlasdust[indx].h_alpha[0]*alog(10.0)*10^(0.4*Aha)*Aha_err

    y = alog10(hb/ha)
    yerr = im_compute_error(hb,hb_err,ha,ha_err,/log)

    hbnfgs = nfgsdust[indxnfgs].h_beta[0]*10^(0.4*Ahb)
    hbnfgs_err = nfgsdust[indxnfgs].h_beta[1]*10^(0.4*Ahb) + $
      nfgsdust[indxnfgs].h_beta[0]*alog(10.0)*10^(0.4*Ahb)*Ahb_err

    hanfgs = nfgsdust[indxnfgs].h_alpha[0]*10^(0.4*Aha)
    hanfgs_err = nfgsdust[indxnfgs].h_alpha[1]*10^(0.4*Aha) + $
      nfgsdust[indxnfgs].h_alpha[0]*alog(10.0)*10^(0.4*Aha)*Aha_err

    ynfgs = alog10(hbnfgs/hanfgs)
    yerrnfgs = im_compute_error(hbnfgs,hbnfgs_err,hanfgs,hanfgs_err,/log)

    stats = im_stats([y,ynfgs],/verbose)
    xstr = strtrim(string(stats.median,format='(F12.2)'),2)+'\pm'+strtrim(string(stats.sig68mean,format='(F12.2)'),2)
        
    atlas1d_lineplot, x, y, xerr, yerr, plottype=3, postscript=postscript, $
      xtitle='', ytitle='', xrange=xrange, yrange=yrange, legendtype=0, $
      charsize=charsize, ytickname=replicate(' ',10), position=pos[*,3], /noerase, $
      xnfgs=xnfgs, ynfgs=ynfgs, xerrnfgs=xerrnfgs, yerrnfgs=yerrnfgs, xtickname=replicate(' ',10)
    djs_oplot, !x.crange, alog10(1.0/HaHb)*[1,1], line=2, thick=postthick

    legend, '(d)', /left, /top, box=0, charsize=charsize, charthick=postthick

    label = ['Absorption-Corrected','Mean A(H\beta), Mean A(H\alpha)']
    legend, textoidl(xstr), /right, /top, box=0, charsize=charsize_3, charthick=postthick
    legend, textoidl(label), /left, /bottom, box=0, charsize=charsize_3, $
      charthick=postthick
    
; ###############################################################
; Panel 5: Absorption-Corrected, A(Hb)=0, Individual A(Ha)
; ###############################################################

    hb = atlasdust[indx].h_beta[0]
    hb_err = atlasdust[indx].h_beta[1]

    ha = atlasnodust[indx].h_alpha[0]
    ha_err = atlasnodust[indx].h_alpha[1]

    y = alog10(hb/ha)
    yerr = im_compute_error(hb,hb_err,ha,ha_err,/log)

    hbnfgs = nfgsdust[indxnfgs].h_beta[0]
    hbnfgs_err = nfgsdust[indxnfgs].h_beta[1]
    
    hanfgs = nfgsnodust[indxnfgs].h_alpha[0]
    hanfgs_err = nfgsnodust[indxnfgs].h_alpha[1]

    ynfgs = alog10(hbnfgs/hanfgs)
    yerrnfgs = im_compute_error(hbnfgs,hbnfgs_err,hanfgs,hanfgs_err,/log)
       
    stats = im_stats([y,ynfgs],/verbose)
    xstr = strtrim(string(stats.median,format='(F12.2)'),2)+'\pm'+strtrim(string(stats.sig68mean,format='(F12.2)'),2)

    atlas1d_lineplot, x, y, xerr, yerr, plottype=3, postscript=postscript, $
      xtitle=xtitle, ytitle=ytitle, xrange=xrange, yrange=yrange, legendtype=0, $
      charsize=charsize, position=pos[*,4], /noerase, xnfgs=xnfgs, ynfgs=ynfgs, $
      xerrnfgs=xerrnfgs, yerrnfgs=yerrnfgs
    djs_oplot, !x.crange, alog10(1.0/HaHb)*[1,1], line=2, thick=postthick

    legend, '(e)', /left, /top, box=0, charsize=charsize, charthick=postthick

    label = ['Absorption-Corrected','A(H\beta) = 0, Individual A(H\alpha)']
    legend, textoidl(xstr), /right, /top, box=0, charsize=charsize_3, charthick=postthick
    legend, textoidl(label), /left, /bottom, box=0, charsize=charsize_3, $
      charthick=postthick

    
; ###############################################################
; Panel 6: Absorption-Corrected, Mean A(Hb), Individual A(Ha)
; ###############################################################

    hb = atlasdust[indx].h_beta[0]*10^(0.4*Ahb)
    hb_err = atlasdust[indx].h_beta[1]*10^(0.4*Ahb) + $
      atlasdust[indx].h_beta[0]*alog(10.0)*10^(0.4*Ahb)*Ahb_err

    ha = atlasnodust[indx].h_alpha[0]
    ha_err = atlasnodust[indx].h_alpha[1]

    y = alog10(hb/ha)
    yerr = im_compute_error(hb,hb_err,ha,ha_err,/log)

    hbnfgs = nfgsdust[indxnfgs].h_beta[0]*10^(0.4*Ahb)
    hbnfgs_err = nfgsdust[indxnfgs].h_beta[1]*10^(0.4*Ahb) + $
      nfgsdust[indxnfgs].h_beta[0]*alog(10.0)*10^(0.4*Ahb)*Ahb_err

    hanfgs = nfgsnodust[indxnfgs].h_alpha[0]
    hanfgs_err = nfgsnodust[indxnfgs].h_alpha[1]

    ynfgs = alog10(hbnfgs/hanfgs)
    yerrnfgs = im_compute_error(hbnfgs,hbnfgs_err,hanfgs,hanfgs_err,/log)
       
    stats = im_stats([y,ynfgs],/verbose)
    xstr = strtrim(string(stats.median,format='(F12.2)'),2)+'\pm'+strtrim(string(stats.sig68mean,format='(F12.2)'),2)

    atlas1d_lineplot, x, y, xerr, yerr, plottype=3, postscript=postscript, $
      xtitle=xtitle, ytitle='', xrange=xrange, yrange=yrange, legendtype=0, $
      charsize=charsize, ytickname=replicate(' ',10), position=pos[*,5], /noerase, $
      xnfgs=xnfgs, ynfgs=ynfgs, xerrnfgs=xerrnfgs, yerrnfgs=yerrnfgs
    djs_oplot, !x.crange, alog10(1.0/HaHb)*[1,1], line=2, thick=postthick

    legend, '(f)', /left, /top, box=0, charsize=charsize, charthick=postthick

    label = ['Absorption-Corrected','Mean A(H\beta), Individual A(H\alpha)']
    legend, textoidl(xstr), /right, /top, box=0, charsize=charsize_3, charthick=postthick
    legend, textoidl(label), /left, /bottom, box=0, charsize=charsize_3, $
      charthick=postthick

    im_openclose, postscript=postscript, /close

; ------------------------------------------------------------
; 6-panel H-beta plot - L(B) vs Hb/Ha
; ------------------------------------------------------------
    
    psname = 'sdss_LB_vs_hbha_multipanel'
    im_openclose, pspath+psname, postscript=postscript, /encapsulated

    indx = where((sdssdust.b_lum gt -900) and (sdssnodust.ebv_hahb_err gt 0.0) and $
      (sdssdust.h_beta_ew_uncor[0] ge ewcut),nindx)

    x = sdssdust[indx].b_lum
    xabs = sdssdust[indx].m_b
    xerr = sdssdust[indx].b_lum_err
    
    Aha = 1.0
    Aha_err = 0.1

    Ahb = Aha*k_lambda(4861,/odonnell)/k_lambda(6563,/odonnell)
    Ahb_err = 0.1*k_lambda(4861,/odonnell)/k_lambda(6563,/odonnell)

    xrange = LBrange
    yrange = [-2.4,0.8]

    xtitle = 'log L(B) [L(B)'+sunsymbol()+']'
    ytitle = 'log (H\beta/H\alpha)'

    pagemaker, nx=2, ny=3, yspace=0, xspace=0, xmargin=[1.25,0.35], $
      ymargin=[0.85,1.05], position=pos, /normal

; ##########################################################
; Panel 1: No absorption correction, A(Hb)=0, A(Ha)=0
; ##########################################################

    hb = sdssdust[indx].h_beta_uncor[0]
    hb_err = sdssdust[indx].h_beta_uncor[1]
    
    ha = sdssdust[indx].h_alpha[0]
    ha_err = sdssdust[indx].h_alpha[1]
    
    y = alog10(hb/ha)
    yerr = im_compute_error(hb,hb_err,ha,ha_err,/log)

    stats = im_stats(y,/verbose)
    xstr = strtrim(string(stats.median,format='(F12.2)'),2)+'\pm'+strtrim(string(stats.sig68mean,format='(F12.2)'),2)

    sdss_lineplot, x, y, xerr, yerr, plottype=3, postscript=postscript, $
      xtitle='', ytitle=ytitle, xrange=xrange, yrange=yrange, legendtype=0, $
      charsize=charsize, xstyle=11, xtickname=replicate(' ',10), position=pos[*,0]
    axis, /xaxis, xrange=interpol(xabs,x,!x.crange), xthick=postthick, xstyle=1, $
      charsize=charsize, charthick=postthick
    xyouts, mean(!x.crange), !y.crange[1]+0.16*(!y.crange[1]-!y.crange[0]), textoidl('M_{B} [mag]'), $
      align=0.5, charsize=charsize, charthick=postthick
    djs_oplot, !x.crange, alog10(1.0/HaHb)*[1,1], line=2, thick=postthick

    legend, '(a)', /left, /top, box=0, charsize=charsize, charthick=postthick
    legend, textoidl(xstr), /right, /top, box=0, charsize=charsize_3, charthick=postthick

    label = ['No Absorption Correction','A(H\beta) = 0, A(H\alpha) = 0']
    legend, textoidl(label), /left, /bottom, box=0, charsize=charsize_3, $
      charthick=postthick

; ############################################################
; Panel 2: Absorption-Corrected, A(Hb)=0, A(Ha)=0
; ############################################################

    hb = sdssdust[indx].h_beta[0]
    hb_err = sdssdust[indx].h_beta[1]
    
    ha = sdssdust[indx].h_alpha[0]
    ha_err = sdssdust[indx].h_alpha[1]
    
    y = alog10(hb/ha)
    yerr = im_compute_error(hb,hb_err,ha,ha_err,/log)

    stats = im_stats(y,/verbose)
    xstr = strtrim(string(stats.median,format='(F12.2)'),2)+'\pm'+strtrim(string(stats.sig68mean,format='(F12.2)'),2)

    sdss_lineplot, x, y, xerr, yerr, plottype=3, postscript=postscript, $
      xtitle='', ytitle='', xrange=xrange, yrange=yrange, legendtype=0, $
      charsize=charsize, xstyle=11, xtickname=replicate(' ',10), ytickname=replicate(' ',10), $
      position=pos[*,1], /noerase
    axis, /xaxis, xrange=interpol(xabs,x,!x.crange), xthick=postthick, xstyle=1, $
      charsize=charsize, charthick=postthick
    xyouts, mean(!x.crange), !y.crange[1]+0.16*(!y.crange[1]-!y.crange[0]), textoidl('M_{B} [mag]'), $
      align=0.5, charsize=charsize, charthick=postthick
    djs_oplot, !x.crange, alog10(1.0/HaHb)*[1,1], line=2, thick=postthick

    legend, '(b)', /left, /top, box=0, charsize=charsize, charthick=postthick
    legend, textoidl(xstr), /right, /top, box=0, charsize=charsize_3, charthick=postthick

    label = ['Absorption-Corrected','A(H\beta) = 0, A(H\alpha) = 0']
    legend, textoidl(label), /left, /bottom, box=0, charsize=charsize_3, $
      charthick=postthick
    
; ############################################################
; Panel 3: Absorption-Corrected, A(Hb)=0, Mean A(Ha)
; ############################################################

    hb = sdssdust[indx].h_beta[0]
    hb_err = sdssdust[indx].h_beta[1]
    
    ha = sdssdust[indx].h_alpha[0]*10^(0.4*Aha)
    ha_err = sdssdust[indx].h_alpha[1]*10^(0.4*Aha) + $
      sdssdust[indx].h_alpha[0]*alog(10.0)*10^(0.4*Aha)*Aha_err

    y = alog10(hb/ha)
    yerr = im_compute_error(hb,hb_err,ha,ha_err,/log)

    stats = im_stats(y,/verbose)
    xstr = strtrim(string(stats.median,format='(F12.2)'),2)+'\pm'+strtrim(string(stats.sig68mean,format='(F12.2)'),2)
    
    sdss_lineplot, x, y, xerr, yerr, plottype=3, postscript=postscript, $
      xtitle='', ytitle=ytitle, xrange=xrange, yrange=yrange, legendtype=0, $
      charsize=charsize, position=pos[*,2], /noerase, xtickname=replicate(' ',10)
    djs_oplot, !x.crange, alog10(1.0/HaHb)*[1,1], line=2, thick=postthick

    legend, '(c)', /left, /top, box=0, charsize=charsize, charthick=postthick
    legend, textoidl(xstr), /right, /top, box=0, charsize=charsize_3, charthick=postthick

    label = ['Absorption-Corrected','A(H\beta) = 0, Mean A(H\alpha)']
    legend, textoidl(label), /left, /bottom, box=0, charsize=charsize_3, $
      charthick=postthick
    
; ###############################################################
; Panel 4: Absorption-Corrected, A(Hb)=0, Individual A(Ha)
; ###############################################################

    hb = sdssdust[indx].h_beta[0]*10^(0.4*Ahb)
    hb_err = sdssdust[indx].h_beta[1]*10^(0.4*Ahb) + $
      sdssdust[indx].h_beta[0]*alog(10.0)*10^(0.4*Ahb)*Ahb_err

    ha = sdssdust[indx].h_alpha[0]*10^(0.4*Aha)
    ha_err = sdssdust[indx].h_alpha[1]*10^(0.4*Aha) + $
      sdssdust[indx].h_alpha[0]*alog(10.0)*10^(0.4*Aha)*Aha_err

    y = alog10(hb/ha)
    yerr = im_compute_error(hb,hb_err,ha,ha_err,/log)

    stats = im_stats(y,/verbose)
    xstr = strtrim(string(stats.median,format='(F12.2)'),2)+'\pm'+strtrim(string(stats.sig68mean,format='(F12.2)'),2)

    sdss_lineplot, x, y, xerr, yerr, plottype=3, postscript=postscript, $
      xtitle='', ytitle='', xrange=xrange, yrange=yrange, legendtype=0, $
      charsize=charsize, ytickname=replicate(' ',10), position=pos[*,3], /noerase
    djs_oplot, !x.crange, alog10(1.0/HaHb)*[1,1], line=2, thick=postthick

    legend, '(d)', /left, /top, box=0, charsize=charsize, charthick=postthick
    legend, textoidl(xstr), /right, /top, box=0, charsize=charsize_3, charthick=postthick

    label = ['Absorption-Corrected','Mean A(H\beta), Mean A(H\alpha)']
    legend, textoidl(label), /left, /bottom, box=0, charsize=charsize_3, $
      charthick=postthick

; ############################################################
; Panel 5: Absorption-Corrected, Mean A(Hb), Mean A(Ha)
; ############################################################

    hb = sdssdust[indx].h_beta[0]
    hb_err = sdssdust[indx].h_beta[1]

    ha = sdssnodust[indx].h_alpha[0]
    ha_err = sdssnodust[indx].h_alpha[1]

    y = alog10(hb/ha)
    yerr = im_compute_error(hb,hb_err,ha,ha_err,/log)

    stats = im_stats(y,/verbose)
    xstr = strtrim(string(stats.median,format='(F12.2)'),2)+'\pm'+strtrim(string(stats.sig68mean,format='(F12.2)'),2)

    sdss_lineplot, x, y, xerr, yerr, plottype=3, postscript=postscript, $
      xtitle=xtitle, ytitle=ytitle, xrange=xrange, yrange=yrange, legendtype=0, $
      charsize=charsize, position=pos[*,4], /noerase
    djs_oplot, !x.crange, alog10(1.0/HaHb)*[1,1], line=2, thick=postthick

    legend, '(e)', /left, /top, box=0, charsize=charsize, charthick=postthick
    legend, textoidl(xstr), /right, /top, box=0, charsize=charsize_3, charthick=postthick

    label = ['Absorption-Corrected','A(H\beta) = 0, Individual A(H\alpha)']
    legend, textoidl(label), /left, /bottom, box=0, charsize=charsize_3, $
      charthick=postthick
    
; ###############################################################
; Panel 6: Absorption-Corrected, Mean A(Hb), Individual A(Ha)
; ###############################################################

    hb = sdssdust[indx].h_beta[0]*10^(0.4*Ahb)
    hb_err = sdssdust[indx].h_beta[1]*10^(0.4*Ahb) + $
      sdssdust[indx].h_beta[0]*alog(10.0)*10^(0.4*Ahb)*Ahb_err

    ha = sdssnodust[indx].h_alpha[0]
    ha_err = sdssnodust[indx].h_alpha[1]

    y = alog10(hb/ha)
    yerr = im_compute_error(hb,hb_err,ha,ha_err,/log)

    stats = im_stats(y,/verbose)
    xstr = strtrim(string(stats.median,format='(F12.2)'),2)+'\pm'+strtrim(string(stats.sig68mean,format='(F12.2)'),2)
    
    sdss_lineplot, x, y, xerr, yerr, plottype=3, postscript=postscript, $
      xtitle=xtitle, ytitle='', xrange=xrange, yrange=yrange, legendtype=0, $
      charsize=charsize, ytickname=replicate(' ',10), position=pos[*,5], /noerase
    djs_oplot, !x.crange, alog10(1.0/HaHb)*[1,1], line=2, thick=postthick

    legend, '(f)', /left, /top, box=0, charsize=charsize, charthick=postthick
    legend, textoidl(xstr), /right, /top, box=0, charsize=charsize_3, charthick=postthick

    label = ['Absorption-Corrected','Mean A(H\beta), Individual A(H\alpha)']
    legend, textoidl(label), /left, /bottom, box=0, charsize=charsize_3, $
      charthick=postthick

    im_openclose, postscript=postscript, /close

; ------------------------------------------------------------
; 3-panel H-alpha plot - A(Ha) vs L(Ha)/L(IR)
; ------------------------------------------------------------
    
    psname = 'AHa_vs_LHa_LIR_multipanel'
    im_openclose, pspath+psname, postscript=postscript

    indx = where((atlasdust.h_alpha[0]/atlasdust.h_alpha[1] gt snrcut) and $   
      (atlasnodust.h_alpha[0]/atlasnodust.h_alpha[1] gt snrcut) and $
      (atlasnodust.ir_flux gt -900) and (atlasnodust.ebv_hahb_err gt 0.0),nindx)

    kha = k_lambda(6563,/odonnel)    

    x = atlasnodust[indx].ebv_hahb*kha
    xerr = atlasnodust[indx].ebv_hahb_err*kha

    lir = atlasdust[indx].ir_flux ; [erg/s/cm2]
    lir_err = sqrt(atlasdust[indx].ir_flux_err^2 + (atlasdust[indx].ir_flux*0.15)^2)

    if keyword_set(nfgs) then begin

       indxnfgs = where((nfgsdust.h_alpha[0]/nfgsdust.h_alpha[1] gt snrcut) and $      
         (nfgsnodust.h_alpha[0]/nfgsnodust.h_alpha[1] gt snrcut) and $
         (nfgsdust.ir_flux gt -900) and (nfgsnodust.ebv_hahb_err gt 0.0),nindx)

       xnfgs = nfgsnodust[indxnfgs].ebv_hahb*kha
       xerrnfgs = nfgsnodust[indxnfgs].ebv_hahb_err*kha
       
       lirnfgs = nfgsdust[indxnfgs].ir_flux
       lirnfgs_err = sqrt(nfgsdust[indxnfgs].ir_flux_err^2 + (nfgsdust[indxnfgs].ir_flux*0.15)^2)

    endif
    
    xrange = AHarange
    yrange = haLIRrange

    xtitle = 'A(H\alpha)'
    ytitle = 'log [L(H\alpha)/L(IR)]'

    lhalir = alog10(4.5D-44/7.9D-42)
    
    pagemaker, nx=1, ny=3, yspace=0, xspace=0, xmargin=[1.1,0.2], $
      ymargin=[0.75,1.0], position=pos, /normal

; ##########################################################
; Panel 1: No reddening correction
; ##########################################################

    ha = atlasdust[indx].h_alpha[0]
    ha_err = atlasdust[indx].h_alpha[1]
    
    y = alog10(ha/lir)
    yerr = im_compute_error(ha,ha_err,lir,lir_err,/log)

    if keyword_set(nfgs) then begin

       hanfgs = nfgsdust[indxnfgs].h_alpha[0]
       hanfgs_err = nfgsdust[indxnfgs].h_alpha[1]
       
       ynfgs = alog10(hanfgs/lirnfgs)
       yerrnfgs = im_compute_error(hanfgs,hanfgs_err,lirnfgs,lirnfgs_err,/log)
       
    endif

    stats = im_stats(y,/verbose)

    atlas1d_lineplot, x, y, xerr, yerr, plottype=3, postscript=postscript, $
      xtitle='', ytitle='', xrange=xrange, yrange=yrange, legendtype=0, $
      charsize=charsize, position=pos[*,0], xnfgs=xnfgs, ynfgs=ynfgs, $
      xerrnfgs=xerrnfgs, yerrnfgs=yerrnfgs, xtickname=replicate(' ',10), yminor=3
    djs_oplot, !x.crange, lhalir*[1,1], line=2, thick=postthick

    legend, '(a)', /left, /top, box=0, charsize=charsize, charthick=postthick

    label = ['A(H\alpha) = 0']
    legend, textoidl(label), /left, /bottom, box=0, charsize=charsize_3, $
      charthick=postthick

; ############################################################
; Panel 2: A(Ha) = 1
; ############################################################

    Aha = 1.0
    Aha_err = 0.1
    
    ha = atlasdust[indx].h_alpha[0]*10^(0.4*Aha)
    ha_err = atlasdust[indx].h_alpha[1]*10^(0.4*Aha) + $
      atlasdust[indx].h_alpha[0]*alog(10.0)*10^(0.4*Aha)*Aha_err

    y = alog10(ha/lir)
    yerr = im_compute_error(ha,ha_err,lir,lir_err,/log)

    if keyword_set(nfgs) then begin

       hanfgs = nfgsdust[indxnfgs].h_alpha[0]*10^(0.4*Aha)
       hanfgs_err = nfgsdust[indxnfgs].h_alpha[1]*10^(0.4*Aha) + $
         nfgsdust[indxnfgs].h_alpha[0]*alog(10.0)*10^(0.4*Aha)*Aha_err

       ynfgs = alog10(hanfgs/lirnfgs)
       yerrnfgs = im_compute_error(hanfgs,hanfgs_err,lirnfgs,lirnfgs_err,/log)
       
    endif

    stats = im_stats(y,/verbose,/no_head)

    atlas1d_lineplot, x, y, xerr, yerr, plottype=3, postscript=postscript, $
      xtitle='', ytitle=ytitle, xrange=xrange, yrange=yrange, legendtype=0, $
      charsize=charsize, xtickname=replicate(' ',10), yminor=3, $
      position=pos[*,1], /noerase, xnfgs=xnfgs, ynfgs=ynfgs, xerrnfgs=xerrnfgs, yerrnfgs=yerrnfgs
    djs_oplot, !x.crange, lhalir*[1,1], line=2, thick=postthick

    legend, '(b)', /left, /top, box=0, charsize=charsize, charthick=postthick

    label = ['Mean A(H\alpha)']
    legend, textoidl(label), /left, /bottom, box=0, charsize=charsize_3, $
      charthick=postthick
    
; ############################################################
; Panel 3: Good Reddening Correction
; ############################################################

    ha = atlasnodust[indx].h_alpha[0]
    ha_err = atlasnodust[indx].h_alpha[1]

    y = alog10(ha/lir)
    yerr = im_compute_error(ha,ha_err,lir,lir_err,/log)

    if keyword_set(nfgs) then begin

       hanfgs = nfgsnodust[indxnfgs].h_alpha[0]
       hanfgs_err = nfgsnodust[indxnfgs].h_alpha[1]

       ynfgs = alog10(hanfgs/lirnfgs)
       yerrnfgs = im_compute_error(hanfgs,hanfgs_err,lirnfgs,lirnfgs_err,/log)
       
    endif

    stats = im_stats(y,/verbose,/no_head)

    atlas1d_lineplot, x, y, xerr, yerr, plottype=3, postscript=postscript, $
      xtitle=xtitle, ytitle='', xrange=xrange, yrange=yrange, legendtype=0, $
      charsize=charsize, position=pos[*,2], /noerase, xnfgs=xnfgs, ynfgs=ynfgs, $
      xerrnfgs=xerrnfgs, yerrnfgs=yerrnfgs, yminor=3
    djs_oplot, !x.crange, lhalir*[1,1], line=2, thick=postthick

    legend, '(c)', /left, /top, box=0, charsize=charsize, charthick=postthick

    label = ['Individual A(H\alpha)']
    legend, textoidl(label), /left, /bottom, box=0, charsize=charsize_3, $
      charthick=postthick

    im_openclose, postscript=postscript, /close

; ------------------------------------------------------------
; A(Ha) vs L(Ha)/L(IR)
; ------------------------------------------------------------
    
    psname = 'AHa_vs_LHa_LIR'
    im_openclose, pspath+psname, postscript=postscript, square=0, ysize=4.9

    indx = where((atlasnodust.h_alpha[0]/atlasnodust.h_alpha[1] gt snrcut) and $   
      (atlasnodust.ebv_hahb_err gt 0.0) and (atlasnodust.ir_flux gt -900),nindx)

    kha = k_lambda(6563,/odonnel)    

    x = atlasnodust[indx].ebv_hahb*kha
    xerr = atlasnodust[indx].ebv_hahb_err*kha

    lir = atlasdust[indx].ir_flux         ; [erg/s/cm2]
    lir_err = sqrt(atlasdust[indx].ir_flux_err^2 + (atlasdust[indx].ir_flux*0.15)^2)

    ha = atlasnodust[indx].h_alpha[0]
    ha_err = atlasnodust[indx].h_alpha[1]

    y = alog10(ha/lir)
    yerr = im_compute_error(ha,ha_err,lir,lir_err,/log)

    ybig = y
    xbig = x
    
    if keyword_set(nfgs) then begin

       indxnfgs = where((nfgsnodust.h_alpha[0]/nfgsnodust.h_alpha[1] gt snrcut) and $
         (nfgsnodust.ebv_hahb_err gt 0.0) and (nfgsdust.ir_flux gt -900),nindxnfgs)

       xnfgs = nfgsnodust[indxnfgs].ebv_hahb*kha
       xerrnfgs = nfgsnodust[indxnfgs].ebv_hahb_err*kha
       
       lirnfgs = nfgsdust[indxnfgs].ir_flux ; [erg/s/cm2]
       lirnfgs_err = sqrt(nfgsdust[indxnfgs].ir_flux_err^2 + (nfgsdust[indxnfgs].ir_flux*0.15)^2)

       hanfgs = nfgsnodust[indxnfgs].h_alpha[0]
       hanfgs_err = nfgsnodust[indxnfgs].h_alpha[1]

       ynfgs = alog10(hanfgs/lirnfgs)
       yerrnfgs = im_compute_error(hanfgs,hanfgs_err,lirnfgs,lirnfgs_err,/log)

       xbig = [xbig,ynfgs]
       ybig = [ybig,ynfgs]

    endif
    
    xrange = AHarange
    yrange = haLIRrange

    xtitle = 'A(H\alpha)'
    ytitle = 'log [L(H\alpha)/L(IR)]'

    lhalir = alog10(4.5D-44/7.9D-42)
    
    pagemaker, nx=1, ny=1, ypage=4.9, xspace=0, xmargin=[1.1,0.2], $
      height=3.08, ymargin=[0.75,1.0], position=pos, /normal

    stats = im_stats(y,/verbose)
    rcor = r_correlate(xbig,ybig,zd=zd,probd=probd)
    splog, 'Spearman rank test: ', rcor, probd

    atlas1d_lineplot, x, y, xerr, yerr, plottype=3, postscript=postscript, $
      xtitle=xtitle, ytitle=ytitle, xrange=xrange, yrange=yrange, legendtype=0, $
      position=pos, xnfgs=xnfgs, ynfgs=ynfgs, $
      xerrnfgs=xerrnfgs, yerrnfgs=yerrnfgs, yminor=3

    label = ['Individual A(H\alpha)']
    legend, textoidl(label), /left, /bottom, box=0, charsize=charsize_3, $
      charthick=postthick

    im_openclose, postscript=postscript, /close

; ------------------------------------------------------------
; 3-panel H-alpha plot - L(FIR) vs L(Ha)/L(FIR)
; ------------------------------------------------------------
    
    psname = 'LFIR_vs_LHa_LFIR_multipanel'
    im_openclose, pspath+psname, postscript=postscript

    indx = where((atlasdust.h_alpha[0]/atlasdust.h_alpha[1] gt snrcut) and $   
      (atlasdust.fir_lum gt -900) and (atlasnodust.h_alpha[0]/atlasnodust.h_alpha[1] gt snrcut) and $
      (atlasnodust.fir_flux gt -900) and (atlasnodust.ebv_hahb_err gt 0.0),nindx)

    x = atlasdust[indx].fir_lum
    xerr = atlasdust[indx].fir_lum_err
    
    lfir = atlasdust[indx].fir_flux ; [erg/s/cm2]
    lfir_err = atlasdust[indx].fir_flux_err

    if keyword_set(nfgs) then begin

       indxnfgs = where((nfgsdust.h_alpha[0]/nfgsdust.h_alpha[1] gt snrcut) and $      
         (nfgsdust.fir_lum gt -900) and (nfgsnodust.h_alpha[0]/nfgsnodust.h_alpha[1] gt snrcut) and $
         (nfgsdust.fir_flux gt -900) and (nfgsnodust.ebv_hahb_err gt 0.0),nindx)

       xnfgs = nfgsnodust[indxnfgs].fir_lum
       xerrnfgs = nfgsnodust[indxnfgs].fir_lum_err
       
       lfirnfgs = nfgsdust[indxnfgs].fir_flux
       lfirnfgs_err = nfgsdust[indxnfgs].fir_flux_err

    endif
    
    xrange = LIRrange
    yrange = haLIRrange

    xtitle = 'log [L(FIR)/L'+sunsymbol()+']'
    ytitle = 'log [L(H\alpha)/L(FIR)]'

    pagemaker, nx=1, ny=3, yspace=0, xspace=0, xmargin=[1.1,0.2], $
      ymargin=[0.75,1.0], position=pos, /normal

; ##########################################################
; Panel 1: No reddening correction
; ##########################################################

    ha = atlasdust[indx].h_alpha[0]
    ha_err = atlasdust[indx].h_alpha[1]
    
    y = alog10(ha/lfir)
    yerr = im_compute_error(ha,ha_err,lfir,lfir_err,/log)

    if keyword_set(nfgs) then begin

       hanfgs = nfgsdust[indxnfgs].h_alpha[0]
       hanfgs_err = nfgsdust[indxnfgs].h_alpha[1]
       
       ynfgs = alog10(hanfgs/lfirnfgs)
       yerrnfgs = im_compute_error(hanfgs,hanfgs_err,lfirnfgs,lfirnfgs_err,/log)
       
    endif

    stats = im_stats(y,/verbose)

    atlas1d_lineplot, x, y, xerr, yerr, plottype=3, postscript=postscript, $
      xtitle='', ytitle='', xrange=xrange, yrange=yrange, legendtype=0, $
      charsize=charsize, position=pos[*,0], xnfgs=xnfgs, ynfgs=ynfgs, $
      xerrnfgs=xerrnfgs, yerrnfgs=yerrnfgs, xtickname=replicate(' ',10), yminor=3

    legend, '(a)', /left, /top, box=0, charsize=charsize, charthick=postthick

    label = ['A(H\alpha) = 0']
    legend, textoidl(label), /left, /bottom, box=0, charsize=charsize_3, $
      charthick=postthick

; ############################################################
; Panel 2: A(Ha) = 1
; ############################################################

    Aha = 1.0
    Aha_err = 0.1
    
    ha = atlasdust[indx].h_alpha[0]*10^(0.4*Aha)
    ha_err = atlasdust[indx].h_alpha[1]*10^(0.4*Aha) + $
      atlasdust[indx].h_alpha[0]*alog(10.0)*10^(0.4*Aha)*Aha_err

    y = alog10(ha/lfir)
    yerr = im_compute_error(ha,ha_err,lfir,lfir_err,/log)

    if keyword_set(nfgs) then begin

       hanfgs = nfgsdust[indxnfgs].h_alpha[0]*10^(0.4*Aha)
       hanfgs_err = nfgsdust[indxnfgs].h_alpha[1]*10^(0.4*Aha) + $
         nfgsdust[indxnfgs].h_alpha[0]*alog(10.0)*10^(0.4*Aha)*Aha_err

       ynfgs = alog10(hanfgs/lfirnfgs)
       yerrnfgs = im_compute_error(hanfgs,hanfgs_err,lfirnfgs,lfirnfgs_err,/log)
       
    endif

    stats = im_stats(y,/verbose,/no_head)

    atlas1d_lineplot, x, y, xerr, yerr, plottype=3, postscript=postscript, $
      xtitle='', ytitle=ytitle, xrange=xrange, yrange=yrange, legendtype=0, $
      charsize=charsize, xtickname=replicate(' ',10), yminor=3, $
      position=pos[*,1], /noerase, xnfgs=xnfgs, ynfgs=ynfgs, xerrnfgs=xerrnfgs, yerrnfgs=yerrnfgs

    legend, '(b)', /left, /top, box=0, charsize=charsize, charthick=postthick

    label = ['Mean A(H\alpha)']
    legend, textoidl(label), /left, /bottom, box=0, charsize=charsize_3, $
      charthick=postthick
    
; ############################################################
; Panel 3: Good Reddening Correction
; ############################################################

    ha = atlasnodust[indx].h_alpha[0]
    ha_err = atlasnodust[indx].h_alpha[1]

    y = alog10(ha/lfir)
    yerr = im_compute_error(ha,ha_err,lfir,lfir_err,/log)

    if keyword_set(nfgs) then begin

       hanfgs = nfgsnodust[indxnfgs].h_alpha[0]
       hanfgs_err = nfgsnodust[indxnfgs].h_alpha[1]

       ynfgs = alog10(hanfgs/lfirnfgs)
       yerrnfgs = im_compute_error(hanfgs,hanfgs_err,lfirnfgs,lfirnfgs_err,/log)
       
    endif

    stats = im_stats(y,/verbose,/no_head)

    atlas1d_lineplot, x, y, xerr, yerr, plottype=3, postscript=postscript, $
      xtitle=xtitle, ytitle='', xrange=xrange, yrange=yrange, legendtype=0, $
      charsize=charsize, position=pos[*,2], /noerase, xnfgs=xnfgs, ynfgs=ynfgs, $
      xerrnfgs=xerrnfgs, yerrnfgs=yerrnfgs, yminor=3

    legend, '(c)', /left, /top, box=0, charsize=charsize, charthick=postthick

    label = ['Individual A(H\alpha)']
    legend, textoidl(label), /left, /bottom, box=0, charsize=charsize_3, $
      charthick=postthick

    im_openclose, postscript=postscript, /close

;;; ------------------------------------------------------------
;;; 4-panel H-beta plot - L(B) vs Hb/Ha
;;; ------------------------------------------------------------
;;    
;;    psname = 'LB_vs_hbha_multipanel'
;;    im_openclose, pspath+psname, postscript=postscript
;;
;;    indx = where((atlasdust.h_beta_uncor[0]/atlasdust.h_beta_uncor[1] gt snrcut) and $      
;;      (atlasdust.b_lum_obs gt -900) and (atlasnodust.h_alpha[0]/atlasnodust.h_alpha[1] gt snrcut) and $
;;      (atlasnodust.ebv_hahb_err gt 0.0),nindx)
;;
;;    x = atlasdust[indx].b_lum_obs
;;    xabs = atlasdust[indx].m_b_obs
;;    xerr = atlasdust[indx].b_lum_obs_err
;;    
;;    ha = atlasnodust[indx].h_alpha[0]
;;    ha_err = atlasnodust[indx].h_alpha[1]
;;
;;    if keyword_set(nfgs) then begin
;;
;;       indxnfgs = where((nfgsdust.h_beta_uncor[0]/nfgsdust.h_beta_uncor[1] gt snrcut) and $      
;;         (nfgsdust.b_lum_obs gt -900) and (nfgsnodust.ebv_hahb_err gt 0.0) and $
;;         (nfgsnodust.h_alpha[0]/nfgsnodust.h_alpha[1] gt snrcut),nindx)
;;
;;       xnfgs = nfgsnodust[indxnfgs].b_lum_obs
;;       xerrnfgs = nfgsnodust[indxnfgs].b_lum_obs_err
;;       
;;       hanfgs = nfgsnodust[indxnfgs].h_alpha[0]
;;       hanfgs_err = nfgsnodust[indxnfgs].h_alpha[1]
;;
;;    endif
;;    
;;    hahb_true = atlasdust[indx].h_alpha[0]/atlasdust[indx].h_beta[0]
;;    hahb_uncor = atlasdust[indx].h_alpha_uncor[0]/atlasdust[indx].h_beta_uncor[0]
;;
;;    hb_ew_abs     = 3.5 ; EW [mean absorption correction]
;;    hb_ew_abs_err = 1.0
;;
;;    xrange = LBrange
;;    yrange = [-2.8,0.45]
;;
;;    xtitle = 'log L(B) [L(B)'+sunsymbol()+']'
;;    ytitle = 'log (H\beta/H\alpha)'
;;
;;    pagemaker, nx=2, ny=2, yspace=0, xspace=0, xmargin=[1.25,0.35], $
;;      ymargin=[0.9,1.1], position=pos, /normal
;;
;;; ##########################################################
;;; Panel 1: No absorption correction, no reddening correction
;;; ##########################################################
;;
;;    hb = atlasdust[indx].h_beta_uncor[0]
;;    hb_err = atlasdust[indx].h_beta_uncor[1]
;;    
;;    y = alog10(hb/ha)
;;    yerr = im_compute_error(hb,hb_err,ha,ha_err,/log)
;;
;;    if keyword_set(nfgs) then begin
;;
;;       hbnfgs = nfgsdust[indxnfgs].h_beta_uncor[0]
;;       hbnfgs_err = nfgsdust[indxnfgs].h_beta_uncor[1]
;;       
;;       ynfgs = alog10(hbnfgs/hanfgs)
;;       yerrnfgs = im_compute_error(hbnfgs,hbnfgs_err,hanfgs,hanfgs_err,/log)
;;       
;;    endif
;;
;;    stats = im_stats(y,/verbose)
;;
;;    atlas1d_lineplot, x, y, xerr, yerr, plottype=3, postscript=postscript, $
;;      xtitle='', ytitle=ytitle, xrange=xrange, yrange=yrange, legendtype=0, $
;;      charsize=charsize, xstyle=11, xtickname=replicate(' ',10), position=pos[*,0], $
;;      xnfgs=xnfgs, ynfgs=ynfgs, xerrnfgs=xerrnfgs, yerrnfgs=yerrnfgs
;;    axis, /xaxis, xrange=interpol(xabs,x,!x.crange), xthick=postthick, charthick=postthick, $
;;      charsize=charsize, xtitle = textoidl('M_{B} [mag]'), xsty=1
;;    djs_oplot, !x.crange, alog10(1.0/HaHb)*[1,1], line=2, thick=postthick
;;
;;    legend, '(a)', /left, /top, box=0, charsize=charsize, charthick=postthick
;;
;;    label = ['No Abs. Correction','A(H\beta) = 0','A(H\alpha) = 0']
;;;   label = ['No Abs. Correction','No Dust Correction']
;;    legend, textoidl(label), /left, /bottom, box=0, charsize=charsize, $
;;      charthick=postthick
;;
;;; ############################################################
;;; Panel 2: Mean absorption correction, no reddening correction
;;; ############################################################
;;
;;    hb = atlasdust[indx].h_beta_uncor[0] + hb_ew_abs*atlasdust[indx].babs_h_beta_continuum[0]
;;    hb_err = sqrt(atlasdust[indx].h_beta[1]^2 + (hb_ew_abs_err*atlasdust[indx].babs_h_beta_continuum[0])^2)
;;    
;;    y = alog10(hb/ha)
;;    yerr = im_compute_error(hb,hb_err,ha,ha_err,/log)
;;
;;    if keyword_set(nfgs) then begin
;;
;;       hbnfgs = nfgsdust[indxnfgs].h_beta_uncor[0] + hb_ew_abs*nfgsdust[indxnfgs].babs_h_beta_continuum[0]
;;       hbnfgs_err = nfgsdust[indxnfgs].h_beta[1]
;;       
;;       ynfgs = alog10(hbnfgs/hanfgs)
;;       yerrnfgs = im_compute_error(hbnfgs,hbnfgs_err,hanfgs,hanfgs_err,/log)
;;
;;    endif
;;
;;    stats = im_stats(y,/verbose,/no_head)
;;
;;    atlas1d_lineplot, x, y, xerr, yerr, plottype=3, postscript=postscript, $
;;      xtitle='', ytitle='', xrange=xrange, yrange=yrange, legendtype=0, $
;;      charsize=charsize, xstyle=11, xtickname=replicate(' ',10), ytickname=replicate(' ',10), $
;;      position=pos[*,1], /noerase, xnfgs=xnfgs, ynfgs=ynfgs, xerrnfgs=xerrnfgs, yerrnfgs=yerrnfgs
;;    axis, /xaxis, xrange=interpol(xabs,x,!x.crange), xthick=postthick, charthick=postthick, $
;;      charsize=charsize, xtitle = textoidl('M_{B} [mag]'), xsty=1
;;    djs_oplot, !x.crange, alog10(1.0/HaHb)*[1,1], line=2, thick=postthick
;;
;;    legend, '(b)', /left, /top, box=0, charsize=charsize, charthick=postthick
;;
;;    label = ['Mean Abs. Correction','No Dust Correction']
;;    legend, textoidl(label), /left, /bottom, box=0, charsize=charsize, $
;;      charthick=postthick
;;    
;;; ############################################################
;;; Panel 3: Good absorption correction, no reddening correction
;;; ############################################################
;;
;;    hb = atlasdust[indx].h_beta[0]
;;    hb_err = atlasdust[indx].h_beta[1]
;;    
;;    y = alog10(hb/ha)
;;    yerr = im_compute_error(hb,hb_err,ha,ha_err,/log)
;;
;;    if keyword_set(nfgs) then begin
;;
;;       hbnfgs = nfgsdust[indxnfgs].h_beta[0]
;;       hbnfgs_err = nfgsdust[indxnfgs].h_beta[1]
;;    
;;       ynfgs = alog10(hbnfgs/hanfgs)
;;       yerrnfgs = im_compute_error(hbnfgs,hbnfgs_err,hanfgs,hanfgs_err,/log)
;;
;;    endif
;;
;;    stats = im_stats(y,/verbose,/no_head)
;;    
;;    atlas1d_lineplot, x, y, xerr, yerr, plottype=3, postscript=postscript, $
;;      xtitle=xtitle, ytitle=ytitle, xrange=xrange, yrange=yrange, legendtype=0, $
;;      charsize=charsize, position=pos[*,2], /noerase, xnfgs=xnfgs, ynfgs=ynfgs, $
;;      xerrnfgs=xerrnfgs, yerrnfgs=yerrnfgs
;;    djs_oplot, !x.crange, alog10(1.0/HaHb)*[1,1], line=2, thick=postthick
;;
;;    legend, '(c)', /left, /top, box=0, charsize=charsize, charthick=postthick
;;
;;    label = ['Individual Abs. Correction','No Dust Correction']
;;    legend, textoidl(label), /left, /bottom, box=0, charsize=charsize, $
;;      charthick=postthick
;;    
;;; ###############################################################
;;; Panel 4: Mean absorption correction, Mean reddening correction
;;; ###############################################################
;;
;;    Aha = 1.0
;;    Aha_err = 0.1
;;
;;;   yha = atlasdust[indx].h_alpha[0]*10^(0.4*Aha)
;;;   yha_err = atlasdust[indx].h_alpha[1]*10^(0.4*Aha) + $
;;;     atlasdust[indx].h_alpha[0]*alog(10.0)*10^(0.4*Aha)*Aha_err
;;
;;    Ahb = Aha*k_lambda(4861,/odonnell)/k_lambda(6563,/odonnell)
;;    Ahb_err = 0.1*k_lambda(4861,/odonnell)/k_lambda(6563,/odonnell)
;;
;;    mhb = atlasdust[indx].h_beta_uncor[0] + hb_ew_abs*atlasdust[indx].babs_h_beta_continuum[0]
;;    mhb_err = sqrt(atlasdust[indx].h_beta[1]^2 + (hb_ew_abs_err*atlasdust[indx].babs_h_beta_continuum[0])^2)
;;
;;    hb = mhb*10^(0.4*Ahb)
;;    hb_err = mhb_err*10^(0.4*Ahb) + mhb*alog(10.0)*10^(0.4*Ahb)*Ahb_err
;;    
;;    y = alog10(hb/ha)
;;    yerr = im_compute_error(hb,hb_err,ha,ha_err,/log)
;;
;;    if keyword_set(nfgs) then begin
;;
;;       mhbnfgs = nfgsdust[indxnfgs].h_beta_uncor[0] + hb_ew_abs*nfgsdust[indxnfgs].babs_h_beta_continuum[0]
;;       mhbnfgs_err = nfgsdust[indxnfgs].h_beta[1]
;;       
;;       hbnfgs = mhbnfgs*10^(0.4*Ahb)
;;       hbnfgs_err = mhbnfgs_err*10^(0.4*Ahb) + mhbnfgs*alog(10.0)*10^(0.4*Ahb)*Ahb_err
;;
;;;      yhanfgs = nfgsdust[indxnfgs].h_alpha[0]*10^(0.4*Aha)
;;;      yhanfgs_err = nfgsdust[indxnfgs].h_alpha[1]*10^(0.4*Aha) + $
;;;        nfgsdust[indxnfgs].h_alpha[0]*alog(10.0)*10^(0.4*Aha)*Aha_err
;;
;;       ynfgs = alog10(hbnfgs/hanfgs)
;;       yerrnfgs = im_compute_error(hbnfgs,hbnfgs_err,hanfgs,hanfgs_err,/log)
;;
;;    endif
;;       
;;    stats = im_stats(y,/verbose,/no_head)
;;
;;    atlas1d_lineplot, x, y, xerr, yerr, plottype=3, postscript=postscript, $
;;      xtitle=xtitle, ytitle='', xrange=xrange, yrange=yrange, legendtype=0, $
;;      charsize=charsize, ytickname=replicate(' ',10), position=pos[*,3], /noerase, $
;;      xnfgs=xnfgs, ynfgs=ynfgs, xerrnfgs=xerrnfgs, yerrnfgs=yerrnfgs      
;;    djs_oplot, !x.crange, alog10(1.0/HaHb)*[1,1], line=2, thick=postthick
;;
;;    legend, '(d)', /left, /top, box=0, charsize=charsize, charthick=postthick
;;
;;    label = ['Mean Abs. Correction','Mean Dust Correction']
;;    legend, textoidl(label), /left, /bottom, box=0, charsize=charsize, $
;;      charthick=postthick
;;
;;    im_openclose, postscript=postscript, /close

; ------------------------------------------------------------
; High-Redshift E(Hb-Ha) vs [O II]/Ha_obs
; ------------------------------------------------------------

    psname = 'highz_ebv_vs_oiihb_obs'
    im_openclose, pspath+psname, postscript=postscript

    cut = where(atlasnodust.ebv_hahb_err gt 0.0)
    lineratio, atlasdust[cut], '', '', 'OII_3727', 'H_BETA', $
      dum1, dum2, y, yerr, index=indx, nindex=nindx, snrcut=snrcut, $
      xsyserr=syserr, ysyserr=syserr

    x = atlasnodust[cut[indx]].ebv_hahb
    xerr = atlasnodust[cut[indx]].ebv_hahb_err

    xbig = x
    ybig = y
    
    if keyword_set(nfgs) then begin

       cutnfgs = where(nfgsnodust.ebv_hahb_err gt 0.0)
       lineratio, nfgsdust[cutnfgs], '', '', 'OII_3727', 'H_BETA', $
         dum1, dum2, ynfgs, yerrnfgs, index=indxnfgs, snrcut=snrcut, $
      xsyserr=nfgssyserr, ysyserr=nfgssyserr
       xnfgs = nfgsnodust[cutnfgs[indxnfgs]].ebv_hahb
       xerrnfgs = nfgsnodust[cutnfgs[indxnfgs]].ebv_hahb_err

       xbig = [xbig,xnfgs]
       ybig = [ybig,ynfgs]
       
    endif

; ####################
; Liang et al. (2004)    
; ####################
    
    cutliang = where(liangnodust.ebv_hbhg_err gt 0.0)
    lineratio, liang04[cutliang], '', '', 'OII_3727', 'H_BETA', dum1, dum2, $
      yliang, yerrliang, index=indxliang, nindex=nindxliang, snrcut=snrcut_highz

    xliang = liangnodust[cutliang[indxliang]].ebv_hbhg
    xerrliang = liangnodust[cutliang[indxliang]].ebv_hbhg_err

;;; ####################
;;; H02 et al. (2002)    
;;; ####################
;;    
;;    cuth02 = where(h02nodust.ebv_hbhg_err gt 0.0)
;;    lineratio, h02[cuth02], '', '', 'OII_3727', 'H_BETA', dum1, dum2, $
;;      yh02, yerrh02, index=indxh02, nindex=nindxh02, snrcut=snrcut_highz
;;
;;    xh02 = h02nodust[cuth02[indxh02]].ebv_hbhg
;;    xerrh02 = h02nodust[cuth02[indxh02]].ebv_hbhg_err

    xtitle = 'E(B-V) [mag]'
    ytitle = 'log ([O II]/H\beta)_{obs}'

    xrange = [-0.1,1.5]
    yrange = oiihbrange

    djs_plot, [0], [0], /nodata, xthick=postthick, ythick=postthick, charsize=2.0, $
      charthick=postthick, thick=postthick, xtitle=xtitle, ytitle=ytitle, $
      xstyle=3, ystyle=3, xrange=xrange, yrange=yrange

    plotsym, 0, 1.0, /fill
    djs_oplot, xbig, ybig, ps=8, color='grey'

    plotsym, 8, 1.5, /fill
    oploterror, xliang, yliang, xerrliang, yerrliang, ps=8, thick=postthick, $
      color=djs_icolor('orange'), errcolor=djs_icolor('orange'), /nohat, errthick=0.5
    
;   plotsym, 8, 1.5, /fill
;   oploterror, xh02, yh02, xerrh02, yerrh02, ps=8, thick=postthick, $
;     color=djs_icolor('dark green'), errcolor=djs_icolor('dark green'), /nohat, errthick=0.5

    im_openclose, postscript=postscript, /close    

; ------------------------------------------------------------
; This plot was Rob's idea.  Essentially, this plot provides a
; visualization of the various SFR calibrations by representing them
; as Gaussian functions with different mean values and standard
; deviations. 
; ------------------------------------------------------------
    
    psname = 'visualize_sfr_calibrations'
    im_openclose, pspath+psname, postscript=postscript

    xrange = [-1.0,1.0]
    yrange = [0,0.04]

    xtitle = 'L/L(H\alpha)'
    ytitle = 'Probability'

    ratioaxis = findgen((xrange[1]-xrange[0])/0.01+1)*0.01+xrange[0]

    djs_plot, [0], [0], /nodata, xrange=xrange, yrange=yrange, xtitle=xtitle, $
      ytitle=ytitle, charsize=2.0, charthick=postthick, xthick=postthick, $
      ythick=postthick, ytickname=replicate(' ',10)

; --------------------
; [O II]/Ha observed    
; --------------------
    
    oii_obs_median = 0.01 & oii_obs_sigma = 0.12
    oii_obs = gauss1(ratioaxis,[oii_obs_median,oii_obs_sigma,1.0])
    oii_obs = oii_obs/total(oii_obs)

    djs_oplot, ratioaxis, oii_obs, line=0
    
; --------------------
; [O II]/Ha corrected - low metallicity
; --------------------
    
    oii_cor_lowZ_median = -0.05 & oii_cor_lowZ_sigma = 0.22
    oii_cor_lowZ = gauss1(ratioaxis,[oii_cor_lowZ_median,oii_cor_lowZ_sigma,1.0])
    oii_cor_lowZ = oii_cor_lowZ/total(oii_cor_lowZ)
    
    djs_oplot, ratioaxis, oii_cor_lowZ, line=0

; --------------------
; [O II]/Ha corrected - intermediate metallicity
; --------------------
    
    oii_cor_interZ_median = 0.5 & oii_cor_interZ_sigma = 0.15
    oii_cor_interZ = gauss1(ratioaxis,[oii_cor_interZ_median,oii_cor_interZ_sigma,1.0])
    oii_cor_interZ = oii_cor_interZ/total(oii_cor_interZ)
    
    djs_oplot, ratioaxis, oii_cor_interZ, line=0

; --------------------
; [O II]/Ha corrected - high metallicity
; --------------------
    
    oii_cor_highZ_median = 0.1 & oii_cor_highZ_sigma = 0.2
    oii_cor_highZ = gauss1(ratioaxis,[oii_cor_highZ_median,oii_cor_highZ_sigma,1.0])
    oii_cor_highZ = oii_cor_highZ/total(oii_cor_highZ)

    djs_oplot, ratioaxis, oii_cor_highZ, line=0
    
; --------------------
; Ha corrected
; --------------------
    
    ha_cor_median = -0.1 & ha_cor_sigma = 0.20
    ha_cor = gauss1(ratioaxis,[ha_cor_median,ha_cor_sigma,1.0])
    ha_cor = ha_cor/total(ha_cor)

    djs_oplot, ratioaxis, ha_cor, line=0
    
    im_openclose, postscript=postscript, /close

;;;; ------------------------------------------------------------
;;;; 3-panel [O II] plot - E(B-V)
;;;; ------------------------------------------------------------
;;;
;;;    psname = 'oiiha_vs_ebv_multipanel'
;;;    im_openclose, pspath+psname, postscript=postscript
;;;
;;;    indx = where((atlasdust.oii_3727[0]/atlasdust.oii_3727[1] gt snrcut) and $
;;;      (atlasnodust.oii_3727[0]/atlasnodust.oii_3727[1] gt snrcut) and $
;;;      (atlasdust.h_alpha[0]/atlasdust.h_alpha[1] gt snrcut) and $
;;;      (atlasnodust.h_alpha[0]/atlasnodust.h_alpha[1] gt snrcut) and $
;;;      (atlasnodust.ehbha_err gt 0.0),nindx)
;;;
;;;    x = atlasnodust[indx].ehbha
;;;    xerr = atlasnodust[indx].ehbha_err
;;;    xebv = x/(k_lambda(4861,/odonnell)-k_lambda(6563,/odonnell))
;;;
;;;    if keyword_set(nfgs) then begin
;;;
;;;       indxnfgs = where((nfgsdust.oii_3727[0]/nfgsdust.oii_3727[1] gt snrcut) and $
;;;         (nfgsnodust.oii_3727[0]/nfgsnodust.oii_3727[1] gt snrcut) and $
;;;         (nfgsdust.h_alpha[0]/nfgsdust.h_alpha[1] gt snrcut) and $
;;;         (nfgsnodust.h_alpha[0]/nfgsnodust.h_alpha[1] gt snrcut) and $
;;;         (nfgsnodust.ehbha_err gt 0.0),nindxnfgs)
;;;
;;;       xnfgs = nfgsnodust[indxnfgs].ehbha
;;;       xerrnfgs = nfgsnodust[indxnfgs].ehbha_err
;;;
;;;    endif
;;;
;;;    xtitle = 'E(H\beta-H\alpha) [mag]'
;;;    ytitle = 'log ([O II]/H\alpha)'
;;;
;;;    xrange = ehbharange
;;;    yrange = oiihacorrange
;;;
;;;    pagemaker, nx=3, ny=1, yspace=0, xspace=0, xmargin=[1.25,0.35], $
;;;      ymargin=[1.1,1.1], width=2, height=2, position=pos, /normal
;;;
;;;; ##########################################################
;;;; Panel 1: [O II] observed, Ha observed
;;;; ##########################################################
;;;
;;;    y1 = atlasdust[indx].oii_3727[0]
;;;    y1err = atlasdust[indx].oii_3727[1]
;;;
;;;    y2 = atlasdust[indx].h_alpha[0]
;;;    y2err = atlasdust[indx].h_alpha[1]
;;;
;;;    y = alog10(y1/y2)
;;;    yerr = im_compute_error(y1,y1err,y2,y2err,/log)
;;;
;;;    xbig = x
;;;    ybig = y
;;;    
;;;    if keyword_set(nfgs) then begin
;;;       
;;;       y1nfgs = nfgsdust[indxnfgs].oii_3727[0]
;;;       y1errnfgs = nfgsdust[indxnfgs].oii_3727[1]
;;;
;;;       y2nfgs = nfgsdust[indxnfgs].h_alpha[0]
;;;       y2errnfgs = nfgsdust[indxnfgs].h_alpha[1]
;;;
;;;       ynfgs = alog10(y1nfgs/y2nfgs)
;;;       yerrnfgs = im_compute_error(y1nfgs,y1errnfgs,y2nfgs,y2errnfgs,/log)
;;;
;;;       ybig = [ybig,ynfgs]
;;;       xbig = [xbig,xnfgs]
;;;
;;;    endif
;;;
;;;    stats = im_stats(ybig,/verbose)
;;;    rcor = r_correlate(xbig,ybig,zd=zd,probd=probd)
;;;    splog, 'Spearman rank for E(Hb-Ha) vs [O II]/Ha observed: ', rcor, probd
;;;
;;;    atlas1d_lineplot, x, y, xerr, yerr, plottype=3, postscript=postscript, $
;;;      ytitle=ytitle, xtitle=xtitle, xrange=xrange, yrange=yrange, legendtype=0, xminor=3, yminor=3, $
;;;      xstyle=11, /right, /top, xnfgs=xnfgs, ynfgs=ynfgs, xerrnfgs=xerrnfgs, $
;;;      yerrnfgs=yerrnfgs, position=pos[*,0], charsize=1.5
;;;    axis, /xaxis, xrange=interpol(xebv,x,!x.crange), xthick=postthick, charthick=postthick, $
;;;      charsize=1.5, xtitle = 'E(B-V) [mag]', xsty=1
;;;
;;;    legend, '(a)', /left, /top, box=0, charsize=charsize, charthick=postthick
;;;
;;;    label = ['[O II] Observed','H\alpha Observed']
;;;    legend, textoidl(label), /left, /bottom, box=0, charsize=charsize, $
;;;      charthick=postthick
;;;
;;;; overplot reddening relations
;;;
;;;    ratio = 0.0
;;;
;;;    x_ebv = findgen((1.1-(-0.5))/0.1)*0.1+(-0.5)
;;;    yodonnell =  -0.4*x_ebv*(k_lambda(3727.0,/odonnell)-k_lambda(6563.0,/odonnell)) + ratio
;;;    ysmc = -0.4*x_ebv*(k_lambda(3727,/smc)-k_lambda(6563,/smc)) + ratio
;;;
;;;    djs_oplot, x_ebv, yodonnell, line=0, thick=postthick
;;;    djs_oplot, x_ebv, ysmc, line=2, thick=postthick
;;;
;;;    label = ['Milky Way','SMC Bar']
;;;    linestyle = [0,2]
;;;    legend, label, /right, /top, box=0, linestyle=linestyle, $
;;;      charsize=charsize, charthick=postthick, thick=postthick, clear=postscript
;;;
;;;; ##########################################################
;;;; Panel 2: [O II] observed, Ha corrected
;;;; ##########################################################
;;;
;;;    y1 = atlasdust[indx].oii_3727[0]
;;;    y1err = atlasdust[indx].oii_3727[1]
;;;
;;;    y2 = atlasnodust[indx].h_alpha[0]
;;;    y2err = atlasnodust[indx].h_alpha[1]
;;;
;;;    y = alog10(y1/y2)
;;;    yerr = im_compute_error(y1,y1err,y2,y2err,/log)
;;;
;;;    xbig = x
;;;    ybig = y
;;;    
;;;    if keyword_set(nfgs) then begin
;;;       
;;;       y1nfgs = nfgsdust[indxnfgs].oii_3727[0]
;;;       y1errnfgs = nfgsdust[indxnfgs].oii_3727[1]
;;;
;;;       y2nfgs = nfgsnodust[indxnfgs].h_alpha[0]
;;;       y2errnfgs = nfgsnodust[indxnfgs].h_alpha[1]
;;;
;;;       ynfgs = alog10(y1nfgs/y2nfgs)
;;;       yerrnfgs = im_compute_error(y1nfgs,y1errnfgs,y2nfgs,y2errnfgs,/log)
;;;
;;;       ybig = [ybig,ynfgs]
;;;       xbig = [xbig,xnfgs]
;;;
;;;    endif
;;;
;;;    stats = im_stats(ybig,/verbose)
;;;
;;;    atlas1d_lineplot, x, y, xerr, yerr, plottype=3, postscript=postscript, $
;;;      xrange=xrange, yrange=yrange, legendtype=0, /noerase, xtitle=xtitle, ytickname=replicate(' ',10), $
;;;      /right, /top, xnfgs=xnfgs, ynfgs=ynfgs, xerrnfgs=xerrnfgs, xminor=3, yminor=3, $
;;;      yerrnfgs=yerrnfgs, position=pos[*,1], charsize=1.5, ytitle=''
;;;    
;;;    legend, '(b)', /left, /top, box=0, charsize=charsize, charthick=postthick
;;;
;;;    label = ['[O II] Observed','H\alpha Corrected']
;;;    legend, textoidl(label), /left, /bottom, box=0, charsize=charsize, $
;;;      charthick=postthick
;;;
;;;; overplot a histogram
;;;
;;;    binsize = 0.1
;;;    plothist, ybig, bin=binsize, xbin, ybin, /noplot
;;;
;;;    yhistrange = minmax(ybin)*[1.0,1.05]
;;;
;;;; ##########################################################
;;;; Panel 3: [O II] corrected, Ha corrected
;;;; ##########################################################
;;;
;;;    y1 = atlasnodust[indx].oii_3727[0]
;;;    y1err = atlasnodust[indx].oii_3727[1]
;;;
;;;    y2 = atlasnodust[indx].h_alpha[0]
;;;    y2err = atlasnodust[indx].h_alpha[1]
;;;
;;;    y = alog10(y1/y2)
;;;    yerr = im_compute_error(y1,y1err,y2,y2err,/log)
;;;
;;;    xbig = x
;;;    ybig = y
;;;    
;;;    if keyword_set(nfgs) then begin
;;;       
;;;       y1nfgs = nfgsnodust[indxnfgs].oii_3727[0]
;;;       y1errnfgs = nfgsnodust[indxnfgs].oii_3727[1]
;;;
;;;       y2nfgs = nfgsnodust[indxnfgs].h_alpha[0]
;;;       y2errnfgs = nfgsnodust[indxnfgs].h_alpha[1]
;;;
;;;       ynfgs = alog10(y1nfgs/y2nfgs)
;;;       yerrnfgs = im_compute_error(y1nfgs,y1errnfgs,y2nfgs,y2errnfgs,/log)
;;;
;;;       ybig = [ybig,ynfgs]
;;;       xbig = [xbig,xnfgs]
;;;
;;;    endif
;;;
;;;    stats = im_stats(ybig,/verbose)
;;;    rcor = r_correlate(xbig,ybig,zd=zd,probd=probd)
;;;    splog, 'Spearman rank for E(Hb-Ha) vs [O II]/Ha corrected: ', rcor, probd
;;;
;;;    atlas1d_lineplot, x, y, xerr, yerr, plottype=3, postscript=postscript, $
;;;      xtitle=xtitle, xrange=xrange, yrange=yrange, legendtype=0, /noerase, $
;;;      xstyle=3, /right, /top, xnfgs=xnfgs, ynfgs=ynfgs, xerrnfgs=xerrnfgs, $
;;;      yerrnfgs=yerrnfgs, position=pos[*,2], charsize=1.5, xminor=3, yminor=3, $
;;;      ytickname=replicate(' ',10)
;;;    
;;;    legend, '(c)', /left, /top, box=0, charsize=charsize, charthick=postthick
;;;
;;;    label = ['[O II] Corrected','H\alpha Corrected']
;;;    legend, textoidl(label), /left, /bottom, box=0, charsize=charsize, $
;;;      charthick=postthick
;;;
;;;    im_openclose, postscript=postscript, /close
    
; ------------------------------------------------------------
; 4-panel [O II] plot - E(B-V)
; ------------------------------------------------------------

    psname = 'oiiha_vs_ebv_4panel'
    im_openclose, pspath+psname, postscript=postscript

    indx = where((atlasdust.oii_3727[0]/atlasdust.oii_3727[1] gt snrcut) and $
      (atlasnodust.oii_3727[0]/atlasnodust.oii_3727[1] gt snrcut) and $
      (atlasdust.h_alpha[0]/atlasdust.h_alpha[1] gt snrcut) and $
      (atlasnodust.h_alpha[0]/atlasnodust.h_alpha[1] gt snrcut) and $
      (atlasnodust.ehbha_err gt 0.0),nindx)

    x = atlasnodust[indx].ehbha
    xerr = atlasnodust[indx].ehbha_err
    xebv = x/(k_lambda(4861,/odonnell)-k_lambda(6563,/odonnell))

    if keyword_set(nfgs) then begin

       indxnfgs = where((nfgsdust.oii_3727[0]/nfgsdust.oii_3727[1] gt snrcut) and $
         (nfgsnodust.oii_3727[0]/nfgsnodust.oii_3727[1] gt snrcut) and $
         (nfgsdust.h_alpha[0]/nfgsdust.h_alpha[1] gt snrcut) and $
         (nfgsnodust.h_alpha[0]/nfgsnodust.h_alpha[1] gt snrcut) and $
         (nfgsnodust.ehbha_err gt 0.0),nindxnfgs)

       xnfgs = nfgsnodust[indxnfgs].ehbha
       xerrnfgs = nfgsnodust[indxnfgs].ehbha_err

    endif

    xtitle = 'E(H\beta-H\alpha) [mag]'
    ytitle = 'log ([O II]/H\alpha)'

    xrange = ehbharange
    yrange = [-2.2,1.1] ; oiihacorrange

    pagemaker, nx=2, ny=2, yspace=0, xspace=0, xmargin=[1.25,0.35], $
      ymargin=[0.8,1.1], position=pos, /normal

; ##########################################################
; Panel 1: [O II] observed, Ha observed
; ##########################################################

    y1 = atlasdust[indx].oii_3727[0]
    y1err = atlasdust[indx].oii_3727[1]

    y2 = atlasdust[indx].h_alpha[0]
    y2err = atlasdust[indx].h_alpha[1]

    y = alog10(y1/y2)
    yerr = im_compute_error(y1,y1err,y2,y2err,/log)

    xbig = x
    ybig = y
    
    if keyword_set(nfgs) then begin
       
       y1nfgs = nfgsdust[indxnfgs].oii_3727[0]
       y1errnfgs = nfgsdust[indxnfgs].oii_3727[1]

       y2nfgs = nfgsdust[indxnfgs].h_alpha[0]
       y2errnfgs = nfgsdust[indxnfgs].h_alpha[1]

       ynfgs = alog10(y1nfgs/y2nfgs)
       yerrnfgs = im_compute_error(y1nfgs,y1errnfgs,y2nfgs,y2errnfgs,/log)

       ybig = [ybig,ynfgs]
       xbig = [xbig,xnfgs]

    endif

    stats = im_stats(ybig,/verbose)
    rcor = r_correlate(xbig,ybig,zd=zd,probd=probd)
    splog, 'Spearman rank for E(Hb-Ha) vs [O II]/Ha observed: ', rcor, probd

    atlas1d_lineplot, x, y, xerr, yerr, plottype=3, postscript=postscript, $
      ytitle=ytitle, xrange=xrange, yrange=yrange, legendtype=0, $
      xstyle=11, /right, /top, xnfgs=xnfgs, ynfgs=ynfgs, xerrnfgs=xerrnfgs, $
      yerrnfgs=yerrnfgs, position=pos[*,0], charsize=charsize, xtickname=replicate(' ',10)
    axis, /xaxis, xrange=interpol(xebv,x,!x.crange), xthick=postthick, charthick=postthick, $
      charsize=charsize, xtitle = 'E(B-V) [mag]', xsty=1

    legend, '(a)', /right, /top, box=0, charsize=charsize, charthick=postthick

    label = ['[O II] Observed','H\alpha Observed']
    legend, textoidl(label), /left, /bottom, box=0, charsize=charsize, $
      charthick=postthick

; overplot reddening relations

    ratio = 0.0

    x_ebv = findgen((1.1-(-0.5))/0.1)*0.1+(-0.5)
    yodonnell =  -0.4*x_ebv*(k_lambda(3727.0,/odonnell)-k_lambda(6563.0,/odonnell)) + ratio
    ysmc = -0.4*x_ebv*(k_lambda(3727,/smc)-k_lambda(6563,/smc)) + ratio

    djs_oplot, x_ebv, yodonnell, line=0, thick=postthick
    djs_oplot, x_ebv, ysmc, line=2, thick=postthick

    label = ['Milky Way','SMC Bar']
    linestyle = [0,2]
    legend, label, /left, /top, box=0, linestyle=linestyle, $
      charsize=charsize_3, charthick=postthick, thick=postthick, clear=postscript

; overplot a histogram

    binsize = 0.1
    plothist, ybig, bin=binsize, xbin, ybin, /noplot, /halfbin

    yhistrange = minmax(ybin)*[1.0,1.05]

;   djs_plot, [0], [0], /nodata, /noerase, xthick=postthick, ythick=postthick, charsize=0.8, $
;     charthick=postthick, thick=postthick, xtitle=textoidl('log ([O II]/H\alpha)'), ytitle='Number', $
;     xstyle=1, ystyle=1, xrange=yrange, yrange=yhistrange, /normal, $
;     position=[0.2,0.55,0.35,0.7]
;   plothist, ybig, bin=binsize, /overplot, thick=postthick, /fill, /fline, /halfbin, $
;     forientation=45, fcolor=djs_icolor('grey'), color=djs_icolor('grey')

; ##########################################################
; Panel 2: [O II] observed, Ha corrected
; ##########################################################

    y1 = atlasdust[indx].oii_3727[0]
    y1err = atlasdust[indx].oii_3727[1]

    y2 = atlasnodust[indx].h_alpha[0]
    y2err = atlasnodust[indx].h_alpha[1]

    y = alog10(y1/y2)
    yerr = im_compute_error(y1,y1err,y2,y2err,/log)

    xbig = x
    ybig = y
    
    if keyword_set(nfgs) then begin
       
       y1nfgs = nfgsdust[indxnfgs].oii_3727[0]
       y1errnfgs = nfgsdust[indxnfgs].oii_3727[1]

       y2nfgs = nfgsnodust[indxnfgs].h_alpha[0]
       y2errnfgs = nfgsnodust[indxnfgs].h_alpha[1]

       ynfgs = alog10(y1nfgs/y2nfgs)
       yerrnfgs = im_compute_error(y1nfgs,y1errnfgs,y2nfgs,y2errnfgs,/log)

       ybig = [ybig,ynfgs]
       xbig = [xbig,xnfgs]

    endif

    stats = im_stats(ybig,/verbose)

    atlas1d_lineplot, x, y, xerr, yerr, plottype=3, postscript=postscript, $
      xrange=xrange, yrange=yrange, legendtype=0, /noerase, xtickname=replicate(' ',10), $
      xstyle=11, /right, /top, xnfgs=xnfgs, ynfgs=ynfgs, xerrnfgs=xerrnfgs, $
      yerrnfgs=yerrnfgs, position=pos[*,1], charsize=charsize, ytickname=replicate(' ',10)
    axis, /xaxis, xrange=interpol(xebv,x,!x.crange), xthick=postthick, charthick=postthick, $
      charsize=charsize, xtitle = 'E(B-V) [mag]', xsty=1
    
    legend, '(b)', /right, /top, box=0, charsize=charsize, charthick=postthick

    label = ['[O II] Observed','H\alpha Corrected']
    legend, textoidl(label), /left, /bottom, box=0, charsize=charsize, $
      charthick=postthick

; overplot a histogram

    binsize = 0.1
    plothist, ybig, bin=binsize, xbin, ybin, /noplot, /halfbin

    yhistrange = minmax(ybin)*[1.0,1.05]

;   djs_plot, [0], [0], /nodata, /noerase, xthick=postthick, ythick=postthick, charsize=1.0, $
;     charthick=postthick, thick=postthick, xtitle=textoidl('log ([O II]/H\alpha)'), ytitle='Number', $
;     xstyle=1, ystyle=1, xrange=yrange, yrange=yhistrange, /normal, $
;     position=[0.7,0.55,0.85,0.7]
;   plothist, ybig, bin=binsize, /overplot, thick=postthick, /fill, /fline, /halfbin, $
;     forientation=45, fcolor=djs_icolor('grey'), color=djs_icolor('grey')

; ##########################################################
; Panel 3: [O II] mean corrected, Ha corrected
; ##########################################################

    Aha = 1.0
    Aha_err = 0.1
    
    Aoii = Aha*k_lambda(3727,/odonnel)/k_lambda(6563,/odonnel)
    Aoii_err = Aha_err*k_lambda(3727,/odonnel)/k_lambda(6563,/odonnel)

    y1 = atlasdust[indx].oii_3727[0]*10^(0.4*Aoii)
    y1err = atlasdust[indx].oii_3727[1]*10^(0.4*Aoii) + $
      atlasdust[indx].oii_3727[0]*alog(10.0)*10^(0.4*Aoii)*Aoii_err

    y2 = atlasnodust[indx].h_alpha[0]
    y2err = atlasnodust[indx].h_alpha[1]

    y = alog10(y1/y2)
    yerr = im_compute_error(y1,y1err,y2,y2err,/log)

    xbig = x
    ybig = y
    
    stats = im_stats(ybig,/verbose)

    if keyword_set(nfgs) then begin
       
       y1nfgs = nfgsdust[indxnfgs].oii_3727[0]*10^(0.4*Aoii)
       y1errnfgs = nfgsdust[indxnfgs].oii_3727[1]*10^(0.4*Aoii) + $
         nfgsdust[indxnfgs].oii_3727[0]*alog(10.0)*10^(0.4*Aoii)*Aoii_err

       y2nfgs = nfgsnodust[indxnfgs].h_alpha[0]
       y2errnfgs = nfgsnodust[indxnfgs].h_alpha[1]

       ynfgs = alog10(y1nfgs/y2nfgs)
       yerrnfgs = im_compute_error(y1nfgs,y1errnfgs,y2nfgs,y2errnfgs,/log)

       ybig = [ybig,ynfgs]
       xbig = [xbig,xnfgs]

    endif

    atlas1d_lineplot, x, y, xerr, yerr, plottype=3, postscript=postscript, $
      xtitle=xtitle, ytitle=ytitle, xrange=xrange, yrange=yrange, legendtype=0, /noerase, $
      xstyle=3, /right, /top, xnfgs=xnfgs, ynfgs=ynfgs, xerrnfgs=xerrnfgs, $
      yerrnfgs=yerrnfgs, position=pos[*,2], charsize=charsize
    
    legend, '(c)', /right, /top, box=0, charsize=charsize, charthick=postthick

    label = ['[O II] Mean Corrected','H\alpha Corrected']
    legend, textoidl(label), /left, /bottom, box=0, charsize=charsize, $
      charthick=postthick

; overplot a histogram

    binsize = 0.1
    plothist, ybig, bin=binsize, xbin, ybin, /noplot, /halfbin

    yhistrange = minmax(ybin)*[1.0,1.05]

;   djs_plot, [0], [0], /nodata, /noerase, xthick=postthick, ythick=postthick, charsize=1.0, $
;     charthick=postthick, thick=postthick, xtitle=textoidl('log ([O II]/H\alpha)'), ytitle='Number', $
;     xstyle=1, ystyle=1, xrange=yrange, yrange=yhistrange, /normal, $
;     position=[0.2,0.55,0.35,0.7]
;   plothist, ybig, bin=binsize, /overplot, thick=postthick, /fill, /fline, /halfbin, $
;     forientation=45, fcolor=djs_icolor('grey'), color=djs_icolor('grey')

; ##########################################################
; Panel 4: [O II] corrected, Ha corrected
; ##########################################################

    y1 = atlasnodust[indx].oii_3727[0]
    y1err = atlasnodust[indx].oii_3727[1]

    y2 = atlasnodust[indx].h_alpha[0]
    y2err = atlasnodust[indx].h_alpha[1]

    y = alog10(y1/y2)
    yerr = im_compute_error(y1,y1err,y2,y2err,/log)

    xbig = x
    ybig = y
    
    if keyword_set(nfgs) then begin
       
       y1nfgs = nfgsnodust[indxnfgs].oii_3727[0]
       y1errnfgs = nfgsnodust[indxnfgs].oii_3727[1]

       y2nfgs = nfgsnodust[indxnfgs].h_alpha[0]
       y2errnfgs = nfgsnodust[indxnfgs].h_alpha[1]

       ynfgs = alog10(y1nfgs/y2nfgs)
       yerrnfgs = im_compute_error(y1nfgs,y1errnfgs,y2nfgs,y2errnfgs,/log)

       ybig = [ybig,ynfgs]
       xbig = [xbig,xnfgs]

    endif

    stats = im_stats(ybig,/verbose)
    rcor = r_correlate(xbig,ybig,zd=zd,probd=probd)
    splog, 'Spearman rank for E(Hb-Ha) vs [O II]/Ha corrected: ', rcor, probd

    atlas1d_lineplot, x, y, xerr, yerr, plottype=3, postscript=postscript, $
      xtitle=xtitle, xrange=xrange, yrange=yrange, legendtype=0, /noerase, $
      xstyle=3, /right, /top, xnfgs=xnfgs, ynfgs=ynfgs, xerrnfgs=xerrnfgs, $
      yerrnfgs=yerrnfgs, position=pos[*,3], charsize=charsize, ytickname=replicate(' ',10)
    
    legend, '(d)', /right, /top, box=0, charsize=charsize, charthick=postthick

    label = ['[O II] Corrected','H\alpha Corrected']
    legend, textoidl(label), /left, /bottom, box=0, charsize=charsize, $
      charthick=postthick

; overplot a histogram

    binsize = 0.1
    plothist, ybig, bin=binsize, xbin, ybin, /noplot, /halfbin

    yhistrange = minmax(ybin)*[1.0,1.05]

;   djs_plot, [0], [0], /nodata, /noerase, xthick=postthick, ythick=postthick, charsize=1.0, $
;     charthick=postthick, thick=postthick, xtitle=textoidl('log ([O II]/H\alpha)'), ytitle='Number', $
;     xstyle=1, ystyle=1, xrange=yrange, yrange=yhistrange, /normal, $
;     position=[0.7,0.55,0.85,0.7]
;   plothist, ybig, bin=binsize, /overplot, thick=postthick, /fill, /fline, /halfbin, $
;     forientation=45, fcolor=djs_icolor('grey'), color=djs_icolor('grey')

    im_openclose, postscript=postscript, /close
    
; ------------------------------------------------------------
; 2-panel f(3727) - E(B-V)
; ------------------------------------------------------------
    
    psname = 'f3727_ehbha_multipanel'
    im_openclose, pspath+psname, postscript=postscript

    indx = where((atlasnodust.h_alpha[0]/atlasnodust.h_alpha[1] gt snrcut) and $   
      (atlasdust.oii_3727_continuum[0]/atlasdust.oii_3727_continuum[1] gt snrcut) and $
      (atlasnodust.oii_3727_continuum[0]/atlasnodust.oii_3727_continuum[1] gt snrcut) and $
      (atlasnodust.ebv_hahb_err gt 0.0),nindx)

    x = atlasnodust[indx].ehbha
    xerr = atlasnodust[indx].ehbha_err
    xebv = x/(k_lambda(4861,/odonnell)-k_lambda(6563,/odonnell))
    
    ha = atlasnodust[indx].h_alpha[0]
    ha_err = atlasnodust[indx].h_alpha[1]
    
    if keyword_set(nfgs) then begin

       indxnfgs = where((nfgsnodust.h_alpha[0]/nfgsnodust.h_alpha[1] gt snrcut) and $   
         (nfgsdust.oii_3727_continuum[0]/nfgsdust.oii_3727_continuum[1] gt snrcut) and $
         (nfgsnodust.oii_3727_continuum[0]/nfgsnodust.oii_3727_continuum[1] gt snrcut) and $
         (nfgsnodust.ebv_hahb_err gt 0.0),nindxnfgs)

       xnfgs = nfgsnodust[indxnfgs].ehbha
       xerrnfgs = nfgsnodust[indxnfgs].ehbha_err
    
       hanfgs = nfgsnodust[indxnfgs].h_alpha[0]
       hanfgs_err = nfgsnodust[indxnfgs].h_alpha[1]
    
    endif
    
    xrange = ehbharange
    yrange = f3727harange

    xtitle = 'E(H\beta-H\alpha) [mag]'

    pagemaker, nx=1, ny=2, yspace=0, xspace=0, xmargin=[1.25,0.35], $
      ymargin=[1.1,1.1], position=pos, /normal

; ##########################################################
; Panel 1: No reddening correction
; ##########################################################

    f3727 = atlasdust[indx].oii_3727_continuum[0]
    f3727_err = atlasdust[indx].oii_3727_continuum[1]
    
    y = alog10(f3727/ha)
    yerr = im_compute_error(f3727,f3727_err,ha,ha_err,/log)

    if keyword_set(nfgs) then begin

       f3727nfgs = nfgsdust[indxnfgs].oii_3727_continuum[0]
       f3727nfgs_err = nfgsdust[indxnfgs].oii_3727_continuum[1]
       
       ynfgs = alog10(f3727nfgs/hanfgs)
       yerrnfgs = im_compute_error(f3727nfgs,f3727nfgs_err,hanfgs,hanfgs_err,/log)
       
    endif

    stats = im_stats(y,/verbose)

    ytitle = 'log (F_{\lambda3727}_{obs}/H\alpha_{cor})'

    atlas1d_lineplot, x, y, xerr, yerr, plottype=3, postscript=postscript, $
      xtitle='', ytitle=ytitle, xrange=xrange, yrange=yrange, legendtype=0, $
      charsize=charsize, position=pos[*,0], xnfgs=xnfgs, ynfgs=ynfgs, $
      xerrnfgs=xerrnfgs, yerrnfgs=yerrnfgs, xtickname=replicate(' ',10), $
      xstyle=11, ymargin=[4,3]
    axis, /xaxis, xrange=interpol(xebv,x,!x.crange), xthick=postthick, $
      charthick=postthick, charsize=charsize, xtitle = 'E(B-V) [mag]', xsty=1

    legend, '(a)', /left, /top, box=0, charsize=charsize, charthick=postthick

    label = ['No Reddening Correction']
    legend, textoidl(label), /left, /bottom, box=0, charsize=charsize, $
      charthick=postthick

;;; ############################################################
;;; Panel 2: A(Ha) = 1
;;; ############################################################
;;
;;    A3727 = 1.0*k_lambda(3727.0,/odonnell)/k_lambda(6563.0,/odonnell)
;;    A3727_err = 0.1*k_lambda(3727.0,/odonnell)/k_lambda(6563.0,/odonnell)
;;    
;;    f3727 = atlasdust[indx].oii_3727_continuum[0]*10^(0.4*A3727)
;;    f3727_err = atlasdust[indx].oii_3727_continuum[1]*10^(0.4*A3727) + $
;;      atlasdust[indx].oii_3727_continuum[0]*alog(10.0)*10^(0.4*A3727)*A3727_err
;;
;;    y = alog10(f3727/ha)
;;    yerr = im_compute_error(f3727,f3727_err,ha,ha_err,/log)
;;
;;    if keyword_set(nfgs) then begin
;;
;;       f3727nfgs = nfgsdust[indxnfgs].oii_3727_continuum[0]*10^(0.4*A3727)
;;       f3727nfgs_err = nfgsdust[indxnfgs].oii_3727_continuum[1]*10^(0.4*A3727) + $
;;         nfgsdust[indxnfgs].oii_3727_continuum[0]*alog(10.0)*10^(0.4*A3727)*A3727_err
;;
;;       ynfgs = alog10(f3727nfgs/hanfgs)
;;       yerrnfgs = im_compute_error(f3727nfgs,f3727nfgs_err,hanfgs,hanfgs_err,/log)
;;       
;;    endif
;;
;;    stats = im_stats(y,/verbose)
;;
;;    atlas1d_lineplot, x, y, xerr, yerr, plottype=3, postscript=postscript, $
;;      xtitle='', ytitle=ytitle, xrange=xrange, yrange=yrange, legendtype=0, $
;;      charsize=charsize, xtickname=replicate(' ',10), $
;;      position=pos[*,1], /noerase, xnfgs=xnfgs, ynfgs=ynfgs, xerrnfgs=xerrnfgs, yerrnfgs=yerrnfgs
;;
;;    legend, '(b)', /left, /top, box=0, charsize=charsize, charthick=postthick
;;
;;    label = ['A(\lambda3727) = 1.85 \pm 0.2']
;;    legend, textoidl(label), /left, /bottom, box=0, charsize=charsize, $
;;      charthick=postthick
    
; ############################################################
; Panel 2: Good Reddening Correction
; ############################################################

    f3727 = atlasnodust[indx].oii_3727_continuum[0]
    f3727_err = atlasnodust[indx].oii_3727_continuum[1]

    y = alog10(f3727/ha)
    yerr = im_compute_error(f3727,f3727_err,ha,ha_err,/log)

    if keyword_set(nfgs) then begin

       f3727nfgs = nfgsnodust[indxnfgs].oii_3727_continuum[0]
       f3727nfgs_err = nfgsnodust[indxnfgs].oii_3727_continuum[1]

       ynfgs = alog10(f3727nfgs/hanfgs)
       yerrnfgs = im_compute_error(f3727nfgs,f3727nfgs_err,hanfgs,hanfgs_err,/log)
       
    endif

    stats = im_stats(y,/verbose,/no_head)

    ytitle = 'log (F_{\lambda3727}_{cor}/H\alpha_{cor})'

    atlas1d_lineplot, x, y, xerr, yerr, plottype=3, postscript=postscript, $
      xtitle=xtitle, ytitle=ytitle, xrange=xrange, yrange=yrange, legendtype=0, $
      charsize=charsize, position=pos[*,1], /noerase, xnfgs=xnfgs, ynfgs=ynfgs, $
      xerrnfgs=xerrnfgs, yerrnfgs=yerrnfgs

    legend, '(c)', /left, /top, box=0, charsize=charsize, charthick=postthick

    label = ['Individual Reddening Correction']
    legend, textoidl(label), /left, /bottom, box=0, charsize=charsize, $
      charthick=postthick

    im_openclose, postscript=postscript, /close

; ------------------------------------------------------------
; E(Hb-Ha) vs f3727/Ha_obs
; ------------------------------------------------------------

    psname = 'ehbha_vs_f3727ha_obs'
    im_openclose, pspath+psname, postscript=postscript

    cut = where(atlasnodust.ehbha_err gt 0.0)
    lineratio, atlasdust[cut], '', '', 'OII_3727_CONTINUUM', 'H_ALPHA', $
      dum1, dum2, y, yerr, index=indx, nindex=nindx, snrcut=snrcut

    x = atlasnodust[cut[indx]].ehbha
    xerr = atlasnodust[cut[indx]].ehbha_err
    xebv = x/(k_lambda(4861,/odonnell)-k_lambda(6563,/odonnell))

    ybig = y
    xbig = x
    
    if keyword_set(kenn92) then begin
       cutkenn92 = where(kenn92nodust.ehbha_err gt 0.0)
       lineratio, kenn92dust[cutkenn92], '', '', 'OII_3727_CONTINUUM', 'H_ALPHA', $
         dum1, dum2, ykenn92, yerrkenn92, index=indxkenn92, snrcut=snrcut
       xkenn92 = kenn92nodust[cutkenn92[indxkenn92]].ehbha
       xerrkenn92 = kenn92nodust[cutkenn92[indxkenn92]].ehbha_err
    endif

    if keyword_set(nfgs) then begin
       cutnfgs = where(nfgsnodust.ehbha_err gt 0.0)
       lineratio, nfgsdust[cutnfgs], '', '', 'OII_3727_CONTINUUM', 'H_ALPHA', $
         dum1, dum2, ynfgs, yerrnfgs, index=indxnfgs, snrcut=snrcut, $
      xsyserr=nfgssyserr, ysyserr=nfgssyserr
       xnfgs = nfgsnodust[cutnfgs[indxnfgs]].ehbha
       xerrnfgs = nfgsnodust[cutnfgs[indxnfgs]].ehbha_err

       ybig = [y,ynfgs]
       xbig = [x,xnfgs]
    endif

    xtitle = 'E(H\beta-H\alpha) [mag]'
    ytitle = 'log (F_{\lambda3727}/H\alpha)_{obs}'

    xrange = ehbharange
    yrange = f3727harange

    atlas1d_lineplot, x, y, xerr, yerr, plottype=3, postscript=postscript, $
      xtitle=xtitle, ytitle=ytitle, xrange=xrange, yrange=yrange, legendtype=0, $
      xstyle=11, ymargin=[4,3], /right, /top, $
;     xkenn92=xkenn92, ykenn92=ykenn92, xerrkenn92=xerrkenn92, yerrkenn92=yerrkenn92, $
      xnfgs=xnfgs, ynfgs=ynfgs, xerrnfgs=xerrnfgs, yerrnfgs=yerrnfgs
    axis, /xaxis, xrange=interpol(xebv,x,!x.crange), xthick=postthick, $
      charthick=postthick, charsize=2.0, xtitle = 'E(B-V) [mag]', xsty=1

    im_openclose, postscript=postscript, /close    
    
; ------------------------------------------------------------
; E(Hb-Ha) vs f3727/Ha_cor
; ------------------------------------------------------------

    psname = 'ehbha_vs_f3727ha_cor'
    im_openclose, pspath+psname, postscript=postscript

    lineratio, atlasnodust, '', '', 'OII_3727_CONTINUUM', 'H_ALPHA', $
      dum1, dum2, y, yerr, index=indx, nindex=nindx, snrcut=snrcut

    x = atlasnodust[indx].ehbha
    xerr = atlasnodust[indx].ehbha_err
    xebv = x/(k_lambda(4861,/odonnell)-k_lambda(6563,/odonnell))

    ybig = y
    xbig = x
    
    if keyword_set(kenn92) then begin
       lineratio, kenn92nodust, '', '', 'OII_3727_CONTINUUM', 'H_ALPHA', $
         dum1, dum2, ykenn92, yerrkenn92, index=indxkenn92, snrcut=snrcut
       xkenn92 = kenn92nodust[indxkenn92].ehbha
       xerrkenn92 = kenn92nodust[indxkenn92].ehbha_err
    endif

    if keyword_set(nfgs) then begin
       lineratio, nfgsnodust, '', '', 'OII_3727_CONTINUUM', 'H_ALPHA', $
         dum1, dum2, ynfgs, yerrnfgs, index=indxnfgs, snrcut=snrcut, $
      xsyserr=nfgssyserr, ysyserr=nfgssyserr
       xnfgs = nfgsnodust[indxnfgs].ehbha
       xerrnfgs = nfgsnodust[indxnfgs].ehbha_err

       ybig = [y,ynfgs]
       xbig = [x,xnfgs]
    endif

    xtitle = 'E(H\beta-H\alpha) [mag]'
    ytitle = 'log (F_{\lambda3727}/H\alpha)_{cor}'

    xrange = ehbharange
    yrange = f3727harange

    atlas1d_lineplot, x, y, xerr, yerr, plottype=3, postscript=postscript, $
      xtitle=xtitle, ytitle=ytitle, xrange=xrange, yrange=yrange, legendtype=0, $
      xstyle=11, ymargin=[4,3], /right, /top, $
;     xkenn92=xkenn92, ykenn92=ykenn92, xerrkenn92=xerrkenn92, yerrkenn92=yerrkenn92, $
      xnfgs=xnfgs, ynfgs=ynfgs, xerrnfgs=xerrnfgs, yerrnfgs=yerrnfgs
    axis, /xaxis, xrange=interpol(xebv,x,!x.crange), xthick=postthick, $
      charthick=postthick, charsize=2.0, xtitle = 'E(B-V) [mag]', xsty=1

    im_openclose, postscript=postscript, /close    
    
; ------------------------------------------------------------
; E(Hb-Ha) vs f3727_obs/Ha_cor
; ------------------------------------------------------------

    psname = 'ehbha_vs_f3727obs_hacor'
    im_openclose, pspath+psname, postscript=postscript

    indx = where((atlasnodust.ehbha_err gt 0.0) and $
      (atlasdust.oii_3727_continuum[0]/atlasdust.oii_3727_continuum[1] gt snrcut) and $
      (atlasnodust.h_alpha[0]/atlasnodust.h_alpha[1] gt snrcut),nindx)

    x = atlasnodust[indx].ehbha
    xerr = atlasnodust[indx].ehbha_err
    xebv = x/(k_lambda(4861,/odonnell)-k_lambda(6563,/odonnell))

    y = alog10(atlasdust[indx].oii_3727_continuum[0]/atlasnodust[indx].h_alpha[0])
    yerr = im_compute_error(atlasdust[indx].oii_3727_continuum[0],atlasdust[indx].oii_3727_continuum[1],$
      atlasnodust[indx].h_alpha[0],atlasnodust[indx].h_alpha[1],/log)
    
    if keyword_set(nfgs) then begin
       indxnfgs = where((nfgsnodust.ehbha_err gt 0.0) and $
         (nfgsdust.oii_3727_continuum[0]/nfgsdust.oii_3727_continuum[1] gt snrcut) and $
         (nfgsnodust.h_alpha[0]/nfgsnodust.h_alpha[1] gt snrcut),nindxnfgs)

       xnfgs = nfgsnodust[indxnfgs].ehbha
       xnfgserr = nfgsnodust[indxnfgs].ehbha_err

       ynfgs = alog10(nfgsdust[indxnfgs].oii_3727_continuum[0]/nfgsnodust[indxnfgs].h_alpha[0])
       yerrnfgs = im_compute_error(nfgsdust[indxnfgs].oii_3727_continuum[0],nfgsdust[indxnfgs].oii_3727_continuum[1],$
         nfgsnodust[indxnfgs].h_alpha[0],nfgsnodust[indxnfgs].h_alpha[1],/log)
    endif

    xtitle = 'E(H\beta-H\alpha) [mag]'
    ytitle = 'log (F_{\lambda3727 obs}/H\alpha_{cor})'

    xrange = ehbharange
    yrange = f3727hacorrange

    atlas1d_lineplot, x, y, xerr, yerr, plottype=3, postscript=postscript, $
      xtitle=xtitle, ytitle=ytitle, xrange=xrange, yrange=yrange, legendtype=1, $
      xstyle=11, ymargin=[4,3], /right, /top, xnfgs=xnfgs, ynfgs=ynfgs, xerrnfgs=xerrnfgs, $
      yerrnfgs=yerrnfgs
    axis, /xaxis, xrange=interpol(xebv,x,!x.crange), xthick=postthick, charthick=postthick, $
      charsize=2.0, xtitle = 'E(B-V) [mag]', xsty=1

    im_openclose, postscript=postscript, /close    
    
; ------------------------------------------------------------
; 12+log(O/H) [NIIHA] vs [O II]_cor/Ha_cor
; ------------------------------------------------------------

    psname = '12oh_NIIHA_vs_oiiha_cor'
    im_openclose, pspath+psname, postscript=postscript

    cut = where(atlasnodust.zstrong_12oh_niiha gt -900)
    lineratio, atlasnodust[cut], '', '', 'OII_3727', 'H_ALPHA', $
      dum1, dum2, y, yerr, index=indx, nindex=nindx, snrcut=snrcut

    x = atlasnodust[cut[indx]].zstrong_12oh_niiha
    xerr = atlasnodust[cut[indx]].zstrong_12oh_niiha_err

    xbig = x
    ybig = y

    if keyword_set(nfgs) then begin

       cutnfgs = where(nfgsnodust.zstrong_12oh_niiha gt -900)
       lineratio, nfgsnodust[cutnfgs], '', '', 'OII_3727', 'H_ALPHA', $
         dum1, dum2, ynfgs, yerrnfgs, index=indxnfgs, snrcut=snrcut

       xnfgs = nfgsnodust[cutnfgs[indxnfgs]].zstrong_12oh_niiha
       xerrnfgs = nfgsnodust[cutnfgs[indxnfgs]].zstrong_12oh_niiha_err

       xbig = [xbig,xnfgs]
       ybig = [ybig,ynfgs]

    endif    

    xtitle = '12 + log (O/H) [N_{2}]'
    ytitle = 'log ([O II]/H\alpha)_{cor}'

    xrange = ohrange
    yrange = oiiharange

; HII regions

    if keyword_set(hiiregions) then begin
       good = where((hii.zstrong_12oh_niiha gt -900.0) and (hii.oii_h_alpha gt -900.0))
       xregion = hii[good].zstrong_12oh_niiha & xerrregion = hii[good].zstrong_12oh_niiha_err
       yregion = hii[good].oii_h_alpha & yerrregion = hii[good].oii_h_alpha_err
    endif

    atlas1d_lineplot, x, y, xerr, yerr, plottype=3, postscript=postscript, $
      xtitle=xtitle, ytitle=ytitle, xrange=xrange, yrange=yrange, legendtype=0, $
      /right, /top, xregion=xregion, yregion=yregion, xerrregion=xerrregion, $
      yerrregion=yerrregion, xnfgs=xnfgs, ynfgs=ynfgs, xerrnfgs=xerrnfgs, yerrnfgs=yerrnfgs

    im_openclose, postscript=postscript, /close    

; ------------------------------------------------------------
; 3-panel H-alpha plot - Ha & Ha_cor
; ------------------------------------------------------------
    
    psname = 'haha_cor_ewoii_multipanel'
    im_openclose, pspath+psname, postscript=postscript

    indx = where((atlasdust.h_alpha[0]/atlasdust.h_alpha[1] gt snrcut) and $   
      (atlasdust.oii_3727_ew[0]/atlasdust.oii_3727_ew[1] gt snrcut) and $
      (atlasnodust.h_alpha[0]/atlasnodust.h_alpha[1] gt snrcut) and $
      (atlasnodust.ebv_hahb_err gt 0.0),nindx)

    x = alog10(atlasdust[indx].oii_3727_ew[0])
    xerr = atlasdust[indx].oii_3727_ew[1]/atlasdust[indx].oii_3727_ew[0]/alog(10)
    
    ha = atlasnodust[indx].h_alpha[0]
    ha_err = atlasnodust[indx].h_alpha[1]

    if keyword_set(nfgs) then begin

       indxnfgs = where((nfgsdust.h_alpha[0]/nfgsdust.h_alpha[1] gt snrcut) and $      
         (nfgsdust.oii_3727_ew[0]/nfgsdust.oii_3727_ew[1] gt snrcut) and $
         (nfgsnodust.h_alpha[0]/nfgsnodust.h_alpha[1] gt snrcut) and (nfgsnodust.ebv_hahb_err gt 0.0),nindx)

       xnfgs = alog10(nfgsnodust[indxnfgs].oii_3727_ew[0])
       xerrnfgs = nfgsnodust[indxnfgs].oii_3727_ew[1]/nfgsnodust[indxnfgs].oii_3727_ew[0]/alog(10)
       
       hanfgs = nfgsnodust[indxnfgs].h_alpha[0]
       hanfgs_err = nfgsnodust[indxnfgs].h_alpha[1]

    endif
    
    xrange = [0.3,2.3] ; ewoiirange
    yrange = hahacorrange

    xtitle = 'log EW([O II])'
    ytitle = 'log (H\alpha/H\alpha_{cor})'

    pagemaker, nx=1, ny=3, yspace=0, xspace=0, xmargin=[1.25,0.35], $
      ymargin=[0.3,1.1], position=pos, /normal

; ##########################################################
; Panel 1: No reddening correction
; ##########################################################

    yha = atlasdust[indx].h_alpha[0]
    yha_err = atlasdust[indx].h_alpha[1]
    
    y = alog10(yha/ha)
    yerr = im_compute_error(yha,yha_err,ha,ha_err,/log)

    if keyword_set(nfgs) then begin

       yhanfgs = nfgsdust[indxnfgs].h_alpha[0]
       yhanfgs_err = nfgsdust[indxnfgs].h_alpha[1]
       
       ynfgs = alog10(yhanfgs/hanfgs)
       yerrnfgs = im_compute_error(yhanfgs,yhanfgs_err,hanfgs,hanfgs_err,/log)
       
    endif

    stats = im_stats(y,/verbose)

    atlas1d_lineplot, x, y, xerr, yerr, plottype=3, postscript=postscript, $
      xtitle='', ytitle=ytitle, xrange=xrange, yrange=yrange, legendtype=0, $
      charsize=charsize, position=pos[*,0], xnfgs=xnfgs, ynfgs=ynfgs, xerrnfgs=xerrnfgs, $
      yerrnfgs=yerrnfgs, xtickname=replicate(' ',10)
    djs_oplot, !x.crange, [0,0], line=2, thick=postthick

    legend, '(a)', /left, /top, box=0, charsize=charsize, charthick=postthick

    label = ['A(H\alpha) = 0']
    legend, textoidl(label), /left, /bottom, box=0, charsize=charsize, $
      charthick=postthick

; ############################################################
; Panel 2: A(Ha) = 1
; ############################################################

    Aha = 1.0
    Aha_err = 0.1
    
    yha = atlasdust[indx].h_alpha[0]*10^(0.4*Aha)
    yha_err = atlasdust[indx].h_alpha[1]*10^(0.4*Aha) + $
      atlasdust[indx].h_alpha[0]*alog(10.0)*10^(0.4*Aha)*Aha_err

    y = alog10(yha/ha)
    yerr = im_compute_error(yha,yha_err,ha,ha_err,/log)

    if keyword_set(nfgs) then begin

       yhanfgs = nfgsdust[indxnfgs].h_alpha[0]*10^(0.4*Aha)
       yhanfgs_err = nfgsdust[indxnfgs].h_alpha[1]*10^(0.4*Aha) + $
         nfgsdust[indxnfgs].h_alpha[0]*alog(10.0)*10^(0.4*Aha)*Aha_err

       ynfgs = alog10(yhanfgs/hanfgs)
       yerrnfgs = im_compute_error(yhanfgs,yhanfgs_err,hanfgs,hanfgs_err,/log)
       
    endif

    stats = im_stats(y,/verbose)

    atlas1d_lineplot, x, y, xerr, yerr, plottype=3, postscript=postscript, $
      xtitle='', ytitle=ytitle, xrange=xrange, yrange=yrange, legendtype=0, $
      charsize=charsize, xtickname=replicate(' ',10), $
      position=pos[*,1], /noerase, xnfgs=xnfgs, ynfgs=ynfgs, xerrnfgs=xerrnfgs, yerrnfgs=yerrnfgs
    djs_oplot, !x.crange, [0,0], line=2, thick=postthick

    legend, '(b)', /left, /top, box=0, charsize=charsize, charthick=postthick

    label = ['A(H\alpha) = 1 \pm 0.1']
    legend, textoidl(label), /left, /bottom, box=0, charsize=charsize, $
      charthick=postthick
    
; ############################################################
; Panel 3: A(Ha) = 0.55
; ############################################################

    Aha = 0.55
    Aha_err = 0.45
    
    yha = atlasdust[indx].h_alpha[0]*10^(0.4*Aha)
    yha_err = atlasdust[indx].h_alpha[1]*10^(0.4*Aha) + $
      atlasdust[indx].h_alpha[0]*alog(10.0)*10^(0.4*Aha)*Aha_err

    y = alog10(yha/ha)
    yerr = im_compute_error(yha,yha_err,ha,ha_err,/log)

    if keyword_set(nfgs) then begin

       yhanfgs = nfgsdust[indxnfgs].h_alpha[0]*10^(0.4*Aha)
       yhanfgs_err = nfgsdust[indxnfgs].h_alpha[1]*10^(0.4*Aha) + $
         nfgsdust[indxnfgs].h_alpha[0]*alog(10.0)*10^(0.4*Aha)*Aha_err

       ynfgs = alog10(yhanfgs/hanfgs)
       yerrnfgs = im_compute_error(yhanfgs,yhanfgs_err,hanfgs,hanfgs_err,/log)
       
    endif

    stats = im_stats(y,/verbose)

    atlas1d_lineplot, x, y, xerr, yerr, plottype=3, postscript=postscript, $
      xtitle=xtitle, ytitle=ytitle, xrange=xrange, yrange=yrange, legendtype=0, $
      charsize=charsize, position=pos[*,2], /noerase, xnfgs=xnfgs, ynfgs=ynfgs, $
      xerrnfgs=xerrnfgs, yerrnfgs=yerrnfgs
    djs_oplot, !x.crange, [0,0], line=2, thick=postthick

    legend, '(b)', /left, /top, box=0, charsize=charsize, charthick=postthick

    label = ['A(H\alpha) = 0.55 \pm 0.45']
    legend, textoidl(label), /left, /bottom, box=0, charsize=charsize, $
      charthick=postthick

    im_openclose, postscript=postscript, /close

; ------------------------------------------------------------
; 4-panel H-beta plot - E(B-V) vs Hb/Ha
; ------------------------------------------------------------
    
    psname = 'ebv_vs_hbha_multipanel'
    im_openclose, pspath+psname, postscript=postscript

    indx = where((atlasdust.h_beta_uncor[0]/atlasdust.h_beta_uncor[1] gt snrcut) and $      
      (atlasnodust.h_alpha[0]/atlasnodust.h_alpha[1] gt snrcut) and $
      (atlasnodust.ebv_hahb_err gt 0.0)$
      ,nindx)

    x = atlasnodust[indx].ebv_hahb
    xerr = atlasnodust[indx].ebv_hahb_err
    
    ha = atlasnodust[indx].h_alpha[0]
    ha_err = atlasnodust[indx].h_alpha[1]

    if keyword_set(nfgs) then begin

       indxnfgs = where((nfgsdust.h_beta_uncor[0]/nfgsdust.h_beta_uncor[1] gt snrcut) and $      
         (nfgsdust.oii_3727_ew[0]/nfgsdust.oii_3727_ew[1] gt snrcut) and $
         (nfgsnodust.h_alpha[0]/nfgsnodust.h_alpha[1] gt snrcut) and (nfgsnodust.ebv_hahb_err gt 0.0),nindx)

       xnfgs = nfgsnodust[indxnfgs].ebv_hahb
       xerrnfgs = nfgsnodust[indxnfgs].ebv_hahb_err
       
       hanfgs = nfgsnodust[indxnfgs].h_alpha[0]
       hanfgs_err = nfgsnodust[indxnfgs].h_alpha[1]

    endif
    
    hahb_true = atlasdust[indx].h_alpha[0]/atlasdust[indx].h_beta[0]
    hahb_uncor = atlasdust[indx].h_alpha_uncor[0]/atlasdust[indx].h_beta_uncor[0]

    xrange = [-0.1,0.9]
    yrange = [-3.2,0.4]

    xtitle = 'E(B-V) [H\alpha/H\beta]'
    ytitle = 'log (H\beta/H\alpha_{cor})'

    pagemaker, nx=2, ny=2, yspace=0, xspace=0, xmargin=[1.25,0.35], $
      ymargin=[0.3,1.1], position=pos, /normal

; ##########################################################
; Panel 1: No absorption correction, no reddening correction
; ##########################################################

    hb = atlasdust[indx].h_beta_uncor[0]
    hb_err = atlasdust[indx].h_beta_uncor[1]
    
    y = alog10(hb/ha)
    yerr = im_compute_error(hb,hb_err,ha,ha_err,/log)

    if keyword_set(nfgs) then begin

       hbnfgs = nfgsdust[indxnfgs].h_beta_uncor[0]
       hbnfgs_err = nfgsdust[indxnfgs].h_beta_uncor[1]
       
       ynfgs = alog10(hbnfgs/hanfgs)
       yerrnfgs = im_compute_error(hbnfgs,hbnfgs_err,hanfgs,hanfgs_err,/log)
       
    endif

    stats = im_stats(y,/verbose)

    atlas1d_lineplot, x, y, xerr, yerr, plottype=3, postscript=postscript, $
      xtitle='', ytitle=ytitle, xrange=xrange, yrange=yrange, legendtype=0, $
      charsize=charsize, xtickname=replicate(' ',10), position=pos[*,0], $
      xnfgs=xnfgs, ynfgs=ynfgs, xerrnfgs=xerrnfgs, yerrnfgs=yerrnfgs
    djs_oplot, !x.crange, alog10(1.0/HaHb)*[1,1], line=2, thick=postthick

    legend, '(a)', /left, /top, box=0, charsize=charsize, charthick=postthick
    
    label = ['No Absorption','No Reddening']
    legend, textoidl(label), /left, /bottom, box=0, charsize=charsize, $
      spacing=1.7, charthick=postthick

; ############################################################
; Panel 2: Mean absorption correction, no reddening correction
; ############################################################

    hb_ew_abs     = 3.5 ; EW [mean absorption]
    hb_ew_abs_err = 1.0

    hb = atlasdust[indx].h_beta_uncor[0] + hb_ew_abs*atlasdust[indx].babs_h_beta_continuum[0]
    hb_err = sqrt(atlasdust[indx].h_beta[1]^2 + (hb_ew_abs_err*atlasdust[indx].babs_h_beta_continuum[0])^2)
    
    y = alog10(hb/ha)
    yerr = im_compute_error(hb,hb_err,ha,ha_err,/log)

    if keyword_set(nfgs) then begin

       hbnfgs = nfgsdust[indxnfgs].h_beta_uncor[0] + hb_ew_abs*nfgsdust[indxnfgs].babs_h_beta_continuum[0]
       hbnfgs_err = nfgsdust[indxnfgs].h_beta[1]
       
       ynfgs = alog10(hbnfgs/hanfgs)
       yerrnfgs = im_compute_error(hbnfgs,hbnfgs_err,hanfgs,hanfgs_err,/log)

    endif

    stats = im_stats(y,/verbose,/no_head)

    atlas1d_lineplot, x, y, xerr, yerr, plottype=3, postscript=postscript, $
      xtitle='', ytitle='', xrange=xrange, yrange=yrange, legendtype=0, $
      charsize=charsize, xtickname=replicate(' ',10), ytickname=replicate(' ',10), $
      position=pos[*,1], /noerase, xnfgs=xnfgs, ynfgs=ynfgs, xerrnfgs=xerrnfgs, yerrnfgs=yerrnfgs
    djs_oplot, !x.crange, alog10(1.0/HaHb)*[1,1], line=2, thick=postthick

    legend, '(b)', /left, /top, box=0, charsize=charsize, charthick=postthick

    label = ['Mean Absorption','No Reddening']
    legend, textoidl(label), /left, /bottom, box=0, charsize=charsize, $
      spacing=1.7, charthick=postthick
    
; ############################################################
; Panel 3: Good absorption correction, no reddening correction
; ############################################################

    hb = atlasdust[indx].h_beta[0]
    hb_err = atlasdust[indx].h_beta[1]
    
    y = alog10(hb/ha)
    yerr = im_compute_error(hb,hb_err,ha,ha_err,/log)

    if keyword_set(nfgs) then begin

       hbnfgs = nfgsdust[indxnfgs].h_beta[0]
       hbnfgs_err = nfgsdust[indxnfgs].h_beta[1]
    
       ynfgs = alog10(hbnfgs/hanfgs)
       yerrnfgs = im_compute_error(hbnfgs,hbnfgs_err,hanfgs,hanfgs_err,/log)

    endif

    stats = im_stats(y,/verbose,/no_head)
    
    atlas1d_lineplot, x, y, xerr, yerr, plottype=3, postscript=postscript, $
      xtitle=xtitle, ytitle=ytitle, xrange=xrange, yrange=yrange, legendtype=0, $
      charsize=charsize, position=pos[*,2], /noerase, xnfgs=xnfgs, ynfgs=ynfgs, $
      xerrnfgs=xerrnfgs, yerrnfgs=yerrnfgs
    djs_oplot, !x.crange, alog10(1.0/HaHb)*[1,1], line=2, thick=postthick

    legend, '(c)', /left, /top, box=0, charsize=charsize, charthick=postthick

    label = ['Individual Absorption','No Reddening']
    legend, textoidl(label), /left, /bottom, box=0, charsize=charsize, $
      spacing=1.7, charthick=postthick
    
; ###############################################################
; Panel 4: Good absorption correction, Good reddening correction
; ###############################################################

    hb = atlasnodust[indx].h_beta[0]
    hb_err = atlasnodust[indx].h_beta[1]
    
    y = alog10(hb/ha)
    yerr = im_compute_error(hb,hb_err,ha,ha_err,/log)

    if keyword_set(nfgs) then begin

       hbnfgs = nfgsnodust[indxnfgs].h_beta[0]
       hbnfgs_err = nfgsnodust[indxnfgs].h_beta[1]

       ynfgs = alog10(hbnfgs/hanfgs)
       yerrnfgs = im_compute_error(hbnfgs,hbnfgs_err,hanfgs,hanfgs_err,/log)

    endif
       
    stats = im_stats(y,/verbose,/no_head)

    atlas1d_lineplot, x, y, xerr, yerr, plottype=3, postscript=postscript, $
      xtitle=xtitle, ytitle='', xrange=xrange, yrange=yrange, legendtype=0, $
      charsize=charsize, ytickname=replicate(' ',10), position=pos[*,3], /noerase, $
      xnfgs=xnfgs, ynfgs=ynfgs, xerrnfgs=xerrnfgs, yerrnfgs=yerrnfgs      
    djs_oplot, !x.crange, alog10(1.0/HaHb)*[1,1], line=2, thick=postthick

    legend, '(d)', /left, /top, box=0, charsize=charsize, charthick=postthick

    label = ['Individual Absorption','Individual Reddening']
    legend, textoidl(label), /left, /bottom, box=0, charsize=charsize, $
      spacing=1.7, charthick=postthick

    im_openclose, postscript=postscript, /close

; ------------------------------------------------------------
; [N II]/Ha vs Hb_obs/Ha_cor 
; ------------------------------------------------------------
    
    psname = 'niiha_vs_hbobs_hacor'
    im_openclose, pspath+psname, postscript=postscript

    cut = where((atlasdust.h_beta[0]/atlasdust.h_beta[1] gt snrcut))
    lineratio, atlasnodust[cut], 'NII_6584', 'H_ALPHA', '', '', $
      x, xerr, dum1, dum2, index=indx, nindex=nindx, snrcut=snrcut

    ha = atlasnodust[cut[indx]].h_alpha[0]
    ha_err = atlasnodust[cut[indx]].h_alpha[1]

    hb = atlasdust[cut[indx]].h_beta[0]
    hb_err = atlasdust[cut[indx]].h_beta[1]

    y = alog10(hb/ha)
    yerr = im_compute_error(hb,hb_err,ha,ha_err,/log)

    if keyword_set(nfgs) then begin

       cutnfgs = where(nfgsdust.h_beta[0]/nfgsdust.h_beta[1] gt snrcut)
       lineratio, nfgsnodust[cutnfgs], 'NII_6584', 'H_ALPHA', '', '', $
         xnfgs, xerrnfgs, dum1, dum2, index=indxnfgs, nindex=nindxnfgs, snrcut=snrcut

       hanfgs = nfgsnodust[cutnfgs[indxnfgs]].h_alpha[0]
       hanfgs_err = nfgsnodust[cutnfgs[indxnfgs]].h_alpha[1]

       hbnfgs = nfgsdust[cutnfgs[indxnfgs]].h_beta[0]
       hbnfgs_err = nfgsdust[cutnfgs[indxnfgs]].h_beta[1]

       ynfgs = alog10(hbnfgs/hanfgs)
       yerrnfgs = im_compute_error(hbnfgs,hbnfgs_err,hanfgs,hanfgs_err,/log)

    endif
    
    xrange = niiharange
    yrange = hbhacorrange

    xtitle = 'log ([N II] \lambda6584/H\alpha)_{cor}'
    ytitle = 'log (H\beta_{obs}/H\alpha_{cor})'

    atlas1d_lineplot, x, y, xerr, yerr, plottype=3, postscript=postscript, $
      xtitle=xtitle, ytitle=ytitle, xrange=xrange, yrange=yrange, legendtype=1, $
      /left, /bottom, xnfgs=xnfgs, ynfgs=ynfgs, xerrnfgs=xerrnfgs, yerrnfgs=yerrnfgs
    djs_oplot, !x.crange, alog10(1.0/HaHb)*[1,1], line=2, thick=postthick
    
    im_openclose, postscript=postscript, /close

; ------------------------------------------------------------
; [O III]/[O II] vs Hb_obs/Ha_cor 
; ------------------------------------------------------------
    
    psname = 'oiiioii_vs_hbobs_hacor'
    im_openclose, pspath+psname, postscript=postscript

    cut = where((atlasdust.h_beta[0]/atlasdust.h_beta[1] gt snrcut))
    lineratio, atlasnodust[cut], 'OIII_5007', 'OII_3727', '', '', $
      x, xerr, dum1, dum2, index=indx, nindex=nindx, snrcut=snrcut

    ha = atlasnodust[cut[indx]].h_alpha[0]
    ha_err = atlasnodust[cut[indx]].h_alpha[1]

    hb = atlasdust[cut[indx]].h_beta[0]
    hb_err = atlasdust[cut[indx]].h_beta[1]

    y = alog10(hb/ha)
    yerr = im_compute_error(hb,hb_err,ha,ha_err,/log)

    if keyword_set(nfgs) then begin

       cutnfgs = where(nfgsdust.h_beta[0]/nfgsdust.h_beta[1] gt snrcut)
       lineratio, nfgsnodust[cutnfgs], 'OIII_5007', 'OII_3727', '', '', $
         xnfgs, xerrnfgs, dum1, dum2, index=indxnfgs, nindex=nindxnfgs, snrcut=snrcut

       hanfgs = nfgsnodust[cutnfgs[indxnfgs]].h_alpha[0]
       hanfgs_err = nfgsnodust[cutnfgs[indxnfgs]].h_alpha[1]

       hbnfgs = nfgsdust[cutnfgs[indxnfgs]].h_beta[0]
       hbnfgs_err = nfgsdust[cutnfgs[indxnfgs]].h_beta[1]

       ynfgs = alog10(hbnfgs/hanfgs)
       yerrnfgs = im_compute_error(hbnfgs,hbnfgs_err,hanfgs,hanfgs_err,/log)

    endif
    
    xrange = oiiioiirange
    yrange = hbhacorrange

    xtitle = 'log ([O III] \lambda5007/[O II])_{cor}'
    ytitle = 'log (H\beta_{obs}/H\alpha_{cor})'

    atlas1d_lineplot, x, y, xerr, yerr, plottype=3, postscript=postscript, $
      xtitle=xtitle, ytitle=ytitle, xrange=xrange, yrange=yrange, legendtype=1, $
      /left, /bottom, xnfgs=xnfgs, ynfgs=ynfgs, xerrnfgs=xerrnfgs, yerrnfgs=yerrnfgs
    djs_oplot, !x.crange, alog10(1.0/HaHb)*[1,1], line=2, thick=postthick
    
    im_openclose, postscript=postscript, /close

; ------------------------------------------------------------
; R23 vs Hb_obs/Ha_cor 
; ------------------------------------------------------------
    
    psname = 'R23_vs_hbobs_hacor'
    im_openclose, pspath+psname, postscript=postscript

    cut = where((atlasdust.h_beta[0]/atlasdust.h_beta[1] gt snrcut))
    lineratio, atlasnodust[cut], ['OIII_4959','OIII_5007', 'OII_3727'], 'H_BETA', '', '', $
      x, xerr, dum1, dum2, index=indx, nindex=nindx, snrcut=snrcut

    ha = atlasnodust[cut[indx]].h_alpha[0]
    ha_err = atlasnodust[cut[indx]].h_alpha[1]

    hb = atlasdust[cut[indx]].h_beta[0]
    hb_err = atlasdust[cut[indx]].h_beta[1]

    y = alog10(hb/ha)
    yerr = im_compute_error(hb,hb_err,ha,ha_err,/log)

    if keyword_set(nfgs) then begin

       cutnfgs = where(nfgsdust.h_beta[0]/nfgsdust.h_beta[1] gt snrcut)
       lineratio, nfgsnodust[cutnfgs], ['OIII_4959','OIII_5007', 'OII_3727'], 'H_BETA', '', '', $
         xnfgs, xerrnfgs, dum1, dum2, index=indxnfgs, nindex=nindxnfgs, snrcut=snrcut

       hanfgs = nfgsnodust[cutnfgs[indxnfgs]].h_alpha[0]
       hanfgs_err = nfgsnodust[cutnfgs[indxnfgs]].h_alpha[1]

       hbnfgs = nfgsdust[cutnfgs[indxnfgs]].h_beta[0]
       hbnfgs_err = nfgsdust[cutnfgs[indxnfgs]].h_beta[1]

       ynfgs = alog10(hbnfgs/hanfgs)
       yerrnfgs = im_compute_error(hbnfgs,hbnfgs_err,hanfgs,hanfgs_err,/log)

    endif
    
    xtitle = 'log R_{23}'
    ytitle = 'log (H\beta_{obs}/H\alpha_{cor})'

    xrange = R23range
    yrange = hbhacorrange

    atlas1d_lineplot, x, y, xerr, yerr, plottype=3, postscript=postscript, $
      xtitle=xtitle, ytitle=ytitle, xrange=xrange, yrange=yrange, legendtype=1, $
      /left, /bottom, xnfgs=xnfgs, ynfgs=ynfgs, xerrnfgs=xerrnfgs, yerrnfgs=yerrnfgs
    djs_oplot, !x.crange, alog10(1.0/HaHb)*[1,1], line=2, thick=postthick
    
    im_openclose, postscript=postscript, /close

; ------------------------------------------------------------
; L(B) distribution
; ------------------------------------------------------------

    psname = 'histogram_LB'
    im_openclose, pspath+psname, postscript=postscript

    indx = where(atlasdust.b_lum_obs gt -900.0)
    x = atlasdust[indx].b_lum_obs
    xabs = atlasdust[indx].m_b_obs

    if keyword_set(nfgs) then begin
       indxnfgs = where(nfgsdust.b_lum_obs gt -900.0)
       xnfgs = nfgsdust[indxnfgs].b_lum_obs
    endif

    xtitle = 'log L(B) [L(B)'+sunsymbol()+']'

    binsize = 0.25

    plothist, x, bin=binsize, xbin, ybin, /noplot, /halfbin
    if keyword_set(nfgs) then plothist, xnfgs, bin=binsize, xbinnfgs, ybinnfgs, /noplot, /halfbin

    xrange = LBrange
    yrange = minmax(ybin)*[1.0,1.3]

    djs_plot, [0], [0], /nodata, xthick=postthick, ythick=postthick, charsize=2.0, $
      charthick=postthick, thick=postthick, xtitle=textoidl(xtitle), ytitle='Number', $
      xstyle=11, ystyle=1, xrange=xrange, yrange=yrange
    axis, /xaxis, xrange=interpol(xabs,x,!x.crange), xthick=postthick, charthick=postthick, $
      charsize=2.0, xtitle = textoidl('M_{B} [mag]'), xsty=1

    plothist, x, bin=binsize, /overplot, color=djs_icolor('blue'), thick=postthick, /halfbin
    if keyword_set(nfgs) then plothist, xnfgs, bin=binsize, thick=postthick, line=2, /halfbin, $
      /overplot, /fill, /fline, forientation=45, fcolor=djs_icolor('red'), $
      color=djs_icolor('red')

    legend, ['K/M Atlas','NFGS'], /right, /top, charthick=postthick, $
      charsize=2.0, box=0, line=[0,2], thick=postthick, clear=postscript, $
      color=djs_icolor(['blue','red'])
    legend, '(a)', /left, /top, charthick=postthick, charsize=2.0, box=0, clear=postscript

    im_openclose, postscript=postscript, /close    

; ------------------------------------------------------------
; B-V distribution
; ------------------------------------------------------------

    psname = 'histogram_BV'
    im_openclose, pspath+psname, postscript=postscript

    indx = where((atlasdust.B gt -900.0) and (atlasdust.V gt -900.0),nindx)
    x = atlasdust[indx].B-atlasdust[indx].V
    
    if keyword_set(nfgs) then begin
       indxnfgs = where((nfgsdust.B gt -900.0) and (nfgsdust.V gt -900.0))
       xnfgs = nfgsdust[indxnfgs].B-nfgsdust[indxnfgs].V
    endif

    xtitle = 'B-V'

    binsize = 0.05
    
    plothist, x, bin=binsize, xbin, ybin, /noplot, /halfbin
    if keyword_set(nfgs) then plothist, xnfgs, bin=binsize, xbinnfgs, ybinnfgs, /noplot, /halfbin

    xrange = BVrange
    yrange = minmax(ybin)*[1.0,1.1]

    djs_plot, [0], [0], /nodata, xthick=postthick, ythick=postthick, charsize=2.0, $
      charthick=postthick, thick=postthick, xtitle=textoidl(xtitle), ytitle='Number', $
      xstyle=1, ystyle=1, xrange=xrange, yrange=yrange
    plothist, x, bin=binsize, /overplot, color=djs_icolor('blue'), thick=postthick, /halfbin
    if keyword_set(nfgs) then plothist, xnfgs, bin=binsize, thick=postthick, line=2, /halfbin, $
      /overplot, /fill, /fline, forientation=45, fcolor=djs_icolor('red'), $
      color=djs_icolor('red')

    legend, '(b)', /left, /top, charthick=postthick, charsize=2.0, box=0, clear=postscript

    im_openclose, postscript=postscript, /close    

; ------------------------------------------------------------
; L(IR) distribution
; ------------------------------------------------------------

    psname = 'histogram_L_IR'
    im_openclose, pspath+psname, postscript=postscript

    indx = where(atlasdust.ir_lum gt -900.0)
    x = atlasdust[indx].ir_lum
    
    if keyword_set(nfgs) then begin
       indxnfgs = where(nfgsdust.ir_lum gt -900.0)
       xnfgs = nfgsdust[indxnfgs].ir_lum
    endif

    xtitle = 'log L_{IR} [L'+sunsymbol()+']'

    binsize = 0.3
    
    plothist, x, bin=binsize, xbin, ybin, /noplot, /halfbin
    if keyword_set(nfgs) then plothist, xnfgs, bin=binsize, xbinnfgs, ybinnfgs, /noplot, /halfbin

    xrange = LIRrange
    yrange = minmax(ybin)*[1.0,1.1]

    djs_plot, [0], [0], /nodata, xthick=postthick, ythick=postthick, charsize=2.0, $
      charthick=postthick, thick=postthick, xtitle=textoidl(xtitle), ytitle='Number', $
      xstyle=1, ystyle=1, xrange=xrange, yrange=yrange
    plothist, x, bin=binsize, /overplot, color=djs_icolor('blue'), thick=postthick, /halfbin
    if keyword_set(nfgs) then plothist, xnfgs, bin=binsize, thick=postthick, line=2, /halfbin, $
      /overplot, /fill, /fline, forientation=45, fcolor=djs_icolor('red'), $
      color=djs_icolor('red')

    legend, '(c)', /left, /top, charthick=postthick, charsize=2.0, box=0, clear=postscript

    im_openclose, postscript=postscript, /close    

; ------------------------------------------------------------
; E(B-V) distribution
; ------------------------------------------------------------

    psname = 'histogram_ebv'
    im_openclose, pspath+psname, postscript=postscript

    indx = where(atlasnodust.ebv_hahb_err gt 0.0)
    x = atlasnodust[indx].ebv_hahb
    
    if keyword_set(nfgs) then begin
       indxnfgs = where(nfgsnodust.ebv_hahb_err gt 0.0)
       xnfgs = nfgsnodust[indxnfgs].ebv_hahb
    endif

    xtitle = 'E(B-V) [mag]'

    binsize = 0.05
    
    plothist, x, bin=binsize, xbin, ybin, /noplot, /halfbin
    if keyword_set(nfgs) then plothist, xnfgs, bin=binsize, xbinnfgs, ybinnfgs, /noplot, /halfbin

    xrange = ebvrange
    yrange = minmax(ybin)*[1.0,1.15]

    djs_plot, [0], [0], /nodata, xthick=postthick, ythick=postthick, charsize=2.0, $
      charthick=postthick, thick=postthick, xtitle=textoidl(xtitle), ytitle='Number', $
      xstyle=1, ystyle=1, xrange=xrange, yrange=yrange
    plothist, x, bin=binsize, /overplot, color=djs_icolor('blue'), thick=postthick, /halfbin
    if keyword_set(nfgs) then plothist, xnfgs, bin=binsize, thick=postthick, line=2, /halfbin, $
      /overplot, /fill, /fline, forientation=45, fcolor=djs_icolor('red'), $
      color=djs_icolor('red')

    legend, '(d)', /left, /top, charthick=postthick, charsize=2.0, box=0, clear=postscript

    im_openclose, postscript=postscript, /close    

; ------------------------------------------------------------
; L(IR)/L(B) distribution
; ------------------------------------------------------------

    psname = 'histogram_LIR_LB'
    im_openclose, pspath+psname, postscript=postscript

    indx = where(atlasdust.L_IR_L_B gt -900.0)
    x = alog10(atlasdust[indx].L_IR_L_B)
    
    if keyword_set(nfgs) then begin
       indxnfgs = where(nfgsdust.L_IR_L_B gt -900.0)
       xnfgs = alog10(nfgsdust[indxnfgs].L_IR_L_B)
    endif

    xtitle = 'log (L_{IR}/L_{B})'

    binsize = 0.2
    
    plothist, x, bin=binsize, xbin, ybin, /noplot, /halfbin
    if keyword_set(nfgs) then plothist, xnfgs, bin=binsize, xbinnfgs, ybinnfgs, /noplot, /halfbin

    xrange = iroptrange
    yrange = minmax(ybin)*[1.0,1.1]

    djs_plot, [0], [0], /nodata, xthick=postthick, ythick=postthick, charsize=2.0, $
      charthick=postthick, thick=postthick, xtitle=textoidl(xtitle), ytitle='Number', $
      xstyle=1, ystyle=1, xrange=xrange, yrange=yrange
    plothist, x, bin=binsize, /overplot, color=djs_icolor('blue'), thick=postthick, /halfbin
    if keyword_set(nfgs) then plothist, xnfgs, bin=binsize, thick=postthick, line=2, /halfbin, $
      /overplot, /fill, /fline, forientation=45, fcolor=djs_icolor('red'), $
      color=djs_icolor('red')

;   legend, '(b)', /left, /top, charthick=postthick, charsize=2.0, box=0, clear=postscript
    
    im_openclose, postscript=postscript, /close    

; ------------------------------------------------------------
; O32 distribution
; ------------------------------------------------------------

    psname = 'histogram_O32'
    im_openclose, pspath+psname, postscript=postscript

    lineratio, atlasnodust, ['OIII_4959','OIII_5007'],'OII_3727', '', '', $
      x, xerr, dum1, dum2, index=indx, nindex=nindx, snrcut=snrcut
    
    if keyword_set(nfgs) then begin
       lineratio, nfgsnodust, ['OIII_4959','OIII_5007'],'OII_3727', '', '', $
         xnfgs, xerrnfgs, dum1, dum2, index=indxnfgs, nindex=nindxnfgs, snrcut=snrcut
    endif

    xtitle = 'log O_{32}'

    binsize = 0.15
    
    plothist, x, bin=binsize, xbin, ybin, /noplot, /halfbin
    if keyword_set(nfgs) then plothist, xnfgs, bin=binsize, xbinnfgs, ybinnfgs, /noplot, /halfbin

    xrange = O32range
    yrange = minmax(ybin)*[1.0,1.1]

    djs_plot, [0], [0], /nodata, xthick=postthick, ythick=postthick, charsize=2.0, $
      charthick=postthick, thick=postthick, xtitle=textoidl(xtitle), ytitle='Number', $
      xstyle=1, ystyle=1, xrange=xrange, yrange=yrange
    plothist, x, bin=binsize, /overplot, color=djs_icolor('blue'), thick=postthick, /halfbin

    if keyword_set(nfgs) then plothist, xnfgs, bin=binsize, thick=postthick, line=2, /halfbin, $
      /overplot, /fill, /fline, forientation=45, fcolor=djs_icolor('red'), $
      color=djs_icolor('red')

;   if keyword_set(hiiregions) then begin
;      good = where(regions.oiii_5007_oii gt -900.0)
;      xregion = regions[good].oiii_5007_oii
;      plothist, xregion, xbinregion, ybinregion, bin=0.15, /noplot, /halfbin
;      ybinregion = round(max(ybin)*ybinregion/float(max(ybinregion)))
;      
;      oplot, xbinregion, ybinregion, ps=10, line=0
;      polyfill, xbinregion, ybinregion, /line_fill, linestyle=1, ps=10
;
;   endif
    
    im_openclose, postscript=postscript, /close    
           
; ------------------------------------------------------------
; R23 distribution
; ------------------------------------------------------------

    psname = 'histogram_R23'
    im_openclose, pspath+psname, postscript=postscript

    lineratio, atlasnodust, ['OII_3727','OIII_5007','OIII_4959'], 'H_BETA', $
      '', '', x, xerr, dum1, dum2, index=indx, nindex=nindx, snrcut=snrcut
    
    if keyword_set(nfgs) then begin
       lineratio, nfgsnodust, ['OII_3727','OIII_5007','OIII_4959'], 'H_BETA', $
         '', '', xnfgs, xerrnfgs, dum1, dum2, index=indxnfgs, snrcut=snrcut
    endif

;   xtitle = 'log [([O II] + [O III])/H\beta]'
    xtitle = 'log R_{23}'

    binsize = 0.1
    
    plothist, x, bin=binsize, xbin, ybin, /noplot, /halfbin
    if keyword_set(nfgs) then plothist, xnfgs, bin=binsize, xbinnfgs, ybinnfgs, /noplot, /halfbin

    xrange = R23range
    yrange = minmax(ybin)*[1.0,1.1]

    djs_plot, [0], [0], /nodata, xthick=postthick, ythick=postthick, charsize=2.0, $
      charthick=postthick, thick=postthick, xtitle=textoidl(xtitle), ytitle='Number', $
      xstyle=1, ystyle=1, xrange=xrange, yrange=yrange
    plothist, x, bin=binsize, /overplot, color=djs_icolor('blue'), thick=postthick, /halfbin
    if keyword_set(nfgs) then plothist, xnfgs, bin=binsize, thick=postthick, line=2, /halfbin, $
       /overplot, /fill, /fline, forientation=45, fcolor=djs_icolor('red'), $
      color=djs_icolor('red')

;   legend, '(d)', /left, /top, charthick=postthick, charsize=2.0, box=0, clear=postscript

    im_openclose, postscript=postscript, /close    

; ------------------------------------------------------------
; E(Hb-Ha) vs [O II]/Ha_obs
; ------------------------------------------------------------

    psname = 'ehbha_vs_oiiha_obs'
    im_openclose, pspath+psname, postscript=postscript

    cut = where(atlasnodust.ehbha_err gt 0.0)
    lineratio, atlasdust[cut], '', '', 'OII_3727', 'H_ALPHA', $
      dum1, dum2, y, yerr, index=indx, nindex=nindx, snrcut=snrcut

    x = atlasnodust[cut[indx]].ehbha
    xerr = atlasnodust[cut[indx]].ehbha_err
    xebv = x/(k_lambda(4861,/odonnell)-k_lambda(6563,/odonnell))

    ybig = y
    xbig = x
    
    if keyword_set(kenn92) then begin
       cutkenn92 = where(kenn92nodust.ehbha_err gt 0.0)
       lineratio, kenn92dust[cutkenn92], '', '', 'OII_3727', 'H_ALPHA', $
         dum1, dum2, ykenn92, yerrkenn92, index=indxkenn92, snrcut=snrcut
       xkenn92 = kenn92nodust[cutkenn92[indxkenn92]].ehbha
       xerrkenn92 = kenn92nodust[cutkenn92[indxkenn92]].ehbha_err
    endif

    if keyword_set(nfgs) then begin
       cutnfgs = where(nfgsnodust.ehbha_err gt 0.0)
       lineratio, nfgsdust[cutnfgs], '', '', 'OII_3727', 'H_ALPHA', $
         dum1, dum2, ynfgs, yerrnfgs, index=indxnfgs, snrcut=snrcut
       xnfgs = nfgsnodust[cutnfgs[indxnfgs]].ehbha
       xerrnfgs = nfgsnodust[cutnfgs[indxnfgs]].ehbha_err

       ybig = [y,ynfgs]
       xbig = [x,xnfgs]
    endif

    xtitle = 'E(H\beta-H\alpha) [mag]'
    ytitle = 'log ([O II]/H\alpha)_{obs}'

    xrange = ehbharange
    yrange = oiihacorrange
;   yrange = oiiharange
;   yrange = [-1.0,0.4]

    atlas1d_lineplot, x, y, xerr, yerr, plottype=3, postscript=postscript, $
      xtitle=xtitle, ytitle=ytitle, xrange=xrange, yrange=yrange, legendtype=0, $
      xstyle=11, ymargin=[4,3], /right, /top, $
;     xkenn92=xkenn92, ykenn92=ykenn92, xerrkenn92=xerrkenn92, yerrkenn92=yerrkenn92, $
      xnfgs=xnfgs, ynfgs=ynfgs, xerrnfgs=xerrnfgs, yerrnfgs=yerrnfgs
    axis, /xaxis, xrange=interpol(xebv,x,!x.crange), xthick=postthick, $
      charthick=postthick, charsize=2.0, xtitle = 'E(B-V) [mag]', xsty=1

    rcor = r_correlate(xbig,ybig,zd=zd,probd=probd)
    splog, 'Spearman rank for E(Hb-Ha) vs [O II]/Ha observed: ', rcor, probd

; label a low-metallicity galaxy    
    
;   xyouts, 0.05, -0.91, 'I Zw 18', charsize=charsize, charthick=postthick, /data

; fit a line

    if keyword_set(nfgs) then begin
       xpts = [x,xnfgs]
       ypts = [y,ynfgs]
    endif else begin
       xpts = x
       ypts = y
    endelse
    
    sixlin, xpts, ypts, a, siga, b, sigb
    ii = 2                      ; Ordinary Least Squares Bisector
    coeff = [a[ii],b[ii]]

    axis = findgen((1.2-(-0.1))/0.01)*0.01-0.1
    yfit = poly(axis,coeff)    
;   djs_oplot, axis, yfit, thick=postthick, line=0

; overplot reddening relations: CCM (with new coefficients from
; O'Donnell), Calzetti (2001), Charlot & Fall (2001), and SMC.  assume
; an intrinsic ratio log [O II]/Ha 

    ratio = 0.0
;   ratio = interpol(yfit,axis,0.0)
    splog, 'Ratio = '+string(10^ratio,format='(F0.0)')
;   ratio = 0.046

    x_ebv = findgen((1.1-(-0.5))/0.1)*0.1+(-0.5)
    yodonnell =  -0.4*x_ebv*(k_lambda(3727.0,/odonnell)-k_lambda(6563.0,/odonnell)) + ratio
;   yccm =  -0.4*x_ebv*(k_lambda(3727.0,/ccm)-k_lambda(6563.0,/ccm)) + ratio
;   ycalzetti = -0.4*x_ebv*(k_lambda(3727,/calzetti)-k_lambda(6563,/calzetti)) + ratio
    ycharlot = -0.4*x_ebv*(k_lambda(3727,/charlot)-k_lambda(6563,/charlot)) + ratio
    ysmc = -0.4*x_ebv*(k_lambda(3727,/smc)-k_lambda(6563,/smc)) + ratio

    djs_oplot, x_ebv, yodonnell, line=0, thick=postthick
;   djs_oplot, x_ebv, yccm, line=2, thick=postthick
;   djs_oplot, x_ebv, ycalzetti, line=0, thick=3.0
    djs_oplot, x_ebv, ycharlot, line=1, thick=3.0
    djs_oplot, x_ebv, ysmc, line=2, thick=postthick

    label = ['Milky Way','SMC Bar','Charlot & Fall']
;   label = ['Bisector Fit','Milky Way','SMC Bar']n
;   linestyle = [0,2,1]
    linestyle = [0,2,1]
    legend, label, /right, /top, box=0, linestyle=linestyle, $
      charsize=charsize_3, charthick=postthick, thick=postthick, clear=postscript

; overplot a histogram

    binsize = 0.1
    plothist, ybig, bin=binsize, xbin, ybin, /noplot, /halfbin

    yhistrange = minmax(ybin)*[1.0,1.05]

    djs_plot, [0], [0], /nodata, /noerase, xthick=postthick, ythick=postthick, charsize=1.0, $
      charthick=postthick, thick=postthick, xtitle=textoidl('log ([O II]/H\alpha)'), ytitle='Number', $
      xstyle=1, ystyle=1, xrange=yrange, yrange=yhistrange, /normal, $
      position=[0.35,0.25,0.6,0.45] ; [0.6,0.22,0.85,0.42]
    plothist, ybig, bin=binsize, /overplot, thick=postthick, /fill, /fline, /halfbin, $
      forientation=45, fcolor=djs_icolor('grey'), color=djs_icolor('grey')

    im_openclose, postscript=postscript, /close    

; ------------------------------------------------------------
; E(Hb-Ha) vs [O II]_cor/Ha_cor
; ------------------------------------------------------------

    psname = 'ehbha_vs_oiiha_cor'
    im_openclose, pspath+psname, postscript=postscript

    lineratio, atlasnodust, '', '', 'OII_3727', 'H_ALPHA', $
      dum1, dum2, y, yerr, index=indx, nindex=nindx, snrcut=snrcut

    x = atlasnodust[indx].ehbha
    xerr = atlasnodust[indx].ehbha_err
    xebv = x/(k_lambda(4861,/odonnell)-k_lambda(6563,/odonnell))

    xbig = x
    ybig = y
    
    if keyword_set(kenn92) then begin
       lineratio, kenn92nodust, '', '', 'OII_3727', 'H_ALPHA', $
         dum1, dum2, ykenn92, yerrkenn92, index=indxkenn92, snrcut=snrcut
       xkenn92 = kenn92nodust[indxkenn92].ehbha
       xerrkenn92 = kenn92nodust[indxkenn92].ehbha_err
    endif

    if keyword_set(nfgs) then begin
       lineratio, nfgsnodust, '', '', 'OII_3727', 'H_ALPHA', $
         dum1, dum2, ynfgs, yerrnfgs, index=indxnfgs, snrcut=snrcut
       xnfgs = nfgsnodust[indxnfgs].ehbha
       xerrnfgs = nfgsnodust[indxnfgs].ehbha_err

       xbig = [xbig,xnfgs]
       ybig = [ybig,ynfgs]
    endif

    ystats = im_stats(ybig,sigrej=3.0,mask=mask,/verbose)
    out = where(mask eq 0L,nout)
    splog, '3-sigma outliers: '+string(nout,format='(I0)')+'/'+string(n_elements(xbig),format='(I0)')
    
    xtitle = 'E(H\beta-H\alpha) [mag]'
    ytitle = 'log ([O II]/H\alpha)_{cor}'

    xrange = ehbharange
    yrange = oiihacorrange
;   yrange = oiiharange
;   yrange = [-1.0,0.4]

    rcor = r_correlate(xbig,ybig,zd=zd,probd=probd)
;   rcor = r_correlate(xbig[where(mask)],ybig[where(mask)],zd=zd,probd=probd)
    splog, 'Spearman rank for E(Hb-Ha) vs [O II]/Ha corrected: ', rcor, probd, zd

; suppress external data on this plot

    atlas1d_lineplot, x, y, xerr, yerr, plottype=3, postscript=postscript, $
      xtitle=xtitle, ytitle=ytitle, xrange=xrange, yrange=yrange, legendtype=0, $
      xstyle=11, ymargin=[4,3], /left, /top, $
;     xkenn92=xkenn92, ykenn92=ykenn92, xerrkenn92=xerrkenn92, yerrkenn92=yerrkenn92, $
      xnfgs=xnfgs, ynfgs=ynfgs, xerrnfgs=xerrnfgs, yerrnfgs=yerrnfgs, /save
    axis, /xaxis, xrange=interpol(xebv,x,!x.crange), xthick=postthick, charthick=postthick, $
      charsize=2.0, xtitle = 'E(B-V) [mag]', xsty=1

;   djs_oplot, !x.crange, ystats.median*[1,1], line=0, thick=postthick
;   djs_oplot, !x.crange, ystats.median+ystats.sigma*[1,1], line=2, thick=postthick
;   djs_oplot, !x.crange, ystats.median-ystats.sigma*[1,1], line=2, thick=postthick

; overplot a histogram

    binsize = 0.1
    plothist, ybig, bin=binsize, xbin, ybin, /noplot, /halfbin

    yhistrange = minmax(ybin)*[1.0,1.05]

    djs_plot, [0], [0], /nodata, /noerase, xthick=postthick, ythick=postthick, charsize=1.0, $
      charthick=postthick, thick=postthick, xtitle=textoidl('log ([O II]/H\alpha)'), ytitle='Number', $
      xstyle=1, ystyle=1, xrange=yrange, yrange=yhistrange, /normal, $
      position=[0.35,0.25,0.6,0.45] ; [0.6,0.22,0.85,0.42]
    plothist, ybig, bin=binsize, /overplot, thick=postthick, /fill, /fline, /halfbin, $
      forientation=45, fcolor=djs_icolor('grey'), color=djs_icolor('grey')

    im_openclose, postscript=postscript, /close    

; ------------------------------------------------------------
; E(Hb-Ha) vs [O II]_obs/Ha_cor
; ------------------------------------------------------------

    psname = 'ehbha_vs_oiiobs_hacor'
    im_openclose, pspath+psname, postscript=postscript

    indx = where((atlasnodust.ehbha_err gt 0.0) and $
      (atlasdust.oii_3727[0]/atlasdust.oii_3727[1] gt snrcut) and $
      (atlasnodust.h_alpha[0]/atlasnodust.h_alpha[1] gt snrcut),nindx)

    x = atlasnodust[indx].ehbha
    xerr = atlasnodust[indx].ehbha_err
    xebv = x/(k_lambda(4861,/odonnell)-k_lambda(6563,/odonnell))

    y = alog10(atlasdust[indx].oii_3727[0]/atlasnodust[indx].h_alpha[0])
    yerr = im_compute_error(atlasdust[indx].oii_3727[0],atlasdust[indx].oii_3727[1],$
      atlasnodust[indx].h_alpha[0],atlasnodust[indx].h_alpha[1],/log)

    xbig = x
    ybig = y
    
    if keyword_set(nfgs) then begin

       indxnfgs = where((nfgsnodust.ehbha_err gt 0.0) and $
         (nfgsdust.oii_3727[0]/nfgsdust.oii_3727[1] gt snrcut) and $
         (nfgsnodust.h_alpha[0]/nfgsnodust.h_alpha[1] gt snrcut),nindxnfgs)

       xnfgs = nfgsnodust[indxnfgs].ehbha
       xnfgserr = nfgsnodust[indxnfgs].ehbha_err

       ynfgs = alog10(nfgsdust[indxnfgs].oii_3727[0]/nfgsnodust[indxnfgs].h_alpha[0])
       yerrnfgs = im_compute_error(nfgsdust[indxnfgs].oii_3727[0],nfgsdust[indxnfgs].oii_3727[1],$
         nfgsnodust[indxnfgs].h_alpha[0],nfgsnodust[indxnfgs].h_alpha[1],/log)

       xbig = [xbig,xnfgs]
       ybig = [ybig,ynfgs]

    endif

    xtitle = 'E(H\beta-H\alpha) [mag]'
    ytitle = 'log ([O II]_{obs}/H\alpha_{cor})'

    xrange = ehbharange
    yrange = oiihacorrange

    atlas1d_lineplot, x, y, xerr, yerr, plottype=3, postscript=postscript, $
      xtitle=xtitle, ytitle=ytitle, xrange=xrange, yrange=yrange, legendtype=1, $
      xstyle=11, ymargin=[4,3], /right, /top, xnfgs=xnfgs, ynfgs=ynfgs, xerrnfgs=xerrnfgs, $
      yerrnfgs=yerrnfgs
    axis, /xaxis, xrange=interpol(xebv,x,!x.crange), xthick=postthick, charthick=postthick, $
      charsize=2.0, xtitle = 'E(B-V) [mag]', xsty=1

; overplot a histogram

    binsize = 0.1
    plothist, ybig, bin=binsize, xbin, ybin, /noplot, /halfbin

    yhistrange = minmax(ybin)*[1.0,1.05]

    djs_plot, [0], [0], /nodata, /noerase, xthick=postthick, ythick=postthick, charsize=1.0, $
      charthick=postthick, thick=postthick, xtitle=textoidl('log ([O II]/H\alpha)'), ytitle='Number', $
      xstyle=1, ystyle=1, xrange=yrange, yrange=yhistrange, /normal, $
      position=[0.35,0.25,0.6,0.45] ;[0.65,0.45,0.9,0.65]
    plothist, ybig, bin=binsize, /overplot, thick=postthick, /fill, /fline, /halfbin, $
      forientation=45, fcolor=djs_icolor('grey'), color=djs_icolor('grey')

    im_openclose, postscript=postscript, /close    

; ------------------------------------------------------------
; E(Hb-Ha) vs [O II]_obs/Ha_meancor
; ------------------------------------------------------------

; apply the standard 1 mag extinction correction at H-alpha    
    
    psname = 'ehbha_vs_oiiobs_hameancor'
    im_openclose, pspath+psname, postscript=postscript

    indx = where((atlasnodust.ehbha_err gt 0.0) and $
      (atlasdust.oii_3727[0]/atlasdust.oii_3727[1] gt snrcut) and $
      (atlasdust.h_alpha[0]/atlasdust.h_alpha[1] gt snrcut),nindx)

    x = atlasnodust[indx].ehbha
    xerr = atlasnodust[indx].ehbha_err
    xebv = x/(k_lambda(4861,/odonnell)-k_lambda(6563,/odonnell))

    y1 = atlasdust[indx].oii_3727[0]
    y1err = atlasdust[indx].oii_3727[1]

    y2 = atlasdust[indx].h_alpha[0]*10^(0.4)    ; 1 mag correction
    y2err = atlasdust[indx].h_alpha[1]*10^(0.4)
    
    y = alog10(y1/y2)
    yerr = im_compute_error(y1,y1err,y2,y2err,/log)

    xbig = x
    ybig = y
    
    if keyword_set(nfgs) then begin

       indxnfgs = where((nfgsnodust.ehbha_err gt 0.0) and $
         (nfgsdust.oii_3727[0]/nfgsdust.oii_3727[1] gt snrcut) and $
         (nfgsdust.h_alpha[0]/nfgsdust.h_alpha[1] gt snrcut),nindxnfgs)

       xnfgs = nfgsnodust[indxnfgs].ehbha
       xnfgserr = nfgsnodust[indxnfgs].ehbha_err

       y1nfgs = nfgsdust[indxnfgs].oii_3727[0]
       y1errnfgs = nfgsdust[indxnfgs].oii_3727[1]

       y2nfgs = nfgsdust[indxnfgs].h_alpha[0]*10^(0.4) ; 1 mag correction
       y2errnfgs = nfgsdust[indxnfgs].h_alpha[1]*10^(0.4)
       
       ynfgs = alog10(y1nfgs/y2nfgs)
       yerrnfgs = im_compute_error(y1nfgs,y1errnfgs,y2nfgs,y2errnfgs,/log)

       xbig = [xbig,xnfgs]
       ybig = [ybig,ynfgs]

    endif

    xtitle = 'E(H\beta-H\alpha) [mag]'
    ytitle = 'log ([O II]_{obs}/H\alpha_{cor})'

    xrange = ehbharange
    yrange = oiihacorrange

    atlas1d_lineplot, x, y, xerr, yerr, plottype=3, postscript=postscript, $
      xtitle=xtitle, ytitle=ytitle, xrange=xrange, yrange=yrange, legendtype=1, $
      xstyle=11, ymargin=[4,3], /right, /top, xnfgs=xnfgs, ynfgs=ynfgs, xerrnfgs=xerrnfgs, $
      yerrnfgs=yerrnfgs
    axis, /xaxis, xrange=interpol(xebv,x,!x.crange), xthick=postthick, charthick=postthick, $
      charsize=2.0, xtitle = 'E(B-V) [mag]', xsty=1

; overplot a histogram

    binsize = 0.1
    plothist, ybig, bin=binsize, xbin, ybin, /noplot, /halfbin

    yhistrange = minmax(ybin)*[1.0,1.05]

    djs_plot, [0], [0], /nodata, /noerase, xthick=postthick, ythick=postthick, charsize=1.0, $
      charthick=postthick, thick=postthick, xtitle=textoidl('log ([O II]/H\alpha)'), ytitle='Number', $
      xstyle=1, ystyle=1, xrange=yrange, yrange=yhistrange, /normal, $
      position=[0.35,0.25,0.6,0.45] ;[0.65,0.45,0.9,0.65]
    plothist, ybig, bin=binsize, /overplot, thick=postthick, /fill, /fline, /halfbin, $
      forientation=45, fcolor=djs_icolor('grey'), color=djs_icolor('grey')

    im_openclose, postscript=postscript, /close    

; ------------------------------------------------------------
; [N II]/[O III] vs [O II]_cor/Ha_cor
; ------------------------------------------------------------

    psname = 'niioiii_vs_oiiha_cor'
    im_openclose, pspath+psname, postscript=postscript

    lineratio, atlasnodust, 'NII_6584', 'OIII_5007', 'OII_3727', 'H_ALPHA', $
      x, xerr, y, yerr, index=indx, nindex=nindx, snrcut=snrcut

    if keyword_set(kenn92) then begin
       lineratio, kenn92nodust, 'NII_6584', 'OIII_5007', 'OII_3727', 'H_ALPHA', $
         xkenn92, xerrkenn92, ykenn92, yerrkenn92, index=indxkenn92, snrcut=snrcut
    endif

    if keyword_set(nfgs) then begin
       lineratio, nfgsnodust, 'NII_6584', 'OIII_5007', 'OII_3727', 'H_ALPHA', $
         xnfgs, xerrnfgs, ynfgs, yerrnfgs, index=indxnfgs, snrcut=snrcut
    endif

    xtitle = 'log ([N II] \lambda6584/[O III] \lambda5007)_{cor}'
    ytitle = 'log ([O II]/H\alpha)_{cor}'

    xrange = niioiiirange
    yrange = oiiharange

    if keyword_set(hiiregions) then begin
       good = where((hii.nii_6584_oiii_5007 gt -900.0) and (hii.oii_h_alpha gt -900.0))
       xregion = hii[good].nii_6584_oiii_5007 & xerrregion = hii[good].nii_6584_oiii_5007_err
       yregion = hii[good].oii_h_alpha & yerrregion = hii[good].oii_h_alpha_err
    endif
       
    atlas1d_lineplot, x, y, xerr, yerr, plottype=3, postscript=postscript, $
      xtitle=xtitle, ytitle=ytitle, xrange=xrange, yrange=yrange, legendtype=0, $
      /left, /top, $
      xregion=xregion, yregion=yregion, xerrregion=xerrregion, yerrregion=yerrregion, $
;     xkenn92=xkenn92, ykenn92=ykenn92, xerrkenn92=xerrkenn92, yerrkenn92=yerrkenn92, $
      xnfgs=xnfgs, ynfgs=ynfgs, xerrnfgs=xerrnfgs, yerrnfgs=yerrnfgs

;   if keyword_set(kewley_grids) then begin
;      plot_kewley_grids, plotnumber=12, model=3, labeltype=4, /noZgrid, $
;        /overplot, Umax=-2.5, Umin=-3.5, Zmax=2.9, Zmin=0.15, postscript=postscript
;      plot_kewley_grids, plotnumber=12, model=3, labeltype=2, /noUgrid, $
;        /overplot, Umax=-2.5, Umin=-3.5, Zmax=2.9, Zmin=0.15, postscript=postscript
;   endif

    im_openclose, postscript=postscript, /close    

; ------------------------------------------------------------
; [O III]/[O II] vs [O II]/Ha_cor
; ------------------------------------------------------------

    psname = 'oiiioii_vs_oiiha_cor'
    im_openclose, pspath+psname, postscript=postscript

    lineratio, atlasnodust, 'OIII_5007', 'OII_3727', 'OII_3727', 'H_ALPHA', $
      x, xerr, y, yerr, index=indx, nindex=nindx, snrcut=snrcut

    if keyword_set(kenn92) then begin
       lineratio, kenn92nodust, 'OIII_5007', 'OII_3727', 'OII_3727', 'H_ALPHA', $
         xkenn92, xerrkenn92, ykenn92, yerrkenn92, index=indxkenn92, snrcut=snrcut
    endif    

    if keyword_set(nfgs) then begin
       lineratio, nfgsnodust, 'OIII_5007', 'OII_3727', 'OII_3727', 'H_ALPHA', $
         xnfgs, xerrnfgs, ynfgs, yerrnfgs, index=indxnfgs, snrcut=snrcut
    endif    
    
    xtitle = 'log ([O III] \lambda5007/[O II])_{cor}'
    ytitle = 'log ([O II]/H\alpha)_{cor}'

    xrange = oiiioiirange
    yrange = oiiharange

; HII regions

    if keyword_set(hiiregions) then begin
       good = where((hii.oii_h_alpha gt -900.0) and (hii.oiii_5007_oii gt -900.0))
       xregion = hii[good].oiii_5007_oii & xerrregion = hii[good].oiii_5007_oii_err
       yregion = hii[good].oii_h_alpha & yerrregion = hii[good].oii_h_alpha_err
    endif

    atlas1d_lineplot, x, y, xerr, yerr, plottype=3, postscript=postscript, $
      xtitle=xtitle, ytitle=ytitle, xrange=xrange, yrange=yrange, legendtype=1, $
      /right, /top, $
      xregion=xregion, yregion=yregion, xerrregion=xerrregion, yerrregion=yerrregion, $
;     xkenn92=xkenn92, ykenn92=ykenn92, xerrkenn92=xerrkenn92, yerrkenn92=yerrkenn92, $
      xnfgs=xnfgs, ynfgs=ynfgs, xerrnfgs=xerrnfgs, yerrnfgs=yerrnfgs
    
    if keyword_set(kewley_grids) then plot_kewley_grids, plotnumber=4, model=3, $
      labeltype=1, /overplot, postscript=postscript

    im_openclose, postscript=postscript, /close    

; ------------------------------------------------------------
; [O III]/Hb vs [O II]/Ha_cor
; ------------------------------------------------------------

    psname = 'oiiihb_vs_oiiha_cor'
    im_openclose, pspath+psname, postscript=postscript

    lineratio, atlasnodust, 'OIII_5007', 'H_BETA', 'OII_3727', 'H_ALPHA', $
      x, xerr, y, yerr, index=indx, nindex=nindx, snrcut=snrcut

    if keyword_set(kenn92) then begin
       lineratio, kenn92nodust, 'OIII_5007', 'H_BETA', 'OII_3727', 'H_ALPHA', $
         xkenn92, xerrkenn92, ykenn92, yerrkenn92, index=indxkenn92, snrcut=snrcut
    endif    

    if keyword_set(nfgs) then begin
       lineratio, nfgsnodust, 'OIII_5007', 'H_BETA', 'OII_3727', 'H_ALPHA', $
         xnfgs, xerrnfgs, ynfgs, yerrnfgs, index=indxnfgs, snrcut=snrcut
    endif    
    
    xtitle = 'log ([O III] \lambda5007/H\beta)_{cor}'
    ytitle = 'log ([O II]/H\alpha)_{cor}'

    xrange = oiiihbrange
    yrange = oiiharange

; HII regions

    if keyword_set(hiiregions) then begin
       good = where((hii.oii_h_alpha gt -900.0) and (hii.oiii_5007_h_beta gt -900.0))
       xregion = hii[good].oiii_5007_h_beta & xerrregion = hii[good].oiii_5007_h_beta_err
       yregion = hii[good].oii_h_alpha & yerrregion = hii[good].oii_h_alpha_err
    endif

    atlas1d_lineplot, x, y, xerr, yerr, plottype=3, postscript=postscript, $
      xtitle=xtitle, ytitle=ytitle, xrange=xrange, yrange=yrange, legendtype=1, $
      /right, /top, $
      xregion=xregion, yregion=yregion, xerrregion=xerrregion, yerrregion=yerrregion, $
;     xkenn92=xkenn92, ykenn92=ykenn92, xerrkenn92=xerrkenn92, yerrkenn92=yerrkenn92, $
      xnfgs=xnfgs, ynfgs=ynfgs, xerrnfgs=xerrnfgs, yerrnfgs=yerrnfgs
    
    if keyword_set(kewley_grids) then plot_kewley_grids, plotnumber=13, model=3, $
      labeltype=1, /overplot, postscript=postscript

    im_openclose, postscript=postscript, /close    

; ------------------------------------------------------------
; [N II]/[O II] vs [O II]/Ha_cor
; ------------------------------------------------------------

    psname = 'niioii_vs_oiiha_cor'
    im_openclose, pspath+psname, postscript=postscript

    lineratio, atlasnodust, 'NII_6584', 'OII_3727', 'OII_3727', 'H_ALPHA', $
      x, xerr, y, yerr, index=indx, nindex=nindx, snrcut=snrcut

    if keyword_set(kenn92) then begin
       lineratio, kenn92nodust, 'NII_6584', 'OII_3727', 'OII_3727', 'H_ALPHA', $
         xkenn92, xerrkenn92, ykenn92, yerrkenn92, index=indxkenn92, snrcut=snrcut
    endif    

    if keyword_set(nfgs) then begin
       lineratio, nfgsnodust, 'NII_6584', 'OII_3727', 'OII_3727', 'H_ALPHA', $
         xnfgs, xerrnfgs, ynfgs, yerrnfgs, index=indxnfgs, snrcut=snrcut
    endif    

    xtitle = 'log ([N II] \lambda6584/[O II])_{cor}'
    ytitle = 'log ([O II]/H\alpha)_{cor}'

    xrange = niioiirange
    yrange = oiiharange

; HII regions

    if keyword_set(hiiregions) then begin
       good = where((hii.nii_6584_oii gt -900.0) and (hii.oii_h_alpha gt -900.0))
       xregion = hii[good].nii_6584_oii & xerrregion = hii[good].nii_6584_oii_err
       yregion = hii[good].oii_h_alpha & yerrregion = hii[good].oii_h_alpha_err
    endif

    atlas1d_lineplot, x, y, xerr, yerr, plottype=3, postscript=postscript, $
      xtitle=xtitle, ytitle=ytitle, xrange=xrange, yrange=yrange, legendtype=1, $
      /right, /top, $
      xregion=xregion, yregion=yregion, xerrregion=xerrregion, yerrregion=yerrregion, $
;     xkenn92=xkenn92, ykenn92=ykenn92, xerrkenn92=xerrkenn92, yerrkenn92=yerrkenn92, $
      xnfgs=xnfgs, ynfgs=ynfgs, xerrnfgs=xerrnfgs, yerrnfgs=yerrnfgs

    if keyword_set(kewley_grids) then plot_kewley_grids, plotnumber=14, model=3, $
      labeltype=1, Zmax=2.0, /overplot, postscript=postscript

; show that assuming [N II]/Ha = constant gives you a 1/x relation! 

    nii_ha = -0.408
    nii_ha_var = 0.1
    
    xaxis = findgen(((+0.2)-(-0.9))/0.1)*0.1+(-0.9)
    yline = poly(xaxis,[nii_ha,-1.0])
    yleft = poly(xaxis,[nii_ha-nii_ha_var,-1.0])
    yrght = poly(xaxis,[nii_ha+nii_ha_var,-1.0])

;   oplot, xaxis, yline, linestyle=0, thick=postthick
;   oplot, xaxis, yleft, linestyle=2, thick=postthick
;   oplot, xaxis, yrght, linestyle=2, thick=postthick

;   legend, textoidl(['log [N II]/H\alpha = '+string(nii_ha,format='(F5.2)')+$
;     '\pm'+string(nii_ha_var,format='(F4.2)')]), $
;     /left, /top, box=0, charsize=1.8, charthick=postthick
    
    im_openclose, postscript=postscript, /close    

; ------------------------------------------------------------
; [N II]/[S II] vs [O II]_cor/Ha_cor
; ------------------------------------------------------------

    psname = 'niisii_vs_oiiha_cor'
    im_openclose, pspath+psname, postscript=postscript

    lineratio, atlasnodust, 'NII_6584', ['SII_6716','SII_6731'], 'OII_3727', 'H_ALPHA', $
      x, xerr, y, yerr, index=indx, nindex=nindx, snrcut=snrcut
    
    if keyword_set(kenn92) then begin
       lineratio, kenn92nodust, 'NII_6584', ['SII_6716','SII_6731'], 'OII_3727', 'H_ALPHA', $
         xkenn92, xerrkenn92, ykenn92, yerrkenn92, index=indxkenn92, snrcut=snrcut
    endif
       
    if keyword_set(nfgs) then begin
       lineratio, nfgsnodust, 'NII_6584', ['SII_6716','SII_6731'], 'OII_3727', 'H_ALPHA', $
         xnfgs, xerrnfgs, ynfgs, yerrnfgs, index=indxnfgs, snrcut=snrcut
    endif
       
    xtitle = 'log ([N II] \lambda6584/[S II] \lambda\lambda6716,6731)_{cor}'
    ytitle = 'log ([O II]/H\alpha)_{cor}'

    xrange = niisiirange
    yrange = oiiharange

; HII regions

    if keyword_set(hiiregions) then begin
       good = where((hii.nii_6584_sii gt -900.0) and (hii.oii_h_alpha gt -900.0))
       xregion = hii[good].nii_6584_sii & xerrregion = hii[good].nii_6584_sii_err
       yregion = hii[good].oii_h_alpha & yerrregion = hii[good].oii_h_alpha_err
    endif

    atlas1d_lineplot, x, y, xerr, yerr, plottype=3, postscript=postscript, $
      xtitle=xtitle, ytitle=ytitle, xrange=xrange, yrange=yrange, legendtype=1, $
      /right, /top, $
      xregion=xregion, yregion=yregion, xerrregion=xerrregion, yerrregion=yerrregion, $
;     xkenn92=xkenn92, ykenn92=ykenn92, xerrkenn92=xerrkenn92, yerrkenn92=yerrkenn92, $
      xnfgs=xnfgs, ynfgs=ynfgs, xerrnfgs=xerrnfgs, yerrnfgs=yerrnfgs
    
    if keyword_set(kewley_grids) then plot_kewley_grids, plotnumber=5, labeltype=1, $
      model=3, /overplot, postscript=postscript
    
    im_openclose, postscript=postscript, /close    

; ------------------------------------------------------------
; EW(Ha) vs EW([O II])
; ------------------------------------------------------------

    psname = 'ewha_vs_ewoii'
    im_openclose, pspath+psname, postscript=postscript

    lineratio, atlasdust, 'H_ALPHA_EW', '', 'OII_3727_EW', '', x, xerr, $
      y, yerr, index=indx, nindex=nindx, snrcut=snrcut

    if keyword_set(kenn92) then begin
       lineratio, kenndust, 'H_ALPHA_EW', '', 'OII_3727_EW', '', $
         xkenn92, xerrkenn92, ykenn92, yerrkenn92, index=indxkenn92, snrcut=snrcut
    endif

    if keyword_set(nfgs) then begin
       lineratio, nfgsnodust, 'H_ALPHA_EW', '', 'OII_3727_EW', '', xnfgs, xerrnfgs, $
         ynfgs, yerrnfgs, index=indxnfgs, nindex=nindxnfgs, snrcut=snrcut
    endif

    xtitle = 'log EW(H\alpha)  ['+angstrom()+']' 
    ytitle = 'log EW([O II])  ['+angstrom()+']'
    
    xrange = ewharange
    yrange = ewharange

    atlas1d_lineplot, x, y, xerr, yerr, plottype=3, postscript=postscript, $
      xtitle=xtitle, ytitle=ytitle, xrange=xrange, yrange=yrange, legendtype=1, $
;     xkenn92=xkenn92, ykenn92=ykenn92, xerrkenn92=xerrkenn92, yerrkenn92=yerrkenn92, $
      xnfgs=xnfgs, ynfgs=ynfgs, xerrnfgs=xerrnfgs, yerrnfgs=yerrnfgs

    im_openclose, postscript=postscript, /close        
    
; ------------------------------------------------------------
; 41-50 vs EW(OII)
; ------------------------------------------------------------

    psname = '41_50_vs_ewoii'
    im_openclose, pspath+psname, postscript=postscript

    cut = where(atlasdust.c41_50[1] gt 0.0)
    lineratio, atlasdust[cut], '', '', 'OII_3727_EW', '', $
      dum1, dum2, y, yerr, index=indx, nindex=nindx, snrcut=snrcut

    if keyword_set(kenn92) then begin
       cutkenn92 = where(kenn92dust.c41_50[1] gt 0.0)
       lineratio, kenn92dust[cut], '', '', 'OII_3727_EW', '', $
         dum1, dum2, ykenn92, yerrkenn92, index=indxkenn92
       xkenn92 = kenn92dust[cutkenn92[indxkenn92]].c41_50[0]
       xerrkenn92 = kenn92dust[cutkenn92[indxkenn92]].c41_50[1]
    endif
       
    if keyword_set(nfgs) then begin
       cutnfgs = where(nfgsdust.c41_50[1] gt 0.0)
       lineratio, nfgsdust[cutnfgs], '', '', 'OII_3727_EW', '', $
         dum1, dum2, ynfgs, yerrnfgs, index=indxnfgs
       xnfgs = nfgsdust[cutnfgs[indxnfgs]].c41_50[0]
       xerrnfgs = nfgsdust[cutnfgs[indxnfgs]].c41_50[1]
    endif

    x = atlasdust[cut[indx]].c41_50[0]
    xerr = atlasdust[cut[indx]].c41_50[1]

    xtitle = '41 - 50' 
    ytitle = 'log EW([O II])'

    xrange = c4150range
    yrange = ewoiirange

    atlas1d_lineplot, x, y, xerr, yerr, plottype=3, postscript=postscript, $
      xtitle=xtitle, ytitle=ytitle, xrange=xrange, yrange=yrange, legendtype=1, $
      /right, /top, $
;     xkenn92=xkenn92, ykenn92=ykenn92, xerrkenn92=xerrkenn92, yerrkenn92=yerrkenn92, $
      xnfgs=xnfgs, ynfgs=ynfgs, xerrnfgs=xerrnfgs, yerrnfgs=yerrnfgs

    im_openclose, postscript=postscript, /close    

; ------------------------------------------------------------    
; D(4000) vs EW([O II])
; ------------------------------------------------------------    

    psname = 'D4000_vs_ewoii'
    im_openclose, pspath+psname, postscript=postscript

    cut = where(atlasdust.D4000_narrow[1] gt 0.0)
    lineratio, atlasdust[cut], '', '', 'OII_3727_EW', '', $
      dum1, dum2, y, yerr, index=indx, nindex=nindx, snrcut=snrcut

    if keyword_set(kenn92) then begin
       cutkenn92 = where(kenn92dust.D4000_narrow[1] gt 0.0)
       lineratio, kenn92dust[cut], '', '', 'OII_3727_EW', '', $
         dum1, dum2, ykenn92, yerrkenn92, index=indxkenn92
       xkenn92 = kenn92dust[cutkenn92[indxkenn92]].D4000_narrow[0]
       xerrkenn92 = kenn92dust[cutkenn92[indxkenn92]].D4000_narrow[1]
    endif
       
    if keyword_set(nfgs) then begin
       cutnfgs = where(nfgsdust.D4000_narrow[1] gt 0.0)
       lineratio, nfgsdust[cutnfgs], '', '', 'OII_3727_EW', '', $
         dum1, dum2, ynfgs, yerrnfgs, index=indxnfgs
       xnfgs = nfgsdust[cutnfgs[indxnfgs]].D4000_narrow[0]
       xerrnfgs = nfgsdust[cutnfgs[indxnfgs]].D4000_narrow[1]
    endif

    x = atlasdust[cut[indx]].D4000_narrow[0]
    xerr = atlasdust[cut[indx]].D4000_narrow[1]

    xtitle = 'D_{n}(4000)'
    ytitle = 'log EW([O II])'

    xrange = D4000range
    yrange = ewoiirange

    atlas1d_lineplot, x, y, xerr, yerr, plottype=3, postscript=postscript, $
      xtitle=xtitle, ytitle=ytitle, xrange=xrange, yrange=yrange, legendtype=1, $
      /right, /top, $
;     xkenn92=xkenn92, ykenn92=ykenn92, xerrkenn92=xerrkenn92, yerrkenn92=yerrkenn92, $
      xnfgs=xnfgs, ynfgs=ynfgs, xerrnfgs=xerrnfgs, yerrnfgs=yerrnfgs

    im_openclose, postscript=postscript, /close    

; ------------------------------------------------------------
; L(IR)/L(B) vs [O II]/Ha_obs
; ------------------------------------------------------------

    psname = 'LIR_LB_vs_oiiha_obs'
    im_openclose, pspath+psname, postscript=postscript

    cut = where(atlasdust.L_IR_L_B gt -900.0)
    lineratio, atlasdust[cut], '', '', 'OII_3727', 'H_ALPHA', $
      dum1, dum2, y, yerr, index=indx, nindex=nindx, snrcut=snrcut

    x = alog10(atlasdust[cut[indx]].L_IR_L_B)
    xerr = atlasdust[cut[indx]].L_IR_L_B_err/atlasdust[cut[indx]].L_IR_L_B/alog(10.0)
    
    if keyword_set(nfgs) then begin
       cutnfgs = where(nfgsdust.L_IR_L_B gt -900.0)
       lineratio, nfgsdust[cutnfgs], '', '', 'OII_3727', 'H_ALPHA', $
         dum1, dum2, ynfgs, yerrnfgs, index=indxnfgs, snrcut=snrcut
       xnfgs = alog10(nfgsdust[cutnfgs[indxnfgs]].L_IR_L_B)
       xerrnfgs = nfgsdust[cutnfgs[indxnfgs]].L_IR_L_B_err/nfgsdust[cutnfgs[indxnfgs]].L_IR_L_B/alog(10.0)
    endif

    xtitle = 'log (L_{IR}/L_{B})'
    ytitle = 'log ([O II]/H\alpha)_{obs}'

    xrange = iroptrange
    yrange = oiiharange

    atlas1d_lineplot, x, y, xerr, yerr, plottype=3, postscript=postscript, $
      xtitle=xtitle, ytitle=ytitle, xrange=xrange, yrange=yrange, legendtype=1, $
      /right, /top, xnfgs=xnfgs, ynfgs=ynfgs, xerrnfgs=xerrnfgs, yerrnfgs=yerrnfgs

    im_openclose, postscript=postscript, /close    

; ------------------------------------------------------------
; L(IR)/L(B) vs [O II]/Ha_cor
; ------------------------------------------------------------

    psname = 'L_IR_L_B_vs_oiiha_cor'
    im_openclose, pspath+psname, postscript=postscript

    cut = where(atlasdust.L_IR_L_B gt -900.0)
    lineratio, atlasnodust[cut], '', '', 'OII_3727', 'H_ALPHA', $
      dum1, dum2, y, yerr, index=indx, nindex=nindx, snrcut=snrcut

    x = alog10(atlasdust[cut[indx]].L_IR_L_B)
    xerr = atlasdust[cut[indx]].L_IR_L_B_err/atlasdust[cut[indx]].L_IR_L_B/alog(10.0)

    if keyword_set(nfgs) then begin
       cutnfgs = where(nfgsdust.L_IR_L_B gt -900.0)
       lineratio, nfgsnodust[cutnfgs], '', '', 'OII_3727', 'H_ALPHA', $
         dum1, dum2, ynfgs, yerrnfgs, index=indxnfgs, snrcut=snrcut
       xnfgs = alog10(nfgsdust[cutnfgs[indxnfgs]].L_IR_L_B)
       xerrnfgs = nfgsdust[cutnfgs[indxnfgs]].L_IR_L_B_err/nfgsdust[cutnfgs[indxnfgs]].L_IR_L_B/alog(10.0)
    endif
    
    xtitle = 'log (L_{IR}/L_{B})'
    ytitle = 'log ([O II]/H\alpha)_{cor}'

    xrange = iroptrange
    yrange = oiiharange

    atlas1d_lineplot, x, y, xerr, yerr, plottype=3, postscript=postscript, $
      xtitle=xtitle, ytitle=ytitle, xrange=xrange, yrange=yrange, legendtype=1, $
      /right, /top, xnfgs=xnfgs, ynfgs=ynfgs, xerrnfgs=xerrnfgs, yerrnfgs=yerrnfgs

    im_openclose, postscript=postscript, /close    

; ------------------------------------------------------------
; L(B) vs [O II]_obs/Ha_obs
; ------------------------------------------------------------
    
    psname = 'LB_vs_oiiha_obs'
    im_openclose, pspath+psname, postscript=postscript

    cut = where(atlasdust.b_lum_obs gt -900.0)
    lineratio, atlasdust[cut], '', '', 'OII_3727', 'H_ALPHA', $
      dum1, dum2, y, yerr, index=indx, nindex=nindx, snrcut=snrcut

    x = atlasdust[cut[indx]].b_lum_obs
    xerr = atlasdust[cut[indx]].b_lum_obs_err
    xabs = atlasdust[cut[indx]].m_b_obs

    if keyword_set(nfgs) then begin
       cutnfgs = where(nfgsdust.b_lum_obs gt -900.0)
       lineratio, nfgsdust[cutnfgs], '', '', 'OII_3727', 'H_ALPHA', $
         dum1, dum2, ynfgs, yerrnfgs, index=indxnfgs, snrcut=snrcut
       xnfgs = nfgsdust[cutnfgs[indxnfgs]].b_lum_obs
       xerrnfgs = nfgsdust[cutnfgs[indxnfgs]].b_lum_obs_err
    endif

    xtitle = 'log L(B) [L(B)'+sunsymbol()+']'
    ytitle = 'log ([O II]/H\alpha)_{obs}'

    xrange = LBrange
    yrange = oiiharange

    atlas1d_lineplot, x, y, xerr, yerr, plottype=3, postscript=postscript, $
      xtitle=xtitle, ytitle=ytitle, xrange=xrange, yrange=yrange, legendtype=1, $
      xstyle=11, ymargin=[4,3], /right, /top, xnfgs=xnfgs, ynfgs=ynfgs, xerrnfgs=xerrnfgs, $
      yerrnfgs=yerrnfgs
    axis, /xaxis, xrange=interpol(xabs,x,!x.crange), xthick=postthick, charthick=postthick, $
      charsize=2.0, xtitle = textoidl('M_{B} [mag]'), xsty=1

; overplot the NFGS relation and the Aragon-Salamanca relation

    nfgs_yint = 0.24   ; y intercept
    nfgs_slope = 0.087 ; slope
    
;   oplot, xabs, poly(xabs,[nfgs_yint,nfgs_slope]), line=0, thick=postthick
    
; overplot histograms
    
    w1 = where(xabs ge -19.2)
    w2 = where(xabs lt -19.2)

; bin 1    
    
    im_plothist, y[w1], xbin, ybin, bin=0.2, /noplot, /fraction, /halfbin
;   im_plothist, y[w1], bin=0.2, xsty=3, ysty=1, charsize=1.0, charthick=postthick, /halfbin, $
;     /noerase, xthick=postthick, ythick=postthick, position=[0.32,0.23,0.6,0.48], $
;     /normal, ytitle='Fraction', xtitle=textoidl(ytitle), thick=postthick, $
;     /fraction, yrange=minmax(ybin)*[0,1.1], $
;     /fill, forientation=45, fcolor=djs_icolor('grey')
    
; bin 2

    im_plothist, y[w2], xbin, ybin, bin=0.2, /noplot, /fraction, /halfbin
;   im_plothist, y[w2], bin=0.2, /overplot, thick=postthick, /fraction, line=2, /halfbin

    im_openclose, postscript=postscript, /close    
    
; ------------------------------------------------------------
; L(B) vs [O II]_cor/Ha_cor
; ------------------------------------------------------------

    psname = 'LB_vs_oiiha_cor'
    im_openclose, pspath+psname, postscript=postscript

    cut = where(atlasdust.b_lum_obs gt -900.0)
    lineratio, atlasnodust[cut], '', '', 'OII_3727', 'H_ALPHA', $
      dum1, dum2, y, yerr, index=indx, nindex=nindx, snrcut=snrcut

    x = atlasdust[cut[indx]].b_lum_obs
    xerr = atlasdust[cut[indx]].b_lum_obs_err
    xabs = atlasdust[cut[indx]].m_b_obs

    if keyword_set(nfgs) then begin
       cutnfgs = where(nfgsdust.b_lum_obs gt -900.0)
       lineratio, nfgsnodust[cutnfgs], '', '', 'OII_3727', 'H_ALPHA', $
         dum1, dum2, ynfgs, yerrnfgs, index=indxnfgs, snrcut=snrcut
       xnfgs = nfgsdust[cutnfgs[indxnfgs]].b_lum_obs
       xerrnfgs = nfgsdust[cutnfgs[indxnfgs]].b_lum_obs_err
    endif

    ystats = im_stats(y)

    xtitle = 'log L(B) [L(B)'+sunsymbol()+']'
    ytitle = 'log ([O II]/H\alpha)_{cor}'

    xrange = LBrange
    yrange = oiiharange

    atlas1d_lineplot, x, y, xerr, yerr, plottype=3, postscript=postscript, $
      xtitle=xtitle, ytitle=ytitle, xrange=xrange, yrange=yrange, legendtype=1, $
      xstyle=11, ymargin=[4,3], /right, /top, xnfgs=xnfgs, ynfgs=ynfgs, xerrnfgs=xerrnfgs, $
      yerrnfgs=yerrnfgs
    axis, /xaxis, xrange=interpol(xabs,x,!x.crange), xthick=postthick, charthick=postthick, $
      charsize=2.0, xtitle = textoidl('M_{B} [mag]'), xsty=1

; overplot histograms
    
    w1 = where(xabs ge -19.2)
    w2 = where(xabs lt -19.2)
 
; bin 1    
    
    im_plothist, y[w1], xbin, ybin, bin=0.2, /noplot, /fraction, /halfbin
;   im_plothist, y[w1], bin=0.2, xsty=3, ysty=1, charsize=1.0, charthick=postthick, /halfbin, $
;     /noerase, xthick=postthick, ythick=postthick, position=[0.32,0.23,0.6,0.48], $
;     /normal, ytitle='Fraction', xtitle=textoidl(ytitle), thick=postthick, $
;     /fraction, yrange=minmax(ybin)*[0,1.1], $
;     /fill, forientation=45, fcolor=djs_icolor('grey')
    
; bin 2
 
    im_plothist, y[w2], xbin, ybin, bin=0.2, /noplot, /fraction, /halfbin
;   im_plothist, y[w2], bin=0.2, /overplot, thick=postthick, /fraction, line=2, /halfbin

    im_openclose, postscript=postscript, /close    

; ------------------------------------------------------------
; L(B) vs [O II]_obs/Ha_cor
; ------------------------------------------------------------

    psname = 'LB_vs_oiiobs_hacor'
    im_openclose, pspath+psname, postscript=postscript

    indx = where((atlasdust.b_lum_obs gt -900.0) and (atlasdust.oii_3727[0]/atlasdust.oii_3727[1] gt snrcut) and $
      (atlasnodust.h_alpha[0]/atlasnodust.h_alpha[1] gt snrcut),nindx)

    x = atlasdust[indx].b_lum_obs
    xerr = atlasdust[indx].b_lum_obs_err
    xabs = atlasdust[indx].m_b_obs

    y = alog10(atlasdust[indx].oii_3727[0]/atlasnodust[indx].h_alpha[0])
    yerr = im_compute_error(atlasdust[indx].oii_3727[0],atlasdust[indx].oii_3727[1],$
      atlasnodust[indx].h_alpha[0],atlasnodust[indx].h_alpha[1],/log)
    
    if keyword_set(nfgs) then begin

       indxnfgs = where((nfgsdust.b_lum_obs gt -900.0) (nfgsdust.oii_3727[0]/nfgsdust.oii_3727[1] gt snrcut) and $
         (nfgsnodust.h_alpha[0]/nfgsnodust.h_alpha[1] gt snrcut),nindx)

       xnfgs = nfgsdust[indxnfgs].b_lum_obs
       xnfgserr = nfgsdust[indxnfgs].b_lum_obs_err

       ynfgs = alog10(nfgsdust[indxnfgs].oii_3727[0]/nfgsnodust[indxnfgs].h_alpha[0])
       yerrnfgs = im_compute_error(nfgsdust[indxnfgs].oii_3727[0],nfgsdust[indxnfgs].oii_3727[1],$
         nfgsnodust[indxnfgs].h_alpha[0],nfgsnodust[indxnfgs].h_alpha[1],/log)

    endif

    xtitle = 'log L(B) [L(B)'+sunsymbol()+']'
    ytitle = 'log ([O II]_{obs}/H\alpha_{cor})'

    xrange = LBrange
    yrange = oiihacorrange

    atlas1d_lineplot, x, y, xerr, yerr, plottype=3, postscript=postscript, $
      xtitle=xtitle, ytitle=ytitle, xrange=xrange, yrange=yrange, legendtype=1, $
      xstyle=11, ymargin=[4,3], /right, /top, xnfgs=xnfgs, ynfgs=ynfgs, xerrnfgs=xerrnfgs, $
      yerrnfgs=yerrnfgs
    axis, /xaxis, xrange=interpol(xabs,x,!x.crange), xthick=postthick, charthick=postthick, $
      charsize=2.0, xtitle = textoidl('M_{B} [mag]'), xsty=1

; overplot histograms
    
    w1 = where(xabs ge -19.2)
    w2 = where(xabs lt -19.2)

; bin 1    

    im_plothist, y[w1], xbin, ybin, bin=0.2, /noplot, /fraction, /halfbin
;   im_plothist, y[w1], bin=0.2, xsty=3, ysty=1, charsize=1.0, charthick=postthick, /halfbin, $
;     /noerase, xthick=postthick, ythick=postthick, position=[0.32,0.23,0.6,0.48], $
;     /normal, ytitle='Fraction', xtitle=textoidl(ytitle), thick=postthick, $
;     /fraction, yrange=minmax(ybin)*[0,1.1], $
;     /fill, forientation=45, fcolor=djs_icolor('grey')
    
; bin 2

    im_plothist, y[w2], xbin, ybin, bin=0.2, /noplot, /fraction, /halfbin
;   im_plothist, y[w2], bin=0.2, /overplot, thick=postthick, /fraction, line=2, /halfbin

    im_openclose, postscript=postscript, /close    

; ------------------------------------------------------------
; L_Ks vs [O II]/Ha_obs
; ------------------------------------------------------------
    
    psname = 'L_Ks_vs_oiiha_obs'
    im_openclose, pspath+psname, postscript=postscript

    cut = where(atlasdust.twomass_Ks_lum gt -900.0)
    lineratio, atlasdust[cut], '', '', 'OII_3727', 'H_ALPHA', $
      dum1, dum2, y, yerr, index=indx, nindex=nindx, snrcut=snrcut

    x = atlasdust[cut[indx]].twomass_Ks_lum
    xerr = atlasdust[cut[indx]].twomass_Ks_lum_err
    xabs = atlasdust[cut[indx]].twomass_M_Ks

    if keyword_set(nfgs) then begin
       cutnfgs = where(nfgsdust.twomass_Ks_lum gt -900.0)
       lineratio, nfgsdust[cutnfgs], '', '', 'OII_3727', 'H_ALPHA', $
         dum1, dum2, ynfgs, yerrnfgs, index=indxnfgs, snrcut=snrcut
       xnfgs = nfgsdust[cutnfgs[indxnfgs]].twomass_Ks_lum
       xerrnfgs = nfgsdust[cutnfgs[indxnfgs]].twomass_Ks_lum_err
    endif

    xtitle = 'log L_{K_{s}} [L'+sunsymbol()+']'
    ytitle = 'log ([O II]/H\alpha)_{obs}'

    xrange = LKsrange
    yrange = oiiharange

    atlas1d_lineplot, x, y, xerr, yerr, plottype=3, postscript=postscript, $
      xtitle=xtitle, ytitle=ytitle, xrange=xrange, yrange=yrange, legendtype=1, $
      xstyle=11, ymargin=[4,3], /right, /top, xnfgs=xnfgs, ynfgs=ynfgs, $
      xerrnfgs=xerrnfgs, yerrnfgs=yerrnfgs
    axis, /xaxis, xrange=interpol(xabs,x,!x.crange), xthick=postthick, $
      charthick=postthick, charsize=2.0, xtitle=textoidl('M_{Ks}'), xsty=1

    im_openclose, postscript=postscript, /close    

; ------------------------------------------------------------
; L_Ks vs [O II]/Ha_cor
; ------------------------------------------------------------
    
    psname = 'L_Ks_vs_oiiha_cor'
    im_openclose, pspath+psname, postscript=postscript

    cut = where(atlasdust.twomass_Ks_lum gt -900.0)
    lineratio, atlasnodust[cut], '', '', 'OII_3727', 'H_ALPHA', $
      dum1, dum2, y, yerr, index=indx, nindex=nindx, snrcut=snrcut

    x = atlasdust[cut[indx]].twomass_Ks_lum
    xerr = atlasdust[cut[indx]].twomass_Ks_lum_err
    xabs = atlasdust[cut[indx]].twomass_M_Ks

    if keyword_set(nfgs) then begin
       cutnfgs = where(nfgsdust.twomass_Ks_lum gt -900.0)
       lineratio, nfgsnodust[cutnfgs], '', '', 'OII_3727', 'H_ALPHA', $
         dum1, dum2, ynfgs, yerrnfgs, index=indxnfgs, snrcut=snrcut
       xnfgs = nfgsdust[cutnfgs[indxnfgs]].twomass_Ks_lum
       xerrnfgs = nfgsdust[cutnfgs[indxnfgs]].twomass_Ks_lum_err
    endif

    xtitle = 'log L_{K_{s}} [L'+sunsymbol()+']'
    ytitle = 'log ([O II]/H\alpha)_{cor}'

    xrange = LKsrange
    yrange = oiiharange

    atlas1d_lineplot, x, y, xerr, yerr, plottype=3, postscript=postscript, $
      xtitle=xtitle, ytitle=ytitle, xrange=xrange, yrange=yrange, legendtype=1, $
      xstyle=11, ymargin=[4,3], /right, /top, xnfgs=xnfgs, ynfgs=ynfgs, $
      xerrnfgs=xerrnfgs, yerrnfgs=yerrnfgs
    axis, /xaxis, xrange=interpol(xabs,x,!x.crange), xthick=postthick, $
      charthick=postthick, charsize=2.0, xtitle=textoidl('M_{Ks}'), xsty=1

    im_openclose, postscript=postscript, /close    

; ------------------------------------------------------------
; L_Ks vs [O II]_obs/Ha_cor
; ------------------------------------------------------------

    psname = 'L_Ks_vs_oiiobs_hacor'
    im_openclose, pspath+psname, postscript=postscript

    indx = where((atlasdust.twomass_Ks_lum gt -900.0) and (atlasdust.oii_3727[0]/atlasdust.oii_3727[1] gt snrcut) and $
      (atlasnodust.h_alpha[0]/atlasnodust.h_alpha[1] gt snrcut),nindx)

    x = atlasdust[indx].twomass_Ks_lum
    xerr = atlasdust[indx].twomass_Ks_lum_err
    xabs = atlasdust[indx].twomass_M_Ks

    y = alog10(atlasdust[indx].oii_3727[0]/atlasnodust[indx].h_alpha[0])
    yerr = im_compute_error(atlasdust[indx].oii_3727[0],atlasdust[indx].oii_3727[1],$
      atlasnodust[indx].h_alpha[0],atlasnodust[indx].h_alpha[1],/log)
    
    if keyword_set(nfgs) then begin
       indxnfgs = where((nfgsdust.twomass_Ks_lum gt -900.0) and (nfgsdust.oii_3727[0]/nfgsdust.oii_3727[1] gt snrcut) and $
         (nfgsnodust.h_alpha[0]/nfgsnodust.h_alpha[1] gt snrcut),nindx)
       xnfgs = nfgsdust[indxnfgs].twomass_Ks_lum
       xnfgserr = nfgsdust[indxnfgs].twomass_Ks_lum_err
       ynfgs = alog10(nfgsdust[indxnfgs].oii_3727[0]/nfgsnodust[indxnfgs].h_alpha[0])
       yerrnfgs = im_compute_error(nfgsdust[indxnfgs].oii_3727[0],nfgsdust[indxnfgs].oii_3727[1],$
         nfgsnodust[indxnfgs].h_alpha[0],nfgsnodust[indxnfgs].h_alpha[1],/log)
    endif

    xtitle = 'log L_{K_{s}} [L'+sunsymbol()+']'
    ytitle = 'log ([O II]_{obs}/H\alpha_{cor})'

    xrange = LKsrange
    yrange = oiihacorrange

    atlas1d_lineplot, x, y, xerr, yerr, plottype=3, postscript=postscript, $
      xtitle=xtitle, ytitle=ytitle, xrange=xrange, yrange=yrange, legendtype=1, $
      xstyle=11, ymargin=[4,3], /right, /top, xnfgs=xnfgs, ynfgs=ynfgs, $
      xerrnfgs=xerrnfgs, yerrnfgs=yerrnfgs
    axis, /xaxis, xrange=interpol(xabs,x,!x.crange), xthick=postthick, $
      charthick=postthick, charsize=2.0, xtitle=textoidl('M_{Ks}'), xsty=1

    im_openclose, postscript=postscript, /close    

; ------------------------------------------------------------
; Mass vs [O II]/Ha_obs
; ------------------------------------------------------------
    
    psname = 'mass_vs_oiiha_obs'
    im_openclose, pspath+psname, postscript=postscript

    cut = where(atlasdust.mass_gr_r gt -900.0)
    lineratio, atlasdust[cut], '', '', 'OII_3727', 'H_ALPHA', $
      dum1, dum2, y, yerr, index=indx, nindex=nindx, snrcut=snrcut

    x = atlasdust[cut[indx]].mass_gr_r
    xerr = atlasdust[cut[indx]].mass_gr_r_err

    if keyword_set(nfgs) then begin
       cutnfgs = where(nfgsdust.mass_gr_r gt -900.0)
       lineratio, nfgsdust[cutnfgs], '', '', 'OII_3727', 'H_ALPHA', $
         dum1, dum2, ynfgs, yerrnfgs, index=indxnfgs, snrcut=snrcut
       xnfgs = nfgsdust[cutnfgs[indxnfgs]].mass_gr_r
       xerrnfgs = nfgsdust[cutnfgs[indxnfgs]].mass_gr_r_err
    endif

    xtitle = 'log (M/M'+sunsymbol()+')'
    ytitle = 'log ([O II]/H\alpha)_{obs}'

    xrange = massrange
    yrange = oiiharange

    atlas1d_lineplot, x, y, xerr, yerr, plottype=3, postscript=postscript, $
      xtitle=xtitle, ytitle=ytitle, xrange=xrange, yrange=yrange, legendtype=1, $
      /right, /top, xnfgs=xnfgs, ynfgs=ynfgs, xerrnfgs=xerrnfgs, yerrnfgs=yerrnfgs

    im_openclose, postscript=postscript, /close    

; ------------------------------------------------------------
; Mass vs [O II]/Ha_cor
; ------------------------------------------------------------
    
    psname = 'mass_vs_oiiha_cor'
    im_openclose, pspath+psname, postscript=postscript

    cut = where(atlasdust.mass_gr_r gt -900.0)
    lineratio, atlasnodust[cut], '', '', 'OII_3727', 'H_ALPHA', $
      dum1, dum2, y, yerr, index=indx, nindex=nindx, snrcut=snrcut

    x = atlasdust[cut[indx]].mass_gr_r
    xerr = atlasdust[cut[indx]].mass_gr_r_err

    if keyword_set(nfgs) then begin
       cutnfgs = where(nfgsdust.mass_gr_r gt -900.0)
       lineratio, nfgsnodust[cutnfgs], '', '', 'OII_3727', 'H_ALPHA', $
         dum1, dum2, ynfgs, yerrnfgs, index=indxnfgs, snrcut=snrcut
       xnfgs = nfgsdust[cutnfgs[indxnfgs]].mass_gr_r
       xerrnfgs = nfgsdust[cutnfgs[indxnfgs]].mass_gr_r_err
    endif

    xtitle = 'log M [M'+sunsymbol()+']'
    ytitle = 'log ([O II]/H\alpha)_{cor}'

    xrange = massrange
    yrange = oiiharange

    atlas1d_lineplot, x, y, xerr, yerr, plottype=3, postscript=postscript, $
      xtitle=xtitle, ytitle=ytitle, xrange=xrange, yrange=yrange, legendtype=1, $
      /right, /top, xnfgs=xnfgs, ynfgs=ynfgs, xerrnfgs=xerrnfgs, yerrnfgs=yerrnfgs

    im_openclose, postscript=postscript, /close    

; ------------------------------------------------------------
; L_Ks vs [O II]_obs/Ha_cor
; ------------------------------------------------------------

    psname = 'L_Ks_vs_oiiobs_hacor'
    im_openclose, pspath+psname, postscript=postscript

    indx = where((atlasdust.mass_gr_r gt -900.0) and (atlasdust.oii_3727[0]/atlasdust.oii_3727[1] gt snrcut) and $
      (atlasnodust.h_alpha[0]/atlasnodust.h_alpha[1] gt snrcut),nindx)

    x = atlasdust[indx].mass_gr_r
    xerr = atlasdust[indx].mass_gr_r_err

    y = alog10(atlasdust[indx].oii_3727[0]/atlasnodust[indx].h_alpha[0])
    yerr = im_compute_error(atlasdust[indx].oii_3727[0],atlasdust[indx].oii_3727[1],$
      atlasnodust[indx].h_alpha[0],atlasnodust[indx].h_alpha[1],/log)
    
    if keyword_set(nfgs) then begin
       indxnfgs = where((nfgsdust.mass_gr_r gt -900.0) and (nfgsdust.oii_3727[0]/nfgsdust.oii_3727[1] gt snrcut) and $
         (nfgsnodust.h_alpha[0]/nfgsnodust.h_alpha[1] gt snrcut),nindx)
       xnfgs = nfgsdust[indxnfgs].mass_gr_r
       xnfgserr = nfgsdust[indxnfgs].mass_gr_r_err
       ynfgs = alog10(nfgsdust[indxnfgs].oii_3727[0]/nfgsnodust[indxnfgs].h_alpha[0])
       yerrnfgs = im_compute_error(nfgsdust[indxnfgs].oii_3727[0],nfgsdust[indxnfgs].oii_3727[1],$
         nfgsnodust[indxnfgs].h_alpha[0],nfgsnodust[indxnfgs].h_alpha[1],/log)
    endif

    xtitle = 'log M [M'+sunsymbol()+']'
    ytitle = 'log ([O II]_{obs}/H\alpha_{cor})'

    xrange = massrange
    yrange = oiihacorrange

    atlas1d_lineplot, x, y, xerr, yerr, plottype=3, postscript=postscript, $
      xtitle=xtitle, ytitle=ytitle, xrange=xrange, yrange=yrange, legendtype=1, $
      /right, /top, xnfgs=xnfgs, ynfgs=ynfgs, xerrnfgs=xerrnfgs, yerrnfgs=yerrnfgs

    im_openclose, postscript=postscript, /close    

; ------------------------------------------------------------
; J-Ks vs [O II]/Ha_cor
; ------------------------------------------------------------

    psname = 'JKs_vs_oiiha_cor'
    im_openclose, pspath+psname, postscript=postscript

    cut = where((atlasdust.twomass_J gt -900.0) and (atlasdust.twomass_Ks gt -900.0))
    lineratio, atlasnodust[cut], '', '', 'OII_3727', 'H_ALPHA', $
      dum1, dum2, y, yerr, index=indx, nindex=nindx, snrcut=snrcut

    x = atlasdust[cut[indx]].twomass_J - atlasdust[cut[indx]].twomass_Ks
    xerr = sqrt(atlasdust[cut[indx]].twomass_J_err^2 + atlasdust[cut[indx]].twomass_Ks_err^2)
    
    if keyword_set(nfgs) then begin
       cutnfgs = where((nfgsdust.twomass_J gt -900.0) and (nfgsdust.twomass_Ks gt -900.0))
       lineratio, nfgsnodust[cutnfgs], '', '', 'OII_3727', 'H_ALPHA', $
         dum1, dum2, ynfgs, yerrnfgs, index=indxnfgs, snrcut=snrcut
       xnfgs = nfgsdust[cutnfgs[indxnfgs]].twomass_J - nfgsdust[cutnfgs[indxnfgs]].twomass_Ks
       xerrnfgs = sqrt(nfgsdust[cutnfgs[indxnfgs]].twomass_J_err^2.0 + nfgsdust[cutnfgs[indxnfgs]].twomass_Ks_err^2.0)
    endif

    xtitle = 'J - K_{s}' 
    ytitle = 'log ([O II]/H\alpha)_{cor}'

    xrange = JKsrange
    yrange = oiiharange

    atlas1d_lineplot, x, y, xerr, yerr, plottype=3, postscript=postscript, $
      xtitle=xtitle, ytitle=ytitle, xrange=xrange, yrange=yrange, legendtype=1, $
      /right, /top, xnfgs=xnfgs, ynfgs=ynfgs, xerrnfgs=xerrnfgs, yerrnfgs=yerrnfgs
    
    im_openclose, postscript=postscript, /close    
    
; ------------------------------------------------------------
; EW([O II]) vs [O II]/Ha_cor
; ------------------------------------------------------------

    psname = 'ewoii_vs_oiiha_cor'
    im_openclose, pspath+psname, postscript=postscript

    lineratio, atlasnodust, 'OII_3727_EW', '', 'OII_3727', 'H_ALPHA', $
      x, xerr, y, yerr, index=indx, nindex=nindx, snrcut=snrcut

    if keyword_set(kenn92) then begin
       lineratio, kenn92nodust, 'OII_3727_EW', '', 'OII_3727', 'H_ALPHA', $
         xkenn92, xerrkenn92, ykenn92, yerrkenn92, index=indxkenn92, snrcut=snrcut
    endif

    if keyword_set(nfgs) then begin
       lineratio, nfgsnodust, 'OII_3727_EW', '', 'OII_3727', 'H_ALPHA', $
         xnfgs, xerrnfgs, ynfgs, yerrnfgs, index=indxnfgs, snrcut=snrcut
    endif

    xtitle = 'log EW([O II])'
    ytitle = 'log ([O II]/H\alpha)_{cor}'

    xrange = ewoiirange
    yrange = oiiharange

    atlas1d_lineplot, x, y, xerr, yerr, plottype=3, postscript=postscript, $
      xtitle=xtitle, ytitle=ytitle, xrange=xrange, yrange=yrange, legendtype=1, $
      /right, /top, $
;     xkenn92=xkenn92, ykenn92=ykenn92, xerrkenn92=xerrkenn92, yerrkenn92=yerrkenn92, $
      xnfgs=xnfgs, ynfgs=ynfgs, xerrnfgs=xerrnfgs, yerrnfgs=yerrnfgs

    im_openclose, postscript=postscript, /close    

; ------------------------------------------------------------
; EW([O II]) vs [O II]_obs/Ha_cor
; ------------------------------------------------------------

    psname = 'ewoii_vs_oiiobs_hacor'
    im_openclose, pspath+psname, postscript=postscript

    lineratio, atlasnodust, 'OII_3727_EW', '', 'OII_3727', 'H_ALPHA', $
      x, xerr, dum1, dum2, index=indx, nindex=nindx, snrcut=snrcut

    y = alog10(atlasdust[indx].oii_3727[0]/atlasnodust[indx].h_alpha[0])
    yerr = im_compute_error(atlasdust[indx].oii_3727[0],atlasdust[indx].oii_3727[1],$
      atlasnodust[indx].h_alpha[0],atlasnodust[indx].h_alpha[1],/log)
    
    if keyword_set(kenn92) then begin

       lineratio, kenn92nodust, 'OII_3727_EW', '', 'OII_3727', 'H_ALPHA', $
         xkenn92, xerrkenn92, dum1, dum2, index=indxkenn92, snrcut=snrcut
       ykenn92 = alog10(kenn92dust[indxkenn92].oii_3727[0]/kenn92nodust[indxkenn92].h_alpha[0])
       yerrkenn92 = im_compute_error(kenn92dust[indxkenn92].oii_3727[0],kenn92nodust[indxkenn92].oii_3727[1],$
         kenn92nodust[indxkenn92].h_alpha[0],kenn92nodust[indxkenn92].h_alpha[1],/log)
       
    endif

    if keyword_set(nfgs) then begin

       lineratio, nfgsnodust, 'OII_3727_EW', '', 'OII_3727', 'H_ALPHA', $
         xnfgs, xerrnfgs, dum1, dum2, index=indxnfgs, snrcut=snrcut
       ynfgs = alog10(nfgsdust[indxnfgs].oii_3727[0]/nfgsnodust[indxnfgs].h_alpha[0])
       yerrnfgs = im_compute_error(nfgsdust[indxnfgs].oii_3727[0],nfgsnodust[indxnfgs].oii_3727[1],$
         nfgsnodust[indxnfgs].h_alpha[0],nfgsnodust[indxnfgs].h_alpha[1],/log)
       
    endif

    xtitle = 'log EW([O II])'
    ytitle = 'log ([O II]_{obs}/H\alpha_{cor})'

    xrange = ewoiirange
    yrange = oiihacorrange

    atlas1d_lineplot, x, y, xerr, yerr, plottype=3, postscript=postscript, $
      xtitle=xtitle, ytitle=ytitle, xrange=xrange, yrange=yrange, legendtype=1, $
      /right, /top, $
;     xkenn92=xkenn92, ykenn92=ykenn92, xerrkenn92=xerrkenn92, yerrkenn92=yerrkenn92, $
      xnfgs=xnfgs, ynfgs=ynfgs, xerrnfgs=xerrnfgs, yerrnfgs=yerrnfgs

    sixlin, [x,xnfgs], [y,ynfgs], a, siga, b, sigb ; Ordinary Least Squares Bisector
    coeff = [a[2],b[2]]

;   coeff = linfit([x,xnfgs],[y,ynfgs])
    
    xaxis = findgen((2.0-0.6)/0.05+1)*0.05+0.6
;   djs_oplot, xaxis, poly(xaxis,coeff), line=2, thick=postthick
    
    im_openclose, postscript=postscript, /close    

; --------------------------------------------------    
; GENERATE THE WEB PAGE
; --------------------------------------------------    
    
    if keyword_set(postscript) then begin

       im_ps2html, htmlbase, html_path=html_path, cleanpng=0, npscols=2, _extra=extra
       
    endif

; ------------------------------------------------------------
; Abundances vs [O II]_cor/Ha_cor
; ------------------------------------------------------------

    psname = 'abundances_vs_oiicor_hacor'
    im_openclose, pspath+psname, postscript=postscript

    plotsym, 0, 1.2

    pagemaker, nx=2, ny=2, position=pos, /normal, xspace=0.0, $
      xmargin=[1.5,0.1], ymargin=[0.5,1.2], yspace=0.0

    ytitle = 'log ([O II]/H\alpha)_{cor}'
    xtitle = '12 + log (O/H)'

    xrange = abundrange
    yrange = oiiharange
    pcharsize = 1.4
    lcharsize = 1.4

; --------------------------------------------------    
;   cut = where(atlasnodust.zstrong_12oh_P01_melbourne gt -900.0)
    cut = where(atlasnodust.zstrong_12oh_CL01A gt -900.0)
    lineratio, atlasnodust[cut], '', '', 'OII_3727', 'H_ALPHA', $
      dum1, dum2, y, yerr, index=indx, nindex=nindx, snrcut=snrcut
;   x = atlasnodust[cut[indx]].zstrong_12oh_P01_melbourne
    x = atlasnodust[cut[indx]].zstrong_12oh_CL01A

    djs_plot, [0], [0], xrange=xrange, yrange=yrange, xtickname=replicate(' ',10), $
      charsize=pcharsize, charthick=postthick, xthick=postthick, ythick=postthick, $
      xsty=3, ysty=3, position=pos[*,0], ytitle=ytitle
    djs_oplot, x, y, ps=8

    legend, 'CL01', /left, /bottom, box=0, charsize=lcharsize, charthick=postthick
;   legend, 'P01', /left, /bottom, box=0, charsize=lcharsize, charthick=postthick
; --------------------------------------------------    
    cut = where(atlasnodust.zstrong_12oh_KD02_NIIOII gt -900.0)
    lineratio, atlasnodust[cut], '', '', 'OII_3727', 'H_ALPHA', $
      dum1, dum2, y, yerr, index=indx, nindex=nindx, snrcut=snrcut
    x = atlasnodust[cut[indx]].zstrong_12oh_KD02_NIIOII

    djs_plot, [0], [0], xrange=xrange, yrange=yrange, charsize=pcharsize, charthick=postthick, $
      xthick=postthick, ythick=postthick, xsty=3, ysty=3, xtickname=replicate(' ',10), $
      ytickname=replicate(' ',10), position=pos[*,1], /noerase
    djs_oplot, x, y, ps=8

    legend, 'KD02 - [N II]/[O II]', /left, /bottom, box=0, charsize=lcharsize, charthick=postthick

; --------------------------------------------------    
    cut = where(atlasnodust.zstrong_12oh_M91_CONTINI gt -900.0)
    lineratio, atlasnodust[cut], '', '', 'OII_3727', 'H_ALPHA', $
      dum1, dum2, y, yerr, index=indx, nindex=nindx, snrcut=snrcut
    x = atlasnodust[cut[indx]].zstrong_12oh_M91_CONTINI

    djs_plot, [0], [0], xrange=xrange, yrange=yrange, xtitle=xtitle, charsize=pcharsize, $
      charthick=postthick, xthick=postthick, ythick=postthick, xsty=3, ysty=3, $
      position=pos[*,2], /noerase, ytitle=ytitle
    djs_oplot, x, y, ps=8
    legend, 'M91 - Contini', /left, /bottom, box=0, charsize=lcharsize, charthick=postthick

; --------------------------------------------------    
    cut = where(atlasnodust.zstrong_12oh_ZKH94 gt -900.0)
    lineratio, atlasnodust[cut], '', '', 'OII_3727', 'H_ALPHA', $
      dum1, dum2, y, yerr, index=indx, nindex=nindx, snrcut=snrcut
    x = atlasnodust[cut[indx]].zstrong_12oh_ZKH94

    djs_plot, [0], [0], xrange=xrange, yrange=yrange, xtitle=xtitle, $
      charsize=pcharsize, charthick=postthick, xthick=postthick, ythick=postthick, xsty=3, ysty=3, $
      position=pos[*,3], /noerase, ytickname=replicate(' ',10)
    djs_oplot, x, y, ps=8
    legend, 'ZKH94', /left, /bottom, box=0, charsize=lcharsize, charthick=postthick
; --------------------------------------------------    

    im_openclose, postscript=postscript, /close    

; ---------------------------------------------------------------------------    
; 12+log (O/H) - CL01 vs [O II]_cor/Ha_cor
; ---------------------------------------------------------------------------    

    psname = '12oh_CL01_vs_oiicor_hacor'
    im_openclose, pspath+psname, postscript=postscript

    cut = where(atlasnodust.zstrong_12oh_CL01A gt -900.0)
    lineratio, atlasnodust[cut], '', '', 'OII_3727', 'H_ALPHA', $
      dum1, dum2, y, yerr, index=indx, nindex=nindx, snrcut=snrcut

    x = atlasnodust[cut[indx]].zstrong_12oh_CL01A
;   xerr = atlasnodust[cut[indx]].zstrong_12oh_CL01A_err
    xerr = x*0.0+0.1

    if keyword_set(nfgs) then begin
       cutnfgs = where(nfgsnodust.zstrong_12oh_CL01A gt -900.0)
       lineratio, nfgsnodust[cutnfgs], '', '', 'OII_3727', 'H_ALPHA', $
         dum1, dum2, ynfgs, yerrnfgs, index=indxnfgs, snrcut=snrcut
       xnfgs = nfgsnodust[cutnfgs[indxnfgs]].zstrong_12oh_CL01A
       xerrnfgs = xnfgs*0.0+0.1
;      xerrnfgs = nfgsnodust[cutnfgs[indxnfgs]].zstrong_12oh_CL01A_err
    endif
    
    if keyword_set(hiiregions) then begin
       good = where((hii.oii_h_alpha gt -900.0) and (hii.zstrong_12oh_CL01A gt -900.0))
       xregion = hii[good].zstrong_12oh_CL01A & xregion = xregion*0.0
       yregion = hii[good].oii_h_alpha & yerrregion = hii[good].oii_h_alpha_err
    endif

    ytitle = 'log ([O II]/H\alpha)_{cor}'
    xtitle = '12 + log (O/H) [CL01 - Case A]'

    xrange = abundrange
    yrange = oiiharange
    
    atlas1d_lineplot, x, y, xerr, yerr, plottype=3, postscript=postscript, $
      xtitle=xtitle, ytitle=ytitle, xrange=xrange, yrange=yrange, legendtype=1, $
      /right, /top, xregion=xregion, yregion=yregion, xerrregion=xerrregion, $
      yerrregion=yerrregion, xnfgs=xnfgs, ynfgs=ynfgs, xerrnfgs=xerrnfgs, yerrnfgs=yerrnfgs

; overplot the Kewley result using ZKH94

    ohaxis = findgen((9.4-8.4)/0.01)*0.01+8.4 ; only applicable over this range!
    kewley_oiiha = alog10(16.73-1.75*ohaxis - (Zsun_old-Zsun_new))
;   kewley_oiiha = alog10(16.73-1.75*ohaxis)

    djs_oplot, ohaxis, kewley_oiiha, line=2, thick=postthick
    legend, 'KGJ04', /left, /top, box=0, charsize=charsize, $
      charthick=postthick, line=2, thick=postthick

    im_openclose, postscript=postscript, /close    

; ---------------------------------------------------------------------------    
; 12+log (O/H) - ZKH94 vs [O II]_cor/Ha_cor
; ---------------------------------------------------------------------------    

    psname = '12oh_ZKH94_vs_oiicor_hacor'
    im_openclose, pspath+psname, postscript=postscript

    cut = where(atlasnodust.zstrong_12oh_ZKH94 gt -900.0)
    lineratio, atlasnodust[cut], '', '', 'OII_3727', 'H_ALPHA', $
      dum1, dum2, y, yerr, index=indx, nindex=nindx, snrcut=snrcut

    x = atlasnodust[cut[indx]].zstrong_12oh_ZKH94
;   xerr = atlasnodust[cut[indx]].zstrong_12oh_ZKH94_err
    xerr = x*0.0+0.1

    if keyword_set(nfgs) then begin
       cutnfgs = where(nfgsnodust.zstrong_12oh_ZKH94 gt -900.0)
       lineratio, nfgsnodust[cutnfgs], '', '', 'OII_3727', 'H_ALPHA', $
         dum1, dum2, ynfgs, yerrnfgs, index=indxnfgs, snrcut=snrcut
       xnfgs = nfgsnodust[cutnfgs[indxnfgs]].zstrong_12oh_ZKH94
       xerrnfgs = xnfgs*0.0+0.1
;      xerrnfgs = nfgsnodust[cutnfgs[indx]].zstrong_12oh_ZKH94_err
    endif
    
    if keyword_set(hiiregions) then begin
       good = where((hii.oii_h_alpha gt -900.0) and (hii.zstrong_12oh_ZKH94 gt -900.0))
       xregion = hii[good].zstrong_12oh_ZKH94 & xregion = xregion*0.0
       yregion = hii[good].oii_h_alpha & yerrregion = hii[good].oii_h_alpha_err
    endif

    ytitle = 'log ([O II]/H\alpha)_{cor}'
    xtitle = '12 + log (O/H) [ZKH94]'

    xrange = abundrange
    yrange = oiiharange
    
    atlas1d_lineplot, x, y, xerr, yerr, plottype=3, postscript=postscript, $
      xtitle=xtitle, ytitle=ytitle, xrange=xrange, yrange=yrange, legendtype=1, $
      /right, /top, xregion=xregion, yregion=yregion, xerrregion=xerrregion, $
      yerrregion=yerrregion, xnfgs=xnfgs, ynfgs=ynfgs, xerrnfgs=xerrnfgs, yerrnfgs=yerrnfgs

; overplot the Kewley result using ZKH94

    ohaxis = findgen((9.4-8.4)/0.01)*0.01+8.4 ; only applicable over this range!
    kewley_oiiha = alog10(16.73-1.75*ohaxis - (Zsun_old-Zsun_new))
;   kewley_oiiha = alog10(16.73-1.75*ohaxis)

    djs_oplot, ohaxis, kewley_oiiha, line=2, thick=postthick

    legend, 'KGJ04', /left, /top, box=0, charsize=charsize, $
      charthick=postthick, line=2, thick=postthick

    im_openclose, postscript=postscript, /close    

; ------------------------------------------------------------
; example spectra with various [O II]/Ha ratios
; ------------------------------------------------------------

    psname = 'spectra_oii_examples'
    im_openclose, pspath+psname, postscript=postscript
    
    spectra_oii_examples, atlasdust, atlasnodust, postscript=postscript

    im_openclose, postscript=postscript, /close    

;; ------------------------------------------------------------
;; [O II]/Ha histogram
;; ------------------------------------------------------------
;
;    psname = 'oii_histogram'
;    im_openclose, pspath+psname, postscript=postscript
;
;; exclude all galaxies with M_B > -M_B_cut
;    
;    cut = where((atlasdust.m_b_obs gt -900) and (atlasdust.m_b_obs lt M_B_cut))
;    lineratio, atlasdust[cut], 'OII_3727', 'H_ALPHA', '', '', x, xerr
;    lineratio, atlasnodust[cut], 'OII_3727', 'H_ALPHA', '', '', xnodust, xerrnodust
;
;    xstats = im_stats(x)
;    xnoduststats = im_stats(xnodust)
;
;    fmt = '(F6.3)'
;    xstr = ' ('+string(xstats.median,format=fmt)+' +/- '+$
;      string(xstats.sigma,format=fmt)+')'
;    xnoduststr = ' ('+string(xnoduststats.median,format=fmt)+' +/- '+$
;      string(xnoduststats.sigma,format=fmt)+')'
;
;; figure out the xyrange and plot it
;
;    im_plothist, x, xbin, ybin, bin=0.1, /noplot, /fraction, /halfbin
;    im_plothist, xnodust, xbindust, ybindust, bin=0.1, /noplot, /fraction, /halfbin
;    
;    yrange = [0.0,max(ybin)>max(ybindust)]*[1.0,1.2]
;    xrange = [-1.2,1.2] ; one point above this XRANGE
;
;    xtitle = textoidl('log [O II] - log H\alpha')
;    ytitle = 'Fraction'
;
;    im_plothist, x, bin=0.1, xsty=3, ysty=1, xthick=postthick, ythick=postthick, /halfbin, $
;      charsize=2.0, charthick=postthick, xrange=xrange, yrange=yrange, line=3, $
;      xtitle=xtitle, ytitle=ytitle, thick=postthick, /fraction
;    im_plothist, xnodust, bin=0.1, /overplot, color=djs_icolor('red'), /halfbin, $
;      thick=postthick, line=0, /fraction
;
;; legend
;
;    legend, ['Observed '+xstr,'Corrected '+xnoduststr], line=[3,0], $
;      color=djs_icolor(['default','red']), /left, /top, box=0, $
;      charsize=charsize, charthick=postthick, thick=4.0
;
;    if not keyword_set(postscript) then cc = get_kbrd(1)

;; ------------------------------------------------------------
;; redshift range of nebular emission lines
;; ------------------------------------------------------------
;
;    lines = ['[O II]','H\beta','[O III]','H\alpha']
;    waves = [3727.0,4861.0,5007.0,6563.0]
;    wavepos = [3727.0,4800.0,5250.0,6563.0]
;
;    zmin = 0.0
;    zmax = 2.0
;
;    wavemin = 3600.0
;    wavemax = 9900.0
;    xwavemin = 2500.0
;    xwavemax = 8500.0
;    
;    djs_plot, [0], [0], /nodata, xrange=[xwavemin,xwavemax], yrange=[zmin,zmax], $
;      xsty=1, ysty=1, xthick=postthick, ythick=postthick, $
;      charsize=2.0, charthick=postthick, xtitle='Wavelength ['+angstrom()+']', $
;      ytitle='Redshift', ymargin=[5,3]
;    for i = 0L, n_elements(lines)-1L do xyouts, wavepos[i], 2.07, textoidl(lines[i]), $
;      orientation=90, align=0.0, charsize=charsize, charthick=postthick, /data
;    for i = 0L, n_elements(lines)-1L do oplot, waves[i]*[1,1], !y.crange, $
;      line=2, thick=postthick
;
;    waveaxis = findgen((xwavemax-xwavemin)/10.0+1)*10.0+xwavemin
;
;    x = [3727.0,6563.0]
;
;    y1 = wavemax/x-1.0
;    y2 = wavemin/x-1.0
;
;    slope1 = (y1[1]-y1[0])/(x[1]-x[0])
;    intercept1 = y1[0]-slope1*x[0]
;    slope2 = (y2[1]-y2[0])/(x[1]-x[0])
;    intercept2 = y2[0]-slope2*x[0]
;    
;    line1 = poly(waveaxis,[intercept1,slope1])
;    line2 = poly(waveaxis,[intercept2,slope2])
;    oplot, waveaxis, line1, line=0, thick=postthick
;;   oplot, waveaxis, line2
;    
;;   polyfill, [xwavemin,xwavemin,wavemin,wavemax,interpol(waveaxis,line1,zmax)],$
;;     [zmax,interpol(line0,waveaxis,xwavemin),zmin,zmin,zmax], $
;;     /line_fill, color=djs_icolor('red'), linestyle=1
;    
