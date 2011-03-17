pro plot_sfrm_nuvmr, ps=ps
; jm10feb05ucsd - make various plots that rely on the NUV-r color

    if keyword_set(ps) then suffix = '.ps' else suffix = '.eps'

    sfrmpath = ages_path(/projects)+'sfrm/'
    paperpath = ages_path(/papers)+'sfrm/'

    parent = read_sfrm_sample()
    sdss = read_sfrm_sample(/sdss)
    zbins = sfrm_zbins(nzbins)

    levels = [0.01,0.05,0.1,0.3,0.6,0.9,0.95,0.99]
;   levels = [0.5,0.75,0.9,0.95]
;   cannotation = strtrim(levels,2)

; specify the models    
    modelspath = getenv('ISEDFIT_SFHGRID_DIR')+'/basemodels/bc03/'
    sspfits = 'chab_Z0.02_tau_00.0Gyr.fits.gz'
;   taufits = 'chab_Z0.004_tau_100Gyr.fits.gz'
;   taufits = 'chab_Z0.02_tau_07.0Gyr.fits.gz'
;   taufits = 'chab_Z0.004_tau_07.0Gyr.fits.gz'
    taufits = 'chab_Z0.02_tau_100Gyr.fits.gz'

    ssp = mrdfits(modelspath+sspfits,1,/silent)
    tau = mrdfits(modelspath+taufits,1,/silent)
    wave = ssp.wave

; choose the formation redshift and the ages to plot    
    minage = 0.1
    maxage = 13.5
    age = range(minage,maxage,10,/log)
    nmodel = n_elements(age)
    
; interpolate the models at the appropriate ages
    sspflux = interpolate(ssp.flux,findex(ssp.age/1D9,age))
    tauflux = interpolate(tau.flux,findex(tau.age/1D9,age))
    models = {$
      ssp_nuvmr: 0.0, ssp_rmj: 0.0, ssp_umv: 0.0, ssp_vmj: 0.0, $
      tau_nuvmr: 0.0, tau_rmj: 0.0, tau_umv: 0.0, tau_vmj: 0.0}
    models = replicate(models,nmodel)

    lambda = k_lambda_to_edges(wave)
    ff = ['galex_NUV','bessell_U','bessell_V','sdss_r0','twomass_J']+'.par'
    for ii = 0, nmodel-1 do begin
; ssp
       mm = reform(k_project_filters(lambda,sspflux[*,ii],$
         band_shift=0.1,filterlist=ff))
       models[ii].ssp_nuvmr = -2.5*alog10(mm[0]/mm[3])
       models[ii].ssp_rmj = -2.5*alog10(mm[3]/mm[4])
       models[ii].ssp_umv = -2.5*alog10(mm[1]/mm[2])
       models[ii].ssp_vmj = -2.5*alog10(mm[2]/mm[4])
; tau
       mm = reform(k_project_filters(lambda,tauflux[*,ii],$
         band_shift=0.1,filterlist=ff))
       models[ii].tau_nuvmr = -2.5*alog10(mm[0]/mm[3])
       models[ii].tau_rmj = -2.5*alog10(mm[3]/mm[4])
       models[ii].tau_umv = -2.5*alog10(mm[1]/mm[2])
       models[ii].tau_vmj = -2.5*alog10(mm[2]/mm[4])
    endfor

; --------------------------------------------------
; SDSS - NUV-r vs M_0.1r

    psfile = paperpath+'sdss_nuvmr_mr'+suffix
    im_plotconfig, 0, pos, psfile=psfile

    xrange = [-17.6,-23.8]
    yrange = [0.01,6.7]
    xtitle = textoidl('M_{0.1r}')
    ytitle = textoidl('^{0.1}(NUV - r)')

    these = where((total(sdss.k_maggies[0:1] gt 0.0,1) ge 1),nthese)

    mr = sdss[these].k_ugriz_absmag_01[2]
    nuvmr = sdss[these].k_galex_absmag_01[1]-sdss[these].k_ugriz_absmag_01[2]
    hogg_scatterplot, mr, nuvmr, position=pos, xsty=1, ysty=1, $
      yrange=yrange, xrange=xrange, xtitle=xtitle, ytitle=ytitle, $
      /outliers, outsymsize=0.1, outcolor=fsc_color('slate grey',101), $
      /internal, levels=levels, cannotation=cannotation, $
      smooth_contours=1.6
;   djs_oplot, -17.6*[1,1], !y.crange, line=5, thick=5

; selection of quiescent galaxies plus labels
;   oplot_select_quiescent
;   xyouts, 0.4, 5.8, 'Quiescent', align=0.5
;   xyouts, 1.4, 0.6, 'Star-Forming', align=0.5
    
    im_plotconfig, /psclose

; --------------------------------------------------
; NUV-r vs r-J for the SDSS + models

; #########################
; NUV-r vs r-J plot
    psfile = paperpath+'nuvmr_rmj_models'+suffix
    im_plotconfig, 0, pos, psfile=psfile

    xrange = [0.0,1.9]
    yrange = [0.01,6.7]
    xtitle = textoidl('^{0.1}(r - J)')
    ytitle = textoidl('^{0.1}(NUV - r)')
    
; sdss; require NUV and at least *one* of JHKs observed photometry 
    nuvmr = sdss.k_galex_absmag_01[1]-sdss.k_ugriz_absmag_01[2]
    rmj = sdss.k_ugriz_absmag_01[2]-sdss.k_ubvrijhk_absmag_01[5]
    gg = where((sdss.k_maggies[1] gt 0) and (total(sdss.k_maggies[7:9] gt 0.0,1) ge 1),ngg)
    rmj = rmj[gg] & nuvmr = nuvmr[gg]
    hogg_scatterplot, rmj, nuvmr, position=pos, xsty=1, ysty=1, $
      yrange=yrange, xrange=xrange, xtitle=xtitle, ytitle=ytitle, $
      /outliers, outsymsize=0.1, outcolor=fsc_color('slate grey',101), $
      /internal, levels=levels, cannotation=cannotation, $
      smooth_contours=1.6       ;, xnpix=30, ynpix=30

; overplot the models    
    djs_oplot, models.ssp_rmj, models.ssp_nuvmr, thick=8.0, $
      color='red', psym=-symcat(9,thick=8), symsize=2.0
    djs_oplot, models.tau_rmj, models.tau_nuvmr, thick=8.0, $
      color='blue', psym=-symcat(6,thick=8), symsize=2.0

; selection of quiescent galaxies plus labels
    oplot_select_quiescent
    xyouts, 0.4, 5.8, 'Quiescent', align=0.5
    xyouts, 1.4, 0.6, 'Star-Forming', align=0.5
    
; overlay a reddening vector
    nuvmr_true = 1.8
    rmj_true = 1.3

    ebv = 0.3
    nuvmr_red = nuvmr_true + 0.4*ebv*(k_lambda(k_lambda_eff(filterlist='galex_NUV.par',band_shift=0.0),/calz)-$
      k_lambda(k_lambda_eff(filterlist='sdss_r0.par',band_shift=0.0),/calz))
    rmj_red = rmj_true + 0.4*ebv*(k_lambda(k_lambda_eff(filterlist='sdss_r0.par',band_shift=0.0),/calz)-$
      k_lambda(k_lambda_eff(filterlist='twomass_J.par',band_shift=0.0),/calz))

    arrow, rmj_true, nuvmr_true, rmj_red, nuvmr_red, /data, $
      hsize=-0.25, hthick=8, thick=8
    xyouts, rmj_true+0.05, nuvmr_true-0.15, 'E(B-V) = '+string(ebv,format='(F3.1)'), $
      charsize=1.3, align=0.0

; ages
;   djs_plot, [0], [0], /nodata, position=pos, $
;     xrange=xrange, yrange=yrange, xsty=1, ysty=1, $
;     xtitle=xtitle, ytitle=ytitle
;    these = where((parent.z ge zbins[0].zlo) and $
;      (parent.z lt zbins[0].zup),ngal)
;    nuvmr = parent[these].k_galex_absmag_01[1]-parent[these].k_ugriz_absmag_01[2]
;    rmj = parent[these].k_ugriz_absmag_01[2]-parent[these].k_ubvrijhk_absmag_01[5]
;    djs_oplot, rmj, nuvmr, psym=symcat(6,thick=2), symsize=0.3, color='orange'
    
    im_plotconfig, /psclose

; #########################
; U-V vs V-J plot
    psfile = paperpath+'umv_vmj_models'+suffix
    im_plotconfig, 0, pos, psfile=psfile

    xrange = [0.0,2.5]
    yrange = [-0.2,2.9]
    xtitle = textoidl('V - J')
    ytitle = textoidl('U - V')
    
; sdss; require NUV and at least *one* of JHKs observed photometry 
    umv = sdss.k_ubvrijhk_absmag_00[0]-sdss.k_ubvrijhk_absmag_00[2]
    vmj = sdss.k_ubvrijhk_absmag_00[2]-sdss.k_ubvrijhk_absmag_00[5]
    gg = where((total(sdss.k_maggies[7:9] gt 0.0,1) ge 1),ngg)
    vmj = vmj[gg] & umv = umv[gg]
    hogg_scatterplot, vmj, umv, position=pos, xsty=1, ysty=1, $
      yrange=yrange, xrange=xrange, xtitle=xtitle, ytitle=ytitle, $
      /outliers, outsymsize=0.1, outcolor=djs_icolor('grey'), /internal, $
      levels=[0.5,0.75,0.9];, xnpix=30, ynpix=30

; overplot the models    
    djs_oplot, models.ssp_vmj, models.ssp_umv, thick=8.0, $
      color='red', psym=-symcat(9,thick=8), symsize=2.0
    djs_oplot, models.tau_vmj, models.tau_umv, thick=8.0, $
      color='blue', psym=-symcat(6,thick=8), symsize=2.0

;; selection of quiescent galaxies plus labels
;    oplot_select_quiescent
;    xyouts, 0.4, 5.8, 'Quiescent', align=0.5
;    xyouts, 1.4, 0.6, 'Star-Forming', align=0.5
;    
;; overlay a reddening vector
;    nuvmr_true = 1.8
;    rmj_true = 1.3
;
;    ebv = 0.3
;    nuvmr_red = nuvmr_true + 0.4*ebv*(k_lambda(k_lambda_eff(filterlist='galex_NUV.par',band_shift=0.0),/calz)-$
;      k_lambda(k_lambda_eff(filterlist='sdss_r0.par',band_shift=0.0),/calz))
;    rmj_red = rmj_true + 0.4*ebv*(k_lambda(k_lambda_eff(filterlist='sdss_r0.par',band_shift=0.0),/calz)-$
;      k_lambda(k_lambda_eff(filterlist='twomass_J.par',band_shift=0.0),/calz))
;
;    arrow, rmj_true, nuvmr_true, rmj_red, nuvmr_red, /data, $
;      hsize=-0.25, hthick=8, thick=8
;    xyouts, rmj_true+0.05, nuvmr_true-0.15, 'E(B-V) = '+string(ebv,format='(F3.1)'), $
;      charsize=1.3, align=0.0

    im_plotconfig, /psclose

stop    
    
; --------------------------------------------------
; NUV-r vs r-J in six redshift bins
    psfile = paperpath+'nuvmr_vs_rmj'+suffix
    im_plotconfig, 7, pos, psfile=psfile, height=3.0*[1,1,1]

    xrange = [0.0,1.9]
    yrange = [0.01,6.7]
    xtitle = textoidl('^{0.1}r-J')    

    mjaxis = im_array(-1.0,2.0,0.02)
    
    for ii = 0, nzbins-1 do begin
       these = where((parent.z ge zbins[ii].zlo) and $
         (parent.z lt zbins[ii].zup),ngal)
       mips = where((parent[these].phot_mips24 gt 0.0),comp=notmips)
       
; NUV-r and r-J       
       nuvmr = parent[these].galex_absmag[1]-parent[these].ugriz_absmag[2]
       rmj = parent[these].ugriz_absmag[2]-(parent[these].ubvrijhk_absmag[5]+jv2ab)

       qq = select_quiescent(nuvmr,rmj,active=aa)
       
       if odd(ii) then begin
          ytickname = replicate(' ',10)
          delvarx, ytitle
       endif else begin
          delvarx, ytickname
          ytitle = textoidl('NUV-^{0.1}r')
       endelse
       if (ii lt 4) then xtickname = replicate(' ',10) else $
         delvarx, xtickname
;      if (ii eq 2) then ytitle = 'M_{NUV}-M_{0.1r}' else delvarx, ytitle
       
      hogg_scatterplot, rmj, nuvmr, noerase=(ii gt 0), $
        position=pos[*,ii], xsty=1, ysty=1, $
        yrange=yrange, xrange=xrange, xtickname=xtickname, $
        ytickname=ytickname, xtitle='', ytitle=ytitle, $
        /outliers, outcolor=djs_icolor('grey'), /internal, $
        levels=[0.5,0.75,0.9], xnpix=15, ynpix=15
      djs_oplot, rmj[mips], nuvmr[mips], psym=symcat(6,thick=2), $
        color='red', symsize=0.3
      oplot_select_quiescent
;     djs_oplot, mjaxis, poly(mjaxis,[1.7,2.1]), line=5, thick=3

;      djs_plot, [0], [0], /nodata, noerase=(ii gt 0), $
;        position=pos[*,ii], xsty=1, ysty=1, $
;        yrange=yrange, xrange=xrange, xtickname=xtickname, $
;        ytickname=ytickname, xtitle='', ytitle=ytitle
;      djs_oplot, rmj[qq], nuvmr[qq], psym=symcat(6), symsize=0.6, color='red'
;      djs_oplot, rmj[aa], nuvmr[aa], psym=symcat(6), symsize=0.6, color='blue'
;;     djs_oplot, rmj, nuvmr, psym=symcat(6), symsize=0.6

       legend, 'z='+strtrim(string(zbins[ii].zlo,format='(F12.2)'),2)+$
         '-'+strtrim(string(zbins[ii].zup,format='(F12.2)'),2), $
         /right, /bottom, box=0, charsize=1.4, margin=0
    endfor
    xyouts, pos[0,5], pos[1,5]-0.07, xtitle, align=0.5, /norm
    
    im_plotconfig, /psclose
    
return
end
    
