pro vimos_calibrate, quicksex=quicksex
; jm09mar23nyu - calibrate the VIMOS imaging using the SDSS
; jm09aug05ucsd - new calibration with updated Pickles plots; SDSS
;   photometry now includes K-correct AB corrections  

    common vimos_sdss_cat, allsdss, allbcat, allvcat, allrcat

    sexpath = sg1120_path(/sex)
    mosaicpath = vimos_path(/mosaics)
    
; build a quick SE catalog of each mosaic for calibration; detect off
; the chi2 image so that the gr catalogs have the same number of
; sources 
    if keyword_set(quicksex) then begin
       imagelist = mosaicpath+'sg1120_'+['B','V','R']+'.fits'
       weightlist = repstr(imagelist,'.fits','.weight.fits')
       catlist = repstr(imagelist,'.fits','.calib.cat')
       detect_imagename = mosaicpath+'sg1120_BVR_chi2.fits'
       detect_weightname = mosaicpath+'sg1120_BVR_chi2.weight.fits'
       
; initialize the SE configuration parameters
       config = init_sex_config(nimage)
       config.parameters_name = sg1120_path(/sex)+'sg1120.sex.param.mosaic'
       config.catalog_type = 'ASCII_HEAD'
       config.detect_thresh = 3.0
       config.analysis_thresh = 3.0

       config.weight_type = 'MAP_WEIGHT'
       config.weight_gain = 'Y'
       config.seeing_fwhm = 0.8 ; update this!

       pixscale = 0.205         ; pixel scale [arcsec/pixel]
       photaper_arcsec = findgen(5)+1.0 ; [diameter,arcsec]
       photaper_pixel = photaper_arcsec/pixscale ; [diameter,pixel]
       photaper = strjoin(strtrim(string(photaper_pixel,format='(F10.2)'),2),',')

       config.phot_apertures = photaper
       config.phot_fluxfrac = '0.2,0.5,0.9'
       config.pixel_scale = pixscale

       t0 = systime(1)
       for ii = 0, n_elements(imagelist)-1 do begin
          config.catalog_name = catlist[ii]
          config.weight_image = detect_weightname+','+weightlist[ii]
          im_sex, imagelist[ii], config, detect_imagelist=detect_imagename
          splog, 'Total time to generate SE catalogs = ', $
            (systime(1)-t0)/60.0, ' minutes'
; write a region file
;         radecregname = repstr(catlist[ii],'.cat','.reg')
;         splog, 'Witing DS9 region file'+radecregname
;         cat = rsex(catlist[ii])
;         write_ds9_regionfile, cat.xwin_world, cat.ywin_world, $
;           filename=radecregname, color='red', symbol='box'
       endfor
       return
    endif

; read the catalogs
    if (n_elements(allsdss) eq 0L) then begin
       splog, 'Reading '+sexpath+'sg1120_sdss_dr7.fits.gz'
       allsdss = mrdfits(sexpath+'sg1120_sdss_dr7.fits.gz',1)
    endif
    if (n_elements(allbcat) eq 0L) then begin
       splog, 'Reading '+mosaicpath+'sg1120_B.calib.cat'
       allbcat = rsex(mosaicpath+'sg1120_B.calib.cat')
    endif
    if (n_elements(allvcat) eq 0L) then begin
       splog, 'Reading '+mosaicpath+'sg1120_V.calib.cat'
       allvcat = rsex(mosaicpath+'sg1120_V.calib.cat')
    endif
    if (n_elements(allrcat) eq 0L) then begin
       splog, 'Reading '+mosaicpath+'sg1120_R.calib.cat'
       allrcat = rsex(mosaicpath+'sg1120_R.calib.cat')
    endif

; identify stars in the SDSS catalog; this needs to match
; LDSS3_CALIBRATE!
    ww = where($
      (allsdss.g ge 14.5) and (allsdss.g le 22.0) and $
      (allsdss.r ge 14.5) and (allsdss.r le 21.0) and $
      (allsdss.radius_g le 2.0) and (allsdss.radius_g ge 1.1) and $
      (allsdss.radius_r le 2.0) and (allsdss.radius_r ge 1.1) and $
      (allsdss.type eq 6))
    sdss1 = allsdss[ww]

; clean up the VIMOS photometry    
    good = where((allbcat.flags eq 0) and (allvcat.flags eq 0) and $
      (allrcat.flags eq 0),ngood)
    bcat1 = allbcat[good]
    vcat1 = allvcat[good]
    rcat1 = allrcat[good]
    
; spherematch the SDSS stars against the R-band catalog and write a
; region file
    spherematch, rcat1.xwin_world, rcat1.ywin_world, $
      sdss1.ra, sdss1.dec, 1.0/3600.0, m1, m2
    bcat = bcat1[m1]
    vcat = vcat1[m1]
    rcat = rcat1[m1]
    sdss = sdss1[m2]

    write_ds9_regionfile, rcat1.xwin_world, rcat1.ywin_world, filename=mosaicpath+$
      'vimos_calibrate.all.reg', color='red', symbol='circle'
    write_ds9_regionfile, rcat.xwin_world, rcat.ywin_world, filename=mosaicpath+$
      'vimos_calibrate.stars.reg', color='blue', symbol='circle'

; correct for Galactic extinction (not negligible!); for the SDSS
; magnitudes use the extinction values obtained from CAS (see
; BUILD_SG1120_SDSS_REFCAT) 
    euler, rcat.xwin_world, rcat.ywin_world, gl, gb, 1
    ebv = dust_getval(gl,gb,/interp)

    bdust = ebv*k_lambda(k_lambda_eff(filterlist='vimos_B.par'),/odonnell)
    vdust = ebv*k_lambda(k_lambda_eff(filterlist='vimos_V.par'),/odonnell)
    rdust = ebv*k_lambda(k_lambda_eff(filterlist='vimos_R.par'),/odonnell)

; convert the SDSS (AB) magnitudes to BVR (Vega) using the color
; transformations from Lupton+05; see
; http://www.sdss.org/dr7/algorithms/sdssUBVRITransform.html
    
;    B = u - 0.8116*(u - g) + 0.1313;  sigma = 0.0095
;    B = g + 0.3130*(g - r) + 0.2271;  sigma = 0.0107
;  
;    V = g - 0.2906*(u - g) + 0.0885;  sigma = 0.0129
;    V = g - 0.5784*(g - r) - 0.0038;  sigma = 0.0054
;  
;    R = r - 0.1837*(g - r) - 0.0971;  sigma = 0.0106
;    R = r - 0.2936*(r - i) - 0.1439;  sigma = 0.0072
;  
;    I = r - 1.2444*(r - i) - 0.3820;  sigma = 0.0078
;    I = i - 0.3780*(i - z)  -0.3974;  sigma = 0.006

    uu = sdss.u - sdss.udust
    gg = sdss.g - sdss.gdust
    rr = sdss.r - sdss.rdust
    ii = sdss.i - sdss.idust

    uu_err = sdss.uerr
    gg_err = sdss.gerr
    rr_err = sdss.rerr
    ii_err = sdss.ierr

    ug = uu-gg
    gr = gg-rr
    ri = rr-ii
    
    ug_err = sqrt(uu_err^2+gg_err^2)
    gr_err = sqrt(gg_err^2+rr_err^2)
    ri_err = sqrt(rr_err^2+ii_err^2)

    filterlist = 'vimos_'+['B','V','R','I']+'.par'
;   filterlist = 'bessell_'+['B','V','R','I']+'.par'
    v2ab = k_vega2ab(filterlist=filterlist,/kurucz,/silent)

;   sdss_bmag = (uu - 0.8116*ug + 0.1313) + v2ab[0];  sigma = 0.0095
    sdss_bmag = (gg + 0.3130*gr + 0.2271) + v2ab[0] ;  sigma = 0.0107
    sdss_vmag = (gg - 0.5784*gr - 0.0038) + v2ab[1];-0.05 ;  sigma = 0.0054
;   sdss_vmag = (gg - 0.2906*ug + 0.0885) + v2ab[1] ;  sigma = 0.0129
    sdss_rmag = (rr - 0.2936*ri - 0.1439) + v2ab[2] ;  sigma = 0.0072
    sdss_imag = (rr - 1.2444*ri - 0.3820) + v2ab[3] ;  sigma = 0.0078

    sdss_bmag_err = sqrt(gg_err^2 + (0.3130*gr_err)^2)
    sdss_vmag_err = sqrt(gg_err^2 + (0.5784*gr_err)^2)
    sdss_rmag_err = sqrt(rr_err^2 + (0.2936*ri_err)^2)
    sdss_imag_err = sqrt(rr_err^2 + (1.2444*ri_err)^2)

; VIMOS photometry
    instr_bmag = bcat.mag_auto - bdust; + v2ab[0] ; AB
    instr_vmag = vcat.mag_auto - vdust; + v2ab[1] ; AB
    instr_rmag = rcat.mag_auto - rdust; + v2ab[2] ; AB

; compute the half-light radii    
    pixscale = 0.205
    bhalf1 = pixscale*bcat1.flux_radius1*bcat1.awin_image
    vhalf1 = pixscale*vcat1.flux_radius1*vcat1.awin_image
    rhalf1 = pixscale*rcat1.flux_radius1*rcat1.awin_image

    bhalf = pixscale*bcat.flux_radius1*bcat.awin_image
    vhalf = pixscale*vcat.flux_radius1*vcat.awin_image
    rhalf = pixscale*rcat.flux_radius1*rcat.awin_image
    
    fit = where($
      (instr_bmag ge 18.5) and (instr_bmag le 23.0) and $
      (instr_vmag ge 18.5) and (instr_vmag le 22.0) and $
      (instr_rmag ge 18.5) and (instr_rmag le 21.0) and $
      (bhalf le 0.8) and (vhalf le 0.6) and (rhalf le 0.6),nfit)

; region file of the objects used in the fitting    
    write_ds9_regionfile, rcat[fit].xwin_world, rcat[fit].ywin_world, $
      filename=mosaicpath+'vimos_calibrate.fitstars.reg', $
      color='red', symbol='circle'
    
; make a quick plot
    psfile = mosaicpath+'qaplot_vimos_calibrate.ps'
    im_plotconfig, 0, pos1, psfile=psfile, xmargin=[1.2,0.2], charsize=1.8
; ### SDSS stars    
; r-band
    djs_plot, allsdss.r, allsdss.radius_r, psym=6, symsize=0.2, yrange=[0,8], $
      xrange=[13,27], xsty=1, ysty=1, xtitle='r_{SDSS} (AB mag)', $
      ytitle='Petrosian radius (arcsec)', thick=3.0, $
      position=pos1
;   legend, 'SDSS/DR7', /left, /top, box=0, charsize=2.0
    plotsym, 0, 0.8, /fill
    djs_oplot, sdss.r, sdss.radius_r, ps=8, color='red'
    plotsym, 8, 1.1, fill=0, thick=5.0
    djs_oplot, sdss[fit].r, sdss[fit].radius_r, ps=8, color='blue'
; g-band
    djs_plot, allsdss.g, allsdss.radius_g, psym=6, symsize=0.2, yrange=[0,8], $
      xrange=[13,27], xsty=1, ysty=1, xtitle='g_{SDSS} (AB mag)', $
      ytitle='Petrosian radius (arcsec)', thick=3.0, $
      position=pos1
;   legend, 'SDSS/DR7', /left, /top, box=0, charsize=2.0
    plotsym, 0, 0.8, /fill
    djs_oplot, sdss1.g, sdss1.radius_g, ps=8, color='red'
    plotsym, 8, 1.1, fill=0, thick=5.0
    djs_oplot, sdss[fit].g, sdss[fit].radius_g, ps=8, color='blue'
; ### B-band
    djs_plot, [0], [0], /nodata, /xsty, /ysty, xthick=4.0, ythick=4.0, $
      charthick=4.0, xrange=[14.0,29.0], yrange=[0,8.0], $
      xtitle='B_{VIMOS} (instrumental mag)', ytitle='B-band Half-light Radius (arcsec)', $
      position=pos1
    djs_oplot, bcat1.mag_auto, bhalf1, ps=6, sym=0.4, thick=5.0
; objects identified as stars in the SDSS with non-zero SE flags
    plotsym, 0, 0.8, /fill
    djs_oplot, instr_bmag, bhalf, psym=8, color='red'
    plotsym, 8, 1.1, fill=0, thick=5.0
    djs_oplot, instr_bmag[fit], bhalf[fit], psym=8, color='blue'
; ### V-band
    djs_plot, [0], [0], /nodata, /xsty, /ysty, xthick=4.0, ythick=4.0, $
      charthick=4.0, xrange=[14.0,29.0], yrange=[0,8.0], $
      xtitle='V_{VIMOS} (instrumental mag)', ytitle='V-band Half-light Radius (arcsec)', $
      position=pos1
    djs_oplot, vcat1.mag_auto, vhalf1, ps=6, sym=0.4, thick=5.0
; objects identified as stars in the SDSS with non-zero SE flags
    plotsym, 0, 0.8, /fill
    djs_oplot, instr_vmag, vhalf, psym=8, color='red'
    plotsym, 8, 1.1, fill=0, thick=5.0
    djs_oplot, instr_vmag[fit], vhalf[fit], psym=8, color='blue'
; ### R-band
    djs_plot, [0], [0], /nodata, /xsty, /ysty, xthick=4.0, ythick=4.0, $
      charthick=4.0, xrange=[14.0,29.0], yrange=[0,8.0], $
      xtitle='R_{VIMOS} (instrumental mag)', ytitle='R-band Half-light Radius (arcsec)', $
      position=pos1
    djs_oplot, rcat1.mag_auto, rhalf1, ps=6, sym=0.4, thick=5.0
; stars
    plotsym, 0, 0.8, /fill
    djs_oplot, instr_rmag, rhalf, psym=8, color='red'
    plotsym, 8, 1.1, fill=0, thick=5.0
    djs_oplot, instr_rmag[fit], rhalf[fit], psym=8, color='blue'
;   im_plotconfig, /psclose

; make a plot showing the rough transformation
    im_plotconfig, 6, pos2, xmargin=[1.5,0.2], height=[4.5,2.5], $
      yspace=1.1, charsize=1.8
    im_plotfaves, /postscript
    magaxis = im_array(10.0,30.0,0.1)
; B-band
    plotsym, 0, 1.0, /fill
    djs_plot, instr_bmag, sdss_bmag, psym=8, xsty=1, ysty=1, $
      xtitle='B_{VIMOS} (instrumental mag)', $
      ytitle='B_{SDSS} (AB mag)', xrange=[15,24.0], $
      yrange=[15,24.0], color='blue', position=pos2[*,0]
    plotsym, 8, 1.2, fill=0, thick=5.0
    djs_oplot, instr_bmag[fit], sdss_bmag[fit], psym=8, color='red'
    bcoeff = im_linefit(instr_bmag[fit],sdss_bmag[fit],coeff_fixed=[0,1],$
      coeff_guess=[-1.0,1.0],yerr=sdss_bmag_err[fit],chi2=bchi2,yfit=bfit,$
      coeff_err=bcoeff_err)
    bresid = djsig(sdss_bmag[fit]-bfit)
    djs_oplot, magaxis, poly(magaxis,bcoeff), line=0, thick=4
    im_legend, ['B_{AB} = '+strtrim(string(bcoeff[0],format='(F12.3)'),2)+$
      '+B_{VIMOS}*'+strtrim(string(bcoeff[1],format='(F12.2)'),2),$
      '\sigma_{B} = '+string(bresid,format='(F5.3)')], /right, /bottom, box=0
    im_legend, ['N = '+string(nfit,format='(I0)'),'Extinction Corrected'], $
      /left, /top, box=0
; residuals
    plotsym, 0, 1.0, /fill
    djs_plot, instr_bmag-instr_vmag, sdss_bmag-instr_bmag, /noerase, psym=8, xsty=1, ysty=1, $
      ytitle='B_{SDSS}-B_{VIMOS}', $
      xtitle='B_{VIMOS} - V_{VIMOS} (instrumental mag)', xrange=[-0.1,1.9], $
      yrange=[-0.8,0.2], color='blue', position=pos2[*,1]
    plotsym, 8, 1.2, fill=0, thick=5.0
    djs_oplot, instr_bmag[fit]-instr_vmag[fit], sdss_bmag[fit]-instr_bmag[fit], psym=8, color='red'
    djs_oplot, !x.crange, bcoeff[0]*[1,1], line=2, thick=5.0
;   im_legend, ['\sigma_{B} = '+string(bresid,format='(F5.3)')], $
;     /left, /top, box=0
; V-band
    plotsym, 0, 1.0, /fill
    djs_plot, instr_vmag, sdss_vmag, psym=8, xsty=1, ysty=1, $
      xtitle='V_{VIMOS} (instrumental mag)', $
      ytitle='V_{SDSS} (AB mag)', xrange=[15,22.5], $
      yrange=[15,22.5], color='blue', position=pos2[*,0]
    plotsym, 8, 1.2, fill=0, thick=5.0
    djs_oplot, instr_vmag[fit], sdss_vmag[fit], psym=8, color='red'
    vcoeff = im_linefit(instr_vmag[fit],sdss_vmag[fit],coeff_fixed=[0,1],$
      coeff_guess=[-1.0,1.0],yerr=sdss_vmag_err[fit],chi2=vchi2,yfit=vfit,$
      coeff_err=vcoeff_err)
    vresid = djsig(sdss_vmag[fit]-vfit)
    djs_oplot, magaxis, poly(magaxis,vcoeff), line=0, thick=4
    im_legend, ['V_{AB} = '+strtrim(string(vcoeff[0],format='(F12.3)'),2)+$
      '+V_{VIMOS}*'+strtrim(string(vcoeff[1],format='(F12.2)'),2),$
      '\sigma_{V} = '+string(vresid,format='(F5.3)')], /right, /bottom, box=0
    im_legend, ['N = '+string(nfit,format='(I0)'),'Extinction Corrected'], /left, /top, box=0
; residuals
    plotsym, 0, 1.0, /fill
    djs_plot, instr_vmag-instr_rmag, sdss_vmag-instr_vmag, /noerase, psym=8, xsty=1, ysty=1, $
      ytitle='V_{SDSS}-V_{VIMOS}', $
      xtitle='V_{VIMOS} - R_{VIMOS} (instrumental mag)', xrange=[-0.2,1.2], $
      yrange=[-0.7,0.4], color='blue', position=pos2[*,1]
    plotsym, 8, 1.2, fill=0, thick=5.0
    djs_oplot, instr_vmag[fit]-instr_rmag[fit], sdss_vmag[fit]-instr_vmag[fit], psym=8, color='red'
    djs_oplot, !x.crange, vcoeff[0]*[1,1], line=2, thick=5.0
;   im_legend, ['\sigma_{V} = '+string(vresid,format='(F5.3)')], $
;     /left, /top, box=0
; R-band
    plotsym, 0, 1.0, /fill
    djs_plot, instr_rmag, sdss_rmag, psym=8, xsty=1, ysty=1, $
      xtitle='R_{VIMOS} (instrumental mag)', $
      ytitle='R_{SDSS} (AB mag)', xrange=[15,22.5], $
      yrange=[15,22.5], color='blue', position=pos2[*,0]
    plotsym, 8, 1.2, fill=0, thick=5.0
    djs_oplot, instr_rmag[fit], sdss_rmag[fit], psym=8, color='red'
    rcoeff = im_linefit(instr_rmag[fit],sdss_rmag[fit],coeff_fixed=[0,1],$
      coeff_guess=[-1.0,1.0],yerr=sdss_rmag_err[fit],chi2=rchi2,yfit=rfit,$
      coeff_err=rcoeff_err)
    rresid = djsig(sdss_rmag[fit]-rfit)
    djs_oplot, magaxis, poly(magaxis,rcoeff), line=0, thick=4
    im_legend, ['R_{AB} = '+strtrim(string(rcoeff[0],format='(F12.3)'),2)+$
      '+R_{VIMOS}*'+strtrim(string(rcoeff[1],format='(F12.2)'),2),$
      '\sigma_{R} = '+string(rresid,format='(F5.3)')], /right, /bottom, box=0
    im_legend, ['N = '+string(nfit,format='(I0)'),'Extinction Corrected'], /left, /top, box=0
; residuals
    plotsym, 0, 1.0, /fill
    djs_plot, instr_vmag-instr_rmag, sdss_rmag-instr_rmag, /noerase, psym=8, xsty=1, ysty=1, $
      ytitle='R_{SDSS}-R_{VIMOS}', $
      xtitle='V_{VIMOS} - R_{VIMOS} (instrumental mag)', xrange=[-0.2,1.2], $
      yrange=[-0.5,0.6], color='blue', position=pos2[*,1]
    plotsym, 8, 1.2, fill=0, thick=5.0
    djs_oplot, instr_vmag[fit]-instr_rmag[fit], sdss_rmag[fit]-instr_rmag[fit], psym=8, color='red'
    djs_oplot, !x.crange, rcoeff[0]*[1,1], line=2, thick=5.0
;   im_legend, ['\sigma_{R} = '+string(rresid,format='(F5.3)')], $
;     /left, /top, box=0

    im_plotconfig, /psclose
;   spawn, 'ps2pdf13 '+psfile+' '+repstr(psfile,'.ps','.pdf'), /sh
    
; final zeropoints:
    bzero = bcoeff[0]
    vzero = vcoeff[0]
    rzero = rcoeff[0]
    splog, '##################################################'
    splog, 'Final zeropoints (BVR):'
    print, bzero, vzero, rzero
    
; plot the stellar locus and compare against Pickles+98; just keep the
; dwarfs! 
    p98 = read_98pickles()
;   keep = where(p98.feh eq 0.0)
    keep = where(strmatch(p98.type,'*V*'))
    p98 = p98[keep]
    lambda = k_lambda_to_edges(p98[0].wave)
; test
    vv1 = fltarr(n_elements(p98))
    vv2 = fltarr(n_elements(p98))
    for ii = 0L, n_elements(p98)-1L do vv1[ii] = $
      -2.5*alog10(k_project_filters(lambda,p98[ii].flux,$
      filterlist='vimos_V.par',/silent))
    for ii = 0L, n_elements(p98)-1L do vv2[ii] = $
      -2.5*alog10(k_project_filters(lambda,p98[ii].flux,$
      filterlist='bessell_V.par',/silent))

; BVR
    filterlist = 'vimos_'+['B','V','R','I']+'.par'
    v2ab = k_vega2ab(filterlist=filterlist,/kurucz,/silent)
    maggies = fltarr(n_elements(filterlist),n_elements(p98))
    for ii = 0L, n_elements(p98)-1L do maggies[*,ii] = $
      k_project_filters(lambda,p98[ii].flux,filterlist=filterlist,/silent)
    bv_p98 = -2.5*alog10(maggies[0,*]/maggies[1,*]) ; AB
    br_p98 = -2.5*alog10(maggies[0,*]/maggies[2,*]) ; AB
    vr_p98 = -2.5*alog10(maggies[1,*]/maggies[2,*]) ; AB
    vi_p98 = -2.5*alog10(maggies[1,*]/maggies[3,*]) ; AB
    ri_p98 = -2.5*alog10(maggies[2,*]/maggies[3,*]) ; AB
; gri
    filterlist = 'sdss_'+['g0','r0','i0']+'.par'
    maggies = fltarr(n_elements(filterlist),n_elements(p98))
    for ii = 0L, n_elements(p98)-1L do maggies[*,ii] = $
      k_project_filters(lambda,p98[ii].flux,filterlist=filterlist,/silent)
    gr_sdss_p98 = -2.5*alog10(maggies[0,*]/maggies[1,*]) ; AB
    ri_sdss_p98 = -2.5*alog10(maggies[1,*]/maggies[2,*]) ; AB
    
    gr_sdss = gg[fit]-rr[fit]
    ri_sdss = rr[fit]-ii[fit]

    bv_sdss = sdss_bmag[fit]-sdss_vmag[fit]
    br_sdss = sdss_bmag[fit]-sdss_rmag[fit]
    vr_sdss = sdss_vmag[fit]-sdss_rmag[fit]
    vi_sdss = sdss_vmag[fit]-sdss_imag[fit]
    ri_sdss = sdss_rmag[fit]-sdss_imag[fit]

    bv_vimos = (instr_bmag[fit]+bzero)-(instr_vmag[fit]+vzero)
    vr_vimos = (instr_vmag[fit]+vzero)-(instr_rmag[fit]+rzero)
    
    im_plotconfig, 0, psfile=mosaicpath+'qaplot_vimos_pickles.ps', $
      charsize=1.8
; B-V vs V-R - SDSS & VIMOS
    djs_plot, [0], [0], /nodata, /xsty, /ysty, xrange=[-0.5,1.3], $
      yrange=[-0.6,2.3], xtitle='V - R (AB mag)', ytitle='B - V (AB mag)'
    djs_oplot, vr_vimos, bv_vimos, psym=6, color='red', sym=2.0, thick=4.0
    djs_oplot, vr_sdss, bv_sdss, psym=6, color='blue', sym=2.0, thick=4.0
;   djs_oplot, p98.vr+(v2ab[1]-v2ab[2]), p98.bv+(v2ab[0]-v2ab[1]), $
;     psym=7, thick=3.0, symsize=0.8 ; AB colors from the Pickles paper
    djs_oplot, vr_p98, bv_p98, psym=7, thick=4.0, symsize=1.5
    im_legend, ['VIMOS','SDSS','Pickles+98 (synth)'], /left, /top, $
      box=0, psym=[6,6,7], color=['red','blue','']
;   im_legend, ['VIMOS','SDSS','Pickles+98 (synth)',$
;     'Pickles+98 (published)'], /left, /top, $
;     box=0, psym=[6,6,6,7], color=['red','blue','','']
; g-r vs r-i - SDSS
    djs_plot, [0], [0], /nodata, /xsty, /ysty, xrange=[-0.5,1.3], $
      yrange=[-0.6,2.1], xtitle='r - i (AB mag)', ytitle='g - r (AB mag)'
    djs_oplot, ri_sdss, gr_sdss, psym=6, color='blue', sym=2.0, thick=4.0
    djs_oplot, ri_sdss_p98, gr_sdss_p98, psym=7, thick=4.0, symsize=1.5
    im_legend, ['SDSS','Pickles+98 (synth)'], /left, /top, $
      box=0, psym=[6,7], color=['blue','']
; B-V vs V-R - SDSS
    djs_plot, [0], [0], /nodata, /xsty, /ysty, xrange=[-0.5,1.3], $
      yrange=[-0.6,2.3], xtitle='V - R (AB mag)', ytitle='B - V (AB mag)'
    djs_oplot, vr_sdss, bv_sdss, psym=6, color='blue', sym=2.0, thick=4.0
    djs_oplot, vr_p98, bv_p98, psym=7, thick=4.0, symsize=1.5
    im_legend, ['SDSS','Pickles+98 (synth)'], /left, /top, $
      box=0, psym=[6,7], color=['blue','']
; B-R vs R-I - SDSS
    djs_plot, [0], [0], /nodata, /xsty, /ysty, xrange=[-0.6,1.7], $
      yrange=[-0.2,2.9], xtitle='R - I (AB mag)', ytitle='B - R (AB mag)'
    djs_oplot, ri_sdss, br_sdss, psym=6, color='blue', sym=2.0, thick=4.0
    djs_oplot, ri_p98, br_p98, psym=7, thick=4.0, symsize=1.5
    im_legend, ['SDSS','Pickles+98 (synth)'], /left, /top, $
      box=0, psym=[6,7], color=['blue','']
; B-V vs V-I - SDSS
    djs_plot, [0], [0], /nodata, /xsty, /ysty, xrange=[-0.7,2.0], $
      yrange=[-0.5,2.2], xtitle='V - I (AB mag)', ytitle='B - V (AB mag)'
    djs_oplot, vi_sdss, bv_sdss, psym=6, color='blue', sym=2.0, thick=4.0
    djs_oplot, vi_p98, bv_p98, psym=7, thick=4.0, symsize=1.5
    im_legend, ['SDSS','Pickles+98 (synth)'], /left, /top, $
      box=0, psym=[6,7], color=['blue','']
; V-R vs R-I - SDSS
    djs_plot, [0], [0], /nodata, /xsty, /ysty, xrange=[-0.6,1.7], $
      yrange=[-0.5,1.4], xtitle='R - I (AB mag)', ytitle='V - R (AB mag)'
    djs_oplot, ri_sdss, vr_sdss, psym=6, color='blue', sym=2.0, thick=4.0
    djs_oplot, ri_p98, vr_p98, psym=7, thick=4.0, symsize=1.5
    im_legend, ['SDSS','Pickles+98 (synth)'], /left, /top, $
      box=0, psym=[6,7], color=['blue','']

; B-V vs ;V-I, and V-R vs R-I

    im_plotconfig, /psclose
    
stop

return
end
