function zpt_mpfunc, x, p, color=color
; zeropoint function for MPFIT    
    return, x + p[0] + p[1]*color
end

function sg1120_zpt_mpfit, sdssmag, mag, magerr, color=color
; fit for the zeropoint
    parinfo1 = replicate({value: 1.0, fixed: 0},2)
    coeff = mpfitfun('zpt_mpfunc',sdssmag,mag,magerr,$
      parinfo=parinfo1,functargs={color: color},$
      yfit=bfit,/quiet)
return, coeff
end
    
pro vimos_sdss_zeropoints
; jm09aug03ucsd - 

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

; correct for Galactic extinction (not negligible!)    
    euler, rcat.xwin_world, rcat.ywin_world, gl, gb, 1
    ebv = dust_getval(gl,gb,/interp)

    sdss_udust = ebv*k_lambda(k_lambda_eff(filterlist='sdss_u0.par'),/odonnell)
    sdss_gdust = ebv*k_lambda(k_lambda_eff(filterlist='sdss_g0.par'),/odonnell)
    sdss_rdust = ebv*k_lambda(k_lambda_eff(filterlist='sdss_r0.par'),/odonnell)
    sdss_idust = ebv*k_lambda(k_lambda_eff(filterlist='sdss_i0.par'),/odonnell)

    bdust = ebv*k_lambda(k_lambda_eff(filterlist='vimos_B.par'),/odonnell)
    vdust = ebv*k_lambda(k_lambda_eff(filterlist='vimos_V.par'),/odonnell)
    rdust = ebv*k_lambda(k_lambda_eff(filterlist='vimos_R.par'),/odonnell)

;   udust = 0.0
;   gdust = 0.0
;   rdust = 0.0
;   idust = 0.0
;   bdust = 0.0
;   vdust = 0.0
;   bigrdust = 0.0
    
; convert the SDSS gr magnitudes to BVR using the color
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

    uu = sdss.u-sdss_udust
    gg = sdss.g-sdss_gdust
    rr = sdss.r-sdss_rdust
    ii = sdss.i-sdss_idust

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

;;   sdss_bmag = (uu - 0.8116*ug + 0.1313) + v2ab[0];  sigma = 0.0095
;    sdss_bmag = (gg + 0.3130*gr + 0.2271) + v2ab[0] ;  sigma = 0.0107
;    sdss_vmag = (gg - 0.5784*gr - 0.0038) + v2ab[1];-0.05 ;  sigma = 0.0054
;;   sdss_vmag = (gg - 0.2906*ug + 0.0885) + v2ab[1] ;  sigma = 0.0129
;    sdss_rmag = (rr - 0.2936*ri - 0.1439) + v2ab[2] ;  sigma = 0.0072
;    sdss_imag = (rr - 1.2444*ri - 0.3820) + v2ab[3] ;  sigma = 0.0078
;
;    sdss_bmag_err = sqrt(gg_err^2 + (0.3130*gr_err)^2)
;    sdss_vmag_err = sqrt(gg_err^2 + (0.5784*gr_err)^2)
;    sdss_rmag_err = sqrt(rr_err^2 + (0.2936*ri_err)^2)
;    sdss_imag_err = sqrt(rr_err^2 + (1.2444*ri_err)^2)

; VIMOS photometry
    instr_bmag = bcat.mag_auto - bdust; + v2ab[0] ; AB
    instr_vmag = vcat.mag_auto - vdust; + v2ab[1] ; AB
    instr_rmag = rcat.mag_auto - rdust; + v2ab[2] ; AB

    minerror = 0.05
    instr_bmag_err = sqrt(bcat.magerr_auto^2.0 + minerror^2.0)
    instr_vmag_err = sqrt(vcat.magerr_auto^2.0 + minerror^2.0)
    instr_rmag_err = sqrt(rcat.magerr_auto^2.0 + minerror^2.0)

; compute the half-light radii    
    pixscale = 0.205
    bhalf1 = pixscale*bcat1.flux_radius1*bcat1.awin_image
    vhalf1 = pixscale*vcat1.flux_radius1*vcat1.awin_image
    rhalf1 = pixscale*rcat1.flux_radius1*rcat1.awin_image

    bhalf = pixscale*bcat.flux_radius1*bcat.awin_image
    vhalf = pixscale*vcat.flux_radius1*vcat.awin_image
    rhalf = pixscale*rcat.flux_radius1*rcat.awin_image
    
    fit = where($
      (gg ge 14.5) and (gg le 21.0) and $
      (rr ge 14.5) and (rr le 21.0) and $
      (instr_bmag ge 18.5) and (instr_bmag le 23.0) and $
      (instr_vmag ge 18.5) and (instr_vmag le 22.0) and $
      (instr_rmag ge 18.5) and (instr_rmag le 21.0) and $
      (bhalf le 0.8) and (vhalf le 0.6) and (rhalf le 0.6),nfit)
    
    bfit = primus_zpt_mpfit(gg[fit],instr_bmag[fit],$
      instr_bmag_err[fit],color=gr[fit])
    vfit = primus_zpt_mpfit(gg[fit],instr_vmag[fit],$
      instr_vmag_err[fit],color=gr[fit])
    rfit = primus_zpt_mpfit(rr[fit],instr_vmag[fit],$
      instr_rmag_err[fit],color=ri[fit])
    splog, 'B ', bfit
    splog, 'V ', vfit
    splog, 'R ', rfit

    djs_plot, gg, instr_bmag, psym=6, xsty=3, ysty=3, xr=[14,22], yr=[15,23]
    djs_oplot, gg[fit], instr_bmag[fit], psym=6, color='red'
    

stop    
    
    allcat = mrdfits(analysis_path+'sg1120_parent_refR2.0chi2_v2.0.fits.gz',1)
    psfile = mosaicpath+'vimos_sdss_zeropoints.ps'

; select stars
    these = where($
      (allcat.phot_b gt 0.0) and (allcat.phot_b lt 20.0) and $
      (allcat.phot_v gt 0.0) and $
      (allcat.phot_r gt 0.0) and $
      (allcat.phot_sdssg ge 14.5) and (allcat.phot_sdssg le 22.0) and $
      (allcat.phot_sdssr ge 14.5) and (allcat.phot_sdssr le 21.0) and $
;     (allcat.radius_g le 2.0) and (allsdss.radius_g ge 1.1) and $
;     (allcat.radius_r le 2.0) and (allsdss.radius_r ge 1.1) and $
      (allcat.sdss_type eq 6))
    cat = allcat[these]

;   cat.phot_b = cat.phot_b - 0.298
;   cat.phot_v = cat.phot_v - 0.125
;   cat.phot_r = cat.phot_r + 0.085
    
    minerror = 0.05
    cat.phot_b_err = sqrt(cat.phot_b_err^2.0+minerror^2.0)
    cat.phot_v_err = sqrt(cat.phot_v_err^2.0+minerror^2.0)
    cat.phot_r_err = sqrt(cat.phot_r_err^2.0+minerror^2.0)
    
    bfit = primus_zpt_mpfit(cat.phot_sdssg,cat.phot_b,$
      cat.phot_b_err,color=cat.phot_sdssg-cat.phot_sdssr)
    vfit = primus_zpt_mpfit(cat.phot_sdssg,cat.phot_v,$
      cat.phot_v_err,color=cat.phot_sdssg-cat.phot_sdssr)
    rfit = primus_zpt_mpfit(cat.phot_sdssr,cat.phot_r,$
      cat.phot_r_err,color=cat.phot_sdssr-cat.phot_sdssi)
    splog, 'B ', bfit
    splog, 'V ', vfit
    splog, 'R ', rfit

    plot, cat.phot_sdssg, cat.phot_b, psym=6, xsty=3, ysty=3, $
      xr=[14,20.5], yr=[15,21]
    
stop    
    
; spherematch, convert to maggies, and identify stars    
    spherematch, sg1120.ra, sg1120.dec, sdss1.ra, $
      sdss1.dec, 1.0/3600.0, m1, m2

    sdss = zpt_sdsscat2mag(sdss1[m2])
    deep = zpt_sg1120cat2mag(sg1120[m1])

    these = where((deep.b gt 0.0) and (deep.r gt 0.0) and $
      (deep.i gt 0.0) and (deep.badflag eq 0) and $
      (sdss.r ge 19.0) and (sdss.r le 20.5))
    sdss = sdss[these]
    deep = deep[these]

; fit for the color terms and zeropoints
    bfit = primus_zpt_mpfit(sdss.g,deep.b,$
      deep.berr,color=sdss.g-sdss.r)
    rfit = primus_zpt_mpfit(sdss.r,deep.r,$
      deep.rerr,color=sdss.g-sdss.r)
    ifit = primus_zpt_mpfit(sdss.i,deep.i,$
      deep.ierr,color=sdss.r-sdss.i)

; finally make a QAplot
    zpt_scatterplot, deep.b, sdss.g, sdss.g-sdss.r, $
      deep.b-(sdss.g+bfit[0]), xtitle1='B (SG1120)', $
      ytitle1='g (SDSS)', xtitle2='g - r', $
      ytitle2='Residuals'

stop    
    
return
end
    
