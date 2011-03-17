pro ldss3_calibrate, quicksex=quicksex
; jm09mar23nyu - calibrate the LDSS3 imaging using the SDSS
; jm09aug05ucsd - slightly revised to match changes made to
;   VIMOS_CALIBRATE; SDSS photometry now includes K-correct AB
;   corrections 

    common ldss3_sdss_cat, allsdss, allgcat, allrcat

    sexpath = sg1120_path(/sex)
    mosaicpath = ldss3_path(/mosaics)
    
; build a quick SE catalog of each mosaic for calibration; detect off
; the chi2 image so that the gr catalogs have the same number of
; sources 
    if keyword_set(quicksex) then begin
       imagelist = mosaicpath+'sg1120_'+['gprime','rprime']+'.fits'
       weightlist = repstr(imagelist,'.fits','.weight.fits')
       catlist = repstr(imagelist,'.fits','.calib.cat')
       detect_imagename = mosaicpath+'sg1120_gr_chi2.fits'
       detect_weightname = mosaicpath+'sg1120_gr_chi2.weight.fits'
       
; initialize the SE configuration parameters
       config = init_sex_config(nimage)
       config.parameters_name = sg1120_path(/sex)+'sg1120.sex.param.mosaic'
       config.catalog_type = 'ASCII_HEAD'
       config.detect_thresh = 3.0
       config.analysis_thresh = 3.0

       config.weight_type = 'MAP_WEIGHT'
       config.weight_gain = 'Y'
       config.seeing_fwhm = 1.0 ; update this!

       pixscale = 0.188         ; pixel scale [arcsec/pixel]
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
    if (n_elements(allgcat) eq 0L) then begin
       splog, 'Reading '+mosaicpath+'sg1120_gprime.calib.cat'
       allgcat = rsex(mosaicpath+'sg1120_gprime.calib.cat')
    endif
    if (n_elements(allrcat) eq 0L) then begin
       splog, 'Reading '+mosaicpath+'sg1120_rprime.calib.cat'
       allrcat = rsex(mosaicpath+'sg1120_rprime.calib.cat')
    endif

; identify stars in the SDSS catalog
    ww = where($
      (allsdss.g ge 14.5) and (allsdss.g le 22.0) and $
      (allsdss.r ge 14.5) and (allsdss.r le 21.0) and $
      (allsdss.radius_g le 2.0) and (allsdss.radius_g ge 1.1) and $
      (allsdss.radius_r le 2.0) and (allsdss.radius_r ge 1.1) and $
      (allsdss.type eq 6))
    sdss1 = allsdss[ww]

; clean up the LDSS3 photometry    
    good = where((allrcat.flags eq 0) and (allgcat.flags eq 0),ngood)
    gcat1 = allgcat[good]
    rcat1 = allrcat[good]

; correct everything for Galactic extinction; for the SDSS magnitudes
; use the extinction values obtained from CAS (see
; BUILD_SG1120_SDSS_REFCAT)
    euler, rcat1.xwin_world, rcat1.ywin_world, gl, gb, 1
    ebv = dust_getval(gl,gb,/interp)

    gdust = ebv*k_lambda(k_lambda_eff(filterlist='sdss_g0.par'),/odonnell)
    rdust = ebv*k_lambda(k_lambda_eff(filterlist='sdss_r0.par'),/odonnell)
    gcat1.mag_auto = gcat1.mag_auto - gdust
    rcat1.mag_auto = rcat1.mag_auto - rdust

    sdss1.g = sdss1.g - sdss1.gdust
    sdss1.r = sdss1.r - sdss1.rdust
    
; spherematch the SDSS stars against the r-band catalog and write a
; region file
    spherematch, rcat1.xwin_world, rcat1.ywin_world, $
      sdss1.ra, sdss1.dec, 1.0/3600.0, m1, m2
    gcat = gcat1[m1]
    rcat = rcat1[m1]
    sdss = sdss1[m2]

    write_ds9_regionfile, rcat1.xwin_world, rcat1.ywin_world, filename=mosaicpath+$
      'ldss3_calibrate.all.reg', color='red', symbol='circle'
    write_ds9_regionfile, rcat.xwin_world, rcat.ywin_world, filename=mosaicpath+$
      'ldss3_calibrate.stars.reg', color='blue', symbol='circle'

; now do it!    
    sdss_gmag = sdss.g
    sdss_rmag = sdss.r
    sdss_gmag_err = sdss.gerr
    sdss_rmag_err = sdss.rerr

    instr_gmag = gcat.mag_auto
    instr_rmag = rcat.mag_auto
    instr_gmag_err = gcat.magerr_auto
    instr_rmag_err = rcat.magerr_auto

    pixscale = 0.188
    ghalf1 = pixscale*gcat1.flux_radius1*gcat1.awin_image
    rhalf1 = pixscale*rcat1.flux_radius1*rcat1.awin_image
    ghalf = pixscale*gcat.flux_radius1*gcat.awin_image
    rhalf = pixscale*rcat.flux_radius1*rcat.awin_image
    
    fit = where($
      (instr_gmag ge 18.0) and (instr_gmag le 23.0) and $
      (instr_rmag ge 18.7) and (instr_rmag le 21.0) and $
      (ghalf le 2.0) and (rhalf le 2.0),nfit)

; region file of the objects used in the fitting    
    write_ds9_regionfile, rcat[fit].xwin_world, rcat[fit].ywin_world, $
      filename=mosaicpath+'ldss3_calibrate.fitstars.reg', $
      color='red', symbol='circle'
    
; make a quick plot
    psfile = mosaicpath+'qaplot_ldss3_calibrate.ps'
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
; ### LDSS3 ###
; g-band
    djs_plot, [0], [0], /nodata, /xsty, /ysty, xthick=4.0, ythick=4.0, $
      charthick=4.0, charsize=2.0, xrange=[14.0,29.0], yrange=[0,8.0], $
      xtitle='g_{LDSS3} (instrumental mag)', ytitle='g-band Half-light Radius (arcsec)', $
      position=pos1
    djs_oplot, gcat1.mag_auto, ghalf1, ps=6, sym=0.4, thick=5.0
    plotsym, 0, 0.8, /fill
    djs_oplot, instr_gmag, ghalf, psym=8, color='red'
    plotsym, 8, 1.1, fill=0, thick=5.0
    djs_oplot, instr_gmag[fit], ghalf[fit], psym=8, color='blue'
; r-band
    djs_plot, [0], [0], /nodata, /xsty, /ysty, xthick=4.0, ythick=4.0, $
      charthick=4.0, charsize=2.0, xrange=[14.0,29.0], yrange=[0,8.0], $
      xtitle='r_{LDSS3} (instrumental mag)', ytitle='r-band Half-light Radius (arcsec)', $
      position=pos1
    djs_oplot, rcat1.mag_auto, rhalf1, ps=6, sym=0.4, thick=5.0
    plotsym, 0, 0.8, /fill
    djs_oplot, instr_rmag, rhalf, psym=8, color='red'
    plotsym, 8, 1.1, fill=0, thick=5.0
    djs_oplot, instr_rmag[fit], rhalf[fit], psym=8, color='blue'
;   im_plotconfig, /psclose

; solve for the zero-point and a color term
;   functargs = {color: gmag[fit]-rmag[fit]}
;   parinfo = replicate({value: 0.0, fixed: 0},2)
;   parinfo[0].value = -10.0 & parinfo[0].fixed = 0 ; zeropoint
;   parinfo[1].value = 0.01  & parinfo[1].fixed = 0 ; color term
;   param = mpfitfun('zptfunc',rmag[fit],instr_rmag[fit],instr_rmag_err[fit],$
;     parinfo=parinfo,functargs=functargs,perror=perror,quiet=0,$
;     dof=dof,bestnorm=bestnorm,yfit=best)
;   niceprint, best, instr_rmag[fit]

; make a plot showing the rough transformation
    im_plotconfig, 6, pos2, xmargin=[1.5,0.2], height=[4.5,2.5], $
      yspace=1.1, charsize=1.8
    im_plotfaves, /postscript
    magaxis = im_array(10.0,30.0,0.1)
; g-band
    plotsym, 0, 1.0, /fill
    djs_plot, instr_gmag, sdss_gmag, psym=8, xsty=1, ysty=1, $
      xtitle='g_{LDSS3} (instrumental mag)', $
      ytitle='g_{SDSS} (AB mag)', xrange=[14,23.0], $
      yrange=[14,23.0], color='blue', position=pos2[*,0]
    plotsym, 8, 1.2, fill=0, thick=5.0
    djs_oplot, instr_gmag[fit], sdss_gmag[fit], psym=8, color='red'
    gcoeff = im_linefit(instr_gmag[fit],sdss_gmag[fit],coeff_fixed=[0,1],$
      coeff_guess=[-1.0,1.0],yerr=sdss_gmag_err[fit],chi2=gchi2,yfit=gfit,$
      coeff_err=gcoeff_err)
    gresid = djsig(sdss_gmag[fit]-gfit)
    djs_oplot, magaxis, poly(magaxis,gcoeff), line=0, thick=4
    im_legend, ['g_{AB} = '+strtrim(string(gcoeff[0],format='(F12.3)'),2)+$
      '+g_{LDSS3}*'+strtrim(string(gcoeff[1],format='(F12.2)'),2),$
      '\sigma_{g} = '+string(gresid,format='(F5.3)')], /right, /bottom, box=0
    im_legend, ['N = '+string(nfit,format='(I0)')], $
      /left, /top, box=0
; residuals
    plotsym, 0, 1.0, /fill
    djs_plot, instr_gmag-instr_rmag, sdss_gmag-instr_gmag, /noerase, psym=8, xsty=1, ysty=1, $
      ytitle='g_{SDSS}-g_{LDSS3}', $
      xtitle='g_{LDSS3} - r_{LDSS3} (instrumental mag)', xrange=[-0.2,1.7], $
      yrange=[-0.2,0.8], color='blue', position=pos2[*,1]
    plotsym, 8, 1.2, fill=0, thick=5.0
    djs_oplot, instr_gmag[fit]-instr_rmag[fit], sdss_gmag[fit]-instr_gmag[fit], psym=8, color='red'
    djs_oplot, !x.crange, gcoeff[0]*[1,1], line=2, thick=5.0
;   im_legend, ['\sigma_{g} = '+string(bresid,format='(F5.3)')], $
;     /left, /top, box=0
; r-band
    plotsym, 0, 1.0, /fill
    djs_plot, instr_rmag, sdss_rmag, psym=8, xsty=1, ysty=1, $
      xtitle='r_{LDSS3} (instrumental mag)', $
      ytitle='r_{SDSS} (AB mag)', xrange=[14,21.5], $
      yrange=[14,21.5], color='blue', position=pos2[*,0]
    plotsym, 8, 1.2, fill=0, thick=5.0
    djs_oplot, instr_rmag[fit], sdss_rmag[fit], psym=8, color='red'
    rcoeff = im_linefit(instr_rmag[fit],sdss_rmag[fit],coeff_fixed=[0,1],$
      coeff_guess=[-1.0,1.0],yerr=sdss_rmag_err[fit],chi2=rchi2,yfit=rfit,$
      coeff_err=rcoeff_err)
    rresid = djsig(sdss_rmag[fit]-rfit)
    djs_oplot, magaxis, poly(magaxis,rcoeff), line=0, thick=4
    im_legend, ['r_{AB} = '+strtrim(string(rcoeff[0],format='(F12.3)'),2)+$
      '+r_{LDSS3}*'+strtrim(string(rcoeff[1],format='(F12.2)'),2),$
      '\sigma_{r} = '+string(rresid,format='(F5.3)')], /right, /bottom, box=0
    im_legend, ['N = '+string(nfit,format='(I0)')], /left, /top, box=0
; residuals
    plotsym, 0, 1.0, /fill
    djs_plot, instr_gmag-instr_rmag, sdss_rmag-instr_rmag, /noerase, psym=8, xsty=1, ysty=1, $
      ytitle='r_{SDSS}-r_{LDSS3}', $
      xtitle='g_{LDSS3} - r_{LDSS3} (instrumental mag)', xrange=[-0.2,1.6], $
      yrange=[-0.2,0.7], color='blue', position=pos2[*,1]
    plotsym, 8, 1.2, fill=0, thick=5.0
    djs_oplot, instr_gmag[fit]-instr_rmag[fit], sdss_rmag[fit]-instr_rmag[fit], psym=8, color='red'
    djs_oplot, !x.crange, rcoeff[0]*[1,1], line=2, thick=5.0
;   im_legend, ['\sigma_{r} = '+string(rresid,format='(F5.3)')], $
;     /left, /top, box=0

    im_plotconfig, /psclose
    spawn, 'ps2pdf13 '+psfile+' '+repstr(psfile,'.ps','.pdf'), /sh
    
;   niceprint, gcoeff, gcoeff_err, rcoeff, rcoeff_err
    
    
stop    

return
end
    






;;function read_sdss3, sdssfile
;;
;;    data = djs_readlines(sdssfile)
;;    keep = where(strmatch(data,'#*') eq 0B)
;;    data = data[keep]
;;    ndata = n_elements(data)
;;    
;;    sdss = {star: '', ra: 0.0D, dec: 0.0D, g: 0.0D, gerr: 0.0D, r: 0.0D, rerr: 0.0D}
;;    sdss = replicate(sdss,ndata)
;;
;;    for ii = 0L, ndata-1L do begin
;;       sdss[ii].star = strmid(data[ii],5,18)
;;       sdss[ii].ra   = strmid(data[ii],25,10)
;;       sdss[ii].dec  = strmid(data[ii],35,10)
;;       sdss[ii].g    = strmid(data[ii],71,6)
;;       sdss[ii].gerr = strmid(data[ii],78,5)
;;       sdss[ii].r    = strmid(data[ii],84,6)
;;       sdss[ii].rerr = strmid(data[ii],91,5)
;;    
;;    endfor
;;
;;return, sdss
;;end
;;
;;pro ldss3_calibrate, findsdss=findsdss
;;; jm07jan27nyu - calibrate the LDSS3 images
;;;   (11:25:48.43, -02:44:48.4), 10.2'x11.2'
;;
;;    datapath = ldss3_path(/feb06)+'standards/'
;;    
;;    sdssfile = datapath+'sdss3.dat'
;;    if keyword_set(findsdss) then $
;;      spawn, "findsdss3 -c '11:25:48,-02:44:48.4', -b 11.0 -lc 3 -m 1000 > "+sdssfile
;;
;;    sdss = read_sdss3(sdssfile)
;;;   struct_print, sdss[0:10]
;;
;;    result1 = {image: '', mag_auto: 0.0D, magerr_auto: 0.0D, rmag: 0.0D, $
;;      rmagerr: 0.0D, gmag: 0.0D, gmagerr: 0.0D, airmass: 0.0D, exptime: 0.0D}
;;    
;;; loop through each SE catalog    
;;
;;    headlist = file_search(datapath+'ra.????_sdss_1_*_c?.head',count=nimage)
;;    catlist = repstr(headlist,'.head','.cat')
;;    imagelist = repstr(headlist,'.head','.fits')
;;
;;    for ii = 0L, nimage-1L do begin
;;       
;;       hdr = headfits(imagelist[ii])
;;       cat = mrdfits(catlist[ii],2,/silent)
;;       astrhdr = djs_readlines(headlist[ii])
;;       extast, astrhdr, astr
;;
;;       airmass = sxpar(hdr,'AIRMASS')
;;       exptime = sxpar(hdr,'EXPTIME')
;;       
;;       xy2ad, cat.xwin_image, cat.ywin_image, astr, aa, dd
;;       spherematch, aa, dd, sdss.ra, sdss.dec, 4.0/3600.0, $
;;         catmatch, sdssmatch, distance12, maxmatch=1
;;       sortindex = uniq(sdssmatch,sort(sdssmatch))
;;       sdssmatch = sdssmatch[sortindex]
;;       catmatch = catmatch[sortindex]
;;       nsdssmatch = n_elements(sdssmatch)
;;
;;       splog, 'Image '+file_basename(imagelist[ii])+': matched '+$
;;         string(nsdssmatch,format='(I0)')+' SDSS stars.'
;;;      struct_print, sdss[sdssmatch]
;;
;;       result = replicate(result1,nsdssmatch)
;;       result.image       = imagelist[ii]
;;       result.mag_auto    = cat[catmatch].mag_auto
;;       result.magerr_auto = cat[catmatch].magerr_auto
;;       result.rmag        = sdss[sdssmatch].r
;;       result.rmagerr     = sdss[sdssmatch].rerr
;;       result.gmag        = sdss[sdssmatch].g
;;       result.gmagerr     = sdss[sdssmatch].gerr
;;       result.airmass     = airmass
;;       result.exptime     = exptime
;;
;;       if strmatch(imagelist[ii],'*ra.????_sdss_1_r_c?*') then begin       
;;          if (n_elements(bigresult_r) eq 0L) then bigresult_r = result else $
;;            bigresult_r = [bigresult_r,result]
;;       endif else begin
;;          if (n_elements(bigresult_g) eq 0L) then bigresult_g = result else $
;;            bigresult_g = [bigresult_g,result]
;;       endelse
;;       
;;    endfor       
;;
;;stop    
;;    
;;; calibrate the r-band magnitudes    
;;        
;;    parinfo = replicate({value: 0.1D},3)
;;    expr = 'X + P[0] + P[1]*X + P[2]*(PRIVATE.GMAG-PRIVATE.RMAG)'
;;;   expr = 'X + P[0] + P[1]*X + P[2]*(PRIVATE.GMAG-PRIVATE.RMAG) + P[3]*PRIVATE.AIRMASS'
;;    params = mpfitexpr(expr,bigresult_r.mag_auto,bigresult_r.rmag,parinfo=parinfo,functargs=bigresult_r)
;;
;;    plot, bigresult_r.mag_auto, bigresult_r.rmag, ps=4, xsty=3, ysty=3
;;    plot, mpevalexpr(expr,bigresult_r.mag_auto,params,functargs=bigresult_r), $
;;      bigresult_r.rmag, ps=4, xsty=3, ysty=3, xr=[18,25.5], yr=[18,25.5]
;;       
;;stop       
;;
;;return
;;end
    
;;
;;
;;; g-band
;;    plotsym, 0, 1.0, /fill
;;    djs_plot, instr_gmag, gmag, psym=8, xsty=1, ysty=1, $
;;      charsize=2.0, xtitle='g_{LDSS3} (instrumental mag)', $
;;      ytitle='SDSS g (AB mag)', xrange=[14,23.5], $
;;      yrange=[14,23.5], color='blue', position=pos2[*,0]
;;    plotsym, 8, 1.2, fill=0, thick=5.0
;;    djs_oplot, instr_gmag[fit], gmag[fit], psym=8, color='red'
;;    gcoeff = im_linefit(instr_gmag[fit],gmag[fit],coeff_fixed=[0,1],$
;;      coeff_guess=[-1.0,1.0],yerr=gmag_err[fit],chi2=gchi2,yfit=gfit,$
;;      coeff_err=gcoeff_err)
;;    gresid = djsig(gmag[fit]-gfit)
;;;   gcoeff = robust_linefit(instr_gmag[fit],gmag[fit],ff,gsig,/bisect)
;;    djs_oplot, magaxis, poly(magaxis,gcoeff), line=0, thick=4
;;    im_legend, ['g_{AB} = '+strtrim(string(gcoeff[0],format='(F12.3)'),2)+$
;;      '+g_{LDSS3}*'+strtrim(string(gcoeff[1],format='(F12.2)'),2),$
;;      '\sigma_{g} = '+string(gresid,format='(F5.3)')],/right, /bottom, box=0
;;    im_legend, 'N = '+string(nfit,format='(I0)'), /left, /top, box=0
;;;   cc = get_kbrd(1)
;;; r-band
;;    plotsym, 0, 1.0, /fill
;;    djs_plot, instr_rmag, rmag, psym=8, xsty=1, ysty=1, $
;;      charsize=2.0, xtitle='LDSS3 r (instrumental mag)', $
;;      ytitle='SDSS r (AB mag)', xrange=[14,21.5], $
;;      yrange=[14,21.5], color='blue'
;;;   plotsym, 0, 1.2, fill=1, thick=6.0
;;;   djs_oplot, instr_rmag[fit], rmag[fit], psym=8
;;    plotsym, 8, 1.2, fill=0, thick=5.0
;;    djs_oplot, instr_rmag[fit], rmag[fit], psym=8, color='red'
;;    rcoeff = im_linefit(instr_rmag[fit],rmag[fit],coeff_fixed=[0,1],$
;;      coeff_guess=[-1.0,1.0],yerr=rmag_err[fit],chi2=rchi2,yfit=rfit,$
;;      coeff_err=rcoeff_err)
;;    rresid = djsig(rmag[fit]-rfit)
;;;   rcoeff = robust_linefit(instr_rmag[fit],rmag[fit],ff,rsig,/bisect)
;;    djs_oplot, magaxis, poly(magaxis,rcoeff), line=0, thick=4
;;    im_legend, ['r_{AB} = '+strtrim(string(rcoeff[0],format='(F12.3)'),2)+$
;;      '+r_{LDSS3}*'+strtrim(string(rcoeff[1],format='(F12.2)'),2),$
;;      '\sigma_{r} = '+string(rresid,format='(F5.3)')],/right, /bottom, box=0
;;    im_legend, 'N = '+string(nfit,format='(I0)'), /left, /top, box=0
