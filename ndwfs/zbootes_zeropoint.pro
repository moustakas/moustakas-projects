pro zbootes_zeropoint
; jm09sep04ucsd - calibrate the ZBOOTES photometry
    
    common calibrate_zbootes, zbootes1, sdss1

    analysis_path = ages_path(/analysis)
    photodir = getenv('RESEARCHPATH')+'/data/ndwfs/'

; read and match the ZBOOTES and SDSS catalogs    
    if (n_elements(zbootes1) eq 0L) then begin
       splog, 'Reading '+analysis_path+'zbootes-cat.fits.gz'
       zbootes1 = mrdfits(analysis_path+'zbootes-cat.fits.gz',1)
    endif

    if (n_elements(sdss1) eq 0) then begin
       photofile = photodir+'NDWFS_SDSS_stars.fits.gz'
       sdss1 = mrdfits(photofile,1)
    endif

    sdss_filterlist = sdss_filterlist()
    zbootes_filterlist = zbootes_filterlist()

; spherematch, convert to maggies, and select a fiducial sample of
; objects with good multiband photometry
    spherematch, zbootes1.alpha_j2000, zbootes1.delta_j2000, $
      sdss1.ra, sdss1.dec, 1.0/3600.0, m1, m2

    sdss = ndwfs_zpt_sdsscat2mag(sdss1[m2],maggies=sdssmaggies,$
      ivarmaggies=sdssivarmaggies)
    zbootes = zpt_zbootescat2mag(zbootes1[m1])

    these = where((zbootes.z gt 0.0) and $
      (sdss.u gt 0.0) and (sdss.g gt 0.0) and $
      (sdss.r gt 0.0) and (sdss.i gt 0.0) and (sdss.z gt 0.0) and $
      (sdss.r ge 18.0) and (sdss.r le 20.5) and $
      (sdss.g-sdss.r le 1.1))
    sdss = sdss[these]
    zbootes = zbootes[these]
    sdssmaggies = sdssmaggies[*,these]
    sdssivarmaggies = sdssivarmaggies[*,these]

; fit Kurucz models to all the stars and compare the synthesized vs
; observed z-band photometry

    kall = fit_kurucz_models(sdssmaggies,sdssivarmaggies,$
      filterlist=sdss_filterlist)
    these = where(kall.kurucz_chi2min lt 3.0,nstar)
    kthese = kall[these]

    zsynth = fltarr(nstar)
    for ii = 0L, nstar-1L do zsynth[ii] = reform(k_project_filters($
      k_lambda_to_edges(kthese[ii].lambda),kthese[ii].spec,$
      filterlist=zbootes_filterlist))

    psfile = photodir+'zbootes_kurucz_zeropoint.ps'
    im_plotconfig, 6, pos1, psfile=psfile, xmargin=[1.4,0.3], $
      height=[4.5,3.0], width=6.8

    mrange = [17.0,21.5]
    rrange = 0.39*[-1,1]

    band = ['z']
    for ii = 0, n_elements(band)-1 do begin
       magtag = tag_indx(zbootes[0],band[ii])
       xx = zbootes[these].(magtag)
       yy = reform(-2.5*alog10(zsynth))
       djs_plot, xx, yy, position=pos1[*,0], $
         xsty=1, ysty=1, xrange=mrange, yrange=mrange, $
         psym=4, xtitle='', ytitle=band[ii]+' (SDSS synthesized, AB mag)', $
         xtickname=replicate(' ',10)
       djs_oplot, !x.crange, !y.crange, line=0, thick=3, color='red'
;      im_legend, band[ii]+'_{SDSS}-'+band[ii]+'_{ZBOOTES} = '+$
       im_legend, '<\Delta'+band[ii]+'> = '+$
         im_string_stats(yy-xx,sigrej=3.0), /left, /top, box=0
       djs_plot, xx, yy-xx, position=pos1[*,1], /noerase, xsty=1, ysty=1, $
         xrange=mrange, yrange=rrange, psym=4, xtitle=band[ii]+' (ZBOOTES, AB mag)', $
         ytitle='Residuals (AB mag)'
       djs_oplot, !x.crange, [0,0], line=0, thick=3, color='red'
    endfor
       
    im_plotconfig, psfile=psfile, /psclose, /pdf, /pskeep
    
return
end
