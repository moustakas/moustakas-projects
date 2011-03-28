pro bootes_kurucz_zeropoints
; jm11mar23ucsd - 
    
    common calibrate_bootes, sdss1, bootes1

    photodir = getenv('IM_ARCHIVE_DIR')+'/bootes/'
    iopath = ages_path(/mycat)+'kurucz/'
    psfile = iopath+'qa_bootes_kurucz_zeropoints.ps'
    
; read and match the Bootes and SDSS catalogs    
    if (n_elements(sdss1) eq 0L) then begin
       photofile = photodir+'sdss_stars_dr72.fits.gz'
       sdss1 = mrdfits(photofile,1)
    endif

    if (n_elements(bootes1) eq 0L) then begin
       photofile = photodir+'bootes_sdss_stars_dr72.fits.gz'
       if (file_test(photofile) eq 0) then begin
          ver = '2010b'
          iband = mrdfits(photodir+ver+'/bootes_I.fits.gz',1)

          spherematch, iband.alpha_j2000, iband.delta_j2000, $
            sdss1.ra, sdss1.dec, 1.0/3600.0, m1, m2
          srt = sort(m1) & m1 = m1[srt] & m2 = m2[srt]

          bootes1 = iband[m1]
          tags = tag_names(bootes1)
          bootes1 = im_struct_trimtags(temporary(bootes1),$
            select=tags,newtags='I_'+tags)
          
          bands = ['u','Bw','R','z']
          for jj = 0, n_elements(bands)-1 do begin
             phot1 = mrdfits(photodir+ver+'/bootes_'+bands[jj]+'.fits.gz',rows=m1,1)
             tags = tag_names(phot1)
             phot1 = im_struct_trimtags(phot1,select=tags,$
               newtags=bands[jj]+'_'+tags)
             bootes1 = struct_addtags(temporary(bootes1),phot1)
          endfor
          im_mwrfits, bootes1, repstr(photofile,'.gz',''), /clobber
       endif else bootes1 = mrdfits(photofile,1)   
    endif

; spherematch, convert to maggies, and select a fiducial sample of
; objects with good multiband photometry
    spherematch, bootes1.i_alpha_j2000, bootes1.i_delta_j2000, $
      sdss1.ra, sdss1.dec, 1.0/3600.0, m1, m2
    good = where(sdss1[m2].psfflux[2] gt 0.0)
    m1 = m1[good] & m2 = m2[good]
    primus_sdss_to_maggies, sdssmaggies, sdssivarmaggies, calib=sdss1[m2], flux='psf'
    rmag = -2.5*alog10(sdssmaggies[2,*])

    cut1 = where((rmag gt 19.0) and (rmag lt 20.5) and $
      (bootes1.i_flag_duplicate eq 0) and (bootes1.i_flag_subfield eq 1) and $
      (bootes1.i_class_star gt 0.95) and $
      (bootes1.u_mag_psf gt 0.0) and (bootes1.u_mag_psf lt 90.0) and $
      (bootes1.bw_mag_psf gt 0.0) and (bootes1.bw_mag_psf lt 90.0) and $
      (bootes1.r_mag_psf gt 0.0) and (bootes1.r_mag_psf lt 90.0) and $
      (bootes1.i_mag_psf gt 0.0) and (bootes1.i_mag_psf lt 90.0) and $
      (bootes1.z_mag_psf gt 0.0) and (bootes1.z_mag_psf lt 90.0),ncut1)

; fit Kurucz models to all the stars and compare the synthesized vs
; observed photometry
    splog, 'Fitting '+string(ncut1,format='(I0)')+' objects...'
    kall = fit_kurucz_models(sdssmaggies[*,cut1],sdssivarmaggies[*,cut1],$
      filterlist=sdss_filterlist())
    cut2 = where(kall.kurucz_chi2min lt 3.0,nstar)
    splog, 'N = ', nstar

    kthese = kall[cut2]
    sdss = sdss1[m2[cut1[cut2]]]
    bootes = bootes1[m1[cut1[cut2]]]
    sdssmaggies = sdssmaggies[*,cut1[cut2]]
    sdssivarmaggies = sdssivarmaggies[*,cut1[cut2]]

; pick the bands    
;   band = ['u*','g','r','i*','z']
;   mrange = [[19.5,23.0],[19.0,21.5],[19.0,21.0],[18.5,21.0],[18.5,21.0]]
    band = ['U','B_{W}','R','I','z']
    mrange = [[19.5,23.5],[19.0,22.0],[18.5,21.0],[18.0,21.0],[18.5,21.0]]
    rrange = 0.99*[-1,1]
    nband = n_elements(band)

    bootes_to_maggies, bootes, maggies, ivarmaggies, /psf, $
      filterlist=synthmaggies_filterlist, /nozpoffset 
    indx = lindgen(5)
    maggies = maggies[indx,*]   ; exclude JHKs + irac
    ivarmaggies = ivarmaggies[indx,*]
    synthmaggies_filterlist = synthmaggies_filterlist[indx]

; synthesize magnitudes
    synthmaggies = dblarr(nband,nstar)
    for ii = 0L, nstar-1L do synthmaggies[*,ii] = reform(k_project_filters($
      k_lambda_to_edges(kthese[ii].lambda),kthese[ii].spec,$
      filterlist=synthmaggies_filterlist))

; iteratively solve for the best zeropoint offsets
    mask = intarr(nband,nstar)+1 ; all are good
    for iter = 0, 9 do begin
       norm = total(mask*ivarmaggies*synthmaggies^2,2,/double)
       dscale = total(mask*ivarmaggies*maggies*synthmaggies,2,/double)/$
         (norm+(norm eq 0))*(norm ne 0)
       chi2 = total(mask*ivarmaggies*(maggies-rebin(reform(dscale,nband,1),nband,nstar)*$
         synthmaggies)^2,/double,1)/total(maggies gt 0,1)
       mask = mask and rebin(reform(chi2 lt weighted_quantile(chi2,quant=0.95),1,nstar),nband,nstar)
       dzpt = +2.5*alog10(dscale)
    endfor

; make the plot    
    im_plotconfig, 6, pos1, psfile=psfile, xmargin=[1.4,0.3], $
      height=[4.5,3.0], width=6.8
    for ii = 0, nband-1 do begin
       good = where((maggies[ii,*] gt 0.0),ngood)
       xx = reform(-2.5*alog10(maggies[ii,good]))
       yy = reform(-2.5*alog10(synthmaggies[ii,good]))
       djs_plot, xx, yy, position=pos1[*,0], xsty=1, ysty=1, $
         xrange=mrange[*,ii], yrange=mrange[*,ii], $
         xtitle='', ytitle=band[ii]+' (SDSS synthesized, AB mag)', $
         xtickname=replicate(' ',10), psym=6
       djs_oplot, !x.crange, !y.crange, line=0, thick=3, color='red'
       im_legend, ['\Delta'+'Zpt = '+strtrim(string(dzpt[ii],format='(F12.3)'),2),$
         '<\Delta'+band[ii]+'> = '+im_string_stats(yy-xx,sigrej=2.5,ndec=3)], $
         /left, /top, box=0, margin=0
       djs_plot, xx, yy-xx, position=pos1[*,1], /noerase, xsty=1, ysty=1, $
         xrange=mrange[*,ii], yrange=rrange, psym=6, $
         xtitle=band[ii]+' (BOOTES/bootes, AB mag)', $
         ytitle='Residuals (AB mag)'
       djs_oplot, !x.crange, [0,0], line=0, thick=3, color='red'
       djs_oplot, !x.crange, dzpt[ii]*[1,1], line=0, thick=3, color='blue'
    endfor
    im_plotconfig, psfile=psfile, /psclose, /gzip
;   spawn, 'rsync -auv '+psfile+'.gz ~/', /sh

    niceprint, dzpt, synthmaggies_filterlist
    
stop    
    
return
end
