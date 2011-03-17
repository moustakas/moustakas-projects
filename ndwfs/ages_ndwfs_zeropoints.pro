pro zpt_hogg_scatterplot, mag, deltamag, band, input=input
; make a scatterplot
    case band of
       'u': xrange = [17.5,26]
       'g': xrange = [15.5,24.5]
       'r': xrange = [15,22.5]
       'i': xrange = [14.5,22]
       'z': xrange = [14,22]
       'Bw': xrange = [16.5,25.5]
       'R': xrange = [15,22]
       'I': xrange = [15.0,21.0]
       'zBootes': xrange = [14,21.5]
       'J': xrange = [14.5,21.0]
       'K': xrange = [14.5,21.0]
       'Ks': xrange = [14.5,21.0]
    endcase
    
    hogg_scatterplot, mag, deltamag, xsty=3, ysty=1, /internal, $
;     xrange=[floor(min(mag)),ceil(max(mag))], $
      xrange=xrange, $
      yrange=[-1,1], /outliers, $
      xtitle=textoidl(band+' (observed, AB mag)'), $
      ytitle=textoidl('\Delta'+band+' (observed minus synthesized, AB mag)')
    djs_oplot, !x.crange, [0,0], line=0, thick=5
    med = im_medxbin(mag,deltamag,0.2,minpts=10)
    djs_oplot, med.binctr, med.medy, line=0, thick=5, color='blue'
    djs_oplot, med.binctr, med.sigy84, line=5, thick=5, color='blue'
    djs_oplot, med.binctr, med.sigy16, line=5, thick=5, color='blue'
    im_legend, 'Input: '+input, /left, /top, box=0
    im_legend, '<\Delta'+band+'> = '+im_string_stats(deltamag,sigrej=3.0), $
      /left, /bottom, box=0

return
end
    
pro ages_ndwfs_zeropoints, recompute=recompute
; jm09apr30nyu - check the relative and absolute zeropoints of the
;   NDWFS photometry by using the SDSS photometry of the AGES galaxies
; jm09sep03ucsd - rewritten

    common ages_ndwfs_zeropoints, sdss1, ndwfs1, kcorr1, zband1

    analysis_path = ages_path(/analysis)
    photodir = getenv('RESEARCHPATH')+'/data/ndwfs/'

; SDSS galaxy photometry    
    if (n_elements(sdss1) eq 0L) then begin
       splog, 'Reading '+analysis_path+'ages.sdss.phot.dr72.fits.gz'
       sdss1 = mrdfits(analysis_path+'ages.sdss.phot.dr72.fits.gz',1)
    endif    

; NDWFS/DR3 catalogs    
    if (n_elements(ndwfs1) eq 0L) then begin
       bands = ['Bw','R','I','K']
       tags = ['ndwfs_field','ndwfs_name','alpha_j2000','delta_j2000','flux_radius',$
         'class_star','mag_auto','magerr_auto','mag_aper*','magerr_aper*',$
         'flags','imaflags_iso','flag_duplicate','flag_splitmatch']
       files = analysis_path+'catalog.ndwfs'+['b','r','i','k']+'.fits.gz'
       for jj = 0, n_elements(bands)-1 do begin
          splog, 'Reading '+files[jj]
          cat1 = mrdfits(files[jj],1)
          cat1 = struct_trimtags(cat1,select=bands[jj]+'_'+tags)
          if (jj eq 0) then ndwfs1 = cat1 else ndwfs1 = $
            struct_addtags(temporary(ndwfs1),temporary(cat1))
       endfor
    endif

; zBootes
    if (n_elements(zband1) eq 0L) then begin
       splog, 'Reading '+analysis_path+'zbootes-cat.fits.gz'
       allzband1 = mrdfits(analysis_path+'zbootes-cat.fits.gz',1)

       searchrad = 1.0          ; [arcsec]
       spherematch, ndwfs1.i_alpha_j2000, ndwfs1.i_delta_j2000, $
         allzband1.alpha_j2000, allzband1.delta_j2000, searchrad/3600.0, m1, m2
       zband1 = im_empty_structure(allzband1[0],empty_value=-999.0,$
         ncopies=n_elements(ndwfs1))
       zband1[m1] = allzband1[m2]
    endif

; Eisenstein's k-corrections (for the redshifts)
    if (n_elements(kcorr1) eq 0L) then begin
       splog, 'Reading '+analysis_path+'catalog.kcorr.v3.fits.gz'
       kcorr1 = mrdfits(analysis_path+'catalog.kcorr.v3.fits.gz',1)
    endif

; select a well-defined sample
    cut1 = where((sdss1.ra gt -900.0) and (kcorr1.kcorr_z gt 0.03) and $
      (kcorr1.kcorr_z lt 0.5))
    sdss = ndwfs_zpt_sdsscat2mag(sdss1[cut1],flux='petro')
    ndwfs = zpt_ndwfscat2mag(ndwfs1[cut1])
    cut2 = where((sdss.u gt 0.0) and (sdss.g gt 0.0) and $
      (sdss.r gt 0.0) and (sdss.i gt 0.0) and (sdss.z gt 0.0) and $
      (sdss.r lt 21.0) and (ndwfs.bw_flags eq 0) and $
      (ndwfs.r_flags eq 0) and (ndwfs.i_flags eq 0))
    keep = cut1[cut2]
    nkeep = n_elements(keep)
    splog, 'Number of objects', nkeep

    sdss = sdss1[keep]
    ndwfs = ndwfs1[keep]
    zband = zband1[keep]
    kcorr = kcorr1[keep]

    sdss_filterlist = sdss_filterlist()
    out_filterlist = [ndwfs_filterlist(),zbootes_filterlist()]

; fit the SDSS photometry
    sdsscat = ndwfs_zpt_sdsscat2mag(sdss,maggies=maggies,$
      ivarmaggies=ivarmaggies,flux='petro')

    kk = im_kcorrect(kcorr.kcorr_z,maggies,ivarmaggies,$
      sdss_filterlist,out_filterlist,chi2=chi2,coeffs=coeffs,$
      bestmaggies=bestmaggies,synth_outmaggies_obs=maggies_synth,$
      mass=mass,/silent,band_shift=0.0,/reset)
    these = where(chi2 lt 3.0,nthese)
    splog, 'Final number of objects', nthese

; make a K-correct QAplot    
    res = replicate({z: 0.0, chi2: 0.0, mass: 0.0, maggies: fltarr(5), $
      ivarmaggies: fltarr(5), bestmaggies: fltarr(5), $
      coeffs: fltarr(5)},nkeep)
    res.z = kcorr.kcorr_z
    res.chi2 = chi2
    res.mass = alog10(mass)
    res.maggies = maggies
    res.ivarmaggies = ivarmaggies
    res.bestmaggies = bestmaggies
    res.coeffs = coeffs
    psfile = photodir+'ndwfs_ages_kcorr_qaplot.ps'
;   kcorrect_qaplot, res, sdss_filterlist, psfile=psfile, /clobber

; make the final zeropoint plot    
    psfile = photodir+'ndwfs_ages_zeropoints.ps'
    im_plotconfig, 6, pos1, psfile=psfile, xmargin=[1.4,0.3], $
      height=[4.5,3.0], width=6.8

    mrange = [[18.5,24.0],[17.5,22.0],[17.0,21.5],$
      [16.0,21.0],[16.0,21.0]]
    rrange = 0.89*[-1,1]

    ndwfscat = zpt_ndwfscat2mag(ndwfs)
    zbandcat = zpt_zbootescat2mag(zband)
    ndwfscat = struct_addtags(ndwfscat,zbandcat)

    band = ['Bw','R','I','K','z']
    for ii = 0, n_elements(band)-1 do begin
       magtag = tag_indx(ndwfscat[0],band[ii])
       xx = ndwfscat[these].(magtag) ; observed, AB
       yy = reform(-2.5*alog10(maggies_synth[ii,these])) ; synthesized, AB
       good = where(xx gt 0.0)
       xx = xx[good] & yy = yy[good]
       djs_plot, xx, yy, position=pos1[*,0], $
         xsty=1, ysty=1, xrange=mrange[*,ii], yrange=mrange[*,ii], $
         psym=4, xtitle='', ytitle=band[ii]+' (SDSS/AGES synthesized, AB mag)', $
         xtickname=replicate(' ',10)
       djs_oplot, !x.crange, !y.crange, line=0, thick=3, color='red'
       im_legend, '<\Delta'+band[ii]+'> = '+$
         im_string_stats(yy-xx,sigrej=3.0), /left, /top, box=0
       djs_plot, xx, yy-xx, position=pos1[*,1], /noerase, xsty=1, ysty=1, $
         xrange=mrange[*,ii], yrange=rrange, psym=4, xtitle=band[ii]+' (NDWFS, AB mag)', $
         ytitle='Residuals (mag)'
       djs_oplot, !x.crange, [0,0], line=0, thick=3, color='red'
    endfor
       
    im_plotconfig, psfile=psfile, /psclose, /pdf, /pskeep
    
stop    
    
return
end
    


;;    ages1 = mrdfits(analysis_path+'ages_merged_catalogs.fits.gz',1,/silent)
;;    keep = where(ages1.main_flag and ages1.sdss_match and $
;;      (ages1.z gt 0.05) and (ages1.z lt 0.75),ngal)
;;    ages = ages1[keep]
;;
;;; only keep relatively bright SDSS galaxies    
;;    sdss_to_maggies, mm, vv, calibobj=ages, flux='model'
;;    keep = where(-2.5*alog10(mm[2,*]) lt 20.5)
;;;   plot, ages[keep].z, -2.5*alog10(mm[2,keep]), ps=4, xsty=3, ysty=3
;;    ages = ages[keep]
;;
;;; now do K-corrections; start by converting to maggies    
;;    ages_to_maggies, ages, ndwfsmaggies, ndwfsmaggies_ivar, $
;;      filterlist=ndwfsfilters
;;    ages_to_maggies, ages, sdssmaggies, sdssmaggies_ivar, $
;;      filterlist=sdssfilters, /sdss
;;;   ages_to_maggies, ages, twomassmaggies, twomassmaggies_ivar, $
;;;     filterlist=twomassfilters, /twomass
;;    
;;    outfile = analysis_path+'ndwfs_zeropoints_data.fits'
;;    if keyword_set(recompute) then begin
;;       splog, 'Computing K-corrections...'
;;; test1: compute ugriz and BwRIzJKsK from the observed ugriz photometry 
;;       kk = im_kcorrect(ages.z,sdssmaggies,sdssmaggies_ivar,sdssfilters,$
;;         band_shift=0.0,chi2=chi2,out_filterlist=ndwfsfilters,$
;;         coeffs=coeffs1,/silent,synth_inmaggies=synth_sdssmaggies1,$
;;         synth_outmaggies=synth_ndwfsmaggies1,mass=mass1,/reset)
;;; test2: compute ugriz and BwRIzJKsK from the observed BwRI photometry 
;;       in_ndwfsmaggies_ivar = ndwfsmaggies_ivar
;;       in_ndwfsmaggies_ivar[3:6,*] = 0.0
;;       kk = im_kcorrect(ages.z,ndwfsmaggies,in_ndwfsmaggies_ivar,ndwfsfilters,$
;;         band_shift=0.0,chi2=chi2,out_filterlist=sdssfilters,$
;;         coeffs=coeffs2,/silent,synth_inmaggies=synth_ndwfsmaggies2,$
;;         synth_outmaggies=synth_sdssmaggies2,mass=mass2,/reset)
;;; test3: compute ugriz and BwRIzJKsK from the observed BwRIz photometry 
;;       in_ndwfsmaggies_ivar = ndwfsmaggies_ivar
;;       in_ndwfsmaggies_ivar[4:6,*] = 0.0
;;       kk = im_kcorrect(ages.z,ndwfsmaggies,in_ndwfsmaggies_ivar,ndwfsfilters,$
;;         band_shift=0.0,chi2=chi2,out_filterlist=sdssfilters,$
;;         coeffs=coeffs3,/silent,synth_inmaggies=synth_ndwfsmaggies3,$
;;         synth_outmaggies=synth_sdssmaggies3,mass=mass3,/reset)
;;; test4: compute ugriz and BwRIzJKsK from the observed BwRIzK photometry 
;;       in_ndwfsmaggies_ivar = ndwfsmaggies_ivar
;;       in_ndwfsmaggies_ivar[4:5,*] = 0.0
;;       kk = im_kcorrect(ages.z,ndwfsmaggies,in_ndwfsmaggies_ivar,ndwfsfilters,$
;;         band_shift=0.0,chi2=chi2,out_filterlist=sdssfilters,$
;;         coeffs=coeffs4,/silent,synth_inmaggies=synth_ndwfsmaggies4,$
;;         synth_outmaggies=synth_sdssmaggies4,mass=mass4,/reset)
;;; write out
;;       out = {ngal: ngal, $
;;         synth_sdssmaggies1: synth_sdssmaggies1, $
;;         synth_sdssmaggies2: synth_sdssmaggies2, $
;;         synth_sdssmaggies3: synth_sdssmaggies3, $
;;         synth_sdssmaggies4: synth_sdssmaggies4, $
;;         synth_ndwfsmaggies1: synth_ndwfsmaggies1, $
;;         synth_ndwfsmaggies2: synth_ndwfsmaggies2, $
;;         synth_ndwfsmaggies3: synth_ndwfsmaggies3, $
;;         synth_ndwfsmaggies4: synth_ndwfsmaggies4, $
;;         mass1: mass1, mass2: mass2, mass3: mass3, mass4: mass4}
;;       im_mwrfits, out, outfile
;;    endif else begin
;;       out = mrdfits(outfile+'.gz',1)
;;       if (ngal ne out.ngal) then begin
;;          splog, 'Call this routine with /RECOMPUTE'
;;          return
;;       endif
;;       synth_sdssmaggies1  = out.synth_sdssmaggies1
;;       synth_sdssmaggies2  = out.synth_sdssmaggies2
;;       synth_sdssmaggies3  = out.synth_sdssmaggies3
;;       synth_sdssmaggies4  = out.synth_sdssmaggies4
;;       synth_ndwfsmaggies1 = out.synth_ndwfsmaggies1
;;       synth_ndwfsmaggies2 = out.synth_ndwfsmaggies2
;;       synth_ndwfsmaggies3 = out.synth_ndwfsmaggies3
;;       synth_ndwfsmaggies4 = out.synth_ndwfsmaggies4
;;       mass1 = out.mass1 & mass2 = out.mass2
;;       mass3 = out.mass3 & mass4 = out.mass4
;;    endelse
;;       
;;    ssband = ['u','g','r','i','z']
;;;   nnband = ['Bw','R','I','zBootes']
;;;   nnband = ['Bw','R','I','zBootes','K']
;;    nnband = ['Bw','R','I','z','J','Ks','K']
;;
;;; ###test1###
;;    psfile = analysis_path+'ndwfs_zeropoints_ugriz.ps'
;;    im_plotconfig, 0, psfile=psfile
;;; compare ugriz observed/synthesized from ugriz input
;;    for ii = 0, n_elements(ssband)-1 do begin
;;       good = where((sdssmaggies[ii,*] gt 0.0) and (synth_sdssmaggies1[ii,*] gt 0.0))
;;       mag = -2.5*alog10(sdssmaggies[ii,good])
;;       deltamag = -2.5*alog10(sdssmaggies[ii,good]/synth_sdssmaggies1[ii,good])
;;       zpt_hogg_scatterplot, mag, deltamag, ssband[ii], input='ugriz'
;;    endfor
;;; compare BwRIzJKsK observed/synthesized from ugriz input
;;    for ii = 0, n_elements(nnband)-1 do begin
;;       good = where((ndwfsmaggies[ii,*] gt 0.0) and (synth_ndwfsmaggies1[ii,*] gt 0.0))
;;       mag = -2.5*alog10(ndwfsmaggies[ii,good])
;;       deltamag = -2.5*alog10(ndwfsmaggies[ii,good]/synth_ndwfsmaggies1[ii,good])
;;       zpt_hogg_scatterplot, mag, deltamag, nnband[ii], input='ugriz'
;;    endfor
;;    im_plotconfig, psfile=psfile, /psclose, /gzip
;;
;;; ###test2###
;;    psfile = analysis_path+'ndwfs_zeropoints_BwRI.ps'
;;    im_plotconfig, 0, psfile=psfile
;;; compare ugriz observed/synthesized from BwRI input
;;    for ii = 0, n_elements(ssband)-1 do begin
;;       good = where((sdssmaggies[ii,*] gt 0.0) and (synth_sdssmaggies2[ii,*] gt 0.0))
;;       mag = -2.5*alog10(sdssmaggies[ii,good])
;;       deltamag = -2.5*alog10(sdssmaggies[ii,good]/synth_sdssmaggies2[ii,good])
;;       zpt_hogg_scatterplot, mag, deltamag, ssband[ii], input='BwRI'
;;    endfor
;;; compare BwRIzJKsK observed/synthesized from BwRIzJKsK input
;;    for ii = 0, n_elements(nnband)-1 do begin
;;       good = where((ndwfsmaggies[ii,*] gt 0.0) and (synth_ndwfsmaggies2[ii,*] gt 0.0))
;;       mag = -2.5*alog10(ndwfsmaggies[ii,good])
;;       deltamag = -2.5*alog10(ndwfsmaggies[ii,good]/synth_ndwfsmaggies2[ii,good])
;;       zpt_hogg_scatterplot, mag, deltamag, nnband[ii], input='BwRI'
;;    endfor
;;    im_plotconfig, psfile=psfile, /psclose, /gzip
;;
;;; ###test3###
;;    psfile = analysis_path+'ndwfs_zeropoints_BwRIz.ps'
;;    im_plotconfig, 0, psfile=psfile
;;; compare ugriz observed/synthesized from BwRIzJKsK input
;;    for ii = 0, n_elements(ssband)-1 do begin
;;       good = where((sdssmaggies[ii,*] gt 0.0) and (synth_sdssmaggies3[ii,*] gt 0.0))
;;       mag = -2.5*alog10(sdssmaggies[ii,good])
;;       deltamag = -2.5*alog10(sdssmaggies[ii,good]/synth_sdssmaggies3[ii,good])
;;       zpt_hogg_scatterplot, mag, deltamag, ssband[ii], input='BwRIz'
;;    endfor
;;; compare BwRIzJKsK observed/synthesized from BwRIzJKsK input
;;    for ii = 0, n_elements(nnband)-1 do begin
;;       good = where((ndwfsmaggies[ii,*] gt 0.0) and (synth_ndwfsmaggies3[ii,*] gt 0.0))
;;       mag = -2.5*alog10(ndwfsmaggies[ii,good])
;;       deltamag = -2.5*alog10(ndwfsmaggies[ii,good]/synth_ndwfsmaggies3[ii,good])
;;       zpt_hogg_scatterplot, mag, deltamag, nnband[ii], input='BwRIz'
;;    endfor
;;    im_plotconfig, psfile=psfile, /psclose, /gzip
;;
;;; ###test4###
;;    psfile = analysis_path+'ndwfs_zeropoints_BwRIzK.ps'
;;    im_plotconfig, 0, psfile=psfile
;;; compare ugriz observed/synthesized from BwRIzJKsK input
;;    for ii = 0, n_elements(ssband)-1 do begin
;;       good = where((sdssmaggies[ii,*] gt 0.0) and (synth_sdssmaggies4[ii,*] gt 0.0))
;;       mag = -2.5*alog10(sdssmaggies[ii,good])
;;       deltamag = -2.5*alog10(sdssmaggies[ii,good]/synth_sdssmaggies4[ii,good])
;;       zpt_hogg_scatterplot, mag, deltamag, ssband[ii], input='BwRIzK'
;;    endfor
;;; compare BwRIzJKsK observed/synthesized from BwRIzK input
;;    for ii = 0, n_elements(nnband)-1 do begin
;;       good = where((ndwfsmaggies[ii,*] gt 0.0) and (synth_ndwfsmaggies4[ii,*] gt 0.0))
;;       mag = -2.5*alog10(ndwfsmaggies[ii,good])
;;       deltamag = -2.5*alog10(ndwfsmaggies[ii,good]/synth_ndwfsmaggies4[ii,good])
;;       zpt_hogg_scatterplot, mag, deltamag, nnband[ii], input='BwRIzK'
;;    endfor
;;    im_plotconfig, psfile=psfile, /psclose, /gzip
;;
;;; mass comparison    
;;    psfile = analysis_path+'ndwfs_zeropoints_mass.ps'
;;    im_plotconfig, 0, psfile=psfile, charsize=2.0
;;    hogg_scatterplot, alog10(mass1), alog10(mass1/mass2), xrange=[7.5,12.5], $
;;      yrange=[-1.0,1.0], /xsty, /ysty, /outliers, /internal, $
;;      xtitle=textoidl('log_{10} (M_{ugriz}/M_{\odot})'), $
;;      ytitle=textoidl('log_{10} (M_{ugriz}/M_{BwRI})')
;;    djs_oplot, !x.crange, [0,0], line=0, thick=5
;;    im_legend, '<\Delta(M)> = '+im_string_stats(alog10(mass1/mass2),$
;;      sigrej=3.0), /left, /top, box=0
;;
;;    hogg_scatterplot, alog10(mass1), alog10(mass1/mass3), xrange=[7.5,12.5], $
;;      yrange=[-1.0,1.0], /xsty, /ysty, /outliers, /internal, $
;;      xtitle=textoidl('log_{10} (M_{ugriz}/M_{\odot})'), $
;;      ytitle=textoidl('log_{10} (M_{ugriz}/M_{BwRIz})')
;;    djs_oplot, !x.crange, [0,0], line=0, thick=5
;;    im_legend, '<\Delta(M)> = '+im_string_stats(alog10(mass1/mass3),$
;;      sigrej=3.0), /left, /top, box=0
;;
;;    hogg_scatterplot, alog10(mass1), alog10(mass1/mass4), xrange=[7.5,12.5], $
;;      yrange=[-1.0,1.0], /xsty, /ysty, /outliers, /internal, $
;;      xtitle=textoidl('log_{10} (M_{ugriz}/M_{\odot})'), $
;;      ytitle=textoidl('log_{10} (M_{ugriz}/M_{BwRIzK})')
;;    djs_oplot, !x.crange, [0,0], line=0, thick=5
;;    im_legend, '<\Delta(M)> = '+im_string_stats(alog10(mass1/mass4),$
;;      sigrej=3.0), /left, /top, box=0
;;    im_plotconfig, psfile=psfile, /psclose, /gzip
;;    
;;; ###########################################################################    
;;; intercompare the NDWFS K, FLAMEX JKs, and 2MASS JKs photometry; note
;;; that AGES_NDWFS_KBAND contains a more general comparison
;;    psfile = analysis_path+'ndwfs_zeropoints_kband.ps'
;;    im_plotconfig, 6, pos1, psfile=psfile, xmargin=[1.4,0.2], $
;;      height=[4.5,3.0]
;;    indx1 = where((ages1.phot_k gt 0.0) and (ages1.k_m_ext gt 0.0) and $
;;      (ages1.phot_flamj gt 0.0) and (ages1.phot_flamk gt 0.0))
;;;   indx2 = where((ages1.phot_flamj gt 0.0) and (ages1.phot_flamk gt 0.0) and $
;;;     (ages1.k_m_ext gt 0.0))
;;
;;    mrange1 = [11.0,15.0]
;;    mrange2 = [12.0,16.2]
;;    rrange = 0.95*[-1,1]
;;    
;;; mean K-Ks color of galaxies (see KMINUSKS); in AB magnitudes K-Ks is
;;; +0.07 across a range of star-formation histories (assuming an age of
;;; ~12 Gyr)
;;    kmks = +0.10 ; [Vega]
;;    strkmks = string(kmks,format='(F4.2)')
;;    
;;; NDWFS vs 2MASS
;;    xx = ages1[indx1].k_m_ext
;;    yy = ages1[indx1].phot_k-kmks
;;    djs_plot, xx, yy, position=pos1[*,0], $
;;      xsty=1, ysty=1, xrange=mrange1, yrange=mrange1, psym=4, $
;;      xtitle='', ytitle='K-'+strkmks+' (NDWFS, Vega mag)', xtickname=replicate(' ',10)
;;    djs_oplot, !x.crange, !y.crange, line=0, thick=3
;;    im_legend, '<\Delta'+'K_{s}> = '+im_string_stats(yy-xx), $
;;      /left, /top, box=0
;;    im_legend, '<K-K_{s}>_{Vega}=+'+strkmks+' for galaxies', /right, /bottom, box=0, margin=0
;;    djs_plot, xx, yy-xx, $
;;      position=pos1[*,1], /noerase, xsty=1, ysty=1, xrange=mrange1, $
;;      yrange=rrange, psym=4, xtitle='K_{s} (2MASS, Vega mag)', $
;;      ytitle='\Delta'+'K_{s} (Vega mag)'
;;    djs_oplot, !x.crange, [0,0], line=0, thick=3
;;;   kaxis = im_array(10,20,0.05)
;;;   djs_oplot, kaxis, kaxis+0.07, line=5, thick=3, color='red' ; <K-Ks>=+kmks for galaxies
;;; K-FLAMEX vs 2MASS
;;    xx = ages1[indx1].k_m_ext
;;    yy = ages1[indx1].phot_flamk
;;    djs_plot, xx, yy, position=pos1[*,0], $
;;      xsty=1, ysty=1, xrange=mrange1, yrange=mrange1, psym=4, $
;;      xtitle='', ytitle='K_{s} (FLAMEX, Vega mag)', xtickname=replicate(' ',10)
;;    djs_oplot, !x.crange, !y.crange, line=0, thick=3
;;    im_legend, '<\Delta'+'K_{s}> = '+im_string_stats(yy-xx), $
;;      /left, /top, box=0
;;    djs_plot, xx, yy-xx, $
;;      position=pos1[*,1], /noerase, xsty=1, ysty=1, xrange=mrange1, $
;;      yrange=rrange, psym=4, xtitle='K_{s} (2MASS, Vega mag)', $
;;      ytitle='\Delta'+'K_{s} (Vega mag)'
;;    djs_oplot, !x.crange, [0,0], line=0, thick=3
;;; J-FLAMEX vs 2MASS
;;    xx = ages1[indx1].j_m_ext
;;    yy = ages1[indx1].phot_flamj
;;    djs_plot, xx, yy, position=pos1[*,0], $
;;      xsty=1, ysty=1, xrange=mrange2, yrange=mrange2, psym=4, $
;;      xtitle='', ytitle='J (FLAMEX, Vega mag)', xtickname=replicate(' ',10)
;;    djs_oplot, !x.crange, !y.crange, line=0, thick=3
;;    im_legend, '<\Delta'+'J> = '+im_string_stats(yy-xx), $
;;      /left, /top, box=0
;;    djs_plot, xx, yy-xx, $
;;      position=pos1[*,1], /noerase, xsty=1, ysty=1, xrange=mrange2, $
;;      yrange=rrange, psym=4, xtitle='J (2MASS, Vega mag)', $
;;      ytitle='\Delta'+'J (Vega mag)'
;;    djs_oplot, !x.crange, [0,0], line=0, thick=3
;;; FLAMEX vs NDWFS
;;    xx = ages1[indx1].phot_flamk
;;    yy = ages1[indx1].phot_k-kmks
;;    djs_plot, xx, yy, position=pos1[*,0], $
;;      xsty=1, ysty=1, xrange=mrange1, yrange=mrange1, psym=4, $
;;      xtitle='', ytitle='K-'+strkmks+' (NDWFS, Vega mag)', xtickname=replicate(' ',10)
;;    djs_oplot, !x.crange, !y.crange, line=0, thick=3
;;    im_legend, '<\Delta'+'K_{s}> = '+im_string_stats(yy-xx), $
;;      /left, /top, box=0
;;    im_legend, '<K-K_{s}>_{Vega}=+'+strkmks+' for galaxies', /right, /bottom, box=0, margin=0
;;    djs_plot, xx, yy-xx, $
;;      position=pos1[*,1], /noerase, xsty=1, ysty=1, xrange=mrange1, $
;;      yrange=rrange, psym=4, xtitle='K_{s} (FLAMEX, Vega mag)', $
;;      ytitle='\Delta'+'K_{s} (Vega mag)'
;;    djs_oplot, !x.crange, [0,0], line=0, thick=3
;;
;;    im_plotconfig, psfile=psfile, /psclose, /gzip
;;
