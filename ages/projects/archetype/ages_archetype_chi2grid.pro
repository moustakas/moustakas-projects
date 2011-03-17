function archetype_sample
; define the sample for the archetype project

    ages1 = read_ages(/ispec)

;    good = where((ages1.pass ne 106) and (ages1.pass ne 110) and (ages1.pass ne 209) and $
;      (ages1.pass ne 310) and (ages1.pass ne 311) and (ages1.main_flag eq 1) and $
;      (ages1.z gt 0.01) and (ages1.z lt 0.4),ngood)
;;     (ages1.z gt 0.21) and (ages1.z lt 0.212),ngood)
;    ages = ages1[good]

    badplates = [104, 105, 106, 110, 113, 115, $
      201, 203, 205, 206, 207, 208, 209, 212, 213, $
      214, 215, 301, 302, 303, 304, 306, 307, 308, $
      309, 310, 311, 312, 313, 314, 315]

    good = lindgen(n_elements(ages1))
    good[where_array(badplates,ages1.pass)] = -1L
    good = good[where(good ne -1L)]
    ages = ages1[good]
    if arg_present(ancdata) then ancdata = ancdata1[good]

    zmin = 0.05 & zmax = 0.35
    snrcut_cont = 10.0
    snrcut_line = 5.0
    ewerrcut = 3.0
    chi2cut = 10.0
    
    these = where(ages.class eq 'GALAXY' and $
      (ages.z ge zmin) and (ages.z le zmax) and $
      (ages.continuum_snr ge snrcut_cont) and $
      ((ages.oii_3727_chi2 lt chi2cut) and ((ages.oii_3727_ew[1] lt ewerrcut) or $
      (ages.oii_3727_ew[0]/ages.oii_3727_ew[1] gt snrcut_line))) and $
      ((ages.h_beta_chi2 lt chi2cut) and ((ages.h_beta_ew[1] lt ewerrcut) or $
      (ages.h_beta_ew[0]/ages.h_beta_ew[1] gt snrcut_line))) and $
      ((ages.h_alpha_chi2 lt chi2cut) and ((ages.h_alpha_ew[1] lt ewerrcut) or $
      (ages.h_alpha_ew[0]/ages.h_alpha_ew[1] gt snrcut_line))),nthese)

;   ages_display_spectrum, ages[these], /postscript, $
;     psname='primus_ages_quality.ps'

; v2.0 crap spectra
    crap = 'ages_'+['101/121','101/164','101/257','108/046','404/081','406/086',$
      '406/092','417/040','420/120','609/131','611/097','713/116','722/194']
    good = lindgen(nthese)
    good[where_array(crap,ages[these].galaxy)] = -1L
    good = good[where(good ne -1L)]
    ages = ages[these[good]]

;   ages = ages[0:20]
    
return, ages
end

pro ages_archetype_chi2grid, ss, out
; jm08sep19nyu - 

; read the quality sample and the observed spectra 

    ages = archetype_sample()
    ngalaxy = n_elements(ages)
    if (n_elements(ss) eq 0L) then ss = read_ages_specfit(ages.galaxy,/obswave)

; define a common wavelength array and interpolate
    
    allwave = reform(ss[*,0,*])
    allflux = reform(ss[*,1,*])
    allmodels = reform(ss[*,2,*]+ss[*,3,*])
    allivar = reform(ss[*,5,*])

    minwave = 0.0 & maxwave = 1E4
    for ii = 0L, ngalaxy-1L do begin
       good = where((allwave[*,ii] gt 0.0) and (allivar[*,ii] gt 0.0))
       minwave = min(allwave[good,ii])>minwave
       maxwave = max(allwave[good,ii])<maxwave
    endfor

    good = where((allwave gt 0.0) and (allivar gt 0.0))
    dwave = median(allwave[good]-shift(allwave[good],1))
    wave = findgen((maxwave-minwave)/dwave+1)*dwave+minwave
    npix = n_elements(wave)

    flux = fltarr(npix,ngalaxy)
    ivar = flux*0.0
    models = flux*0.0
    
    for ii = 0L, ngalaxy-1L do begin
       good = where((allwave[*,ii] gt 0.0) and (allivar[*,ii] gt 0.0))
       linterp, allwave[good,ii], allflux[good,ii], wave, flux1
       linterp, allwave[good,ii], allivar[good,ii], wave, ivar1
       linterp, allwave[good,ii], allmodels[good,ii], wave, models1
       flux[*,ii] = flux1 & ivar[*,ii] = ivar1 & models[*,ii] = models1
    endfor

; now construct the chi2 grid, and then subtract out the chi2 value
; along the diagonal

    out = {ages_id: 0L, chi2: fltarr(ngalaxy), scale: fltarr(ngalaxy)}
    out = replicate(out,ngalaxy)
    out.ages_id = ages.ages_id

    for ii = 0L, ngalaxy-1L do begin
; build the mask
       mask = iemission_mask(wave,z=ages[ii].z,/telluric,/sky);,/nebular,/qso)
       good = where((ivar[*,ii]*mask gt 0.0) and (finite(ivar[*,ii]) eq 1),ngood)
       rejmask = mask*0.0+1.0
       djs_iterstat, flux[good,ii]-models[good,ii], sigrej=3.0, mask=rejmask1
       rejmask[good] = rejmask1
       mask = (mask and rejmask) and (wave gt ages[ii].continuum_minwave) and $
         (wave lt ages[ii].continuum_maxwave)
       djs_plot, wave, flux[*,ii], ps=10, xsty=3, ysty=3
       w1 = where(mask eq 0) & djs_oplot, wave[w1], flux[w1,ii], ps=10, color='red'
       w2 = where(rejmask eq 0) & djs_oplot, wave[w2], flux[w2,ii], ps=10, color='green'
stop
       
       for jj = 0L, ngalaxy-1L do begin
          out[ii].scale[jj] = total(mask*ivar[*,ii]*flux[*,ii]*models[*,jj])/$
            total(ivar[*,ii]*models[*,jj]^2.0)
          out[ii].chi2[jj] = total(mask*ivar[*,ii]*(flux[*,ii]-out[ii].scale[jj]*models[*,jj])^2.0)
;         djs_plot, wave, flux[*,ii], ps=10, xsty=3, ysty=3
;         djs_oplot, wave, out[ii].scale[ii]*models[*,ii], ps=10, color='red'
;         djs_oplot, wave, out[ii].scale[jj]*models[*,jj], ps=10, color='green'
       endfor

       out[ii].chi2 = out[ii].chi2 - out[ii].chi2[ii]
;      dum = out[ii].chi2 - out[ii].chi2[ii]
;      neg = where(dum lt 0.0) & if neg[0] ne -1 then stop
;      print, ages[ii].continuum_chi2, out[ii].chi2[ii], ages[neg].continuum_chi2, out[ii].chi2[neg]
;      djs_plot, wave, models[*,ii], ps=10, xsty=3, ysty=3
;      djs_oplot, wave, (out[ii].scale[neg])[0]*models[*,neg], ps=10, color='red'
;      ages_display_spectrum, ages[[ii,neg]], plottype=2

    endfor

    mwrfits, out, 'ages_archetype_chi2grid.fits', /create
    
return
end
    
