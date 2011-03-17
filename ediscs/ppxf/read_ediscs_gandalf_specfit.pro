function read_ediscs_gandalf_specfit, specdata, solar=solar, $
  observed=observed, linear=linear, silent=silent
; jm10apr28ucsd - based on READ_EDISCS_SPECFIT; default is to return
;   *rest-frame* ln-binned spectra; use /OBSERVED to convert to the
;   observed frame, and /LINEAR to convert to linear wavelength

    nobj = n_elements(specdata)
    if (nobj eq 0) then begin
       doc_library, 'read_ediscs_gandalf_specfit'
       return, -1
    endif

    version = ediscs_version(/ppxf_specfit)    
    specfitpath = ediscs_path(/ppxf)

    if keyword_set(solar) then suffix = 'solar_' else suffix = ''
    specfitfile = specfitpath+'ediscs_specfit_'+suffix+version+'.fits.gz'
    
; read the SPECFIT file       
    if (file_test(specfitfile,/reg) eq 0) then begin
       splog, 'SPECFITFILE '+specfitfile+' not found'
       return, -1
    endif
    specfit1 = mrdfits(specfitfile,1,/silent)

    dwave = replicate(ediscs_pixel_size(),nobj) ; [A/pixel] (see /linear)
    run34 = where(strtrim(specdata.run,2) eq 'run34',nrun34)
    if (nrun34 ne 0L) then dwave[run34] = ediscs_pixel_size(/run34)

; pull out the relevant spectra
    match, strtrim(specdata.ediscs_id,2), strtrim(specfit1.ediscs_id,2), m1, m2
    if (n_elements(m1) ne nobj) then message, 'One or more spectra not found!'
    srt = sort(m1) & m1 = m1[srt] & m2 = m2[srt]
;   niceprint, specfit1[m2].galaxy, specdata.galaxy
    specfit = specfit1[m2]
    
; convert to linear wavelength and/or the observed frame, if desired
    if keyword_set(observed) then begin
       splog, 'Shifting to the observed frame'
       for ii = 0L, nobj-1L do begin
          zfactor = (1.0+specfit[ii].zabs)
          specfit[ii].wave = specfit[ii].wave + alog(zfactor)*(specfit[ii].wave gt 0.0)
          specfit[ii].flux = specfit[ii].flux/zfactor
          specfit[ii].ferr = specfit[ii].ferr/zfactor
          specfit[ii].linefit = specfit[ii].linefit/zfactor
          specfit[ii].continuum = specfit[ii].continuum/zfactor
          specfit[ii].smooth_continuum = specfit[ii].smooth_continuum/zfactor
       endfor       
    endif

    if keyword_set(linear) then begin
       splog, 'Rebinning linearly in wavelength'
       ln_specfit = specfit
       specfit.wave = 0.0
       specfit.flux = 0.0
       specfit.ferr = 0.0
       specfit.linefit = 0.0
       specfit.continuum = 0.0
       specfit.smooth_continuum = 0.0
       for ii = 0L, nobj-1L do begin
          good = where(ln_specfit[ii].wave gt 0.0)
          wave = exp(ln_specfit[ii].wave[good])*1.0D
          flux = ln_specfit[ii].flux[good]*1.0D
          ferr = ln_specfit[ii].ferr[good]*1.0D
          linefit = ln_specfit[ii].linefit[good]*1.0D
          continuum = ln_specfit[ii].continuum[good]*1.0D
          smooth_continuum = ln_specfit[ii].smooth_continuum[good]*1.0D

; call XREBIN() directly to avoid additional overhead; note that
; OUTVAR can *sometimes* have a handful of negative pixels, but
; we don't care here
          outflux = im_linear_rebin(wave,flux,dwave=dwave[ii],$
            minwave=min(wave),maxwave=max(wave),outwave=outwave,$
            borders=borders,newborders=newborders)
          outvar = xrebin(borders,ferr^2.0,newborders,/splinf)
;         zero = where(outvar le 0.0)
;         if (zero[0] ne -1) then message, 'Fix me!'
          outferr = sqrt(abs(outvar))

          outlinefit = xrebin(borders,linefit,newborders,/splinf)
          outcontinuum = xrebin(borders,continuum,newborders,/splinf)
          outsmooth_continuum = xrebin(borders,smooth_continuum,newborders,/splinf)

          npix = n_elements(outflux)
          specfit[ii].wave[0:npix-1] = outwave
          specfit[ii].flux[0:npix-1] = outflux
          specfit[ii].ferr[0:npix-1] = outferr
          specfit[ii].linefit[0:npix-1] = outlinefit
          specfit[ii].continuum[0:npix-1] = outcontinuum
          specfit[ii].smooth_continuum[0:npix-1] = outsmooth_continuum
       endfor
    endif

return, specfit
end
