function read_ages_gandalf_specfit, specdata1, unfluxed=unfluxed, $
  observed=observed, linear=linear, solar=solar, silent=silent
; jm09nov30ucsd - based on READ_AGES_SPECFIT; default is to return
;   *rest-frame* ln-binned spectra; use /OBSERVED to convert to the
;   observed frame, and /LINEAR to convert to linear wavelength

    nobj = n_elements(specdata1)
    if (nobj eq 0) then begin
       print, 'specfit = read_ages_gandalf_specfit(galaxy,_extra=extra)'
       return, -1L
    endif

    version = ages_version(/ppxf_specfit)    
    specfitpath = ages_path(/spec1d)+'fluxed/tweak/'+version+'/' 

    pass = fix(specdata1.pass)
    upass = pass[uniq(pass,sort(pass))]
    npass = n_elements(upass)

    dwave = ages_pixel_size() ; [A/pixel] (see /linear)
    
    for ipass = 0, npass-1 do begin
       print, format='("Reading pass ",I0,"/",I0,A10,$)', $
         ipass+1, npass, string(13b)

; read the SPECFIT file       
       if keyword_set(solar) then suffix = '_solar' else suffix = ''
       specfitfile = specfitpath+'specfit_'+$
         string(upass[ipass],format='(I3.3)')+suffix+'.fits.gz'
       if (file_test(specfitfile,/reg) eq 0) then begin
          splog, 'SPECFITFILE '+specfitfile+' not found'
          continue
       endif
       specfit1 = mrdfits(specfitfile,1,/silent)

; pull out the relevant spectra
       inindx = where(upass[ipass] eq pass,nin)
       specdata = specdata1[inindx]

       match, specdata.ages_id, specfit1.ages_id, m1, m2
       if (n_elements(m1) ne nin) then message, 'One or more spectra not found!'
       srt = sort(m1) & m1 = m1[srt] & m2 = m2[srt]
;      niceprint, specfit1[m2].ages_id, specdata.ages_id
       specfit = specfit1[m2]

       if (ipass eq 0) then out_specfit = specfit else $
         out_specfit = [temporary(out_specfit),specfit]
    endfor

; convert to linear wavelength and/or the observed frame, if desired
    if keyword_set(observed) then begin
       splog, 'Shifting to the observed frame'
       for ii = 0L, nobj-1L do begin
          zfactor = (1.0+out_specfit[ii].zabs)
          out_specfit[ii].wave = out_specfit[ii].wave + alog(zfactor)*(out_specfit[ii].wave gt 0.0)
          out_specfit[ii].flux = out_specfit[ii].flux/zfactor
          out_specfit[ii].ferr = out_specfit[ii].ferr/zfactor
          out_specfit[ii].linefit = out_specfit[ii].linefit/zfactor
          out_specfit[ii].continuum = out_specfit[ii].continuum/zfactor
          out_specfit[ii].smooth_continuum = out_specfit[ii].smooth_continuum/zfactor
       endfor       
    endif

    if keyword_set(linear) then begin
       splog, 'Rebinning linearly in wavelength'
       ln_specfit = out_specfit
       out_specfit.wave = 0.0
       out_specfit.flux = 0.0
       out_specfit.ferr = 0.0
       out_specfit.linefit = 0.0
       out_specfit.continuum = 0.0
       out_specfit.smooth_continuum = 0.0
       for ii = 0L, nobj-1L do begin
          good = where(ln_specfit[ii].wave gt 0.0)
          wave = exp(ln_specfit[ii].wave[good])*1.0D
          flux = ln_specfit[ii].flux[good]*1.0D
          ferr = ln_specfit[ii].ferr[good]*1.0D
          linefit = ln_specfit[ii].linefit[good]*1.0D
          continuum = ln_specfit[ii].continuum[good]*1.0D
          smooth_continuum = ln_specfit[ii].smooth_continuum[good]*1.0D

; x_specrebin is very slow compared to im_linear_rebin()
;         outwave = im_array(min(wave),max(wave),1.2)
;         x_specrebin, wave, flux, outwave, outflux, $
;           var=ferr^2, nwvar=outvar, /silent
;         zero = where(outvar le 0.0)
;         if (zero[0] ne -1) then message, 'Fix me!'
;         outferr = sqrt(outvar)
;         x_specrebin, wave, linefit, outwave, outlinefit, /silent
;         x_specrebin, wave, continuum, outwave, outcontinuum, /silent
;         x_specrebin, wave, smooth_continuum, outwave, outsmooth_continuum, /silent

; call XREBIN() directly to avoid additional overhead; note that
; OUTVAR can *sometimes* have a handful of negative pixels, but
; we don't care here
          outflux = im_linear_rebin(wave,flux,dwave=dwave,$
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
          out_specfit[ii].wave[0:npix-1] = outwave
          out_specfit[ii].flux[0:npix-1] = outflux
          out_specfit[ii].ferr[0:npix-1] = outferr
          out_specfit[ii].linefit[0:npix-1] = outlinefit
          out_specfit[ii].continuum[0:npix-1] = outcontinuum
          out_specfit[ii].smooth_continuum[0:npix-1] = outsmooth_continuum
       endfor
    endif

; resort the output
    match, specdata1.ages_id, out_specfit.ages_id, m1, m2
    srt = sort(m1)
    out_specfit = out_specfit[m2[srt]]
    
return, out_specfit
end
