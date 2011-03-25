pro rebuild_sings_models
; jm11mar21ucsd - build the model fits to the SINGS galaxies
; for M. Brown's efforts to build new photoz templates

;   nuc = read_sings(/nuclear)
;   d20 = read_sings(/drift20)

    d56 = read_sings_gandalf(/drift56,/solar)
    d56fit = read_sings_gandalf_specfit(d56,/drift56,/solar,/linear)
    sd56 = restore_sings_ppxf_bestfit(d56.continuum_coeff,/solar,$
      ebv=d56.continuum_ebv,bestwave=modelwave,bestage=modelage)

    ngal = n_elements(d56)
    npix = n_elements(modelwave)
    
    out = {galaxy: d56.galaxy, ra: dblarr(ngal), dec: dblarr(ngal), $
      zabs: d56.zabs, snr: d56.continuum_snr, wave: modelwave, $
      flux: fltarr(npix,ngal), isdata: intarr(npix,ngal)}
    out.ra = 15D*im_hms2dec(d56.ra)
    out.dec = im_hms2dec(d56.dec)
    for ii = 0, ngal-1 do begin
       good = where(d56fit[ii].wave gt 0.0)
       linterp, d56fit[ii].wave[good], d56fit[ii].flux[good], modelwave, flux1, missing=0
       out.flux[*,ii] = sd56[*,ii]*(flux1 eq 0) + flux1
       out.isdata[*,ii] = flux1 ne 0
    endfor

; write out
    outfile = sings_path(/proj)+'photoztemplates/sings_drift56_models.fits'
    im_mwrfits, out, outfile, /clobber

return
end
    
