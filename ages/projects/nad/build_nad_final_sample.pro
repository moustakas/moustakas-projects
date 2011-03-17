pro build_nad_final_sample
; jm11mar01ucsd

    nadpath = ages_path(/proj)+'nad/'

    ind = mrdfits(nadpath+'nad_indices.fits.gz',1)
    ppxf = mrdfits(nadpath+'nad_ppxf.fits.gz',1)
    phot = mrdfits(nadpath+'nad_phot.fits.gz',1)

    ww = where(ind.lick_nad[0] gt 20,ngal)
    
    specfit = read_ages_gandalf_specfit(ppxf[ww],/solar)
    for ii = 0, ngal-1 do begin

       good = where(specfit[ii].wave gt 0.0,npix)
       wave = exp(specfit[ii].wave[good])
       flux = specfit[ii].flux[good]
       plot, wave, flux, xsty=3, ysty=3, xrange=[5850,5950], psym=10
       cc = get_kbrd(1)
    endfor

    
    
return
end

