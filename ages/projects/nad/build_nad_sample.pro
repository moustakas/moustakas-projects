pro build_nad_sample
; jm11mar01ucsd

    nadpath = ages_path(/proj)+'nad/'

    ppxf = read_ages(/ppxf)
    phot = read_ages(/phot)
    phot = phot[ppxf.ages_id]

    rejplates = [$
      104,$                     ; average fluxing vector
      106,110,209,310,311,$     ; unfluxed
      604]                      ; very crummy 

    nadwave = 5890.0 ; [Angstroms]
    nadwidth = 35.0  ; [Angstroms]
    snrcut = 1.0
    
    parent = where((ppxf.minwave lt nadwave-20.0) and (ppxf.maxwave gt nadwave+20.0) and $
      (phot.phot_mips24 gt -900.0) and (phot.pass ne rejplates[0]) and $
      (phot.pass ne rejplates[1]) and (phot.pass ne rejplates[2]) and $
      (phot.pass ne rejplates[3]) and (phot.pass ne rejplates[4]) and $
      (phot.pass ne rejplates[5]) and (phot.pass ne rejplates[6]),ngal)
;   parent = parent[0:10]
;   ngal = n_elements(parent)
    
    nadspec = ppxf[parent]
    nadphot = phot[parent]

;   qaplot_ages_gandalf_specfit, nadspec[0:10], /solar, $
;     psfile=nadpath+'qa_nadspec.ps'

    specfit = read_ages_gandalf_specfit(nadspec,/solar)
    for ii = 0, ngal-1 do begin
       print, ii, ngal, string(13B), format='(I0,"/",I0,A10,$)'
       good = where(specfit[ii].wave gt 0.0,npix)
       wave = exp(specfit[ii].wave[good])
       flux = specfit[ii].flux[good]
       ferr = specfit[ii].ferr[good]
 
       ind1 = spectral_indices(wave,flux,ivar=(ferr lt 1E5)/ferr^2,$
         indexpath=nadpath,indexfile='nadindex.txt',/silent)
;      nad = im_splotew(wave,flux,ferr,nadwave,boxwidth=nadwidth,$
;        /absorption,doplot=debug,/silent)
       if (ii eq 0) then ind = ind1 else ind = [ind,ind1]
    endfor

    keep = where(ind.lick_nad[0]/ind.lick_nad[1] gt snrcut,nkeep)
    
; write out    
    im_mwrfits, nadspec[keep], nadpath+'nad_ppxf.fits', /clobber
    im_mwrfits, nadphot[keep], nadpath+'nad_phot.fits', /clobber
    im_mwrfits, ind[keep], nadpath+'nad_indices.fits', /clobber
    
stop    
    
return
end
    

    
