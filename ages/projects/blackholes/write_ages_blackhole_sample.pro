pro write_ages_blackhole_sample, allagesdust1, allagesancillary1, write=write
; jm06oct25nyu - written, based on WRITE_AGES_MZ_SAMPLE; find objects
;   with broad lines, from which we can derive the black-hole mass

    outpath = ages_path(/projects)+'blackholes/'

    if (n_elements(allagesdust1) eq 0L) then allagesdust1 = read_ages(ancillary=allagesancillary1)
    ngalaxy = n_elements(allagesdust1)

    allagesdust = allagesdust1
    allagesancillary = allagesancillary1

    snrcut = 4.0 ; determine this using simulations
    fwhm2sig = 2.35
    broad_sigcut = 500.0
    
; broad Mg II *and* H-beta
    
    broad = where((allagesdust.mgii_2800[0]/allagesdust.mgii_2800[1] gt snrcut) and $
      (allagesdust.mgii_2800[0] gt 0.0) and (allagesdust.mgii_2800_sigma[0] ge broad_sigcut/fwhm2sig) and $
      (allagesdust.h_beta[0]/allagesdust.h_beta[1] gt snrcut) and $
      (allagesdust.h_beta_ew[0] gt 0.0) and (allagesdust.h_beta_sigma[0] ge broad_sigcut/fwhm2sig),nbroad)

    splog, 'Black hole sample: '+string(nbroad,format='(I0)')+'/'+string(ngalaxy,format='(I0)')+' ('+$
      string(round(100.0*nbroad/ngalaxy),format='(I0)')+'%).'
    
    zstats = im_stats(allagesancillary[broad].z)
    splog, '   Redshift      : ['+strtrim(string(zstats.min,format='(F12.2)'),2)+'-'+strtrim(string(zstats.max,format='(F12.2)'),2)+'] '+$
      strtrim(string(zstats.median,format='(F12.2)'),2)+' ('+strtrim(string(zstats.mean,format='(F12.2)'),2)+'+/-'+$
      strtrim(string(zstats.sigma,format='(F12.2)'),2)+')'
    good = where(allagesancillary[broad].infiber_r gt -900.0)

    good = where(allagesancillary[broad].infiber_r gt -900.0)
    fstats = im_stats(allagesancillary[broad[good]].infiber_r)
    splog, '   Light fraction: ['+strtrim(string(fstats.min,format='(F12.2)'),2)+'-'+strtrim(string(fstats.max,format='(F12.2)'),2)+'] '+$
      strtrim(string(fstats.median,format='(F12.2)'),2)+' ('+strtrim(string(fstats.mean,format='(F12.2)'),2)+'+/-'+$
      strtrim(string(fstats.sigma,format='(F12.2)'),2)+')'

    nirphot = where((allagesancillary[broad].phot_k gt -900.0) or ((allagesancillary[broad].phot_flamj gt -900.0) and $
      (allagesancillary[broad].phot_flamk gt -900.0)),nnirphot)
    splog, '   NIR photometry: '+string(nnirphot,format='(I0)')+'/'+string(nbroad,format='(I0)')+' ('+$
      string(round(100.0*nnirphot/nbroad),format='(I0)')+'%).'
    print

;   ages_display_spectrum, allagesdust[broad[0:20]], allagesancillary[broad[0:20]], specfit=specbroad

; write out
    
    if keyword_set(write) then begin

       splog, 'Writing '+outpath+'ages_blackhole_speclinefit.fits.gz'
       mwrfits, allagesdust[broad], outpath+'ages_blackhole_speclinefit.fits', /create
       spawn, ['gzip -f '+outpath+'ages_blackhole_speclinefit.fits'], /sh
 
       splog, 'Writing '+outpath+'ages_blackhole_ancillary.fits.gz'
       mwrfits, allagesancillary[broad], outpath+'ages_blackhole_ancillary.fits', /create
       spawn, ['gzip -f '+outpath+'ages_blackhole_ancillary.fits'], /sh

    endif

return
end
