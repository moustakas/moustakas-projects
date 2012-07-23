pro get_fabrice_spec1d
; jm12mar18ucsd - send the AGES/MZ 1d spectra to Fabrice Lamareille

    pp = read_mz_sample(/parent)
    ss = read_ages_gandalf_specfit(pp)

    im_mwrfits, ss, '~/moustakas_ages_spectra.fits', /clobber
    

return
end
    
