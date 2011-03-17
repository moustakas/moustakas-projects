function ages_restore_speclinefit, ispec, templateinfo=templateinfo
; jm08aug29nyu - restore the SPECLINEFIT results; by default the
;   models are returned at the BC03 wavelength spacing; need to write
;   a wrapper to restore at the wavelength spacing of the AGES spectra
;   themselves 
    
    if (n_elements(ispec) eq 0L) then ispec = read_ages(/ispec)
    model = irestore_speclinefit(ispec,templateinfo=templateinfo,$
      templatedir=ages_path(/specfit),/repair_4959)

return, model
end
