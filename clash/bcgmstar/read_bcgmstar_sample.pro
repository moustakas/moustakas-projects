function read_bcgmstar_sample, zsort=zsort
; jm13oct19siena - read the sample
;   sample = rsex(getenv('CLASH_DIR')+'/clash_sample.sex')
    sample = rsex(bcgmstar_path(/propath)+'bcgmstar_sample.sex')

    splog, 'HACK!!!!!'
    sample = sample[0]
    if keyword_set(zsort) then sample = sample[sort(sample.z)]
    
;   struct_print, sample

return, sample
end
