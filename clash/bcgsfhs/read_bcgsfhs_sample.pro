function read_bcgsfhs_sample, zsort=zsort
; jm13oct19siena - read the sample
;   sample = rsex(getenv('CLASH_DIR')+'/clash_sample.sex')
    sample = rsex(bcgsfhs_path(/propath)+'bcgsfhs_sample.sex')
    if keyword_set(zsort) then sample = sample[sort(sample.z)]
    
;   struct_print, sample

return, sample
end
