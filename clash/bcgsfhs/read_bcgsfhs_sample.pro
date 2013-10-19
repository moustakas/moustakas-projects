function read_bcgsfhs_sample, noa2261=noa2261
; jm13oct19siena - read the sample
;   sample = rsex(getenv('CLASH_DIR')+'/clash_sample.sex')
    sample = rsex(bcgsfhs_path(/propath)+'bcgsfhs_sample.sex')
    struct_print, sample

return, sample
end
