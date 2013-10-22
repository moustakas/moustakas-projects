function read_bcgsfhs_sample, noa2261=noa2261
; jm13oct19siena - read the sample
;   sample = rsex(getenv('CLASH_DIR')+'/clash_sample.sex')
    sample = rsex(bcgsfhs_path(/propath)+'bcgsfhs_sample.sex')

    if keyword_set(noa2261) then begin
       splog, 'IGNORING A2261!!!'
       keep = where(strtrim(sample.shortname,2) ne 'a2261')
       sample = sample[keep]
    endif
    struct_print, sample

return, sample
end
