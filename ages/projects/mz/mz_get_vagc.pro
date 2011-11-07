function mz_get_vagc, sample=sample, letter=letter, $
  poststr=poststr, zminmax=zminmax
; jm10sep10ucsd - define the VAGC sample

; POSTSTR/35: -17>Mr>-24; 0.033<z<0.25
    zminmax = [0.05,0.2]
;   zminmax = [0.033,0.25]
    sample = 'dr72'
    letter = 'bsafe'
    poststr = '35'
    vagc_sample = sample+letter+poststr

return, vagc_sample
end
    
