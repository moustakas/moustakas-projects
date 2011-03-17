pro batch_parse_ages_gandalf_specfit, firstpass, nprocess
; jm09nov18ucsd - batch wrapper script on PARSE_AGES_GANDALF_SPECFIT

    if (n_elements(nprocess) eq 0) then nprocess = 8
    
; assign a unique process number to each pass    
    allpass1 = ages_allpasses(/fluxed)
    nallpass = n_elements(allpass1)
;   index = indgen(n_elements(allpass1))+1
;   prockeep = where(((index mod nprocess) eq process),npass)
;   prockeep = process*indgen(
    if ((nallpass mod nprocess) ne 0) then message, 'This might not work'
    npass = nallpass/nprocess
    dopass = npass*(firstpass-1)+indgen(npass)
    allpass = allpass1[dopass]
    print, allpass
    
; now do each pass
;   for ipass = 5, 5 do begin
    for ipass = 0, npass-1 do begin
       splog, 'Pass '+allpass[ipass]
       parse_ages_gandalf_specfit, allpass[ipass]
    endfor

return
end
