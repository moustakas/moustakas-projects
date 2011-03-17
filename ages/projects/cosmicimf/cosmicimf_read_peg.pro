function cosmicimf_read_peg, sfhfile, age=age, sfr_const=sfr_const
; wrapper to read the Pegase models; scale the constant SFR models

    peg = im_read_peg(sfhfile)
    if (n_elements(age) ne 0) then begin
       get_element, peg.age, age, indx ; [Myr]
       peg = peg[indx]
    endif

; scale the constant SFR models back to 1 M_sun/yr
    if strmatch(sfhfile,'*const*',/fold) then begin
       factor = 1.0/(sfr_const*1D-6)
       peg.flux = peg.flux*factor
       peg.lineflux = peg.lineflux*factor
       peg.sfr = peg.sfr*factor
; ...and probably more       
    endif

return, peg
end

