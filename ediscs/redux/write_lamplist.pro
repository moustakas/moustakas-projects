pro write_lamplist
; jm04jun24uofa
; use Bo's lamp list to generate an ispec-compatible line list

    readcol, 'arc_lines.dat', bowave, boelement, format='F,A', $
      comment='#', /silent

    lamp = iread_lamplist('HeNeArHgCd')

    get_element, lamp.lambda, bowave, indx
    blamp = lamp[indx]

    bad = where(abs(blamp.lambda-bowave) gt 0.05,comp=good)
    blamp = blamp[good]
    bowave = bowave[good]
    boelement = boelement[good]
    
    niceprint, blamp.lambda, bowave, blamp.element, boelement

; write out

    openw, lun, 'GRIS_600RI_lamplist.dat', /get_lun
    printf, lun, '# Line list for the GRIS_600RI grism.'
    struct_print, blamp, lun=lun, /no_head
    free_lun, lun
    
return
end
    
