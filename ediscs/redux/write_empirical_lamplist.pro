pro write_empirical_lamplist
; jm04jun24uofa

    readcol, 'initial_arc_solution.dat', wave, intensity, pixel, $
      format='F,A', comment='#', /silent
    nline = n_elements(wave)
    
; write out

    openw, lun, 'GRIS_600RI_lamplist.dat', /get_lun
    printf, lun, '# Line list for the GRIS_600RI grism.'
    for i = 0L, nline-1L do printf, lun, wave[i], intensity[i], 'GOOD', $
      'El', format='(F7.2,3x,F12.2,3x,A4,3x,A7)'
    free_lun, lun
    
return
end
    
