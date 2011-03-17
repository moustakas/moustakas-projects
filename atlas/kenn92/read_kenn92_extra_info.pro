function read_kenn92_extra_info
; jm05jul24uofa - read the KENN92_EXTRA_INFO.TXT file, store everything
;                 in a structure, and return

    analysis_path = kenn92_path(/analysis)

    readcol, analysis_path+'kenn92_extra_info.txt', galaxy, $
      nicegalaxy, altgalaxy, morph_type, agnflag, lowresflag, $
      comments, format='A,A,A,A,L,L,A', delimiter='|', $
      comment='#', /silent
    ngalaxy = n_elements(galaxy)

    info = {$
      galaxy:      '', $
      nicegalaxy:  '', $
      altgalaxy:   '', $
      morph_type:  '', $
      agnflag:     0L, $
      lowresflag:  0L}
    info = replicate(info,ngalaxy)

    info.galaxy      = galaxy
    info.nicegalaxy  = nicegalaxy
    info.altgalaxy   = altgalaxy
    info.morph_type  = morph_type
    info.agnflag     = agnflag
    info.lowresflag  = lowresflag

return, info
end
