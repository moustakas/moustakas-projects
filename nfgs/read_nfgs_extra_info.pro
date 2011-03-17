function read_nfgs_extra_info
; jm05jul24uofa - read the NFGS_EXTRA_INFO.TXT file, store everything
;                 in a structure, and return

    analysis_path = nfgs_path(/analysis)

    readcol, analysis_path+'nfgs_extra_info.txt', galaxy, $
      nicegalaxy, nedgalaxy, ledagalaxy, morph_type, morph_t, $
      int_agnflag, format='A,A,A,A,A,F,L,L', $
      delimiter='|', comment='#', /silent
    ngalaxy = n_elements(galaxy)

    info = {$
      galaxy:      '', $
      nicegalaxy:  '', $
      nedgalaxy:   '', $
      ledagalaxy:  '', $
      morph_type:  '', $
      morph_t:    0.0, $
      int_agnflag: 0L}
    info = replicate(info,ngalaxy)

    info.galaxy      = galaxy
    info.nicegalaxy  = nicegalaxy
    info.nedgalaxy   = nedgalaxy
    info.ledagalaxy  = ledagalaxy
    info.morph_type  = morph_type
    info.morph_t     = morph_t
    info.int_agnflag = int_agnflag

return, info
end
