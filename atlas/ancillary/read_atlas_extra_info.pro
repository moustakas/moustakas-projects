function read_atlas_extra_info
; jm05jul21uofa - read the READ_ATLAS_EXTRA_INFO.TXT file, store
;                 everything in a structure, and return

    analysis_path = atlas_path(/analysis)

    readcol, analysis_path+'atlas_extra_info.txt', galaxy, $
      nicegalaxy, altgalaxy, nedgalaxy, ledagalaxy, morph_type, $
      morph_t, redshift, dmaj, dmin, drift_code, int_agnflag, nuc_agnflag, $
      comments, format='A,A,A,A,A,A,F,F,F,F,L,L,L,A', $
      delimiter='|', comment='#', /silent
    ngalaxy = n_elements(galaxy)

    info = {$
      galaxy:      '', $
      nicegalaxy:  '', $
      altgalaxy:   '', $
      nedgalaxy:   '', $
      ledagalaxy:  '', $
      morph_type:  '', $
      morph_t:    0.0, $
      redshift:   0.0, $
      dmaj:       0.0, $
      dmin:       0.0, $
      drift_code:  0L, $
      int_agnflag: 0L, $
      nuc_agnflag: 0L, $
      comments:    ''}
    info = replicate(info,ngalaxy)

    info.galaxy      = galaxy
    info.nicegalaxy  = nicegalaxy
    info.altgalaxy   = altgalaxy
    info.nedgalaxy   = nedgalaxy
    info.ledagalaxy  = ledagalaxy
    info.morph_type  = morph_type
    info.morph_t     = morph_t
    info.redshift    = redshift
    info.dmaj        = dmaj
    info.dmin        = dmin
    info.drift_code  = drift_code 
    info.int_agnflag = int_agnflag
    info.nuc_agnflag = nuc_agnflag
    info.comments    = comments

return, info
end
