function read_sings_extra_info
; jm05jul25uofa - read the SINGS_EXTRA_INFO.TXT file, store everything
;                 in a structure, and return

    analysis_path = sings_path(/analysis)

    readcol, analysis_path+'sings_extra_info.txt', galaxy, $
      nicegalaxy, altgalaxy, drift20_comments, drift56_comments, $
      nuclear_comments, drift20_agnflag, drift56_agnflag, $
      nuclear_agnflag, format='A,A,A,A,A,A,L,L,L', $
      delimiter='|', comment='#', /silent
    ngalaxy = n_elements(galaxy)

    moreinfo = rsex(analysis_path+'sings_ancillary_data.sex')
    
    info = {$
      galaxy:           '', $
      nicegalaxy:       '', $
      altgalaxy:        '', $
      drift20_comments: '', $
      drift56_comments: '', $
      nuclear_comments: '', $
      drift20_agnflag:  0L, $
      drift56_agnflag:  0L, $
      nuclear_agnflag:  0L}
    info = replicate(info,ngalaxy)

    info.galaxy           = galaxy
    info.nicegalaxy       = nicegalaxy
    info.altgalaxy        = altgalaxy
    info.drift20_agnflag  = drift20_agnflag
    info.drift56_agnflag  = drift56_agnflag
    info.nuclear_agnflag  = nuclear_agnflag
    info.drift20_comments = drift20_comments
    info.drift56_comments = drift56_comments
    info.nuclear_comments = nuclear_comments

    info = struct_addtags(info,struct_trimtags(moreinfo,$
      select=['TYPE','CLASS','DISTANCE','DMAJ','DMIN']))
    
return, info
end
