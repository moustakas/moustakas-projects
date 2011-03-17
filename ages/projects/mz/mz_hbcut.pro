function mz_hbcut, sdss=sdss
; H-beta flux limit for the MZ sample
    if keyword_set(sdss) then hbcut = 1D-16 else $
      hbcut = 3D-17
return, hbcut
end
    
