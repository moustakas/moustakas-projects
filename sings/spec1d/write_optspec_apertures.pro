pro write_optspec_apertures
; jm08jan15nyu - write out the spectroscopic apertures so that D. Dale
;                can give me optical light fractions
    
    s = sings_read_info()
    ss = struct_trimtags(s,select=['GALAXY','RA','DEC','DRIFT56','drift20','nuclear','DRIFT56_SCAN',$
      'DRIFT56_AP','DRIFT56_POSANGLE','drift20_scan','drift20_ap','drift20_posangle','nuclear_scan',$
      'nuclear_ap','nuclear_posangle'])

; I need to go back and throw away some of the dwarf spectra as being
; crappy; plus I'm not really sure where we were pointing, so
; the light fractions won't be reliable; also the radial strip
; spectrum of NGC3034 is pretty useless, so toss it


    


return
end
