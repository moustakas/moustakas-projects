pro atlas_write_nedlist, write=write
; jm02mar22uofa
; jm05jul21uofa - rewritten to accomodate the new "extra info" text
;                 file, which contains NED- and LEDA-compatible names 

; also generate the LEDA input list
    
    nedpath = atlas_path(/ned)

    info = read_atlas_extra_info()    
    ngalaxy = n_elements(info)
    splog, 'There are '+string(ngalaxy,format='(I0)')+' unique galaxies in the ATLAS.'

    if keyword_set(write) then begin
       
; write the list of unique galaxy names to be sent to NED and LEDA 

       openw, lun1, nedpath+'atlas_ned_input.txt', /get_lun
       openw, lun2, nedpath+'atlas_leda_input.txt', /get_lun

       for j = 0L, ngalaxy-1L do begin
          printf, lun1, info[j].nedgalaxy, format='(A0)'
          printf, lun2, info[j].ledagalaxy, format='(A0)'
       endfor

       free_lun, lun1, lun2

    endif
       
return
end
