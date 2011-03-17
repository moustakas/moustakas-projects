pro nfgs_write_nedlist, write=write
; jm04jan21uofa
; jm05may17uofa - also write out a list for LEDA
; jm05jul24uofa - updated

    nedpath = nfgs_path(/ned)

    info = read_nfgs_extra_info()    
    ngalaxy = n_elements(info)
    splog, 'There are '+string(ngalaxy,format='(I0)')+' unique galaxies in the NFGS.'

    if keyword_set(write) then begin
       
; write the list of unique galaxy names to be sent to NED and LEDA 

       openw, lun1, nedpath+'nfgs_ned_input.txt', /get_lun
       openw, lun2, nedpath+'nfgs_leda_input.txt', /get_lun

       for j = 0L, ngalaxy-1L do begin
          printf, lun1, info[j].nedgalaxy, format='(A0)'
          printf, lun2, info[j].ledagalaxy, format='(A0)'
       endfor

       free_lun, lun1, lun2

    endif
       
return
end
