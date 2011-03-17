pro write_highzea_wrapper, cleanfits=cleanfits, overwrite=overwrite, wfits=wfits
; jm06junuofa - based on WRITE_ATLAS2D_WRAPPER
    
    datapath = highzea_path()
    outpath = highzea_path(/spec2d)
    
    path = datapath+[$
      '04dec/',$
      '06feb/',$
      '06mar/',$
      '06may/']  ; SINGS+11HUGS
    np = n_elements(path)

    if keyword_set(cleanfits) then begin
       if (file_test(outpath,/directory) ne 0L) then begin
          splog, 'Would you like to remove all FITS files from '+outpath+' [Y/N]?'
          cc = get_kbrd(1)
          if (strupcase(cc) eq 'Y') then spawn, ['/bin/rm -f '+outpath+'*.fits'], /sh
       endif
    endif
    
    for i = 0L, np-1L do begin
       write_atlas2d, 'speclist.txt', datapath=path[i], outpath=outpath, $
         /upcase, overwrite=overwrite, wfits=wfits
    endfor

stop

return
end    
