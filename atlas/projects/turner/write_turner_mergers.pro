pro write_turner_mergers, cleanfits=cleanfits, overwrite=overwrite, wfits=wfits
; jm05jul01uofa - written
    
    datapath = atlas_path(/spec2dturner)
    outpath = atlas_path(/spec2dturner)+'spec2d/'
    
    path = datapath+[$
      '94nov/',$
      '95mar/',$
      '95oct/',$
      '96apr/',$
      '97apr/']
    np = n_elements(path)

    if keyword_set(cleanfits) then begin
       if (file_test(outpath,/directory) ne 0L) then begin
          splog, 'Would you like to remove all FITS files from '+outpath+' [Y/N]?'
          cc = get_kbrd(1)
          if (strupcase(cc) eq 'Y') then spawn, ['/bin/rm -f '+outpath+'*.fits'], /sh
       endif
       if (file_test(outpath+'rectified/',/directory) ne 0L) then begin
          splog, 'Would you like to remove all FITS files from '+outpath+'rectified/'+' [Y/N]?'
          cc = get_kbrd(1)
          if (strupcase(cc) eq 'Y') then spawn, ['/bin/rm -f '+outpath+'rectified/'+'*.fits'], /sh
       endif
    endif
    
    for i = 0L, np-1L do begin

       write_atlas2d, 'mergerlist.txt', datapath=path[i], outpath=outpath, $
         overwrite=overwrite, wfits=wfits
       write_atlas2d, 'mergerlist.txt', datapath=path[i], outpath=outpath, $
         overwrite=overwrite, /rectified, wfits=wfits

    endfor

stop

return
end    
