pro write_atlas2d_wrapper, atlas=atlas, sings=sings, hii=hii, $
  cleanfits=cleanfits, overwrite=overwrite, wfits=wfits
; jm02oct25uofa
; jm03dec25uofa - minor updates
; jm04may01uofa - update the absolute error budget
; jm05jun29uofa - updated for iSPEC2d v2.0 reductions 
; jm06jan24uofa - added STITCHED support for the SINGS spectra; added
;                 05APR directory
; jm06jul31uofa - added 06MAR and 06MAY directories, and support for
;                 the SINGS/HII observations
; jm08jan17nyu - do not distinguish between rectified and
;                non-rectified spectra; *only* write out the rectified
;                spectra!; added support for STITCHED atlas spectra
;                (specifically, NGC3344 and NGC5194 in 99may)
    
; send the spectral atlas, SINGS, GTO starbursts, and miscellaneous
; spectra to individual output directories

    datapath = atlas_path(/spec2datlas)
    singspath = sings_path(/spec2d)
    hiipath = sings_hii_path(/spec2d)
    atlaspath = atlas_path(/atlas2d)
    
    path = datapath+[$
      '94nov/',$ ; Turner
      '95mar/',$ ; Turner
      '95oct/',$ ; Turner
      '96apr/',$ ; Turner
      '97apr/',$ ; Turner
      '98mar/',$
      '98apr/',$
      '98jun/',$
      '98oct/',$
      '99apr/',$
      '99may/',$
      '99nov/',$
      '00apr/',$
      '01nov/',$ ; SINGS+ATLAS2D
      '01dec/',$ ; SINGS+ATLAS2D
      '02feb/',$ ; SINGS+ATLAS2D
      '02apr/',$
      '02may/',$ ; SINGS+ATLAS2D
      '03may/',$ ; ATLAS2D+GTO
      '05apr/',$ ; SINGS+11HUGS
      '06mar/',$ ; SINGS+11HUGS
      '06may/']  ; SINGS+11HUGS
;   path = datapath+'05apr/'
    np = n_elements(path)

    if keyword_set(atlas) then begin
       
       if keyword_set(cleanfits) then begin
          if (file_test(atlaspath,/directory) ne 0L) then begin
             splog, 'Would you like to remove all FITS files from '+atlaspath+' [Y/N]?'
             cc = get_kbrd(1)
             if (strupcase(cc) eq 'Y') then spawn, ['/bin/rm -f '+atlaspath+'*.fits'], /sh
          endif
;         if (file_test(atlaspath+'rectified/',/directory) ne 0L) then begin
;            splog, 'Would you like to remove all FITS files from '+atlaspath+'rectified/'+' [Y/N]?'
;            cc = get_kbrd(1)
;            if (strupcase(cc) eq 'Y') then spawn, ['/bin/rm -f '+atlaspath+'rectified/'+'*.fits'], /sh
;         endif
       endif
       
       for i = 0L, np-1L do begin

          write_atlas2d, 'atlaslist.txt', datapath=path[i], outpath=atlaspath, $
            overwrite=overwrite, /rectified, wfits=wfits
          write_atlas2d, 'atlaslist_stitched.txt', datapath=path[i], outpath=atlaspath, $
            overwrite=overwrite, /stitched, wfits=wfits

;         write_atlas2d, 'atlaslist.txt', datapath=path[i], outpath=atlaspath, $
;           overwrite=overwrite, wfits=wfits
;         write_atlas2d, 'atlaslist.txt', datapath=path[i], outpath=atlaspath+'rectified/', $
;           overwrite=overwrite, /rectified, wfits=wfits

       endfor

    endif

; SINGS    
    
    if keyword_set(sings) then begin
       
       if keyword_set(cleanfits) then begin
          if (file_test(singspath,/directory) ne 0L) then begin
             splog, 'Would you like to remove all FITS files from '+singspath+' [Y/N]?'
             cc = get_kbrd(1)
             if (strupcase(cc) eq 'Y') then spawn, ['/bin/rm -f '+singspath+'*.fits'], /sh
          endif
;         if (file_test(singspath+'rectified/',/directory) ne 0L) then begin
;            splog, 'Would you like to remove all FITS files from '+singspath+'rectified/'+' [Y/N]?'
;            cc = get_kbrd(1)
;            if (strupcase(cc) eq 'Y') then spawn, ['/bin/rm -f '+singspath+'rectified/'+'*.fits'], /sh
;         endif
       endif
       
;      for i = 20L, np-1L do begin
       for i = 0L, np-1L do begin

          write_atlas2d, 'singslist.txt', datapath=path[i], outpath=singspath, $
            overwrite=overwrite, /rectified, wfits=wfits
          write_atlas2d, 'singslist_stitched.txt', datapath=path[i], outpath=singspath, $
            overwrite=overwrite, /stitched, wfits=wfits

;         write_atlas2d, 'singslist.txt', datapath=path[i], outpath=singspath, $
;           overwrite=overwrite, wfits=wfits
;         write_atlas2d, 'singslist.txt', datapath=path[i], outpath=singspath+'rectified/', $
;           overwrite=overwrite, /rectified, wfits=wfits

       endfor

    endif
       
; SINGS/HII
    
    if keyword_set(hii) then begin
       
       if keyword_set(cleanfits) then begin
          if (file_test(hiipath,/directory) ne 0L) then begin
             splog, 'Would you like to remove all FITS files from '+hiipath+' [Y/N]?'
             cc = get_kbrd(1)
             if (strupcase(cc) eq 'Y') then spawn, ['/bin/rm -f '+hiipath+'*.fits'], /sh
          endif
       endif
       
;      for i = 20L, np-1L do begin
       for i = 0L, np-1L do begin

          write_atlas2d, 'hiilist.txt', datapath=path[i], outpath=hiipath, $
            overwrite=overwrite, wfits=wfits, /rectified

       endfor

    endif
       
; these observations are OBSOLETE    
    
;; write the NFGS and the stellar template observations (2002 June)
;
;    write_atlas2d, 'stellar_template_list.txt', datapath=datapath+'02jun/', $
;      outpath=atlas_path(/templates2d), /overwrite, wfits=wfits
;
;    write_atlas2d, 'nfgslist.txt', datapath=datapath+'02jun/', $
;      outpath=atlas_path(/analysis1d)+'nfgs/spec2d/', /overwrite, wfits=wfits

;      write_atlas2d, 'gtolist.txt', datapath=path[i], $
;        outpath=gtopath, overwrite=overwrite, wfits=wfits
;     
;      write_atlas2d, 'speclist.txt', datapath=path[i], $
;        outpath=projectpath+'misc/', /overwrite, wfits=wfits

return
end    
