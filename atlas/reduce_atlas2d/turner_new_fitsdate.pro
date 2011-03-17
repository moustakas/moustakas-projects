pro turner_new_fitsdate, root=root, datapath=datapath, update=update
; jm05jun28uofa - convert the old date format in Anne's thesis spectra
;                 (DD/MM/YY) to the new FITS standard (YYYY-MM-DD),
;                 and optionally update the header

    if (n_elements(root) eq 0L) then root = ''
    if (n_elements(datapath) eq 0L) then datapath = cwd()

    pushd, datapath
    flist = file_search(root+'*.fits',count=fcount)
    popd

    if (fcount eq 0L) then begin
       print, 'No FITS files found in path '+datapath
       return
    endif

    for icount = 0L, fcount-1L do begin
       h = headfits(datapath+flist[icount])
       date = sxpar(h,'DATE-OBS',count=datecount)
       if (datecount eq 1L) then begin
          splitdate = strtrim(strsplit(date,'/',/extract),2)
          newdate = '19'+splitdate[2]+'-'+splitdate[1]+'-'+splitdate[0]
          if keyword_set(update) then begin
             sxaddpar, h, 'DATE-OBS', newdate
             djs_modfits, datapath+flist[icount], 0, h
          endif else begin
             splog, 'To be updated: '+flist[icount]+' '+strtrim(date,2)+' --> '+newdate
          endelse
       endif
    endfor
    
return
end
    
