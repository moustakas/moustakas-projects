function read_tkrs, zcat=zcat, kcorr=kcorr, silent=silent
; jm08apr24nyu - written

    common tkrs_read, tkrs_zcat, tkrs_kcorr
    
    path = tkrs_path()

    if keyword_set(zcat) then begin
       thisfile = path+'tkrs_zcat.fits.gz'
       if (size(tkrs_zcat,/type) ne 8L) then begin
          if (not keyword_set(silent)) then splog, 'Reading '+thisfile
          tkrs_zcat = mrdfits(thisfile,1,silent=0)
       endif else if (not keyword_set(silent)) then splog, 'Restoring '+file_basename(thisfile)
       return, tkrs_zcat
    endif

    if keyword_set(kcorr) then begin
       thisfile = path+'tkrs_kcorr.fits.gz'
       if (size(tkrs_kcorr,/type) ne 8L) then begin
          if (not keyword_set(silent)) then splog, 'Reading '+thisfile
          tkrs_kcorr = mrdfits(thisfile,1,silent=0)
       endif else if (not keyword_set(silent)) then splog, 'Restoring '+file_basename(thisfile)
       return, tkrs_kcorr
    endif

end
