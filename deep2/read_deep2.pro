function read_deep2, kcorr=kcorr, ispec=ispec, silent=silent
; jm07dec04nyu
; jm08mar31nyu - major update

    common deep2_read, deep2_kcorr, deep2_ispec
    
    catpath = deep2_path(/catalogs)

    if keyword_set(kcorr) then begin
       thisfile = catpath+'deep2_kcorr_'+deep2_version(/kcorr)+'.fits.gz'
       if (size(deep2_kcorr,/type) ne 8L) then begin
          if (not keyword_set(silent)) then splog, 'Reading '+thisfile
          deep2_kcorr = mrdfits(thisfile,1,silent=0)
       endif else if (not keyword_set(silent)) then splog, 'Restoring '+file_basename(thisfile)
       return, deep2_kcorr
    endif
    
    if keyword_set(ispec) then begin
       thisfile = catpath+'deep2_specdata_'+deep2_version(/ispec)+'.fits.gz'
       if (size(deep2_ispec,/type) ne 8L) then begin
          if (not keyword_set(silent)) then splog, 'Reading '+thisfile
          deep2_ispec = mrdfits(thisfile,1,silent=0)
       endif else if (not keyword_set(silent)) then splog, 'Restoring '+file_basename(thisfile)
       return, deep2_ispec
    endif
end
