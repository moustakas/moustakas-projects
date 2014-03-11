function read_deep2, kcorr=kcorr, ispec=ispec, ppxf=ppxf, fixoii=fixoii, silent=silent
; jm07dec04nyu
; jm08mar31nyu - major update

    common deep2_read, deep2_kcorr, deep2_ispec, deep2_ppxf, deep2_ppxf_fixoii
    
    catpath = deep2_path(/catalogs)
    if keyword_set(kcorr) then begin
       thisfile = catpath+'kcorr.dr4.goodspec1d.Q34.fits.gz'
       if (size(deep2_kcorr,/type) ne 8L) then begin
          if (not keyword_set(silent)) then splog, 'Reading '+thisfile
          deep2_kcorr = mrdfits(thisfile,1,silent=0)
       endif else if (not keyword_set(silent)) then splog, 'Restoring '+file_basename(thisfile)
       return, deep2_kcorr
    endif
    
    if keyword_set(ispec) then begin
       thisfile = catpath+'speclinefit.dr4.goodspec1d.Q34.fits.gz'
       if (size(deep2_ispec,/type) ne 8L) then begin
          if (not keyword_set(silent)) then splog, 'Reading '+thisfile
          deep2_ispec = mrdfits(thisfile,1,silent=0)
       endif else if (not keyword_set(silent)) then splog, 'Restoring '+file_basename(thisfile)
       return, deep2_ispec
    endif

    if keyword_set(ppxf) and keyword_set(fixoii) eq 0 then begin
       thisfile = catpath+'deep2.ppxf.specdata.dr4_v1.0.fits.gz'
       if (size(deep2_ppxf,/type) ne 8L) then begin
          if (not keyword_set(silent)) then splog, 'Reading '+thisfile
          deep2_ppxf = mrdfits(thisfile,1,silent=0)
       endif else if (not keyword_set(silent)) then splog, 'Restoring '+file_basename(thisfile)
       return, deep2_ppxf
    endif

    if keyword_set(ppxf) and keyword_set(fixoii) then begin
       thisfile = catpath+'deep2.ppxf.specdata.fixoii.dr4_v1.0.fits.gz'
       if (size(deep2_ppxf_fixoii,/type) ne 8L) then begin
          if (not keyword_set(silent)) then splog, 'Reading '+thisfile
          deep2_ppxf_fixoii = mrdfits(thisfile,1,silent=0)
       endif else if (not keyword_set(silent)) then splog, 'Restoring '+file_basename(thisfile)
       return, deep2_ppxf_fixoii
    endif
    
end
