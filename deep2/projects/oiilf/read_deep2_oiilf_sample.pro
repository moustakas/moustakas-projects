function read_deep2_oiilf_sample, parent_kcorr=parent_kcorr, kcorr=kcorr, ispec=ispec
; jm08apr12nyu - written, based on READ_DEEP2_MZ_SAMPLE()

    common deep2_oiilf, deep2_parent_kcorr, deep2_kcorr, deep2_ispec

    datapath = deep2_path(/projects)+'oiilf/'

; k-corrections
    
    if keyword_set(parent_kcorr) then begin
       thisfile = datapath+'deep2_oiilf_parent_kcorr.fits'
       if (size(deep2_parent_kcorr,/type) ne 8L) then begin
          if (not keyword_set(silent)) then splog, 'Reading '+thisfile+'.'
          deep2_parent_kcorr = mrdfits(thisfile,1,silent=0)
       endif else if (not keyword_set(silent)) then splog, 'Restoring '+file_basename(thisfile)+'.'
       return, deep2_parent_kcorr
    endif
    
    if keyword_set(kcorr) then begin
       thisfile = datapath+'deep2_oiilf_kcorr.fits'
       if (size(deep2_kcorr,/type) ne 8L) then begin
          if (not keyword_set(silent)) then splog, 'Reading '+thisfile+'.'
          deep2_kcorr = mrdfits(thisfile,1,silent=0)
       endif else if (not keyword_set(silent)) then splog, 'Restoring '+file_basename(thisfile)+'.'
       return, deep2_kcorr
    endif
    
; emission-line data    
    
    if keyword_set(ispec) then begin
       thisfile = datapath+'deep2_oiilf_ispec.fits'
       if (size(deep2_ispec,/type) ne 8L) then begin
          if (not keyword_set(silent)) then splog, 'Reading '+thisfile+'.'
          deep2_ispec = mrdfits(thisfile,1,silent=0)
       endif else if (not keyword_set(silent)) then splog, 'Restoring '+file_basename(thisfile)+'.'
       return, deep2_ispec
    endif

end
