function read_deep2_mz_sample, parent_kcorr=parent_kcorr, mz_kcorr=mz_kcorr, $
  mz_ispec=mz_ispec, mz_log12oh=mz_log12oh
; jm08apr07nyu - written, based on READ_AGES_MZ_SAMPLE

    common deep2_mz, deep2_parent_kcorr, deep2_mz_kcorr, $
      deep2_mz_ispec, deep2_mz_log12oh

    datapath = deep2_path(/projects)+'mz/'

; k-corrections
    
    if keyword_set(parent_kcorr) then begin
       thisfile = datapath+'deep2_parent_kcorr.fits'
       if (size(deep2_parent_kcorr,/type) ne 8L) then begin
          if (not keyword_set(silent)) then splog, 'Reading '+thisfile+'.'
          deep2_parent_kcorr = mrdfits(thisfile,1,silent=0)
       endif else if (not keyword_set(silent)) then splog, 'Restoring '+file_basename(thisfile)+'.'
       return, deep2_parent_kcorr
    endif
    
    if keyword_set(mz_kcorr) then begin
       thisfile = datapath+'deep2_mz_kcorr.fits'
       if (size(deep2_mz_kcorr,/type) ne 8L) then begin
          if (not keyword_set(silent)) then splog, 'Reading '+thisfile+'.'
          deep2_mz_kcorr = mrdfits(thisfile,1,silent=0)
       endif else if (not keyword_set(silent)) then splog, 'Restoring '+file_basename(thisfile)+'.'
       return, deep2_mz_kcorr
    endif
    
; emission-line data    
    
    if keyword_set(mz_ispec) then begin
       thisfile = datapath+'deep2_mz_ispec.fits'
       if (size(deep2_mz_ispec,/type) ne 8L) then begin
          if (not keyword_set(silent)) then splog, 'Reading '+thisfile+'.'
          deep2_mz_ispec = mrdfits(thisfile,1,silent=0)
       endif else if (not keyword_set(silent)) then splog, 'Restoring '+file_basename(thisfile)+'.'
       return, deep2_mz_ispec
    endif

; oxygen abundances

    if keyword_set(mz_log12oh) then begin
       thisfile = datapath+'deep2_mz_log12oh.fits'
       if (size(deep2_mz_log12oh,/type) ne 8L) then begin
          if (not keyword_set(silent)) then splog, 'Reading '+thisfile+'.'
          deep2_mz_log12oh = mrdfits(thisfile,1,silent=0)
       endif else if (not keyword_set(silent)) then splog, 'Restoring '+file_basename(thisfile)+'.'
       return, deep2_mz_log12oh
    endif
    
end
