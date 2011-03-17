function read_ediscs_mz_sample, parent_ancillary=parent_ancillary, mz_ancillary=mz_ancillary, $
  mz_ispec=mz_ispec, mz_log12oh=mz_log12oh
; jm08apr07nyu - written, based on READ_AGES_MZ_SAMPLE

    common ediscs_mz, ediscs_parent_ancillary, ediscs_mz_ancillary, $
      ediscs_mz_ispec, ediscs_mz_log12oh

    datapath = ediscs_path(/projects)+'mz/'

; k-corrections
    
    if keyword_set(parent_ancillary) then begin
       thisfile = datapath+'ediscs_parent_ancillary.fits'
       if (size(ediscs_parent_ancillary,/type) ne 8L) then begin
          if (not keyword_set(silent)) then splog, 'Reading '+thisfile+'.'
          ediscs_parent_ancillary = mrdfits(thisfile,1,silent=0)
       endif else if (not keyword_set(silent)) then splog, 'Restoring '+file_basename(thisfile)+'.'
       return, ediscs_parent_ancillary
    endif
    
    if keyword_set(mz_ancillary) then begin
       thisfile = datapath+'ediscs_mz_ancillary.fits'
       if (size(ediscs_mz_ancillary,/type) ne 8L) then begin
          if (not keyword_set(silent)) then splog, 'Reading '+thisfile+'.'
          ediscs_mz_ancillary = mrdfits(thisfile,1,silent=0)
       endif else if (not keyword_set(silent)) then splog, 'Restoring '+file_basename(thisfile)+'.'
       return, ediscs_mz_ancillary
    endif
    
; emission-line data    
    
    if keyword_set(mz_ispec) then begin
       thisfile = datapath+'ediscs_mz_ispec.fits'
       if (size(ediscs_mz_ispec,/type) ne 8L) then begin
          if (not keyword_set(silent)) then splog, 'Reading '+thisfile+'.'
          ediscs_mz_ispec = mrdfits(thisfile,1,silent=0)
       endif else if (not keyword_set(silent)) then splog, 'Restoring '+file_basename(thisfile)+'.'
       return, ediscs_mz_ispec
    endif

; oxygen abundances

    if keyword_set(mz_log12oh) then begin
       thisfile = datapath+'ediscs_mz_log12oh.fits'
       if (size(ediscs_mz_log12oh,/type) ne 8L) then begin
          if (not keyword_set(silent)) then splog, 'Reading '+thisfile+'.'
          ediscs_mz_log12oh = mrdfits(thisfile,1,silent=0)
       endif else if (not keyword_set(silent)) then splog, 'Restoring '+file_basename(thisfile)+'.'
       return, ediscs_mz_log12oh
    endif
    
end
