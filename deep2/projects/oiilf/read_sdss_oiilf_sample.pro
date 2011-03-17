function read_sdss_oiilf_sample, parent_kcorr=parent_kcorr, kcorr=kcorr, ispec=ispec
; jm08apr12nyu - written, based on READ_SDSS_MZ_SAMPLE()

    common sdss_oiilf, sdss_parent_kcorr, sdss_kcorr, sdss_ispec

    datapath = deep2_path(/projects)+'oiilf/'

; k-corrections
    
    if keyword_set(parent_kcorr) then begin
       thisfile = datapath+'sdss_oiilf_parent_kcorr.fits'
       if (size(sdss_parent_kcorr,/type) ne 8L) then begin
          if (not keyword_set(silent)) then splog, 'Reading '+thisfile+'.'
          sdss_parent_kcorr = mrdfits(thisfile,1,silent=0)
       endif else if (not keyword_set(silent)) then splog, 'Restoring '+file_basename(thisfile)+'.'
       return, sdss_parent_kcorr
    endif
    
    if keyword_set(kcorr) then begin
       thisfile = datapath+'sdss_oiilf_kcorr.fits'
       if (size(sdss_kcorr,/type) ne 8L) then begin
          if (not keyword_set(silent)) then splog, 'Reading '+thisfile+'.'
          sdss_kcorr = mrdfits(thisfile,1,silent=0)
       endif else if (not keyword_set(silent)) then splog, 'Restoring '+file_basename(thisfile)+'.'
       return, sdss_kcorr
    endif
    
; emission-line data    
    
    if keyword_set(ispec) then begin
       thisfile = datapath+'sdss_oiilf_ispec.fits'
       if (size(sdss_ispec,/type) ne 8L) then begin
          if (not keyword_set(silent)) then splog, 'Reading '+thisfile+'.'
          sdss_ispec = mrdfits(thisfile,1,silent=0)
       endif else if (not keyword_set(silent)) then splog, 'Restoring '+file_basename(thisfile)+'.'
       return, sdss_ispec
    endif

end
