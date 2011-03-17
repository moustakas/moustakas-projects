function read_oiine_sample, sdss=sdss, silent=silent, $
  ancillary=ancillary, ispec=ispec
; jm09apr01 - read the output from BUILD_OIINE_SAMPLE

    common oiine_sample, sdss_ancillary, sdss_ispec, $
      deep2_ancillary, deep2_ispec
    
    datapath = deep2_path(/projects)+'oiine/'
    if keyword_set(sdss) then begin ; SDSS
; ancillary
       if keyword_set(ancillary) then begin
          thisfile = datapath+'sdss_oiine_ancillary.fits.gz'
          if (size(sdss_ancillary,/type) ne 8L) then begin
             if (not keyword_set(silent)) then splog, 'Reading '+thisfile
             sdss_ancillary = mrdfits(thisfile,1,silent=0)
          endif else if (not keyword_set(silent)) then splog, 'Restoring '+file_basename(thisfile)
          return, sdss_ancillary
       endif
; ispec
       if keyword_set(ispec) then begin
          thisfile = datapath+'sdss_oiine_ispec.fits.gz'
          if (size(sdss_ispec,/type) ne 8L) then begin
             if (not keyword_set(silent)) then splog, 'Reading '+thisfile
             sdss_ispec = mrdfits(thisfile,1,silent=0)
          endif else if (not keyword_set(silent)) then splog, 'Restoring '+file_basename(thisfile)
          return, sdss_ispec
       endif
    endif else begin            ; DEEP2
; ancillary
       if keyword_set(ancillary) then begin
          thisfile = datapath+'deep_oiine_ancillary.fits.gz'
          if (size(deep2_ancillary,/type) ne 8L) then begin
             if (not keyword_set(silent)) then splog, 'Reading '+thisfile
             deep2_ancillary = mrdfits(thisfile,1,silent=0)
          endif else if (not keyword_set(silent)) then splog, 'Restoring '+file_basename(thisfile)
          return, deep2_ancillary
       endif
; ispec
       if keyword_set(ispec) then begin
          thisfile = datapath+'deep_oiine_ispec.fits.gz'
          if (size(deep2_ispec,/type) ne 8L) then begin
             if (not keyword_set(silent)) then splog, 'Reading '+thisfile
             deep2_ispec = mrdfits(thisfile,1,silent=0)
          endif else if (not keyword_set(silent)) then splog, 'Restoring '+file_basename(thisfile)
          return, deep2_ispec
       endif
    endelse 
    
end
