function read_sings, nuclear=nuclear, drift20=drift20,drift56=drift56, $
  version=version, silent=silent
; jm04sep02uofa - written
; jm05jul26uofa - updated
; jm08oct20nyu - updated to the latest data model

    common sings_read, sings_nuclear, sings_drift20, sings_drift56, ispec_version
    
    specfitpath = sings_path(/specfit)
    version = sings_version(/specfit)
    
    if keyword_set(nuclear) then begin
       thisfile = specfitpath+'sings_nuclear_speclinefit_'+version+'.fits.gz'
       if (size(sings_nuclear,/type) ne 8L) then begin
          if (not keyword_set(silent)) then splog, 'Reading '+thisfile
          sings_nuclear = mrdfits(thisfile,1,silent=0)
       endif else if (not keyword_set(silent)) then splog, 'Restoring '+file_basename(thisfile)
       return, sings_nuclear
    endif

    if keyword_set(drift20) then begin
       thisfile = specfitpath+'sings_drift20_speclinefit_'+version+'.fits.gz'
       if (size(sings_drift20,/type) ne 8L) then begin
          if (not keyword_set(silent)) then splog, 'Reading '+thisfile
          sings_drift20 = mrdfits(thisfile,1,silent=0)
       endif else if (not keyword_set(silent)) then splog, 'Restoring '+file_basename(thisfile)
       return, sings_drift20
    endif

    if keyword_set(drift56) then begin
       thisfile = specfitpath+'sings_drift56_speclinefit_'+version+'.fits.gz'
       if (size(sings_drift56,/type) ne 8L) then begin
          if (not keyword_set(silent)) then splog, 'Reading '+thisfile
          sings_drift56 = mrdfits(thisfile,1,silent=0)
       endif else if (not keyword_set(silent)) then splog, 'Restoring '+file_basename(thisfile)
       return, sings_drift56
    endif

end
    
