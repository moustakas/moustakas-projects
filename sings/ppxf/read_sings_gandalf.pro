function read_sings_gandalf, ppxf=ppxf, nuclear=nuclear, $
  drift20=drift20, drift56=drift56, solar=solar, silent=silent
; jm10mar04ucsd - see SINGS_GANDALF_SPECFIT

;   common sings_read_gandalf, nuclear_ppxf, drift20_ppxf, drift56_ppxf
    
    if (keyword_set(nuclear) eq 0) and (keyword_set(drift20) eq 0) and $
      (keyword_set(drift56) eq 0) then begin
       splog, 'Either /NUCLEAR, DRIFT20, or DRIFT56 must be set!'
       return, -1
    endif

    version = sings_version(/ppxf_specfit)    
    specfitpath = sings_path(/ppxf)

    if keyword_set(solar) then suffix = 'solar_' else suffix = ''

    if keyword_set(nuclear) then begin
       thisfile = specfitpath+'sings_specdata_'+suffix+'nuclear_'+version+'.fits.gz'
       if (size(nuclear_ppxf,/type) ne 8) then begin
          if (not keyword_set(silent)) then splog, 'Reading '+thisfile
          nuclear_ppxf = mrdfits(thisfile,1,silent=0)
       endif else if (not keyword_set(silent)) then splog, 'Restoring '+file_basename(thisfile)
       return, nuclear_ppxf
    endif
    
    if keyword_set(drift20) then begin
       thisfile = specfitpath+'sings_specdata_'+suffix+'drift20_'+version+'.fits.gz'
       if (size(drift20_ppxf,/type) ne 8) then begin
          if (not keyword_set(silent)) then splog, 'Reading '+thisfile
          drift20_ppxf = mrdfits(thisfile,1,silent=0)
       endif else if (not keyword_set(silent)) then splog, 'Restoring '+file_basename(thisfile)
       return, drift20_ppxf
    endif
    
    if keyword_set(drift56) then begin
       thisfile = specfitpath+'sings_specdata_'+suffix+'drift56_'+version+'.fits.gz'
       if (size(drift56_ppxf,/type) ne 8) then begin
          if (not keyword_set(silent)) then splog, 'Reading '+thisfile
          drift56_ppxf = mrdfits(thisfile,1,silent=0)
       endif else if (not keyword_set(silent)) then splog, 'Restoring '+file_basename(thisfile)
       return, drift56_ppxf
    endif
    
end
