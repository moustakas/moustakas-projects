function read_atlas_gandalf, ppxf=ppxf, nuclear=nuclear, solar=solar, silent=silent
; jm10mar04ucsd - see ATLAS_GANDALF_SPECFIT

;   common atlas_read_gandalf, nuclear_ppxf, drift_ppxf
    
    version = atlas_version(/ppxf_specfit)    
    specfitpath = atlas_path(/ppxf)

    if keyword_set(solar) then suffix = 'solar_' else suffix = ''

    if keyword_set(nuclear) then begin
       thisfile = specfitpath+'atlas_specdata_'+suffix+'nuclear_'+version+'.fits.gz'
       if (size(nuclear_ppxf,/type) ne 8) then begin
          if (not keyword_set(silent)) then splog, 'Reading '+thisfile
          nuclear_ppxf = mrdfits(thisfile,1,silent=0)
       endif else if (not keyword_set(silent)) then splog, 'Restoring '+file_basename(thisfile)
       return, nuclear_ppxf
    endif else begin
       thisfile = specfitpath+'atlas_specdata_'+suffix+'drift_'+version+'.fits.gz'
       if (size(drift_ppxf,/type) ne 8) then begin
          if (not keyword_set(silent)) then splog, 'Reading '+thisfile
          drift_ppxf = mrdfits(thisfile,1,silent=0)
       endif else if (not keyword_set(silent)) then splog, 'Restoring '+file_basename(thisfile)
       return, drift_ppxf
    endelse 

end
