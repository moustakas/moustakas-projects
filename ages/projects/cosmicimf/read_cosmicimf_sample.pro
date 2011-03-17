function read_cosmicimf_sample, sfrs=sfrs, mass=mass, pegase=pegase, $
  witt=witt, silent=silent
; jm10mar16ucsd - read the output from BUILD_COSMICIMF_SAMPLE

;   common cosmicimf_parent, ages_parent, ages_sfrs, $
;     ages_mass, pegase_sfrs, witt_dust
    datapath = ages_path(/projects)+'cosmicimf/'

; star-formation rates    
    if keyword_set(sfrs) then begin
       thisfile = datapath+'ages_cosmicimf_sfrs.fits.gz'
       if (size(ages_sfrs,/type) ne 8) then begin
          if (keyword_set(silent) eq 0) then splog, 'Reading '+thisfile
          ages_sfrs = mrdfits(thisfile,1,silent=0)
       endif else if (keyword_set(silent) eq 0) then $
         splog, 'Restoring '+file_basename(thisfile)
       return, ages_sfrs
    endif

; stellar masses
    if keyword_set(mass) then begin
       thisfile = datapath+'ages_cosmicimf_mass.fits.gz'
       if (size(ages_mass,/type) ne 8) then begin
          if (keyword_set(silent) eq 0) then splog, 'Reading '+thisfile
          ages_mass = mrdfits(thisfile,1,silent=0)
       endif else if (keyword_set(silent) eq 0) then $
         splog, 'Restoring '+file_basename(thisfile)
       return, ages_mass
    endif

; Pegase model results    
    if keyword_set(pegase) then begin
       thisfile = datapath+'pegase_sfrs.fits.gz'
       if (size(pegase_sfrs,/type) ne 8) then begin
          if (keyword_set(silent) eq 0) then splog, 'Reading '+thisfile
          pegase_sfrs = mrdfits(thisfile,1,silent=0)
       endif else if (keyword_set(silent) eq 0) then $
         splog, 'Restoring '+file_basename(thisfile)
       return, pegase_sfrs
    endif

; Witt & Gordon dust models
    if keyword_set(witt) then begin
       thisfile = datapath+'witt_dust.fits.gz'
       if (size(witt_dust,/type) ne 8) then begin
          if (keyword_set(silent) eq 0) then splog, 'Reading '+thisfile
          witt_dust = mrdfits(thisfile,1,silent=0)
       endif else if (keyword_set(silent) eq 0) then $
         splog, 'Restoring '+file_basename(thisfile)
       return, witt_dust
    endif

; by default, read the parent sample    
    thisfile = datapath+'ages_cosmicimf.fits.gz'
    if (size(ages_parent,/type) ne 8) then begin
       if (keyword_set(silent) eq 0) then splog, 'Reading '+thisfile
       ages_parent = mrdfits(thisfile,1,silent=0)
    endif else if (keyword_set(silent) eq 0) then $
      splog, 'Restoring '+file_basename(thisfile)
    return, ages_parent

end
