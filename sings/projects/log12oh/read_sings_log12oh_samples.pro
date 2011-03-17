function read_sings_log12oh_samples, nuclear=nuclear, drift20=drift20, drift56=drift56, $
  nodust_nuclear=nodust_nuclear, nodust_drift20=nodust_drift20, nodust_drift56=nodust_drift56, $
  sdss=sdss, parent_hii=parent_hii, hii=hii
; jm06feb10uofa - see WRITE_SINGS_LOG12OH_SAMPLES
; jm08feb07nyu - added JUSTNODUST keyword
; jm08oct20nyu - major update

    common read_sings_log12oh, log12oh_sdss

;   common read_sings_log12oh, log12oh_nuclear, log12oh_drift20, log12oh_drift56, $
;     log12oh_nodust_nuclear, log12oh_nodust_drift20, log12oh_nodust_drift56, $
;     log12oh_sdss, log12oh_parent_hii, log12oh_hii
    
    version = sings_log12oh_version()
    specfitpath = sings_path(/projects)+'log12oh/'

    if keyword_set(sdss) then begin
       thisfile = specfitpath+'sdss_log12oh_'+version+'.fits.gz' 
       if (size(log12oh_sdss,/type) ne 8L) then begin
          if (not keyword_set(silent)) then splog, 'Reading '+thisfile
          log12oh_sdss = mrdfits(thisfile,1,silent=0)
       endif else if (not keyword_set(silent)) then splog, 'Restoring '+file_basename(thisfile)
       return, log12oh_sdss
    endif

; this file corresponds to the output from sings_log12oh    
    if keyword_set(hii) then begin
       thisfile = specfitpath+'sings_hiiregions_'+version+'.fits.gz' 
       if (size(log12oh_hii,/type) ne 8L) then begin
          if (not keyword_set(silent)) then splog, 'Reading '+thisfile
          log12oh_hii = mrdfits(thisfile,1,silent=0)
       endif else if (not keyword_set(silent)) then splog, 'Restoring '+file_basename(thisfile)
       return, log12oh_hii
    endif

; this is the parent sample written out by write_sings_log12oh_hiiregions
    if keyword_set(parent_hii) then begin
       thisfile = specfitpath+'sings_hiiregions_'+version+'.fits.gz' 
       if (size(log12oh_parent_hii,/type) ne 8L) then begin
          if (not keyword_set(silent)) then splog, 'Reading '+thisfile
          log12oh_parent_hii = mrdfits(thisfile,1,silent=0)
       endif else if (not keyword_set(silent)) then splog, 'Restoring '+file_basename(thisfile)
       return, log12oh_parent_hii
    endif

    if keyword_set(nuclear) then begin
       thisfile = specfitpath+'sings_log12oh_nuclear_speclinefit_'+version+'.fits.gz' 
       if (size(log12oh_nuclear,/type) ne 8L) then begin
          if (not keyword_set(silent)) then splog, 'Reading '+thisfile
          log12oh_nuclear = mrdfits(thisfile,1,silent=0)
       endif else if (not keyword_set(silent)) then splog, 'Restoring '+file_basename(thisfile)
       return, log12oh_nuclear
    endif

    if keyword_set(nodust_nuclear) then begin
       thisfile = specfitpath+'sings_log12oh_nuclear_speclinefit_'+version+'_nodust.fits.gz' 
       if (size(log12oh_nodust_nuclear,/type) ne 8L) then begin
          if (not keyword_set(silent)) then splog, 'Reading '+thisfile
          log12oh_nodust_nuclear = mrdfits(thisfile,1,silent=0)
       endif else if (not keyword_set(silent)) then splog, 'Restoring '+file_basename(thisfile)
       return, log12oh_nodust_nuclear
    endif

    if keyword_set(drift20) then begin
       thisfile = specfitpath+'sings_log12oh_drift20_speclinefit_'+version+'.fits.gz' 
       if (size(log12oh_drift20,/type) ne 8L) then begin
          if (not keyword_set(silent)) then splog, 'Reading '+thisfile
          log12oh_drift20 = mrdfits(thisfile,1,silent=0)
       endif else if (not keyword_set(silent)) then splog, 'Restoring '+file_basename(thisfile)
       return, log12oh_drift20
    endif

    if keyword_set(nodust_drift20) then begin
       thisfile = specfitpath+'sings_log12oh_drift20_speclinefit_'+version+'_nodust.fits.gz' 
       if (size(log12oh_nodust_drift20,/type) ne 8L) then begin
          if (not keyword_set(silent)) then splog, 'Reading '+thisfile
          log12oh_nodust_drift20 = mrdfits(thisfile,1,silent=0)
       endif else if (not keyword_set(silent)) then splog, 'Restoring '+file_basename(thisfile)
       return, log12oh_nodust_drift20
    endif

    if keyword_set(drift56) then begin
       thisfile = specfitpath+'sings_log12oh_drift56_speclinefit_'+version+'.fits.gz' 
       if (size(log12oh_drift56,/type) ne 8L) then begin
          if (not keyword_set(silent)) then splog, 'Reading '+thisfile
          log12oh_drift56 = mrdfits(thisfile,1,silent=0)
       endif else if (not keyword_set(silent)) then splog, 'Restoring '+file_basename(thisfile)
       return, log12oh_drift56
    endif

    if keyword_set(nodust_drift56) then begin
       thisfile = specfitpath+'sings_log12oh_drift56_speclinefit_'+version+'_nodust.fits.gz' 
       if (size(log12oh_nodust_drift56,/type) ne 8L) then begin
          if (not keyword_set(silent)) then splog, 'Reading '+thisfile
          log12oh_nodust_drift56 = mrdfits(thisfile,1,silent=0)
       endif else if (not keyword_set(silent)) then splog, 'Restoring '+file_basename(thisfile)
       return, log12oh_nodust_drift56
    endif

end
    
