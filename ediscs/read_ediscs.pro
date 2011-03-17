function read_ediscs, spec1dinfo=spec1dinfo, photo=photo, $
  specfit=specfit, ppxf=ppxf, kcorr=kcorr, synthmags=synthmags, $
  silent=silent
; jm06sep27nyu
; jm08mar31nyu - major update

    common ediscs_read, ediscs_photo, ediscs_kcorr, ediscs_spec1dinfo, $
      ediscs_specfit, ediscs_ppxf
    
    mycatpath = ediscs_path(/mycatalogs)
    specfitpath = ediscs_path(/specfit)
    ppxfpath = ediscs_path(/ppxf)

; see UNPACK_EDISCS_SPEC1D
    if keyword_set(spec1dinfo) then begin
       thisfile = mycatpath+'ediscs_spec1d_info.fits.gz'
       if (size(ediscs_spec1dinfo,/type) ne 8) then begin
          if (keyword_set(silent) eq 0) then splog, 'Reading '+thisfile
          ediscs_spec1dinfo = mrdfits(thisfile,1,silent=0)
       endif else if (keyword_set(silent) eq 0) then splog, 'Restoring '+file_basename(thisfile)
       return, ediscs_spec1dinfo
    endif

; see BUILD_EDISCS_PHOTOMETRY    
    if keyword_set(photo) then begin
       thisfile = mycatpath+'ediscs_photometry.v23.fits.gz'
       if (size(ediscs_photo,/type) ne 8) then begin
          if (keyword_set(silent) eq 0) then splog, 'Reading '+thisfile
          ediscs_photo = mrdfits(thisfile,1,silent=0)
       endif else if (keyword_set(silent) eq 0) then splog, 'Restoring '+file_basename(thisfile)
       return, ediscs_photo
    endif

; see EDISCS_GANDALF_SPECFIT    
    if keyword_set(ppxf) then begin
       if (size(ediscs_ppxf,/type) ne 8) then $
         ediscs_ppxf = read_ediscs_gandalf(silent=silent,_extra=extra) else $
           if (keyword_set(silent) eq 0) then splog, 'Restoring PPXF file'
       return, ediscs_ppxf
    endif

; see BUILD_EDISCS_KCORRECT
    if keyword_set(kcorr) then begin
       thisfile = mycatpath+'ediscs_kcorrect_'+ediscs_version(/kcorr)+'.fits.gz'
       if (size(ediscs_kcorr,/type) ne 8) then begin
          if (keyword_set(silent) eq 0) then splog, 'Reading '+thisfile
          ediscs_kcorr = mrdfits(thisfile,1,silent=0)
       endif else if (keyword_set(silent) eq 0) then splog, 'Restoring '+file_basename(thisfile)
       return, ediscs_kcorr
    endif

; see EDISCS_SYNTHMAGS
    if keyword_set(synthmags) then begin
       thisfile = mycatpath+'ediscs_synthmags_'+ediscs_version(/specfit)+'.fits.gz'
       if (size(ediscs_synthmags,/type) ne 8) then begin
          if (keyword_set(silent) eq 0) then splog, 'Reading '+thisfile
          ediscs_synthmags = mrdfits(thisfile,1,silent=0)
       endif else if (keyword_set(silent) eq 0) then splog, 'Restoring '+file_basename(thisfile)
       return, ediscs_synthmags
    endif

; obsolete!!    
    if keyword_set(specfit) then begin
       thisfile = specfitpath+'ediscs_specdata_'+ediscs_version(/specfit)+'.fits.gz'
       if (size(ediscs_specfit,/type) ne 8) then begin
          if (keyword_set(silent) eq 0) then splog, 'Reading '+thisfile
          ediscs_specfit = mrdfits(thisfile,1,silent=0)
       endif else if (keyword_set(silent) eq 0) then splog, 'Restoring '+file_basename(thisfile)
       return, ediscs_specfit
    endif

end
