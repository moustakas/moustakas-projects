;+
; NAME:
;   READ_AGES()
;
; PURPOSE:
;   Read all the principal AGES catalogs (all containing 39943 rows).
;
; INPUTS: 
;
; KEYWORD PARAMETERS: 
;
; OUTPUTS: 
;
; COMMENTS:
;
; MODIFICATION HISTORY:
;   J. Moustakas, 2010 May 01, UCSD - long previous history!
;
; Copyright (C) 2010, John Moustakas
; 
; This program is free software; you can redistribute it and/or modify 
; it under the terms of the GNU General Public License as published by 
; the Free Software Foundation; either version 2 of the License, or
; (at your option) any later version. 
; 
; This program is distributed in the hope that it will be useful, but 
; WITHOUT ANY WARRANTY; without even the implied warranty of
; MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
; General Public License for more details. 
;-

function read_ages, photo=photo, kcorr=kcorr, ppxf=ppxf, $
  Zmulti=Zmulti, unfluxed=unfluxed, silent=silent

    common ages_read, ages_photo, ages_kcorr, ages_ppxf_zmulti, $
      ages_ppxf_solar, ages_ppxf_unfluxed
    
    mycatpath = ages_path(/mycatalogs)
    ppxfpath = ages_path(/ppxf)
    
; merged photometry
    if keyword_set(photo) then begin
       vv = ages_version(/photo)
       thisfile = mycatpath+'ages_photometry_'+vv+'.fits.gz'
       if (size(ages_photo,/type) ne 8) then begin
          if (keyword_set(silent) eq 0) then splog, 'Reading '+thisfile
          ages_photo = mrdfits(thisfile,1,silent=0)
       endif else if (keyword_set(silent) eq 0) then splog, 'Restoring '+file_basename(thisfile)
       return, ages_photo
    endif

; K-corrections    
    if keyword_set(kcorr) then begin
       vv = ages_version(/kcorr)
       thisfile = mycatpath+'ages_kcorrect_'+vv+'.fits.gz'
       if (size(ages_kcorr,/type) ne 8) then begin
          if (keyword_set(silent) eq 0) then splog, 'Reading '+thisfile
          ages_kcorr = mrdfits(thisfile,1,silent=0)
       endif else if (keyword_set(silent) eq 0) then splog, 'Restoring '+file_basename(thisfile)
       return, ages_kcorr
    endif
    
; see AGES_GANDALF_SPECFIT    
    if keyword_set(ppxf) then begin
       if keyword_set(unfluxed) then begin
          if (size(ages_ppxf_unfluxed,/type) ne 8) then $
            ages_ppxf_unfluxed = read_ages_gandalf(silent=silent,/unfluxed,_extra=extra) else $
              if (keyword_set(silent) eq 0) then splog, 'Restoring unfluxed PPXF file'
          return, ages_ppxf_unfluxed
       endif else begin
          if keyword_set(Zmulti) then begin
             if (size(ages_ppxf_zmulti,/type) ne 8) then $
               ages_ppxf_zmulti = read_ages_gandalf(silent=silent,_extra=extra) else $
                 if (keyword_set(silent) eq 0) then splog, 'Restoring PPXF/Zmulti file'
             return, ages_ppxf_zmulti
          endif else begin
             if (size(ages_ppxf_solar,/type) ne 8) then $
               ages_ppxf_solar = read_ages_gandalf(silent=silent,/solar,_extra=extra) else $
                 if (keyword_set(silent) eq 0) then splog, 'Restoring PPXF/Zsolar file'
             return, ages_ppxf_solar
          endelse
       endelse
    endif

;    if keyword_set(isedfit) then begin
;       thisfile = isedfitpath+'BwRIzK_salp_sfhgrid02_isedfit.fits.gz'
;       if (size(ages_isedfit,/type) ne 8L) then begin
;          if (not keyword_set(silent)) then splog, 'Reading '+thisfile
;          ages_isedfit = mrdfits(thisfile,1,silent=0)
;       endif else if (not keyword_set(silent)) then splog, 'Restoring '+file_basename(thisfile)
;       return, ages_isedfit
;    endif
;    
;    
;    if keyword_set(ispec) then begin ; this needs to be last to not interfere with the keywords
;       thisfile = specfitpath+'ages_specdata_ispec_tweak_'+ages_version(/ispec)+'.fits.gz'
;       if (size(ages_ispec,/type) ne 8L) then begin
;          if (not keyword_set(silent)) then splog, 'Reading '+thisfile
;          ages_ispec = mrdfits(thisfile,1,silent=0)
;       endif else if (not keyword_set(silent)) then splog, 'Restoring '+file_basename(thisfile)
;       return, ages_ispec
;    endif
    
end
