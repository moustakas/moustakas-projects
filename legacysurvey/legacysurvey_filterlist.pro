;+
; NAME:
;   LEGACYSURVEY_FILTERLIST()
; PURPOSE:
;   Return the DECam-Legacy filters.
; MODIFICATION HISTORY:
;   J. Moustakas, 2014 Sep 15, Siena
;-

function legacysurvey_filterlist, north=north
    if keyword_set(north) then begin
       filt = ['bass_'+['g','r'], 'mzls_z']+'.par'
    endif else begin
       filt = 'decam_'+['g','r','z']+'.par'
    endelse
return, filt
end
