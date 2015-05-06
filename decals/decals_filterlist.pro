;+
; NAME:
;   DECALS_FILTERLIST()
; PURPOSE:
;   Return the DECam-Legacy filters.
; MODIFICATION HISTORY:
;   J. Moustakas, 2014 Sep 15, Siena
;-

function decals_filterlist, grz=grz
    filt = 'decam_'+['u','g','r','i','z','Y']+'.par'
    if keyword_set(grz) then filt = 'decam_'+['g','r','z']+'.par'
;   return, ['decam_'+['g','r','z']+'.par',sdss_filterlist()]
return, filt
end
