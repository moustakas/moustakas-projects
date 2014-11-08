;+
; NAME:
;   DECAMLEGACY_FILTERLIST()
; PURPOSE:
;   Return the DECam-Legacy filters.
; MODIFICATION HISTORY:
;   J. Moustakas, 2014 Sep 15, Siena
;-

function decamlegacy_filterlist
    return, 'decam_'+['g','r','z']+'.par'
;   return, ['decam_'+['g','r','z']+'.par',sdss_filterlist()]
end
