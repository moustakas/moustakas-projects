function sings_log12oh_version
; jm06apr16uofa
; v5.3 -
; v6.0 - major update; spectra totally re-reduced (new/better error
;        spectra), updated and expanded HII-region database; improved
;        sample selection
; v7.0 - spectra refitted one final time using the new version of
;        ispec (e.g., better continuum fitting and deblending of
;        the Broad-line AGN)
; v7.1 - updated HII region structure (v2.0)
; v8.0 - version for submitted paper; updated HII-region database; new
;   SDSS sample; new PPXF fitting of the SINGS spectra
; v9.0 - version for resubmitted paper; general QAchecking and small
;   tweaks in how the gradients and average abundances are computed 
; v10.0 - new MONTE_LOG12OH_KK04 code
    
    version = 'v10.0'

return, version
end
