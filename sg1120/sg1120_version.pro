function sg1120_version, kcorrect=kcorrect, parent_catalog=parent_catalog
; jm08jun09nyu - written
;
; kcorrect
;   v1.0 - earlier effort
;   v1.1 - properly implemented VIMOS filters; updated the
;          k-correction code to use SG1120_KCORRECT
;   v2.0 - new VIMOS and LDSS3 reductions and zeropoints 
;   v3.0 - added Ks-band photometry and iterative tweaking of the
;     zeropoints 
; parent catalog (see BUILD_SG1120_PARENT_CATALOG)
;   v1.0 - earlier effort
;   v2.0 - based on a complete rereduction of the VIMOS, LDSS3, and
;     FLAMINGOS imaging; includes HST and SDSS ancillary data
;   v3.0 - include Ks-band photometry; rederived zeropoints 

    version = -1.0
    if keyword_set(kcorrect) then version = 'v3.0'
    if keyword_set(parent_catalog) then version = 'v3.0'

return, version
end
