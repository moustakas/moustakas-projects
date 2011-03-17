function ediscs_version, kcorrect=kcorrect, $ ; ancillary=ancillary, $
  specfit=specfit, ppxf_templates=ppxf_templates, $
  ppxf_specfit=ppxf_specfit
; jm08feb05nyu 
; jm08mar24nyu - changed
;
; kcorrect 
;   v3.0 - 
;   v3.1 - now I'm used properly scaled totals magnitudes (see
;     BUILD_EDISCS_PHOTOMETRY)  
;   v4.0 - code rewritten but using the same photometry
; ancillary
;   v1.0 - earlier effort
;   v2.0 - mainly updated k-corrections
;   v2.1 - additional catalog entries (e.g., cluster properties, etc.)
;   v3.0 - use IM_KCORRECT()
;   v3.1 - the content and organization of the structure has changed,
;     and the v3.1 K-corrections have been adopted
;   v3.2 - incorporates the v2.2 synthesized magnitudes
; specfit
;   v1.0 - earlier effort
;   v1.1 - fixed spectral index bug; new BC03 templates
;   v2.0 - rewritten ispec code
;   v2.1 - exclusively solar-metallicity models; also fit H-delta in
;     emission  
;   v2.2 - do not tweak from the EDisCS redshift

;   ppxf_templates 
;     v1.0 - see SINGS_GANDALF_SPECFIT
;   ppxf_specfit
;     v1.0 - see SINGS_GANDALF_SPECFIT

    version = -1.0
    if keyword_set(kcorrect) then version = 'v4.0'
;   if keyword_set(ancillary) then version = 'v3.2'
    if keyword_set(specfit) then version = 'v2.2'

    if keyword_set(ppxf_templates) then version = 'v1.0'
    if keyword_set(ppxf_specfit) then version = 'v1.0'

return, version
end
