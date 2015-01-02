function atlas_version, ancillary=ancillary, specfit=specfit, $
  templates=templates, ppxf_templates=ppxf_templates, $
  ppxf_specfit=ppxf_specfit
; jm08jan18nyu
; jm08apr07nyu - renamed atlas_version() and v3.1 implemented 
;
; ancillary
;    v1.0 - original
;    v1.1 - fixed some coordinates and the redshift for
;           UGC09425NW=UGC09425NED01 
; specfit
;    v1.0 - original thesis data
;    v2.0 - bug! (overestimated errors due to use of SYSERR); better
;           treatment of broad-line AGN
;    v3.0 - new 2D reductions (distortion maps, better sky subtraction,
;           slightly modified error maps); more BC03 templates
;    v3.1 - bug fix in how the three spectra of NGC3690/ARP299 were
;           being fitted 
; templates 
;    v1.0 - not really original
;    v1.1 - minor tweaks to the BC03 info structure; write out a
;           limited wavelength range
;
    
    if keyword_set(ancillary) then version = 'v1.1'
    if keyword_set(specfit) then version = 'v3.1'
    if keyword_set(templates) then version = 'v1.1'

    if keyword_set(ppxf_templates) then version = 'v1.0'
    if keyword_set(ppxf_specfit) then version = 'v1.0'
    if keyword_set(ppxf_ancillary) then version = 'v1.0'

return, version
end
