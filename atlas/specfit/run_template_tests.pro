pro run_template_tests, object
; jm03dec4uofa
; study the completeness of our choice of templates; fit the entire
; atlas with 4, 7, 14, and 20 templates and compare the distribution
; of continuum chi2 values and the fluxes and equivalent widths; for
; 347 galaxies and 7 templates ATLAS1D_SPECFIT takes approximately 70
; minutes to complete.  4, 14, and 20 templates should take,
; respectively, 40, 140, and 200 minutes, or 6.3 hours total
;
; the code to generate the templates for this exercise is:
;
;    write_bc03_templates, dlogage=0.18, tpath=atlas_path(/specfit)+'template_tests/', suffix='20', /write
;    write_bc03_templates, dlogage=0.26, tpath=atlas_path(/specfit)+'template_tests/', suffix='14', /write
;    write_bc03_templates, agegrid=[10.0,100,1000,10000], tpath=atlas_path(/specfit)+'template_tests/', suffix='04', /write

    datapath = atlas_path(/atlas1d)
    atlas = read_integrated()

    if n_elements(object) ne 0L then begin

       doit = match_string(object,atlas.galaxy,index=match) ; match by galaxy name
       if match[0] eq -1L then begin
          splog, 'Object '+object+' not found!'
          return
       endif
       atlas = atlas[match] 
       
    endif

    speclist = atlas.driftfile

; fit only with the solar templates    
    
    eigendir = atlas_path(/specfit)+'template_tests/'
    eigenspec = 'BC03_Z02'+['_04','_14','_20']+'_templates.fits'
    suffix = 'integrated_atlas_Z02'+['_04','_14','_20']
    neigenspec = n_elements(eigenspec)

    for k = 0L, neigenspec-1L do begin

       specdata = atlas1d_specfit(speclist,datapath=datapath,eigendir=eigendir,$
         eigenspec=eigenspec[k],suffix=suffix[k],/postscript,/write)

    endfor

; parse the results

    for k = 0L, neigenspec-1L do atlas1d_parse_specfit, root=suffix[k], datapath=eigendir

stop    
    
return
end    
