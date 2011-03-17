pro read_model_grids, filters, zarray, filterflux, colorflux, kcorrection, templates=templates
;+
; NAME:
;	READ_MODEL_GRIDS		
;
; PURPOSE:
;	Extract the broadband flux, the color, and the k-correction
;	for a subset of the known filters, given set of templates.
;
; INPUTS:
;	filters - filter observations to read
;
; OPTIONAL INPUTS:
;	templates - SEDs to read (default to sirtf.templates)
;	
; OUTPUTS:
;	zarray - redshift array (interval at which the models were
;                redshifted) 
;
; COMMON BLOCKS:
;	sirtf_simulations
;
; PROCEDURES USED:
;	CMRESTORE
;
; MODIFICATION HISTORY:
;	J. Moustakas, 2001 Sep 20, U of A
;-

    common sirtf_simulations

    if not keyword_set(templates) then templates = sirtf.templates
    path = filepath('',root_dir=getenv('SIRTFZ_DIR'),subdirectory='lib')

    datafile = path+strlowcase(templates)+'.idlsave'
    cmrestore, datafile, modelgrid, /quiet

; redshift array
    
    obands = filter_match(filters,sirtf.bandcube.bandnames)
    filterflux = modelgrid.filterflux[*,*,obands]
    colorflux = modelgrid.colorflux[*,*,obands]
    kcorrection = modelgrid.kcorrection[*,*,obands]
    zarray = modelgrid.zarray

return
end
