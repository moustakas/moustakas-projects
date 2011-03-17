pro zgalaxy, catalog, interp=interp, pdz=pdz, zmin=zmin, zmax=zmax, templates=templates
;+
; NAME:
;	ZGALAXY
;
; PURPOSE:
;	Wrapper to extract the photometric redshift of a galaxy.
;
; CALLING SEQUENCE:
;	zgalaxy, catalog, [interp=], [pdz=], [zmin=], [zmax=], [templates=]
;
; INPUTS:
;	catalog  - root name of the photometric catalog [the assumed
;                  extension is .cat]
;
; OPTIONAL INPUTS:
;	interp    - number of linear interpolation points between SED
;                   type [default 0]
;	pdz       - confidence interval [0,1] within which to compute
;                   the redshift uncertainty [default 0.95]
;	zmin/zmax - restrict the redshift search window to [zmin,zmax]
;	templates - SEDs with which to compute the photometric
;                   redshift [default sirtf.templates]
;
; KEYWORD PARAMETERS:
;
; OUTPUTS:
;	An ASCII file (catalog.dat) and an IDL save set
;	(catalog.idlsave) are written to the results subdirectory.
;	The IDL save set contains all the information in the ASCII
;	file in addition to the p(z) for each galaxy. 
;
; COMMON BLOCKS:
;	sirtf_simulations
;
; PROCEDURES USED:
;	READ_CATALOG, FILTER_MATCH, READ_MODEL_GRIDS, PHOTOZ,
;	GET_ELEMENT, CMSAVE
;
; MODIFICATION HISTORY:
;	J. Moustakas, 2001 May 18, U of A, written
;	jm01aug2uofa, documented
;-

    common sirtf_simulations

    if n_params() eq 0L then begin
       print, 'Syntax - zgalaxy, catalog, [interp=], [pdz=], [zmin=], [zmax=], [templates=]'
       return
    endif

; read in the photometric catalog and information file

    read_catalog, catalog+'.cat', oflux, filters, nsources, catinfo

; match the observed filters to the filter database
    
    obands = filter_match(filters,sirtf.bandcube.bandnames)
    nbands = n_elements(obands)
    
; restore the model grids, restricting the redshift range if specified

    if not keyword_set(templates) then templates = sirtf.templates
    read_model_grids, filters, zarray, mflux, templates=templates

    if keyword_set(zmin) then get_element, zarray, zmin, zstart else zstart = 0L
    if keyword_set(zmax) then get_element, zarray, zmax, zend else zend = size(zarray,/n_elements)-1L

    zarray = zarray[zstart:zend]
    mflux = double(mflux[*,zstart:zend,*])

    mfluxsz = size(mflux,/dimension)
    nseds = mfluxsz[0]
    nz = mfluxsz[1]

; interpolate between SEDs
    
    if n_elements(interp) then interp = round(interp) else interp = 0

    npts = interp*(nseds-1)+nseds                    ; number of points in the interpolated flux
    intgrid = (findgen(npts)/(npts-1.0))*(nseds-1.0) ; fractional indices of interpolation
    
    mflux = interpolate(mflux,intgrid,lindgen(nz),lindgen(nbands),/grid)

; generate an array of structures for the results

    template = {zphot_mode  : fltarr(1), $
                zphot_mean  : fltarr(1), $
                zphot_median: fltarr(1), $
                sigmaz      : fltarr(2), $
                chi2_nu     : fltarr(1), $
                zpdf        : dblarr(nz),$
                constant    : fltarr(1), $
                zindx       : lonarr(1), $
                tindx       : lonarr(1)}
    pzs = replicate(template,nsources)    

; write the results to a file

    rpath = filepath('',root_dir=getenv('SIRTFZ_DIR'),subdirectory='results')
    openw, lun, rpath+catalog+'.dat', /get_lun
    printf, lun, '# SIRTFz analysis of '+catalog+' catalog on '+strmid(systime(),4,20)
    printf, lun, '# Source, [S/N] > 3, zphot (mode), zphot (mean), zphot (median), -sigma(z), +sigma(z), '+$
      'chi2_nu, constant, type, zspec'

    for i = 0L, nsources-1L do begin ; loop on each galaxy in the catalog

       print, format='("Analyzing galaxy ",i0," of ",i0,".",a1,$)', i+1, nsources, string(13b)

       flux = reform(oflux[0,*,i])
       ferror = reform(oflux[1,*,i])
       snhigh = where(flux/ferror gt 10.0,sncount) ; [S/N]

       photoz, flux, ferror, mflux, zarray, pdz=pdz, template
       pzs[i] = template

; write out

       printf, lun, i+1L, sncount, pzs[i].zphot_mode, pzs[i].zphot_mean, pzs[i].zphot_median, $
         pzs[i].sigmaz, pzs[i].chi2_nu, pzs[i].constant, pzs[i].tindx, catinfo.zspec[i], $
         format='(I5,I6,6G15.7,G14.6,F7.2,F14.6)'

    endfor
    print
    free_lun, lun

; save the data structure

    print, 'Saving the data structure.'
    cmsave, filename=rpath+catalog+'.idlsave', pzs, zarray
    
return
end
