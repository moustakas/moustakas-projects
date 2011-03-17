pro make_catalog, filters, fref, ncatalog=ncatalog, templates=templates, $
        fmin=fmin, fmax=fmax, zmin=zmin, zmax=zmax, alpha=alpha, mag=mag
;+
; NAME:
;	MAKE_CATALOG
;
; PURPOSE:
;	Generate a photometric catalog with Gaussian noise.
;
; CALLING SEQUENCE:
;	make_catalog, 
;
; INPUTS:
;	filters  - string array specifying the filters in which to
;                  generate a catalog
;	fref     - reference filter (zero-indexed)
;
; OPTIONAL INPUTS:
;	ncatalog   - number of simulated galaxies to generate (default 1000)
;	templates  - SED templates to use (Devriendt, CWW)
;	fmin       - minimum apparent flux in fref in mJy to accept in
;                    the catalog (default: 5-sigma in the fiducial band)
;	fmax       - maximum apparent flux in fref in mJy to accept in
;                    the catalog (default: 10 Jy)
;	zmin       - minimum redshift
;	zmax       - maximum redshift
;	
; KEYWORD PARAMETERS:
;	mag      - fmin and fmax are in apparent magnitudes
;
; OUTPUTS:
;
;
; OPTIONAL OUTPUTS:
;
;
; COMMON BLOCKS:
;	sirtf_simulations
;
; COMMENTS:
;	The mag keyword has not been implemented yet.
;
; PROCEDURES USED:
;
;
; MODIFICATION HISTORY:
;	J. Moustakas, 2001 August 3, U of A, written
;-

    common sirtf_simulations
    common cosmology

; ----------------------------------------------------------------------

    filters = *sirtf.filters
    fref = 4L

;   filters = sirtf.bandcube[14:17].bandnames ; WFPC2 filters
;   fref = 3L

; ----------------------------------------------------------------------

;   if n_params() lt 2L then message, 'Please specify the filters and the reference filter.'
    if not keyword_set(ncatalog) then ncatalog = 10L
    if not keyword_set(templates) then templates = sirtf.templates

    obands = filter_match(filters,sirtf.bandcube.bandnames)
    nbands = size(filters,/n_elements)

    if not keyword_set(fmin) then fmin = 5.0*sirtf.bandcube[obands[fref]].flimit ; 5-sigma
    if not keyword_set(fmax) then fmax = 10.0*1E3                                ; 10,000 mJy
    if not keyword_set(alpha) then alpha = sirtf.evolution.alpha

    read_model_grids, filters, mflux, zarray, templates=templates

    if keyword_set(zmin) then get_element, zarray, zmin, zstart else zstart = 0L
    if keyword_set(zmax) then get_element, zarray, zmax, zend else zend = size(zarray,/n_elements)-1L

    zarray = zarray[zstart:zend]
    nz = n_elements(zarray)
    mflux = mflux[*,zstart:zend,*]

    path = filepath('',root_dir=getenv('SIRTFZ_DIR'),subdirectory='catalogs')
    file = 'simcat.cat'
    print, 'Opening '+path+file
    
    openw, lun, path+file, /get_lun
    printf, lun, '# '+strmid(systime(),4,20)+' simulated catalog: ncat = '+strn(ncatalog)
    header = ['# '+filters[0]+',',filters[1:nbands-2L]+',',filters[nbands-1L]]
    printf, lun, header

    nsources = 0L
    while nsources lt ncatalog do begin ; generate the catalog

       z = (zarray[floor(randomu(seed,1)*nz)])[0] ; random redshift (could put in a z-distribution)
       l_evol = (1.0+z)^alpha                     ; luminosity evolution

       get_element, zarray, z, zindx
       good = where((l_evol*mflux[*,zindx,fref] gt fmin) and $ ; above the limiting flux
                    (l_evol*mflux[*,zindx,fref] lt fmax),ngood)

       if ngood ne 0L then begin

          sindx = floor(randomu(seed,1)*ngood) ; random SED type

          flux = reform(l_evol*mflux[good[sindx],zindx,*])
          ferror = sqrt(sirtf.bandcube[obands].flimit^2 + (0.05*flux)^2) ; flux error
          dflux = flux + randomn(seed,nbands)*ferror                     ; Gaussian-perturbed flux

          print, format='("Writing source ",I0," of ",I0,".",A1,$)', nsources, ncatalog, string(13b)
          printf, lun, dflux, ferror, z, format='('+strn(nbands)+'E15.5,'+strn(nbands)+'E15.5,3x,D0.0)'
          
          nsources = nsources + 1L

       endif
       
    endwhile
    print
    free_lun, lun

return
end
