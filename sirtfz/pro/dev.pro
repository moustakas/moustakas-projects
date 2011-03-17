pro dev, alpha=alpha, sigma=sigma
; jm01aug6uofa
; write out a text file indicating if an SED is a sigma detection in
; the user-selected filters
    
    common sirtf_simulations, sirtf
    if not keyword_set(sirtf) then sirtfz_initialize, templates='DEVRIENDT'
    if not keyword_set(sigma) then sigma=3.0
    
    rpath = filepath('',root_dir=getenv('SIRTFZ_DIR'),subdirectory='results')
    filters = *sirtf.filters
    obands = filter_match(filters,sirtf.bandcube.bandnames)
    nbands = size(filters,/n_elements)
    
    yesno = strarr(nbands)
    
    read_model_grids, filters, mflux, zarray, templates=sirtf.templates

    dz = 1.0
    z = (findgen(5/dz)+1.0)*dz
    nz = size(z,/n_elements)

    get_element, zarray, z, zindx

;   sindx = [0,8,14]
    sindx = lindgen(n_elements(sirtf.sedcube))

    if keyword_set(alpha) then fname = 'choose_filters_evol.dat' else fname = 'choose_filters.dat'
    openw, lun, rpath+fname, /get_lun
    printf, lun, '# ', sirtf.bandcube[obands].bandnames
    if keyword_set(alpha) then printf, lun, '#  Evolution: '+strn(alpha,length=3) else $
      printf, lun, '#  Evolution: None'
    printf, lun, ' '
    for i = 0L, n_elements(sindx)-1L do begin

       for j = 0L, nz-1L do begin

          flux = mflux[sindx[i],zindx[j],*]
          if keyword_set(alpha) then flux = flux * (1.0+z[j])^alpha

          yes = where(flux gt sigma*sirtf.bandcube[obands].flimit eq 1B,nryes,comp=no,ncomp=nrno)
          if nryes ne 0L then yesno[yes] = 'YES'
          if nrno ne 0L then yesno[no] = 'NO'

;         printf, lun, sindx[i], sirtf.sedcube[sindx[i]].lum60_sun, z[j], reform(mflux[sindx[i],zindx[j],*]), $
;           yesno, format='(I5,F10.5,F10.5,'+strn(nbands)+'G15.5,'+strn(nbands)+'A7)'
          printf, lun, sindx[i], sirtf.sedcube[sindx[i]].lum60_sun, z[j], yesno, $
            format='(I2,F10.5,F10.5,'+strn(nbands)+'A7)'
          
       endfor

       printf, lun, ' '

    endfor
    free_lun, lun

return
end
