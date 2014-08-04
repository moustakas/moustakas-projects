pro redmapper_photoid_extract
; jm14aug01siena - parse photoid so that I can run
; build_redmapper_unwise on NERSC without running into memory issues

    path = redmapper_path(version=ver)
    cat = mrdfits(path+'dr8_run_redmapper_'+ver+$
      '_lgt5_catalog_members.fits.gz',1,columns='photoid')
    ngal = n_elements(cat)
    out = replicate({run: 0L, rerun: 0, camcol: 0, field: 0, id: 0L, file: ''},ngal)

    photoid_extract, cat.photoid, run, rerun, camcol, field, id
;   out.photoid = cat.photoid
    out.run = run
    out.rerun = rerun
    out.camcol = camcol
    out.field = field
    out.id = id
    out.file = strtrim(run,2)+'/'+strtrim(camcol,2)+'/'+$
      'photoWiseForced-'+string(run,format='(I6.6)')+$
      '-'+strtrim(camcol,2)+'-'+string(field,format='(I4.4)')+'.fits'

    im_mwrfits, out, path+'redmapper_'+ver+'_photoid.fits', /clobber

stop    
    
return
end
    
