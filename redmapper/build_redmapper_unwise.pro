pro build_redmapper_unwise, query=query, parse=parse
; jm13mar28siena - match the REDMAPPER/v5.10 catalog to
; Dustin's unWISE catalogs

    path = '/global/u2/i/ioannis/'
;   path = redmapper_path(version=ver)
    unwisepath = '/global/project/projectdirs/desi/users/dstn/sdssv4-pobj/301/'
    
    cat = mrdfits(path+'dr8_run_redmapper_'+ver+$
      '_lgt5_catalog_members.fits.gz',1)
    ngal = n_elements(cat)

; extract     
    photoid_extract, cat.photoid, run, rerun, camcol, field, id

    allfile = unwisepath+strtrim(run,2)+'/'+strtrim(camcol,2)+'/'+$
      'photoWiseForced-'+string(run,format='(I6.6)')+$
      '-'+strtrim(camcol,2)+'-'+string(field,format='(I4.4)')+'.fits'
    ufile = allfile[uniq(allfile,sort(allfile))]
    nfile = n_elements(ufile)
    
    for ii = 0L, ufile-1 do begin

stop       
;      if n_elements(out) eq 0L then out = 
    endfor

return
end    
