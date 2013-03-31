pro redmapper_write_unique

    path = getenv('REDMAPPER_DATA')
    all = mrdfits(path+'/dr8_run_redmapper_v5.2_lgt20_catalog_members.fits.gz',1)
    nall = n_elements(all)

    

    
; write out the catalog of unique objecs
    indx = lindgen(nall)
    uu = uniq(all.photoid,sort(all.photoid))
    uindx = indx[uu]

    outfile = path+'/dr8_run_redmapper_v5.2_lgt20_catalog_members_uniq.fits'
    im_mwrfits, all[uu], outfile, /clobber

stop    
    
return
end
    

    
