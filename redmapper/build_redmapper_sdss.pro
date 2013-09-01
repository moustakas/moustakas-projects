pro build_redmapper_sdss
; jm13mar28siena - build the SDSS catalog for the sample using
; CasJobs; be sure to write the output file from CasJobs as
; 'redmapper_'+ver+'_sdss.fits' 

    path = redmapper_path(version=ver)
    cat = mrdfits(path+'dr8_run_redmapper_'+ver+$
      '_lgt20_catalog_members.fits.gz',1)
    ngal = n_elements(cat)

    out = struct_addtags(replicate({redmapper_id: 0L, casid: long64(0)},ngal),cat)
    out.redmapper_id = lindgen(ngal)
    out.casid = photoid2casid(cat.photoid)

    outfile = '~/tmp/redmapper_'+ver+'_casjobs.dat'
;   outfile = path+'redmapper_'+ver+'_casjobs.dat'
    openw, lun, outfile, /get_lun
    printf, lun, '# redmapper_id casid ra dec'
    struct_print, struct_trimtags(out,select=['redmapper_id','casid','ra','dec']), $
      lun=lun, ddigit=12, /no_head
    free_lun, lun

return
end
    
