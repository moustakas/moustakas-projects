pro write_redmapper_casjobs_input
; jm13mar28siena
;
    ver = 'v5.2'
    outpath = getenv('REDMAPPER_DATA')+'/catalogs/'
    galexpath = outpath

    cat = mrdfits(outpath+'dr8_run_redmapper_'+ver+$
      '_lgt20_catalog_members.fits.gz',1)
    ngal = n_elements(cat)

    out = struct_addtags(replicate({redmapper_id: 0L, casid: long64(0)},ngal),cat)
    out.redmapper_id = lindgen(ngal)
    out.casid = photoid2casid(cat.photoid)

    outfile = galexpath+'redmapper_'+ver+'_casjobs.dat'
    openw, lun, outfile, /get_lun
    printf, lun, '# redmapper_id casid ra dec'
    struct_print, struct_trimtags(out,select=['redmapper_id','casid','ra','dec']), $
      lun=lun, ddigit=12, /no_head
    free_lun, lun

return
end
    
