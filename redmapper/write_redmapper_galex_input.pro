pro write_redmapper_galex_input, gr=gr
; jm13jan19siena - based on WRITE_VAGC_GALEX_CASJOBS_INPUT; modified
; to run on all of REDMAPPER/DR9
;
    if (n_elements(gr) eq 0) then gr = 'gr6'
    ver = 'v5.2'
    
    outpath = getenv('REDMAPPER_DATA')+'/catalogs/'
    galexpath = outpath

    cat = hogg_mrdfits(outpath+'dr8_run_redmapper_'+ver+$
      '_lgt20_catalog_members.fits.gz',1,$
      columns=['ra','dec'],nrow=100000L)
    ngal = n_elements(cat)

    out = struct_addtags(replicate({redmapper_id: 0L},ngal),cat)
    out.redmapper_id = lindgen(ngal)

    outfile = galexpath+'redmapper_'+ver+'_galex_'+gr+'_casjobs.dat'
    openw, lun, outfile, /get_lun
    printf, lun, '# redmapper_id ra dec'
    struct_print, out, lun=lun, ddigit=12, /no_head
    free_lun, lun

return
end
    
